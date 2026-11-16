#include "phz.hpp"

namespace phz {
static int nComp(const VI& sh) {
    return std::accumulate(sh.begin(), sh.end(), 1, std::multiplies<int>());
}

void System::initMixed() {
    req(mdl->dim == 2 || mdl->dim == 3, "Mixed flow requires 2D or 3D");
    req(bool(mdl->pde), "Mixed flow requires a PDE");
    auto& pd = *mdl->pde;
    std::string vel, press;
    for (auto s : pd.fields) {
        auto sh = pd.ns.at(s)->sh;
        if (sh == VI{mdl->dim}) {
            req(vel.empty(), "One velocity field required");
            vel = s;
        } else if (sh.empty()) {
            req(press.empty(), "One pressure field required");
            press = s;
        } else
            throw Err("Expected velocity vector and scalar pressure");
    }
    req(!vel.empty() && !press.empty() && pd.fields.size() == 2,
        "Mixed flow requires velocity and pressure");
    mdl->fields = {vel, press};
    mesh = std::make_unique<Mesh>(*mdl);
    req(mesh->geo == mesh->dim, "Mixed flow requires full-dimensional mesh");
    int d = mdl->dim, q = mdl->disc.value("q", 4);
    auto u = std::make_shared<Space>(*mesh, 2, q), p = std::make_shared<Space>(*mesh, 1, q);
    spaces.assign(d, u);
    spaces.push_back(p);
    for (int k = 0; k < d; ++k)
        compNames.push_back(vel + "_" + std::to_string(k));
    compNames.push_back(press);
    offset = {0};
    for (auto s : spaces)
        offset.push_back(offset.back() + s->n);
    size = offset.back();
    n = u->n;
    nf = d + 1;
    req(size <= 30000, "Mixed state exceeds 30000 coefficients");
    x = u->x;
    w = u->w;
    b = u->b;
    g = u->g;
    mass = u->mass;
    stiff = u->stiff;
    massv = mass * V::Ones(n);
    std::map<std::string, E> substitutions;
    J src = mdl->input.value("src", J::object()), co = mdl->input.value("coef", J::object());
    req(src.size() == co.size(), "All source coefficients require expressions");
    Parser parser(pd.ns, d);
    for (auto it = co.begin(); it != co.end(); ++it) {
        E e = parser.parse(src.at(it.key()));
        req(e->sh == it.value().get<VI>() &&
                !has(e, {"fld", "coef", "lap", "div", "grad", "intg", "dt"}),
            "Source must be local prescribed data");
        auto values = components(e, d), symbols = components(pd.ns.at(it.key()), d);
        for (int i = 0; i < values.size(); ++i) {
            substitutions[symbols[i]->name] = values[i];
            for (int k = 0; k < d; ++k)
                substitutions[symbols[i]->name + "$" + std::to_string(k)] = dx(values[i], k);
        }
    }
    auto strong = components(pd.eq.at(vel), d);
    E continuity = components(pd.eq.at(press), d)[0];
    std::map<std::string, double> zeroCtx = mdl->par;
    for (int k = 0; k < 3; ++k)
        zeroCtx[std::string(1, "xyz"[k])] = 0;
    zeroCtx["t"] = 0;
    for (auto s : compNames) {
        zeroCtx[s] = 0;
        for (int k = 0; k < 4; ++k)
            zeroCtx[s + "$" + std::to_string(k)] = 0;
    }
    for (int f = 0; f < nf; ++f) {
        req(eval(diff(continuity, compNames[f]), mdl->par) == 0,
            "Continuity must equal div(velocity)");
        for (int k = 0; k < d; ++k)
            req(eval(diff(continuity, compNames[f] + "$" + std::to_string(k)), mdl->par) ==
                    (f == k ? 1 : 0),
                "Continuity must equal div(velocity)");
    }
    req(eval(continuity, zeroCtx) == 0, "Continuity has extra source");
    int time = -1;
    for (int f = 0; f < d; ++f) {
        E e = subst(strong[f], substitutions);
        double a = eval(diff(e, compNames[f] + "$3"), mdl->par);
        req(a == 0 || a == 1, "Momentum requires unit time derivative");
        if (time < 0)
            time = int(a);
        req(time == int(a), "Inconsistent velocity time derivatives");
        req(eval(diff(e, press), mdl->par) == 0, "Extra pressure dependence");
        for (int k = 0; k < d; ++k)
            req(eval(diff(e, press + "$" + std::to_string(k)), mdl->par) == (f == k ? 1 : 0),
                "Momentum requires +grad(pressure) in residual");
    }
    transient = time == 1;
    for (auto field : mdl->fields) {
        E v = sym("v_" + field, pd.ns.at(field)->sh, "test");
        auto wk = weak(pd.eq.at(field), v, d);
        E ex = subst(components(wk.first, d).at(0), substitutions);
        std::function<void(E)> check = [&](E e) {
            req(std::count(e->name.begin(), e->name.end(), '$') <= 1 &&
                    !has(e, {"intg", "lap", "div", "grad", "dt"}),
                "Mixed weak form requires local values and first derivatives");
            for (auto a : e->a)
                check(a);
        };
        check(ex);
        VI rows = field == vel ? VI{} : VI{d};
        if (field == vel)
            for (int f = 0; f < d; ++f)
                rows.push_back(f);
        for (int f : rows)
            for (int k = 0; k <= d; ++k) {
                std::string test = "v_" + compNames[f] + (k ? "$" + std::to_string(k - 1) : "");
                E r = diff(ex, test);
                if (zero(r))
                    continue;
                terms.push_back({f, k, r});
                for (int c = 0; c < nf; ++c)
                    for (int l = 0; l <= d + 1; ++l) {
                        std::string symbol =
                            compNames[c] + (l ? "$" + std::to_string(l == d + 1 ? 3 : l - 1) : "");
                        E j = diff(r, symbol);
                        if (!zero(j))
                            jterms.push_back({f, k, c, l, j});
                    }
            }
    }
    if (!mdl->bc.contains(vel))
        mdl->bc[vel] = J{{"type", "dirichlet"}, {"val", 0}};
    initBoundary();
    std::vector<Tr> dt;
    for (int f = 0; f < d; ++f) {
        SM wg = u->g[f];
        for (int k = 0; k < wg.outerSize(); ++k)
            for (SM::InnerIterator it(wg, k); it; ++it)
                it.valueRef() *= w[it.row()];
        SM block = p->b.transpose() * wg;
        for (int k = 0; k < block.outerSize(); ++k)
            for (SM::InnerIterator it(block, k); it; ++it)
                dt.emplace_back(it.row(), offset[f] + int(it.col()), it.value());
    }
    div.resize(p->n, offset[d]);
    div.setFromTriplets(dt.begin(), dt.end());
    gauge = V::Zero(size);
    gauge.tail(p->n) = p->mass * V::Ones(p->n) / w.sum();
    VI freeV;
    for (int i : active)
        if (i < offset[d])
            freeV.push_back(i);
    VI pressIds(p->n);
    std::iota(pressIds.begin(), pressIds.end(), 0);
    if (p->n <= 256) {
        M dv = block(div, pressIds, freeV);
        Eigen::FullPivLU<M> rank(dv);
        rank.setThreshold(1e-10);
        req(rank.rank() == p->n - 1, "Additional pressure nullspace: refine mesh");
    }
    u0 = V::Zero(size);
    J init = mdl->input.value("ini", J::object());
    double t0 = mdl->run.value("t", J::array({0, 1}))[0].get<double>();
    for (auto field : mdl->fields) {
        E e = mdl->parse(init.value(field, J(0)));
        VI sh = pd.ns.at(field)->sh;
        if (zero(e) && e->sh.empty()) {
            for (int f = field == vel ? 0 : d; f < (field == vel ? d : d + 1); ++f)
                u0.segment(offset[f], spaces[f]->n).setZero();
        } else {
            req(e->sh == sh && !has(e, {"fld", "coef", "lap", "grad", "div", "dt", "intg"}),
                "Invalid mixed initial data");
            auto a = components(e, d);
            for (int k = 0; k < a.size(); ++k) {
                int f = field == vel ? k : d;
                u0.segment(offset[f], spaces[f]->n) = at(a[k], spaces[f]->p, t0);
            }
        }
    }
    u0 = lift(t0, u0);
    u0.tail(p->n).array() -= gauge.dot(u0);
}

V System::projectVelocity(const V& u) const {
    int d = mdl->dim, np = spaces.back()->n, nv = offset[d];
    VI ids;
    for (int i : active)
        if (i < nv)
            ids.push_back(i);
    VI pr(np);
    std::iota(pr.begin(), pr.end(), 0);
    SM dv = block(div, pr, ids);
    std::vector<Tr> trip;
    std::map<int, int> index;
    for (int i = 0; i < ids.size(); ++i)
        index[ids[i]] = i;
    for (int f = 0; f < d; ++f) {
        auto& m = spaces[f]->mass;
        for (int k = 0; k < m.outerSize(); ++k)
            for (SM::InnerIterator it(m, k); it; ++it)
                if (index.count(offset[f] + int(it.row())) &&
                    index.count(offset[f] + int(it.col())))
                    trip.emplace_back(index[offset[f] + int(it.row())],
                                      index[offset[f] + int(it.col())],
                                      it.value());
    }
    int nf = int(ids.size());
    for (int k = 0; k < dv.outerSize(); ++k)
        for (SM::InnerIterator it(dv, k); it; ++it) {
            trip.emplace_back(nf + it.row(), it.col(), it.value());
            trip.emplace_back(it.col(), nf + it.row(), it.value());
        }
    for (int i = 0; i < np; ++i) {
        trip.emplace_back(nf + i, nf + np, gauge[nv + i]);
        trip.emplace_back(nf + np, nf + i, gauge[nv + i]);
    }
    SM a(nf + np + 1, nf + np + 1);
    a.setFromTriplets(trip.begin(), trip.end());
    V r = V::Zero(a.rows());
    r.segment(nf, np) = -div * u.head(nv);
    V step = solve(a, r), out = u;
    put(out, ids, pick(out, ids) + step.head(nf));
    req((div * out.head(nv)).lpNorm<Eigen::Infinity>() < 1e-8,
        "Boundary flux incompatible with incompressibility");
    return out;
}
}
