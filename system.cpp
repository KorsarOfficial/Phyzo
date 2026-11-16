#include "phz.hpp"

namespace phz {
static bool constantExpr(E e, const Model& m) {
    if (e->op == "num")
        return true;
    if (!e->name.empty())
        return e->name == "pi" || m.par.count(e->name);
    if (has(e, {"lap", "divg", "intg", "grad", "div", "dt"}))
        return false;
    for (auto a : e->a)
        if (!constantExpr(a, m))
            return false;
    return true;
}

static bool polynomial(E e, const Model& m) {
    if (constantExpr(e, m))
        return true;
    if (e->op == "fld")
        return true;
    if (e->op == "add" || e->op == "mul" || e->op == "lap" || e->op == "divg") {
        for (auto a : e->a)
            if (!polynomial(a, m))
                return false;
        return true;
    }
    return e->op == "pow" && e->a[1]->n >= 0 && e->a[1]->n <= 4 &&
           std::floor(e->a[1]->n) == e->a[1]->n && polynomial(e->a[0], m);
}

static SM weighted(const SM& b, const V& w) {
    SM a = b;
    for (int k = 0; k < a.outerSize(); ++k)
        for (SM::InnerIterator it(a, k); it; ++it)
            it.valueRef() *= w[it.row()];
    return a;
}

static void boundaryForm(E e, const Model& m) {
    bool seen = false;
    std::function<void(E)> walk = [&](E v) {
        if (!has(v, {"lap", "divg"}))
            return;
        if (v->op == "add") {
            walk(v->a[0]);
            walk(v->a[1]);
            return;
        }
        if (v->op == "mul") {
            for (int i = 0; i < 2; ++i)
                if (constantExpr(v->a[i], m)) {
                    walk(v->a[1 - i]);
                    return;
                }
            throw Err("Spatially varying outer coefficient is not a conservative flux");
        }
        req(v->op == "lap" || v->op == "divg",
            "Boundary requires a conservative second-order operator");
        for (auto a : v->a)
            req(!has(a, {"lap", "divg"}),
                "Nested operators require separate intermediate-field boundary conditions");
        seen = true;
    };
    walk(e);
    req(seen, "Boundary requires a spatial operator");
}

static bool sameBoundary(E a, E b, const Model& m, const M& point) {
    std::map<std::string, E> values;
    for (auto [name, value] : m.par)
        values[name] = num(value);
    values["pi"] = num(pi);
    for (int i = 0; i < 3; ++i)
        values[std::string(1, "xyz"[i])] = num(i < point.cols() ? point(0, i) : 0);
    a = subst(a, values);
    b = subst(b, values);
    if (str(a) == str(b))
        return true;
    using P = std::vector<double>;
    std::function<P(E)> poly = [&](E e) -> P {
        if (e->op == "num")
            return {e->n};
        if (e->name == "t")
            return {0, 1};
        req(e->op == "add" || e->op == "mul" || e->op == "pow",
            "Cannot establish symbolic boundary agreement");
        P x = poly(e->a[0]);
        if (e->op == "add") {
            P y = poly(e->a[1]);
            x.resize(std::max(x.size(), y.size()));
            for (int i = 0; i < y.size(); ++i)
                x[i] += y[i];
            return x;
        }
        auto product = [](P x, P y) {
            req(x.size() + y.size() <= 32, "Boundary polynomial too large");
            P z(x.size() + y.size() - 1);
            for (int i = 0; i < x.size(); ++i)
                for (int j = 0; j < y.size(); ++j)
                    z[i + j] += x[i] * y[j];
            return z;
        };
        if (e->op == "mul")
            return product(x, poly(e->a[1]));
        double n = e->a[1]->n;
        req(n >= 0 && n <= 12 && n == std::floor(n),
            "Cannot establish symbolic boundary agreement");
        P z{1};
        for (int i = 0; i < n; ++i)
            z = product(z, x);
        return z;
    };
    try {
        P x = poly(a), y = poly(b);
        x.resize(std::max(x.size(), y.size()));
        y.resize(x.size());
        for (int i = 0; i < x.size(); ++i)
            if (std::abs(x[i] - y[i]) > 1e-12 * std::max({1., std::abs(x[i]), std::abs(y[i])}))
                return false;
        return true;
    } catch (const Err&) {
        return false;
    }
}

static std::unique_ptr<LU> factor(const SM& m) {
    auto lu = std::make_unique<LU>();
    lu->compute(m);
    req(lu->info() == Eigen::Success, "Singular mass or constraint matrix");
    return lu;
}

V System::at(E e, const M& pts, double t) const {
    return mdl->at(e, pts, t);
}

System::System(J j) : mdl(std::make_shared<Model>(std::move(j))) {
    mixed = mdl->type == "th";
    if (mixed)
        initMixed();
    else
        initScalar();
}

void System::initScalar() {
    nf = int(mdl->fields.size());
    if (mdl->type == "fe") {
        mesh = std::make_unique<Mesh>(*mdl);
        auto sp = std::make_shared<Space>(*mesh, 1, mdl->disc.value("q", 4));
        spaces.assign(nf, sp);
        x = sp->x;
        w = sp->w;
        b = sp->b;
        g = sp->g;
        mass = sp->mass;
        stiff = sp->stiff;
        n = sp->n;
        one = V::Ones(n);
    } else {
        VI counts = mdl->disc.value("n", VI(mdl->dim, 4));
        n = 1;
        for (int k : counts)
            n *= k;
        V len(mdl->dim);
        for (int i = 0; i < mdl->dim; ++i)
            len[i] = mdl->par.at("Lx" + std::to_string(i + 1));
        for (int i = 0; i < n; ++i) {
            int r = i;
            VI a(mdl->dim);
            for (int k = mdl->dim - 1; k >= 0; --k) {
                a[k] = r % counts[k];
                r /= counts[k];
            }
            modes.push_back(a);
        }
        int q = mdl->disc.value("q",
                                std::max(16, 3 * *std::max_element(counts.begin(), counts.end()))),
            nq = 1;
        for (int k = 0; k < mdl->dim; ++k)
            nq *= q;
        req(int64_t(nq) * n <= 8000000, "Spectral quadrature memory limit");
        auto [z, wt] = gauss(q);
        x.resize(nq, mdl->dim);
        w.resize(nq);
        for (int i = 0; i < nq; ++i) {
            int r = i;
            double ww = 1;
            for (int k = mdl->dim - 1; k >= 0; --k) {
                int v = r % q;
                r /= q;
                x(i, k) = z[v] * len[k];
                ww *= wt[v] * len[k];
            }
            w[i] = ww;
        }
        M aB(nq, n);
        std::vector<M> aG(mdl->dim, M(nq, n));
        std::vector<Tr> mt, kt;
        for (int k = 0; k < n; ++k) {
            double nm = 1, ev = 0;
            for (int d = 0; d < mdl->dim; ++d) {
                nm *= len[d] * (modes[k][d] ? 0.5 : 1);
                ev += std::pow(pi * modes[k][d] / len[d], 2);
            }
            mt.emplace_back(k, k, nm);
            kt.emplace_back(k, k, nm * ev);
            for (int i = 0; i < nq; ++i) {
                double v = 1;
                for (int d = 0; d < mdl->dim; ++d)
                    v *= std::cos(pi * modes[k][d] * x(i, d) / len[d]);
                aB(i, k) = v;
                for (int ax = 0; ax < mdl->dim; ++ax) {
                    v = 1;
                    for (int d = 0; d < mdl->dim; ++d) {
                        double f = pi * modes[k][d] / len[d];
                        v *= ax == d ? -f * std::sin(f * x(i, d)) : std::cos(f * x(i, d));
                    }
                    aG[ax](i, k) = v;
                }
            }
        }
        b = sparse(aB);
        for (auto a : aG)
            g.push_back(sparse(a));
        mass.resize(n, n);
        mass.setFromTriplets(mt.begin(), mt.end());
        stiff.resize(n, n);
        stiff.setFromTriplets(kt.begin(), kt.end());
        one = V::Zero(n);
        one[0] = 1;
    }
    size = n * nf;
    for (int i = 0; i <= nf; ++i)
        offset.push_back(i * n);
    compNames = mdl->fields;
    massv = mass * one;
    massLU = factor(mass);
    initBoundary();
    u0 = V::Zero(size);
    J ini = mdl->input.value("ini", J::object());
    for (int f = 0; f < nf; ++f) {
        J in = ini.value(mdl->fields[f], J(0));
        V u;
        if (in.is_array()) {
            u = vec(in);
            req(u.size() == n, "Initial coefficient count mismatch");
        } else {
            E e = mdl->parse(in);
            req(!has(e, {"fld", "lap", "divg", "intg"}), "Initial data must be prescribed");
            u = mdl->type == "fe"
                    ? at(e, mesh->p, mdl->run.value("t", J::array({0, 1}))[0].get<double>())
                    : project(at(e, x, mdl->run.value("t", J::array({0, 1}))[0].get<double>()));
        }
        u0.segment(f * n, n) = u;
    }
    u0 = lift(mdl->run.value("t", J::array({0, 1}))[0].get<double>(), u0);
    if (!mdl->alg.empty()) {
        std::vector<Tr> t;
        for (int k = 0; k < stiff.outerSize(); ++k)
            for (SM::InnerIterator it(stiff, k); it; ++it)
                t.emplace_back(it.row(), it.col(), it.value());
        for (int i = 0; i < n; ++i) {
            t.emplace_back(i, n, massv[i]);
            t.emplace_back(n, i, massv[i]);
        }
        SM a(n + 1, n + 1);
        a.setFromTriplets(t.begin(), t.end());
        poissonLU = factor(a);
        for (auto a : mdl->alg)
            req(mdl->constant(a.value("a", J(1))) > 0, "Poisson coefficient must be positive");
    }
    exactSp = mdl->type == "sp" && size <= 128 && mdl->alg.empty();
    for (auto [s, e] : mdl->eq)
        exactSp &= polynomial(e, *mdl);
    if (exactSp)
        exportCpp();
    if (mdl->type == "fe" && mdl->alg.empty()) {
        try {
            for (int f = 0; f < nf; ++f) {
                E v = sym("v_" + compNames[f], {}, "test");
                auto wk = weak(mdl->eq.at(compNames[f]), v, mesh->geo);
                E e = components(wk.first, mesh->geo).at(0);
                req(!has(e, {"intg", "lap", "divg", "dt"}), "Nonlocal weak expression");
                std::function<void(E)> check = [&](E a) {
                    req(std::count(a->name.begin(), a->name.end(), '$') <= 1,
                        "Higher derivative weak expression");
                    for (auto c : a->a)
                        check(c);
                };
                check(e);
                for (int k = 0; k <= mesh->geo; ++k) {
                    std::string test = "v_" + compNames[f] + (k ? "$" + std::to_string(k - 1) : "");
                    E r = diff(e, test);
                    if (zero(r))
                        continue;
                    terms.push_back({f, k, r});
                    for (int c = 0; c < nf; ++c)
                        for (int l = 0; l <= mesh->geo; ++l) {
                            E j = diff(r, compNames[c] + (l ? "$" + std::to_string(l - 1) : ""));
                            if (!zero(j))
                                jterms.push_back({f, k, c, l, j});
                        }
                }
            }
        } catch (const Err&) {
            terms.clear();
            jterms.clear();
        }
    }
}

void System::initBoundary() {
    int nc = int(compNames.size());
    bcs.resize(nc);
    active.clear();
    freeLU.clear();
    for (int f = 0; f < nc; ++f) {
        int nn = mixed ? spaces[f]->n : n;
        std::string name = mixed ? mdl->fields[f < mdl->dim ? 0 : 1] : compNames[f];
        J patches = mdl->bc.value(name, J::array());
        if (patches.is_object())
            patches = J::array({patches});
        req(patches.is_array(), "Boundary must be object or list");
        Boundary& bc = bcs[f];
        std::map<int, E> fixed;
        std::set<VI> occupied;
        for (auto patch : patches) {
            keys(patch, "type on val coef", "boundary");
            std::string ty = patch.at("type"), tag = patch.value("on", "all");
            req(ty == "dirichlet" || ty == "neumann" || ty == "robin", "Invalid boundary type");
            E val = mdl->parse(patch.value("val", J(0)));
            req(!has(val, {"fld", "grad", "lap", "divg", "intg", "dt"}),
                "Boundary data must be prescribed");
            if (mixed) {
                if (f == mdl->dim)
                    throw Err("Pressure boundary conditions are not supported");
                req(ty == "dirichlet" && tag == "all",
                    "Mixed velocity requires full Dirichlet boundary");
                if (!zero(val)) {
                    req(val->sh == VI{mdl->dim}, "Velocity boundary shape mismatch");
                    val = components(val, mdl->dim)[f];
                } else
                    val = num(0);
            }
            if (!mesh) {
                req(ty == "neumann" && zero(val) && tag == "all",
                    "Cosine basis requires homogeneous Neumann boundary");
                continue;
            }
            req(mesh->tags.count(tag), "Unknown boundary tag: " + tag);
            IM faces = mesh->tags.at(tag);
            req(faces.rows() > 0 || (ty == "neumann" && zero(val) && tag == "all"),
                "Empty boundary tag");
            for (int i = 0; i < faces.rows(); ++i) {
                VI a(mesh->dim);
                for (int k = 0; k < mesh->dim; ++k)
                    a[k] = faces(i, k);
                req(occupied.insert(a).second, "Overlapping boundary facets");
            }
            if (!mixed && !(ty == "neumann" && zero(val)))
                boundaryForm(mdl->eq.at(compNames[f]), *mdl);
            if (ty == "dirichlet") {
                VI ids = spaces[f]->boundary(*mesh, tag);
                for (int id : ids) {
                    if (fixed.count(id)) {
                        M pt = spaces[f]->p.row(id);
                        req(sameBoundary(fixed[id], val, *mdl, pt),
                            "Conflicting or symbolically unverified Dirichlet data");
                    }
                    fixed[id] = val;
                }
                continue;
            }
            if (ty == "neumann" && zero(val))
                continue;
            req(!mixed, "Mixed natural boundaries unsupported");
            req(has(mdl->eq.at(compNames[f]), {"lap", "divg"}),
                "Nonzero flux requires conservative PDE");
            auto [l, wt] = simplexQuad(mesh->dim - 1, mdl->disc.value("q", 4));
            Boundary::Patch pp;
            int nq = int(l.rows());
            pp.x.resize(nq * faces.rows(), mesh->geo);
            pp.w.resize(pp.x.rows());
            std::vector<Tr> trip;
            for (int e = 0; e < faces.rows(); ++e) {
                M pts(mesh->dim, mesh->geo);
                for (int k = 0; k < mesh->dim; ++k)
                    pts.row(k) = mesh->p.row(faces(e, k));
                double vol = 1;
                if (mesh->dim > 1) {
                    M ed = (pts.bottomRows(mesh->dim - 1).rowwise() - pts.row(0)).transpose();
                    vol = std::sqrt((ed.transpose() * ed).determinant()) / (mesh->dim == 3 ? 2 : 1);
                }
                pp.x.middleRows(e * nq, nq) = l * pts;
                pp.w.segment(e * nq, nq) = vol * wt;
                for (int i = 0; i < nq; ++i)
                    for (int k = 0; k < mesh->dim; ++k)
                        trip.emplace_back(e * nq + i, faces(e, k), l(i, k));
            }
            pp.b.resize(pp.x.rows(), n);
            pp.b.setFromTriplets(trip.begin(), trip.end());
            pp.val = val;
            pp.coef = ty == "robin" ? mdl->parse(patch.at("coef")) : E();
            if (pp.coef)
                req((at(pp.coef, pp.x, 0).array() >= 0).all(),
                    "Robin coefficient must be nonnegative");
            bc.natural.push_back(std::move(pp));
        }
        for (auto [i, e] : fixed) {
            bc.fix.push_back(i);
            bc.val.push_back(e);
            bc.dt.push_back(diff(e, "t"));
        }
        for (int i = 0; i < nn; ++i)
            if (!fixed.count(i)) {
                bc.free.push_back(i);
                active.push_back(offset[f] + i);
            }
        SM mm = mixed ? spaces[f]->mass : mass;
        if (!bc.free.empty())
            freeLU.push_back(factor(block(mm, bc.free, bc.free)));
        else
            freeLU.push_back(nullptr);
    }
}

V System::lift(double t, const V& u) const {
    req(u.size() == size && u.allFinite(), "State must be a finite vector of the declared size");
    V a = u;
    for (int f = 0; f < bcs.size(); ++f) {
        auto& bc = bcs[f];
        for (int k = 0; k < bc.fix.size(); ++k) {
            M pt = spaces[f]->p.row(bc.fix[k]);
            a[offset[f] + bc.fix[k]] = at(bc.val[k], pt, t)[0];
        }
    }
    return a;
}

V System::boundaryDt(int f, double t) const {
    int nn = mixed ? spaces[f]->n : n;
    V a = V::Zero(nn);
    auto& bc = bcs[f];
    for (int k = 0; k < bc.fix.size(); ++k) {
        M pt = spaces[f]->p.row(bc.fix[k]);
        a[bc.fix[k]] = at(bc.dt[k], pt, t)[0];
    }
    return a;
}

V System::boundaryLoad(int f, double t, const V& u) const {
    V a = V::Zero(n);
    for (auto& p : bcs[f].natural) {
        V val = at(p.val, p.x, t);
        if (p.coef) {
            V k = at(p.coef, p.x, t);
            req((k.array() >= 0).all(), "Negative Robin coefficient");
            val = k.array() * (p.b * u - val).array();
        }
        a -= p.b.transpose() * (p.w.array() * val.array()).matrix();
    }
    return a;
}

V System::project(const V& q) const {
    req(q.size() == w.size(), "Quadrature size mismatch");
    V r = massLU->solve(b.transpose() * (w.array() * q.array()).matrix());
    req(r.allFinite(), "Projection failed");
    return r;
}

struct Context {
    const System& s;
    double t;
    std::map<std::string, V> fld;
    std::map<E, V, std::owner_less<E>> dc, dq;

    Context(const System& ss, double tt, const V& u, bool neutral) : s(ss), t(tt) {
        V v = s.lift(t, u);
        for (int f = 0; f < s.nf; ++f)
            fld[s.compNames[f]] = v.segment(f * s.n, s.n);
        for (auto a : s.mdl->alg) {
            E ex = s.mdl->parse(a.at("src"));
            V src = co(ex);
            double charge = s.massv.dot(src);
            if (neutral)
                req(std::abs(charge) <= 1e-9 * std::max(s.w.dot((s.b * src).cwiseAbs()), s.w.sum()),
                    "Neumann Poisson source is not neutral");
            V rhs(s.n + 1);
            rhs.head(s.n) = s.mass * src / s.mdl->constant(a.value("a", J(1)));
            rhs[s.n] = 0;
            V pot = s.poissonLU->solve(rhs);
            req(pot.allFinite(), "Poisson solve failed");
            fld[a.at("nm").get<std::string>()] = pot.head(s.n);
        }
    }

    V qu(E e) {
        if (dq.count(e))
            return dq.at(e);
        V r;
        if (e->op == "num")
            r = V::Constant(s.w.size(), e->n);
        else if (!e->name.empty()) {
            if (fld.count(e->name))
                r = s.b * fld.at(e->name);
            else
                r = s.at(e, s.x, t);
        } else if (e->op == "add")
            r = qu(e->a[0]) + qu(e->a[1]);
        else if (e->op == "mul")
            r = qu(e->a[0]).array() * qu(e->a[1]).array();
        else if (e->op == "pow")
            r = qu(e->a[0]).array().pow(e->a[1]->n);
        else if (e->op == "lap" || e->op == "divg")
            r = s.b * co(e);
        else if (e->op == "intg")
            r = V::Constant(s.w.size(), s.w.dot(qu(e->a[0])));
        else {
            V a = qu(e->a[0]);
            if (e->op == "exp")
                r = a.array().exp();
            else if (e->op == "log")
                r = a.array().log();
            else if (e->op == "sin")
                r = a.array().sin();
            else if (e->op == "cos")
                r = a.array().cos();
            else if (e->op == "sqrt")
                r = a.array().sqrt();
            else
                throw Err("Unsupported quadrature expression: " + e->op);
        }
        req(r.allFinite(), "Expression left the finite real domain");
        dq[e] = r;
        return r;
    }

    V co(E e) {
        if (dc.count(e))
            return dc.at(e);
        V r;
        if (fld.count(e->name))
            r = fld.at(e->name);
        else if (constantExpr(e, *s.mdl))
            r = eval(e, s.mdl->par) * s.one;
        else if (e->op == "add")
            r = co(e->a[0]) + co(e->a[1]);
        else if (e->op == "lap")
            r = -s.massLU->solve(s.stiff * co(e->a[0]));
        else if (e->op == "divg") {
            V mob = qu(e->a[0]), phi = co(e->a[1]), rhs = V::Zero(s.n);
            for (auto& g : s.g)
                rhs += g.transpose() * (s.w.array() * mob.array() * (g * phi).array()).matrix();
            r = -s.massLU->solve(rhs);
        } else if (e->op == "mul" && constantExpr(e->a[0], *s.mdl))
            r = eval(e->a[0], s.mdl->par) * co(e->a[1]);
        else if (e->op == "mul" && constantExpr(e->a[1], *s.mdl))
            r = eval(e->a[1], s.mdl->par) * co(e->a[0]);
        else
            r = s.project(qu(e));
        req(r.allFinite(), "Non-finite projected operator");
        dc[e] = r;
        return r;
    }
};

V System::rhs(double t, const V& u, bool neutral) const {
    checkpoint();
    req(!mixed, "Mixed system requires constrained integration");
    req(u.size() == size && u.allFinite(), "Invalid state");
    if (exactSp && !spRhs.empty()) {
        auto d = mdl->par;
        for (int i = 0; i < size; ++i)
            d["u[" + std::to_string(i) + "]"] = u[i];
        V out(size);
        for (int i = 0; i < size; ++i)
            out[i] = eval(spRhs[i], d);
        return out;
    }
    Context ctx(*this, t, u, neutral);
    V r(size);
    for (int f = 0; f < nf; ++f) {
        V v = ctx.co(mdl->eq.at(compNames[f]));
        if (!bcs[f].fix.empty() || !bcs[f].natural.empty()) {
            V dt = boundaryDt(f, t),
              load = mass * v + boundaryLoad(f, t, ctx.fld.at(compNames[f])) - mass * dt;
            r.segment(f * n, n) = dt;
            if (!bcs[f].free.empty()) {
                V free = freeLU[f]->solve(pick(load, bcs[f].free));
                for (int i = 0; i < bcs[f].free.size(); ++i)
                    r[f * n + bcs[f].free[i]] = free[i];
            }
        } else
            r.segment(f * n, n) = v;
    }
    req(r.allFinite(), "Non-finite RHS");
    return r;
}

M System::jacobian(double t, const V& u) const {
    return jac(
        [&](const V& v) {
            return rhs(t, v);
        },
        u);
}

std::vector<V> System::algebraic(double t, const V& u) const {
    if (mixed || mdl->alg.empty())
        return {};
    Context ctx(*this, t, u, true);
    std::vector<V> a;
    for (auto d : mdl->alg)
        a.push_back(ctx.fld.at(d.at("nm").get<std::string>()));
    return a;
}

J poissonAudit(const System& s, double t, const V& u) {
    J out = J::array();
    if (s.mixed || s.mdl->alg.empty())
        return out;
    Context ctx(s, t, u, true);
    for (auto a : s.mdl->alg) {
        V phi = ctx.fld.at(a.at("nm").get<std::string>()), src = ctx.co(s.mdl->parse(a.at("src")));
        V residual = s.mdl->constant(a.value("a", J(1))) * (s.stiff * phi) - s.mass * src;
        out.push_back(J{{"field", a["nm"]},
                        {"weak_residual_l2_vector", residual.norm()},
                        {"source_integral", s.massv.dot(src)},
                        {"mean", s.massv.dot(phi)}});
    }
    return out;
}

double System::energy(double t, const V& u) const {
    if (mixed) {
        double en = 0;
        for (int f = 0; f < mdl->dim; ++f) {
            V v = u.segment(offset[f], spaces[f]->n);
            en += .5 * v.dot(spaces[f]->mass * v);
        }
        return en;
    }
    req(!mdl->en.is_null(), "No energy functional supplied");
    Context c(*this, t, u, true);
    double en = w.dot(c.qu(mdl->parse(mdl->en.at("f"))));
    for (int f = 0; f < nf; ++f) {
        double k = mdl->constant(mdl->en.value("kap", J::object()).value(compNames[f], J(0)));
        V v = u.segment(f * n, n);
        en += .5 * k * v.dot(stiff * v);
    }
    return en;
}

V System::norm(const V& u, bool dual) const {
    V out = V::Zero(mixed ? 2 : nf);
    for (int f = 0; f < compNames.size(); ++f) {
        int nn = mixed ? spaces[f]->n : n;
        V a = u.segment(offset[f], nn);
        double v;
        if (dual) {
            a = pick(a, bcs[f].free);
            v = a.size() ? a.dot(freeLU[f]->solve(a)) : 0;
        } else
            v = a.dot((mixed ? spaces[f]->mass : mass) * a);
        out[mixed ? (f < mdl->dim ? 0 : 1) : f] += std::max(0., v);
    }
    return out.array().sqrt();
}

std::pair<V, SM>
System::assemble(double t, const V& input, const V* old, double h, bool derivative) const {
    req(!terms.empty(), "Local sparse assembly unavailable");
    checkpoint();
    V u = lift(t, input);
    int nc = int(compNames.size()), ng = mesh->geo, nq = int(w.size());
    std::vector<std::vector<SM>> ops(nc);
    std::vector<std::vector<V>> vals(nc);
    for (int f = 0; f < nc; ++f) {
        auto s = spaces[f];
        ops[f].push_back(s->b);
        ops[f].insert(ops[f].end(), s->g.begin(), s->g.end());
        V a = u.segment(offset[f], s->n);
        for (auto& o : ops[f])
            vals[f].push_back(o * a);
        V vt = V::Zero(nq);
        if (old && h > 0)
            vt = s->b * ((a - old->segment(offset[f], s->n)) / h);
        vals[f].push_back(vt);
    }
    std::vector<V> rr(terms.size(), V(nq)), jj(derivative ? jterms.size() : 0, V(nq));

    struct Instr {
        int op = 0, a = 0, b = 0;
        double c = 0;
        const double* input = nullptr;
    };

    std::vector<Instr> program;
    std::map<std::string, const double*> inputs;
    for (int k = 0; k < ng; ++k)
        inputs[std::string(1, "xyz"[k])] = x.col(k).data();
    for (int f = 0; f < nc; ++f) {
        inputs[compNames[f]] = vals[f][0].data();
        for (int k = 0; k < ng; ++k)
            inputs[compNames[f] + "$" + std::to_string(k)] = vals[f][k + 1].data();
        inputs[compNames[f] + "$3"] = vals[f].back().data();
    }
    std::map<E, int, std::owner_less<E>> memo;
    std::function<int(E)> compile = [&](E e) {
        if (memo.count(e))
            return memo.at(e);
        Instr ins;
        if (e->op == "num")
            ins.c = e->n;
        else if (!e->name.empty()) {
            if (inputs.count(e->name)) {
                ins.op = 1;
                ins.input = inputs.at(e->name);
            } else if (e->name == "t")
                ins.c = t;
            else if (e->name == "pi")
                ins.c = pi;
            else if (e->name == "x" || e->name == "y" || e->name == "z")
                ins.c = 0;
            else {
                req(mdl->par.count(e->name), "Unresolved assembly symbol: " + e->name);
                ins.c = mdl->par.at(e->name);
            }
        } else {
            ins.a = compile(e->a[0]);
            if (e->op == "add" || e->op == "mul") {
                ins.b = compile(e->a[1]);
                ins.op = e->op == "add" ? 2 : 3;
            } else if (e->op == "pow") {
                ins.op = 4;
                ins.c = e->a[1]->n;
            } else if (e->op == "sin")
                ins.op = 5;
            else if (e->op == "cos")
                ins.op = 6;
            else if (e->op == "exp")
                ins.op = 7;
            else if (e->op == "log")
                ins.op = 8;
            else if (e->op == "sqrt")
                ins.op = 9;
            else
                throw Err("Unsupported assembly instruction: " + e->op);
        }
        int index = int(program.size());
        req(index < 200000, "Assembly expression program too large");
        program.push_back(ins);
        memo[e] = index;
        return index;
    };
    VI ri, ji;
    for (auto term : terms)
        ri.push_back(compile(term.e));
    if (derivative)
        for (auto term : jterms)
            ji.push_back(compile(term.e));
    std::vector<double> value(program.size());
    for (int i = 0; i < nq; ++i) {
        if (i % 1024 == 0)
            checkpoint();
        for (int k = 0; k < program.size(); ++k) {
            auto ins = program[k];
            double a = value[ins.a], b = value[ins.b], v = 0;
            switch (ins.op) {
                case 0:
                    v = ins.c;
                    break;
                case 1:
                    v = ins.input[i];
                    break;
                case 2:
                    v = a + b;
                    break;
                case 3:
                    v = a * b;
                    break;
                case 4:
                    v = std::pow(a, ins.c);
                    break;
                case 5:
                    v = std::sin(a);
                    break;
                case 6:
                    v = std::cos(a);
                    break;
                case 7:
                    v = std::exp(a);
                    break;
                case 8:
                    v = std::log(a);
                    break;
                case 9:
                    v = std::sqrt(a);
                    break;
            }
            if (!std::isfinite(v))
                throw Err("Assembly expression left the finite real domain");
            value[k] = v;
        }
        for (int j = 0; j < ri.size(); ++j)
            rr[j][i] = value[ri[j]];
        for (int j = 0; j < ji.size(); ++j)
            jj[j][i] = value[ji[j]];
    }
    V r = V::Zero(size);
    std::vector<Tr> jt;
    for (int j = 0; j < terms.size(); ++j) {
        auto term = terms[j];
        r.segment(offset[term.row], spaces[term.row]->n) +=
            ops[term.row][term.k].transpose() * (w.array() * rr[j].array()).matrix();
    }
    if (derivative)
        for (int j = 0; j < jterms.size(); ++j) {
            auto term = jterms[j];
            if (term.l == ng + 1 && !old)
                continue;
            double fac = term.l == ng + 1 ? 1 / h : 1;
            int l = term.l == ng + 1 ? 0 : term.l;
            SM block = ops[term.row][term.k].transpose() *
                       weighted(ops[term.col][l], (fac * w.array() * jj[j].array()).matrix());
            for (int k = 0; k < block.outerSize(); ++k)
                for (SM::InnerIterator it(block, k); it; ++it)
                    jt.emplace_back(
                        offset[term.row] + it.row(), offset[term.col] + it.col(), it.value());
        }
    if (!mixed)
        for (int f = 0; f < nf; ++f) {
            r.segment(offset[f], n) += boundaryLoad(f, t, u.segment(offset[f], n));
            if (derivative)
                for (auto& p : bcs[f].natural)
                    if (p.coef) {
                        V k = at(p.coef, p.x, t);
                        SM a = p.b.transpose() * weighted(p.b, (-p.w.array() * k.array()).matrix());
                        for (int col = 0; col < a.outerSize(); ++col)
                            for (SM::InnerIterator it(a, col); it; ++it)
                                jt.emplace_back(
                                    offset[f] + it.row(), offset[f] + it.col(), it.value());
                    }
        }
    SM jm(size, size);
    if (derivative)
        jm.setFromTriplets(jt.begin(), jt.end());
    return {r, jm};
}

V System::values(const V& u, const M& pts, int f) const {
    if (!mesh) {
        V out = V::Zero(pts.rows());
        for (int j = 0; j < n; ++j)
            for (int i = 0; i < pts.rows(); ++i) {
                double v = u[offset[f] + j];
                for (int k = 0; k < mdl->dim; ++k)
                    v *= std::cos(pi * modes[j][k] * pts(i, k) /
                                  mdl->par.at("Lx" + std::to_string(k + 1)));
                out[i] += v;
            }
        return out;
    }
    auto sp = spaces[f];
    int ne = int(mesh->c.rows());
    std::vector<M> inverse(ne), edge(ne);
    M origins(ne, mesh->geo), lo(ne, mesh->geo), hi(ne, mesh->geo);
    for (int e = 0; e < ne; ++e) {
        M points(mesh->dim + 1, mesh->geo);
        for (int k = 0; k <= mesh->dim; ++k)
            points.row(k) = mesh->p.row(mesh->c(e, k));
        origins.row(e) = points.row(0);
        lo.row(e) = points.colwise().minCoeff();
        hi.row(e) = points.colwise().maxCoeff();
        edge[e] = (points.bottomRows(mesh->dim).rowwise() - points.row(0)).transpose();
        inverse[e] = edge[e].completeOrthogonalDecomposition().pseudoInverse();
    }
    V out(pts.rows());
    for (int i = 0; i < pts.rows(); ++i) {
        if (i % 1024 == 0)
            checkpoint();
        bool found = false;
        for (int e = 0; e < ne; ++e) {
            double tol = 1e-9 * std::max(1., edge[e].norm());
            if ((pts.row(i).array() < lo.row(e).array() - tol).any() ||
                (pts.row(i).array() > hi.row(e).array() + tol).any())
                continue;
            V r = inverse[e] * (pts.row(i) - origins.row(e)).transpose();
            V l(mesh->dim + 1);
            l[0] = 1 - r.sum();
            l.tail(mesh->dim) = r;
            if (l.minCoeff() >= -1e-9 &&
                (edge[e] * r - (pts.row(i) - origins.row(e)).transpose()).norm() < tol) {
                M ll = l.transpose();
                auto bs = sp->basisAt(ll).first;
                double v = 0;
                for (int k = 0; k < sp->c.cols(); ++k)
                    v += bs(0, k) * u[offset[f] + sp->c(e, k)];
                out[i] = v;
                found = true;
                break;
            }
        }
        req(found, "Evaluation point is outside the mesh");
    }
    return out;
}
}
