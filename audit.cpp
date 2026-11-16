#include "phz.hpp"

namespace phz {
J poissonAudit(const System&, double, const V&);
double deadline(double);

void Result::diagnose() {
    auto& s = *sys;
    auto& m = *s.mdl;
    double tol = m.chk.value("tol", 1e-6);
    auto addCheck = [&](std::string name, double val, double lim, std::string scope) {
        checks.push_back(J{{"name", name},
                           {"value", val},
                           {"limit", lim},
                           {"status", val <= lim ? "pass" : "fail"},
                           {"scope", scope}});
    };
    J mass = J::array(), minimum = J::array(), maximum = J::array(), l2 = J::array(),
      en = J::array(), alg = J::array(), divergence = J::array(), cont = J::array(),
      mean = J::array(), bcerr = J::array(), weakr = J::array();
    J exact = m.chk.value("exact", J::object());
    int nc = int(s.compNames.size());
    std::vector<double> maxerr(s.mixed ? 2 : nc, 0), minval(nc, INFINITY);
    double massErr = 0, divMax = 0, weakMax = 0;
    for (int it = 0; it < u.size(); ++it) {
        checkpoint();
        J mm = J::array(), mi = J::array(), ma = J::array(), ee = J::array();
        std::vector<double> err(nc, 0);
        double dir = 0;
        for (int f = 0; f < nc; ++f) {
            int n = s.mixed ? s.spaces[f]->n : s.n;
            V a = u[it].segment(s.offset[f], n);
            const SM& b = s.mixed ? s.spaces[f]->b : s.b;
            const SM& massM = s.mixed ? s.spaces[f]->mass : s.mass;
            V q = b * a;
            V one = s.mixed ? V(V::Ones(n)) : s.one;
            mm.push_back(one.dot(massM * a));
            double mn = q.minCoeff(), mx = q.maxCoeff();
            if (s.mesh) {
                mn = std::min(mn, a.minCoeff());
                mx = std::max(mx, a.maxCoeff());
            }
            mi.push_back(mn);
            ma.push_back(mx);
            minval[f] = std::min(minval[f], mn);
            std::string field = s.mixed ? m.fields[f < m.dim ? 0 : 1] : s.compNames[f];
            if (exact.contains(field)) {
                E e = m.parse(exact[field]);
                if (s.mixed && f < m.dim) {
                    if (zero(e))
                        e = num(0);
                    else {
                        req(e->sh == VI{m.dim}, "Exact velocity shape mismatch");
                        e = components(e, m.dim)[f];
                    }
                }
                V diff = q - s.at(e, s.x, t[it]);
                err[f] = std::sqrt(std::max(0., s.w.dot(diff.array().square().matrix())));
            }
            auto& boundary = s.bcs[f];
            for (int k = 0; k < boundary.fix.size(); ++k) {
                M pt = s.spaces[f]->p.row(boundary.fix[k]);
                dir = std::max(dir,
                               std::abs(a[boundary.fix[k]] - s.at(boundary.val[k], pt, t[it])[0]));
            }
        }
        if (s.mixed) {
            double ev = 0;
            for (int f = 0; f < m.dim; ++f)
                ev += err[f] * err[f];
            ee.push_back(std::sqrt(ev));
            ee.push_back(err.back());
            maxerr[0] = std::max(maxerr[0], std::sqrt(ev));
            maxerr[1] = std::max(maxerr[1], err.back());
            int nv = s.offset[m.dim];
            V dv = V::Zero(s.w.size());
            for (int f = 0; f < m.dim; ++f)
                dv += s.spaces[f]->g[f] * u[it].segment(s.offset[f], s.spaces[f]->n);
            double dn = std::sqrt(s.w.dot(dv.array().square().matrix())),
                   cn = (s.div * u[it].head(nv)).norm();
            divergence.push_back(dn);
            cont.push_back(cn);
            mean.push_back(s.gauge.dot(u[it]));
            divMax = std::max(divMax, dn);
            double residual = 0;
            if (stationary)
                residual = pick(s.assemble(t[it], u[it], nullptr, 0, false).first, s.active).norm();
            weakr.push_back(stationary ? J(residual) : J());
            weakMax = std::max(weakMax, residual);
        } else {
            for (int f = 0; f < nc; ++f) {
                ee.push_back(err[f]);
                maxerr[f] = std::max(maxerr[f], err[f]);
            }
            J ag = J::array();
            for (auto v : s.algebraic(t[it], u[it]))
                ag.push_back(arr(v));
            alg.push_back(ag);
            if (!s.mdl->alg.empty())
                diag["poisson_residuals"].push_back(poissonAudit(s, t[it], u[it]));
        }
        mass.push_back(mm);
        minimum.push_back(mi);
        maximum.push_back(ma);
        l2.push_back(ee);
        bcerr.push_back(dir);
        if (s.mixed || !m.en.is_null())
            en.push_back(s.energy(t[it], u[it]));
    }
    diag["mass"] = mass;
    diag["min"] = minimum;
    diag["max"] = maximum;
    diag["l2"] = l2;
    diag["dirichlet_residual_max"] = bcerr;
    if (!alg.empty())
        diag["algebraic_fields"] = alg;
    if (!en.empty())
        diag["energy"] = en;
    if (s.mixed) {
        diag["div_l2"] = divergence;
        diag["continuity_weak"] = cont;
        diag["pressure_mean"] = mean;
        diag["stationary_weak_residual"] = weakr;
        double cn = 0, pn = 0, bn = 0;
        for (auto x : cont)
            cn = std::max(cn, x.get<double>());
        for (auto x : mean)
            pn = std::max(pn, std::abs(x.get<double>()));
        for (auto x : bcerr)
            bn = std::max(bn, x.get<double>());
        addCheck("continuity", cn, 1e-8, "Weak discrete divergence moments");
        addCheck("pressure_gauge", pn, 1e-8, "Mean-zero pressure");
        addCheck("boundary", bn, 1e-10, "Prescribed velocity at boundary DOFs");
        if (m.chk.contains("div"))
            addCheck("divergence_l2",
                     divMax,
                     m.chk["div"],
                     "Pointwise divergence integrated by quadrature");
    }
    for (auto field : m.chk.value("mass", J::array())) {
        std::string f = field;
        auto pos = std::find(s.compNames.begin(), s.compNames.end(), f);
        req(pos != s.compNames.end(), "Mass check for unknown field");
        int k = int(pos - s.compNames.begin());
        if (stationary)
            checks.push_back(J{{"name", "mass:" + f},
                               {"status", "not_run"},
                               {"scope", "One equilibrium cannot test time conservation"}});
        else {
            double err = 0, initial = mass[0][k];
            for (auto a : mass)
                err = std::max(err, std::abs(a[k].get<double>() - initial));
            addCheck("mass:" + f,
                     err,
                     tol * std::max(1., std::abs(initial)),
                     "Conservation between saved times");
        }
    }
    for (auto field : m.chk.value("nonneg", J::array())) {
        std::string f = field;
        auto pos = std::find(s.compNames.begin(), s.compNames.end(), f);
        req(pos != s.compNames.end(), "Nonnegative check for unknown field");
        addCheck("nonnegative:" + f,
                 std::max(0., -minval[pos - s.compNames.begin()]),
                 tol,
                 "Quadrature points and mesh DOFs");
    }
    for (int f = 0; f < m.fields.size(); ++f)
        if (exact.contains(m.fields[f]))
            addCheck("exact:" + m.fields[f],
                     maxerr[f],
                     tol,
                     "L2 error against supplied analytical expression");
    if (m.chk.value("energy", false)) {
        if (en.size() < 2)
            checks.push_back(J{{"name", "energy"},
                               {"status", "not_run"},
                               {"scope", "One equilibrium cannot test energy decay"}});
        else {
            double growth = 0;
            for (int i = 1; i < en.size(); ++i)
                growth = std::max(growth, en[i].get<double>() - en[i - 1].get<double>());
            addCheck("energy",
                     growth,
                     tol * std::max(1., std::abs(en[0].get<double>())),
                     "Non-increase at saved times");
        }
    }

    if (stationary && !s.mixed) {
        J masses = m.stat.value("mass", J::object());
        for (auto it = masses.begin(); it != masses.end(); ++it) {
            int f = int(std::find(s.compNames.begin(), s.compNames.end(), it.key()) -
                        s.compNames.begin());
            double target = it.value().is_string() ? s.massv.dot(s.u0.segment(f * s.n, s.n))
                                                   : number(it.value(), "mass target");
            addCheck("constraint_mass:" + it.key(),
                     std::abs(mass[0][f].get<double>() - target),
                     std::max(1., std::abs(target)) * m.stat.value("tol", 1e-9),
                     "Prescribed stationary integral");
        }
    }
    if (scheme.contains("heat") && scheme["heat"].value("dt_status", "") == "fail")
        checks.push_back(J{{"name", "heat_step_stability"},
                           {"status", "fail"},
                           {"scope", "Step exceeds derived negative-axis stability interval"}});
    if (checks.empty())
        checks.push_back(J{{"name", "validation"},
                           {"status", "not_run"},
                           {"scope", "No physical or reference checks requested"}});
    info["assurance"] =
        "Numerical verification within stated model assumptions; no formal proof or experimental validation";
}

bool Result::passed() const {
    for (auto c : checks)
        if (c.value("status", "") == "fail")
            return false;
    return true;
}

J quadratureConfig(const Model& m) {
    if (!m.disc.contains("aq"))
        return J();
    J a = m.disc["aq"];
    keys(a, "atol rtol qmax step confirm timeout", "disc.aq");
    bool sp = m.type == "sp";
    VI counts = m.disc.value("n", VI(m.dim, 4));
    int q = m.disc.value(
            "q", sp ? std::max(16, 3 * *std::max_element(counts.begin(), counts.end())) : 4),
        lo = sp               ? std::max(4, 2 * *std::max_element(counts.begin(), counts.end()))
             : m.type == "th" ? 4
             : m.dim == 3     ? 3
                              : 2,
        hi = sp ? 100 : 20;
    J out{{"atol", 1e-8},
          {"rtol", 1e-5},
          {"qmax", sp ? 100 : 10},
          {"step", sp ? 4 : 2},
          {"confirm", 2},
          {"timeout", 300.},
          {"qstart", q}};
    for (auto it = a.begin(); it != a.end(); ++it)
        out[it.key()] = it.value();
    for (auto s : {"atol", "rtol", "timeout"})
        req(number(out[s], s) >=
                (std::string(s) == "rtol" ? 0 : std::numeric_limits<double>::min()),
            "Invalid quadrature tolerance");
    integer(q, lo, hi, "qstart");
    int maxq = integer(out["qmax"], lo, hi, "qmax"), step = integer(out["step"], 1, hi, "qstep"),
        confirm = integer(out["confirm"], 2, 4, "qconfirm");
    req(int(std::ceil(double(maxq - q) / step)) >= confirm,
        "qmax cannot satisfy consecutive confirmations");
    return out;
}

Result execute(J j, bool st) {
    auto m = std::make_shared<Model>(j);
    J cfg = quadratureConfig(*m);
    if (cfg.is_null())
        return st ? stationary(std::make_shared<System>(j))
                  : integrate(std::make_shared<System>(j));
    J original = j, history = J::array();
    j["disc"].erase("aq");
    int q = cfg["qstart"], qmax = cfg["qmax"], step = cfg["step"], required = cfg["confirm"],
        good = 0;
    double atol = cfg["atol"], rtol = cfg["rtol"], limit = cfg["timeout"];
    auto started = std::chrono::steady_clock::now();

    struct Deadline {
        double old;

        ~Deadline() {
            deadline(old);
        }
    } guard{deadline(std::chrono::duration<double>(started.time_since_epoch()).count() + limit)};

    std::unique_ptr<Result> previous;
    for (;;) {
        checkpoint();
        req(std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count() <
                limit,
            "Quadrature control timeout");
        j["disc"]["q"] = q;
        auto s = std::make_shared<System>(j);
        if (st && !previous && !s->mixed && j.value("stat", J::object()).contains("mass")) {
            J masses = j["stat"]["mass"];
            for (auto it = masses.begin(); it != masses.end(); ++it)
                if (it.value().is_string() && it.value() == "initial") {
                    int f = int(std::find(s->compNames.begin(), s->compNames.end(), it.key()) -
                                s->compNames.begin());
                    req(f < s->nf, "Unknown mass constraint field");
                    j["stat"]["mass"][it.key()] = s->massv.dot(s->u0.segment(f * s->n, s->n));
                }
        }
        Result now = st ? stationary(s) : integrate(s);
        J row{{"q", q}, {"max_ratio", 0}, {"consecutive_passes", 0}, {"metrics", J::object()}};
        double worst = 0;
        if (previous) {
            auto& old = *previous;
            req(old.sys->size == s->size && old.sys->active == s->active && old.t == now.t,
                "Quadrature adaptation changed approximation space or times");
            auto metric = [&](std::string name, V diff, V scale) {
                V ratio = diff.array() / (atol + rtol * scale.array());
                req(ratio.allFinite(), "Non-finite quadrature metric");
                double r = ratio.maxCoeff();
                worst = std::max(worst, r);
                if (!row["metrics"].contains(name) ||
                    r > row["metrics"][name]["max_ratio"].get<double>())
                    row["metrics"][name] = J{{"absolute", arr(diff)},
                                             {"scale", arr(scale)},
                                             {"ratio", arr(ratio)},
                                             {"max_ratio", r}};
            };
            auto opval = [](const System& s, double t, const V& u) {
                return s.mixed ? s.assemble(t, u, nullptr, 0, false).first : s.rhs(t, u);
            };
            for (auto init : {old.sys->u0, s->u0}) {
                V a = opval(*old.sys, now.t[0], init), b = opval(*s, now.t[0], init);
                metric("initial_operator",
                       s->norm(b - a, s->mixed),
                       s->norm(a, s->mixed).cwiseMax(s->norm(b, s->mixed)));
            }
            for (int i = 0; i < now.u.size(); ++i) {
                metric("solution_l2",
                       s->norm(now.u[i] - old.u[i]),
                       s->norm(now.u[i]).cwiseMax(s->norm(old.u[i])));
                for (auto state : {old.u[i], now.u[i]}) {
                    V a = opval(*old.sys, now.t[i], state), b = opval(*s, now.t[i], state);
                    metric("operator",
                           s->norm(b - a, s->mixed),
                           s->norm(a, s->mixed).cwiseMax(s->norm(b, s->mixed)));
                    if (!s->mixed && !s->mdl->alg.empty()) {
                        auto aa = old.sys->algebraic(now.t[i], state),
                             ab = s->algebraic(now.t[i], state);
                        V diff(aa.size()), scale(aa.size());
                        for (int k = 0; k < aa.size(); ++k) {
                            V d = ab[k] - aa[k];
                            diff[k] = std::sqrt(std::max(0., d.dot(s->mass * d)));
                            scale[k] = std::sqrt(
                                std::max(aa[k].dot(s->mass * aa[k]), ab[k].dot(s->mass * ab[k])));
                        }
                        metric("algebraic_l2", diff, scale);
                    }
                    if (!s->mdl->en.is_null() || s->mixed) {
                        double a = old.sys->energy(now.t[i], state), b = s->energy(now.t[i], state);
                        V dif(1), sc(1);
                        dif << std::abs(a - b);
                        sc << std::max(std::abs(a), std::abs(b));
                        metric("energy", dif, sc);
                    }
                }
            }
            good = worst <= 1 ? good + 1 : 0;
            row["max_ratio"] = worst;
            row["consecutive_passes"] = good;
            row["saved_times"] = now.t.size();
        }
        history.push_back(row);
        bool pass = good >= required, done = pass || q == qmax;
        if (done) {
            now.info["quadrature_control"] = J{
                {"status", pass ? "pass" : "fail"},
                {"selected_q", q},
                {"reason", pass ? "consecutive_confirmations" : "qmax_reached"},
                {"history", history},
                {"config", cfg},
                {"scope",
                 "Same-space quadrature consistency; not a mesh convergence or physical validation certificate"}};
            now.info["quadrature_request"] = original;
            now.checks.push_back(J{{"name", "quadrature_control"},
                                   {"status", pass ? "pass" : "fail"},
                                   {"value", worst},
                                   {"limit", 1},
                                   {"achieved_consecutive", good}});
            return now;
        }
        previous = std::make_unique<Result>(std::move(now));
        q = std::min(qmax, q + step);
    }
}
}
