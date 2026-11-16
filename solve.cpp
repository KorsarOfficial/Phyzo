#include "phz.hpp"
#include "coeff.hpp"
#include <cvode/cvode.h>
#include <cvode/cvode_ls.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>

namespace phz {
struct Watch {
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    double limit;

    explicit Watch(double x) : limit(x) {
    }

    void check() const {
        checkpoint();
        if (std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count() >=
            limit)
            throw Stop("Calculation timeout");
    }

    double seconds() const {
        return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    }
};

static V full(const System& s, double t, const V& a) {
    V u = s.u0;
    put(u, s.active, a);
    return s.lift(t, u);
}

static SM saddle(const SM& j, const M& c) {
    int n = int(j.rows()), k = int(c.rows());
    std::vector<Tr> a;
    for (int col = 0; col < j.outerSize(); ++col)
        for (SM::InnerIterator it(j, col); it; ++it)
            a.emplace_back(it.row(), it.col(), it.value());
    for (int r = 0; r < k; ++r)
        for (int col = 0; col < n; ++col)
            if (c(r, col)) {
                a.emplace_back(n + r, col, c(r, col));
                a.emplace_back(col, n + r, c(r, col));
            }
    SM out(n + k, n + k);
    out.setFromTriplets(a.begin(), a.end());
    return out;
}

static V newton(System& s, double t, V u, const V* old, double h, J& trace, J& lin, Watch& watch) {
    int n = int(s.active.size());
    bool useWeak = s.mixed || (!s.terms.empty() && s.mdl->stat.value("solver", "auto") != "dense");
    req(useWeak || s.mdl->stat.value("solver", "auto") != "sparse",
        "Requested sparse assembly is unavailable for this operator");
    req(useWeak || n <= 256, "Dense stationary solver is limited to 256 free coefficients");
    double tol = s.mdl->stat.value("tol", s.mixed ? 1e-10 : 1e-9);
    int maxit = s.mdl->stat.value("max_it", 40), maxls = s.mdl->stat.value("max_ls", 24);
    M constraints(0, n);
    V target(0);
    if (s.mixed) {
        constraints.resize(1, n);
        constraints.row(0) = pick(s.gauge, s.active);
        target = V::Zero(1);
    } else {
        J masses = s.mdl->stat.value("mass", J::object());
        req(masses.is_object(), "stat.mass must be an object");
        constraints = M::Zero(masses.size(), n);
        target.resize(masses.size());
        int row = 0;
        for (auto it = masses.begin(); it != masses.end(); ++it) {
            auto pos = std::find(s.compNames.begin(), s.compNames.end(), it.key());
            req(pos != s.compNames.end(), "Mass target for unknown field");
            int f = int(pos - s.compNames.begin());
            V mf = V::Zero(s.size);
            mf.segment(f * s.n, s.n) = s.massv;
            double tg = it.value().is_string() && it.value().get<std::string>() == "initial"
                            ? mf.dot(s.u0)
                            : number(it.value(), "mass target");
            V cf = pick(mf, s.active);
            double norm = cf.norm();
            req(norm > 0, "Mass constraint has no free coefficients");
            V zero = full(s, t, V::Zero(n));
            constraints.row(row) = cf / norm;
            target[row] = (tg - mf.dot(zero)) / norm;
            ++row;
        }
    }
    V v = pick(s.lift(t, u), s.active);
    if (constraints.rows())
        v += constraints.transpose() *
             (constraints * constraints.transpose()).ldlt().solve(target - constraints * v);
    double alpha = 0;
    auto residual = [&](const V& x, bool derivative) {
        watch.check();
        V uu = full(s, t, x);
        if (useWeak)
            return s.assemble(t, uu, old, h, derivative);
        V r = s.rhs(t, uu);
        SM j;
        if (derivative)
            j = sparse(s.jacobian(t, uu));
        return std::make_pair(r, j);
    };
    for (int it = 0; it <= maxit; ++it) {
        auto [all, j] = residual(v, true);
        V r = pick(all, s.active);
        V rr = r;
        if (!s.mixed && useWeak) {
            int at = 0;
            for (int f = 0; f < s.nf; ++f) {
                int sz = int(s.bcs[f].free.size());
                if (sz) {
                    rr.segment(at, sz) = s.freeLU[f]->solve(r.segment(at, sz));
                    at += sz;
                }
            }
        }
        double fn = rr.size() ? rr.lpNorm<Eigen::Infinity>() : 0;
        trace.push_back(J::array({t, h, it, fn, r.norm(), alpha}));
        if (fn <= tol) {
            u = full(s, t, v);
            if (constraints.rows())
                req((constraints * v - target).lpNorm<Eigen::Infinity>() <= std::max(1e-10, tol),
                    "Mass/gauge constraint lost accuracy");
            if (!s.mixed && !s.active.empty())
                req(pick(s.rhs(t, u), s.active).lpNorm<Eigen::Infinity>() <= tol,
                    "Full PDE residual exceeds tolerance");
            return u;
        }
        req(it < maxit, "Newton iteration limit; residual=" + txt(fn));
        SM jf = block(j, s.active, s.active);
        V rhs(n + constraints.rows());
        rhs.head(n) = -r;
        if (constraints.rows())
            rhs.tail(constraints.rows()) = target - constraints * v;
        SM a = saddle(jf, constraints);
        V step = solve(a, rhs);
        lin.push_back(J::array({it, (a * step - rhs).norm()}));
        double before = r.norm();
        bool accepted = false;
        for (int ls = 0; ls < maxls; ++ls) {
            alpha = std::ldexp(1., -ls);
            V trial = v + alpha * step.head(n);
            try {
                V rnew = pick(residual(trial, false).first, s.active);
                if (rnew.norm() < (1 - 1e-4 * alpha) * before) {
                    v = trial;
                    accepted = true;
                    break;
                }
            } catch (const Stop&) {
                throw;
            } catch (const Err&) {
            }
        }
        req(accepted, "Newton line search failed; residual=" + txt(fn));
    }
    throw Err("Newton failed");
}

Result stationary(std::shared_ptr<System> s) {
    Watch watch(s->mdl->stat.value("timeout", 120.));
    for (auto [name, e] : s->mdl->eq)
        req(!depends(e, "t"), "Stationary PDE must be autonomous");
    for (auto& bc : s->bcs) {
        for (auto e : bc.val)
            req(!depends(e, "t"), "Stationary boundary must be autonomous");
        for (auto p : bc.natural)
            req(!depends(p.val, "t") && (!p.coef || !depends(p.coef, "t")),
                "Stationary flux must be autonomous");
    }
    if (s->mixed)
        for (auto t : s->terms)
            req(!depends(t.e, "t"), "Stationary forcing must be autonomous");
    Result out;
    out.sys = s;
    out.stationary = true;
    out.t = {0};
    J trace = J::array(), lin = J::array();
    V u = newton(*s, 0, s->u0, nullptr, 0, trace, lin, watch);
    out.u = {u};
    out.diag["stationary_trace"] = trace;
    out.diag["stationary_linear_trace"] = lin;
    J st{{"iterations", trace.size() - 1},
         {"residual", trace.back()[3]},
         {"tol", s->mdl->stat.value("tol", s->mixed ? 1e-10 : 1e-9)},
         {"solver",
          s->terms.empty() || s->mdl->stat.value("solver", "auto") == "dense"
              ? "domain-aware dense Newton"
              : "analytic sparse weak Newton"},
         {"mass_targets", s->mdl->stat.value("mass", J::object())},
         {"stability", "not_run"}};
    if (!s->mixed && s->mdl->stat.value("eig", true)) {
        int n = int(s->active.size());
        if (n == 0)
            st["stability"] = "no_free_coefficients";
        else if (n > 256)
            st["stability"] = "not_computed_size_limit";
        else {
            M jacobian = s->jacobian(0, u);
            M a(n, n);
            for (int i = 0; i < n; ++i)
                for (int j = 0; j < n; ++j)
                    a(i, j) = jacobian(s->active[i], s->active[j]);
            J masses = s->mdl->stat.value("mass", J::object());
            M c = M::Zero(masses.size(), n);
            int k = 0;
            for (auto it = masses.begin(); it != masses.end(); ++it) {
                int f = int(std::find(s->compNames.begin(), s->compNames.end(), it.key()) -
                            s->compNames.begin());
                for (int i = 0; i < n; ++i)
                    if (s->active[i] / s->n == f)
                        c(k, i) = s->massv[s->active[i] % s->n];
                ++k;
            }
            M z = M::Identity(n, n);
            if (c.rows()) {
                Eigen::JacobiSVD<M> svd(c, Eigen::ComputeFullV);
                z = svd.matrixV().rightCols(n - c.rows());
            }
            double leak = c.rows() ? (c * a * z).norm() : 0;
            M ar = z.transpose() * a * z;
            Eigen::EigenSolver<M> es(ar, false);
            req(es.info() == Eigen::Success, "Stationary eigenvalue analysis failed");
            J vals = J::array();
            double top = -INFINITY;
            for (auto e : es.eigenvalues()) {
                vals.push_back(J::array({e.real(), e.imag()}));
                top = std::max(top, e.real());
            }
            double tol = s->mdl->stat.value("eig_tol", 1e-7);
            st["eigenvalues"] = vals;
            st["tangent_leakage"] = leak;
            st["stability"] = leak > 1e-7  ? "inconclusive_noninvariant_constraints"
                              : top < -tol ? "asymptotically_stable"
                              : top > tol  ? "unstable"
                                           : "inconclusive";
        }
    }
    out.info = J{{"kind", "stationary"},
                 {"stationary", st},
                 {"seconds", watch.seconds()},
                 {"mixed_flow", s->mixed}};
    out.diagnose();
    return out;
}

struct CvCtx {
    System* s;
    Watch* watch;
    std::string error;
};

static int cvRhs(double t, N_Vector y, N_Vector dy, void* data) {
    auto& c = *static_cast<CvCtx*>(data);
    try {
        c.watch->check();
        V v = Eigen::Map<V>(NV_DATA_S(y), NV_LENGTH_S(y));
        V r = pick(c.s->rhs(t, full(*c.s, t, v)), c.s->active);
        Eigen::Map<V>(NV_DATA_S(dy), r.size()) = r;
        return 0;
    } catch (const Stop& e) {
        c.error = e.what();
        return -1;
    } catch (const std::exception& e) {
        c.error = e.what();
        return 1;
    }
}

static int
cvJac(double t, N_Vector y, N_Vector, SUNMatrix a, void* data, N_Vector, N_Vector, N_Vector) {
    auto& c = *static_cast<CvCtx*>(data);
    try {
        c.watch->check();
        V u = Eigen::Map<V>(NV_DATA_S(y), NV_LENGTH_S(y));
        M j = jac(
            [&](const V& v) {
                return pick(c.s->rhs(t, full(*c.s, t, v)), c.s->active);
            },
            u);
        for (int i = 0; i < j.rows(); ++i)
            for (int k = 0; k < j.cols(); ++k)
                SM_ELEMENT_D(a, i, k) = j(i, k);
        return 0;
    } catch (const std::exception& e) {
        c.error = e.what();
        return -1;
    }
}

struct CvMem {
    SUNContext ctx = nullptr;
    N_Vector y = nullptr;
    SUNMatrix a = nullptr;
    SUNLinearSolver ls = nullptr;
    void* mem = nullptr;

    ~CvMem() {
        if (mem)
            CVodeFree(&mem);
        if (ls)
            SUNLinSolFree(ls);
        if (a)
            SUNMatDestroy(a);
        if (y)
            N_VDestroy(y);
        if (ctx)
            SUNContext_Free(&ctx);
    }
};

static V radauStep(
    const std::function<V(double, const V&)>& fn, double t, const V& u, double h, double tol) {
    int n = int(u.size());
    double r = std::sqrt(6.);
    M a(3, 3);
    a << (88 - 7 * r) / 360, (296 - 169 * r) / 1800, (-2 + 3 * r) / 225, (296 + 169 * r) / 1800,
        (88 + 7 * r) / 360, (-2 - 3 * r) / 225, (16 - r) / 36, (16 + r) / 36, 1. / 9;
    V c(3);
    c << (4 - r) / 10, (4 + r) / 10, 1;
    V stages(3 * n);
    for (int i = 0; i < 3; ++i)
        stages.segment(i * n, n) = u;
    auto residual = [&](const V& v) {
        V out(3 * n);
        std::vector<V> f;
        for (int k = 0; k < 3; ++k)
            f.push_back(fn(t + c[k] * h, v.segment(k * n, n)));
        for (int i = 0; i < 3; ++i) {
            out.segment(i * n, n) = v.segment(i * n, n) - u;
            for (int k = 0; k < 3; ++k)
                out.segment(i * n, n) -= h * a(i, k) * f[k];
        }
        return out;
    };
    for (int it = 0; it < 12; ++it) {
        V f = residual(stages);
        if (f.lpNorm<Eigen::Infinity>() <= tol)
            return stages.tail(n);
        M j = M::Identity(3 * n, 3 * n);
        for (int k = 0; k < 3; ++k) {
            M jf = jac(
                [&](const V& x) {
                    return fn(t + c[k] * h, x);
                },
                stages.segment(k * n, n));
            for (int i = 0; i < 3; ++i)
                j.block(i * n, k * n, n, n) -= h * a(i, k) * jf;
        }
        V step = j.partialPivLu().solve(-f);
        bool ok = false;
        for (int ls = 0; ls < 12; ++ls) {
            V v = stages + std::ldexp(1., -ls) * step;
            try {
                if (residual(v).norm() < f.norm()) {
                    stages = v;
                    ok = true;
                    break;
                }
            } catch (const Err&) {
            }
        }
        req(ok, "Radau stage line search failed");
    }
    throw Err("Radau stage iteration limit");
}

Result integrate(std::shared_ptr<System> s) {
    auto run = s->mdl->run;
    double t0 = run.value("t", J::array({0, 1}))[0], t1 = run.value("t", J::array({0, 1}))[1],
           rtol = run.value("rtol", 1e-7), atol = run.value("atol", 1e-9),
           maxh = run.value("max_step", INFINITY);
    int nout = run.value("n", 51);
    req(int64_t(nout) * s->size <= 8000000, "Saved states exceed memory limit");
    std::string method = run.value("method", s->mixed ? "BE" : "BDF");
    Watch watch(run.value("timeout", 120.));
    Result out;
    out.sys = s;
    out.info = J{{"kind", "transient"},
                 {"method", method},
                 {"rtol", rtol},
                 {"atol", atol},
                 {"mixed_flow", s->mixed}};
    V current = s->mixed ? s->projectVelocity(s->u0) : s->u0;
    for (int i = 0; i < nout; ++i)
        out.t.push_back(t0 + (t1 - t0) * i / (nout - 1));
    out.u.push_back(current);
    J trace = J::array(), lin = J::array(), stepT = J::array({t0}),
      stepU = J::array({arr(current)});
    int64_t steps = 0, fev = 0;
    auto fn = [&](double t, const V& v) {
        watch.check();
        ++fev;
        return pick(s->rhs(t, full(*s, t, v)), s->active);
    };
    if (s->mixed) {
        req(s->transient, "Stationary mixed PDE has no dt term");
        double t = t0;
        for (int io = 1; io < nout; ++io) {
            int count = std::isfinite(maxh)
                            ? std::max(1, int(std::ceil((out.t[io] - t) / maxh - 1e-12)))
                            : 1;
            req(count <= 100000, "Too many mixed time steps");
            double h = (out.t[io] - t) / count;
            for (int k = 0; k < count; ++k) {
                watch.check();
                V prev = current;
                current = newton(*s, t + h, current, &prev, h, trace, lin, watch);
                t += h;
                ++steps;
                stepT.push_back(t);
                req(stepU.size() * size_t(s->size) <= 8000000,
                    "Accepted-state history exceeds memory limit");
                stepU.push_back(arr(current));
            }
            out.u.push_back(current);
        }
    } else if (s->active.empty()) {
        for (int i = 1; i < nout; ++i)
            out.u.push_back(s->lift(out.t[i], current));
    } else if (method == "BDF") {
        CvMem c;
        req(SUNContext_Create(SUN_COMM_NULL, &c.ctx) == 0, "SUNDIALS context failed");
        int n = int(s->active.size());
        req(n <= 2048,
            "Dense CVODE backend is limited to 2048 free coefficients; use a sparse stationary or mixed solver for larger systems");
        c.y = N_VNew_Serial(n, c.ctx);
        req(c.y != nullptr, "SUNDIALS allocation failed");
        Eigen::Map<V>(NV_DATA_S(c.y), n) = pick(current, s->active);
        c.mem = CVodeCreate(CV_BDF, c.ctx);
        req(c.mem != nullptr, "CVODE allocation failed");
        CvCtx ctx{s.get(), &watch, {}};
        req(CVodeInit(c.mem, cvRhs, t0, c.y) == 0 && CVodeSStolerances(c.mem, rtol, atol) == 0,
            "CVODE initialization failed");
        CVodeSetUserData(c.mem, &ctx);
        CVodeSetMaxNumSteps(c.mem, 100000);
        if (std::isfinite(maxh))
            CVodeSetMaxStep(c.mem, maxh);
        c.a = SUNDenseMatrix(n, n, c.ctx);
        c.ls = SUNLinSol_Dense(c.y, c.a, c.ctx);
        req(c.a && c.ls, "CVODE linear solver allocation failed");
        req(CVodeSetLinearSolver(c.mem, c.ls, c.a) == 0 && CVodeSetJacFn(c.mem, cvJac) == 0,
            "CVODE linear solver setup failed");
        for (int i = 1; i < nout; ++i) {
            double t;
            int code = CVode(c.mem, out.t[i], c.y, &t, CV_NORMAL);
            req(code >= 0, "BDF failed: " + ctx.error + " (" + std::to_string(code) + ")");
            out.u.push_back(full(*s, out.t[i], Eigen::Map<V>(NV_DATA_S(c.y), n)));
        }
        long int nni = 0, ncf = 0;
        CVodeGetNonlinSolvStats(c.mem, &nni, &ncf);
        out.info["nonlinear_iterations"] = nni;
        out.info["nonlinear_convergence_failures"] = ncf;
        out.info["internal_newton_residual"] =
            "CVODE convergence statistics; individual Newton residuals are not recorded";
        long int ns = 0, nf = 0, nj = 0;
        CVodeGetNumSteps(c.mem, &ns);
        CVodeGetNumRhsEvals(c.mem, &nf);
        CVodeGetNumJacEvals(c.mem, &nj);
        steps = ns;
        fev = nf;
        out.info["njev"] = nj;
        out.info["backend"] = "SUNDIALS CVODE 7.4.0 variable-order BDF";
    } else {
        V u = pick(current, s->active);
        double t = t0, h = std::min(std::isfinite(maxh) ? maxh : (t1 - t0) / 10, (t1 - t0) / 10);
        bool fixed = method == "FE" || method == "BE" || method == "RK2" || method == "RK4" ||
                     method == "SSPRK3";
        for (int io = 1; io < nout; ++io) {
            int fixedN = fixed ? (std::isfinite(maxh)
                                      ? std::max(1, int(std::ceil((out.t[io] - t) / maxh - 1e-12)))
                                      : 1)
                               : 0;
            req(fixedN <= 100000, "Too many time steps");
            double fixedH = fixed ? (out.t[io] - t) / fixedN : 0;
            while (t < out.t[io] - 8 * std::numeric_limits<double>::epsilon() *
                                       std::max(1., std::abs(out.t[io]))) {
                watch.check();
                req(++steps <= 100000, "Time step limit");
                h = std::min(fixed ? fixedH : h, out.t[io] - t);
                h = std::min(h, maxh);
                req(t + h > t, "Time step below floating-point resolution");
                V v = u, err = V::Zero(u.size());
                double errNorm = 0;
                bool ok = true;
                try {
                    if (method == "BE") {
                        for (int it = 0; it < 20; ++it) {
                            V r = v - u - h * fn(t + h, v),
                              sc = (atol + rtol * u.cwiseAbs().cwiseMax(v.cwiseAbs()).array())
                                       .matrix();
                            double nr = (r.array() / sc.array()).matrix().norm() /
                                        std::sqrt(double(v.size()));
                            trace.push_back(J::array({t + h, h, it, nr, r.norm()}));
                            if (nr <= 1)
                                break;
                            req(it < 19, "BE Newton iteration limit");
                            M j = M::Identity(v.size(), v.size()) - h * jac(
                                                                            [&](const V& x) {
                                                                                return fn(t + h, x);
                                                                            },
                                                                            v);
                            V dv = j.partialPivLu().solve(-r);
                            bool accepted = false;
                            for (int ls = 0; ls < 15; ++ls) {
                                V trial = v + std::ldexp(1., -ls) * dv;
                                try {
                                    if ((trial - u - h * fn(t + h, trial)).norm() < r.norm()) {
                                        v = trial;
                                        accepted = true;
                                        break;
                                    }
                                } catch (const Err&) {
                                }
                            }
                            req(accepted, "BE line search failed");
                        }
                    } else if (method == "Radau") {
                        double tol = std::max(1e-13, atol * .02);
                        V big = radauStep(fn, t, u, h, tol);
                        V half = radauStep(fn, t, u, h / 2, tol);
                        v = radauStep(fn, t + h / 2, half, h / 2, tol);
                        err = (v - big) / 31;
                    } else if (fixed) {
                        J tab = tableau(method);
                        int ns = int(tab["b"].size());
                        M a = M::Zero(ns, ns);
                        V weights(ns);
                        auto val = [](J j) {
                            if (j.is_number())
                                return j.get<double>();
                            std::string s = j.get<std::string>();
                            size_t k = s.find('/');
                            return k == std::string::npos
                                       ? std::stod(s)
                                       : std::stod(s.substr(0, k)) / std::stod(s.substr(k + 1));
                        };
                        for (int i = 0; i < ns; ++i) {
                            weights[i] = val(tab["b"][i]);
                            for (int j = 0; j < ns; ++j)
                                a(i, j) = val(tab["A"][i][j]);
                        }
                        std::vector<V> k;
                        for (int i = 0; i < ns; ++i) {
                            V stage = u;
                            for (int j = 0; j < i; ++j)
                                stage += h * a(i, j) * k[j];
                            k.push_back(fn(t + h * a.row(i).sum(), stage));
                        }
                        for (int i = 0; i < ns; ++i)
                            v += h * weights[i] * k[i];
                    } else {
                        bool dop = method == "DOP853";
                        int ns = dop ? 12 : 6;
                        std::vector<V> k;
                        for (int i = 0; i < ns; ++i) {
                            V stage = u;
                            for (int j = 0; j < i; ++j)
                                stage += h * (dop ? dpA[i][j] : rkA[i][j]) * k[j];
                            k.push_back(fn(t + h * (dop ? dpC[i] : rkC[i]), stage));
                        }
                        for (int i = 0; i < ns; ++i)
                            v += h * (dop ? dpB[i] : rkB[i]) * k[i];
                        k.push_back(fn(t + h, v));
                        if (dop) {
                            V e3 = V::Zero(u.size()), e5 = e3;
                            for (int i = 0; i <= ns; ++i) {
                                e3 += dpE3[i] * k[i];
                                e5 += dpE5[i] * k[i];
                            }
                            V sc = (atol + rtol * u.cwiseAbs().cwiseMax(v.cwiseAbs()).array())
                                       .matrix();
                            double a = (e5.array() / sc.array()).square().sum(),
                                   b = (e3.array() / sc.array()).square().sum();
                            errNorm = (a + b) > 0 ? h * a / std::sqrt((a + .01 * b) * u.size()) : 0;
                        } else
                            for (int i = 0; i <= ns; ++i)
                                err += h * rkE[i] * k[i];
                    }
                } catch (const Stop&) {
                    throw;
                } catch (const Err&) {
                    if (fixed)
                        throw;
                    ok = false;
                }
                if (ok && method != "DOP853" && !fixed) {
                    V sc = (atol + rtol * u.cwiseAbs().cwiseMax(v.cwiseAbs()).array()).matrix();
                    errNorm =
                        (err.array() / sc.array()).matrix().norm() / std::sqrt(double(u.size()));
                }
                if (ok && (fixed || errNorm <= 1)) {
                    req(v.allFinite(), "Non-finite time state");
                    t += h;
                    u = v;
                    stepT.push_back(t);
                    req(stepU.size() * size_t(s->size) <= 8000000,
                        "Accepted-state history exceeds memory limit");
                    stepU.push_back(arr(full(*s, t, u)));
                }
                if (!fixed) {
                    int p = method == "DOP853" ? 8 : method == "Radau" ? 6 : 5;
                    double fac = !ok ? .25
                                 : errNorm == 0
                                     ? 5
                                     : std::clamp(.9 * std::pow(errNorm, -1. / p), .2, 5.);
                    h *= fac;
                }
            }
            out.u.push_back(full(*s, out.t[io], u));
        }
    }
    out.info["steps"] = steps;
    out.info["nfev"] = fev;
    out.info["seconds"] = watch.seconds();
    out.diag["newton_trace"] = trace;
    out.diag["linear_trace"] = lin;
    if (method != "BDF") {
        out.diag["step_t"] = stepT;
        out.diag["step_u"] = stepU;
        if (!s->mixed) {
            J defect = J::array();
            for (int k = 1; k < stepT.size(); ++k) {
                double a = stepT[k - 1], b = stepT[k], mid = (a + b) / 2;
                V prev = vec(stepU[k - 1]), next = vec(stepU[k]),
                  state = s->lift(mid, (prev + next) / 2);
                V d = (next - prev) / (b - a) - s->rhs(mid, state);
                defect.push_back(J{{"t", mid}, {"l2", arr(s->norm(d))}});
            }
            out.diag["ode_defect"] = defect;
            out.info["residual_definition"] =
                "Slope of piecewise-linear accepted-state interpolant minus RHS at interval midpoints; diagnostic, not an error bound";
        }
    }
    if (method == "BE" || method == "FE" || method == "RK2" || method == "RK4" ||
        method == "SSPRK3") {
        out.scheme = scheme(method);
        try {
            double mh = 0;
            for (int i = 1; i < stepT.size(); ++i)
                mh = std::max(mh, stepT[i].get<double>() - stepT[i - 1].get<double>());
            out.scheme["heat"] = heat(*s->mdl, method, mh);
        } catch (const Err& e) {
            out.scheme["heat"] = J{{"status", "not_applicable"}, {"reason", e.what()}};
        }
    }
    out.diagnose();
    return out;
}
}
