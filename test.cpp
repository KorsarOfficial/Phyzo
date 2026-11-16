#include "phz.hpp"

namespace phz {
struct Random {
    std::mt19937_64 rng;

    explicit Random(uint64_t seed) : rng(seed) {
    }

    double f(double a = -1, double b = 1) {
        return std::uniform_real_distribution<double>(a, b)(rng);
    }

    int i(int a, int b) {
        return std::uniform_int_distribution<int>(a, b)(rng);
    }

    V v(int n, double a = -1, double b = 1) {
        V x(n);
        for (auto& k : x)
            k = f(a, b);
        return x;
    }
};

static void close(const V& a, const V& b, double tol, const std::string& what) {
    req(a.size() == b.size() && a.allFinite() && b.allFinite(), "Invalid " + what + " comparison");
    double err = (a - b).size() ? (a - b).cwiseAbs().maxCoeff() : 0;
    req(err <= tol, "Mismatch " + what + ": " + txt(err) + " > " + txt(tol));
}

static void close(const M& a, const M& b, double tol, const std::string& what) {
    req(a.rows() == b.rows() && a.cols() == b.cols(), "Matrix shape mismatch: " + what);
    V av = Eigen::Map<const V>(a.data(), a.size()), bv = Eigen::Map<const V>(b.data(), b.size());
    close(av, bv, tol, what);
}

static double fac(int n) {
    double a = 1;
    for (int i = 2; i <= n; ++i)
        a *= i;
    return a;
}

static double moment(const VI& powers) {
    double x = fac(int(powers.size()) - 1);
    int sum = 0;
    for (int k : powers) {
        x *= fac(k);
        sum += k;
    }
    return x / fac(int(powers.size()) - 1 + sum);
}

static double polyIntegral(const V& u, int p, VI extra = {}) {
    if (extra.empty())
        extra = VI(u.size());
    VI powers(u.size());
    double sum = 0;
    std::function<void(int, int)> rec = [&](int i, int left) {
        if (i == u.size() - 1) {
            powers[i] = left;
            double c = fac(p);
            VI alpha = extra;
            for (int k = 0; k < u.size(); ++k) {
                c *= std::pow(u[k], powers[k]) / fac(powers[k]);
                alpha[k] += powers[k];
            }
            sum += c * moment(alpha);
            return;
        }
        for (int k = 0; k <= left; ++k) {
            powers[i] = k;
            rec(i + 1, left - k);
        }
    };
    rec(0, p);
    return sum;
}

static std::pair<M, IM> geometry(Random& r, int d, int g = 0) {
    if (!g)
        g = d;
    M a(g, d);
    for (int i = 0; i < g; ++i)
        for (int j = 0; j < d; ++j)
            a(i, j) = r.f();
    Eigen::HouseholderQR<M> qr(a);
    M q = qr.householderQ() * M::Identity(g, d);
    double scale = std::pow(10., r.f(-1, 1));
    for (int k = 0; k < d; ++k)
        q.col(k) *= scale * r.f(.5, 2);
    M pts(d + 1, g);
    V base = r.v(g, -scale, scale);
    pts.row(0) = base;
    for (int k = 0; k < d; ++k)
        pts.row(k + 1) = base + q.col(k);
    IM c(1, d + 1);
    VI order(d + 1);
    std::iota(order.begin(), order.end(), 0);
    std::shuffle(order.begin(), order.end(), r.rng);
    for (int i = 0; i <= d; ++i)
        c(0, i) = order[i];
    return {pts, c};
}

static J base(int d = 1, std::string type = "fe") {
    return J{
        {"dim", d},
        {"fld", {"u"}},
        {"eq", {{"u", "lap(u)"}}},
        {"ini", {{"u", 0}}},
        {"disc", type == "sp" ? J{{"type", "sp"}, {"n", VI(d, 2)}} : J{{"type", "fe"}, {"q", 6}}},
        {"dom", {{"type", "box"}, {"n", VI(d, 2)}}},
        {"stat", {{"tol", 1e-10}, {"eig", false}}},
        {"run",
         {{"method", "DOP853"}, {"t", {0, .05}}, {"n", 3}, {"rtol", 1e-10}, {"atol", 1e-12}}}};
}

static std::string affine(const V& a, double b) {
    std::string s = txt(b);
    for (int k = 0; k < a.size(); ++k)
        s += "+(" + txt(a[k]) + ")*" + std::string(1, "xyz"[k]);
    return s;
}

static std::pair<double, M> reference(const M& p, const IM& c) {
    int d = int(c.cols()) - 1, g = int(p.cols());
    M e(g, d);
    for (int i = 0; i < d; ++i)
        e.col(i) = p.row(c(0, i + 1)) - p.row(c(0, 0));
    Eigen::JacobiSVD<M> svd(e, Eigen::ComputeThinU | Eigen::ComputeThinV);
    double vol = svd.singularValues().prod() / fac(d);
    M inv = svd.solve(M::Identity(g, g)), grad(d + 1, g);
    grad.row(0) = -inv.colwise().sum();
    grad.bottomRows(d) = inv;
    return {vol, grad};
}

static void geometryCase(Random& r, int k) {
    int d = 1 + k % 3, g = r.i(d, 3);
    auto [p, c] = geometry(r, d, g);
    auto [vol, grad] = reference(p, c);
    Mesh m(p, c);
    M mass = M::Zero(d + 1, d + 1), stiff = mass;
    for (int i = 0; i <= d; ++i)
        for (int j = 0; j <= d; ++j) {
            VI a(d + 1);
            ++a[i];
            ++a[j];
            mass(c(0, i), c(0, j)) = vol * moment(a);
            stiff(c(0, i), c(0, j)) = vol * grad.row(i).dot(grad.row(j));
        }
    close(M(m.mass), mass, 2e-10 * vol, "simplex mass");
    close(M(m.stiff), stiff, 2e-10 * std::max(1., stiff.norm()), "simplex stiffness");
    req(std::abs(m.vol.sum() - vol) < 2e-11 * vol, "Simplex measure");
    V zero = V::Zero(m.n);
    close(
        V(m.stiff * V::Ones(m.n)), zero, 2e-10 * std::max(1., stiff.norm()), "constant nullspace");
}

static void basisCase(Random& r, int k) {
    int d = 1 + k % 3, degree = 1 + (k / 3) % 4;
    auto [p, c] = geometry(r, d);
    Mesh m(p, c);
    Space s(m, degree, degree + 2);
    V a = r.v(d, -.2, .2);
    double b = r.f(-.2, .2);
    V u = (s.p * a + V::Constant(s.n, b)).array().pow(degree),
      z = s.x * a + V::Constant(s.x.rows(), b);
    double scale = std::max(1., u.cwiseAbs().maxCoeff());
    close(V(s.b * u), V(z.array().pow(degree)), 2e-9 * scale, "polynomial reproduction");
    for (int i = 0; i < d; ++i)
        close(V(s.g[i] * u),
              V(degree * a[i] * z.array().pow(degree - 1)),
              3e-8 * scale * std::max(1., m.grad[0].norm()),
              "polynomial gradient");
    V at(d + 1);
    for (int i = 0; i <= d; ++i)
        at[i] = p.row(c(0, i)).dot(a) + b;
    double exact = m.vol[0] * polyIntegral(at, 2 * degree), value = u.dot(s.mass * u);
    req(std::abs(value - exact) < 3e-8 * std::max(1., std::abs(exact)), "Polynomial L2 integral");
}

static void nonlinearCase(Random& r, int k) {
    int d = 1 + k % 3, degree = 2 + (k / 3) % 3;
    auto [p, c] = geometry(r, d);
    auto [vol, grad] = reference(p, c);
    V nodal = r.v(d + 1, -1.2, 1.2), u(d + 1);
    for (int i = 0; i <= d; ++i)
        u[i] = nodal[c(0, i)];
    double a = r.f(.1, 2), b = r.f(.1, 2), gamma = r.f(-2, 2);
    J j = base(d);
    j["dom"] = J{{"type", "mesh"}, {"p", arr(p)}, {"c", arr(c)}};
    j["par"] = J{{"a", a}, {"b", b}, {"g", gamma}};
    j["eq"]["u"] = "divg(a+b*u**2,u)+g*u**" + std::to_string(degree);
    System s(j);
    V exact = V::Zero(d + 1);
    M deriv = M::Zero(d + 1, d + 1);
    V dot = grad * (u.transpose() * grad).transpose();
    double mobility = a + b * polyIntegral(u, 2);
    for (int i = 0; i <= d; ++i) {
        VI powers(d + 1);
        ++powers[i];
        exact[c(0, i)] = vol * (gamma * polyIntegral(u, degree, powers) - dot[i] * mobility);
        for (int kk = 0; kk <= d; ++kk) {
            VI pp = powers;
            ++pp[kk];
            VI vj(d + 1);
            ++vj[kk];
            deriv(c(0, i), c(0, kk)) = vol * (gamma * degree * polyIntegral(u, degree - 1, pp) -
                                              grad.row(i).dot(grad.row(kk)) * mobility -
                                              dot[i] * 2 * b * polyIntegral(u, 1, vj));
        }
    }
    auto [actual, jacobian] = s.assemble(0, nodal);
    close(actual, exact, 3e-9 * std::max(1., exact.norm()), "nonlinear weak load");
    close(M(jacobian), deriv, 3e-9 * std::max(1., deriv.norm()), "nonlinear weak Jacobian");
    close(V(s.mass * s.rhs(0, nodal)),
          exact,
          3e-9 * std::max(1., exact.norm()),
          "projected nonlinear load");
}

static void spectralCase(Random& r, int k) {
    int d = 1 + k % 3, degree = 2 + k % 2;
    J j = base(d, "sp");
    j["eq"]["u"] = "0.17*lap(u)+u**" + std::to_string(degree) + "-0.3*u";
    System s(j);
    V u = r.v(s.n, -.6, .6);
    int nq = 1, q = 17;
    for (int i = 0; i < d; ++i)
        nq *= q;
    M b(nq, s.n);
    for (int i = 0; i < nq; ++i) {
        VI ix(d);
        int at = i;
        for (int kk = d - 1; kk >= 0; --kk) {
            ix[kk] = at % q;
            at /= q;
        }
        for (int mode = 0; mode < s.n; ++mode) {
            double v = 1;
            for (int kk = 0; kk < d; ++kk)
                v *= std::cos(pi * s.modes[mode][kk] * (ix[kk] + .5) / q);
            b(i, mode) = v;
        }
    }
    V v = b * u, exact = b.transpose() * (v.array().pow(degree) - .3 * v.array()).matrix() / nq;
    for (int i = 0; i < s.n; ++i) {
        double norm = 1, lambda = 0;
        for (int kk : s.modes[i]) {
            norm *= kk ? .5 : 1;
            lambda += pi * pi * kk * kk;
        }
        exact[i] = exact[i] / norm - .17 * lambda * u[i];
    }
    close(s.rhs(0, u), exact, 2e-10, "dealiased spectral projection");
}

static void timeCase(Random& r, int k) {
    int d = 1 + k % 3;
    J j = base(d, "sp");
    std::string method = std::vector<std::string>{
        "FE", "BE", "RK2", "SSPRK3", "RK4", "BDF", "Radau", "DOP853"}[k % 8];
    j["eq"]["u"] = "0.1*lap(u)";
    int n = 1 << d;
    V u = r.v(n, -1, 1);
    j["ini"]["u"] = arr(u);
    double rho = .1 * pi * pi * d, T = r.f(.05, .4) / rho, h = T / 16;
    j["run"] = J{{"method", method},
                 {"t", {0, T}},
                 {"n", 5},
                 {"max_step", h},
                 {"rtol", 1e-10},
                 {"atol", 1e-12}};
    Result result = execute(j);
    int order = method == "FE" || method == "BE" ? 1
                : method == "RK2"                ? 2
                : method == "SSPRK3"             ? 3
                : method == "RK4"                ? 4
                                                 : 0;
    double tol =
        order ? 2 * std::max(1., u.norm()) * std::pow(rho * h, order) * rho * T + 2e-9 : 3e-7;
    for (int it = 0; it < result.t.size(); ++it) {
        V exact(n);
        for (int i = 0; i < n; ++i) {
            double lam = 0;
            for (int kk : result.sys->modes[i])
                lam += .1 * pi * pi * kk * kk;
            exact[i] = u[i] * std::exp(-lam * result.t[it]);
        }
        close(result.u[it], exact, tol, "time integration " + method);
    }
    req(result.u.back().squaredNorm() <= u.squaredNorm() + 1e-8, "Heat energy grew");
}

static void boundaryCase(Random& r, int k) {
    int d = 1 + k % 3;
    J j = base(d);
    double D = r.f(.1, 1), speed = r.f(), b = r.f();
    V a = r.v(d);
    std::string ini = affine(a, b), exact = ini + "+(" + txt(speed) + ")*t";
    J patches = J::array({J{{"type", "dirichlet"}, {"on", "x0"}, {"val", exact}},
                          J{{"type", "neumann"}, {"on", "x1"}, {"val", -D * a[0]}}});
    for (int i = 1; i < d; ++i) {
        patches.push_back(J{{"type", "robin"},
                            {"on", std::string(1, "xyz"[i]) + "0"},
                            {"coef", 2},
                            {"val", exact + "-(" + txt(D * a[i] / 2) + ")"}});
        patches.push_back(
            J{{"type", "neumann"}, {"on", std::string(1, "xyz"[i]) + "1"}, {"val", -D * a[i]}});
    }
    j["par"] = J{{"D", D}};
    j["eq"]["u"] = "D*lap(u)+(" + txt(speed) + ")";
    j["ini"]["u"] = ini;
    j["bc"]["u"] = patches;
    j["run"]["method"] = "BE";
    j["run"]["max_step"] = .025;
    auto out = execute(j);
    for (int it = 0; it < out.t.size(); ++it) {
        V expected = out.sys->mesh->p * a + V::Constant(out.sys->n, b + speed * out.t[it]);
        close(out.u[it], expected, 3e-8, "mixed boundary manufactured solution");
    }
}

static void equilibriumCase(Random& r, int k) {
    int d = 1 + k % 3;
    J j = base(d);
    double mean = r.f(-2, 2);
    if (k < 30) {
        V a = r.v(d, -.4, .4);
        std::string exact = affine(a, mean);
        j["eq"]["u"] = "lap(u)-u-u**3+(" + exact + ")+(" + exact + ")**3";
        j["bc"]["u"] = J{{"type", "dirichlet"}, {"val", exact}};
        auto out = execute(j, true);
        close(out.u[0],
              V(out.sys->mesh->p * a + V::Constant(out.sys->n, mean)),
              2e-7,
              "nonlinear stationary field");
    } else {
        j["ini"]["u"] = "1+cos(pi*x)";
        j["stat"]["mass"]["u"] = mean;
        auto out = execute(j, true);
        close(out.u[0], V::Constant(out.sys->n, mean), 2e-7, "mass constrained equilibrium");
    }
}

static void flowCase(Random& r, int k) {
    int d = 2 + k % 2, axis = r.i(0, d - 1), cross = (axis + 1) % d;
    double amp = r.f(-.7, .7), nu = r.f(.3, 1.3), offset = r.f(-.2, .2);
    std::string velocity = "a*" + std::string(1, "xyz"[cross]) + "*(1-" +
                           std::string(1, "xyz"[cross]) + ")+b",
                pressure = "-2*nu*a*(" + std::string(1, "xyz"[axis]) + "-0.5)", initial = "[",
                boundary = "[";
    for (int i = 0; i < d; ++i) {
        if (i) {
            initial += ",";
            boundary += ",";
        }
        initial += i == axis ? velocity : "0";
        boundary += i == axis ? (k >= 20 ? "(1+t)*(" + velocity + ")" : velocity) : "0";
    }
    initial += "]";
    boundary += "]";
    J j{{"dim", d},
        {"fld", {{"u", VI{d}}, {"p", VI{}}}},
        {"par", {{"a", amp}, {"b", offset}, {"nu", nu}}},
        {"coef", {{"f", VI{d}}}},
        {"src",
         {{"f",
           k >= 20  ? initial
           : d == 2 ? "[0,0]"
                    : "[0,0,0]"}}},
        {"pde", {{"u", "dt(u)+adv(u,u)=-grad(p)+nu*lap(u)+f"}, {"p", "div(u)=0"}}},
        {"ini", {{"u", k >= 20 ? J(initial) : J(0)}, {"p", k >= 20 ? J(pressure) : J(0)}}},
        {"bc", {{"u", {{"type", "dirichlet"}, {"val", boundary}}}}},
        {"dom", {{"type", "box"}, {"n", VI(d, 2)}}},
        {"disc", {{"type", "th"}, {"q", 4}}},
        {"stat", {{"tol", 1e-10}}},
        {"run", {{"method", "BE"}, {"t", {0, .1}}, {"n", 3}, {"max_step", .05}}}};
    auto out = execute(j, k < 20);
    auto& s = *out.sys;
    for (int it = 0; it < out.t.size(); ++it) {
        double factor = k >= 20 ? 1 + out.t[it] : 1;
        for (int f = 0; f < d; ++f) {
            auto sp = s.spaces[f];
            V ex = V::Zero(sp->n);
            if (f == axis)
                ex = factor *
                     (amp * sp->p.col(cross).array() * (1 - sp->p.col(cross).array()) + offset)
                         .matrix();
            close(V(out.u[it].segment(s.offset[f], sp->n)), ex, 2e-7, "mixed velocity");
        }
        V p = -2 * factor * nu * amp * (s.spaces.back()->p.col(axis).array() - .5).matrix();
        close(V(out.u[it].tail(p.size())), p, 2e-6, "mixed pressure");
    }
}

static void refinementCase(Random& r, int k) {
    if (k < 30) {
        int d = 1 + k % 3;
        auto [p, c] = geometry(r, d);
        Mesh old(p, c);
        auto [fine, parent] = old.refine();
        M a = M::Ones(d + 1, d + 1);
        a.rightCols(d) = p;
        M b = M::Ones(fine.n, d + 1);
        b.rightCols(d) = fine.p;
        M transfer = b * a.inverse();
        close(M(transfer.transpose() * fine.mass * transfer),
              M(old.mass),
              3e-9 * old.vol.sum(),
              "refinement mass pullback");
        close(M(transfer.transpose() * fine.stiff * transfer),
              M(old.stiff),
              3e-9 * std::max(1., old.stiff.norm()),
              "refinement stiffness pullback");
    } else {
        int d = 1 + k % 2, n = r.i(4, 6);
        double amp = r.f(.3, 2);
        std::vector<double> errors;
        for (int mul : {1, 2}) {
            J j = base(d);
            std::string ex = txt(amp) + "*sin(pi*x)" + (d == 2 ? "*sin(pi*y)" : "");
            j["eq"]["u"] = "lap(u)+" + std::to_string(d) + "*pi**2*(" + ex + ")";
            j["bc"]["u"] = J{{"type", "dirichlet"}, {"val", 0}};
            j["dom"]["n"] = VI(d, n * mul);
            auto result = execute(j, true);
            auto& s = *result.sys;
            V exact = amp * (pi * s.x.col(0).array()).sin().matrix();
            if (d == 2)
                exact.array() *= (pi * s.x.col(1).array()).sin();
            V err = s.b * result.u[0] - exact;
            errors.push_back(std::sqrt(s.w.dot(err.array().square().matrix())));
        }
        req(errors[1] < .4 * errors[0], "Mesh refinement did not reduce L2 error by 2.5");
    }
}

static void expressionCase(Random& r, int k) {
    if (k < 40) {
        int a = r.i(1, 12), b = r.i(1, 12);
        J tab{{"A", {{0, 0}, {std::to_string(a) + "/" + std::to_string(b), 0}}},
              {"b",
               {std::to_string(2 * a - b) + "/" + std::to_string(2 * a),
                std::to_string(b) + "/" + std::to_string(2 * a)}}};
        auto out = scheme("custom", tab, 4);
        req(out["order"]["value"] == 2 && out["order"]["exact"].get<bool>(),
            "Exact RK2 family order");
    } else {
        double x = r.f(-2, 2), y = r.f(-2, 2), a = r.f(), b = r.f();
        int p = r.i(2, 6);
        Parser parser({{"x", sym("x")}, {"y", sym("y")}});
        E e = parser.parse(std::string("(") + txt(a) + ")*sin(x)+(" + txt(b) +
                           ")*cos(y)+log(1+x**2)+sqrt(1+y**2)+x*y**" + std::to_string(p) +
                           "+exp(x/5)");
        double actual = eval(e, {{"x", x}, {"y", y}}),
               expected = a * std::sin(x) + b * std::cos(y) + std::log1p(x * x) +
                          std::sqrt(1 + y * y) + x * std::pow(y, p) + std::exp(x / 5);
        req(std::abs(actual - expected) < 2e-12 * std::max(1., std::abs(expected)),
            "Independent elementary expression");
    }
}

static void invalidCase(Random& r, int k) {
    bool rejected = false;
    try {
        if (k < 20) {
            int d = 1 + k % 3;
            auto [p, c] = geometry(r, d);
            if (k % 4 == 0)
                p.row(1) = p.row(0);
            else if (k % 4 == 1)
                p(0, 0) = NAN;
            else if (k % 4 == 2)
                c(0, 0) = d + 2;
            else {
                IM cc(2, d + 1);
                cc.row(0) = c;
                cc.row(1) = c;
                c = cc;
            }
            Mesh m(p, c);
        } else if (k < 40) {
            std::string s = std::vector<std::string>{"u**300",
                                                     "u.__class__",
                                                     "u[1]",
                                                     "f999(u)",
                                                     "1/0",
                                                     "sqrt(-1)",
                                                     "lambda u: u",
                                                     "{u:1}",
                                                     "u if 1 else 0",
                                                     "[u]*[u]"}[k % 10];
            Parser p({{"u", sym("u")}});
            p.parse(s);
        } else {
            J j = base();
            if (k % 5 == 0)
                j["run"]["t"] = {1, 0};
            if (k % 5 == 1)
                j["run"]["n"] = 2.5;
            if (k % 5 == 2)
                j["run"]["rtol"] = -r.f(.1, 1);
            if (k % 5 == 3)
                j["run"]["atol"] = 0;
            if (k % 5 == 4)
                j["run"]["max_step"] = -1;
            Model m(j);
        }
    } catch (const Err&) {
        rejected = true;
    }
    req(rejected, "Invalid input was accepted");
}

static void archiveCase(Random& r, int k, const fs::path& root) {
    J j = base(1 + k % 3, k % 2 ? "sp" : "fe");
    j["eq"]["u"] = "0.1*lap(u)-u";
    j["ini"]["u"] = txt(r.f()) + "+(" + txt(r.f()) + ")*cos(pi*x)";
    J before = j;
    auto result = execute(j);
    fs::path path = result.archive(root);
    req(j == before, "Archive mutated request");
    req(verify(path)["ok"].get<bool>(), "Valid archive rejected");
    auto loaded = restore(path);
    req(loaded.t == result.t, "Archive times changed");
    for (int i = 0; i < result.u.size(); ++i)
        close(loaded.u[i], result.u[i], 0, "archive field roundtrip");
    if (k < 5) {
        auto repeated = replay(path);
        close(repeated.u.back(), result.u.back(), 2e-9, "archive replay");
    }
    fs::path file = path / std::vector<std::string>{
                               "state.json", "model.json", "diag.json", "src/phz.hpp"}[k % 4];
    std::string bytes = read(file);
    bytes[r.i(0, int(bytes.size()) - 1)] ^= 1;
    write(file, bytes);
    req(!verify(path)["ok"].get<bool>(), "Archive corruption not detected");
}

int selftest(const fs::path& out, int count, unsigned seed) {
    req(count >= 1 && count <= 1000, "Test count must be 1..1000");
    std::array<int, 12> sizes{120, 120, 140, 100, 100, 80, 60, 40, 60, 80, 60, 40};
    std::array<std::string, 12> names{"geometry",
                                      "basis",
                                      "nonlinear",
                                      "spectral",
                                      "time",
                                      "boundary",
                                      "equilibrium",
                                      "flow",
                                      "refinement",
                                      "expression",
                                      "invalid",
                                      "archive"};
    J cases = J::array();
    int failed = 0;
    auto start = std::chrono::steady_clock::now();
    fs::path scratch =
        out.parent_path() /
        ("fixtures_" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    for (int ii = 0; ii < count; ++ii) {
        int id = count == 1000 ? ii : int(int64_t(ii) * 1000 / count), family = 0, k = id;
        while (k >= sizes[family])
            k -= sizes[family++];
        Random rng(uint64_t(seed) ^ ((uint64_t(family) + 1) * 0x9e3779b97f4a7c15ULL) ^
                   ((uint64_t(k) + 1) * 0xbf58476d1ce4e5b9ULL));
        auto t = std::chrono::steady_clock::now();
        J row{
            {"id", id}, {"family", names[family]}, {"case", k}, {"seed", seed}, {"status", "pass"}};
        try {
            switch (family) {
                case 0:
                    geometryCase(rng, k);
                    break;
                case 1:
                    basisCase(rng, k);
                    break;
                case 2:
                    nonlinearCase(rng, k);
                    break;
                case 3:
                    spectralCase(rng, k);
                    break;
                case 4:
                    timeCase(rng, k);
                    break;
                case 5:
                    boundaryCase(rng, k);
                    break;
                case 6:
                    equilibriumCase(rng, k);
                    break;
                case 7:
                    flowCase(rng, k);
                    break;
                case 8:
                    refinementCase(rng, k);
                    break;
                case 9:
                    expressionCase(rng, k);
                    break;
                case 10:
                    invalidCase(rng, k);
                    break;
                case 11:
                    archiveCase(rng, k, scratch);
                    break;
            }
        } catch (const std::exception& e) {
            row["status"] = "fail";
            row["error"] = e.what();
            ++failed;
            std::cout << "FAIL " << names[family] << "/" << k << ": " << e.what() << std::endl;
        }
        row["seconds"] =
            std::chrono::duration<double>(std::chrono::steady_clock::now() - t).count();
        cases.push_back(row);
        if ((ii + 1) % 25 == 0 || ii + 1 == count) {
            save(out,
                 J{{"seed", seed},
                   {"completed", ii + 1},
                   {"requested", count},
                   {"passed", ii + 1 - failed},
                   {"failed", failed},
                   {"cases", cases}});
            std::cout << ii + 1 << "/" << count << " cases, " << failed << " failures" << std::endl;
        }
    }
    J report{
        {"seed", seed},
        {"tests", count},
        {"passed", count - failed},
        {"failed", failed},
        {"skipped", 0},
        {"seconds",
         std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count()},
        {"cases", cases},
        {"scope",
         "Seeded native tests with independent mathematical references; not a proof for arbitrary PDEs"}};
    save(out, report);
    std::cout << report["passed"] << " passed, " << failed
              << " failed; report=" << fs::absolute(out).string() << std::endl;
    return failed ? 2 : 0;
}
}
