#include "phz.hpp"
#include <boost/multiprecision/cpp_int.hpp>
#include <regex>

namespace phz {
using Big = boost::multiprecision::cpp_int;
using Rat = boost::multiprecision::cpp_rational;
using Poly = std::vector<Rat>;

static std::string rs(const Rat& r) {
    return r.str();
}

static Rat rational(const J& j) {
    req(j.is_string() || j.is_number(), "Expected rational number");
    std::string s = j.is_string() ? j.get<std::string>() : j.dump();
    req(s.size() <= 256, "Rational literal too long");
    static std::regex rx(
        "[+-]?(?:[0-9]+/[0-9]+|(?:[0-9]+(?:\\.[0-9]*)?|\\.[0-9]+)(?:[eE][+-]?[0-9]+)?)");
    req(std::regex_match(s, rx), "Invalid rational literal");
    int sign = 1;
    if (s[0] == '-' || s[0] == '+') {
        if (s[0] == '-')
            sign = -1;
        s.erase(0, 1);
    }
    auto big = [](const std::string& a) {
        Big n = 0;
        for (char c : a)
            n = n * 10 + (c - '0');
        return n;
    };
    size_t slash = s.find('/');
    Rat r;
    if (slash != std::string::npos) {
        Big den = big(s.substr(slash + 1));
        req(den != 0, "Zero rational denominator");
        r = Rat(big(s.substr(0, slash))) / den;
    } else {
        int exponent = 0;
        size_t e = s.find_first_of("eE");
        if (e != std::string::npos) {
            exponent = std::stoi(s.substr(e + 1));
            req(std::abs(exponent) <= 300, "Rational exponent exceeds limit");
            s.resize(e);
        }
        size_t dot = s.find('.');
        if (dot != std::string::npos) {
            exponent -= int(s.size() - dot - 1);
            s.erase(dot, 1);
        }
        Big n = big(s), scale = 1;
        for (int i = 0; i < std::abs(exponent); ++i)
            scale *= 10;
        r = exponent >= 0 ? Rat(n * scale) : Rat(n) / scale;
    }
    r *= sign;
    Big n = numerator(r), d = denominator(r);
    if (n < 0)
        n = -n;
    req((n == 0 || boost::multiprecision::msb(n) < 1024) && boost::multiprecision::msb(d) < 1024,
        "Rational precision limit");
    return r;
}

static void trim(Poly& p) {
    while (p.size() > 1 && p.back() == 0)
        p.pop_back();
}

static Poly plus(Poly a, const Poly& b, Rat sign = 1) {
    a.resize(std::max(a.size(), b.size()));
    for (int i = 0; i < b.size(); ++i)
        a[i] += sign * b[i];
    trim(a);
    return a;
}

static Poly times(const Poly& a, const Poly& b) {
    Poly c(a.size() + b.size() - 1);
    for (int i = 0; i < a.size(); ++i)
        for (int j = 0; j < b.size(); ++j)
            c[i + j] += a[i] * b[j];
    trim(c);
    return c;
}

static std::pair<Poly, Poly> divrem(Poly a, const Poly& b) {
    req(!(b.size() == 1 && b[0] == 0), "Zero polynomial denominator");
    Poly q(a.size() >= b.size() ? a.size() - b.size() + 1 : 1);
    while (a.size() >= b.size() && !(a.size() == 1 && a[0] == 0)) {
        int k = int(a.size() - b.size());
        Rat c = a.back() / b.back();
        q[k] += c;
        for (int i = 0; i < b.size(); ++i)
            a[i + k] -= c * b[i];
        trim(a);
    }
    trim(q);
    return {q, a};
}

static Poly gcd(Poly a, Poly b) {
    while (!(b.size() == 1 && b[0] == 0)) {
        Poly r = divrem(a, b).second;
        a = b;
        b = r;
    }
    Rat lead = a.back();
    for (auto& x : a)
        x /= lead;
    return a;
}

static Poly det(const std::vector<std::vector<Poly>>& a) {
    VI p(a.size());
    std::iota(p.begin(), p.end(), 0);
    Poly sum{0};
    do {
        int inversions = 0;
        for (int i = 0; i < p.size(); ++i)
            for (int j = i + 1; j < p.size(); ++j)
                inversions += p[i] > p[j];
        Poly v{inversions % 2 ? -1 : 1};
        for (int i = 0; i < p.size(); ++i)
            v = times(v, a[i][p[i]]);
        sum = plus(sum, v);
    } while (std::next_permutation(p.begin(), p.end()));
    return sum;
}

static std::string ps(const Poly& a) {
    std::string s = "0";
    for (int i = 0; i < a.size(); ++i)
        if (a[i] != 0)
            s += "+(" + rs(a[i]) + ")" + (i ? "*z" + (i > 1 ? "**" + std::to_string(i) : "") : "");
    return "(" + s + ")";
}

using P2 = std::map<std::pair<int, int>, Rat>;

static P2 mult(const P2& a, const P2& b) {
    P2 c;
    for (auto [ka, va] : a)
        for (auto [kb, vb] : b)
            c[{ka.first + kb.first, ka.second + kb.second}] += va * vb;
    return c;
}

static P2 absSquare(const Poly& p) {
    P2 re, im;
    for (int n = 0; n < p.size(); ++n) {
        Rat binom = 1;
        for (int j = 0; j <= n; ++j) {
            Rat k = p[n] * binom * ((n - j) % 2 ? -1 : 1);
            if (j % 4 >= 2)
                k = -k;
            (j % 2 ? im : re)[{n - j, j}] += k;
            binom = binom * (n - j) / (j + 1);
        }
    }
    P2 c = mult(re, re);
    for (auto [k, v] : mult(im, im))
        c[k] += v;
    return c;
}

static J cert(const P2& p) {
    J terms = J::array();
    bool ok = true;
    for (auto [k, v] : p)
        if (v != 0) {
            terms.push_back(J{{"x", k.first}, {"y", k.second}, {"coefficient", rs(v)}});
            ok &= v >= 0 && k.second % 2 == 0;
        }
    return J{{"terms", terms}, {"nonnegative_coefficients", ok}};
}

struct Tree {
    VI child;
    int n;
    Big gamma;
};

static std::vector<Tree> trees(int order) {
    std::vector<Tree> all{{{}, 1, Big(1)}};
    for (int n = 2; n <= order; ++n) {
        int last = int(all.size());
        VI ch;
        std::function<void(int, int)> rec = [&](int left, int first) {
            if (left == 0) {
                Big g = n;
                for (int i : ch)
                    g *= all[i].gamma;
                all.push_back({ch, n, g});
                return;
            }
            for (int i = first; i < last; ++i) {
                if (all[i].n > left)
                    break;
                ch.push_back(i);
                rec(left - all[i].n, i);
                ch.pop_back();
            }
        };
        rec(n - 1, 0);
    }
    return all;
}

J tableau(const std::string& s) {
    if (s == "FE")
        return J{{"A", {{0}}}, {"b", {1}}};
    if (s == "BE")
        return J{{"A", {{1}}}, {"b", {1}}};
    if (s == "CN")
        return J{{"A", {{0, 0}, {"1/2", "1/2"}}}, {"b", {"1/2", "1/2"}}};
    if (s == "RK2")
        return J{{"A", {{0, 0}, {1, 0}}}, {"b", {"1/2", "1/2"}}};
    if (s == "SSPRK3")
        return J{{"A", {{0, 0, 0}, {1, 0, 0}, {"1/4", "1/4", 0}}}, {"b", {"1/6", "1/6", "2/3"}}};
    if (s == "RK4")
        return J{{"A", {{0, 0, 0, 0}, {"1/2", 0, 0, 0}, {0, "1/2", 0, 0}, {0, 0, 1, 0}}},
                 {"b", {"1/6", "1/3", "1/3", "1/6"}}};
    if (s == "RADAU2")
        return J{{"A", {{"5/12", "-1/12"}, {"3/4", "1/4"}}}, {"b", {"3/4", "1/4"}}};
    throw Err("No fixed RK tableau for " + s);
}

static double peval(const Poly& p, double x) {
    double a = 0;
    for (int i = int(p.size()) - 1; i >= 0; --i)
        a = a * x + p[i].convert_to<double>();
    return a;
}

static Rat ratAt(const Poly& p, const Rat& x) {
    Rat a = 0;
    for (int i = int(p.size()) - 1; i >= 0; --i)
        a = a * x + p[i];
    return a;
}

static Poly derivative(const Poly& p) {
    Poly q(std::max(size_t(1), p.size() - 1));
    for (int i = 1; i < p.size(); ++i)
        q[i - 1] = i * p[i];
    trim(q);
    return q;
}

static std::vector<Poly> sturm(Poly p) {
    trim(p);
    std::vector<Poly> s{p};
    if (p.size() == 1)
        return s;
    s.push_back(derivative(p));
    while (s.back().size() > 1) {
        Poly r = divrem(s[s.size() - 2], s.back()).second;
        if (r.size() == 1 && r[0] == 0)
            break;
        for (auto& x : r)
            x = -x;
        Rat scale = r.back();
        if (scale < 0)
            scale = -scale;
        for (auto& x : r)
            x /= scale;
        s.push_back(r);
    }
    return s;
}

static int variations(const std::vector<Poly>& s, const Rat& x) {
    int prev = 0, n = 0;
    for (auto p : s) {
        Rat v = ratAt(p, x);
        int sign = v > 0 ? 1 : v < 0 ? -1 : 0;
        if (sign) {
            if (prev && prev != sign)
                ++n;
            prev = sign;
        }
    }
    return n;
}

struct RootInterval {
    Rat lo, hi;
};

static std::vector<RootInterval> positiveRoots(Poly p) {
    trim(p);
    while (p.size() > 1 && p[0] == 0)
        p.erase(p.begin());
    if (p.size() == 1)
        return {};
    p = divrem(p, gcd(p, derivative(p))).first;
    auto chain = sturm(p);
    Rat bound = 2;
    for (int i = 0; i + 1 < p.size(); ++i) {
        Rat a = p[i] / p.back();
        if (a < 0)
            a = -a;
        bound = std::max(bound, Rat(2 + a));
    }
    while (ratAt(p, bound) == 0)
        bound += 1;
    std::vector<RootInterval> roots;
    std::function<void(Rat, Rat, int, int)> split = [&](Rat lo, Rat hi, int count, int depth) {
        if (!count)
            return;
        req(depth < 1024, "Root isolation precision limit");
        checkpoint();
        if (count == 1) {
            if (ratAt(p, hi) == 0) {
                roots.push_back({hi, hi});
                return;
            }
            Rat eps = (1 + hi) / Rat(Big(1) << 80);
            if (hi - lo < eps && ratAt(p, lo) != 0) {
                roots.push_back({lo, hi});
                return;
            }
        }
        Rat mid = (lo + hi) / 2;
        int left = variations(chain, lo) - variations(chain, mid);
        split(lo, mid, left, depth + 1);
        split(mid, hi, count - left, depth + 1);
    };
    split(0, bound, variations(chain, 0) - variations(chain, bound), 0);
    return roots;
}

static J negativeAxis(Poly p, Poly q, Poly stage) {
    for (int i = 1; i < p.size(); i += 2)
        p[i] = -p[i];
    for (int i = 1; i < q.size(); i += 2)
        q[i] = -q[i];
    for (int i = 1; i < stage.size(); i += 2)
        stage[i] = -stage[i];
    Poly gap = plus(times(q, q), times(p, p), -1);
    bool flat = gap.size() == 1 && gap[0] == 0;
    Poly unionPoly = times(stage, q);
    if (!flat)
        unionPoly = times(unionPoly, gap);
    auto roots = positiveRoots(unionPoly);
    Rat sample = roots.empty() ? Rat(1) : roots[0].lo / 2;
    if (sample == 0)
        sample = roots[0].hi / 4;
    J result{{"status", "proved"},
             {"unbounded", true},
             {"radius_approx", J()},
             {"right_closed", false},
             {"certificate", "Sturm isolation with exact rational arithmetic"}};
    if (ratAt(gap, sample) < 0) {
        result["unbounded"] = false;
        result["radius_approx"] = 0;
        result["interval"] = {"0", "0"};
        result["right_closed"] = true;
        return result;
    }
    auto sc = sturm(stage), qc = sturm(q);
    for (int i = 0; i < roots.size(); ++i) {
        auto r = roots[i];
        bool pole = r.lo == r.hi ? (ratAt(stage, r.hi) == 0 || ratAt(q, r.hi) == 0)
                                 : (variations(sc, r.lo) > variations(sc, r.hi) ||
                                    variations(qc, r.lo) > variations(qc, r.hi));
        Rat after = i + 1 < roots.size() ? Rat((r.hi + roots[i + 1].lo) / 2) : Rat(r.hi + 1);
        if (pole || ratAt(gap, after) < 0) {
            result["unbounded"] = false;
            result["radius_approx"] = Rat((r.lo + r.hi) / 2).convert_to<double>();
            result["interval"] = {rs(r.lo), rs(r.hi)};
            result["right_closed"] = !pole;
            result["boundary_reason"] = pole ? "singular stage or stability pole"
                                             : "amplification exceeds one beyond boundary";
            return result;
        }
    }
    return result;
}

J scheme(const std::string& method, J tab, int maxOrder) {
    req(maxOrder >= 1 && maxOrder <= 8, "Tree order must be 1..8");
    if (tab.is_null())
        tab = tableau(method);
    keys(tab, "A b", "tableau");
    req(tab.contains("A") && tab.contains("b") && tab["A"].is_array() && tab["b"].is_array(),
        "Invalid tableau");
    int n = int(tab["b"].size());
    req(n >= 1 && n <= 6 && tab["A"].size() == n, "Tableau requires 1..6 stages");
    std::vector<std::vector<Rat>> a(n, std::vector<Rat>(n));
    std::vector<Rat> b(n);
    for (int i = 0; i < n; ++i) {
        req(tab["A"][i].is_array() && tab["A"][i].size() == n, "Tableau A must be square");
        b[i] = rational(tab["b"][i]);
        for (int j = 0; j < n; ++j)
            a[i][j] = rational(tab["A"][i][j]);
    }
    std::vector<std::vector<Poly>> den(n, std::vector<Poly>(n)), num = den;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            den[i][j] = {i == j ? 1 : 0, -a[i][j]};
            num[i][j] = {i == j ? 1 : 0, b[j] - a[i][j]};
        }
    Poly qStage = det(den), p = det(num), q = qStage, g = gcd(p, q);
    p = divrem(p, g).first;
    q = divrem(q, g).first;
    Rat q0 = q[0];
    for (auto& v : p)
        v /= q0;
    for (auto& v : q)
        v /= q0;
    auto tr = trees(maxOrder);
    std::vector<std::vector<Rat>> weights;
    J cond = J::array();
    int order = maxOrder;
    bool exact = false;
    for (auto tree : tr) {
        std::vector<Rat> w(n, 1);
        for (int ch : tree.child)
            for (int i = 0; i < n; ++i) {
                Rat v = 0;
                for (int j = 0; j < n; ++j)
                    v += a[i][j] * weights[ch][j];
                w[i] *= v;
            }
        Rat lhs = 0;
        for (int i = 0; i < n; ++i)
            lhs += b[i] * w[i];
        Rat rhs = Rat(1) / tree.gamma;
        bool ok = lhs == rhs;
        cond.push_back(J{{"order", tree.n},
                         {"children", tree.child},
                         {"lhs", rs(lhs)},
                         {"rhs", rs(rhs)},
                         {"holds", ok}});
        weights.push_back(w);
        if (!ok) {
            order = std::min(order, tree.n - 1);
            exact = true;
        }
    }
    P2 denAbs = absSquare(q), stageAbs = absSquare(qStage), gap = denAbs;
    for (auto [k, v] : absSquare(p))
        gap[k] -= v;
    J gapCert = cert(gap), denCert = cert(denAbs), stageCert = cert(stageAbs);
    bool proved = gapCert["nonnegative_coefficients"].get<bool>() &&
                  denCert["nonnegative_coefficients"].get<bool>() &&
                  stageCert["nonnegative_coefficients"].get<bool>() && denAbs[{0, 0}] > 0 &&
                  stageAbs[{0, 0}] > 0;
    std::string ast = proved ? "proved" : "unknown";
    J counter = J();
    for (int k = 0; k < 32; ++k) {
        Rat x = -(Big(1) << k);
        auto at = [&](const Poly& p) {
            Rat r = 0;
            for (int i = int(p.size()) - 1; i >= 0; --i)
                r = r * x + p[i];
            return r;
        };
        Rat stage = at(qStage), qp = at(q), pp = at(p);
        if (stage == 0 || qp == 0 || pp * pp > qp * qp) {
            ast = "disproved";
            counter = J{{"z", rs(x)},
                        {"reason",
                         stage == 0 ? "singular stages"
                         : qp == 0  ? "stability pole"
                                    : "amplification exceeds one"}};
            break;
        }
    }
    bool decay = p.size() < q.size();
    std::string lst = ast == "proved"                ? (decay ? "proved" : "disproved")
                      : ast == "disproved" || !decay ? "disproved"
                                                     : "unknown";
    J axis = negativeAxis(p, q, qStage);
    return J{
        {"method", method},
        {"tableau", tab},
        {"stability_function", ps(p) + "/" + ps(q)},
        {"stage_determinant", ps(qStage)},
        {"order", J{{"value", order}, {"exact", exact}, {"conditions", cond}}},
        {"a_stability",
         J{{"status", ast},
           {"modulus_gap", gapCert},
           {"denominator", denCert},
           {"stage_determinant", stageCert},
           {"counterexample", counter}}},
        {"l_stability", J{{"status", lst}, {"limit_zero", decay}}},
        {"negative_real_axis", axis},
        {"assurance",
         J{{"arithmetic", "arbitrary-precision rationals for order and A/L certificates"},
           {"formal_proof_assistant", false},
           {"radius_scope", "Negative-axis boundary has an exact rational isolating interval"}}}};
}

J heat(const Model& m, const std::string& method, double dt) {
    req(m.type != "th" && m.fields.size() == 1, "Heat analysis requires scalar model");
    std::string f = m.fields[0];
    E e = m.eq.at(f);
    double D = 0;
    bool seen = false;
    std::function<void(E, double)> walk = [&](E e, double k) {
        if (zero(e))
            return;
        if (e->op == "add") {
            walk(e->a[0], k);
            walk(e->a[1], k);
            return;
        }
        if (e->op == "mul") {
            for (int i = 0; i < 2; ++i)
                try {
                    double v = eval(e->a[i], m.par);
                    walk(e->a[1 - i], k * v);
                    return;
                } catch (const Err&) {
                }
            throw Err("Heat requires constant diffusivity");
        }
        req(e->op == "lap" && e->a[0]->name == f, "Model is not homogeneous scalar heat equation");
        D += k;
        seen = true;
    };
    walk(e, 1);
    req(seen && D >= 0, "Heat requires D>=0");
    if (m.bc.contains(f)) {
        J bc = m.bc[f];
        if (bc.is_object())
            bc = J::array({bc});
        for (auto b : bc)
            req(b.value("type", "") == "neumann" && m.constant(b.value("val", J(0))) == 0,
                "Heat spectral analysis requires homogeneous Neumann");
    }
    J sch = scheme(method);
    J neg = sch["negative_real_axis"];
    if (m.type == "fe")
        return J{{"status", "conditional"},
                 {"dt_status", "not_evaluated"},
                 {"dt_max", "r/(D*lambda_max(K_ff,M_ff))"},
                 {"assumptions",
                  J::array({"Consistent SPD mass", "No certified assembled eigenvalue bound"})}};
    VI ns = m.disc.value("n", VI(m.dim, 4));
    double rho = 0;
    for (int i = 0; i < m.dim; ++i)
        rho += D * std::pow(pi * (ns[i] - 1) / m.par.at("Lx" + std::to_string(i + 1)), 2);
    double radius = neg["unbounded"].get<bool>()       ? INFINITY
                    : neg["radius_approx"].is_number() ? neg["radius_approx"].get<double>()
                                                       : NAN;
    double mx = rho == 0 ? INFINITY : radius / rho;
    return J{{"status",
              std::isnan(mx) ? "unknown"
              : neg["status"] == "proved" || method == "FE" || method == "RK2"
                  ? "derived"
                  : "numerical_estimate"},
             {"spectral_radius", rho},
             {"dt_max_approx", std::isfinite(mx) ? J(mx) : J()},
             {"unbounded", std::isinf(mx)},
             {"dt", dt < 0 ? J() : J(dt)},
             {"dt_status",
              dt < 0           ? "not_evaluated"
              : std::isnan(mx) ? "unknown"
              : dt <= mx       ? "pass"
                               : "fail"},
             {"scope",
              "Perturbation stability of retained cosine modes; no discretization-error bound"}};
}
}
