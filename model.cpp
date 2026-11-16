#include "phz.hpp"
#include <regex>

namespace phz {
static double parValue(const J& j) {
    return j.is_string() ? eval(Parser().parse(j), {}) : number(j, "parameter");
}

static void nameCheck(const std::string& s) {
    static const std::regex r("[A-Za-z][A-Za-z0-9_]{0,23}");
    req(std::regex_match(s, r) && s.find("__") == std::string::npos, "Invalid name: " + s);
    req(!std::set<std::string>{"x",
                               "y",
                               "z",
                               "t",
                               "pi",
                               "n",
                               "grad",
                               "div",
                               "lap",
                               "dt",
                               "intg",
                               "adv",
                               "inner",
                               "dot",
                               "outer",
                               "eye",
                               "trans",
                               "sym",
                               "sin",
                               "cos",
                               "exp",
                               "log",
                               "sqrt",
                               "divg"}
             .count(s),
        "Reserved name: " + s);
}

PDE::PDE(J j) : input(std::move(j)) {
    dim = integer(input.value("dim", J(2)), 1, 3, "dim");
    J f = input.at("fld");
    if (f.is_array()) {
        for (auto n : f) {
            std::string s = n.get<std::string>();
            fields.push_back(s);
            ns[s] = sym(s, {}, "fld");
        }
    } else {
        req(f.is_object() && !f.empty(), "Invalid fields");
        for (auto i = f.begin(); i != f.end(); ++i) {
            VI sh;
            req(i.value().is_array() && i.value().size() <= 2, "Invalid field shape");
            for (auto v : i.value())
                sh.push_back(integer(v, 1, 3, "shape"));
            fields.push_back(i.key());
            ns[i.key()] = sym(i.key(), sh, "fld");
        }
    }
    for (auto group : {"par", "coef"}) {
        J d = input.value(group, J::object());
        req(d.is_object(), "Invalid declaration");
        for (auto i = d.begin(); i != d.end(); ++i) {
            req(!ns.count(i.key()), "Duplicate name: " + i.key());
            VI sh;
            if (std::string(group) == "coef") {
                req(i.value().is_array(), "Invalid coefficient shape");
                for (auto v : i.value())
                    sh.push_back(integer(v, 1, 3, "shape"));
            } else
                parValue(i.value());
            ns[i.key()] = sym(i.key(), sh, std::string(group) == "coef" ? "coef" : "sym");
        }
    }
    for (auto [s, e] : ns)
        nameCheck(s);
    req(fields.size() >= 1 && fields.size() <= 8, "Expected 1..8 fields");
    for (auto s : {"x", "y", "z", "t"})
        ns[s] = sym(s);
    req(input.contains("pde") && input["pde"].is_object() && input["pde"].size() == fields.size(),
        "One PDE required per field");
    Parser p(ns, dim);
    for (auto s : fields) {
        std::string ex = input["pde"].at(s).get<std::string>();
        size_t k = ex.find('=');
        req(k != std::string::npos && ex.find('=', k + 1) == std::string::npos,
            "PDE requires one equals sign");
        E a = p.parse(ex.substr(0, k)), b = p.parse(ex.substr(k + 1));
        E e = add(a, neg(b));
        if (zero(e) && e->sh.empty())
            e = std::make_shared<Node>(Node{"zero", "", 0, ns[s]->sh, {}});
        req(e->sh == ns[s]->sh, "PDE equation shape mismatch");
        eq[s] = e;
    }
}

J PDE::weakForm() const {
    J out = J::object();
    for (auto s : fields) {
        E v = sym("v_" + s, ns.at(s)->sh, "test");
        auto [w, b] = weak(eq.at(s), v, dim);
        J patches = J::array(), essential = J::array();
        std::set<std::string> regs;
        J bc = input.value("bc", J::object()).value(s, J::array());
        if (bc.is_object())
            bc = J::array({bc});
        Parser p(ns, dim);
        p.ns["n"] = sym("n", {dim});
        for (auto d : bc) {
            keys(d, "type on val coef", "boundary");
            std::string kind = d.at("type"), tag = d.value("on", "all");
            req(!regs.count(tag) && !regs.count("all") && (tag != "all" || regs.empty()),
                "Overlapping symbolic boundary regions");
            regs.insert(tag);
            E val = p.parse(d.value("val", J(0)));
            req(!has(val, {"fld", "dt", "grad", "div", "lap", "intg", "adv"}),
                "Boundary data must be prescribed");
            if (zero(val) && val->sh.empty())
                val = std::make_shared<Node>(Node{"zero", "", 0, v->sh, {}});
            req(val->sh == v->sh, "Boundary shape mismatch");
            E integrand = num(0);
            if (kind == "dirichlet")
                essential.push_back(
                    J{{"region", tag}, {"trial", s + "=" + str(val)}, {"test", v->name + "=0"}});
            else {
                req(kind == "neumann" || kind == "flux" || kind == "traction" || kind == "robin",
                    "Unsupported boundary type");
                req(!zero(b) || (zero(val) && (kind == "neumann" || kind == "flux")),
                    "Flux requires spatial boundary operator");
                if (kind == "robin") {
                    E k = p.parse(d.at("coef"));
                    req(k->sh.empty(), "Robin coefficient must be scalar");
                    val = mul(k, add(ns.at(s), neg(val)));
                }
                if (kind == "traction")
                    val = neg(val);
                integrand = op("inner", {val, v}, dim);
            }
            patches.push_back(J{{"region", tag},
                                {"type", kind},
                                {"integrand", str(integrand)},
                                {"ir", ir(integrand)}});
        }
        if (!regs.count("all") && !zero(b))
            patches.push_back(J{{"region", "unassigned"},
                                {"type", "unassigned"},
                                {"integrand", str(b)},
                                {"ir", ir(b)}});
        out[s] = J{{"input", input["pde"][s]},
                   {"shape", eq.at(s)->sh},
                   {"test", v->name},
                   {"residual", str(eq.at(s))},
                   {"strong_ir", ir(eq.at(s))},
                   {"volume", str(w)},
                   {"volume_ir", ir(w)},
                   {"boundary_raw", str(b)},
                   {"boundary_raw_ir", ir(b)},
                   {"boundary", patches},
                   {"essential", essential}};
    }
    return J{
        {"schema", 1},
        {"dimension", dim},
        {"equations", out},
        {"convention",
         "Integral(volume,Omega)+Integral(boundary,Gamma)=0; grad appends axis, div contracts last axis"},
        {"flux_convention", "positive outward total flux"},
        {"formal_proof_assistant", false}};
}

static E bridge(E e) {
    auto a = e->a;
    if (e->op == "div") {
        E f = a[0];
        if (f->op == "grad")
            return op("lap", {bridge(f->a[0])});
        std::vector<E> factors;
        std::function<void(E)> split = [&](E x) {
            if (x->op == "mul") {
                split(x->a[0]);
                split(x->a[1]);
            } else
                factors.push_back(x);
        };
        split(f);
        E k = num(1), u;
        for (auto x : factors) {
            if (x->op == "grad") {
                req(!u, "Multiple gradients in scalar flux");
                u = bridge(x->a[0]);
            } else
                k = mul(k, bridge(x));
        }
        req(bool(u), "Scalar divergence requires coefficient*grad(field)");
        return op("divg", {k, u});
    }
    if (a.empty())
        return e;
    req(e->op != "grad" && e->op != "adv", "Unsupported scalar numerical operator");
    std::vector<E> b;
    for (auto x : a)
        b.push_back(bridge(x));
    return op(e->op, b);
}

Model::Model(J j) : input(std::move(j)) {
    keys(input,
         "ver nm dim fld alg par coef eq pde src ini bc dom disc run stat chk en asm units",
         "model");
    req(input.value("ver", 1) == 1, "Unsupported model version");
    dim = integer(input.value("dim", J(2)), 1, 3, "dim");
    name = input.value("nm", "experiment");
    disc = input.value("disc", J{{"type", "fe"}});
    keys(disc, "type n q aq", "disc");
    type = disc.value("type", "fe");
    req(type == "fe" || type == "sp" || type == "th", "Unknown discretization");
    dom = input.value("dom", J{{"type", "box"}, {"n", VI(dim, 8)}});
    keys(dom, "type n path p c tags", "dom");
    run = input.value("run", J::object());
    keys(run, "t n method rtol atol max_step timeout", "run");
    J ts = run.value("t", J::array({0, 1}));
    req(ts.is_array() && ts.size() == 2 && number(ts[1], "time") > number(ts[0], "time"),
        "Time interval must increase");
    integer(run.value("n", J(51)), 2, 2000, "run.n");
    for (auto s : {"rtol", "atol", "max_step", "timeout"})
        if (run.contains(s))
            req(number(run[s], s) > 0, "Positive run tolerance/limit required");
    std::string method = run.value("method", type == "th" ? "BE" : "BDF");
    req(std::set<std::string>{"FE", "BE", "RK2", "RK4", "SSPRK3", "BDF", "Radau", "RK45", "DOP853"}
            .count(method),
        "Unsupported integrator");
    if (type == "th")
        req(method == "BE", "Mixed flow requires BE");
    stat = input.value("stat", J::object());
    req(std::set<std::string>{"auto", "dense", "sparse"}.count(stat.value("solver", "auto")),
        "Invalid stationary solver");
    keys(stat, "solver tol max_it max_ls timeout mass eig eig_tol", "stat");
    for (auto s : {"tol", "timeout", "eig_tol"})
        if (stat.contains(s))
            req(number(stat[s], s) > 0, "Positive stationary tolerance required");
    for (auto s : {"max_it", "max_ls"})
        if (stat.contains(s))
            integer(stat[s], 1, 100, s);
    chk = input.value("chk", J::object());
    keys(chk, "mass nonneg energy exact tol div", "chk");
    if (chk.contains("tol"))
        req(number(chk["tol"], "chk.tol") > 0, "Positive check tolerance required");
    bc = input.value("bc", J::object());
    req(bc.is_object(), "bc must be an object");
    en = input.value("en", J());
    alg = input.value("alg", J::array());
    req(alg.is_array(), "alg must be an array");
    J pars = input.value("par", J::object());
    req(pars.is_object(), "par must be an object");
    for (auto it = pars.begin(); it != pars.end(); ++it) {
        nameCheck(it.key());
        par[it.key()] = parValue(it.value());
    }
    for (int k = 0; k < dim; ++k) {
        std::string s = "Lx" + std::to_string(k + 1);
        if (!par.count(s))
            par[s] = 1;
        req(par[s] > 0, "Domain lengths must be positive");
    }
    for (auto [s, f] : par)
        ns[s] = sym(s);
    for (auto s : {"x", "y", "z", "t"})
        ns[s] = sym(s);
    J fld = input.at("fld");
    req(fld.is_array() || fld.is_object(), "Invalid fields");
    if (fld.is_array()) {
        for (auto f : fld)
            fields.push_back(f.get<std::string>());
    } else
        for (auto it = fld.begin(); it != fld.end(); ++it)
            fields.push_back(it.key());
    req(fields.size() >= 1 && fields.size() <= 8, "Expected 1..8 fields");
    std::set<std::string> used;
    for (auto s : fields) {
        nameCheck(s);
        req(!ns.count(s) && used.insert(s).second, "Duplicate field");
        VI sh = fld.is_object() ? fld.at(s).get<VI>() : VI{};
        ns[s] = sym(s, sh, "fld");
        if (type != "th")
            req(sh.empty(), "Scalar solver requires scalar fields");
    }
    for (auto d : alg) {
        keys(d, "nm a src", "alg");
        std::string s = d.at("nm");
        nameCheck(s);
        req(!ns.count(s), "Duplicate algebraic field");
        ns[s] = sym(s, {}, "fld");
    }
    for (auto it = bc.begin(); it != bc.end(); ++it)
        req(used.count(it.key()), "Boundary for unknown field");
    req(int(input.contains("eq")) + int(input.contains("pde")) + int(input.contains("en")) == 1,
        "Specify exactly one of eq,pde,en");
    if (input.contains("pde")) {
        J ip = input;
        ip["par"] = par;
        pde = std::make_unique<PDE>(ip);
        if (type != "th")
            for (auto s : fields) {
                double dt = 0;
                E r = num(0);
                std::function<void(E, E)> split = [&](E e, E k) {
                    if (e->op == "add") {
                        split(e->a[0], k);
                        split(e->a[1], k);
                    } else if (e->op == "mul") {
                        for (int i = 0; i < 2; ++i)
                            if (!has(e->a[i],
                                     {"fld", "coef", "dt", "lap", "grad", "div", "intg"})) {
                                split(e->a[1 - i], mul(k, e->a[i]));
                                return;
                            }
                        r = add(r, mul(k, e));
                    } else if (e->op == "dt") {
                        req(e->a[0]->name == s, "Wrong time derivative field");
                        dt += eval(k, par);
                    } else {
                        req(!has(e, {"dt"}), "Nonlinear time derivative unsupported");
                        r = add(r, mul(k, e));
                    }
                };
                split(pde->eq.at(s), num(1));
                req(std::abs(dt - 1) < 1e-14, "Scalar PDE requires unit dt(field)");
                eq[s] = bridge(neg(r));
            }
    } else if (!en.is_null()) {
        keys(en, "f kap mob kin", "en");
        E f = parse(en.at("f"));
        for (auto s : fields) {
            double k = constant(en.value("kap", J::object()).value(s, J(0))),
                   mob = constant(en.value("mob", J::object()).value(s, J(1)));
            req(k >= 0 && mob > 0, "Require kap>=0 and mob>0");
            E mu = add(diff(f, s), neg(mul(num(k), op("lap", {ns.at(s)}))));
            std::string kin = en.value("kin", J::object()).value(s, "ac");
            req(kin == "ac" || kin == "ch", "Unknown energy dynamics");
            eq[s] = kin == "ch" ? mul(num(mob), op("lap", {mu})) : mul(num(-mob), mu);
        }
    } else {
        J es = input.at("eq");
        req(es.is_object() && es.size() == fields.size(), "One RHS required per field");
        for (auto s : fields)
            eq[s] = parse(es.at(s));
    }
    if (type == "sp") {
        VI nspec = disc.value("n", VI(dim, 4));
        req(nspec.size() == dim, "Mode dimension mismatch");
        int n = 1;
        for (int a : nspec) {
            req(a >= 1 && a <= 24, "Modes must be 1..24");
            n *= a;
        }
        req(n * fields.size() <= 1024, "Spectral state exceeds 1024");
        integer(disc.value("q", J(std::max(16, 3 * *std::max_element(nspec.begin(), nspec.end())))),
                std::max(4, 2 * *std::max_element(nspec.begin(), nspec.end())),
                100,
                "quadrature");
    } else
        integer(disc.value("q", J(4)), type == "th" ? 4 : dim == 3 ? 3 : 2, 20, "quadrature");
    if (disc.contains("aq"))
        quadratureConfig(*this);
}

E Model::parse(const J& j) const {
    return Parser(ns, dim).parse(j);
}

double Model::constant(const J& j) const {
    return eval(parse(j), par);
}

V Model::at(E e, const M& x, double t) const {
    req(e->sh.empty(), "Expected scalar data");
    V a(x.rows());
    auto d = par;
    d["t"] = t;
    for (int i = 0; i < x.rows(); ++i) {
        for (int k = 0; k < 3; ++k)
            d[std::string(1, "xyz"[k])] = k < x.cols() ? x(i, k) : 0;
        a[i] = eval(e, d);
    }
    return a;
}
}
