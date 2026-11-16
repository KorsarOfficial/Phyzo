#include "phz.hpp"
#include <cctype>

namespace phz {
E num(double x) {
    req(std::isfinite(x) && std::abs(x) <= 1e100, "Invalid numeric constant");
    return std::make_shared<Node>(Node{"num", "", x, {}, {}});
}

E sym(const std::string& s, VI sh, const std::string& kind) {
    return std::make_shared<Node>(Node{kind, s, 0, std::move(sh), {}});
}

bool zero(E e) {
    return e->op == "zero" || (e->op == "num" && e->n == 0);
}

bool has(E e, const std::set<std::string>& ops) {
    if (ops.count(e->op))
        return true;
    for (auto a : e->a)
        if (has(a, ops))
            return true;
    return false;
}

bool depends(E e, const std::string& s) {
    if (e->name == s)
        return true;
    for (auto a : e->a)
        if (depends(a, s))
            return true;
    return false;
}

static int count(const VI& s) {
    return std::accumulate(s.begin(), s.end(), 1, std::multiplies<int>());
}

static E shapedZero(VI s) {
    return s.empty() ? num(0) : std::make_shared<Node>(Node{"zero", "", 0, std::move(s), {}});
}

E op(const std::string& s, std::vector<E> a, int dim) {
    VI sh;
    req(dim >= 1 && dim <= 3, "Invalid expression dimension");
    if (s == "eye") {
        req(a.empty(), "eye takes no arguments");
        sh = {dim, dim};
    } else if (s == "arr") {
        req(!a.empty() && a.size() <= 3, "Array size must be 1..3");
        sh = a[0]->sh;
        for (auto e : a)
            req(e->sh == sh, "Array shape mismatch");
        sh.insert(sh.begin(), int(a.size()));
    } else if (s == "add" || s == "mul" || s == "pow" || s == "inner" || s == "dot" ||
               s == "outer" || s == "adv" || s == "divg") {
        req(a.size() == 2, s + " takes two arguments");
        E x = a[0], y = a[1];
        if (s == "add") {
            if (zero(x) && x->sh.empty())
                return y;
            if (zero(y) && y->sh.empty())
                return x;
            req(x->sh == y->sh, "Addition shape mismatch");
            sh = x->sh;
            if (zero(x))
                return y;
            if (zero(y))
                return x;
            if (x->op == "num" && y->op == "num")
                return num(x->n + y->n);
        }
        if (s == "mul") {
            req(x->sh.empty() || y->sh.empty(), "Use dot, inner or outer for tensors");
            sh = x->sh.empty() ? y->sh : x->sh;
            if (zero(x) || zero(y))
                return shapedZero(sh);
            if (x->op == "num" && x->n == 1)
                return y;
            if (y->op == "num" && y->n == 1)
                return x;
            if (x->op == "num" && y->op == "num")
                return num(x->n * y->n);
        }
        if (s == "pow") {
            req(x->sh.empty() && y->op == "num" && std::abs(y->n) <= 1024,
                "Power requires a scalar and a bounded constant exponent");
            if (y->n == 0)
                return num(1);
            if (y->n == 1)
                return x;
            if (x->op == "num")
                return num(std::pow(x->n, y->n));
        }
        if (s == "inner") {
            req(x->sh == y->sh, "Inner product shape mismatch");
            if (zero(x) || zero(y))
                return num(0);
            if (x->sh.empty())
                return mul(x, y);
        }
        if (s == "dot") {
            req(!x->sh.empty() && !y->sh.empty() && x->sh.back() == y->sh.front(),
                "Dot shape mismatch");
            sh = x->sh;
            sh.pop_back();
            sh.insert(sh.end(), y->sh.begin() + 1, y->sh.end());
        }
        if (s == "outer") {
            sh = x->sh;
            sh.insert(sh.end(), y->sh.begin(), y->sh.end());
        }
        if (s == "adv") {
            req(x->sh == VI{dim} && y->sh.size() <= 1, "Invalid advection shape");
            sh = y->sh;
        }
        if (s == "divg") {
            req(x->sh.empty(), "divg requires scalar mobility");
            sh = y->sh;
        }
    } else {
        req(a.size() == 1, s + " takes one argument");
        sh = a[0]->sh;
        if (s == "grad")
            sh.push_back(dim);
        else if (s == "div") {
            req(!sh.empty() && sh.back() == dim, "Divergence shape mismatch");
            sh.pop_back();
        } else if (s == "trans" || s == "sym") {
            req(sh.size() == 2 && (s != "sym" || sh[0] == sh[1]), "Matrix shape mismatch");
            std::swap(sh[0], sh[1]);
        } else if (s == "exp" || s == "log" || s == "sin" || s == "cos" || s == "sqrt") {
            req(sh.empty(), s + " requires scalar");
            if (a[0]->op == "num") {
                double v = a[0]->n;
                return num(s == "exp"   ? std::exp(v)
                           : s == "log" ? std::log(v)
                           : s == "sin" ? std::sin(v)
                           : s == "cos" ? std::cos(v)
                                        : std::sqrt(v));
            }
        } else
            req(s == "lap" || s == "dt" || s == "intg", "Unknown operator: " + s);
    }
    req(sh.size() <= 3 && count(sh) <= 81, "Tensor expression too large");
    return std::make_shared<Node>(Node{s, "", 0, sh, std::move(a)});
}

E add(E a, E b) {
    return op("add", {a, b});
}

E mul(E a, E b) {
    return op("mul", {a, b});
}

E neg(E a) {
    return mul(num(-1), a);
}

E power(E a, double b) {
    return op("pow", {a, num(b)});
}

E diff(E e, const std::string& s) {
    if (e->op == "num" || e->op == "zero")
        return num(0);
    if (!e->name.empty())
        return num(e->name == s ? 1 : 0);
    auto a = e->a;
    E d = diff(a[0], s);
    if (e->op == "add")
        return add(d, diff(a[1], s));
    if (e->op == "mul")
        return add(mul(d, a[1]), mul(a[0], diff(a[1], s)));
    if (e->op == "pow")
        return mul(mul(a[1], power(a[0], a[1]->n - 1)), d);
    if (e->op == "sin")
        return mul(op("cos", {a[0]}), d);
    if (e->op == "cos")
        return neg(mul(op("sin", {a[0]}), d));
    if (e->op == "exp")
        return mul(e, d);
    if (e->op == "log")
        return mul(power(a[0], -1), d);
    if (e->op == "sqrt")
        return mul(mul(num(.5), power(a[0], -.5)), d);
    if (e->op == "lap" || e->op == "dt" || e->op == "intg")
        return op(e->op, {d});
    throw Err("Cannot differentiate unexpanded " + e->op);
}

E subst(E e, const std::map<std::string, E>& d) {
    if (d.count(e->name))
        return d.at(e->name);
    if (e->a.empty())
        return e;
    std::vector<E> a;
    for (auto b : e->a)
        a.push_back(subst(b, d));
    int dim = e->op == "grad"  ? e->sh.back()
              : e->op == "div" ? e->a[0]->sh.back()
              : e->op == "eye" ? e->sh[0]
                               : 3;
    return op(e->op, a, dim);
}

std::string str(E e, bool cpp) {
    if (e->op == "num")
        return txt(e->n);
    if (e->op == "zero")
        return "0";
    if (!e->name.empty())
        return e->name == "pi" && cpp ? "3.141592653589793238462643383279502884" : e->name;
    auto a = e->a;
    std::string s;
    if (e->op == "add" || e->op == "mul")
        return "(" + str(a[0], cpp) + (e->op == "add" ? "+" : "*") + str(a[1], cpp) + ")";
    if (e->op == "pow")
        return cpp ? "std::pow(" + str(a[0], cpp) + "," + str(a[1], cpp) + ")"
                   : "(" + str(a[0]) + "**" + str(a[1]) + ")";
    if (e->op == "arr") {
        s = "[";
        for (auto v : a) {
            if (s != "[")
                s += ",";
            s += str(v, cpp);
        }
        return s + "]";
    }
    s = (cpp && std::set<std::string>{"sin", "cos", "log", "exp", "sqrt"}.count(e->op) ? "std::"
                                                                                       : "") +
        e->op + "(";
    for (int i = 0; i < a.size(); ++i) {
        if (i)
            s += ",";
        s += str(a[i], cpp);
    }
    return s + ")";
}

J ir(E e) {
    J a = J::array();
    for (auto x : e->a)
        a.push_back(ir(x));
    return J{{"op", e->op},
             {"shape", e->sh},
             {"value", e->op == "num" ? txt(e->n) : e->name},
             {"args", a}};
}

void Parser::ws() {
    while (at < s.size() && std::isspace((unsigned char)s[at]))
        ++at;
}

bool Parser::eat(char c) {
    ws();
    if (at < s.size() && s[at] == c) {
        ++at;
        return true;
    }
    return false;
}

E Parser::parse(const J& j) {
    req(j.is_string() || j.is_number(), "Expression must be text or number");
    return parse(j.is_string() ? j.get<std::string>() : txt(number(j, "expression")));
}

E Parser::parse(const std::string& t) {
    s = t;
    at = 0;
    depth = nodes = 0;
    req(!s.empty() && s.size() <= 4096, "Expression length must be 1..4096");
    E e = sum();
    ws();
    req(at == s.size(), "Unexpected syntax at offset " + std::to_string(at));
    return e;
}

E Parser::sum() {
    E e = prod();
    for (;;) {
        if (eat('+'))
            e = add(e, prod());
        else if (eat('-'))
            e = add(e, neg(prod()));
        else
            return e;
    }
}

E Parser::prod() {
    E e = unary();
    for (;;) {
        ws();
        if (at < s.size() && s[at] == '*' && (at + 1 == s.size() || s[at + 1] != '*')) {
            ++at;
            e = mul(e, unary());
        } else if (eat('/'))
            e = mul(e, power(unary(), -1));
        else
            return e;
    }
}

E Parser::unary() {
    req(++depth <= 96, "Expression nesting too deep");
    E e;
    if (eat('-'))
        e = neg(unary());
    else if (eat('+'))
        e = unary();
    else
        e = pow();
    --depth;
    return e;
}

E Parser::pow() {
    E e = atom();
    ws();
    if (at + 1 < s.size() && s.substr(at, 2) == "**") {
        at += 2;
        E exponent = unary();
        req(exponent->op == "num" && std::abs(exponent->n) <= 12,
            "Input exponent must be constant with abs <= 12");
        e = op("pow", {e, exponent}, dim);
    }
    return e;
}

E Parser::atom() {
    req(++nodes <= 2048, "Expression too complex");
    ws();
    req(at < s.size(), "Unexpected end of expression");
    if (eat('(')) {
        E e = sum();
        req(eat(')'), "Missing closing parenthesis");
        return e;
    }
    if (eat('[')) {
        std::vector<E> a;
        a.push_back(sum());
        while (eat(','))
            a.push_back(sum());
        req(eat(']'), "Missing closing bracket");
        return op("arr", a, dim);
    }
    if (std::isdigit((unsigned char)s[at]) || s[at] == '.') {
        char* end = nullptr;
        double n = std::strtod(s.c_str() + at, &end);
        req(end != s.c_str() + at, "Invalid number");
        at = end - s.c_str();
        return num(n);
    }
    req(std::isalpha((unsigned char)s[at]), "Invalid token");
    size_t lo = at++;
    while (at < s.size() && (std::isalnum((unsigned char)s[at]) || s[at] == '_'))
        ++at;
    std::string n = s.substr(lo, at - lo);
    if (eat('(')) {
        std::vector<E> a;
        if (!eat(')')) {
            a.push_back(sum());
            while (eat(','))
                a.push_back(sum());
            req(eat(')'), "Missing function parenthesis");
        }
        return op(n, a, dim);
    }
    if (n == "pi")
        return sym("pi");
    req(ns.count(n), "Unknown symbol: " + n);
    return ns.at(n);
}

double eval(E e, const std::map<std::string, double>& d) {
    double f = 0;
    if (e->op == "num")
        return e->n;
    if (e->op == "zero")
        return 0;
    if (!e->name.empty()) {
        if (e->name == "pi")
            return pi;
        req(d.count(e->name), "Unresolved symbol: " + e->name);
        return d.at(e->name);
    }
    double x = eval(e->a[0], d);
    if (e->op == "add")
        f = x + eval(e->a[1], d);
    else if (e->op == "mul")
        f = x * eval(e->a[1], d);
    else if (e->op == "pow")
        f = std::pow(x, e->a[1]->n);
    else if (e->op == "exp")
        f = std::exp(x);
    else if (e->op == "log")
        f = std::log(x);
    else if (e->op == "sin")
        f = std::sin(x);
    else if (e->op == "cos")
        f = std::cos(x);
    else if (e->op == "sqrt")
        f = std::sqrt(x);
    else
        throw Err("Not a local scalar expression: " + e->op);
    req(std::isfinite(f), "Expression left its finite real domain: " + e->op);
    return f;
}

E dx(E e, int k) {
    if (e->op == "num" || e->op == "zero")
        return num(0);
    if (!e->name.empty()) {
        if (e->op == "fld" || e->op == "test" || e->op == "coef")
            return sym(e->name + "$" + std::to_string(k), {}, e->op);
        return num(e->name == (k == 3 ? "t" : std::string(1, "xyz"[k])) ? 1 : 0);
    }
    auto a = e->a;
    E d = dx(a[0], k);
    if (e->op == "add")
        return add(d, dx(a[1], k));
    if (e->op == "mul")
        return add(mul(d, a[1]), mul(a[0], dx(a[1], k)));
    if (e->op == "pow")
        return mul(mul(a[1], power(a[0], a[1]->n - 1)), d);
    if (e->op == "sin")
        return mul(op("cos", {a[0]}), d);
    if (e->op == "cos")
        return neg(mul(op("sin", {a[0]}), d));
    if (e->op == "exp")
        return mul(e, d);
    if (e->op == "log")
        return mul(power(a[0], -1), d);
    if (e->op == "sqrt")
        return mul(mul(num(.5), power(a[0], -.5)), d);
    throw Err("Cannot spatially differentiate " + e->op);
}

std::vector<E> components(E e, int dim) {
    int n = count(e->sh);
    std::vector<E> out(n, num(0));
    if (e->op == "num")
        return {e};
    if (e->op == "zero")
        return out;
    if (!e->name.empty()) {
        for (int i = 0; i < n; ++i) {
            std::string s = e->name;
            int r = i;
            VI ix(e->sh.size());
            for (int k = int(ix.size()) - 1; k >= 0; --k) {
                ix[k] = r % e->sh[k];
                r /= e->sh[k];
            }
            for (int k : ix)
                s += "_" + std::to_string(k);
            out[i] = sym(s, {}, e->op);
        }
        return out;
    }
    if (e->op == "eye") {
        for (int i = 0; i < dim; ++i)
            out[i * dim + i] = num(1);
        return out;
    }
    std::vector<std::vector<E>> a;
    for (auto x : e->a)
        a.push_back(components(x, dim));
    auto s = e->op;
    if (s == "arr") {
        out.clear();
        for (auto v : a)
            out.insert(out.end(), v.begin(), v.end());
        return out;
    }
    for (int i = 0; i < n; ++i) {
        if (s == "add")
            out[i] = add(a[0][i], a[1][i]);
        else if (s == "mul")
            out[i] = mul(a[0][a[0].size() == 1 ? 0 : i], a[1][a[1].size() == 1 ? 0 : i]);
        else if (s == "grad")
            out[i] = dx(a[0][i / dim], i % dim);
        else if (s == "dt")
            out[i] = dx(a[0][i], 3);
        else if (s == "lap") {
            for (int k = 0; k < dim; ++k)
                out[i] = add(out[i], dx(dx(a[0][i], k), k));
        } else if (s == "div") {
            for (int k = 0; k < dim; ++k)
                out[i] = add(out[i], dx(a[0][i * dim + k], k));
        } else if (s == "adv") {
            for (int k = 0; k < dim; ++k)
                out[i] = add(out[i], mul(a[0][k], dx(a[1][i], k)));
        } else if (s == "inner") {
            for (int k = 0; k < a[0].size(); ++k)
                out[i] = add(out[i], mul(a[0][k], a[1][k]));
        } else if (s == "dot") {
            int kN = e->a[0]->sh.back(), nR = int(a[1].size()) / kN;
            for (int k = 0; k < kN; ++k)
                out[i] = add(out[i], mul(a[0][(i / nR) * kN + k], a[1][k * nR + i % nR]));
        } else if (s == "outer")
            out[i] = mul(a[0][i / a[1].size()], a[1][i % a[1].size()]);
        else if (s == "trans" || s == "sym") {
            int ix = (i % e->sh[1]) * e->sh[0] + i / e->sh[1];
            out[i] = s == "trans" ? a[0][ix] : mul(num(.5), add(a[0][i], a[0][ix]));
        } else if (s == "divg") {
            for (int k = 0; k < dim; ++k)
                out[i] = add(out[i], dx(mul(a[0][0], dx(a[1][i], k)), k));
        } else if (s == "intg")
            out[i] = op(s, {a[0][i]});
        else if (s == "pow")
            out[i] = power(a[0][0], a[1][0]->n);
        else
            out[i] = op(s, {a[0][i]});
    }
    return out;
}

std::pair<E, E> weak(E e, E v, int dim) {
    auto make = [&](const std::string& s, std::vector<E> a) {
        return op(s, a, dim);
    };
    E n = sym("n", {dim});
    if (e->op == "add") {
        auto a = weak(e->a[0], v, dim), b = weak(e->a[1], v, dim);
        return {add(a.first, b.first), add(a.second, b.second)};
    }
    if (e->op == "mul") {
        for (int k = 0; k < 2; ++k) {
            E a = e->a[k], b = e->a[1 - k];
            if (a->sh.empty() && has(b, {"lap", "div", "grad", "divg"})) {
                return weak(b, mul(a, v), dim);
            }
        }
    }
    if (e->op == "lap")
        return {neg(make("inner", {make("grad", {e->a[0]}), make("grad", {v})})),
                make("inner", {make("dot", {make("grad", {e->a[0]}), n}), v})};
    if (e->op == "divg")
        return weak(make("div", {mul(e->a[0], make("grad", {e->a[1]}))}), v, dim);
    if (e->op == "div" && e->a[0]->op != "fld" && e->a[0]->op != "coef" && e->a[0]->op != "sym") {
        E f = e->a[0];
        return {neg(make("inner", {f, make("grad", {v})})),
                make("inner", {make("dot", {f, n}), v})};
    }
    if (e->op == "grad") {
        E f = e->a[0];
        return {neg(make("inner", {f, make("div", {v})})), make("inner", {f, make("dot", {v, n})})};
    }
    return {make("inner", {e, v}), num(0)};
}
}
