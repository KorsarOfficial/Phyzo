#include "phz.hpp"

namespace phz {
using Modes = std::map<VI, E>;

static Modes product(const Modes& a, const Modes& b) {
    req(a.size() * b.size() <= 300000, "Symbolic convolution too large");
    Modes out;
    for (auto [ka, ea] : a)
        for (auto [kb, eb] : b) {
            VI k(ka.size());
            std::function<void(int, E)> rec = [&](int i, E v) {
                if (i == k.size()) {
                    out[k] = out.count(k) ? add(out[k], v) : v;
                    return;
                }
                if (ka[i] == 0 || kb[i] == 0) {
                    k[i] = ka[i] + kb[i];
                    rec(i + 1, v);
                } else {
                    k[i] = ka[i] + kb[i];
                    rec(i + 1, mul(num(.5), v));
                    k[i] = std::abs(ka[i] - kb[i]);
                    rec(i + 1, mul(num(.5), v));
                }
            };
            rec(0, mul(ea, eb));
        }
    return out;
}

std::string System::exportCpp(bool fragment) {
    req(mdl->type == "sp" && exactSp,
        "Explicit export requires autonomous polynomial cosine model, <=128 coefficients");
    if (spRhs.empty()) {
        VI z(mdl->dim, 0);
        std::map<std::string, Modes> fields;
        for (int f = 0; f < nf; ++f)
            for (int k = 0; k < n; ++k)
                fields[compNames[f]][modes[k]] = sym("u[" + std::to_string(f * n + k) + "]");
        std::map<E, Modes, std::owner_less<E>> memo;
        std::function<Modes(E)> convert = [&](E e) -> Modes {
            if (memo.count(e))
                return memo[e];
            Modes out;
            if (fields.count(e->name))
                out = fields.at(e->name);
            else if (!has(e, {"fld", "lap", "divg", "intg"}) && !depends(e, "x") &&
                     !depends(e, "y") && !depends(e, "z") && !depends(e, "t"))
                out[z] = e;
            else if (e->op == "add") {
                out = convert(e->a[0]);
                for (auto [k, v] : convert(e->a[1]))
                    out[k] = out.count(k) ? add(out[k], v) : v;
            } else if (e->op == "mul")
                out = product(convert(e->a[0]), convert(e->a[1]));
            else if (e->op == "pow") {
                req(e->a[1]->n >= 0 && e->a[1]->n <= 4 && e->a[1]->n == std::floor(e->a[1]->n),
                    "Polynomial exponent must be 0..4");
                out[z] = num(1);
                for (int i = 0; i < e->a[1]->n; ++i)
                    out = product(out, convert(e->a[0]));
            } else if (e->op == "lap") {
                out = convert(e->a[0]);
                for (auto& [k, v] : out) {
                    E lam = num(0);
                    for (int i = 0; i < mdl->dim; ++i)
                        lam = add(
                            lam,
                            mul(num(k[i] * k[i]), power(sym("Lx" + std::to_string(i + 1)), -2)));
                    v = neg(mul(mul(power(sym("pi"), 2), lam), v));
                }
            } else if (e->op == "divg") {
                E a = e->a[0], b = e->a[1];
                out = convert(mul(num(.5),
                                  add(add(op("lap", {mul(a, b)}), mul(a, op("lap", {b}))),
                                      neg(mul(b, op("lap", {a}))))));
            } else
                throw Err("Expression does not support exact cosine export");
            memo[e] = out;
            return out;
        };
        for (auto s : compNames) {
            auto out = convert(mdl->eq.at(s));
            for (auto k : modes)
                spRhs.push_back(out.count(k) ? out.at(k) : num(0));
        }
    }
    std::ostringstream out, body;
    std::map<std::string, E> aliases;
    int parameter = 0;
    for (auto [name, value] : mdl->par) {
        std::string key = "fP" + std::to_string(parameter++);
        aliases[name] = sym(key);
        body << "const double " << key << "=dP.contains(\"" << name << "\")?dP.at(\"" << name
             << "\"):" << txt(value) << ";\nif(!std::isfinite(" << key
             << "))throw std::invalid_argument(\"Non-finite parameter\");\n";
    }
    for (int k = 0; k < mdl->dim; ++k)
        body << "if(" << aliases.at("Lx" + std::to_string(k + 1))->name
             << "<=0)throw std::invalid_argument(\"Invalid domain length\");\n";
    body
        << "for(double x:u)if(!std::isfinite(x))throw std::invalid_argument(\"Non-finite state\");\nstd::array<double,"
        << size << "> du_dt{};\n";
    for (int k = 0; k < size; ++k)
        body << "du_dt[" << k << "]=" << str(subst(spRhs[k], aliases), true) << ";\n";
    body
        << "for(double x:du_dt)if(!std::isfinite(x))throw std::domain_error(\"Non-finite RHS\");\n";
    if (fragment)
        return body.str();
    out << "#include <array>\n#include <cmath>\n#include <map>\n#include <string>\n#include <stdexcept>\n\n";
    out << "inline const std::array<const char*," << size << "> ORDER={";
    for (int i = 0; i < size; ++i) {
        if (i)
            out << ",";
        out << "\"" << compNames[i / n];
        for (int k : modes[i % n])
            out << "_" << k;
        out << "\"";
    }
    out << "};\ninline const std::map<std::string,double> DEFAULT_PARAMETERS={";
    int i = 0;
    for (auto [s, v] : mdl->par) {
        if (i++)
            out << ",";
        out << "{\"" << s << "\"," << txt(v) << "}";
    }
    out << "};\n\ninline std::array<double," << size << "> aRhs(double t,const std::array<double,"
        << size << ">& u,const std::map<std::string,double>& dP=DEFAULT_PARAMETERS){\n"
        << body.str() << "return du_dt;\n}\n";
    return out.str();
}
}