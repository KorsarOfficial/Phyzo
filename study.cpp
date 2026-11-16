#include "phz.hpp"

namespace phz {
J study(J j, const fs::path& root, int levels, const std::string& kind, bool st) {
    req(levels >= 2 && levels <= 5, "Refinement levels must be 2..5");
    req(kind == "mesh" || kind == "time" || kind == "quad", "Unknown refinement kind");
    Model m(j);
    req(!m.disc.contains("aq"), "Independent study requires fixed quadrature");
    req(!(st && kind == "time"), "Stationary study has no time step");
    if (kind == "time")
        req(std::set<std::string>{"FE", "BE", "RK2", "RK4", "SSPRK3"}.count(
                m.run.value("method", "BDF")),
            "Time refinement requires a fixed-step integrator");
    if (kind == "mesh" && j.contains("ini"))
        for (auto it = j["ini"].begin(); it != j["ini"].end(); ++it)
            req(!it.value().is_array(), "Mesh refinement requires initial expressions");
    fs::create_directories(root);
    std::unique_ptr<Result> prev;
    J rows = J::array();
    fs::path oldPath;
    double baseH = 0;
    if (kind == "time") {
        double t0 = m.run.value("t", J::array({0, 1}))[0],
               t1 = m.run.value("t", J::array({0, 1}))[1];
        int n = m.run.value("n", 51);
        double interval = (t1 - t0) / (n - 1), h = m.run.value("max_step", interval);
        baseH = interval / std::max(1., std::ceil(interval / h - 1e-12));
    }
    if (kind == "mesh" && m.type == "sp") {
        VI ns = m.disc.value("n", VI(m.dim, 4));
        int maxn = *std::max_element(ns.begin(), ns.end()) + levels - 1;
        int q = m.disc.value("q", std::max(16, 3 * maxn));
        req(q >= 2 * maxn, "Quadrature insufficient for all spectral levels");
        j["disc"]["q"] = q;
    }
    VI ns0 = m.disc.value("n", VI(m.dim, 4));
    int q0 = m.disc.value(
        "q", m.type == "sp" ? std::max(16, 3 * *std::max_element(ns0.begin(), ns0.end())) : 4);
    std::vector<V> previousDiff;
    for (int lv = 0; lv < levels; ++lv) {
        checkpoint();
        if (kind == "time")
            j["run"]["max_step"] = baseH / std::pow(2., lv);
        else if (kind == "quad")
            j["disc"]["q"] = q0 + lv * (m.type == "sp" ? 4 : 2);
        else if (lv) {
            if (m.type == "sp") {
                VI ns = ns0;
                for (int& n : ns)
                    n += lv;
                j["disc"]["n"] = ns;
            } else
                j["dom"] = prev->sys->mesh->refine().first.json();
        }
        Result now = execute(j, st);
        if (!oldPath.empty())
            now.info["parent_run"] = oldPath.string();
        J row{{"level", lv},
              {"coefficients", now.sys->size},
              {"quadrature", now.sys->mdl->disc.value("q", q0)},
              {"checks", now.checks},
              {"solve", st ? "stationary" : "transient"}};
        if (kind == "time")
            row["step"] = j["run"]["max_step"];
        if (now.diag.contains("l2")) {
            row["exact_l2_final"] = now.diag["l2"].back();
        }
        if (prev) {
            std::vector<V> differences;
            V maxdiff = V::Zero(now.sys->mixed ? 2 : now.sys->nf);
            J samples = J::array();
            for (int it = 0; it < now.u.size(); ++it) {
                V norm;
                if (kind != "mesh")
                    norm = now.sys->norm(now.u[it] - prev->u[it]);
                else {
                    auto& s = *now.sys;
                    int nc = int(s.compNames.size());
                    V sq = V::Zero(s.mixed ? 2 : nc);
                    for (int f = 0; f < nc; ++f) {
                        V fine = s.mixed ? s.spaces[f]->b *
                                               now.u[it].segment(s.offset[f], s.spaces[f]->n)
                                         : s.b * now.u[it].segment(f * s.n, s.n);
                        V coarse = prev->sys->values(prev->u[it], s.x, f);
                        double error = s.w.dot((fine - coarse).array().square().matrix());
                        sq[s.mixed ? (f < s.mdl->dim ? 0 : 1) : f] += std::max(0., error);
                    }
                    norm = sq.array().sqrt();
                }
                samples.push_back(arr(norm));
                maxdiff = maxdiff.cwiseMax(norm);
                differences.push_back(norm);
            }
            row["l2_at_t"] = samples;
            row["difference_final"] = arr(differences.back());
            row["difference_max"] = arr(maxdiff);
            if (!previousDiff.empty()) {
                J order = J::array();
                for (int k = 0; k < maxdiff.size(); ++k) {
                    double a = previousDiff.back()[k], b = differences.back()[k];
                    order.push_back((kind == "time" || (kind == "mesh" && m.type != "sp")) &&
                                            a > 0 && b > 0
                                        ? J(std::log2(a / b))
                                        : J());
                }
                row["observed_order_final"] = order;
            }
            previousDiff = differences;
        }
        oldPath = now.archive(root);
        row["run"] = oldPath.string();
        rows.push_back(row);
        if (st && lv == 0 && j.value("stat", J::object()).contains("mass")) {
            J masses = j["stat"]["mass"];
            for (auto it = masses.begin(); it != masses.end(); ++it)
                if (it.value().is_string() && it.value() == "initial") {
                    auto& s = *now.sys;
                    int f = int(std::find(s.compNames.begin(), s.compNames.end(), it.key()) -
                                s.compNames.begin());
                    j["stat"]["mass"][it.key()] = s.massv.dot(s.u0.segment(f * s.n, s.n));
                }
        }
        prev = std::make_unique<Result>(std::move(now));
    }
    J out{{"kind", kind},
          {"solve", st ? "stationary" : "transient"},
          {"levels", rows},
          {"scope",
           "Change one discretization parameter; observed order is numerical evidence, not proof"}};
    save(root / "study.json", out);
    return out;
}
}
