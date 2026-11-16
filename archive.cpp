#include "phz.hpp"
#define WIN32_LEAN_AND_MEAN
#include <windows.h>

namespace phz {
static J sparseJson(const SM& a) {
    J r = J::array(), c = J::array(), v = J::array();
    for (int k = 0; k < a.outerSize(); ++k)
        for (SM::InnerIterator it(a, k); it; ++it) {
            r.push_back(it.row());
            c.push_back(it.col());
            v.push_back(it.value());
        }
    return J{{"rows", a.rows()}, {"cols", a.cols()}, {"row", r}, {"col", c}, {"value", v}};
}

static std::string htmlEscape(std::string s) {
    std::string o;
    for (char c : s) {
        if (c == '&')
            o += "&amp;";
        else if (c == '<')
            o += "&lt;";
        else if (c == '>')
            o += "&gt;";
        else if (c == '\"')
            o += "&quot;";
        else
            o += c;
    }
    return o;
}

static std::string picture(const Result& r) {
    auto& s = *r.sys;
    M pts;
    V vals;
    if (s.mesh) {
        pts = s.mesh->p;
        vals = r.u.back().segment(s.offset[0], pts.rows());
    } else {
        int n = 101;
        pts = M::Zero(n, s.mdl->dim);
        for (int i = 0; i < n; ++i)
            pts(i, 0) = s.mdl->par.at("Lx1") * i / (n - 1);
        vals = s.values(r.u.back(), pts, 0);
    }
    double lo = vals.minCoeff(), hi = vals.maxCoeff(), xmin = pts.col(0).minCoeff(),
           xmax = pts.col(0).maxCoeff(), ymin = pts.cols() > 1 ? pts.col(1).minCoeff() : lo,
           ymax = pts.cols() > 1 ? pts.col(1).maxCoeff() : hi;
    auto xp = [&](double x) {
        return 35 + 730 * (x - xmin) / std::max(1e-15, xmax - xmin);
    };
    auto yp = [&](double y) {
        return 550 - 490 * (y - ymin) / std::max(1e-15, ymax - ymin);
    };
    std::ostringstream o;
    o << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"800\" height=\"600\" viewBox=\"0 0 800 600\"><rect width=\"800\" height=\"600\" fill=\"#111b2a\"/><text x=\"30\" y=\"28\" fill=\"white\" font-family=\"sans-serif\">"
      << htmlEscape(s.mdl->name) << " | " << htmlEscape(s.mdl->fields[0]) << " | t=" << r.t.back()
      << "</text>";
    if (pts.cols() == 1) {
        o << "<polyline fill=\"none\" stroke=\"#5ccef0\" stroke-width=\"2\" points=\"";
        for (int i = 0; i < pts.rows(); ++i)
            o << xp(pts(i, 0)) << "," << yp(vals[i]) << " ";
        o << "\"/>";
    } else
        for (int i = 0; i < pts.rows(); ++i) {
            double a = (vals[i] - lo) / std::max(1e-15, hi - lo);
            int red = int(255 * a), blue = int(255 * (1 - a)),
                green = int(180 * (1 - std::abs(2 * a - 1)));
            o << "<circle cx=\"" << xp(pts(i, 0)) << "\" cy=\"" << yp(pts(i, 1))
              << "\" r=\"4\" fill=\"rgb(" << red << "," << green << "," << blue << ")\"/>";
        }
    o << "<text x=\"30\" y=\"584\" fill=\"white\" font-family=\"sans-serif\">min=" << lo
      << " max=" << hi << (pts.cols() > 2 ? " | projection onto x-y plane" : "") << "</text></svg>";
    return o.str();
}

fs::path Result::archive(const fs::path& root) const {
    fs::create_directories(root);
    SYSTEMTIME time;
    GetSystemTime(&time);
    std::ostringstream stamp;
    stamp << std::setfill('0') << std::setw(4) << time.wYear << std::setw(2) << time.wMonth
          << std::setw(2) << time.wDay << "T" << std::setw(2) << time.wHour << std::setw(2)
          << time.wMinute << std::setw(2) << time.wSecond << "Z_" << std::hex << GetTickCount64()
          << "_" << GetCurrentProcessId();
    std::string id = stamp.str();
    fs::path path = fs::absolute(root) / id, temp = path.parent_path() / ("." + id);
    req(!fs::exists(path) && !fs::exists(temp), "Archive collision");
    fs::create_directory(temp);
    auto& s = *sys;
    J model = s.mdl->input;
    if (s.mesh)
        model["dom"] = s.mesh->json();
    save(temp / "model.json", model);
    J state{{"t", t}, {"u", J::array()}};
    for (auto a : u)
        state["u"].push_back(arr(a));
    save(temp / "state.json", state);
    save(temp / "diag.json", diag);
    save(temp / "checks.json", checks);
    save(temp / "mass.json", sparseJson(s.mass));
    save(temp / "stiff.json", sparseJson(s.stiff));
    J disc{{"method", s.mdl->type},
           {"dimension", s.mdl->dim},
           {"quadrature_order", s.mdl->disc.value("q", s.mdl->type == "sp" ? 16 : 4)},
           {"mass_matrix", "consistent"},
           {"nonlinear_evaluation",
            s.exactSp ? "Exact cosine convolution; projection only at final retained modes"
                      : "Reconstruct in quadrature points, evaluate nonlinearity, integrate"},
           {"basis_order",
            s.mixed               ? J::array({2, 1})
            : s.mdl->type == "fe" ? J(1)
                                  : s.mdl->disc.value("n", J())},
           {"boundary_conditions", s.mdl->bc},
           {"parent_run", info.value("parent_run", J())}};
    if (s.mesh) {
        save(temp / "mesh.json", s.mesh->json());
        disc["mesh_sha256"] = sha256(temp / "mesh.json");
        disc["embedding_dimension"] = s.mesh->geo;
    }
    disc["assumptions"] = s.mdl->input.value("asm", J::array());
    disc["units"] = s.mdl->input.value("units", J::object());
    save(temp / "disc.json", disc);
    save(
        temp / "boundary.json",
        J{{"patches", s.mdl->bc},
          {"flux_convention", "Positive outward total flux q=-D*grad(u).n; Robin q=coef*(u-val)"}});
    J order = J::array();
    for (int f = 0; f < s.compNames.size(); ++f) {
        int n = s.mixed ? s.spaces[f]->n : s.n;
        for (int k = 0; k < n; ++k)
            order.push_back(J{{"field", s.compNames[f]},
                              {"index", s.offset[f] + k},
                              {"basis", s.mesh ? VI{k} : s.modes[k]}});
    }
    save(temp / "order.json", order);
    if (s.mixed) {
        J spaces = J::array();
        for (auto sp : s.spaces)
            spaces.push_back(J{{"degree", sp->pdeg},
                               {"p", arr(sp->p)},
                               {"c", arr(sp->c)},
                               {"mass", sparseJson(sp->mass)},
                               {"stiff", sparseJson(sp->stiff)}});
        save(temp / "spaces.json", spaces);
        save(temp / "div.json", sparseJson(s.div));
    }
    if (!scheme.empty())
        save(temp / "scheme.json", scheme);
    if (info.contains("quadrature_control")) {
        save(temp / "quad.json", info["quadrature_control"]);
        save(temp / "quad_model.json", info["quadrature_request"]);
    }
    J env{{"engine", "Phyzo C++"},
          {"version", "0.9.0-cpp"},
          {"cplusplus", _MSVC_LANG},
          {"compiler", "MSVC " + std::to_string(_MSC_FULL_VER)},
          {"eigen", "3.4.0"},
          {"sundials", "7.4.0"},
          {"json", "3.12.0"},
          {"boost", "1.86.0"},
          {"architecture", "x64"},
          {"source_language", "C++20"},
          {"python_runtime", false},
          {"floating_point", "IEEE754 double; no fast-math"}};
    save(temp / "env.json", env);
    fs::path source = fs::exists(exePath().parent_path() / "src/phz.hpp")
                          ? exePath().parent_path() / "src"
                      : fs::exists("C:/Phyzo/cpp/phz.hpp") ? fs::path("C:/Phyzo/cpp")
                                                           : fs::path("C:/Phyzo");
    fs::create_directory(temp / "src");
    for (auto file : fs::directory_iterator(source))
        if (file.is_regular_file() &&
            (file.path().extension() == ".cpp" || file.path().extension() == ".hpp"))
            fs::copy_file(file.path(), temp / "src" / file.path().filename());
    fs::copy_file(exePath(), temp / "Phyzo.exe");
    for (auto name :
         {"vcruntime140.dll", "vcruntime140_1.dll", "msvcp140.dll", "msvcp140_atomic_wait.dll"}) {
        auto dll = exePath().parent_path() / name;
        if (fs::exists(dll))
            fs::copy_file(dll, temp / name);
    }
    fs::path notices = fs::is_directory(exePath().parent_path() / "licenses")
                           ? exePath().parent_path() / "licenses"
                           : exePath().parent_path().parent_path() / "licenses";
    if (fs::is_directory(notices))
        fs::copy(notices, temp / "licenses", fs::copy_options::recursive);
    write(temp / "field.svg", picture(*this));
    std::ostringstream report;
    report << "# Phyzo C++ calculation\n\nModel: " << s.mdl->name
           << "\n\nKind: " << (stationary ? "stationary" : "transient")
           << "; coefficients: " << s.size << "; states: " << u.size()
           << ".\n\n| Check | Status |\n|---|---|\n";
    for (auto ch : checks)
        report << "| " << ch.value("name", "") << " | " << ch.value("status", "") << " |\n";
    report
        << "\nMetadata: model.json, disc.json, boundary.json, env.json, manifest.json.\n\nReplay uses the sealed native Phyzo.exe in this archive.\n\n![Field](field.svg)\n\nNumerical checks cover this model and sampled states; physical adequacy and formal proof are not inferred.\n";
    write(temp / "report.md", report.str());
    J hashes = J::object();
    for (auto f : fs::recursive_directory_iterator(temp))
        if (f.is_regular_file())
            hashes[fs::relative(f.path(), temp).generic_string()] = sha256(f.path());
    J mf{
        {"schema", 2},
        {"id", id},
        {"info", info},
        {"files", hashes},
        {"integrity_scope",
         "SHA256 detects changes relative to manifest; this is not a signed authenticity certificate"}};
    save(temp / "manifest.json", mf);
    for (int k = 0; k < 5; ++k) {
        try {
            fs::rename(temp, path);
            return path;
        } catch (const fs::filesystem_error&) {
            if (k == 4)
                throw;
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
        }
    }
    throw Err("Archive finalization failed");
}

J verify(const fs::path& p) {
    J mf = load(p / "manifest.json");
    req(mf.is_object() && mf.contains("files") && mf["files"].is_object(), "Invalid manifest");
    int schema = mf.value("schema", 0);
    req(schema == 1 || schema == 2, "Unsupported archive schema");
    std::vector<std::string> mandatory = schema == 2
                                             ? std::vector<std::string>{"model.json",
                                                                        "state.json",
                                                                        "diag.json",
                                                                        "checks.json",
                                                                        "disc.json",
                                                                        "env.json",
                                                                        "order.json",
                                                                        "Phyzo.exe",
                                                                        "src/phz.hpp"}
                                             : std::vector<std::string>{"model.json",
                                                                        "state.npz",
                                                                        "diag.npz",
                                                                        "checks.json",
                                                                        "disc.json",
                                                                        "env.json",
                                                                        "order.json",
                                                                        "mass.npz",
                                                                        "stiff.npz",
                                                                        "report.md",
                                                                        "src/phz/__init__.py"};
    for (auto s : mandatory)
        req(mf["files"].contains(s), "Manifest lacks mandatory entry: " + std::string(s));
    J bad = J::array();
    for (auto it = mf["files"].begin(); it != mf["files"].end(); ++it) {
        fs::path rel(it.key());
        req(!rel.empty() && !rel.is_absolute() && !rel.has_root_path(), "Unsafe manifest path");
        for (auto part : rel)
            req(part != ".." && part != ".", "Unsafe manifest traversal");
        req(it.key().find(':') == std::string::npos, "Unsafe manifest stream");
        fs::path file = p / rel;
        bool ok = fs::is_regular_file(file) && !fs::is_symlink(file);
        if (ok) {
            auto root = fs::weakly_canonical(p), resolved = fs::weakly_canonical(file);
            auto relative = resolved.lexically_relative(root);
            for (auto part : relative)
                if (part == "..")
                    ok = false;
        }
        if (!ok || sha256(file) != it.value().get<std::string>())
            bad.push_back(it.key());
    }
    return J{{"ok", bad.empty()}, {"changed_or_missing", bad}, {"files", mf["files"].size()}};
}

Result restore(const fs::path& p) {
    req(verify(p)["ok"].get<bool>(), "Archive integrity verification failed");
    Result r;
    r.sys = std::make_shared<System>(load(p / "model.json"));
    bool legacy = load(p / "manifest.json").value("schema", 0) == 1;
    J state = legacy ? npz(p / "state.npz") : load(p / "state.json");
    r.t = state.at("t").get<std::vector<double>>();
    for (auto a : state.at("u")) {
        V v = vec(a);
        req(v.size() == r.sys->size, "Archive state shape mismatch");
        r.u.push_back(v);
    }
    req(r.t.size() == r.u.size() && !r.t.empty(), "Archive time/state mismatch");
    r.info = load(p / "manifest.json").at("info");
    r.stationary = r.info.value("kind", "") == "stationary";
    r.diag = legacy ? npz(p / "diag.npz") : load(p / "diag.json");
    r.checks = load(p / "checks.json");
    if (fs::exists(p / "scheme.json"))
        r.scheme = load(p / "scheme.json");
    return r;
}

Result replay(const fs::path& p) {
    req(verify(p)["ok"].get<bool>(), "Archive integrity verification failed");
    J mf = load(p / "manifest.json");
    bool st = mf.at("info").value("kind", "") == "stationary";
    return execute(
        load(p / (mf.at("info").contains("quadrature_control") ? "quad_model.json" : "model.json")),
        st);
}
}
