#include "phz.hpp"
#define WIN32_LEAN_AND_MEAN
#include <windows.h>

namespace phz {
static J exactJson(const fs::path& p) {
    std::string s = read(p), out;
    bool quoted = false, escape = false;
    for (size_t i = 0; i < s.size();) {
        char c = s[i];
        if (quoted) {
            out += c;
            ++i;
            if (escape)
                escape = false;
            else if (c == '\\')
                escape = true;
            else if (c == '\"')
                quoted = false;
        } else if (c == '\"') {
            quoted = true;
            out += c;
            ++i;
        } else if ((c >= '0' && c <= '9') || c == '-') {
            size_t start = i++;
            while (i < s.size() && (std::isdigit((unsigned char)s[i]) || s[i] == '.' ||
                                    s[i] == 'e' || s[i] == 'E' || s[i] == '+' || s[i] == '-'))
                ++i;
            out += '\"';
            out += s.substr(start, i - start);
            out += '\"';
        } else {
            out += c;
            ++i;
        }
    }
    return J::parse(out);
}

static int sealedReplay(const fs::path& path, const fs::path& out) {
    req(verify(path)["ok"].get<bool>(), "Archive integrity verification failed");
    fs::path binary = fs::absolute(path / "Phyzo.exe");
    std::wstring cmd = L"\"" + binary.wstring() + L"\" replay-local \"" +
                       fs::absolute(path).wstring() + L"\" -o \"" + fs::absolute(out).wstring() +
                       L"\"";
    STARTUPINFOW si{};
    si.cb = sizeof(si);
    PROCESS_INFORMATION pi{};
    req(CreateProcessW(
            binary.c_str(), cmd.data(), nullptr, nullptr, FALSE, 0, nullptr, nullptr, &si, &pi) !=
            0,
        "Cannot launch saved native engine");
    WaitForSingleObject(pi.hProcess, INFINITE);
    DWORD code = 1;
    GetExitCodeProcess(pi.hProcess, &code);
    CloseHandle(pi.hThread);
    CloseHandle(pi.hProcess);
    return int(code);
}

int cli(int argc, char** argv) {
    if (argc < 2)
        return gui();
    std::string cmd = argv[1];
    std::vector<std::string> pos;
    std::map<std::string, std::string> options;
    std::set<std::string> flags;
    for (int i = 2; i < argc; ++i) {
        std::string s = argv[i];
        if (s == "--stationary" || s == "--fragment")
            flags.insert(s);
        else if (s == "-o" || s == "--archive" || s == "--method" || s == "--dt" ||
                 s == "--levels" || s == "--kind" || s == "--tableau" || s == "--max-order" ||
                 s == "--count" || s == "--seed") {
            req(i + 1 < argc, "Missing option value: " + s);
            options[s] = argv[++i];
        } else {
            req(s.empty() || s[0] != '-', "Unknown option: " + s);
            pos.push_back(s);
        }
    }
    auto opt = [&](std::string name, std::string fallback) {
        return options.count(name) ? options.at(name) : fallback;
    };
    auto model = [&]() {
        req(pos.size() == 1, "Supply one model name or path");
        return getModel(pos[0]);
    };
    fs::path output = opt("-o", dataRoot().string());
    if (cmd == "help" || cmd == "--help") {
        std::cout
            << "Phyzo 0.9.0 C++20\n\nPhyzo.exe gui\nPhyzo.exe examples\nPhyzo.exe run MODEL [-o RUNS]\nPhyzo.exe stat MODEL [-o SUMMARY] [--archive RUNS]\nPhyzo.exe gen MODEL [-o RHS.cpp] [--fragment]\nPhyzo.exe weak MODEL [-o FORM.json]\nPhyzo.exe study MODEL --kind mesh|quad|time --levels 3 [--stationary] [-o RUNS]\nPhyzo.exe scheme --method FE|BE|CN|RK2|SSPRK3|RK4|RADAU2 [--tableau FILE] [--max-order 6]\nPhyzo.exe analyze MODEL --method FE --dt 0.01\nPhyzo.exe verify RUN\nPhyzo.exe replay RUN [-o RUNS]\nPhyzo.exe show RUN\nPhyzo.exe migrate RUN [-o RUNS]\nPhyzo.exe selftest [--count 1000] [--seed 20261115] [-o REPORT.json]\nPhyzo.exe inspect MODEL [-o REPORT.json]\n\nModels may be embedded names or JSON paths. All computation is native.\nDefault result directory: "
            << dataRoot() << "\n";
        return 0;
    }
    if (cmd == "examples") {
        for (auto [s, t] : examples())
            std::cout << s << "\n";
        return 0;
    }
    if (cmd == "gui")
        return gui(pos.empty() ? "nseq2" : pos[0]);
    if (cmd == "selftest")
        return selftest(options.count("-o") ? output : dataRoot() / "selftest.json",
                        std::stoi(opt("--count", "1000")),
                        std::stoul(opt("--seed", "20261115")));
    if (cmd == "scheme") {
        J tab = options.count("--tableau") ? exactJson(options["--tableau"]) : J();
        J out = scheme(tab.is_null() ? opt("--method", "FE") : "custom",
                       tab,
                       std::stoi(opt("--max-order", "6")));
        auto p = options.count("-o") ? output : dataRoot() / "scheme.json";
        save(p, out);
        std::cout << out.dump(2) << "\n";
        return 0;
    }
    if (cmd == "verify" || cmd == "show" || cmd == "replay" || cmd == "replay-local" ||
        cmd == "migrate") {
        req(pos.size() == 1, "Supply one archive path");
        fs::path p = fs::absolute(pos[0]);
        if (cmd == "verify") {
            J j = verify(p);
            std::cout << j.dump(2) << "\n";
            return j["ok"].get<bool>() ? 0 : 2;
        }
        if (cmd == "show") {
            req(verify(p)["ok"].get<bool>(), "Archive changed");
            std::cout << read(p / "report.md");
            return 0;
        }
        if (cmd == "migrate") {
            Result r = restore(p);
            r.info["imported_from"] = p.string();
            r.info["original_manifest_sha256"] = sha256(p / "manifest.json");
            r.info["migration_scope"] =
                "Preserved original numeric arrays and checks; replay recomputes with native C++ engine";
            auto dst = r.archive(output);
            std::cout << dst.string() << "\n";
            return 0;
        }
        if (cmd == "replay" && load(p / "manifest.json").value("schema", 0) == 2)
            return sealedReplay(p, output);
        Result r = replay(p);
        auto archive = r.archive(output);
        std::cout << archive.string() << "\n" << r.checks.dump(2) << "\n";
        return r.passed() ? 0 : 2;
    }
    if (cmd == "weak") {
        PDE p(model());
        auto target = options.count("-o") ? output : dataRoot() / "weak.json";
        J j = p.weakForm();
        save(target, j);
        std::cout << j.dump(2) << "\n";
        return 0;
    }
    if (cmd == "gen") {
        System s(model());
        auto target = options.count("-o") ? output : dataRoot() / "rhs.cpp";
        req(target.extension() == ".cpp" || target.extension() == ".hpp",
            "Generated equations must use .cpp or .hpp");
        if (!pos.empty() && fs::exists(pos[0]))
            req(fs::weakly_canonical(target) != fs::weakly_canonical(pos[0]),
                "Output must not overwrite model");
        write(target, s.exportCpp(flags.count("--fragment")));
        std::cout << fs::absolute(target).string() << "\n";
        return 0;
    }
    if (cmd == "analyze") {
        Model m(model());
        std::string method = opt("--method", "FE");
        J out = scheme(method);
        out["heat"] = heat(m, method, options.count("--dt") ? std::stod(options["--dt"]) : -1);
        save(options.count("-o") ? output : dataRoot() / "analysis.json", out);
        std::cout << out.dump(2) << "\n";
        return out["heat"].value("dt_status", "") == "fail" ? 2 : 0;
    }
    if (cmd == "study") {
        J out = study(model(),
                      output,
                      std::stoi(opt("--levels", "3")),
                      opt("--kind", "mesh"),
                      flags.count("--stationary"));
        std::cout << out.dump(2) << "\n";
        for (auto r : out["levels"])
            for (auto c : r["checks"])
                if (c.value("status", "") == "fail")
                    return 2;
        return 0;
    }
    if (cmd == "inspect") {
        auto s = std::make_shared<System>(model());
        J out{{"model", s->mdl->input},
              {"u0", arr(s->u0)},
              {"active", s->active},
              {"mass", arr(M(s->mass))},
              {"stiff", arr(M(s->stiff))},
              {"quadrature_points", arr(s->x)},
              {"quadrature_weights", arr(s->w)}};
        if (s->mesh)
            out["mesh"] = s->mesh->json();
        if (s->mixed) {
            auto [r, j] = s->assemble(0, s->u0);
            out["weak_residual"] = arr(r);
            out["weak_jacobian"] = arr(M(j));
        } else {
            out["rhs"] = arr(s->rhs(0, s->u0));
            out["jacobian"] = arr(s->jacobian(0, s->u0));
            if (!s->terms.empty()) {
                auto [r, j] = s->assemble(0, s->u0);
                out["weak_residual"] = arr(r);
                out["weak_jacobian"] = arr(M(j));
            }
        }
        if (options.count("-o"))
            save(output, out);
        else
            std::cout << out.dump(2) << "\n";
        return 0;
    }
    req(cmd == "run" || cmd == "stat", "Unknown command: " + cmd);
    Result r = execute(model(), cmd == "stat");
    auto archive =
        r.archive(cmd == "stat" ? fs::path(opt("--archive", dataRoot().string())) : output);
    if (cmd == "stat") {
        J summary = r.info;
        summary["run"] = archive.string();
        save(options.count("-o") ? output : dataRoot() / "stat.json", summary);
    }
    std::cout << archive.string() << "\n" << r.checks.dump(2) << "\n";
    return r.passed() ? 0 : 2;
}
}

int main(int argc, char** argv) {
    SetConsoleOutputCP(CP_UTF8);
    try {
        return phz::cli(argc, argv);
    } catch (const phz::Stop& e) {
        std::cerr << "Cancelled: " << e.what() << "\n";
        return 3;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
