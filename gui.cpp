#include "phz.hpp"
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <commctrl.h>
#include <commdlg.h>
#include <shellapi.h>
#include <mutex>

namespace phz {
static std::wstring wide(const std::string& s) {
    if (s.empty())
        return {};
    int n = MultiByteToWideChar(CP_UTF8, MB_ERR_INVALID_CHARS, s.data(), int(s.size()), nullptr, 0);
    req(n > 0, "Invalid UTF-8 text");
    std::wstring w(n, 0);
    MultiByteToWideChar(CP_UTF8, 0, s.data(), int(s.size()), w.data(), n);
    return w;
}

static std::string utf8(const std::wstring& s) {
    if (s.empty())
        return {};
    int n = WideCharToMultiByte(CP_UTF8, 0, s.data(), int(s.size()), nullptr, 0, nullptr, nullptr);
    std::string t(n, 0);
    WideCharToMultiByte(CP_UTF8, 0, s.data(), int(s.size()), t.data(), n, nullptr, nullptr);
    return t;
}

static std::string windowText(HWND w) {
    int n = GetWindowTextLengthW(w);
    std::wstring s(n + 1, 0);
    GetWindowTextW(w, s.data(), n + 1);
    s.resize(n);
    return utf8(s);
}

struct Gui {
    HWND window = nullptr, editor = nullptr, log = nullptr, plot = nullptr, models = nullptr,
         fields = nullptr, slider = nullptr, status = nullptr, kind = nullptr, levels = nullptr;
    std::vector<HWND> actions;
    HFONT mono = nullptr, font = nullptr;
    std::thread worker;
    std::mutex mutex;
    std::unique_ptr<Result> result, pending;
    std::string message;
    fs::path archive;
    bool busy = false, closing = false, drag = false;
    int sample = 0, field = 0, lastx = 0, lasty = 0;
    double yaw = .7, pitch = .6, zoom = 1;
    std::vector<std::string> names;
};

static Gui* app = nullptr;

static COLORREF color(double a) {
    a = std::clamp(a, 0., 1.);
    return RGB(int(30 + 220 * a), int(70 + 150 * (1 - std::abs(2 * a - 1))), int(245 - 210 * a));
}

static void label(HWND w, const std::string& s) {
    wchar_t cls[32]{};
    GetClassNameW(w, cls, 32);
    std::string text = s;
    if (_wcsicmp(cls, L"EDIT") == 0) {
        text.clear();
        for (size_t i = 0; i < s.size(); ++i) {
            if (s[i] == '\n' && (i == 0 || s[i - 1] != '\r'))
                text += '\r';
            text += s[i];
        }
    }
    SetWindowTextW(w, wide(text).c_str());
}

static void layout(Gui& g) {
    if (!g.window || !g.editor)
        return;
    RECT r;
    GetClientRect(g.window, &r);
    int w = r.right, h = r.bottom, left = std::max(300, w * 42 / 100), top = 114;
    MoveWindow(g.models, 12, 12, 230, 250, TRUE);
    int x = 252;
    for (int i = 0; i < g.actions.size(); ++i) {
        if (i < 5) {
            MoveWindow(g.actions[i], x, 12, 108, 28, TRUE);
            x += 114;
        } else if (i < 8)
            MoveWindow(g.actions[i], 12 + (i - 5) * 114, 46, 108, 27, TRUE);
        else
            MoveWindow(g.actions[i], 12 + (i - 8) * 114, 80, 108, 27, TRUE);
    }
    MoveWindow(g.kind, 252, 80, 110, 200, TRUE);
    MoveWindow(g.levels, 368, 80, 55, 25, TRUE);
    MoveWindow(g.fields, left + 12, 46, 140, 200, TRUE);
    MoveWindow(g.slider, left + 160, 46, std::max(80, w - left - 178), 26, TRUE);
    MoveWindow(g.editor, 12, top, left - 20, h - top - 45, TRUE);
    MoveWindow(g.plot, left + 8, top, w - left - 20, std::max(150, h - top - 225), TRUE);
    MoveWindow(g.log, left + 8, h - 205, w - left - 20, 160, TRUE);
    MoveWindow(g.status, 12, h - 32, w - 24, 24, TRUE);
}

static void paintPlot(Gui& g, HDC dc, RECT area) {
    HBRUSH bg = CreateSolidBrush(RGB(16, 25, 40));
    FillRect(dc, &area, bg);
    DeleteObject(bg);
    SetBkMode(dc, TRANSPARENT);
    SetTextColor(dc, RGB(220, 230, 245));
    SelectObject(dc, g.font);
    if (!g.result) {
        std::wstring text = L"Choose a model and run a calculation";
        TextOutW(dc, 30, 30, text.c_str(), int(text.size()));
        return;
    }
    auto& r = *g.result;
    auto& s = *r.sys;
    int sample = std::clamp(g.sample, 0, int(r.u.size()) - 1),
        field = std::clamp(g.field, 0, int(s.mdl->fields.size()) - 1);
    M pts;
    V values;
    std::vector<VI> triangles;
    bool spatial3 = false;
    if (s.mesh) {
        pts = s.mesh->p;
        spatial3 = pts.cols() == 3;
        int n = int(pts.rows());
        values = V::Zero(n);
        if (s.mixed && field == 0) {
            for (int f = 0; f < s.mdl->dim; ++f)
                values.array() += r.u[sample].segment(s.offset[f], n).array().square();
            values = values.array().sqrt();
        } else {
            int f = s.mixed ? s.mdl->dim : field;
            values = r.u[sample].segment(s.offset[f], n);
        }
        if (s.mesh->dim == 2)
            for (int i = 0; i < s.mesh->c.rows(); ++i)
                triangles.push_back({s.mesh->c(i, 0), s.mesh->c(i, 1), s.mesh->c(i, 2)});
        if (s.mesh->dim == 3)
            for (int i = 0; i < s.mesh->faces.rows(); ++i)
                triangles.push_back(
                    {s.mesh->faces(i, 0), s.mesh->faces(i, 1), s.mesh->faces(i, 2)});
    } else {
        int d = s.mdl->dim, nx = d == 1 ? 151 : 45, ny = d == 1 ? 1 : 45;
        pts = M::Zero(nx * ny, d);
        for (int i = 0; i < nx; ++i)
            for (int j = 0; j < ny; ++j) {
                pts(i * ny + j, 0) = s.mdl->par.at("Lx1") * i / (nx - 1);
                if (d > 1)
                    pts(i * ny + j, 1) = s.mdl->par.at("Lx2") * j / (ny - 1);
                if (d > 2)
                    pts(i * ny + j, 2) = s.mdl->par.at("Lx3") / 2;
            }
        values = s.values(r.u[sample], pts, field);
        for (int i = 0; i < nx - 1; ++i)
            for (int j = 0; j < ny - 1; ++j) {
                int k = i * ny + j;
                triangles.push_back({k, k + ny, k + 1});
                triangles.push_back({k + 1, k + ny, k + ny + 1});
            }
    }
    double lo = values.minCoeff(), hi = values.maxCoeff();
    M xy(pts.rows(), 2);
    for (int i = 0; i < pts.rows(); ++i) {
        double x = pts(i, 0), y = pts.cols() > 1 ? pts(i, 1) : values[i],
               z = spatial3 ? pts(i, 2) : 0;
        if (spatial3) {
            double xx = std::cos(g.yaw) * x - std::sin(g.yaw) * y;
            double yy = std::sin(g.pitch) * (std::sin(g.yaw) * x + std::cos(g.yaw) * y) -
                        std::cos(g.pitch) * z;
            x = xx;
            y = yy;
        }
        xy(i, 0) = x;
        xy(i, 1) = y;
    }
    double xlo = xy.col(0).minCoeff(), xhi = xy.col(0).maxCoeff(), ylo = xy.col(1).minCoeff(),
           yhi = xy.col(1).maxCoeff();
    int w = area.right, h = area.bottom;
    auto pt = [&](int i) {
        return POINT{LONG(w / 2 + g.zoom * (xy(i, 0) - (xhi + xlo) / 2) /
                                      std::max(1e-14, xhi - xlo) * (w - 100)),
                     LONG(h / 2 - g.zoom * (xy(i, 1) - (yhi + ylo) / 2) /
                                      std::max(1e-14, yhi - ylo) * (h - 100))};
    };
    if (!triangles.empty()) {
        if (spatial3)
            std::sort(triangles.begin(), triangles.end(), [&](const VI& a, const VI& b) {
                double da = 0, db = 0;
                for (int i : a)
                    da += std::cos(g.pitch) *
                              (std::sin(g.yaw) * pts(i, 0) + std::cos(g.yaw) * pts(i, 1)) +
                          std::sin(g.pitch) * pts(i, 2);
                for (int i : b)
                    db += std::cos(g.pitch) *
                              (std::sin(g.yaw) * pts(i, 0) + std::cos(g.yaw) * pts(i, 1)) +
                          std::sin(g.pitch) * pts(i, 2);
                return da < db;
            });
        HPEN pen = CreatePen(PS_SOLID, 1, RGB(35, 48, 62));
        auto oldPen = SelectObject(dc, pen);
        for (auto a : triangles) {
            POINT p[3] = {pt(a[0]), pt(a[1]), pt(a[2])};
            double v = (values[a[0]] + values[a[1]] + values[a[2]]) / 3;
            HBRUSH brush = CreateSolidBrush(color((v - lo) / std::max(1e-14, hi - lo)));
            auto old = SelectObject(dc, brush);
            Polygon(dc, p, 3);
            SelectObject(dc, old);
            DeleteObject(brush);
        }
        SelectObject(dc, oldPen);
        DeleteObject(pen);
    } else {
        HPEN pen = CreatePen(PS_SOLID, 2, RGB(90, 210, 245));
        auto old = SelectObject(dc, pen);
        std::vector<POINT> poly;
        for (int i = 0; i < pts.rows(); ++i)
            poly.push_back(pt(i));
        Polyline(dc, poly.data(), int(poly.size()));
        SelectObject(dc, old);
        DeleteObject(pen);
    }
    std::string title = s.mdl->fields[field] + (s.mixed && field == 0 ? " magnitude" : "") +
                        "   t=" + txt(r.t[sample]);
    std::wstring wt = wide(title);
    TextOutW(dc, 16, 12, wt.c_str(), int(wt.size()));
    std::ostringstream footer;
    footer << std::setprecision(6) << "min=" << lo << "  max=" << hi
           << (spatial3 ? " | boundary projection; drag to rotate"
                        : " | linear display between mesh vertices");
    wt = wide(footer.str());
    TextOutW(dc, 16, h - 25, wt.c_str(), int(wt.size()));
}

static LRESULT CALLBACK plotProc(HWND w, UINT msg, WPARAM wp, LPARAM lp) {
    Gui* g = app;
    if (!g)
        return DefWindowProcW(w, msg, wp, lp);
    if (msg == WM_PAINT) {
        PAINTSTRUCT ps;
        HDC dc = BeginPaint(w, &ps);
        RECT r;
        GetClientRect(w, &r);
        HDC mem = CreateCompatibleDC(dc);
        HBITMAP bm =
            CreateCompatibleBitmap(dc, std::max(1, int(r.right)), std::max(1, int(r.bottom)));
        auto old = SelectObject(mem, bm);
        try {
            paintPlot(*g, mem, r);
        } catch (const std::exception& e) {
            label(g->status, e.what());
        }
        BitBlt(dc, 0, 0, r.right, r.bottom, mem, 0, 0, SRCCOPY);
        SelectObject(mem, old);
        DeleteObject(bm);
        DeleteDC(mem);
        EndPaint(w, &ps);
        return 0;
    }
    if (msg == WM_ERASEBKGND)
        return 1;
    if (msg == WM_LBUTTONDOWN) {
        g->drag = true;
        g->lastx = short(LOWORD(lp));
        g->lasty = short(HIWORD(lp));
        SetCapture(w);
        return 0;
    }
    if (msg == WM_LBUTTONUP) {
        g->drag = false;
        ReleaseCapture();
        return 0;
    }
    if (msg == WM_MOUSEMOVE && g->drag) {
        int x = short(LOWORD(lp)), y = short(HIWORD(lp));
        g->yaw += (x - g->lastx) * .01;
        g->pitch = std::clamp(g->pitch + (y - g->lasty) * .01, -1.5, 1.5);
        g->lastx = x;
        g->lasty = y;
        InvalidateRect(w, nullptr, FALSE);
        return 0;
    }
    if (msg == WM_MOUSEWHEEL) {
        g->zoom = std::clamp(g->zoom * std::pow(1.1, short(HIWORD(wp)) / 120.), .2, 5.);
        InvalidateRect(w, nullptr, FALSE);
        return 0;
    }
    return DefWindowProcW(w, msg, wp, lp);
}

static void startWork(Gui& g, int action) {
    if (g.busy)
        return;
    std::string text = windowText(g.editor);
    J model = J::parse(text);
    std::string kind;
    int selection = int(SendMessageW(g.kind, CB_GETCURSEL, 0, 0));
    kind = selection == 1 ? "quad" : selection == 2 ? "time" : "mesh";
    int levels = std::clamp(std::stoi(windowText(g.levels)), 2, 5);
    if (g.worker.joinable())
        g.worker.join();
    g.busy = true;
    cancelled = false;
    for (int i = 0; i < g.actions.size(); ++i)
        EnableWindow(g.actions[i], i == 4);
    label(g.status, "Calculating...");
    g.worker = std::thread([&g, model, action, kind, levels]() {
        try {
            std::unique_ptr<Result> result;
            std::string msg;
            fs::path archive;
            if (action == 0 || action == 1) {
                result = std::make_unique<Result>(execute(model, action == 1));
                archive = result->archive(dataRoot() / "runs");
                msg = result->checks.dump(2) + "\n\n" + archive.string();
            } else if (action == 2) {
                System s(model);
                fs::path file = dataRoot() / "equations.cpp";
                std::string output = s.exportCpp();
                write(file, output);
                msg = "C++ equations saved: " + file.string() + "\n\n" + output.substr(0, 60000);
            } else if (action == 3) {
                PDE p(model);
                J w = p.weakForm();
                save(dataRoot() / "weak.json", w);
                msg = w.dump(2);
            } else if (action == 7 || action == 9) {
                auto stamp = std::to_string(GetTickCount64());
                J out = study(model, dataRoot() / "studies" / stamp, levels, kind, action == 9);
                msg = out.dump(2);
            } else
                throw Err("Unknown GUI task");
            {
                std::lock_guard lock(g.mutex);
                g.pending = std::move(result);
                g.message = std::move(msg);
                if (!archive.empty())
                    g.archive = archive;
            }
        } catch (const std::exception& e) {
            std::lock_guard lock(g.mutex);
            g.message = "Error: " + std::string(e.what());
        }
        PostMessageW(g.window, WM_APP + 1, 0, 0);
    });
}

static void fileDialog(Gui& g, int mode) {
    wchar_t file[32768] = {0};
    std::wstring root = dataRoot().wstring();
    OPENFILENAMEW of{};
    of.lStructSize = sizeof(of);
    of.hwndOwner = g.window;
    of.lpstrFilter = mode == 2 ? L"Saved run (manifest.json)\0manifest.json\0"
                               : L"PDE models (*.json)\0*.json\0All files\0*.*\0";
    of.lpstrFile = file;
    of.nMaxFile = 32768;
    of.lpstrInitialDir = root.c_str();
    of.lpstrDefExt = L"json";
    of.Flags =
        OFN_EXPLORER | OFN_NOCHANGEDIR | (mode == 1 ? OFN_OVERWRITEPROMPT : OFN_FILEMUSTEXIST);
    if (mode == 1 ? GetSaveFileNameW(&of) : GetOpenFileNameW(&of)) {
        if (mode == 1)
            save(file, J::parse(windowText(g.editor)));
        else if (mode == 0)
            label(g.editor, load(file).dump(2));
        else {
            fs::path path = fs::path(file).parent_path();
            g.pending = std::make_unique<Result>(restore(path));
            g.archive = path;
            g.message =
                "Loaded verified run: " + path.string() + "\n\n" + g.pending->checks.dump(2);
            label(g.editor, g.pending->sys->mdl->input.dump(2));
            PostMessageW(g.window, WM_APP + 1, 0, 0);
        }
    }
}

static LRESULT CALLBACK wndProc(HWND w, UINT msg, WPARAM wp, LPARAM lp) {
    Gui* g = app;
    if (!g)
        return DefWindowProcW(w, msg, wp, lp);
    try {
        if (msg == WM_GETMINMAXINFO) {
            auto mm = reinterpret_cast<MINMAXINFO*>(lp);
            mm->ptMinTrackSize = {920, 650};
            return 0;
        }
        if (msg == WM_SIZE) {
            layout(*g);
            return 0;
        }
        if (msg == WM_HSCROLL && HWND(lp) == g->slider) {
            g->sample = int(SendMessageW(g->slider, TBM_GETPOS, 0, 0));
            InvalidateRect(g->plot, nullptr, FALSE);
            return 0;
        }
        if (msg == WM_COMMAND) {
            int id = LOWORD(wp);
            if (id == 101 && HIWORD(wp) == CBN_SELCHANGE) {
                int i = int(SendMessageW(g->models, CB_GETCURSEL, 0, 0));
                if (i >= 0)
                    label(g->editor, getModel(g->names[i]).dump(2));
                return 0;
            }
            if (id == 102 && HIWORD(wp) == CBN_SELCHANGE) {
                g->field = int(SendMessageW(g->fields, CB_GETCURSEL, 0, 0));
                InvalidateRect(g->plot, nullptr, FALSE);
                return 0;
            }
            if (id >= 200 && id < 210) {
                int action = id - 200;
                if (action == 4) {
                    cancelled = true;
                    label(g->status, "Cancelling...");
                } else if (action == 5)
                    fileDialog(*g, false);
                else if (action == 6)
                    fileDialog(*g, 1);
                else if (action == 8)
                    fileDialog(*g, 2);
                else
                    startWork(*g, action);
                return 0;
            }
        }
        if (msg == WM_APP + 1) {
            if (g->worker.joinable())
                g->worker.join();
            {
                std::lock_guard lock(g->mutex);
                label(g->log, g->message);
                if (g->pending) {
                    g->result = std::move(g->pending);
                    g->sample = int(g->result->u.size()) - 1;
                    SendMessageW(g->fields, CB_RESETCONTENT, 0, 0);
                    for (auto s : g->result->sys->mdl->fields)
                        SendMessageW(g->fields, CB_ADDSTRING, 0, LPARAM(wide(s).c_str()));
                    SendMessageW(g->fields, CB_SETCURSEL, 0, 0);
                    g->field = 0;
                    SendMessageW(g->slider, TBM_SETRANGE, TRUE, MAKELONG(0, g->sample));
                    SendMessageW(g->slider, TBM_SETPOS, TRUE, g->sample);
                }
            }
            g->busy = false;
            for (auto a : g->actions)
                EnableWindow(a, TRUE);
            label(g->status,
                  g->message.rfind("Error:", 0) == 0 ? "Calculation failed; see log"
                  : g->result && !g->result->passed()
                      ? "Calculation finished with failed checks; see log"
                      : "Finished");
            InvalidateRect(g->plot, nullptr, FALSE);
            if (g->closing)
                DestroyWindow(w);
            return 0;
        }
        if (msg == WM_CLOSE) {
            if (g->busy) {
                g->closing = true;
                cancelled = true;
                label(g->status, "Cancelling before close...");
            } else
                DestroyWindow(w);
            return 0;
        }
        if (msg == WM_DESTROY) {
            PostQuitMessage(0);
            return 0;
        }
    } catch (const std::exception& e) {
        label(g->log, "Error: " + std::string(e.what()));
        label(g->status, "Action failed; see log");
    }
    return DefWindowProcW(w, msg, wp, lp);
}

int gui(const std::string& initial) {
    DWORD processes[2];
    if (GetConsoleProcessList(processes, 2) == 1)
        ShowWindow(GetConsoleWindow(), SW_HIDE);
    SetProcessDPIAware();
    INITCOMMONCONTROLSEX controls{sizeof(controls), ICC_BAR_CLASSES};
    InitCommonControlsEx(&controls);
    Gui g;
    app = &g;
    HINSTANCE instance = GetModuleHandleW(nullptr);
    WNDCLASSW wc{};
    wc.hInstance = instance;
    wc.hCursor = LoadCursorW(nullptr, MAKEINTRESOURCEW(32512));
    wc.hbrBackground = HBRUSH(COLOR_BTNFACE + 1);
    wc.lpfnWndProc = wndProc;
    wc.lpszClassName = L"PhyzoNativeWindow";
    RegisterClassW(&wc);
    wc.lpfnWndProc = plotProc;
    wc.lpszClassName = L"PhyzoNativePlot";
    RegisterClassW(&wc);
    g.window = CreateWindowW(L"PhyzoNativeWindow",
                             L"Phyzo - Native C++ PDE Laboratory",
                             WS_OVERLAPPEDWINDOW,
                             CW_USEDEFAULT,
                             CW_USEDEFAULT,
                             1280,
                             850,
                             nullptr,
                             nullptr,
                             instance,
                             nullptr);
    req(g.window != nullptr, "Cannot create native window");
    g.font = CreateFontW(-15,
                         0,
                         0,
                         0,
                         FW_NORMAL,
                         FALSE,
                         FALSE,
                         FALSE,
                         DEFAULT_CHARSET,
                         OUT_DEFAULT_PRECIS,
                         CLIP_DEFAULT_PRECIS,
                         CLEARTYPE_QUALITY,
                         DEFAULT_PITCH,
                         L"Segoe UI");
    g.mono = CreateFontW(-14,
                         0,
                         0,
                         0,
                         FW_NORMAL,
                         FALSE,
                         FALSE,
                         FALSE,
                         DEFAULT_CHARSET,
                         OUT_DEFAULT_PRECIS,
                         CLIP_DEFAULT_PRECIS,
                         CLEARTYPE_QUALITY,
                         FIXED_PITCH,
                         L"Consolas");
    auto control = [&](const wchar_t* cls, const wchar_t* text, DWORD style, int id) {
        HWND w = CreateWindowExW(0,
                                 cls,
                                 text,
                                 WS_CHILD | WS_VISIBLE | style,
                                 0,
                                 0,
                                 1,
                                 1,
                                 g.window,
                                 HMENU(INT_PTR(id)),
                                 instance,
                                 nullptr);
        SendMessageW(w, WM_SETFONT, WPARAM(g.font), TRUE);
        return w;
    };
    g.models = control(L"COMBOBOX", L"", CBS_DROPDOWNLIST | WS_VSCROLL, 101);
    g.fields = control(L"COMBOBOX", L"", CBS_DROPDOWNLIST, 102);
    g.editor = control(L"EDIT",
                       L"",
                       WS_BORDER | WS_VSCROLL | WS_HSCROLL | ES_MULTILINE | ES_AUTOVSCROLL |
                           ES_AUTOHSCROLL,
                       103);
    SendMessageW(g.editor, EM_SETLIMITTEXT, 1000000, 0);
    SendMessageW(g.editor, WM_SETFONT, WPARAM(g.mono), TRUE);
    g.log = control(
        L"EDIT", L"", WS_BORDER | WS_VSCROLL | ES_MULTILINE | ES_READONLY | ES_AUTOVSCROLL, 104);
    SendMessageW(g.log, EM_SETLIMITTEXT, 1000000, 0);
    SendMessageW(g.log, WM_SETFONT, WPARAM(g.mono), TRUE);
    g.plot = control(L"PhyzoNativePlot", L"", WS_BORDER, 105);
    g.slider = control(TRACKBAR_CLASSW, L"", TBS_HORZ, 106);
    g.status = control(L"STATIC", L"Ready - calculations run locally", 0, 107);
    g.kind = control(L"COMBOBOX", L"", CBS_DROPDOWNLIST, 108);
    for (auto s : {L"mesh", L"quad", L"time"})
        SendMessageW(g.kind, CB_ADDSTRING, 0, LPARAM(s));
    SendMessageW(g.kind, CB_SETCURSEL, 0, 0);
    g.levels = control(L"EDIT", L"3", WS_BORDER | ES_NUMBER, 109);
    for (auto text : {L"Calculate",
                      L"Equilibrium",
                      L"Equations",
                      L"Weak form",
                      L"Cancel",
                      L"Load model",
                      L"Save model",
                      L"Study",
                      L"Open run",
                      L"Study eq"})
        g.actions.push_back(control(L"BUTTON", text, BS_PUSHBUTTON, 200 + int(g.actions.size())));
    int selection = 0;
    for (auto [name, text] : examples())
        if (name.rfind("problem_", 0) != 0) {
            if (name == initial)
                selection = int(g.names.size());
            g.names.push_back(name);
            SendMessageW(g.models, CB_ADDSTRING, 0, LPARAM(wide(name).c_str()));
        }
    SendMessageW(g.models, CB_SETCURSEL, selection, 0);
    label(g.editor, getModel(initial).dump(2));
    layout(g);
    ShowWindow(g.window, SW_SHOW);
    UpdateWindow(g.window);
    MSG msg;
    while (GetMessageW(&msg, nullptr, 0, 0) > 0) {
        TranslateMessage(&msg);
        DispatchMessageW(&msg);
    }
    cancelled = true;
    if (g.worker.joinable())
        g.worker.join();
    DeleteObject(g.font);
    DeleteObject(g.mono);
    app = nullptr;
    return 0;
}
}
