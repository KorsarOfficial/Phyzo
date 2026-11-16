#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <windows.h>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
namespace fs = std::filesystem;

static std::wstring quote(const fs::path& p) {
    std::wstring s = p.wstring();
    if (s.find(L'"') != std::wstring::npos)
        throw std::runtime_error("Invalid path");
    return L"\"" + s + L"\"";
}

static int process(const fs::path& exe, std::wstring args) {
    std::wstring line = quote(exe) + L" " + args;
    STARTUPINFOW si{};
    si.cb = sizeof(si);
    PROCESS_INFORMATION pi{};
    if (!CreateProcessW(
            exe.c_str(), line.data(), nullptr, nullptr, FALSE, 0, nullptr, nullptr, &si, &pi))
        throw std::runtime_error("Cannot start build tool");
    WaitForSingleObject(pi.hProcess, INFINITE);
    DWORD code = 1;
    GetExitCodeProcess(pi.hProcess, &code);
    CloseHandle(pi.hProcess);
    CloseHandle(pi.hThread);
    return int(code);
}

int main(int argc, char** argv) {
    try {
        fs::path source = fs::absolute(fs::path(__FILE__)).parent_path();
        if (!fs::exists(source / "phz.hpp"))
            source = fs::current_path();
        if (!fs::exists(source / "phz.hpp")) {
            wchar_t module[32768]{};
            DWORD n = GetModuleFileNameW(nullptr, module, 32768);
            if (n && n < 32768)
                source = fs::path(module).parent_path().parent_path().parent_path() / "Phyzo";
        }
        if (!fs::exists(source / "phz.hpp"))
            throw std::runtime_error("Run builder from the C++ source directory");
        fs::path sdk = source.parent_path() / "Phyzo-native";
        if (const wchar_t* env = _wgetenv(L"PHYZO_SDK"))
            sdk = env;
        sdk = fs::absolute(sdk);
        for (auto name : {"json.hpp",
                          "eigen-3.4.0/Eigen/Dense",
                          "boost_1_86_0/boost/multiprecision/cpp_int.hpp",
                          "sundials/include/cvode/cvode.h",
                          "miniz/miniz.h",
                          "tinyxml2/tinyxml2.h"})
            if (!fs::exists(sdk / "deps" / name))
                throw std::runtime_error(std::string("Missing pinned SDK dependency: ") + name +
                                         "; set PHYZO_SDK to the supplied native SDK");
        fs::path cmake;
        wchar_t path[32768];
        if (SearchPathW(nullptr, L"cmake.exe", nullptr, 32768, path, nullptr))
            cmake = path;
        if (cmake.empty())
            for (auto root : {fs::path(L"C:/Program Files/Microsoft Visual Studio"),
                              fs::path(L"C:/Program Files (x86)/Microsoft Visual Studio")})
                if (fs::exists(root))
                    for (auto year : fs::directory_iterator(root))
                        if (year.is_directory())
                            for (auto edition : fs::directory_iterator(year)) {
                                auto candidate =
                                    edition.path() /
                                    "Common7/IDE/CommonExtensions/Microsoft/CMake/CMake/bin/cmake.exe";
                                if (fs::exists(candidate))
                                    cmake = candidate;
                            }
        if (cmake.empty())
            throw std::runtime_error("Install Visual Studio C++ tools including CMake");
        auto cmakePath = [](const fs::path& p) {
            std::string s = p.generic_string();
            if (s.find_first_of("\";$\n\r") != std::string::npos)
                throw std::runtime_error("Unsupported CMake path character");
            return "\"" + s + "\"";
        };
        fs::path config = sdk / "buildsys", build = sdk / "build/phyzo";
        fs::create_directories(config);
        std::ofstream out(config / "CMakeLists.txt");
        out << "cmake_minimum_required(VERSION 3.24)\nproject(Phyzo LANGUAGES C CXX)\nset(CMAKE_CXX_STANDARD 20)\nfile(GLOB PHZ_SRC CONFIGURE_DEPENDS "
            << cmakePath(source / "*.cpp")
            << ")\nlist(FILTER PHZ_SRC EXCLUDE REGEX \"/build.cpp$\")\nadd_executable(Phyzo ${PHZ_SRC}";
        for (auto s : {"miniz/miniz.c",
                       "miniz/miniz_tdef.c",
                       "miniz/miniz_tinfl.c",
                       "miniz/miniz_zip.c",
                       "tinyxml2/tinyxml2.cpp"})
            out << " " << cmakePath(sdk / "deps" / s);
        out << ")\ntarget_include_directories(Phyzo PRIVATE";
        for (auto s : {"", "eigen-3.4.0", "boost_1_86_0", "sundials/include", "miniz", "tinyxml2"})
            out << " " << cmakePath(sdk / "deps" / s);
        out << ")\ntarget_compile_options(Phyzo PRIVATE /utf-8 /bigobj /W3 /EHsc /MP3)\ntarget_compile_definitions(Phyzo PRIVATE NOMINMAX _CRT_SECURE_NO_WARNINGS)\ntarget_link_directories(Phyzo PRIVATE "
            << cmakePath(sdk / "deps/sundials/lib")
            << ")\ntarget_link_libraries(Phyzo PRIVATE bcrypt user32 gdi32 comctl32 comdlg32 shell32 sundials_cvode_static sundials_nvecserial_static sundials_sunmatrixdense_static sundials_sunlinsoldense_static sundials_core_static)\nset_target_properties(Phyzo PROPERTIES RUNTIME_OUTPUT_DIRECTORY_RELEASE "
            << cmakePath(sdk / "bin") << ")\n";
        out << "file(GLOB PHZ_SNAPSHOT " << cmakePath(source / "*.cpp") << " "
            << cmakePath(source / "*.hpp")
            << ")\nadd_custom_command(TARGET Phyzo POST_BUILD COMMAND ${CMAKE_COMMAND} -E make_directory "
            << cmakePath(sdk / "bin/src")
            << " COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PHZ_SNAPSHOT} "
            << cmakePath(sdk / "bin/src") << " VERBATIM)\n";
        out.close();
        int rc = process(cmake, L"-S " + quote(config) + L" -B " + quote(build) + L" -A x64");
        if (rc)
            return rc;
        rc = process(cmake, L"--build " + quote(build) + L" --config Release --parallel 3");
        if (rc)
            return rc;
        if (argc > 1 && std::string(argv[1]) == "--test")
            return process(sdk / "bin/Phyzo.exe",
                           L"selftest --count 1000 -o " + quote(sdk / "reports/native-1000.json"));
        std::cout << (sdk / "bin/Phyzo.exe").string() << "\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << e.what() << "\n";
        return 1;
    }
}
