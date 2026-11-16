#include "phz.hpp"
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <bcrypt.h>

namespace phz {
std::atomic<bool> cancelled{false};
static thread_local double cutoff = INFINITY;

double deadline(double absolute) {
    double old = cutoff;
    cutoff = absolute;
    return old;
}

void checkpoint() {
    if (cancelled.load())
        throw Stop("Calculation cancelled");
    if (std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch())
            .count() >= cutoff)
        throw Stop("Quadrature control timeout");
}

void req(bool b, const std::string& s) {
    if (!b)
        throw Err(s);
}

std::string txt(double f) {
    req(std::isfinite(f), "Non-finite number");
    std::ostringstream o;
    o << std::setprecision(17) << f;
    return o.str();
}

std::string read(const fs::path& p) {
    std::ifstream f(p, std::ios::binary);
    req(bool(f), "Cannot read " + p.string());
    return {std::istreambuf_iterator<char>(f), {}};
}

void write(const fs::path& p, const std::string& s) {
    if (p.has_parent_path())
        fs::create_directories(p.parent_path());
    std::ofstream f(p, std::ios::binary);
    req(bool(f), "Cannot write " + p.string());
    f.write(s.data(), s.size());
    req(bool(f), "Incomplete write " + p.string());
}

J load(const fs::path& p) {
    return J::parse(read(p));
}

void save(const fs::path& p, const J& j) {
    write(p, j.dump(2));
}

J arr(const V& a) {
    J j = J::array();
    for (auto f : a)
        j.push_back(f);
    return j;
}

J arr(const M& a) {
    J j = J::array();
    for (int i = 0; i < a.rows(); ++i) {
        J r = J::array();
        for (int k = 0; k < a.cols(); ++k)
            r.push_back(a(i, k));
        j.push_back(r);
    }
    return j;
}

J arr(const IM& a) {
    J j = J::array();
    for (int i = 0; i < a.rows(); ++i) {
        J r = J::array();
        for (int k = 0; k < a.cols(); ++k)
            r.push_back(a(i, k));
        j.push_back(r);
    }
    return j;
}

double number(const J& j, const std::string& s) {
    req(j.is_number() && !j.is_boolean(), s + " must be numeric");
    double f = j.get<double>();
    req(std::isfinite(f), s + " must be finite");
    return f;
}

int integer(const J& j, int lo, int hi, const std::string& s) {
    req(j.is_number_integer() && !j.is_boolean(), s + " must be an integer");
    auto n = j.get<int64_t>();
    req(n >= lo && n <= hi, s + " outside permitted range");
    return int(n);
}

V vec(const J& j) {
    req(j.is_array(), "Expected numeric array");
    V a(j.size());
    for (int i = 0; i < a.size(); ++i)
        a[i] = number(j[i], "entry");
    return a;
}

M mat(const J& j) {
    req(j.is_array() && !j.empty() && j[0].is_array(), "Expected matrix");
    M a(j.size(), j[0].size());
    for (int i = 0; i < a.rows(); ++i) {
        req(j[i].size() == a.cols(), "Ragged matrix");
        for (int k = 0; k < a.cols(); ++k)
            a(i, k) = number(j[i][k], "matrix entry");
    }
    return a;
}

IM imat(const J& j) {
    req(j.is_array() && !j.empty() && j[0].is_array(), "Expected integer matrix");
    IM a(j.size(), j[0].size());
    for (int i = 0; i < a.rows(); ++i) {
        req(j[i].size() == a.cols(), "Ragged cells");
        for (int k = 0; k < a.cols(); ++k)
            a(i, k) = integer(j[i][k], 0, 100000000, "cell entry");
    }
    return a;
}

SM sparse(const M& a, double tol) {
    return a.sparseView(tol, 1);
}

SM block(const SM& a, const VI& r, const VI& c) {
    VI ri(a.rows(), -1), ci(a.cols(), -1);
    for (int i = 0; i < r.size(); ++i)
        ri.at(r[i]) = i;
    for (int i = 0; i < c.size(); ++i)
        ci.at(c[i]) = i;
    std::vector<Tr> t;
    for (int k = 0; k < a.outerSize(); ++k)
        for (SM::InnerIterator v(a, k); v; ++v)
            if (ri[v.row()] >= 0 && ci[v.col()] >= 0)
                t.emplace_back(ri[v.row()], ci[v.col()], v.value());
    SM m(r.size(), c.size());
    m.setFromTriplets(t.begin(), t.end());
    return m;
}

V pick(const V& a, const VI& ix) {
    V v(ix.size());
    for (int k = 0; k < ix.size(); ++k)
        v[k] = a[ix[k]];
    return v;
}

void put(V& a, const VI& ix, const V& b) {
    req(ix.size() == b.size(), "Index size mismatch");
    for (int k = 0; k < ix.size(); ++k)
        a[ix[k]] = b[k];
}

V solve(const SM& a, const V& b) {
    req(a.rows() == b.size() && a.cols() == b.size(), "Linear system shape mismatch");
    if (!b.size())
        return V();
    LU lu;
    lu.compute(a);
    req(lu.info() == Eigen::Success, "Singular linear system");
    V x = lu.solve(b);
    req(lu.info() == Eigen::Success && x.allFinite(), "Linear solve failed");
    req((a * x - b).norm() <= 2e-7 * (1 + b.norm() + a.norm() * x.norm()),
        "Linear solve residual too large");
    return x;
}

M jac(const std::function<V(const V&)>& fn, const V& u) {
    V f = fn(u);
    M j(f.size(), u.size());
    for (int i = 0; i < u.size(); ++i) {
        checkpoint();
        double h = std::cbrt(std::numeric_limits<double>::epsilon()) * std::max(1., std::abs(u[i]));
        bool ok = false, side = false;
        V ds;
        for (int k = 0; k < 24; ++k) {
            V a = u, b = u, fa, fb;
            a[i] += h;
            b[i] -= h;
            if (a[i] == u[i] || b[i] == u[i])
                break;
            bool ap = true, bp = true;
            try {
                fa = fn(a);
            } catch (const Stop&) {
                throw;
            } catch (const Err&) {
                ap = false;
            }
            try {
                fb = fn(b);
            } catch (const Stop&) {
                throw;
            } catch (const Err&) {
                bp = false;
            }
            if (ap && bp) {
                j.col(i) = (fa - fb) / (a[i] - b[i]);
                ok = true;
                break;
            }
            if (ap) {
                ds = (fa - f) / (a[i] - u[i]);
                side = true;
            }
            if (bp) {
                ds = (f - fb) / (u[i] - b[i]);
                side = true;
            }
            h /= 4;
        }
        if (!ok) {
            req(side, "No admissible Jacobian perturbation");
            j.col(i) = ds;
        }
    }
    req(j.allFinite(), "Non-finite Jacobian");
    return j;
}

void keys(const J& j, const std::string& s, const std::string& ctx) {
    req(j.is_object(), ctx + " must be an object");
    std::set<std::string> names;
    std::istringstream in(s);
    std::string n;
    while (in >> n)
        names.insert(n);
    for (auto it = j.begin(); it != j.end(); ++it)
        req(names.count(it.key()), "Unknown " + ctx + " key: " + it.key());
}

std::string sha256(const fs::path& p) {
    BCRYPT_ALG_HANDLE h = nullptr;
    req(BCryptOpenAlgorithmProvider(&h, BCRYPT_SHA256_ALGORITHM, nullptr, 0) >= 0,
        "SHA256 unavailable");
    std::string s = read(p);
    std::array<unsigned char, 32> out{};
    auto status = BCryptHash(h,
                             nullptr,
                             0,
                             reinterpret_cast<PUCHAR>(s.data()),
                             ULONG(s.size()),
                             out.data(),
                             ULONG(out.size()));
    BCryptCloseAlgorithmProvider(h, 0);
    req(status >= 0, "SHA256 failed");
    std::ostringstream ss;
    for (auto c : out)
        ss << std::hex << std::setw(2) << std::setfill('0') << int(c);
    return ss.str();
}

fs::path exePath() {
    std::wstring s(32768, L'\0');
    DWORD n = GetModuleFileNameW(nullptr, s.data(), DWORD(s.size()));
    req(n && n < s.size(), "Cannot locate executable");
    s.resize(n);
    return fs::path(s);
}

fs::path dataRoot() {
    const char* s = std::getenv("PHYZO_DATA");
    return s ? fs::path(s) : fs::path("C:/Phyzo-native/data");
}
}
