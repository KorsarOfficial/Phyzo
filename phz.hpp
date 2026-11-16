#pragma once
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <json.hpp>
#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <complex>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <numeric>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace phz {
using J = nlohmann::ordered_json;
using V = Eigen::VectorXd;
using M = Eigen::MatrixXd;
using IM = Eigen::MatrixXi;
using SM = Eigen::SparseMatrix<double>;
using Tr = Eigen::Triplet<double>;
using LU = Eigen::SparseLU<SM>;
using VI = std::vector<int>;
namespace fs = std::filesystem;
constexpr double pi = 3.141592653589793238462643383279502884;

struct Err : std::runtime_error {
    using runtime_error::runtime_error;
};

struct Stop : Err {
    using Err::Err;
};

void req(bool b, const std::string& s);
std::string txt(double f);
std::string read(const fs::path& p);
void write(const fs::path& p, const std::string& s);
J load(const fs::path& p);
J npz(const fs::path& p);
std::pair<M, IM> readMesh(const fs::path& p, int dim);
void save(const fs::path& p, const J& j);
J arr(const V& a);
J arr(const M& a);
J arr(const IM& a);
V vec(const J& j);
M mat(const J& j);
IM imat(const J& j);
SM sparse(const M& a, double fTol = 0);
SM block(const SM& a, const VI& r, const VI& c);
V pick(const V& a, const VI& i);
void put(V& a, const VI& i, const V& b);
V solve(const SM& a, const V& b);
M jac(const std::function<V(const V&)>& fn, const V& u);
void keys(const J& j, const std::string& s, const std::string& ctx);
double number(const J& j, const std::string& s);
int integer(const J& j, int lo, int hi, const std::string& s);
std::string sha256(const fs::path& p);
fs::path exePath();
fs::path dataRoot();
extern std::atomic<bool> cancelled;
void checkpoint();

struct Node;
using E = std::shared_ptr<const Node>;

struct Node {
    std::string op, name;
    double n = 0;
    VI sh;
    std::vector<E> a;
};

E num(double x);
E sym(const std::string& s, VI sh = {}, const std::string& op = "sym");
bool zero(E e);
bool has(E e, const std::set<std::string>& ops);
bool depends(E e, const std::string& s);
E op(const std::string& s, std::vector<E> a, int dim = 1);
E add(E a, E b);
E mul(E a, E b);
E neg(E a);
E power(E a, double b);
E diff(E e, const std::string& s);
E subst(E e, const std::map<std::string, E>& d);
std::string str(E e, bool cpp = false);
J ir(E e);

struct Parser {
    std::string s;
    size_t at = 0;
    int dim = 1, depth = 0, nodes = 0;
    std::map<std::string, E> ns;

    Parser(std::map<std::string, E> d = {}, int n = 1) : dim(n), ns(std::move(d)) {
    }

    E parse(const std::string& t);
    E parse(const J& j);
    E sum();
    E prod();
    E unary();
    E pow();
    E atom();
    bool eat(char c);
    void ws();
};

double eval(E e, const std::map<std::string, double>& d);
E dx(E e, int k);
std::vector<E> components(E e, int dim);
std::pair<E, E> weak(E e, E v, int dim);

struct PDE {
    J input;
    int dim = 2;
    std::map<std::string, E> ns, eq;
    std::vector<std::string> fields;
    explicit PDE(J j);
    J weakForm() const;
};

struct Model {
    J input, dom, disc, run, stat, chk, bc, en, alg;
    int dim = 2;
    std::string type, name;
    std::vector<std::string> fields;
    std::map<std::string, double> par;
    std::map<std::string, E> ns, eq;
    std::unique_ptr<PDE> pde;
    explicit Model(J j);
    E parse(const J& j) const;
    double constant(const J& j) const;
    V at(E e, const M& x, double t) const;
};

struct Mesh {
    int dim = 0, geo = 0, n = 0;
    M p;
    IM c;
    V vol;
    std::vector<M> grad;
    SM mass, stiff;
    IM faces;
    std::map<std::string, IM> tags;
    Mesh() = default;
    Mesh(M aP, IM aC);
    explicit Mesh(const Model& m);
    void init();
    void setTags(const J& j);
    J json() const;
    std::pair<Mesh, VI> refine(bool sphere = false) const;
};

std::pair<V, V> gauss(int q);
std::pair<M, V> simplexQuad(int dim, int q);

struct Space {
    int pdeg = 1, q = 4, n = 0;
    M p, x, basis;
    IM c, alpha;
    V w;
    SM b, mass, stiff;
    std::vector<SM> g;
    explicit Space(const Mesh& mesh, int p = 1, int q = 4);
    std::pair<M, std::vector<M>> basisAt(const M& l) const;
    VI boundary(const Mesh& m, const std::string& tag) const;
};

struct Boundary {
    VI fix, free;
    std::vector<E> val, dt;

    struct Patch {
        M x;
        V w;
        SM b;
        E val, coef;
    };

    std::vector<Patch> natural;
};
struct Result;

struct System {
    std::shared_ptr<Model> mdl;
    std::unique_ptr<Mesh> mesh;
    std::vector<std::shared_ptr<Space>> spaces;
    M x;
    V w, one, massv, u0;
    SM b, mass, stiff;
    std::vector<SM> g;
    std::vector<VI> modes;
    std::vector<Boundary> bcs;
    std::unique_ptr<LU> massLU, poissonLU;
    std::vector<std::unique_ptr<LU>> freeLU;
    int n = 0, nf = 0, size = 0;
    VI active, offset;
    std::vector<std::string> compNames;
    std::vector<E> localWeak, localTests;

    struct Term {
        int row, k;
        E e;
    };

    struct JTerm {
        int row, k, col, l;
        E e;
    };

    std::vector<Term> terms;
    std::vector<JTerm> jterms;
    SM div;
    V gauge;
    bool mixed = false, exactSp = false, transient = true;
    std::vector<E> spRhs;
    explicit System(J j);
    void initScalar();
    void initMixed();
    void initBoundary();
    V at(E e, const M& points, double t) const;
    V lift(double t, const V& u) const;
    V boundaryDt(int f, double t) const;
    V boundaryLoad(int f, double t, const V& u) const;
    V project(const V& q) const;
    V rhs(double t, const V& u, bool neutral = true) const;
    M jacobian(double t, const V& u) const;
    std::pair<V, SM> assemble(
        double t, const V& u, const V* old = nullptr, double h = 0, bool derivative = true) const;
    V projectVelocity(const V& u) const;
    std::vector<V> algebraic(double t, const V& u) const;
    double energy(double t, const V& u) const;
    V norm(const V& u, bool dual = false) const;
    std::string exportCpp(bool fragment = false);
    V values(const V& u, const M& pts, int f = 0) const;
};

struct Result {
    std::shared_ptr<System> sys;
    std::vector<double> t;
    std::vector<V> u;
    J info = J::object(), diag = J::object(), checks = J::array(), scheme = J::object();
    bool stationary = false;
    void diagnose();
    bool passed() const;
    fs::path archive(const fs::path& root) const;
};

Result integrate(std::shared_ptr<System> s);
Result stationary(std::shared_ptr<System> s);
Result execute(J j, bool st = false);
J quadratureConfig(const Model& m);
J study(J j, const fs::path& root, int levels, const std::string& kind, bool st);
J verify(const fs::path& p);
Result restore(const fs::path& p);
Result replay(const fs::path& p);
J tableau(const std::string& method);
J scheme(const std::string& method, J tab = J(), int maxOrder = 6);
J heat(const Model& mdl, const std::string& method, double dt = -1);
const std::map<std::string, std::string>& examples();
const std::map<std::string, std::string>& referenceDocs();
J getModel(const std::string& s);
int selftest(const fs::path& out, int count = 1000, unsigned seed = 20261115);
int gui(const std::string& initial = "nseq2");
int cli(int argc, char** argv);
}
