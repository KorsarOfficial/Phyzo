#include "phz.hpp"

namespace phz {
static double fact(int n) {
    double f = 1;
    for (int i = 2; i <= n; ++i)
        f *= i;
    return f;
}

static IM cells(const std::vector<VI>& v, int d) {
    IM a(v.size(), d);
    for (int i = 0; i < v.size(); ++i)
        for (int j = 0; j < d; ++j)
            a(i, j) = v[i][j];
    return a;
}

static M points(const std::vector<V>& v, int g) {
    M a(v.size(), g);
    for (int i = 0; i < v.size(); ++i)
        a.row(i) = v[i];
    return a;
}

Mesh::Mesh(M aP, IM aC) : p(std::move(aP)), c(std::move(aC)) {
    init();
}

void Mesh::init() {
    dim = int(c.cols()) - 1;
    geo = int(p.cols());
    n = int(p.rows());
    req(dim >= 1 && dim <= 3 && geo >= dim && geo <= 3 && c.rows() > 0 && n > 0,
        "Invalid mesh shape");
    req(n <= 12000 && c.rows() <= 50000, "Mesh exceeds local limits");
    req(p.allFinite(), "Non-finite mesh coordinate");
    std::set<VI> seen;
    std::map<VI, int> fc;
    VI used(n), parent(n);
    std::iota(parent.begin(), parent.end(), 0);
    std::function<int(int)> root = [&](int a) {
        while (parent[a] != a) {
            parent[a] = parent[parent[a]];
            a = parent[a];
        }
        return a;
    };
    for (int e = 0; e < c.rows(); ++e) {
        VI a(dim + 1);
        for (int k = 0; k <= dim; ++k) {
            int i = c(e, k);
            req(i >= 0 && i < n, "Invalid mesh index");
            a[k] = i;
            used[i] = 1;
        }
        VI b = a;
        std::sort(b.begin(), b.end());
        req(std::adjacent_find(b.begin(), b.end()) == b.end(), "Repeated simplex vertex");
        req(seen.insert(b).second, "Duplicate simplex");
        for (int k = 1; k <= dim; ++k)
            parent[root(a[k])] = root(a[0]);
        for (int k = 0; k <= dim; ++k) {
            VI f = b;
            f.erase(f.begin() + k);
            req(++fc[f] <= 2, "Non-manifold mesh");
        }
    }
    req(std::all_of(used.begin(),
                    used.end(),
                    [](int v) {
                        return v == 1;
                    }),
        "Unused mesh nodes");
    for (int i = 0; i < n; ++i)
        req(root(i) == root(0), "Disconnected mesh");
    std::vector<VI> bd;
    for (auto [f, k] : fc)
        if (k == 1)
            bd.push_back(f);
    faces = cells(bd, dim);
    tags.clear();
    tags["all"] = faces;
    for (int ax = 0; ax < geo; ++ax) {
        double lo = p.col(ax).minCoeff(), hi = p.col(ax).maxCoeff(),
               tol = 1e-10 * (hi - lo) + 32 * std::numeric_limits<double>::epsilon() *
                                             std::max(std::abs(lo), std::abs(hi));
        if (hi <= lo)
            continue;
        for (int side = 0; side < 2; ++side) {
            std::vector<VI> sel;
            for (auto f : bd) {
                bool ok = true;
                for (int i : f)
                    ok &= std::abs(p(i, ax) - (side ? hi : lo)) <= tol;
                if (ok)
                    sel.push_back(f);
            }
            tags[std::string(1, "xyz"[ax]) + std::to_string(side)] = cells(sel, dim);
        }
    }
    vol.resize(c.rows());
    grad.clear();
    std::vector<Tr> mt, kt;
    for (int e = 0; e < c.rows(); ++e) {
        M h(geo, dim);
        for (int j = 0; j < dim; ++j)
            h.col(j) = p.row(c(e, j + 1)) - p.row(c(e, 0));
        M gram = h.transpose() * h;
        Eigen::SelfAdjointEigenSolver<M> es(gram);
        req(es.info() == Eigen::Success && es.eigenvalues()[0] > 1e-13 * es.eigenvalues()[dim - 1],
            "Degenerate or ill-conditioned simplex");
        double v = std::sqrt(gram.determinant()) / fact(dim);
        vol[e] = v;
        M r(dim + 1, dim);
        r.row(0).setConstant(-1);
        r.bottomRows(dim).setIdentity();
        M g = r * gram.inverse() * h.transpose();
        grad.push_back(g);
        for (int i = 0; i <= dim; ++i)
            for (int j = 0; j <= dim; ++j) {
                mt.emplace_back(c(e, i), c(e, j), v * (1 + (i == j)) / ((dim + 1) * (dim + 2)));
                kt.emplace_back(c(e, i), c(e, j), v * g.row(i).dot(g.row(j)));
            }
    }
    mass.resize(n, n);
    stiff.resize(n, n);
    mass.setFromTriplets(mt.begin(), mt.end());
    stiff.setFromTriplets(kt.begin(), kt.end());
}

void Mesh::setTags(const J& j) {
    req(j.is_object(), "tags must be an object");
    std::set<VI> all;
    for (int i = 0; i < faces.rows(); ++i) {
        VI f(dim);
        for (int k = 0; k < dim; ++k)
            f[k] = faces(i, k);
        all.insert(f);
    }
    for (auto it = j.begin(); it != j.end(); ++it) {
        IM a = imat(it.value());
        req(a.cols() == dim && !it.key().empty(), "Invalid boundary tag");
        std::set<VI> seen;
        for (int i = 0; i < a.rows(); ++i) {
            VI f(dim);
            for (int k = 0; k < dim; ++k)
                f[k] = a(i, k);
            std::sort(f.begin(), f.end());
            req(all.count(f) && seen.insert(f).second, "Invalid or repeated boundary facet");
            for (int k = 0; k < dim; ++k)
                a(i, k) = f[k];
        }
        if (tags.count(it.key())) {
            std::set<VI> old;
            auto b = tags.at(it.key());
            for (int i = 0; i < b.rows(); ++i) {
                VI f(dim);
                for (int k = 0; k < dim; ++k)
                    f[k] = b(i, k);
                old.insert(f);
            }
            req(old == seen, "Reserved boundary tag");
        }
        tags[it.key()] = a;
    }
}

J Mesh::json() const {
    J tg = J::object();
    for (auto [s, a] : tags)
        if (a.rows())
            tg[s] = arr(a);
    return J{{"type", "mesh"}, {"p", arr(p)}, {"c", arr(c)}, {"tags", tg}};
}

static std::pair<M, IM> importMesh(const fs::path& file, int dim) {
    auto [p0, c] = readMesh(file, dim);
    std::set<int> used;
    for (int i = 0; i < c.rows(); ++i)
        for (int k = 0; k < c.cols(); ++k) {
            req(c(i, k) >= 0 && c(i, k) < p0.rows(), "Invalid imported mesh index");
            used.insert(c(i, k));
        }
    M p(used.size(), p0.cols());
    std::map<int, int> remap;
    int i = 0;
    for (int k : used) {
        remap[k] = i;
        p.row(i++) = p0.row(k);
    }
    for (int i = 0; i < c.rows(); ++i)
        for (int k = 0; k < c.cols(); ++k)
            c(i, k) = remap.at(c(i, k));
    if (p.cols() > dim && p.rightCols(p.cols() - dim).norm() == 0)
        p = p.leftCols(dim).eval();
    return {p, c};
}

Mesh::Mesh(const Model& m) {
    int d = m.dim;
    J dom = m.dom;
    std::string ty = dom.value("type", "box");
    V len(d);
    for (int i = 0; i < d; ++i)
        len[i] = m.par.at("Lx" + std::to_string(i + 1));
    if (ty == "mesh") {
        if (dom.contains("path")) {
            auto [p0, c0] = importMesh(dom["path"].get<std::string>(), d);
            p = p0;
            c = c0;
        } else {
            p = mat(dom.at("p"));
            c = imat(dom.at("c"));
        }
    } else if (ty == "sphere") {
        req(d == 2, "Sphere is topologically two-dimensional");
        int k = integer(dom.value("n", J(2)), 0, 4, "sphere refinement");
        p.resize(6, 3);
        p << 1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1;
        c.resize(8, 3);
        c << 0, 2, 4, 2, 1, 4, 1, 3, 4, 3, 0, 4, 2, 0, 5, 1, 2, 5, 3, 1, 5, 0, 3, 5;
        init();
        for (int i = 0; i < k; ++i)
            *this = refine(true).first;
        p *= len[0];
    } else {
        req(ty == "box" || ty == "lshape" || ty == "ring", "Unknown domain");
        J ns = dom.value("n", J(VI(d, 8)));
        req(ns.is_array() && ns.size() == d, "Domain subdivisions dimension mismatch");
        VI n(d);
        for (int i = 0; i < d; ++i)
            n[i] = integer(ns[i], 1, 80, "domain subdivisions");
        std::vector<V> ps;
        std::vector<VI> cs;
        if (ty == "ring") {
            req(d == 2 && n[1] >= 6, "Ring requires dim=2 and 6 angular cells");
            for (int r = 0; r <= n[0]; ++r)
                for (int t = 0; t < n[1]; ++t) {
                    double rr = .2 + .3 * r / n[0], ang = 2 * pi * t / n[1];
                    V pt(2);
                    pt << len[0] * (.5 + rr * std::cos(ang)), len[1] * (.5 + rr * std::sin(ang));
                    ps.push_back(pt);
                }
            for (int r = 0; r < n[0]; ++r)
                for (int t = 0; t < n[1]; ++t) {
                    int a = r * n[1] + t, b = r * n[1] + (t + 1) % n[1];
                    cs.push_back({a, b, a + n[1]});
                    cs.push_back({b, b + n[1], a + n[1]});
                }
        } else {
            req(ty != "lshape" || (d == 2 && n[0] % 2 == 0 && n[1] % 2 == 0),
                "L-shape requires two even counts");
            int count = 1;
            for (int i : n)
                count *= i + 1;
            req(count <= 24000, "Requested mesh too large");
            auto ix = [&](VI v) {
                int k = 0;
                for (int i = 0; i < d; ++i)
                    k = k * (n[i] + 1) + v[i];
                return k;
            };
            for (int k = 0; k < count; ++k) {
                int t = k;
                V pt(d);
                for (int i = d - 1; i >= 0; --i) {
                    pt[i] = len[i] * (t % (n[i] + 1)) / n[i];
                    t /= n[i] + 1;
                }
                ps.push_back(pt);
            }
            int countC = 1;
            for (int i : n)
                countC *= i;
            for (int k = 0; k < countC; ++k) {
                VI base(d);
                int t = k;
                for (int i = d - 1; i >= 0; --i) {
                    base[i] = t % n[i];
                    t /= n[i];
                }
                if (ty == "lshape" && base[0] >= n[0] / 2 && base[1] >= n[1] / 2)
                    continue;
                VI order(d);
                std::iota(order.begin(), order.end(), 0);
                do {
                    VI pos = base, cell{ix(pos)};
                    for (int ax : order) {
                        ++pos[ax];
                        cell.push_back(ix(pos));
                    }
                    cs.push_back(cell);
                } while (std::next_permutation(order.begin(), order.end()));
            }
        }
        std::set<int> used;
        for (auto a : cs)
            used.insert(a.begin(), a.end());
        std::map<int, int> remap;
        std::vector<V> pp;
        for (int i : used) {
            remap[i] = int(pp.size());
            pp.push_back(ps[i]);
        }
        for (auto& a : cs)
            for (int& i : a)
                i = remap[i];
        p = points(pp, d);
        c = cells(cs, d + 1);
    }
    init();
    req(dim == m.dim, "Mesh topological dimension mismatch");
    if (dom.contains("tags"))
        setTags(dom["tags"]);
}

std::pair<Mesh, VI> Mesh::refine(bool sphere) const {
    std::vector<V> ps;
    for (int i = 0; i < n; ++i)
        ps.push_back(p.row(i));
    std::vector<VI> cs;
    VI parents;
    std::map<std::pair<int, int>, int> mid;
    auto midpoint = [&](int a, int b) {
        if (a > b)
            std::swap(a, b);
        auto k = std::make_pair(a, b);
        if (!mid.count(k)) {
            V x = (ps[a] + ps[b]) / 2;
            if (sphere)
                x.normalize();
            mid[k] = int(ps.size());
            ps.push_back(x);
        }
        return mid.at(k);
    };
    for (int e = 0; e < c.rows(); ++e) {
        int a = c(e, 0), b = c(e, 1);
        std::vector<VI> kids;
        if (dim == 1) {
            int m = midpoint(a, b);
            kids = {{a, m}, {m, b}};
        } else {
            int cc = c(e, 2), ab = midpoint(a, b), ac = midpoint(a, cc), bc = midpoint(b, cc);
            if (dim == 2)
                kids = {{a, ab, ac}, {ab, b, bc}, {ac, bc, cc}, {ab, bc, ac}};
            else {
                int d = c(e, 3), ad = midpoint(a, d), bd = midpoint(b, d), cd = midpoint(cc, d);
                kids = {{a, ab, ac, ad},
                        {b, ab, bc, bd},
                        {cc, ac, bc, cd},
                        {d, ad, bd, cd},
                        {ab, ac, ad, cd},
                        {ab, ac, bc, cd},
                        {ab, bc, bd, cd},
                        {ab, ad, bd, cd}};
            }
        }
        for (auto v : kids) {
            cs.push_back(v);
            parents.push_back(e);
        }
    }
    Mesh out(points(ps, geo), cells(cs, dim + 1));
    J jt = J::object();
    for (auto [tag, fc] : tags) {
        std::vector<VI> ff;
        for (int i = 0; i < fc.rows(); ++i) {
            int a = fc(i, 0);
            if (dim == 1)
                ff.push_back({a});
            else {
                int b = fc(i, 1), ab = midpoint(a, b);
                if (dim == 2) {
                    ff.push_back({a, ab});
                    ff.push_back({ab, b});
                } else {
                    int c = fc(i, 2), ac = midpoint(a, c), bc = midpoint(b, c);
                    ff.push_back({a, ab, ac});
                    ff.push_back({ab, b, bc});
                    ff.push_back({ac, bc, c});
                    ff.push_back({ab, bc, ac});
                }
            }
        }
        if (!ff.empty())
            jt[tag] = arr(cells(ff, dim));
    }
    out.setTags(jt);
    return {std::move(out), parents};
}

std::pair<V, V> gauss(int q) {
    req(q >= 1 && q <= 100, "Invalid Gauss order");
    M a = M::Zero(q, q);
    for (int i = 1; i < q; ++i)
        a(i, i - 1) = a(i - 1, i) = i / std::sqrt(4. * i * i - 1);
    Eigen::SelfAdjointEigenSolver<M> es(a);
    req(es.info() == Eigen::Success, "Quadrature eigensolve failed");
    V x = (es.eigenvalues().array() + 1) / 2;
    V w = es.eigenvectors().row(0).array().square();
    return {x, w};
}

std::pair<M, V> simplexQuad(int d, int q) {
    auto [x, w] = gauss(q);
    int nq = 1;
    for (int i = 0; i < d; ++i)
        nq *= q;
    M l(nq, d + 1);
    V v(nq);
    for (int i = 0; i < nq; ++i) {
        VI ix(d);
        int n = i;
        for (int k = d - 1; k >= 0; --k) {
            ix[k] = n % q;
            n /= q;
        }
        double rem = 1, wt = fact(d);
        for (int k = 0; k < d; ++k) {
            l(i, k + 1) = rem * x[ix[k]];
            wt *= w[ix[k]] * std::pow(1 - x[ix[k]], d - k - 1);
            rem *= 1 - x[ix[k]];
        }
        l(i, 0) = rem;
        v[i] = wt;
    }
    return {l, v};
}

static IM multi(int d, int p) {
    std::vector<VI> a;
    for (int i = 0; i <= d; ++i) {
        VI v(d + 1);
        v[i] = p;
        a.push_back(v);
    }
    VI v(d + 1);
    std::function<void(int, int)> rec = [&](int i, int left) {
        if (i == d) {
            v[i] = left;
            if (std::find(a.begin(), a.end(), v) == a.end())
                a.push_back(v);
            return;
        }
        for (int k = 0; k <= left; ++k) {
            v[i] = k;
            rec(i + 1, left - k);
        }
    };
    rec(0, p);
    return cells(a, d + 1);
}

Space::Space(const Mesh& mesh, int pp, int qq) : pdeg(pp), q(qq) {
    req(pp >= 1 && pp <= 4, "P degree must be 1..4");
    req(qq >= (pp == 1 ? (mesh.dim == 3 ? 3 : 2) : pp + 2) && qq <= 20,
        "Insufficient or excessive quadrature");
    alpha = multi(mesh.dim, pp);
    std::map<std::vector<std::pair<int, int>>, int> nodes;
    std::vector<V> ps;
    for (int i = 0; i < mesh.n; ++i) {
        ps.push_back(mesh.p.row(i));
        nodes[{{i, pp}}] = i;
    }
    c.resize(mesh.c.rows(), alpha.rows());
    for (int e = 0; e < c.rows(); ++e)
        for (int a = 0; a < alpha.rows(); ++a) {
            std::vector<std::pair<int, int>> key;
            V pt = V::Zero(mesh.geo);
            for (int k = 0; k <= mesh.dim; ++k) {
                if (alpha(a, k))
                    key.emplace_back(mesh.c(e, k), alpha(a, k));
                pt += double(alpha(a, k)) / pp * mesh.p.row(mesh.c(e, k)).transpose();
            }
            std::sort(key.begin(), key.end());
            if (!nodes.count(key)) {
                nodes[key] = int(ps.size());
                ps.push_back(pt);
            }
            c(e, a) = nodes.at(key);
        }
    p = points(ps, mesh.geo);
    n = int(p.rows());
    req(n <= 50000, "High-order space exceeds coefficient limit");
    auto [l, wt] = simplexQuad(mesh.dim, q);
    auto [bb, dd] = basisAt(l);
    basis = bb;
    int nq = int(l.rows()), na = int(alpha.rows());
    req(int64_t(nq) * c.rows() * na <= 6000000, "Quadrature exceeds memory limit");
    x.resize(nq * c.rows(), mesh.geo);
    w.resize(x.rows());
    std::vector<Tr> bt;
    std::vector<std::vector<Tr>> gt(mesh.geo);
    for (int e = 0; e < c.rows(); ++e) {
        M pts(mesh.dim + 1, mesh.geo);
        for (int k = 0; k <= mesh.dim; ++k)
            pts.row(k) = mesh.p.row(mesh.c(e, k));
        x.middleRows(e * nq, nq) = l * pts;
        w.segment(e * nq, nq) = mesh.vol[e] * wt;
        for (int k = 0; k < nq; ++k)
            for (int a = 0; a < na; ++a) {
                bt.emplace_back(e * nq + k, c(e, a), bb(k, a));
                for (int ax = 0; ax < mesh.geo; ++ax) {
                    double v = 0;
                    for (int j = 0; j <= mesh.dim; ++j)
                        v += dd[j](k, a) * mesh.grad[e](j, ax);
                    gt[ax].emplace_back(e * nq + k, c(e, a), v);
                }
            }
    }
    b.resize(x.rows(), n);
    b.setFromTriplets(bt.begin(), bt.end());
    g.resize(mesh.geo);
    SM wb = b;
    for (int j = 0; j < wb.outerSize(); ++j)
        for (SM::InnerIterator it(wb, j); it; ++it)
            it.valueRef() *= w[it.row()];
    mass = b.transpose() * wb;
    stiff.resize(n, n);
    for (int ax = 0; ax < mesh.geo; ++ax) {
        g[ax].resize(x.rows(), n);
        g[ax].setFromTriplets(gt[ax].begin(), gt[ax].end());
        SM wg = g[ax];
        for (int j = 0; j < wg.outerSize(); ++j)
            for (SM::InnerIterator it(wg, j); it; ++it)
                it.valueRef() *= w[it.row()];
        stiff += g[ax].transpose() * wg;
    }
    if (pp == 1) {
        mass = mesh.mass;
        stiff = mesh.stiff;
    }
}

std::pair<M, std::vector<M>> Space::basisAt(const M& l) const {
    M b = M::Ones(l.rows(), alpha.rows());
    std::vector<M> d(l.cols(), M::Zero(l.rows(), alpha.rows()));
    for (int a = 0; a < alpha.rows(); ++a)
        for (int q = 0; q < l.rows(); ++q) {
            V f(l.cols()), df(l.cols());
            for (int k = 0; k < l.cols(); ++k) {
                int n = alpha(a, k);
                f[k] = 1;
                df[k] = 0;
                for (int j = 0; j < n; ++j) {
                    double v = (pdeg * l(q, k) - j) / (n - j);
                    df[k] = df[k] * v + f[k] * pdeg / (n - j);
                    f[k] *= v;
                }
            }
            b(q, a) = f.prod();
            for (int k = 0; k < l.cols(); ++k) {
                double v = df[k];
                for (int j = 0; j < l.cols(); ++j)
                    if (j != k)
                        v *= f[j];
                d[k](q, a) = v;
            }
        }
    return {b, d};
}

VI Space::boundary(const Mesh& mesh, const std::string& tag) const {
    req(mesh.tags.count(tag), "Unknown boundary: " + tag);
    std::set<int> ids;
    auto fc = mesh.tags.at(tag);
    std::set<VI> faces;
    for (int i = 0; i < fc.rows(); ++i) {
        VI f(mesh.dim);
        for (int k = 0; k < mesh.dim; ++k)
            f[k] = fc(i, k);
        std::sort(f.begin(), f.end());
        faces.insert(f);
    }
    for (int e = 0; e < mesh.c.rows(); ++e)
        for (int miss = 0; miss <= mesh.dim; ++miss) {
            VI f;
            for (int k = 0; k <= mesh.dim; ++k)
                if (k != miss)
                    f.push_back(mesh.c(e, k));
            std::sort(f.begin(), f.end());
            if (faces.count(f))
                for (int a = 0; a < alpha.rows(); ++a)
                    if (alpha(a, miss) == 0)
                        ids.insert(c(e, a));
        }
    return VI(ids.begin(), ids.end());
}
}
