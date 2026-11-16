#include "phz.hpp"
#include <miniz.h>
#include <tinyxml2.h>
#include <regex>
#include <cstring>

namespace phz {
struct Binary {
    std::string data;
    size_t pos = 0;
    bool little = true;

    explicit Binary(const fs::path& p) : data(read(p)) {
    }

    std::string line() {
        req(pos < data.size(), "Truncated mesh header");
        size_t end = data.find('\n', pos);
        if (end == std::string::npos)
            end = data.size();
        std::string s = data.substr(pos, end - pos);
        pos = std::min(end + 1, data.size());
        if (!s.empty() && s.back() == '\r')
            s.pop_back();
        return s;
    }

    std::string next() {
        std::string s;
        do {
            s = line();
        } while (s.empty());
        return s;
    }

    uint64_t u(int n) {
        req(n >= 1 && n <= 8 && size_t(n) <= data.size() - pos, "Truncated binary mesh");
        uint64_t x = 0;
        for (int i = 0; i < n; ++i)
            x |= uint64_t((unsigned char)data[pos++]) << (8 * (little ? i : n - 1 - i));
        return x;
    }

    double real(int n) {
        uint64_t u0 = u(n);
        if (n == 8) {
            double x;
            std::memcpy(&x, &u0, 8);
            req(std::isfinite(x), "Non-finite mesh point");
            return x;
        }
        uint32_t bits = uint32_t(u0);
        float x;
        std::memcpy(&x, &bits, 4);
        req(std::isfinite(x), "Non-finite mesh point");
        return x;
    }
};

static std::pair<M, IM> meshArrays(const std::vector<V>& p, const std::vector<VI>& c, int d) {
    req(!p.empty() && !c.empty(), "No requested simplex cells");
    M pts(p.size(), 3);
    IM cells(c.size(), d + 1);
    for (int i = 0; i < p.size(); ++i)
        pts.row(i) = p[i];
    for (int i = 0; i < c.size(); ++i)
        for (int k = 0; k <= d; ++k)
            cells(i, k) = c[i][k];
    return {pts, cells};
}

static std::pair<M, IM> gmshBinary(const fs::path& file, int dim) {
    Binary in(file);
    req(in.next() == "$MeshFormat", "Missing Gmsh header");
    std::istringstream header(in.next());
    double version = 0;
    int binary = 0, sz = 0;
    header >> version >> binary >> sz;
    req(binary == 1 && sz == 8 && (version >= 2 && version < 3 || version >= 4.1 && version < 5),
        "Unsupported binary Gmsh header");
    uint64_t endian = in.u(4);
    if (endian == 0x1000000)
        in.little = false;
    else
        req(endian == 1, "Invalid Gmsh byte order");
    req(in.next() == "$EndMeshFormat", "Invalid Gmsh format terminator");
    std::vector<V> pts;
    std::vector<VI> cells;
    std::map<uint64_t, int> ids;
    int target = dim == 1 ? 1 : dim == 2 ? 2 : 4;
    std::map<int, int> sizes{{1, 2},
                             {2, 3},
                             {3, 4},
                             {4, 4},
                             {5, 8},
                             {6, 6},
                             {7, 5},
                             {8, 3},
                             {9, 6},
                             {10, 9},
                             {11, 10},
                             {12, 27},
                             {13, 18},
                             {14, 14},
                             {15, 1},
                             {16, 8},
                             {17, 20},
                             {18, 15},
                             {19, 13},
                             {20, 9},
                             {21, 10},
                             {22, 12},
                             {23, 15},
                             {24, 15},
                             {25, 21},
                             {26, 4},
                             {27, 5},
                             {28, 6},
                             {29, 20},
                             {30, 35},
                             {31, 56}};
    auto node = [&](uint64_t tag) {
        req(!ids.count(tag) && pts.size() < 120000, "Duplicate or excessive Gmsh node");
        ids[tag] = int(pts.size());
        V p(3);
        for (int i = 0; i < 3; ++i)
            p[i] = in.real(8);
        pts.push_back(p);
    };
    auto element = [&](int type, int tags, int bytes) {
        req(sizes.count(type), "Unsupported high-order Gmsh binary element type");
        in.u(bytes);
        for (int i = 0; i < tags; ++i)
            in.u(bytes);
        VI c;
        for (int i = 0; i < sizes.at(type); ++i) {
            uint64_t tag = in.u(bytes);
            if (type == target) {
                req(ids.count(tag), "Unknown binary Gmsh node");
                c.push_back(ids.at(tag));
            }
        }
        if (type == target)
            cells.push_back(c);
    };
    while (in.pos < in.data.size()) {
        std::string section = in.next();
        if (section == "$Nodes") {
            if (version < 3) {
                uint64_t n = std::stoull(in.next());
                req(n > 0 && n <= 120000, "Invalid Gmsh node count");
                for (uint64_t i = 0; i < n; ++i)
                    node(in.u(4));
            } else {
                uint64_t blocks = in.u(sz), n = in.u(sz);
                in.u(sz);
                in.u(sz);
                req(n > 0 && n <= 120000 && blocks <= n, "Invalid Gmsh node blocks");
                for (uint64_t b = 0; b < blocks; ++b) {
                    int d = int(in.u(4));
                    in.u(4);
                    int param = int(in.u(4));
                    uint64_t count = in.u(sz);
                    req(count <= n && d >= 0 && d <= 3 && (param == 0 || param == 1),
                        "Invalid Gmsh node block");
                    std::vector<uint64_t> tags(count);
                    for (auto& t : tags)
                        t = in.u(sz);
                    for (auto tag : tags) {
                        node(tag);
                        if (param)
                            for (int k = 0; k < d; ++k)
                                in.real(8);
                    }
                }
            }
            req(in.next() == "$EndNodes", "Missing Gmsh node terminator");
        } else if (section == "$Elements") {
            if (version < 3) {
                uint64_t n = std::stoull(in.next()), read = 0;
                req(n <= 1000000, "Too many Gmsh elements");
                while (read < n) {
                    int type = int(in.u(4));
                    uint64_t count = in.u(4);
                    int tags = int(in.u(4));
                    req(count > 0 && count <= n - read && tags >= 0 && tags <= 100,
                        "Invalid Gmsh element block");
                    for (uint64_t i = 0; i < count; ++i)
                        element(type, tags, 4);
                    read += count;
                }
            } else {
                uint64_t blocks = in.u(sz), n = in.u(sz);
                in.u(sz);
                in.u(sz);
                req(n <= 1000000 && blocks <= n, "Too many Gmsh elements");
                for (uint64_t b = 0; b < blocks; ++b) {
                    in.u(4);
                    in.u(4);
                    int type = int(in.u(4));
                    uint64_t count = in.u(sz);
                    req(count <= n, "Invalid Gmsh element block");
                    for (uint64_t i = 0; i < count; ++i)
                        element(type, 0, sz);
                }
            }
            req(in.next() == "$EndElements", "Missing Gmsh element terminator");
            break;
        } else if (section == "$Entities") {
            req(version >= 4, "Unexpected binary entities");
            uint64_t counts[4];
            for (auto& n : counts) {
                n = in.u(sz);
                req(n <= 1000000, "Too many Gmsh entities");
            }
            for (int d = 0; d < 4; ++d)
                for (uint64_t i = 0; i < counts[d]; ++i) {
                    in.u(4);
                    for (int k = 0; k < (d == 0 ? 3 : 6); ++k)
                        in.real(8);
                    uint64_t np = in.u(sz);
                    req(np <= 10000, "Too many Gmsh physical groups");
                    for (uint64_t k = 0; k < np; ++k)
                        in.u(4);
                    if (d) {
                        uint64_t nb = in.u(sz);
                        req(nb <= 1000000, "Too many Gmsh bounding entities");
                        for (uint64_t k = 0; k < nb; ++k)
                            in.u(4);
                    }
                }
            req(in.next() == "$EndEntities", "Missing Gmsh entities terminator");
        } else if (section.size() > 1 && section[0] == '$' && section.rfind("$End", 0) != 0) {
            std::string stop = "$End" + section.substr(1);
            while (in.next() != stop) {
            }
        }
    }
    return meshArrays(pts, cells, dim);
}

static std::pair<M, IM> vtkBinary(const fs::path& file, int dim) {
    Binary in(file);
    in.little = false;
    std::string first = in.next();
    double version = 0;
    size_t at = first.rfind(' ');
    if (at != std::string::npos)
        version = std::stod(first.substr(at + 1));
    in.line();
    req(in.next() == "BINARY", "Missing VTK binary marker");
    std::vector<V> pts;
    std::vector<VI> all, cells;
    while (in.pos < in.data.size()) {
        std::istringstream row(in.next());
        std::string key;
        row >> key;
        if (key == "POINTS") {
            int n = 0;
            std::string type;
            row >> n >> type;
            req(n > 0 && n <= 120000 && (type == "float" || type == "double"),
                "Invalid VTK points");
            for (int i = 0; i < n; ++i) {
                V p(3);
                for (int k = 0; k < 3; ++k)
                    p[k] = in.real(type == "double" ? 8 : 4);
                pts.push_back(p);
            }
        } else if (key == "CELLS") {
            int n = 0, total = 0;
            row >> n >> total;
            req(n > 0 && n <= 1000000 && total > 0 && total <= 16000000, "Invalid VTK cells");
            if (version >= 5) {
                std::string header = in.next();
                int bytes = header.find("64") != std::string::npos ? 8 : 4;
                req(header.rfind("OFFSETS", 0) == 0, "Missing VTK offsets");
                std::vector<uint64_t> ends(n);
                for (auto& k : ends)
                    k = in.u(bytes);
                header = in.next();
                bytes = header.find("64") != std::string::npos ? 8 : 4;
                req(header.rfind("CONNECTIVITY", 0) == 0 && ends[0] == 0 && ends.back() == total,
                    "Invalid VTK offsets");
                VI conn(total);
                for (auto& k : conn) {
                    uint64_t v = in.u(bytes);
                    req(v < pts.size(), "Invalid VTK connectivity");
                    k = int(v);
                }
                for (int i = 1; i < n; ++i) {
                    req(ends[i] > ends[i - 1] && ends[i] <= total, "Invalid VTK offset order");
                    all.emplace_back(conn.begin() + ends[i - 1], conn.begin() + ends[i]);
                }
            } else {
                int count = 0;
                for (int i = 0; i < n; ++i) {
                    int size = int(in.u(4));
                    req(size > 0 && size <= 100 && count + size + 1 <= total,
                        "Invalid VTK cell size");
                    VI cell(size);
                    for (auto& k : cell) {
                        uint64_t v = in.u(4);
                        req(v < pts.size(), "Invalid VTK node");
                        k = int(v);
                    }
                    all.push_back(cell);
                    count += size + 1;
                }
                req(count == total, "Invalid VTK cell word count");
            }
        } else if (key == "CELL_TYPES") {
            int n = 0;
            row >> n;
            req(bool(row) && n == all.size(), "VTK type count mismatch");
            for (int i = 0; i < n; ++i) {
                int type = int(in.u(4));
                if (type == (dim == 1 ? 3 : dim == 2 ? 5 : 10)) {
                    req(all[i].size() == dim + 1, "Invalid VTK simplex");
                    cells.push_back(all[i]);
                }
            }
            break;
        }
    }
    return meshArrays(pts, cells, dim);
}

static uint64_t unsignedAt(const std::string& s, size_t at, int n, bool little = true) {
    req(n >= 1 && n <= 8 && at <= s.size() && size_t(n) <= s.size() - at, "Truncated binary array");
    uint64_t x = 0;
    for (int i = 0; i < n; ++i)
        x |= uint64_t((unsigned char)s[at + i]) << (8 * (little ? i : n - 1 - i));
    return x;
}

static V
numeric(const std::string& data, std::string type, bool little = true, bool strict = true) {
    int n = 0;
    bool floating = false, signedInt = false;
    if (type == "Float64" || type == "f8") {
        n = 8;
        floating = true;
    } else if (type == "Float32" || type == "f4") {
        n = 4;
        floating = true;
    } else if (type == "Int64" || type == "i8") {
        n = 8;
        signedInt = true;
    } else if (type == "Int32" || type == "i4") {
        n = 4;
        signedInt = true;
    } else if (type == "Int16" || type == "i2") {
        n = 2;
        signedInt = true;
    } else if (type == "Int8" || type == "i1") {
        n = 1;
        signedInt = true;
    } else if (type == "UInt64" || type == "u8")
        n = 8;
    else if (type == "UInt32" || type == "u4")
        n = 4;
    else if (type == "UInt16" || type == "u2")
        n = 2;
    else if (type == "UInt8" || type == "u1" || type == "b1")
        n = 1;
    else
        throw Err("Unsupported numeric type: " + type);
    req(data.size() % n == 0 && data.size() / n <= 16000000, "Invalid or excessive binary array");
    V out(data.size() / n);
    for (int i = 0; i < out.size(); ++i) {
        uint64_t v = unsignedAt(data, size_t(i) * n, n, little);
        if (floating) {
            if (n == 8) {
                double x;
                std::memcpy(&x, &v, 8);
                out[i] = x;
            } else {
                uint32_t k = uint32_t(v);
                float x;
                std::memcpy(&x, &k, 4);
                out[i] = x;
            }
        } else if (signedInt) {
            if (n < 8 && (v & (uint64_t(1) << (8 * n - 1))))
                v |= (~uint64_t(0)) << (8 * n);
            int64_t k;
            std::memcpy(&k, &v, 8);
            req(k >= -9007199254740992LL && k <= 9007199254740992LL,
                "Integer array exceeds exact binary64 range");
            out[i] = double(k);
        } else {
            req(v <= 9007199254740992ULL, "Integer array exceeds exact binary64 range");
            out[i] = double(v);
        }
    }
    if (strict)
        req(out.allFinite(), "Non-finite imported array");
    return out;
}

static J npy(const std::string& s) {
    req(s.size() >= 10 && s.substr(0, 6) == std::string("\x93NUMPY", 6), "Invalid NPY header");
    int version = (unsigned char)s[6], bytes = version == 1                   ? 2
                                               : version == 2 || version == 3 ? 4
                                                                              : 0;
    req(bytes > 0, "Unsupported NPY version");
    size_t len = size_t(unsignedAt(s, 8, bytes)), offset = 8 + bytes + len;
    req(offset <= s.size() && len <= 65536, "Invalid NPY header length");
    std::string header = s.substr(8 + bytes, len);
    std::smatch m;
    req(std::regex_search(header, m, std::regex("['\"]descr['\"]\\s*:\\s*['\"]([^'\"]+)['\"]")),
        "NPY dtype missing");
    std::string type = m[1];
    req(type.size() >= 3, "Invalid NPY dtype");
    bool little = type[0] != '>';
    type.erase(0, 1);
    V data = numeric(s.substr(offset), type, little, false);
    req(std::regex_search(header, m, std::regex("['\"]shape['\"]\\s*:\\s*\\(([^)]*)\\)")),
        "NPY shape missing");
    std::string shape = m[1];
    std::replace(shape.begin(), shape.end(), ',', ' ');
    std::istringstream in(shape);
    VI dims;
    int n;
    int64_t count = 1;
    while (in >> n) {
        req(n >= 0 && n <= 16000000 && dims.size() < 6, "Invalid NPY shape");
        count *= n;
        req(count <= 16000000, "NPY array too large");
        dims.push_back(n);
    }
    req(count == data.size(), "NPY shape/data mismatch");
    bool fort = header.find("True") != std::string::npos;
    VI ix(dims.size());
    std::function<J(int)> nested = [&](int d) -> J {
        if (d == dims.size()) {
            int at = 0;
            if (fort) {
                int stride = 1;
                for (int k = 0; k < dims.size(); ++k) {
                    at += ix[k] * stride;
                    stride *= dims[k];
                }
            } else
                for (int k = 0; k < dims.size(); ++k)
                    at = at * dims[k] + ix[k];
            return std::isfinite(data[at]) ? J(data[at]) : J();
        }
        J a = J::array();
        for (int i = 0; i < dims[d]; ++i) {
            ix[d] = i;
            a.push_back(nested(d + 1));
        }
        return a;
    };
    return nested(0);
}

J npz(const fs::path& path) {
    std::string bytes = read(path);
    mz_zip_archive zip{};
    req(mz_zip_reader_init_mem(&zip, bytes.data(), bytes.size(), 0) != 0,
        "Cannot read NPZ archive");

    struct End {
        mz_zip_archive* z;

        ~End() {
            mz_zip_reader_end(z);
        }
    } end{&zip};

    J result = J::object();
    uint64_t total = 0;
    for (mz_uint i = 0; i < mz_zip_reader_get_num_files(&zip); ++i) {
        mz_zip_archive_file_stat stat{};
        req(mz_zip_reader_file_stat(&zip, i, &stat) != 0, "Invalid ZIP entry");
        std::string name = stat.m_filename;
        req(name.size() > 4 && name.substr(name.size() - 4) == ".npy" &&
                name.find_first_of("/\\:") == std::string::npos,
            "Invalid NPZ entry name");
        total += stat.m_uncomp_size;
        req(stat.m_uncomp_size <= 128000000 && total <= 512000000, "NPZ exceeds memory limit");
        std::string data(size_t(stat.m_uncomp_size), '\0');
        req(mz_zip_reader_extract_to_mem(&zip, i, data.data(), data.size(), 0) != 0,
            "Corrupt NPZ entry");
        name.resize(name.size() - 4);
        req(!result.contains(name), "Duplicate NPZ array");
        result[name] = npy(data);
    }
    return result;
}

static std::string base64(const std::string& s) {
    std::string out;
    unsigned v = 0;
    int bits = 0;
    for (unsigned char c : s) {
        if (std::isspace(c))
            continue;
        if (c == '=') {
            v = 0;
            bits = 0;
            continue;
        }
        int k = c >= 'A' && c <= 'Z'   ? c - 'A'
                : c >= 'a' && c <= 'z' ? c - 'a' + 26
                : c >= '0' && c <= '9' ? c - '0' + 52
                : c == '+'             ? 62
                : c == '/'             ? 63
                                       : -1;
        req(k >= 0, "Invalid base64 data");
        v = (v << 6) | k;
        bits += 6;
        if (bits >= 8) {
            bits -= 8;
            out += char((v >> bits) & 255);
        }
    }
    return out;
}

static std::pair<M, IM> vtu(const fs::path& file, int dim) {
    std::string text = read(file), appended;
    bool raw = false;
    size_t start = text.find("<AppendedData");
    if (start != std::string::npos) {
        size_t endTag = text.find('>', start), stop = text.rfind("</AppendedData>"),
               under = text.find('_', endTag);
        req(stop != std::string::npos && under < stop, "Invalid appended VTK data");
        raw = text.substr(start, endTag - start).find("raw") != std::string::npos;
        appended = text.substr(under + 1, stop - under - 1);
        text.replace(under + 1, stop - under - 1, "");
    }
    tinyxml2::XMLDocument doc;
    req(doc.Parse(text.data(), text.size()) == tinyxml2::XML_SUCCESS, "Invalid VTU XML");
    auto root = doc.FirstChildElement("VTKFile");
    req(root, "Missing VTKFile");
    bool little = std::string(root->Attribute("byte_order") ? root->Attribute("byte_order")
                                                            : "LittleEndian") == "LittleEndian";
    int hb = std::string(root->Attribute("header_type") ? root->Attribute("header_type")
                                                        : "UInt32") == "UInt64"
                 ? 8
                 : 4;
    std::string compressor = root->Attribute("compressor") ? root->Attribute("compressor") : "";
    req(compressor.empty() || compressor == "vtkZLibDataCompressor",
        "VTU requires uncompressed or zlib data");
    auto array = [&](tinyxml2::XMLElement* a) {
        req(a && a->Attribute("type"), "Missing VTU array");
        std::string fmt = a->Attribute("format") ? a->Attribute("format") : "ascii";
        if (fmt == "ascii") {
            std::istringstream in(a->GetText() ? a->GetText() : "");
            std::vector<double> v;
            double x;
            while (in >> x) {
                req(v.size() < 16000000 && std::isfinite(x), "Invalid VTU numeric array");
                v.push_back(x);
            }
            return Eigen::Map<V>(v.data(), v.size()).eval();
        }
        std::string data;
        if (fmt == "binary")
            data = base64(a->GetText() ? a->GetText() : "");
        else {
            req(fmt == "appended", "Unsupported VTU array encoding");
            uint64_t at = a->Unsigned64Attribute("offset", UINT64_MAX);
            req(at < appended.size(), "Invalid VTU offset");
            data = raw ? appended.substr(size_t(at)) : base64(appended.substr(size_t(at)));
        }
        std::string body;
        if (compressor.empty()) {
            size_t n = size_t(unsignedAt(data, 0, hb, little));
            req(n <= 128000000 && n <= data.size() - hb, "Truncated VTU block");
            body = data.substr(hb, n);
        } else {
            uint64_t nb = unsignedAt(data, 0, hb, little), bs = unsignedAt(data, hb, hb, little),
                     last = unsignedAt(data, 2 * hb, hb, little);
            req(nb <= 100000 && bs <= 128000000 && last <= bs && nb * bs <= 128000000,
                "Invalid VTU compressed layout");
            size_t offset = (3 + size_t(nb)) * hb;
            req(offset <= data.size(), "Truncated VTU compression header");
            for (size_t i = 0; i < nb; ++i) {
                size_t sz = size_t(unsignedAt(data, (3 + i) * hb, hb, little));
                req(sz <= data.size() - offset, "Truncated VTU compressed block");
                mz_ulong want = mz_ulong(i + 1 == nb && last ? last : bs);
                std::string block(want, '\0');
                req(mz_uncompress((unsigned char*)block.data(),
                                  &want,
                                  (const unsigned char*)data.data() + offset,
                                  mz_ulong(sz)) == MZ_OK &&
                        want == block.size(),
                    "Corrupt VTU zlib block");
                body += block;
                offset += sz;
            }
        }
        return numeric(body, a->Attribute("type"), little);
    };
    auto grid = root->FirstChildElement("UnstructuredGrid");
    req(grid, "VTU must contain an unstructured grid");
    std::vector<V> pts;
    std::vector<VI> cc;
    for (auto piece = grid->FirstChildElement("Piece"); piece;
         piece = piece->NextSiblingElement("Piece")) {
        auto pp = piece->FirstChildElement("Points"), cells = piece->FirstChildElement("Cells");
        req(pp && cells, "Missing VTU points or cells");
        auto pa = pp->FirstChildElement("DataArray");
        V p = array(pa);
        int nc = pa->IntAttribute("NumberOfComponents", 3);
        req(nc >= dim && nc <= 3 && p.size() % nc == 0, "VTU point shape");
        int offset = int(pts.size());
        for (int i = 0; i < p.size() / nc; ++i)
            pts.push_back(p.segment(i * nc, nc));
        std::map<std::string, V> arrays;
        for (auto a = cells->FirstChildElement("DataArray"); a;
             a = a->NextSiblingElement("DataArray")) {
            req(a->Attribute("Name"), "Unnamed VTU cell array");
            arrays[a->Attribute("Name")] = array(a);
        }
        req(arrays.count("connectivity") && arrays.count("offsets") && arrays.count("types"),
            "Missing VTU connectivity");
        auto conn = arrays.at("connectivity"), ends = arrays.at("offsets"),
             types = arrays.at("types");
        req(ends.size() == types.size(), "VTU cell count mismatch");
        int prev = 0;
        for (int i = 0; i < ends.size(); ++i) {
            req(ends[i] == std::floor(ends[i]) && ends[i] >= prev && ends[i] <= conn.size(),
                "Invalid VTU cell offset");
            int end = int(ends[i]);
            if (types[i] == (dim == 1 ? 3 : dim == 2 ? 5 : 10)) {
                req(end - prev == dim + 1, "Invalid simplex cell size");
                VI cell;
                for (int j = prev; j < end; ++j) {
                    req(conn[j] == std::floor(conn[j]) && conn[j] >= 0 && conn[j] < p.size() / nc,
                        "Invalid VTU node index");
                    cell.push_back(offset + int(conn[j]));
                }
                cc.push_back(cell);
            }
            prev = end;
        }
    }
    req(!pts.empty() && !cc.empty(), "No requested simplices in VTU");
    M p(pts.size(), pts[0].size());
    for (int i = 0; i < pts.size(); ++i) {
        req(pts[i].size() == p.cols(), "VTU coordinate dimensions differ");
        p.row(i) = pts[i];
    }
    IM c(cc.size(), dim + 1);
    for (int i = 0; i < cc.size(); ++i)
        for (int k = 0; k <= dim; ++k)
            c(i, k) = cc[i][k];
    return {p, c};
}

std::pair<M, IM> readMesh(const fs::path& file, int dim) {
    std::string ext = file.extension().string();
    std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
    if (ext == ".json") {
        auto j = load(file);
        return {mat(j.at("p")), imat(j.at("c"))};
    }
    if (ext == ".vtu")
        return vtu(file, dim);
    if (ext == ".msh") {
        std::istringstream header(read(file));
        std::string marker;
        double version = 0;
        int binary = 0, sz = 0;
        header >> marker >> version >> binary >> sz;
        if (marker == "$MeshFormat" && binary == 1)
            return gmshBinary(file, dim);
    }
    if (ext == ".vtk") {
        std::istringstream header(read(file));
        std::string line;
        std::getline(header, line);
        std::getline(header, line);
        std::getline(header, line);
        if (line.find("BINARY") != std::string::npos)
            return vtkBinary(file, dim);
    }
    std::istringstream in(read(file));
    std::string s;
    std::vector<V> pts;
    std::vector<VI> cs, all;
    std::map<uint64_t, int> ids;
    int wanted = dim == 1 ? 1 : dim == 2 ? 2 : 4;
    double version = 2;
    if (ext == ".msh") {
        while (in >> s) {
            if (s == "$MeshFormat") {
                int binary = 0, sz = 0;
                in >> version >> binary >> sz;
                req(bool(in), "Invalid Gmsh format record");
                req(binary == 0 && (version >= 2 && version < 3 || version >= 4.1 && version < 5),
                    "Gmsh supports ASCII 2.x and 4.1; export binary meshes as ASCII or VTU");
            } else if (s == "$Nodes") {
                uint64_t blocks = 1, n = 0, a, b;
                if (version < 3)
                    in >> n;
                else
                    in >> blocks >> n >> a >> b;
                req(n > 0 && n <= 120000 && blocks <= n, "Invalid Gmsh node counts");
                for (uint64_t block = 0; block < blocks; ++block) {
                    uint64_t count = n;
                    int d = 0, entity = 0, param = 0;
                    if (version >= 4)
                        in >> d >> entity >> param >> count;
                    req(count <= n && d >= 0 && d <= 3, "Invalid Gmsh node block");
                    std::vector<uint64_t> tags(count);
                    if (version >= 4)
                        for (auto& tag : tags)
                            in >> tag;
                    for (uint64_t i = 0; i < count; ++i) {
                        uint64_t tag = tags[i];
                        if (version < 3)
                            in >> tag;
                        V p(3);
                        in >> p[0] >> p[1] >> p[2];
                        if (param)
                            for (int k = 0; k < d; ++k) {
                                double v;
                                in >> v;
                            }
                        req(in.good() && !ids.count(tag), "Malformed or duplicate Gmsh node");
                        ids[tag] = int(pts.size());
                        pts.push_back(p);
                    }
                }
            } else if (s == "$Elements") {
                uint64_t blocks = 1, n = 0, a = 0, b = 0;
                if (version < 3)
                    in >> n;
                else
                    in >> blocks >> n >> a >> b;
                req(n <= 1000000 && blocks <= n, "Invalid Gmsh element counts");
                std::getline(in, s);
                for (uint64_t block = 0; block < blocks; ++block) {
                    uint64_t count = n;
                    int type = 0, d = 0, entity = 0;
                    if (version >= 4) {
                        in >> d >> entity >> type >> count;
                        std::getline(in, s);
                    }
                    req(count <= n, "Invalid Gmsh element block");
                    for (uint64_t i = 0; i < count; ++i) {
                        req(bool(std::getline(in, s)), "Truncated Gmsh element list");
                        std::istringstream row(s);
                        uint64_t id;
                        row >> id;
                        if (version < 3) {
                            int nt = 0;
                            row >> type >> nt;
                            req(bool(row), "Malformed Gmsh element record");
                            req(nt >= 0 && nt <= 100, "Invalid Gmsh element tags");
                            for (int k = 0; k < nt; ++k)
                                row >> id;
                        }
                        if (type == wanted) {
                            VI cell(dim + 1);
                            for (auto& k : cell) {
                                uint64_t node;
                                req(bool(row >> node) && ids.count(node), "Unknown Gmsh node");
                                k = ids.at(node);
                            }
                            cs.push_back(cell);
                        }
                    }
                }
            }
        }
    } else if (ext == ".vtk") {
        while (in >> s) {
            if (s == "BINARY")
                throw Err("Legacy VTK import requires ASCII; use VTU for binary data");
            if (s == "POINTS") {
                int n = 0;
                in >> n >> s;
                req(bool(in), "Malformed VTK point header");
                req(n > 0 && n <= 120000, "Invalid VTK node count");
                for (int i = 0; i < n; ++i) {
                    V p(3);
                    in >> p[0] >> p[1] >> p[2];
                    req(bool(in), "Truncated VTK points");
                    pts.push_back(p);
                }
            } else if (s == "CELLS" || s == "POLYGONS" || s == "LINES") {
                std::string section = s;
                int n = 0, total = 0;
                in >> n >> total;
                req(bool(in), "Malformed VTK cell header");
                req(n > 0 && n <= 1000000 && total >= 0 && total <= 16000000,
                    "Invalid VTK cell count");
                in >> std::ws;
                if (in.peek() == 'O') {
                    std::string token, type;
                    in >> token >> type;
                    req(token == "OFFSETS", "Invalid VTK offsets header");
                    VI offsets(n);
                    for (auto& v : offsets)
                        in >> v;
                    in >> token >> type;
                    req(token == "CONNECTIVITY" && offsets[0] == 0 && offsets.back() == total,
                        "Invalid VTK connectivity header");
                    VI conn(total);
                    for (auto& v : conn)
                        in >> v;
                    req(bool(in), "Truncated VTK connectivity");
                    for (int i = 1; i < n; ++i) {
                        req(offsets[i] > offsets[i - 1] && offsets[i] <= total,
                            "Invalid VTK offset");
                        all.emplace_back(conn.begin() + offsets[i - 1], conn.begin() + offsets[i]);
                    }
                    continue;
                }
                for (int i = 0; i < n; ++i) {
                    int size = 0;
                    in >> size;
                    req(size > 0 && size <= 100, "Invalid VTK cell size");
                    VI cell(size);
                    for (auto& k : cell)
                        in >> k;
                    req(bool(in), "Truncated VTK cells");
                    if (section == "CELLS")
                        all.push_back(cell);
                    else if (size == dim + 1 &&
                             (section == "LINES" && dim == 1 || section == "POLYGONS" && dim == 2))
                        cs.push_back(cell);
                }
            } else if (s == "CELL_TYPES") {
                int n = 0;
                in >> n;
                req(bool(in) && n == all.size(), "VTK cell types mismatch");
                for (int i = 0; i < n; ++i) {
                    int type = 0;
                    in >> type;
                    req(bool(in), "Truncated VTK cell types");
                    if (type == (dim == 1 ? 3 : dim == 2 ? 5 : 10)) {
                        req(all[i].size() == dim + 1, "VTK simplex size mismatch");
                        cs.push_back(all[i]);
                    }
                }
            }
        }
    } else
        throw Err("Mesh formats: JSON, Gmsh ASCII 2.x/4.1, VTK ASCII, VTU ASCII/binary/zlib");
    req(!pts.empty() && !cs.empty(), "No requested simplex cells in mesh");
    M p(pts.size(), 3);
    for (int i = 0; i < pts.size(); ++i)
        p.row(i) = pts[i];
    IM c(cs.size(), dim + 1);
    for (int i = 0; i < cs.size(); ++i)
        for (int j = 0; j <= dim; ++j)
            c(i, j) = cs[i][j];
    return {p, c};
}
}
