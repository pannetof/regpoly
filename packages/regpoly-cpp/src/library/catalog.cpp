#include "catalog.h"

#include "factory.h"
#include "params.h"
#include "transformation.h"

#include <yaml-cpp/yaml.h>

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <iostream>
#include <limits>
#include <regex>
#include <sstream>
#include <sys/stat.h>

namespace regpoly_catalog {

namespace {

const std::regex kIdRe("^[a-z0-9][a-z0-9-]*$");

const std::vector<std::string> kKnownTargets = {
    "tested_generator", "primitive_generator"};

// Walk up the directory tree from `start` looking for a sibling
// `docs/` directory; return its parent (the repo root) or empty if
// not found within 5 levels.
std::filesystem::path find_repo_root(const std::filesystem::path& start) {
    auto here = start;
    for (int i = 0; i < 5; ++i) {
        auto docs = here / "docs";
        std::error_code ec;
        if (std::filesystem::is_directory(docs, ec)
            && here.filename() != "docs") {
            return here;
        }
        if (here == here.parent_path()) break;
        here = here.parent_path();
    }
    return {};
}

double stat_mtime(const std::string& path) {
    struct stat st{};
    if (::stat(path.c_str(), &st) != 0) return 0.0;
    return static_cast<double>(st.st_mtime)
         + static_cast<double>(st.st_mtim.tv_nsec) * 1e-9;
}

// Try parsing a scalar string as an integer, handling decimal AND
// hex (0x... — yaml-cpp leaves these as strings). Mirrors what
// PyYAML's safe_load does: 0xFOO is loaded as int.
bool try_parse_int(const std::string& s, int64_t& out) {
    if (s.empty()) return false;
    try {
        if (s.size() > 2 && s[0] == '0' && (s[1] == 'x' || s[1] == 'X')) {
            // Hex literal. Use stoull so large values (>= 1<<63) keep
            // their bit pattern when cast back to int64_t — matches the
            // Python factory's "store uint64 bitmask as int64".
            uint64_t u = std::stoull(s.substr(2), nullptr, 16);
            out = static_cast<int64_t>(u);
            return true;
        }
        if (s[0] == '-' || std::isdigit(static_cast<unsigned char>(s[0]))) {
            size_t pos = 0;
            int64_t v = std::stoll(s, &pos, 10);
            if (pos == s.size()) { out = v; return true; }
        }
    } catch (...) {}
    return false;
}

// Recursive convert YAML scalar/seq/map → ParamValue / TemperingStep / etc.
ParamValue node_to_param_value(const YAML::Node& n) {
    if (n.IsSequence()) {
        std::vector<int64_t> vec;
        bool all_int = true;
        for (const auto& el : n) {
            int64_t v = 0;
            if (el.IsScalar()) {
                std::string s = el.as<std::string>();
                if (!try_parse_int(s, v)) {
                    // Last attempt: yaml-cpp's own conversion (handles
                    // unquoted decimal ints reliably).
                    try { v = el.as<int64_t>(); } catch (...) { all_int = false; break; }
                }
            } else {
                all_int = false;
                break;
            }
            vec.push_back(v);
        }
        if (all_int) return ParamValue::make_int_list(std::move(vec));
        // Non-int list: stringify the YAML emit so we don't lose data.
        YAML::Emitter em; em << n;
        return ParamValue::make_string(em.c_str());
    }
    if (n.IsScalar()) {
        const auto& tag = n.Tag();
        // yaml-cpp marks pure-bool scalars with tag "!" or "tag:yaml.org,2002:bool"
        // depending on quoting. Try bool first (yes/no/true/false).
        std::string s = n.as<std::string>();
        // booleans
        if (tag == "tag:yaml.org,2002:bool"
            || s == "true" || s == "false"
            || s == "True" || s == "False") {
            try { return ParamValue::make_bool(n.as<bool>()); }
            catch (...) {}
        }
        // ints (decimal); leave 0x… as string for the runtime factory.
        if (!s.empty() && (s[0] == '-' || std::isdigit(s[0]))
            && s.find('.') == std::string::npos
            && s.find('e') == std::string::npos
            && !(s.size() > 1 && s[0] == '0'
                 && (s[1] == 'x' || s[1] == 'X'))) {
            try { return ParamValue::make_int(n.as<int64_t>()); }
            catch (...) {}
        }
        return ParamValue::make_string(s);
    }
    if (n.IsMap()) {
        // Catalog params shouldn't have nested maps; serialize defensively.
        YAML::Emitter em; em << n;
        return ParamValue::make_string(em.c_str());
    }
    return ParamValue::make_string("");
}

ParamMap node_to_param_map(const YAML::Node& n) {
    ParamMap out;
    if (!n || !n.IsMap()) return out;
    for (auto it : n) {
        out.emplace(it.first.as<std::string>(),
                    node_to_param_value(it.second));
    }
    return out;
}

std::string trim_copy(std::string s) {
    auto issp = [](unsigned char c) { return std::isspace(c); };
    auto b = std::find_if_not(s.begin(), s.end(), issp);
    auto e = std::find_if_not(s.rbegin(), s.rend(), issp).base();
    if (b >= e) return "";
    return std::string(b, e);
}

// Push a ParamValue into a runtime Params object so the dry-build
// step can call create_generator / create_transformation. Mirrors
// the Python factory's lenient int/string coercion.
void apply_param(Params& dst, const std::string& key, const ParamValue& v) {
    switch (v.kind) {
        case ParamKind::Int:
            dst.set_int(key, v.int_val);
            break;
        case ParamKind::Bool:
            dst.set_bool(key, v.bool_val);
            break;
        case ParamKind::String: {
            // Try hex/decimal parse; fall back to set_string.
            const auto& s = v.string_val;
            try {
                if (s.size() > 2 && s[0] == '0'
                    && (s[1] == 'x' || s[1] == 'X')) {
                    dst.set_int(key, std::stoll(s.substr(2), nullptr, 16));
                } else {
                    dst.set_int(key, std::stoll(s));
                }
            } catch (...) {
                dst.set_string(key, s);
            }
            break;
        }
        case ParamKind::IntList: {
            // If any entry exceeds INT_MAX (typical for hex bitmasks like
            // SFMT's `msk`), use uint_vec so generators that call
            // get_uint_vec() find it. Mirrors dict_to_params'
            // try-int-then-uint fallback in bindings.cpp.
            bool needs_uint = false;
            for (auto x : v.int_list_val) {
                if (x > static_cast<int64_t>(std::numeric_limits<int>::max())
                    || x < static_cast<int64_t>(std::numeric_limits<int>::min())) {
                    needs_uint = true; break;
                }
            }
            if (needs_uint) {
                std::vector<uint64_t> vu;
                vu.reserve(v.int_list_val.size());
                for (auto x : v.int_list_val) {
                    vu.push_back(static_cast<uint64_t>(x));
                }
                dst.set_uint_vec(key, vu);
            } else {
                std::vector<int> v32;
                v32.reserve(v.int_list_val.size());
                for (auto x : v.int_list_val) v32.push_back(static_cast<int>(x));
                dst.set_int_vec(key, v32);
            }
            break;
        }
    }
}

void parse_authors(const YAML::Node& n,
                   std::vector<Author>& out,
                   std::vector<std::string>& errors) {
    if (!n) {
        errors.emplace_back("authors list is required");
        return;
    }
    if (!n.IsSequence() || n.size() == 0) {
        errors.emplace_back("authors must be a non-empty list");
        return;
    }
    int i = -1;
    for (auto el : n) {
        ++i;
        if (el.IsScalar()) {
            std::string a = trim_copy(el.as<std::string>());
            auto pos = a.find_last_of(" \t");
            Author au;
            if (pos == std::string::npos) au.family = a;
            else { au.given = a.substr(0, pos); au.family = a.substr(pos + 1); }
            out.push_back(au);
            continue;
        }
        if (!el.IsMap()) {
            errors.push_back("authors[" + std::to_string(i)
                             + "] must be a string or mapping");
            continue;
        }
        Author au;
        if (el["family"]) au.family = trim_copy(el["family"].as<std::string>());
        if (el["given"]) au.given = trim_copy(el["given"].as<std::string>());
        if (au.family.empty()) {
            errors.push_back("authors[" + std::to_string(i)
                             + "].family is required");
            continue;
        }
        out.push_back(au);
    }
}

std::vector<Component> parse_components(const YAML::Node& n,
                                        std::vector<std::string>& errors) {
    std::vector<Component> out;
    if (!n.IsSequence() || n.size() == 0) {
        errors.emplace_back("components must be a non-empty list");
        return out;
    }
    int idx = -1;
    for (auto cn : n) {
        ++idx;
        if (!cn.IsMap()) {
            errors.push_back("components[" + std::to_string(idx)
                             + "] must be a mapping");
            continue;
        }
        Component c;
        if (cn["family"]) c.family = trim_copy(cn["family"].as<std::string>());
        if (cn["L"])      c.L      = cn["L"].as<int>(0);
        if (c.family.empty()) {
            errors.push_back("components[" + std::to_string(idx)
                             + "].family is required");
        }
        if (c.L <= 0) {
            errors.push_back("components[" + std::to_string(idx)
                             + "].L must be a positive int");
        }
        if (cn["params"]) {
            if (!cn["params"].IsMap()) {
                errors.push_back("components[" + std::to_string(idx)
                                 + "].params must be a mapping");
            } else {
                c.params = node_to_param_map(cn["params"]);
            }
        }
        if (cn["tempering"]) {
            if (!cn["tempering"].IsSequence()) {
                errors.push_back("components[" + std::to_string(idx)
                                 + "].tempering must be a list");
            } else {
                int tj = -1;
                for (auto tn : cn["tempering"]) {
                    ++tj;
                    if (!tn.IsMap() || !tn["type"]) {
                        errors.push_back("components[" + std::to_string(idx)
                                         + "].tempering[" + std::to_string(tj)
                                         + "] needs a 'type' key");
                        continue;
                    }
                    TemperingStep step;
                    step.type = tn["type"].as<std::string>();
                    for (auto kv : tn) {
                        auto k = kv.first.as<std::string>();
                        if (k == "type") continue;
                        step.params.emplace(k,
                            node_to_param_value(kv.second));
                    }
                    c.tempering.push_back(std::move(step));
                }
            }
        }
        out.push_back(std::move(c));
    }
    return out;
}

std::vector<CatalogGenerator>
parse_generators_node(const YAML::Node& n,
                      std::vector<std::string>& errors,
                      bool allow_empty) {
    std::vector<CatalogGenerator> out;
    if (!n || !n.IsSequence()) {
        if (allow_empty && !n) return out;
        errors.emplace_back("generators must be a non-empty list");
        return out;
    }
    if (n.size() == 0 && !allow_empty) {
        errors.emplace_back("generators must be a non-empty list");
        return out;
    }
    int i = -1;
    for (auto gn : n) {
        ++i;
        if (!gn.IsMap()) {
            errors.push_back("generators[" + std::to_string(i)
                             + "] must be a mapping");
            continue;
        }
        CatalogGenerator g;
        std::vector<std::string>& gerr = g.errors;
        g.id = trim_copy(gn["id"] ? gn["id"].as<std::string>() : "");
        if (!std::regex_match(g.id, kIdRe)) {
            gerr.push_back("id must match [a-z0-9][a-z0-9-]*, got '"
                           + g.id + "'");
        }
        g.display = gn["display"] ? trim_copy(gn["display"].as<std::string>())
                                   : g.id;
        if (g.display.empty()) g.display = g.id;
        g.family = gn["family"] ? trim_copy(gn["family"].as<std::string>()) : "";
        g.target = gn["target"] ? trim_copy(gn["target"].as<std::string>())
                                 : "tested_generator";
        g.combined = gn["combined"] ? gn["combined"].as<bool>() : false;
        g.Lmax = gn["Lmax"] ? gn["Lmax"].as<int>(0) : 0;
        g.starred = gn["starred"] ? gn["starred"].as<bool>() : false;
        g.notes_md = gn["notes"] ? gn["notes"].as<std::string>() : "";
        if (std::find(kKnownTargets.begin(), kKnownTargets.end(), g.target)
            == kKnownTargets.end()) {
            gerr.push_back("target must be one of {tested_generator, "
                           "primitive_generator}, got '" + g.target + "'");
        }
        if (g.Lmax <= 0) {
            gerr.push_back("Lmax must be a positive int, got "
                           + std::to_string(g.Lmax));
        }
        if (g.family.empty()) {
            gerr.push_back("family is required");
        }
        g.components = parse_components(gn["components"], gerr);
        if (g.target == "primitive_generator") {
            if (g.combined) {
                gerr.emplace_back("primitive_generator target cannot be combined");
            }
            if (g.components.size() != 1) {
                gerr.emplace_back("primitive_generator target must have "
                                  "exactly one component");
            } else if (!g.components[0].tempering.empty()) {
                gerr.emplace_back("primitive_generator target must have "
                                  "empty tempering");
            }
        }
        out.push_back(std::move(g));
    }
    return out;
}

void validate_pdf_exists(Paper& paper) {
    if (paper.pdf.empty()) return;
    std::filesystem::path src(paper.source_path);
    auto root = find_repo_root(src.parent_path());
    if (root.empty()) return;
    std::filesystem::path pdf = root / paper.pdf;
    std::error_code ec;
    if (!std::filesystem::is_regular_file(pdf, ec)) {
        paper.errors.push_back("pdf points to missing file: " + paper.pdf);
    }
}

int paper_default_Lmax(const Paper& paper) {
    for (const auto& g : paper.generators) {
        if (g.Lmax) return g.Lmax;
    }
    return 32;
}

void dry_build_generator(const Paper& paper, CatalogGenerator& g) {
    int default_Lmax = paper_default_Lmax(paper);
    (void)default_Lmax;
    int idx = -1;
    for (auto& comp : g.components) {
        ++idx;
        Params gp;
        for (const auto& kv : comp.params) apply_param(gp, kv.first, kv.second);
        std::unique_ptr<Generator> gen;
        try {
            gen = create_generator(comp.family, gp, comp.L);
        } catch (const std::exception& exc) {
            g.errors.push_back("components[" + std::to_string(idx)
                + "]: create_generator(" + comp.family + ", L="
                + std::to_string(comp.L) + ", ...) failed: " + exc.what());
            continue;
        }
        // Inherit `w` from the component's input params for tempering
        // steps that omit it (mirrors Python: Generator.params.get("w")).
        std::optional<int64_t> w_default;
        auto wit = comp.params.find("w");
        if (wit != comp.params.end() && wit->second.kind == ParamKind::Int) {
            w_default = wit->second.int_val;
        }
        int tj = -1;
        for (const auto& step : comp.tempering) {
            ++tj;
            Params tp;
            bool has_w = false;
            for (const auto& kv : step.params) {
                if (kv.first == "w") has_w = true;
                apply_param(tp, kv.first, kv.second);
            }
            if (!has_w && w_default.has_value()) {
                tp.set_int("w", *w_default);
            }
            try {
                auto t = create_transformation(step.type, tp);
                (void)t;
            } catch (const std::exception& exc) {
                g.errors.push_back("components[" + std::to_string(idx)
                    + "].tempering[" + std::to_string(tj) + "] ("
                    + step.type + "): " + exc.what());
            }
        }
    }
}

// FNV-1a 64-bit, returned as 16-char zero-padded hex.
std::string fnv1a_hex(const std::string& canonical) {
    uint64_t h = 14695981039346656037ULL;
    for (unsigned char c : canonical) {
        h ^= static_cast<uint64_t>(c);
        h *= 1099511628211ULL;
    }
    char buf[17];
    std::snprintf(buf, sizeof(buf), "%016lx",
                  static_cast<unsigned long>(h));
    return std::string(buf);
}

void write_canonical_pv(std::ostringstream& os, const ParamValue& v);

void write_canonical_pmap(std::ostringstream& os, const ParamMap& m) {
    os << '{';
    bool first = true;
    for (const auto& kv : m) {        // ordered map → deterministic
        if (!first) os << ',';
        first = false;
        os << '"' << kv.first << "\":";
        write_canonical_pv(os, kv.second);
    }
    os << '}';
}

void write_canonical_pv(std::ostringstream& os, const ParamValue& v) {
    switch (v.kind) {
        case ParamKind::Int:    os << v.int_val; break;
        case ParamKind::Bool:   os << (v.bool_val ? "true" : "false"); break;
        case ParamKind::String: os << '"' << v.string_val << '"'; break;
        case ParamKind::IntList:
            os << '[';
            for (size_t i = 0; i < v.int_list_val.size(); ++i) {
                if (i) os << ',';
                os << v.int_list_val[i];
            }
            os << ']';
            break;
    }
}

}  // namespace

// ── Paper presentation helpers ────────────────────────────────────────

std::string Paper::author_list_short() const {
    if (authors.empty()) return "";
    if (authors.size() == 1) return authors[0].short_name();
    if (authors.size() == 2) {
        return authors[0].short_name() + " and " + authors[1].short_name();
    }
    return authors[0].short_name() + " et al.";
}

std::string Paper::display() const {
    std::string who;
    if (authors.empty())          who = id;
    else if (authors.size() == 1) who = authors[0].short_name();
    else if (authors.size() == 2) who = authors[0].short_name()
                                       + " & " + authors[1].short_name();
    else                          who = authors[0].short_name() + " et al.";
    return who + " " + std::to_string(year);
}

namespace {

std::string format_authors_acm(const std::vector<Author>& authors) {
    auto one = [](const Author& a) -> std::string {
        if (!a.given.empty()) {
            std::string initials;
            std::stringstream ss(a.given);
            std::string tok;
            while (ss >> tok) {
                if (!tok.empty()) {
                    initials += tok[0];
                    initials += '.';
                }
            }
            return a.family + ", " + initials;
        }
        return a.family;
    };
    std::vector<std::string> formatted;
    for (const auto& a : authors) formatted.push_back(one(a));
    if (formatted.size() == 1) return formatted[0];
    if (formatted.size() == 2) return formatted[0] + " and " + formatted[1];
    std::string out;
    for (size_t i = 0; i + 1 < formatted.size(); ++i) {
        if (i) out += ", ";
        out += formatted[i];
    }
    out += ", and " + formatted.back();
    return out;
}

std::string rstrip_dot(std::string s) {
    while (!s.empty() && s.back() == '.') s.pop_back();
    return s;
}

}  // namespace

std::string Paper::acmtrans_citation() const {
    std::vector<std::string> parts;
    if (!authors.empty()) {
        parts.push_back(rstrip_dot(format_authors_acm(authors)) + ".");
    }
    parts.push_back(std::to_string(year) + ".");
    if (!title.empty()) parts.push_back(rstrip_dot(title) + ".");
    std::string venue_str = venue;
    if (!volume.empty()) {
        if (!venue_str.empty()) venue_str += " ";
        venue_str += volume;
        if (!issue.empty()) venue_str += ", " + issue;
    }
    std::string tail = "(" + std::to_string(year) + ")";
    if (!pages.empty()) tail += ", " + pages;
    tail += ".";
    if (!venue_str.empty()) parts.push_back(venue_str + " " + tail);
    if (!doi.empty()) parts.push_back(doi);
    std::string out;
    for (const auto& p : parts) {
        if (p.empty()) continue;
        if (!out.empty()) out += " ";
        out += p;
    }
    return out;
}

// ── Catalog impl ──────────────────────────────────────────────────────

Catalog::Catalog(std::string library_dir)
    : library_dir_(std::move(library_dir)) {}

bool is_cross_check_yaml(const std::string& filename) {
    static const std::string suffix = "_params.yaml";
    if (filename.size() < suffix.size()) return false;
    return std::equal(suffix.rbegin(), suffix.rend(), filename.rbegin());
}

void Catalog::load() {
    papers_.clear();
    generator_index_.clear();
    if (library_dir_.empty()) return;
    std::error_code ec;
    if (!std::filesystem::is_directory(library_dir_, ec)) return;

    std::vector<std::string> paths;
    for (auto& entry : std::filesystem::directory_iterator(library_dir_)) {
        auto name = entry.path().filename().string();
        if (name.size() < 5
            || name.substr(name.size() - 5) != ".yaml") continue;
        if (is_cross_check_yaml(name)) continue;
        paths.push_back(entry.path().string());
    }
    std::sort(paths.begin(), paths.end());
    for (const auto& path : paths) {
        Paper p;
        try {
            p = parse_paper(path);
        } catch (const std::exception& exc) {
            std::cerr << "library: failed to parse " << path
                      << ": " << exc.what() << "\n";
            continue;
        }
        insert(std::move(p));
    }
}

void Catalog::reload_if_stale() {
    if (library_dir_.empty()) {
        papers_.clear();
        generator_index_.clear();
        return;
    }
    std::error_code ec;
    if (!std::filesystem::is_directory(library_dir_, ec)) {
        papers_.clear();
        generator_index_.clear();
        return;
    }
    std::map<std::string, std::string> on_disk;     // stem -> path
    for (auto& entry : std::filesystem::directory_iterator(library_dir_)) {
        auto name = entry.path().filename().string();
        if (name.size() < 5
            || name.substr(name.size() - 5) != ".yaml") continue;
        if (is_cross_check_yaml(name)) continue;
        auto stem = name.substr(0, name.size() - 5);
        on_disk.emplace(stem, entry.path().string());
    }
    // Drop papers whose files vanished.
    std::vector<std::string> dropped;
    for (const auto& kv : papers_) {
        if (!on_disk.count(kv.first)) dropped.push_back(kv.first);
    }
    for (const auto& id : dropped) drop(id);

    // Reload mtime-newer files.
    for (const auto& kv : on_disk) {
        const auto& stem = kv.first;
        const auto& path = kv.second;
        double mtime = stat_mtime(path);
        auto it = papers_.find(stem);
        if (it != papers_.end() && it->second.source_mtime == mtime) continue;
        Paper p;
        try {
            p = parse_paper(path);
        } catch (const std::exception& exc) {
            std::cerr << "library: reload " << path << " failed: "
                      << exc.what() << "\n";
            continue;
        }
        drop(stem);
        insert(std::move(p));
    }
}

void Catalog::insert(Paper paper) {
    if (papers_.count(paper.id)) {
        std::cerr << "library: duplicate paper id '" << paper.id
                  << "' (keeping first)\n";
        return;
    }
    auto pid = paper.id;
    auto [it, _] = papers_.emplace(pid, std::move(paper));
    Paper& stored = it->second;
    for (size_t i = 0; i < stored.generators.size(); ++i) {
        const auto& g = stored.generators[i];
        if (generator_index_.count(g.id)) {
            stored.errors.push_back("generator id '" + g.id
                + "' already exists in paper '"
                + generator_index_[g.id].first + "'");
            continue;
        }
        generator_index_[g.id] = {pid, static_cast<int>(i)};
    }
}

void Catalog::drop(const std::string& paper_id) {
    auto it = papers_.find(paper_id);
    if (it == papers_.end()) return;
    for (const auto& g : it->second.generators) {
        auto gi = generator_index_.find(g.id);
        if (gi != generator_index_.end() && gi->second.first == paper_id) {
            generator_index_.erase(gi);
        }
    }
    papers_.erase(it);
}

std::vector<Paper>
Catalog::papers(const PapersFilter& filter) const {
    std::vector<Paper> out;
    for (const auto& kv : papers_) {
        const Paper& p = kv.second;
        if (!filter.include_invalid && !p.valid()) continue;
        if (filter.starred.has_value() && p.starred != *filter.starred) continue;
        if (filter.tag.has_value()) {
            const auto& tag = *filter.tag;
            if (std::find(p.tags.begin(), p.tags.end(), tag) == p.tags.end()) {
                continue;
            }
        }
        out.push_back(p);
    }
    std::sort(out.begin(), out.end(), [](const Paper& a, const Paper& b) {
        if (a.year != b.year) return a.year < b.year;
        return a.id < b.id;
    });
    return out;
}

std::optional<Paper>
Catalog::paper(const std::string& paper_id) const {
    auto it = papers_.find(paper_id);
    if (it == papers_.end()) return std::nullopt;
    return it->second;
}

std::optional<std::pair<Paper, CatalogGenerator>>
Catalog::generator(const std::string& gen_id) const {
    auto it = generator_index_.find(gen_id);
    if (it == generator_index_.end()) return std::nullopt;
    auto [pid, idx] = it->second;
    auto pit = papers_.find(pid);
    if (pit == papers_.end()) return std::nullopt;
    if (idx >= static_cast<int>(pit->second.generators.size())) return std::nullopt;
    return std::make_pair(pit->second, pit->second.generators[idx]);
}

std::vector<std::pair<Paper, CatalogGenerator>>
Catalog::all_generators(const std::optional<std::string>& family) const {
    std::vector<Paper> sorted = papers({/*starred=*/{}, /*tag=*/{},
                                        /*include_invalid=*/true});
    std::vector<std::pair<Paper, CatalogGenerator>> out;
    for (const auto& p : sorted) {
        for (const auto& g : p.generators) {
            if (family.has_value() && g.family != *family) continue;
            out.emplace_back(p, g);
        }
    }
    return out;
}

// ── Free helpers ──────────────────────────────────────────────────────

std::string config_hash(const std::string& family,
                        const ParamMap& params,
                        const std::vector<TemperingStep>& tempering) {
    std::ostringstream os;
    os << "{\"family\":\"" << family << "\",\"params\":";
    write_canonical_pmap(os, params);
    os << ",\"tempering\":[";
    for (size_t i = 0; i < tempering.size(); ++i) {
        if (i) os << ',';
        os << "{\"type\":\"" << tempering[i].type << "\"";
        for (const auto& kv : tempering[i].params) {
            os << ",\"" << kv.first << "\":";
            write_canonical_pv(os, kv.second);
        }
        os << '}';
    }
    os << "]}";
    return fnv1a_hex(os.str());
}

Paper parse_paper(const std::string& path) {
    Paper paper;
    paper.source_path = path;
    paper.source_mtime = stat_mtime(path);

    YAML::Node raw;
    try {
        raw = YAML::LoadFile(path);
    } catch (const std::exception& exc) {
        throw std::runtime_error(
            std::string("YAML load failed: ") + exc.what());
    }
    if (!raw || !raw.IsMap()) {
        throw std::runtime_error("top-level YAML must be a mapping");
    }

    auto& errors = paper.errors;
    paper.id = trim_copy(raw["id"] ? raw["id"].as<std::string>("") : "");
    if (!std::regex_match(paper.id, kIdRe)) {
        errors.push_back("id must match [a-z0-9][a-z0-9-]*, got '"
                         + paper.id + "'");
    } else {
        std::filesystem::path p(path);
        std::string stem = p.stem().string();
        if (paper.id != stem) {
            errors.push_back("id '" + paper.id
                + "' must equal filename stem '" + stem + "'");
        }
    }

    parse_authors(raw["authors"], paper.authors, errors);

    paper.year = raw["year"] && raw["year"].IsScalar()
        ? raw["year"].as<int>(0) : 0;
    if (paper.year == 0
        && (!raw["year"] || raw["year"].as<std::string>("") != "0")) {
        errors.emplace_back("year must be an integer");
    }

    paper.title = trim_copy(raw["title"] ? raw["title"].as<std::string>("") : "");
    paper.venue = trim_copy(raw["venue"] ? raw["venue"].as<std::string>("") : "");
    if (paper.title.empty()) errors.emplace_back("title is required");
    if (paper.venue.empty()) errors.emplace_back("venue is required");

    paper.doi = trim_copy(raw["doi"] ? raw["doi"].as<std::string>("") : "");
    if (paper.doi.empty()) {
        errors.emplace_back("doi is required");
    } else if (!(paper.doi.rfind("https://doi.org/", 0) == 0
              || paper.doi.rfind("https://www.jstor.org/", 0) == 0
              || paper.doi.rfind("http://www.jstor.org/", 0) == 0)) {
        errors.emplace_back("doi must be a full https URL "
                            "(doi.org or jstor.org stable link)");
    }

    paper.pdf    = trim_copy(raw["pdf"]    ? raw["pdf"].as<std::string>("")    : "");
    paper.bibkey = trim_copy(raw["bibkey"] ? raw["bibkey"].as<std::string>("") : "");
    paper.volume = trim_copy(raw["volume"] ? raw["volume"].as<std::string>("") : "");
    paper.issue  = trim_copy(raw["issue"]  ? raw["issue"].as<std::string>("")  : "");
    paper.pages  = trim_copy(raw["pages"]  ? raw["pages"].as<std::string>("")  : "");
    paper.abstract_md = raw["abstract"] ? raw["abstract"].as<std::string>("") : "";
    paper.notes_md    = raw["notes"]    ? raw["notes"].as<std::string>("")    : "";
    if (raw["tags"] && raw["tags"].IsSequence()) {
        for (auto t : raw["tags"]) paper.tags.push_back(t.as<std::string>(""));
    }
    paper.starred  = raw["starred"]  ? raw["starred"].as<bool>(false)  : false;
    paper.deferred = raw["deferred"] ? raw["deferred"].as<bool>(false) : false;

    paper.generators = parse_generators_node(raw["generators"], errors,
                                             /*allow_empty=*/paper.deferred);

    if (errors.empty()) {
        validate_pdf_exists(paper);
        for (auto& g : paper.generators) {
            if (g.valid()) dry_build_generator(paper, g);
        }
    }

    return paper;
}

}  // namespace regpoly_catalog
