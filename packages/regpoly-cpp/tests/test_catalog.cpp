// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

// Phase 3.2 (TDD): paper-centric Catalog in C++.
//
// Loads the real docs/library/*.yaml fixtures and asserts the
// invariants that test_library.py covers on the Python side.

#include <gtest/gtest.h>

#include <chrono>
#include <filesystem>
#include <fstream>
#include <string>
#include <thread>

#include "catalog.h"

namespace fs = std::filesystem;
using namespace regpoly_catalog;

namespace {

// Walk up from this source file's compile-time __FILE__ to find the
// repo's docs/library directory. Falls back to walking up from CWD
// if that fails (depending on where ctest was invoked from).
fs::path find_library_dir() {
    auto try_paths = std::vector<fs::path>{
        fs::path(__FILE__).parent_path(),
        fs::current_path(),
    };
    for (auto base : try_paths) {
        for (int up = 0; up < 8; ++up) {
            auto cand = base / "docs" / "library";
            if (fs::is_directory(cand)) return cand;
            if (base == base.parent_path()) break;
            base = base.parent_path();
        }
    }
    ADD_FAILURE() << "could not find docs/library directory";
    return {};
}

}  // namespace

TEST(Catalog, LoadsRealLibraryDirectory) {
    Catalog c(find_library_dir().string());
    c.load();

    Catalog::PapersFilter f;
    f.include_invalid = true;
    auto ps = c.papers(f);
    std::vector<std::string> ids;
    for (const auto& p : ps) ids.push_back(p.id);

    auto has = [&](const std::string& s) {
        return std::find(ids.begin(), ids.end(), s) != ids.end();
    };
    EXPECT_TRUE(has("matsumoto-nishimura-1998")) << "MT paper missing";
    EXPECT_TRUE(has("lecuyer-1996"))             << "lecuyer-1996 missing";
    EXPECT_TRUE(has("lecuyer-1999"))             << "lecuyer-1999 missing";
}

TEST(Catalog, AllPapersValidate) {
    Catalog c(find_library_dir().string());
    c.load();

    Catalog::PapersFilter f;
    f.include_invalid = true;
    auto ps = c.papers(f);
    for (const auto& p : ps) {
        if (!p.valid()) {
            std::string msg = p.id + ":";
            for (const auto& e : p.errors) msg += " [" + e + "]";
            for (const auto& g : p.generators) {
                if (!g.valid()) {
                    msg += "\n    " + g.id + ":";
                    for (const auto& ge : g.errors) msg += " [" + ge + "]";
                }
            }
            ADD_FAILURE() << msg;
        }
    }
}

TEST(Catalog, GeneratorLookupIsUnique) {
    Catalog c(find_library_dir().string());
    c.load();
    for (const auto& gid : {"mt19937", "taus88", "lfsr113"}) {
        auto loc = c.generator(gid);
        ASSERT_TRUE(loc.has_value()) << gid << " missing";
        EXPECT_EQ(loc->second.id, gid);
    }
}

TEST(Catalog, PaperDisplayFormat) {
    Catalog c(find_library_dir().string());
    c.load();
    auto p = c.paper("matsumoto-nishimura-1998");
    ASSERT_TRUE(p.has_value());
    EXPECT_EQ(p->display(), "Matsumoto & Nishimura 1998");
    auto q = c.paper("lecuyer-1999");
    ASSERT_TRUE(q.has_value());
    EXPECT_EQ(q->display(), "L'Ecuyer 1999");
}

TEST(Catalog, AcmTransCitationCarriesAllFields) {
    Catalog c(find_library_dir().string());
    c.load();
    auto p = c.paper("lecuyer-1999");
    ASSERT_TRUE(p.has_value());
    auto cite = p->acmtrans_citation();
    EXPECT_NE(cite.find("L'Ecuyer"), std::string::npos);
    EXPECT_NE(cite.find("1999"), std::string::npos);
    EXPECT_NE(cite.find("Tables of Maximally"), std::string::npos);
    EXPECT_NE(cite.find("Mathematics of Computation"), std::string::npos);
    EXPECT_NE(cite.find("261–269"), std::string::npos);
    EXPECT_NE(cite.find("https://www.jstor.org/stable/2585109"),
              std::string::npos);
}

TEST(ConfigHash, StableUnderKeyReorder) {
    ParamMap a;
    a.emplace("k", ParamValue::make_int(31));
    a.emplace("nb_terms", ParamValue::make_int(3));
    a.emplace("quicktaus", ParamValue::make_bool(true));
    a.emplace("poly", ParamValue::make_int_list({0, 6, 31}));
    a.emplace("s", ParamValue::make_int(18));

    // Insert in scrambled order — std::map is ordered by key, so the
    // canonical serialization is identical regardless of insert order.
    ParamMap b;
    b.emplace("s", ParamValue::make_int(18));
    b.emplace("poly", ParamValue::make_int_list({0, 6, 31}));
    b.emplace("nb_terms", ParamValue::make_int(3));
    b.emplace("k", ParamValue::make_int(31));
    b.emplace("quicktaus", ParamValue::make_bool(true));

    EXPECT_EQ(config_hash("Tausworthe", a, {}),
              config_hash("Tausworthe", b, {}));
}

TEST(Catalog, ReloadIfStaleReReadsTouchedFiles) {
    auto src = find_library_dir() / "matsumoto-nishimura-1998.yaml";
    fs::path tmpdir = fs::temp_directory_path() /
        ("regpoly_catalog_reload_" + std::to_string(::getpid()));
    fs::create_directories(tmpdir);
    auto tgt = tmpdir / "matsumoto-nishimura-1998.yaml";

    {
        std::ifstream in(src, std::ios::binary);
        std::ofstream out(tgt, std::ios::binary);
        out << in.rdbuf();
    }

    Catalog c(tmpdir.string());
    c.load();

    auto before_p = c.paper("matsumoto-nishimura-1998");
    ASSERT_TRUE(before_p.has_value());
    double before = before_p->source_mtime;

    std::this_thread::sleep_for(std::chrono::milliseconds(20));

    // Touch the file to bump its mtime.
    auto now = fs::file_time_type::clock::now();
    fs::last_write_time(tgt, now);

    c.reload_if_stale();
    auto after_p = c.paper("matsumoto-nishimura-1998");
    ASSERT_TRUE(after_p.has_value());
    EXPECT_GE(after_p->source_mtime, before);

    fs::remove_all(tmpdir);
}
