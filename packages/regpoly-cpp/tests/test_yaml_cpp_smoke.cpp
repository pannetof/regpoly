// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 Francois Panneton, Ph.D.

// Phase 3.1: yaml-cpp dependency smoke test. Just verifies the
// vendored library is reachable, parses a small document correctly,
// and round-trips a value the way the catalog code will need.

#include <gtest/gtest.h>

#include <yaml-cpp/yaml.h>

#include <string>

TEST(YamlCppSmoke, ParsesScalarsAndSequences) {
    const std::string src =
        "title: Mersenne Twister\n"
        "year: 1998\n"
        "authors:\n"
        "  - Matsumoto\n"
        "  - Nishimura\n"
        "tags: [prng, gf2]\n";

    YAML::Node doc = YAML::Load(src);

    ASSERT_TRUE(doc.IsMap());
    EXPECT_EQ(doc["title"].as<std::string>(), "Mersenne Twister");
    EXPECT_EQ(doc["year"].as<int>(), 1998);

    ASSERT_TRUE(doc["authors"].IsSequence());
    ASSERT_EQ(doc["authors"].size(), 2u);
    EXPECT_EQ(doc["authors"][0].as<std::string>(), "Matsumoto");
    EXPECT_EQ(doc["authors"][1].as<std::string>(), "Nishimura");

    ASSERT_TRUE(doc["tags"].IsSequence());
    ASSERT_EQ(doc["tags"].size(), 2u);
    EXPECT_EQ(doc["tags"][0].as<std::string>(), "prng");
}

TEST(YamlCppSmoke, ParsesNestedMaps) {
    const std::string src =
        "params:\n"
        "  w: 32\n"
        "  r: 624\n"
        "  matrix_a: 0x9908b0df\n"
        "tempering:\n"
        "  - type: tempMK\n"
        "    b: 0x9d2c5680\n";

    YAML::Node doc = YAML::Load(src);

    ASSERT_TRUE(doc["params"].IsMap());
    EXPECT_EQ(doc["params"]["w"].as<int>(), 32);
    EXPECT_EQ(doc["params"]["r"].as<int>(), 624);

    // YAML 1.2 hex-int support — used pervasively by the catalog
    // (tempering masks, MT matrix_a, etc.). yaml-cpp parses 0x... as
    // a string, so the catalog code must convert; document that here.
    auto matrix_a_str = doc["params"]["matrix_a"].as<std::string>();
    EXPECT_EQ(matrix_a_str, "0x9908b0df");

    ASSERT_TRUE(doc["tempering"].IsSequence());
    ASSERT_EQ(doc["tempering"].size(), 1u);
    EXPECT_EQ(doc["tempering"][0]["type"].as<std::string>(), "tempMK");
}

TEST(YamlCppSmoke, EmitsRoundTrippableDocs) {
    YAML::Node out;
    out["paper_id"] = "matsumoto-nishimura-1998";
    out["generators"].push_back("mt19937");
    out["generators"].push_back("mt19937-64");

    YAML::Emitter em;
    em << out;
    std::string yaml = em.c_str();

    YAML::Node round = YAML::Load(yaml);
    EXPECT_EQ(round["paper_id"].as<std::string>(),
              "matsumoto-nishimura-1998");
    ASSERT_EQ(round["generators"].size(), 2u);
    EXPECT_EQ(round["generators"][0].as<std::string>(), "mt19937");
}
