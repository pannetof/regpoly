// regpoly/regpoly.h — public API umbrella header.
//
// Downstream C++ users include this single header to access the
// stable public surface of the regpoly library:
//
//     #include <regpoly/regpoly.h>
//
// then link against the regpoly::regpoly target via find_package.
//
// Phase 1 publishes a small, curated set: Generator and its abstract
// interface, CombinedGenerator, BitVect, BitMatrix, the Transformation
// interface, parameter-spec types, and the lattice-method entry points
// (test_me_lat, test_me_harase, test_me_notprimitive, plus their
// SIMD-aware variant; the StackBase caches TemperOptCache and
// PISCache).
//
// Internal headers (algebra/, lattice/ implementation details,
// transforms/ details) stay private under regpoly/internal/ and are
// not part of the API contract — they may move or change between
// versions without notice.

#pragma once

#include <regpoly/bitvect.h>
#include <regpoly/combined.h>
#include <regpoly/equidistribution_runner.h>
#include <regpoly/gauss.h>
#include <regpoly/generator.h>
#include <regpoly/me_harase.h>
#include <regpoly/me_helpers.h>
#include <regpoly/me_notprimitive.h>
#include <regpoly/me_notprimitive_simd.h>
#include <regpoly/params.h>
#include <regpoly/resolution_sets.h>
#include <regpoly/temper_optimizer.h>
#include <regpoly/transformation.h>
