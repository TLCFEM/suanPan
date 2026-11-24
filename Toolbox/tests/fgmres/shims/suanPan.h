// Minimal stand-in for the application header used by the standalone FGMRES
// test: Armadillo types plus the two logging macros the wrapper needs.
#ifndef TEST_FGMRES_SUANPAN_H
#define TEST_FGMRES_SUANPAN_H

// Keep the test free of BLAS/LAPACK link dependencies.
#define ARMA_DONT_USE_BLAS
#define ARMA_DONT_USE_LAPACK

#include <armadillo/arma>

using namespace arma;

// The production macros use fmt-style "{}" placeholders; test builds simply
// discard the messages (the assertions check return codes, not log text).
#define suanpan_error(...) ((void)0)
#define suanpan_debug(...) ((void)0)

#endif
