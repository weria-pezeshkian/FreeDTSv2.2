#pragma once
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include "../State.h"

static int g_check_failures = 0;

#define CHECK(cond) \
    do { if (!(cond)) { ++g_check_failures; \
        std::cerr << __FILE__ << ":" << __LINE__ << ": CHECK failed: " #cond "\n"; } } while (0)

#define CHECK_NEAR(a, b, tol) \
    do { if (!(std::fabs((a) - (b)) <= (tol))) { ++g_check_failures; \
        std::cerr << __FILE__ << ":" << __LINE__ << ": CHECK_NEAR failed: " #a " vs " #b "\n"; } } while (0)

inline int check_report() {
    if (g_check_failures) std::cerr << g_check_failures << " check(s) failed\n";
    else std::cout << "all checks passed\n";
    return g_check_failures ? 1 : 0;
}

// Input files are looked up in the cwd; run_tests.sh copies tests/data there.
inline State* make_test_state(const char* input, const char* top, int nthreads) {
    std::vector<std::string> args = {"DTS", "-in", input, "-top", top,
                                     "-nt", std::to_string(nthreads), "-seed", "4321"};
    State* pState = new State(args);
    pState->Initialize();
    return pState;
}
