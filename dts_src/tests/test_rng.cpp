#include <set>
#include <omp.h>
#include "check.h"

void test_stream_zero_is_seed() {
    CHECK(RNG::StreamSeed(1234, 0) == 1234);
}

void test_streams_distinct() {
    std::set<int> seeds;
    for (int t = 0; t < 64; t++) {
        int s = RNG::StreamSeed(1234, t);
        CHECK(s >= 0);
        seeds.insert(s);
    }
    CHECK(seeds.size() == 64);
}

void test_thread_rngs_distinct() {
    State* pState = make_test_state("in.dts", "top.top", 4);
    RNG* rngs[4] = {nullptr, nullptr, nullptr, nullptr};
    #pragma omp parallel num_threads(4)
    {
        rngs[omp_get_thread_num()] = pState->GetRandomNumberGenerator();
    }
    std::set<RNG*> distinct(rngs, rngs + 4);
    CHECK(distinct.size() == 4);
    CHECK(rngs[0] == pState->GetRandomNumberGenerator());
    delete pState;
}

int main() {
    test_stream_zero_is_seed();
    test_streams_distinct();
    test_thread_rngs_distinct();
    return check_report();
}
