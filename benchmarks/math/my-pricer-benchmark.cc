//
// Created by peter on 1/25/20.
//


#include <benchmark/benchmark.h>
#include <random>
#include <iterator>
#include <iostream>
#include <functional>

#include <pricer-base.h>
#include <my_pricer.h>


typedef double *__restrict__ __attribute__((aligned(ALIGN_TO))) Real_Ptr;

static void BM_My_Pricer(benchmark::State &state) {


    Real_Ptr ss = (Real_Ptr) aligned_alloc(ALIGN_TO, ((sizeof(FLOAT) * state.range(0)) / ALIGN_TO + 1) * ALIGN_TO);
    Real_Ptr xs = (Real_Ptr) aligned_alloc(ALIGN_TO, ((sizeof(FLOAT) * state.range(0)) / ALIGN_TO + 1) * ALIGN_TO);
    Real_Ptr rs = (Real_Ptr) aligned_alloc(ALIGN_TO, ((sizeof(FLOAT) * state.range(0)) / ALIGN_TO + 1) * ALIGN_TO);
    Real_Ptr sigmas = (Real_Ptr) aligned_alloc(ALIGN_TO, ((sizeof(FLOAT) * state.range(0)) / ALIGN_TO + 1) * ALIGN_TO);
    Real_Ptr ts = (Real_Ptr) aligned_alloc(ALIGN_TO, ((sizeof(FLOAT) * state.range(0)) / ALIGN_TO + 1) * ALIGN_TO);
    Real_Ptr taus = (Real_Ptr) aligned_alloc(ALIGN_TO, ((sizeof(FLOAT) * state.range(0)) / ALIGN_TO + 1) * ALIGN_TO);
    Real_Ptr prices = (Real_Ptr) aligned_alloc(ALIGN_TO, ((sizeof(FLOAT) * state.range(0)) / ALIGN_TO + 1) * ALIGN_TO);
    Uint64_Ptr flags = (Uint64_Ptr) aligned_alloc(ALIGN_TO,
                                                  ((sizeof(UINT64) * state.range(0)) / ALIGN_TO + 1) * ALIGN_TO);

    Real_Ptr tmp1 = (Real_Ptr) aligned_alloc(ALIGN_TO, ((sizeof(FLOAT) * state.range(0)) / ALIGN_TO + 1) * ALIGN_TO);
    Real_Ptr tmp2 = (Real_Ptr) aligned_alloc(ALIGN_TO, ((sizeof(FLOAT) * state.range(0)) / ALIGN_TO + 1) * ALIGN_TO);
    Real_Ptr tmp3 = (Real_Ptr) aligned_alloc(ALIGN_TO, ((sizeof(FLOAT) * state.range(0)) / ALIGN_TO + 1) * ALIGN_TO);
    Real_Ptr tmp4 = (Real_Ptr) aligned_alloc(ALIGN_TO, ((sizeof(FLOAT) * state.range(0)) / ALIGN_TO + 1) * ALIGN_TO);
    Real_Ptr tmp5 = (Real_Ptr) aligned_alloc(ALIGN_TO, ((sizeof(FLOAT) * state.range(0)) / ALIGN_TO + 1) * ALIGN_TO);

    Real_Ptr sigmaA = (Real_Ptr) aligned_alloc(ALIGN_TO, ((sizeof(FLOAT) * state.range(0)) / ALIGN_TO + 1) * ALIGN_TO);
    Real_Ptr sigmaA2T2 = (Real_Ptr) aligned_alloc(ALIGN_TO,
                                                  ((sizeof(FLOAT) * state.range(0)) / ALIGN_TO + 1) * ALIGN_TO);
    Real_Ptr sigmaAsqrtT = (Real_Ptr) aligned_alloc(ALIGN_TO,
                                                    ((sizeof(FLOAT) * state.range(0)) / ALIGN_TO + 1) * ALIGN_TO);
    Real_Ptr emrt = (Real_Ptr) aligned_alloc(ALIGN_TO, ((sizeof(FLOAT) * state.range(0)) / ALIGN_TO + 1) * ALIGN_TO);
    Real_Ptr d1 = (Real_Ptr) aligned_alloc(ALIGN_TO, ((sizeof(FLOAT) * state.range(0)) / ALIGN_TO + 1) * ALIGN_TO);
    Real_Ptr d2 = (Real_Ptr) aligned_alloc(ALIGN_TO, ((sizeof(FLOAT) * state.range(0)) / ALIGN_TO + 1) * ALIGN_TO);

    Real_Ptr d2dx2_prep = (Real_Ptr) aligned_alloc(ALIGN_TO,
                                                   ((sizeof(FLOAT) * state.range(0)) / ALIGN_TO + 1) * ALIGN_TO);


    // First create an instance of an engine.
    std::random_device rnd_device;

    // Specify the engine and distribution.
    std::mt19937 mersenne_engine{rnd_device()};  // Generates random integers

    std::uniform_real_distribution<double> ss_dist{1, 200};
    std::uniform_real_distribution<double> rs_dist{0.01, 1.0};
    std::uniform_real_distribution<double> sigmas_dist{0.01, 1.0};
    std::uniform_real_distribution<double> ts_dist{0.09, 2.0};
    std::uniform_real_distribution<double> taus_dist{0.08, 0.09};
    std::uniform_int_distribution<UINT64> flags_dist{0, 3};

    auto ss_gen = [&ss_dist, &mersenne_engine]() {
        return ss_dist(mersenne_engine);
    };

    auto rs_gen = [&rs_dist, &mersenne_engine]() {
        return rs_dist(mersenne_engine);
    };

    auto sigmas_gen = [&sigmas_dist, &mersenne_engine]() {
        return sigmas_dist(mersenne_engine);
    };

    auto ts_gen = [&ts_dist, &mersenne_engine]() {
        return ts_dist(mersenne_engine);
    };

    auto taus_gen = [&taus_dist, &mersenne_engine]() {
        return taus_dist(mersenne_engine);
    };

    auto flags_gen = [&flags_dist, &mersenne_engine]() {
        return flags_dist(mersenne_engine);
    };


    for (auto _ : state) {

        state.PauseTiming();

        prices[0] = prices[0] / 2;
        std::generate(ss, &ss[state.range(0)], ss_gen);
        std::generate(xs, &xs[state.range(0)], ss_gen);
        std::generate(rs, &rs[state.range(0)], rs_gen);
        std::generate(sigmas, &sigmas[state.range(0)], sigmas_gen);
        std::generate(ts, &ts[state.range(0)], ts_gen);
        std::generate(taus, &taus[state.range(0)], taus_gen);
        std::generate(flags, &flags[state.range(0)], flags_gen);

        for (UINT64 i = 0; i < (UINT64) state.range(0); ++i) {
            ts[i] += taus[i];
        }

        init_my_pricer();

        prepare_my_pricer(
                state.range(0),
                ss,
                sigmas,
                ts,
                taus,
                rs,
                tmp1,
                tmp2,
                tmp3,
                tmp4,
                tmp5,
                sigmaA,
                sigmaA2T2,
                sigmaAsqrtT,
                emrt,
                d2dx2_prep);

        state.ResumeTiming();

        my_pricer(
                state.range(0),
                flags,
                ss,
                xs,
                sigmaA2T2,
                sigmaAsqrtT,
                emrt,
                tmp1,
                tmp2,
                tmp3,
                tmp4,
                d1,
                d2,
                prices);

    }

    free(ss);
    free(xs);
    free(rs);
    free(sigmas);
    free(ts);
    free(taus);
    free(prices);
    free(flags);
    free(tmp1);
    free(tmp2);
    free(tmp3);
    free(tmp4);
    free(tmp5);
    free(sigmaA);
    free(sigmaA2T2);
    free(sigmaAsqrtT);
    free(emrt);
    free(d1);
    free(d2);
    free(d2dx2_prep);


}
// Register the function as a benchmark
BENCHMARK(BM_My_Pricer)->Arg(1 << 20)->Arg(1 << 21)->Arg(1 << 22);
