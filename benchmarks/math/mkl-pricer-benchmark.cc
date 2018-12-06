//
// Created by peter on 12/3/18.
//

#include <benchmark/benchmark.h>
#include <vector>
#include <random>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <functional>
#include <cmath>


#include <src/math/pricers/mkl_pricer.h>
#include <mkl.h>


static void BM_Pricer_MKL(benchmark::State &state) {


    std::vector<double> ss(state.range(0));
    std::vector<double> xs(state.range(0));
    std::vector<double> rs(state.range(0));
    std::vector<double> sigmas(state.range(0));
    std::vector<double> ts(state.range(0));
    std::vector<double> taus(state.range(0));
    std::vector<double> prices(state.range(0));
    std::vector<MKL_INT64> flags(state.range(0));

    std::vector<double> tmp1(state.range(0));
    std::vector<double> tmp2(state.range(0));
    std::vector<double> tmp3(state.range(0));
    std::vector<double> tmp4(state.range(0));
    std::vector<double> tmp5(state.range(0));

    std::vector<double> sigmaA(state.range(0));
    std::vector<double> sigmaA2T2(state.range(0));
    std::vector<double> sigmaAsqrtT(state.range(0));
    std::vector<double> emrt(state.range(0));
    std::vector<double> d1(state.range(0));
    std::vector<double> d2(state.range(0));


    // First create an instance of an engine.
    std::random_device rnd_device;

    // Specify the engine and distribution.
    std::mt19937 mersenne_engine{rnd_device()};  // Generates random integers

    std::uniform_real_distribution<double> ss_dist{1, 200};
    std::uniform_real_distribution<double> rs_dist{0.01, 1.0};
    std::uniform_real_distribution<double> sigmas_dist{0.01, 1.0};
    std::uniform_real_distribution<double> ts_dist{0.09, 2.0};
    std::uniform_real_distribution<double> taus_dist{0.08, 0.09};
    std::uniform_int_distribution<MKL_INT64> flags_dist{0, 3};

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
        std::generate(ss.begin(), ss.end(), ss_gen);
        std::generate(xs.begin(), xs.end(), ss_gen);
        std::generate(rs.begin(), rs.end(), rs_gen);
        std::generate(sigmas.begin(), sigmas.end(), sigmas_gen);
        std::generate(ts.begin(), ts.end(), ts_gen);
        std::generate(taus.begin(), taus.end(), taus_gen);
        std::generate(flags.begin(), flags.end(), flags_gen);

        for (MKL_INT64 i = 0; i < state.range(0); ++i) {
            ts[i] += taus[i];
        }

        init_mkl_pricer();

        prepare_mkl_pricer(
                state.range(0),
                ss.data(),
                sigmas.data(),
                ts.data(),
                taus.data(),
                rs.data(),
                tmp1.data(),
                tmp2.data(),
                tmp3.data(),
                tmp4.data(),
                tmp5.data(),
                sigmaA.data(),
                sigmaA2T2.data(),
                sigmaAsqrtT.data(),
                emrt.data());

        state.ResumeTiming();

        mkl_pricer(
                state.range(0),
                flags.data(),
                ss.data(),
                xs.data(),
                sigmaA2T2.data(),
                sigmaAsqrtT.data(),
                emrt.data(),
                tmp1.data(),
                tmp2.data(),
                tmp3.data(),
                tmp4.data(),
                d1.data(),
                d2.data(),
                prices.data());

    }

}
// Register the function as a benchmark
BENCHMARK(BM_Pricer_MKL)->Arg(1 << 20)->Arg(1 << 21)->Arg(1 << 22);
