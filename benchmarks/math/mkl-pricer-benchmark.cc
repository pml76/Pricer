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


#include <src/math/mkl_pricer.h>
#include <mkl.h>


static void BM_Erf_MKL(benchmark::State &state) {


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



    // First create an instance of an engine.
    std::random_device rnd_device;

    // Specify the engine and distribution.
    std::mt19937 mersenne_engine{rnd_device()};  // Generates random integers

    std::uniform_real_distribution<double> dist{-52, 52};


    auto gen = [&dist, &mersenne_engine]() {

        return dist(mersenne_engine);

    };


    for (auto _ : state) {

        state.PauseTiming();
        dbls.data()[0] = results.data()[0];
        std::generate(dbls.begin(), dbls.end(), gen);

        state.ResumeTiming();

        vdErf(state.range(0), dbls.data(), results.data());

    }

}
// Register the function as a benchmark
BENCHMARK(BM_Erf_MKL)->Arg(1 << 10)->Arg(1 << 11)->Arg(1 << 12)->Arg(1 << 20)->Arg(1 << 21)->Arg(1 << 22);
