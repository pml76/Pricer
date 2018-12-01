//
// Created by peter on 11/22/18.
//


#include <benchmark/benchmark.h>
#include <src/math/erf_wrapper.h>
#include <vector>
#include <random>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <functional>
#include <cmath>

static void BM_Erf_IEEE754(benchmark::State &state) {


    std::vector<double> dbls(state.range(0));

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

        std::generate(dbls.begin(), dbls.end(), gen);

        state.ResumeTiming();

        for (int j = 0; j < state.range(0); ++j) {
            benchmark::DoNotOptimize(ieee754_erf(dbls[j]));
        }

    }

}
// Register the function as a benchmark
BENCHMARK(BM_Erf_IEEE754)->Arg(1 << 10)->Arg(1 << 11)->Arg(1 << 12);


static void BM_Erf_Default(benchmark::State &state) {


    std::vector<double> dbls(state.range(0));

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

        std::generate(dbls.begin(), dbls.end(), gen);

        state.ResumeTiming();

        for (int j = 0; j < state.range(0); ++j) {
            benchmark::DoNotOptimize(erf(dbls[j]));
        }

    }

}
// Register the function as a benchmark
BENCHMARK(BM_Erf_Default)->Arg(1 << 10)->Arg(1 << 11)->Arg(1 << 12);