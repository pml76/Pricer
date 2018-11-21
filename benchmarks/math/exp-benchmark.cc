//
// Created by peter on 11/21/18.
//

#include <benchmark/benchmark.h>
#include <src/math/exp_wrapper.h>
#include <vector>
#include <random>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <functional>

static void BM_IEEE754_Exp(benchmark::State &state) {


    std::vector<double> dbls(1024);

    // First create an instance of an engine.
    std::random_device rnd_device;

    // Specify the engine and distribution.
    std::mt19937 mersenne_engine{rnd_device()};  // Generates random integers

    std::uniform_int_distribution<int> dist{-52, 52};


    auto gen = [&dist, &mersenne_engine]() {

        return dist(mersenne_engine);

    };


    for (auto _ : state) {

        state.PauseTiming();

        std::generate(dbls.begin(), dbls.end(), gen);

        state.ResumeTiming();

        for (int j = 0; j < state.range(1); ++j) {
            ieee754_exp(dbls[j]);
        }

    }

}
// Register the function as a benchmark
BENCHMARK(BM_IEEE754_Exp)->Arg(1024);