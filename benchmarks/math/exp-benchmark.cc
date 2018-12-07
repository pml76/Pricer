//
// Created by peter on 11/21/18.
//

#include <benchmark/benchmark.h>
#include <vector>
#include <random>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <functional>
#include <cmath>

#include <mkl.h>


static void BM_Exp_MKL(benchmark::State &state) {


    vmlSetMode(VML_EP);
    vmlSetMode(VML_FTZDAZ_OFF);
    vmlSetMode(VML_ERRMODE_NOERR);

    std::vector<double> dbls(state.range(0));
    std::vector<double> results(state.range(0));

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

        vdExp(state.range(0), dbls.data(), results.data());

    }

}
// Register the function as a benchmark
BENCHMARK(BM_Exp_MKL)->Arg(1 << 10)->Arg(1 << 11)->Arg(1 << 12);

static void BM_Exp_Default(benchmark::State &state) {


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
            benchmark::DoNotOptimize(exp(dbls[j]));
        }

    }

}
// Register the function as a benchmark
BENCHMARK(BM_Exp_Default)->Arg(1 << 10)->Arg(1 << 11)->Arg(1 << 12);