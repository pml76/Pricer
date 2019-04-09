/*
 *
 * (c) 2019, by Peter Lennartz  // peter.lennartz@gmail.com
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */


#include <benchmark/benchmark.h>
#include <random>
#include <iterator>
#include <iostream>
#include <functional>


#include <math/pricers/sleef_pricer.h>

#include <pricer-dispatcher.h>


typedef double *__restrict__ __attribute__((aligned(ALIGN_TO))) Real_Ptr;

static void BM_Pricer_Sleef(benchmark::State &state) {

    Pricer::pricer_context context( state.range(0));


    // First create an instance of an engine.
    std::random_device rnd_device;

    // Specify the engine and distribution.
    std::mt19937 mersenne_engine{rnd_device()};  // Generates random integers

    std::uniform_real_distribution<double> ss_dist{1, 200};
    std::uniform_real_distribution<double> rs_dist{0.01, 1.0};
    std::uniform_real_distribution<double> sigmas_dist{0.01, 1.0};
    std::uniform_real_distribution<double> ts_dist{0.09, 2.0};
    std::uniform_real_distribution<double> taus_dist{0.08, 0.09};
    std::discrete_distribution<int> flags_dist{{50, 0, 50, 0}};

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
        return flags_dist(mersenne_engine) - 1.;
    };

    for (auto _ : state) {

        state.PauseTiming();

        context.get_prices()[0] = context.get_prices()[0] / 2;
        std::generate(context.get_s(), &context.get_s()[state.range(0)], ss_gen);
        std::generate(context.get_x(), &context.get_x()[state.range(0)], ss_gen);
        std::generate(context.get_r(), &context.get_r()[state.range(0)], rs_gen);
        std::generate(context.get_sigma(), &context.get_sigma()[state.range(0)], sigmas_gen);
        std::generate(context.get_t(), &context.get_t()[state.range(0)], ts_gen);
        std::generate(context.get_tau(), &context.get_tau()[state.range(0)], taus_gen);
        std::generate(context.get_long_short(), &context.get_long_short()[state.range(0)], flags_gen);
        std::generate(context.get_put_call(), &context.get_put_call()[state.range(0)], flags_gen);

        for (int64_t i = 0; i < state.range(0); ++i) {
            context.get_t()[i] += context.get_tau()[i];
        }

        init_tw_pricer();

        prepare_tw_pricer(context);

        state.ResumeTiming();

        tw_pricer(context);

    }


}
// Register the function as a benchmark
BENCHMARK(BM_Pricer_Sleef)->Arg(1 << 20)->Arg(1 << 21)->Arg(1 << 22);


static void BM_Pricer_Full_Sleef(benchmark::State &state) {


    Pricer::d2dx2_pricer_context context( state.range(0));


    // First create an instance of an engine.
    std::random_device rnd_device;

    // Specify the engine and distribution.
    std::mt19937 mersenne_engine{rnd_device()};  // Generates random integers

    std::uniform_real_distribution<double> ss_dist{1, 200};
    std::uniform_real_distribution<double> rs_dist{0.01, 1.0};
    std::uniform_real_distribution<double> sigmas_dist{0.01, 1.0};
    std::uniform_real_distribution<double> ts_dist{0.09, 2.0};
    std::uniform_real_distribution<double> taus_dist{0.08, 0.09};
    std::discrete_distribution<int> flags_dist{{50, 0, 50, 0}};

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
        return flags_dist(mersenne_engine) - 1.;
    };

    for (auto _ : state) {

        state.PauseTiming();

        context.get_prices()[0] = context.get_prices()[0] / 2;
        std::generate(context.get_s(), &context.get_s()[state.range(0)], ss_gen);
        std::generate(context.get_x(), &context.get_x()[state.range(0)], ss_gen);
        std::generate(context.get_r(), &context.get_r()[state.range(0)], rs_gen);
        std::generate(context.get_sigma(), &context.get_sigma()[state.range(0)], sigmas_gen);
        std::generate(context.get_t(), &context.get_t()[state.range(0)], ts_gen);
        std::generate(context.get_tau(), &context.get_tau()[state.range(0)], taus_gen);
        std::generate(context.get_long_short(), &context.get_long_short()[state.range(0)], flags_gen);
        std::generate(context.get_put_call(), &context.get_put_call()[state.range(0)], flags_gen);

        for (int64_t i = 0; i < state.range(0); ++i) {
            context.get_t()[i] += context.get_tau()[i];
        }

        init_tw_pricer();

        prepare_tw_pricer(context);

        state.ResumeTiming();

        full_tw_pricer(context);

    }

}
// Register the function as a benchmark
BENCHMARK(BM_Pricer_Full_Sleef)->Arg(1 << 20)->Arg(1 << 21)->Arg(1 << 22);