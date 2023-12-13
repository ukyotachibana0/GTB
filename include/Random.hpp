/**
 * edit time : 2023-12-12
 * Random class
 */
#pragma once

#include <random>
#include <chrono>
#include <cstdlib>
#include "types.hpp"

namespace gl {
namespace gtb {

class Random
{
private:
    std::default_random_engine generator;

public:
    Random() {
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        generator = std::default_random_engine(seed);
    }

    ~Random() {}

    double nextReal() { // [0, 1]
        std::uniform_real_distribution<double> uni_db;
        return uni_db(generator);
    }

    double nextReal(double bound) { // [0, bound]
        std::uniform_real_distribution<double> uni_db(0, bound);
        return uni_db(generator);
    }

    int nextInt32(int bound) { // [0, bound]
        return (rand() % (bound + 1));
    }

    int_t nextInt(int_t bound) { // [0, bound]
        std::uniform_int_distribution<int_t> uni_int(0, bound);
        return uni_int(generator);
    }

    int_t nextIntExp(double theta) { // [0, \inf): exp distribution
        std::exponential_distribution<> exp_int(theta);
        return std::max(exp_int(generator), .0);
    }

    std::vector<int_t> nextInts(const std::vector<int_t>& w, int_t c, bool chosen) {
        std::discrete_distribution<int_t> dis_int(w.begin(), w.end());
        if (chosen) {
            std::vector<int_t> ans(c);
            for (int i = 0; i < c; i++) { ans[i] = dis_int(generator); }
            return ans;
        } else {
            std::vector<int_t> ans(w.size());
            for (int i = 0; i < c; i++) { ans[dis_int(generator)]++; }
            return ans;
        }
    }

    bool nextBool() {
        return nextInt(1);
    }
}; //! class Random

} //! namespace gtb
} //! namespace gl