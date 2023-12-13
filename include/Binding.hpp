/**
 * edit time : 2023-12-12
 * Binding hpp
 */

#pragma once
#include "headers.hpp"

namespace gl {
namespace gtb {

class Binding {
private:
    int_t min_ts;
    int_t max_ts;

    int_t N_bck;
    double MIN_lambda;
    double MAX_lambda;
    double B;
    double A;

public:
    Binding(int_t mit, int_t mat) : min_ts(mit), max_ts(mat), N_bck(8), MIN_lambda(0.001), MAX_lambda(10.0) {
        assert(mat > mit);
        A = (mat - mit) / 100.0;
        B = - 0.9 * A;
    }

    Binding(int_t mit, int_t mat, int_t n_bck, double mil, double mal) : min_ts(mit), max_ts(mat), N_bck(n_bck), MIN_lambda(mil), MAX_lambda(mal) {
        assert(mat > mit);
        A = (mat - mit) / 100.0;
        B = - 0.9 * A;
    }

    std::vector<std::vector<int_t>> binding(std::vector<int_t>& comm_size);

    static std::vector<std::vector<int_t>> unionWindow(const std::vector<std::vector<int_t>>& window);

    static std::vector<std::vector<int_t>> compleWindow(const std::vector<std::vector<int_t>>& window, int_t e_mit, int_t e_mat);
}; //! class Binding

class Solver {
public:
    double E;           // Expectation
    double m;           // maximum
    double initial_guess;
    double epsilon = 1e-6;

    Solver(double expectation, double ma, double guess) : E(expectation), m(ma), initial_guess(guess) {}
    
    double func(double lambda);

    double binary(double l, double r);

    double getK(double lambda);
}; //! class Solver

} //! namespace gtb
} //! namespace gl