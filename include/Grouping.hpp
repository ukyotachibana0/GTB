/**
 * edit time : 2023-12-12
 * Grouping hpp
 */

#pragma once
#include "headers.hpp"

namespace gl {
namespace gtb {

class Grouping {
private:
    int_t MIN_size;

public:
    Grouping() : MIN_size(2) {}

    Grouping(int_t mis) : MIN_size(mis) {}

    std::vector<std::vector<int_t>> splitCommunity(int_t row, int_t col, int_t k, double lambda);

    std::vector<int_t> splitScalar(int_t n, int_t k, double lambda);

    static std::vector<std::unordered_map<int_t, double>> idenOlComm(int_t n, int_t m, double omega_min, double omega_max);

    static std::vector<std::vector<int_t>> homoOlRange(const std::vector<std::vector<int_t>>& comm_split_psum, std::unordered_map<int_t, double>& ol_comm, int_t i);

    static std::vector<std::vector<std::vector<int_t>>> heteOlRange(const std::vector<std::vector<int_t>>& comm_split_psum, std::unordered_map<int_t, double>& ol_comm, int_t i);
};

} //! namespace gtb
} //! namespace gl