/**
 * edit time : 2023-12-12
 * TEDTimestamp cpp
 */

#include "TEDTimestamp.hpp"

namespace gl {
namespace gtb {

int_t TEDTimestamp::mapRange(int_t mit, int_t mat, int_t ts) {
    return mit + (double)(ts - min_timestamp) / (double)(max_timestamp - min_timestamp) * (mat - mit);
}

void TEDTimestamp::preProcess() {
    int size = ts_range;
    // compute PDF
    std::vector<double> PDF(size);
    double sum = .0;
    for (int i = 0; i < size; i++) { 
        PDF[i] = pdf(i + 1);    // [1, ts_range]
        sum += PDF[i];
    }
    for (int i = 0; i < size; i++) {
        PDF[i] /= sum;
        step = (0 < PDF[i] && PDF[i] < step) ? PDF[i] : step;
    }
    step = (step < MIN_STEP) ? MIN_STEP : step;

    // compute CDF
    std::vector<double> CDF(size);
    std::partial_sum(PDF.begin(), PDF.end(), CDF.begin());
    
    // compute CDF_inv
    int steps = Utility::mathCeil(1.0 / step);
    CDF_inv.resize(steps);
    int_t j = 0;
    for (int i = 0; i < steps; i++) {
        double cur_step = i * step;
        while (j < size && CDF[j] < cur_step) { j++; }
        CDF_inv[i] = j + min_timestamp;
    }
}

int_t TEDTimestamp::genTimestamp(int_t mit, int_t mat) {
    double y = rand.nextReal();
    int_t ts = CDF_inv[y / step];
    return mapRange(mit, mat, ts);
    return ts;
}

int_t TEDTimestamp::genTimestamp(const std::vector<std::vector<int_t>>& window) {
    // get mapped ts
    int_t bound = 0, n = window.size();
    for (auto& win : window) { bound += win[1] - win[0]; }
    int_t ts_map = rand.nextInt(bound);
    // get original ts
    int_t last_sum = 0, sum = window[0][1] - window[0][0];
    for (int_t i = 1; i < n; i++) {
        if (ts_map <= sum) return (ts_map - last_sum) + window[i - 1][0];
        last_sum = sum;
        sum += window[i][1] - window[i][0];
    }
    return (ts_map - last_sum) + window[n - 1][0];
}

}
}