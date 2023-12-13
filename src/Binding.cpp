/**
 * edit time : 2023-12-02
 * Binding cpp
 */

#include "Binding.hpp"
#include "Utility.hpp"
#include "Sampler.hpp"

namespace gl {
namespace gtb {

std::vector<std::vector<int_t>> Binding::binding(std::vector<int_t>& comm_size) {
    int_t n = comm_size.size();
    assert(n >= 2);
    std::vector<int_t> bck_marks = {0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50};
    int_t n_bck = 11;
    double bck_step = 5;

    /* Samplers for each bucket */
    n_bck = bck_marks.size() - 1;
    std::vector<Sampler*> samplers = std::vector<Sampler*>(n_bck, nullptr);
    std::vector<bool> samplers_pre = std::vector<bool>(n_bck, false);
    for (int_t i = 0; i < n_bck; i++) {
        double E_T = B + A * log2((bck_marks[i] + bck_marks[i+1]) / 2);
        // sampler based on power-law distirbution
        std::unordered_map<std::string, double> params;
        Solver solver = Solver(E_T, max_ts - min_ts, 2.1);
        params["lambda"] = solver.binary(MIN_lambda, MAX_lambda);
        params["k"] = solver.getK(params["lambda"]);
        samplers[i] = new SamplerPowerLaw(1, (max_ts - min_ts), params);
        
        #ifdef DEBUG
        std::cout << "[Binding::binding] ---bucket " << i << "--- E_T: " << E_T << ", lambda: " << params["lambda"] << std::endl;
        #endif
    }

    /* Generate time windows for each time-bound community */
    Random rand;
    std::vector<int_t> length(n, 0);
    std::vector<std::vector<int_t>> ans(n, std::vector<int_t>({min_ts, max_ts}));
    for (int_t i = 0; i < n; i++) {
        // time window: length
        int_t i_bck = Utility::mathFloor(comm_size[i] / bck_step);
        i_bck = i_bck < n_bck ? i_bck : n_bck - 1;
        if (!samplers_pre[i_bck]) { samplers[i_bck]->preProcess(); samplers_pre[i_bck] = true; }
        int_t len = samplers[i_bck]->sample();
        // time window: left
        int_t left = rand.nextInt(max_ts - min_ts - len);

        length[i] = len;
        ans[i][0] = left;
        ans[i][1] = left + len;
    }

    #ifdef DEBUG
    std::cout << "[Binding::binding] time windows: "; Utility::show(ans);
    std::cout << "[Binding::binding] window lengths: "; Utility::show(length);
    #endif
    
    return ans;
}

std::vector<std::vector<int_t>> Binding::unionWindow(const std::vector<std::vector<int_t>>& window) {
    std::vector<std::vector<int_t>> ans;
    int_t n = window.size(), i = 0;
    while (i < n - 1) {
        int_t left = window[i][0], right = window[i][1], j = i;
        while (j < n - 1 && window[j + 1][0] <= window[j][1]) { 
            right = window[++j][1];
        }
        ans.push_back({left, right});
        i = j + 1;
    }
    if (i == n - 1) ans.push_back(window.back());

    return ans;
}

std::vector<std::vector<int_t>> Binding::compleWindow(const std::vector<std::vector<int_t>>& window, int_t e_mit, int_t e_mat) {
    std::vector<std::vector<int_t>> ans;
    if (e_mit < window.front()[0]) ans.push_back({e_mit, window.front()[0]});
    int n = window.size();
    for (int i = 0; i < n - 1; ++i) { ans.push_back({window[i][1], window[i + 1][0]}); }
    if (e_mat > window.back()[1]) ans.push_back({window.back()[1], e_mat});

    return ans;
}

double Solver::func(double lambda) {
    return (lambda - 1) / (lambda - 2) * (1 - pow(m, 2-lambda)) / (1 - pow(m, 1-lambda)) - E;
}

double Solver::binary(double l, double r) {
    double mid = .0;
    int_t cnt = 0;
    while (l < r) {
        mid = (l + r) / 2.0;
        cnt ++;
        double ans = func(mid);
        if (fabs(ans) < epsilon) {
            return mid;
        }
        if (ans > 0) l = mid;
        else r = mid;
    }
    return mid;
}

double Solver::getK(double lambda) {
    return (lambda - 1) / (pow(m, 1-lambda) - 1);
}

} //! namespace gtb
} //! namespace gl