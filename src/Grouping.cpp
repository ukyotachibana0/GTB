/**
 * edit time : 2023-12-12
 * Grouping cpp
 */

#include "Grouping.hpp"

#include "Utility.hpp"
#include "Sampler.hpp"
#include "Random.hpp"
#include <assert.h>

namespace gl {
namespace gtb {

std::vector<std::vector<int_t>> Grouping::splitCommunity(int_t row, int_t col, int_t k, double lambda) {
    std::vector<int_t> row_ans = splitScalar(row, k, lambda);

    std::vector<std::vector<int_t>> ans(row_ans.size(), std::vector<int_t>(2));
    if (row == col) {
        for (size_t i = 0; i < row_ans.size(); ++i)
            ans[i][0] = ans[i][1] = row_ans[i];
    } else {
        int_t cum = 0;
        for (int_t i = row_ans.size() - 1; i > 0; --i) {
            ans[i][0] = row_ans[i];
            ans[i][1] = (int_t)(Utility::mathRound(row_ans[i] * 1.0 * col / row * 1.0));
            cum += ans[i][1];
        }
        ans[0][0] = row_ans[0];
        ans[0][1] = col - cum;
    }
    return ans;
}

std::vector<int_t> Grouping::splitScalar(int_t n, int_t k, double lambda) {
    lambda = -fabs(lambda);
    int_t avg_size = n / k;
    int_t step = 0;
    int_t memory = 0;
    int_t lr = 0;
    int_t mem_count = 0;
    std::vector<int_t> al_size(1, 0);
    std::vector<int_t> al_count(1, 0);
    while (true) {
        int_t curr_size = 1 + lr * 2;
        std::vector<int_t> size_list(curr_size);
        size_list[lr] = avg_size;
        step = std::max(avg_size / (lr + k + 2.0), 1.0);
        for (int_t i = 1; i <= lr; ++i) {
            size_list[lr + i] = size_list[lr + i - 1] + step;
            size_list[lr - i] = size_list[lr - i + 1] - step;
        }
        std::vector<double> probability(curr_size);
        for (int_t i = 0; i < curr_size; ++i)
            probability[i] = pow((double)(size_list[i]), lambda);
        double last_p = probability[curr_size - 1];
        std::vector<int_t> count = std::vector<int_t>(curr_size, 0);
        int_t spectrum = 0;
        for (int_t i = 0; i < curr_size; ++i) {
            if (std::isnan(probability[i]) || std::isinf(probability[i])) continue;
            count[i] = (int_t)Utility::mathRound(probability[i] / last_p);
            spectrum += count[i] * size_list[i];
        }
        if (spectrum > n) break;
        memory = spectrum;
        mem_count = 0;
        for (int_t i = 0; i < curr_size; ++i) {
            mem_count += count[i];
            if (i >= al_size.size()) {
                al_size.push_back(size_list[i]);
                al_count.push_back(count[i]);
            } else {
                al_size[i] = size_list[i];
                al_count[i] = count[i];
            }
        }
        lr += 1;
    }
    al_size.back() += n - memory;

    mem_count = 0;
    for (int_t i = 0; i < al_count.size(); ++i)
        mem_count += al_count[i];

    std::vector<int_t> res(mem_count + 2);
    int_t j = 0;
    for (int_t i = al_size.size() - 1; i >= 0; --i) {
        for (int_t p = 0; p < al_count[i]; ++p) {
            if (j < res.size())
                res[j ++] = al_size[i];
        }
    }

    std::vector<int_t> ans(k);
    if (mem_count < k) {
        int_t fission = k / mem_count + ((k % mem_count == 0) ? 0 : 1);
        int_t more = k - mem_count;
        int_t res_i = 0, ans_i = 0;
        while (more > 0) {
            int_t cumu_size = 0;
            int_t per_size = res[res_i] / fission;
            for (int_t i = 0; i < fission - 1; ++i) {
                ans[ans_i] = per_size;
                ans_i ++;
                cumu_size += per_size;
            }
            ans[ans_i] = res[res_i] - cumu_size;
            ans_i ++;
            res_i ++;
            more -= (fission - 1);
        }
        while (ans_i < k && res_i < mem_count) {
            ans[ans_i] = res[res_i];
            ans_i ++;
            res_i ++;
        }
    } else if (mem_count > k) {
        int_t fission = mem_count / k + ((mem_count % k == 0) ? 0 : 1);
        int_t less = mem_count - k;
        int_t res_i = mem_count - 1;
        int_t ans_i = k - 1;
        while (less > 0) {
            ans[ans_i] = res[res_i];
            res_i --;
            for (int_t i = 0; i < fission - 1; ++i) {
                ans[ans_i] += res[res_i];
                res_i --;
            }
            less -= (fission - 1);
            ans_i --;
        }
        while (ans_i >= 0 && res_i >= 0) {
            ans[ans_i] = res[res_i];
            ans_i --;
            res_i --;
        }
    } else {
        std::copy_n(res.begin(), k, ans.begin());
    }

    int_t sum = 0;
    for (int_t i = 0; i < k; i++) { sum += ans[i]; } 
    if (sum != n) ans[0] += n - sum;

    bool less = false;
    for (int_t i = 0; i < k; ++i) { if (ans[i] < MIN_size) {less = true; break;} }
    if (!less) return ans;

    // adjust to [MIN_size, )
    Random rand;
    int_t adjust = 0;
    double alpha = (n - MIN_size * k) / (double)n;
    for (int_t i = 0; i < k; ++i) {
        bool floor = rand.nextBool();
        int_t ans_new = (floor ? ans[i] * alpha : Utility::mathCeil(ans[i] * alpha));
        adjust += ans[i] - ans_new;
        ans[i] = ans_new;
    }
    for (int_t i = 0; i < k; ++i) ans[i] += MIN_size;
    ans[0] += adjust - MIN_size * k;
    if (ans[0] < MIN_size) {    // borrow from others
        ans[1] -= MIN_size - ans[0];
        ans[0] = MIN_size;
    }

    #ifdef DEBUG
    std::cout << "[Grouping::splitScalar] ans: ";
    Utility::show(ans);
    #endif
    
    return ans;
}

std::vector<std::unordered_map<int_t, double>> Grouping::idenOlComm(int_t n, int_t m, double omega_min, double omega_max) {
    assert(m <= n * (n-1) / 2);
    std::vector<std::unordered_map<int_t, double>> ans(n);

    int_t i = 0;
    Random rand;
    while (i < m) {
        std::pair<int_t, int_t> one_pair(rand.nextInt(n - 1), rand.nextInt(n - 1));
        if (one_pair.first == one_pair.second) continue;
        double ol = rand.nextReal(omega_max - omega_min) + omega_min;
        if (ans[one_pair.first].insert({one_pair.second, ol}).second &&
            ans[one_pair.second].insert({one_pair.first, ol}).second) 
            ++i;
    }

    #ifdef DEBUG
    std::cout << "[Grouping::idenOlComm] overlap relationships:" << std::endl;
    for (int i = 0; i < n; i++) {
        std::cout << "\t" << i << ": ";
        for (auto p : ans[i]) std::cout << "(" << p.first << "," << p.second << ") ";
        std::cout << std::endl;
    }
    #endif
    
    return ans;
}

std::vector<std::vector<int_t>> Grouping::homoOlRange(const std::vector<std::vector<int_t>>& comm_split_psum, std::unordered_map<int_t, double>& ol_comm, int_t i) {
    std::vector<std::vector<int_t>> ans;
    for (auto& p : ol_comm) {
        int_t j = p.first;
        if (i > j) {
            int_t j_size = comm_split_psum[j + 1][0] - comm_split_psum[j][0];
            int_t sp_size = j_size - j_size * p.second;
            ans.push_back({ comm_split_psum[j][0] + sp_size, comm_split_psum[j + 1][0] - 1 });
        }
    }
    return ans;
}

std::vector<std::vector<std::vector<int_t>>> Grouping::heteOlRange(const std::vector<std::vector<int_t>>& comm_split_psum, std::unordered_map<int_t, double>& ol_comm, int_t i) {
    std::vector<std::vector<std::vector<int_t>>> ans;
    for (auto& p : ol_comm) {
        int_t j = p.first;
        if (i > j) {
            int_t j_size_s = comm_split_psum[j + 1][0] - comm_split_psum[j][0];
            int_t j_size_t = comm_split_psum[j + 1][1] - comm_split_psum[j][1];
            int_t sp_size_s = j_size_s - j_size_s * p.second;
            int_t sp_size_t = j_size_t - j_size_t * p.second;
            ans.push_back({
                { comm_split_psum[j][0] + sp_size_s, comm_split_psum[j][1] + sp_size_t },
                { comm_split_psum[j + 1][0] - 1, comm_split_psum[j + 1][1] - 1 }
            });
        }
    }
    return ans;
}

} //! namespace gtb
} //! namespace gl