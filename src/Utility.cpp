/**
 * edit time : 2023-12-12
 * Utility
 */

#include "Utility.hpp"

namespace gl {
namespace gtb {

double Utility::mathCeil(double x) {
    double y = (double)((int)x);
    if (y == x) return x;
    if (y >= 0.0) return (y + 1.0);
    else return (y - 1.0);
}

double Utility::mathRound(double x) {
    double y = (double)((int)x);
    double delta = x - y;
    if (delta > 0.5) return (y + 1.0);
    else if (delta >= -0.5 && delta <= 0.5) return y;
    else return y - 1.0;
}

double Utility::mathFloor(double x) {
    return floor(x);
}

double Utility::normPdf(double x) {
    return exp(-x * x / 2.0) / sqrt(2.0 * PI);
}

double Utility::normCdf(double x) {
    if (x < -8.0) return 0.0;
    if (x > 8.0) return 1.0;
    double sum = 0.0, term = x;
    for (int i = 3; sum + term != sum; i += 2) {
        sum += term;
        term *= x * x / i;
    }
    return 0.5 + sum * normPdf(x);
}

void Utility::show(const std::vector<int_t>& list) {
    std::cout << "[";
    for (auto n : list)
        std::cout << n << " ";
    std::cout << "]" << std::endl;
}

void Utility::show(const std::vector<std::vector<int_t>>& list) {
    for (auto l : list) {
        std::cout << "[";
        for (auto n : l)
            std::cout << n << " ";
        std::cout << "], ";
    }
    std::cout << std::endl;
}

int_t Utility::min(const std::vector<int_t>& list) {
    int_t ans = list[0];
    for (auto l : list) { if (l < ans) { ans = l; } }
    return ans;
}

int_t Utility::max(const std::vector<int_t>& list) {
    int_t ans = list[0];
    for (auto l : list) { if (l > ans) { ans = l; } }
    return ans;
}

int_t Utility::max(const std::vector<std::vector<int_t>> list, int col) {
    auto list_col = getColumn(list, col);
    return max(list_col);
}

int_t Utility::max2(const std::vector<int_t>& list) {
    int_t ans = -1, ma = max(list);
    for (auto l : list) { if (l > ans && l < ma) { ans = l; } }
    return ans;
}

std::vector<int_t> Utility::getColumn(const std::vector<std::vector<int_t>>& list, int col) {
    int n = list.size();
    std::vector<int_t> list_col(n);
    for (int i = 0; i < n; i++) { list_col[i] = list[i][col]; }
    return list_col;
}

int_t Utility::max_diff(const std::vector<std::vector<int_t>> list) {
    assert(list[0].size() == 2);
    int_t n = list.size();
    int_t ans = list[0][1] - list[0][0], ans_i = 0;
    for (int_t i = 0; i < n; i++) { 
        if (list[i][1] - list[i][0] > ans) { 
            ans = list[i][1] - list[i][0];  
            ans_i = i;
        }
    }
    return ans_i;
}

int Utility::numOneBitInt(uint32_t x) {
    int ans = 0;
    while (x) {
        ans ++;
        x &= (x - 1);
    }
    return ans;
}

std::string Utility::strSeqNumber(int x, int n_bits) {
    // char *buffer = new char[n_bits + 1];
    char buffer[32];
    std::string format = "%0" + std::to_string(n_bits) + "d";
    int res = snprintf(buffer, sizeof(buffer), format.c_str(), x);
    std::string ans = "";
    if (res >= 0 && res < n_bits + 1)
        ans = buffer;
    return ans;
}

void Utility::randomChoice(std::vector<int>& nums, int K) {
    int N = nums.size();
    if (K >= N) {
        return;
    }
    for (int i = 0; i < K; ++i) {
        int x = (rand() % N) + i;
        std::swap(nums[i], nums[x]);
        N --;
    }
}

void Utility::randomChoice(std::vector<int>& nums, std::vector<int>& accomp, int K) {
    int N = nums.size();
    if (K >= N) {
        return;
    }
    for (int i = 0; i < K; ++i) {
        int x = (rand() % N) + i;
        std::swap(nums[i], nums[x]);
        std::swap(accomp[i], accomp[x]);
        N --;
    }
}


bool ConfChecker::checkJson(nlohmann::json& json_obj) {
    // "graph": string
    bool ans = true;
    if (json_obj.find(schema::json_graph) == json_obj.end()) {
        std::cerr << "[ConfChecker::checkJson] Lack Error: JSON['" << schema::json_graph << "']." << std::endl;
        ans = false;
    }
    if (ans && !json_obj[schema::json_graph].is_string()) {
        std::cerr << "[ConfChecker::checkJson] Type Error: JSON['" << schema::json_graph << "']." << std::endl;
        ans = false;
    }
    // node schema
    // "node": [{"label": string, "amount": number}, ...]
    if (json_obj.find(schema::json_node) == json_obj.end()) {
        std::cerr << "[ConfChecker::checkJson] Lack Error: JSON['" << schema::json_node << "']." << std::endl;
        return false;
    }
    if (!json_obj[schema::json_node].is_array()) {
        std::cerr << "[ConfChecker::checkJson] Type Error: JSON['" << schema::json_node << "'] (List of objects)." << std::endl;
        return false;
    }
    auto& node_schema = json_obj[schema::json_node];
    std::unordered_set<std::string> node_set;
    for (auto& one_node : node_schema) {
        if (one_node.find(schema::json_node_label) == one_node.end()) {
            std::cerr << "[ConfChecker::checkJson] Lack Error: JSON['" << schema::json_node << "'][i]['" << schema::json_node_label << "']." << std::endl;
            ans = false;
            continue;
        }
        if (!one_node[schema::json_node_label].is_string()) {
            std::cerr << "[ConfChecker::checkJson] Type Error: JSON['" << schema::json_node << "'][i]['" << schema::json_node_label << "'] (string)." << std::endl;
            ans = false;
            continue;
        }
        if (node_set.count(one_node[schema::json_node_label])) {
            std::cerr << "[ConfChecker::checkJson] Content Error: JSON['" << schema::json_node << "'][i]['" << schema::json_node_label << "'] duplicated: " << one_node[schema::json_node_label] << std::endl;
            ans = false;
            continue;
        }
        std::string __node_label = one_node[schema::json_node_label];
        node_set.insert(__node_label);
        if (one_node.find(schema::json_node_amount) == one_node.end()) {
            std::cerr << "[ConfChecker::checkJson] Lack Error: JSON['" << schema::json_node << "'][i]['" << schema::json_node_amount << "']." << std::endl;
            ans = false;
            continue;
        }
        if (!one_node[schema::json_node_amount].is_number()) {
            std::cerr << "[ConfChecker::checkJson] Type Error: JSON['" << schema::json_node << "'][i]['" << schema::json_node_amount << "'] (number)." << std::endl;
            ans = false;
            continue;
        }
    }
    // edge schema
    // "edge": [{"label": string, "source": string, "target": string, "amount": number, "out": {}, "in"(optional): {}, "community"(optional): {}, "temporal"(optional): {}, ...]
    if (json_obj.find(schema::json_edge) == json_obj.end()) {
        std::cerr << "[ConfChecker::checkJson] Lack Error: JSON['" << schema::json_edge << "']." << std::endl;
        return false;
    }
    if (!json_obj[schema::json_edge].is_array()) {
        std::cerr << "[ConfChecker::checkJson] Type Error: JSON['" << schema::json_edge << "'] (List of objects)." << std::endl;
        return false;
    }
    auto& edge_schema = json_obj[schema::json_edge];
    std::unordered_set<std::string> edge_set;
    for (auto& one_edge : edge_schema) {
        // label
        if (one_edge.find(schema::json_edge_label) == one_edge.end()) {
            std::cerr << "[ConfChecker::checkJson] Lack Error: JSON['" << schema::json_edge << "'][i]['" << schema::json_edge_label << "']." << std::endl;
            ans = false;
            continue;
        }
        if (!one_edge[schema::json_edge_label].is_string()) {
            std::cerr << "[ConfChecker::checkJson] Type Error: JSON['" << schema::json_edge << "'][i]['" << schema::json_edge_label << "'] (string)." << std::endl;
            ans = false;
            continue;
        }
        if (edge_set.count(one_edge[schema::json_edge_label])) {
            std::cerr << "[ConfChecker::checkJson] Content Error: JSON['" << schema::json_edge << "'][i]['" << schema::json_edge_label << "'] duplicated: " << one_edge[schema::json_edge_label] << std::endl;
            ans = false;
            continue;
        }
        std::string __edge_label = one_edge[schema::json_edge_label];
        edge_set.insert(__edge_label);
        // source
        if (one_edge.find(schema::json_edge_source) == one_edge.end()) {
            std::cerr << "[ConfChecker::checkJson] Lack Error: JSON['" << schema::json_edge << "'][i]['" << schema::json_edge_source << "']." << std::endl;
            ans = false;
            continue;
        }
        if (!one_edge[schema::json_edge_source].is_string()) {
            std::cerr << "[ConfChecker::checkJson] Type Error: JSON['" << schema::json_edge << "'][i]['" << schema::json_edge_source << "'] (string)." << std::endl;
            ans = false;
            continue;
        }
        if (!node_set.count(one_edge[schema::json_edge_source])) {
            std::cerr << "[ConfChecker::checkJson] Content Error: JSON['" << schema::json_edge << "'][i]['" << schema::json_edge_source << "']." << std::endl;
            ans = false;
            continue;
        }
        // target
        if (one_edge.find(schema::json_edge_target) == one_edge.end()) {
            std::cerr << "[ConfChecker::checkJson] Lack Error: JSON['" << schema::json_edge << "'][i]['" << schema::json_edge_target << "']." << std::endl;
            ans = false;
            continue;
        }
        if (!one_edge[schema::json_edge_target].is_string()) {
            std::cerr << "[ConfChecker::checkJson] Type Error: JSON['" << schema::json_edge << "'][i]['" << schema::json_edge_target << "'] (string)." << std::endl;
            ans = false;
            continue;
        }
        if (!node_set.count(one_edge[schema::json_edge_target])) {
            std::cerr << "[ConfChecker::checkJson] Content Error: JSON['" << schema::json_edge << "'][i]['" << schema::json_edge_target << "']." << std::endl;
            ans = false;
            continue;
        }
        // amount
        if (one_edge.find(schema::json_edge_amount) == one_edge.end()) {
            std::cerr << "[ConfChecker::checkJson] Lack Error: JSON['" << schema::json_edge << "'][i]['" << schema::json_edge_amount << "']." << std::endl;
            ans = false;
            continue;
        }
        if (!one_edge[schema::json_edge_amount].is_number()) {
            std::cerr << "[ConfChecker::checkJson] Type Error: JSON['" << schema::json_edge << "'][i]['" << schema::json_edge_amount << "'] (number)." << std::endl;
            ans = false;
            continue;
        }
        // out distribution
        std::string out_info = "[ConfChecker::checkJson] Error: JSON['" + schema::json_edge + "'][i]['" + schema::json_ted_in + "']";
        if (one_edge.find(schema::json_ted_out) == one_edge.end()) {
            std::cerr << out_info << "(lack)." << std::endl;
            ans = false;
            continue;
        }
        auto& out_dist = one_edge[schema::json_ted_out];
        if (!checkJsonDegr(out_dist, out_info)) { ans = false; continue; }
        // in distribution
        std::string in_info = "[ConfChecker::checkJson] Error: JSON['" + schema::json_edge + "'][i]['" + schema::json_ted_in + "']";
        if (one_edge.find(schema::json_ted_in) == one_edge.end()) {
            std::cerr << in_info << "(lack)." << std::endl;
            ans = false;
            continue;
        }
        auto& in_dist = one_edge[schema::json_ted_in];
        if (!checkJsonDegr(in_dist, in_info)) { ans = false; continue; }
        // temporal (optional)
        if (one_edge.find(schema::json_ted_temp) != one_edge.end()) {
            auto& temp = one_edge[schema::json_ted_temp];
            std::string info = "[ConfChecker::checkJson] Error: JSON['" + schema::json_edge + "'][i]['" + schema::json_ted_temp + "']";
            if (!checkJsonTemp(temp, info)) { ans = false; continue; }
        }
        // community (optional)
        if (one_edge.find(schema::json_comm) != one_edge.end()) {
            auto& comm = one_edge[schema::json_comm];
            std::string info = "[ConfChecker::checkJson] Error: JSON['" + schema::json_edge + "'][i]['" + schema::json_comm + "']";
            if (!checkJsonComm(comm, info, true)) { ans = false; continue; }
        }
    }
    // "store-format": string
    if (json_obj.find(schema::json_store_format) != json_obj.end()) {
        if (!json_obj[schema::json_store_format].is_string()) {
            std::cerr << "[ConfChecker::checkJson] Type Error: JSON['" << schema::json_store_format << "'] (store format string)." << std::endl;
            ans = false;
        }
    }
    return ans;
}


bool ConfChecker::checkJsonTemp(nlohmann::json& temp, std::string info) {
    if (!temp.is_object()) {
        std::cerr << info << " (object)." << std::endl;
        return false;
    }
    // type
    if (temp.find(schema::json_dist_type) == temp.end()) {
        std::cerr << info << " (type)." << std::endl;
        return false;
    }
    if (!temp[schema::json_dist_type].is_string()) {
        std::cerr << info << " (type: string)." << std::endl;
        return false;
    }
    std::string temp_type = temp[schema::json_dist_type];
    // min_ts
    if (temp.find(schema::json_ted_min_timestamp) == temp.end()) {
        std::cerr << info << " (min_ts)." << std::endl;
        return false;
    }
    if (!temp[schema::json_ted_min_timestamp].is_number()) {
        std::cerr << info << " (min_ts: number)." << std::endl;
        return false;
    }
    // max_ts
    if (temp.find(schema::json_ted_max_timestamp) == temp.end()) {
        std::cerr << info << " (max_ts)." << std::endl;
        return false;
    }
    if (!temp[schema::json_ted_max_timestamp].is_number()) {
        std::cerr << info << " (max_ts: number)." << std::endl;
        return false;
    }
    if (temp_type == schema::json_dist_PowerLaw || temp_type == schema::json_dist_Exponential) {
        if (temp.find(schema::json_dist_lambda) == temp.end()) {
            std::cerr << info << " (PowerLaw lambda)." << std::endl;
            return false;
        }
        if (!temp[schema::json_dist_lambda].is_number()) {
            std::cerr << info << " (PowerLaw lambda: number)." << std::endl;
            return false;
        }
    } else if (temp_type == schema::json_dist_Normal || temp_type == schema::json_dist_LogNormal) {
        // mu (optional)
        if (temp.find(schema::json_dist_mu) != temp.end() && !temp[schema::json_dist_mu].is_number()) {
            std::cerr << info << " ([Log]Normal mu: number)." << std::endl;
            return false;
        }
        // sigma (required)
        if (temp.find(schema::json_dist_sigma) == temp.end()) {
            std::cerr << info << " ([Log]Normal sigma)." << std::endl;
            return false;
        }
        if (!temp[schema::json_dist_sigma].is_number()) {
            std::cerr << info << " ([Log]Normal sigma: number)." << std::endl;
            return false;
        }
    } else if (temp_type == schema::json_dist_Uniform) {
        //
    } else {
        std::cerr << info << " (unknown distribution or not implement: " << temp_type << ")." << std::endl;
        return false;
    }
    return true;
}

bool ConfChecker::checkJsonComm(nlohmann::json& comm, std::string info, bool required) {
    if (!comm.is_object()) {
        std::cerr << info << " (object)." << std::endl;
        return false;
    }
    // community amount
    if (comm.find(schema::json_comm_amount) == comm.end()) {
        if (required) { std::cerr << info << " (amount)." << std::endl; return false; }
    } else if (!comm[schema::json_comm_amount].is_number()) {
        std::cerr << info << " (amount: number)." << std::endl;
        return false;
    }
    // delta
    if (comm.find(schema::json_comm_delta) == comm.end()) {
        if (required) { std::cerr << info << " (delta)." << std::endl; return false; }
    } else if (!comm[schema::json_comm_delta].is_number()) {
        std::cerr << info << " (delta: number)." << std::endl;
        return false;
    }
    // lambda
    if (comm.find(schema::json_comm_lambda) == comm.end()) {
        if (required) { std::cerr << info << " (lambda)." << std::endl; return false; }
    } else if (!comm[schema::json_comm_lambda].is_number()) {
        std::cerr << info << " (lambda: number)." << std::endl;
        return false;
    }
    // overlap (optional)
    if (comm.find(schema::json_comm_ol) != comm.end()) {
        auto& comm_ol = comm[schema::json_comm_ol];
        // m
        if (comm_ol.find(schema::json_comm_ol_m) == comm_ol.end()) {
            std::cerr << info << "[overlap] (m)." << std::endl; 
            return false;
        } else if (!comm_ol[schema::json_comm_ol_m].is_number()) {
            std::cerr << info << "[overlap] (m: number)." << std::endl; 
            return false;
        }
        // min omega
        if (comm_ol.find(schema::json_comm_ol_min_omega) == comm_ol.end()) {
            std::cerr << info << "[overlap] (min_omega)." << std::endl; 
            return false;
        } else if (!comm_ol[schema::json_comm_ol_min_omega].is_number()) {
            std::cerr << info << "[overlap] (min_omega: number)." << std::endl; 
            return false;
        }
        // max omega
        if (comm_ol.find(schema::json_comm_ol_max_omega) == comm_ol.end()) {
            std::cerr << info << "[overlap] (max_omega)." << std::endl; 
            return false;
        } else if (!comm_ol[schema::json_comm_ol_max_omega].is_number()) {
            std::cerr << info << "[overlap] (max_omega: number)." << std::endl; 
            return false;
        }
    }
    return true;
}

bool ConfChecker::checkJsonDegr(nlohmann::json& dist, std::string info) {
    if (!dist.is_object()) {
        std::cerr << info << " (object)." << std::endl;
        return false;
    }
    if (dist.find(schema::json_dist_type) == dist.end()) {
        std::cerr << info << " (type)." << std::endl;
        return false;
    }
    if (!dist[schema::json_dist_type].is_string()) {
        std::cerr << info << " (type: string)." << std::endl;
        return false;
    }
    std::string dist_type = dist[schema::json_dist_type];
    if (dist.find(schema::json_ted_min_degree) == dist.end()) {
        std::cerr << info << " (min_d)." << std::endl;
        return false;
    }
    if (!dist[schema::json_ted_min_degree].is_number()) {
        std::cerr << info << " (min_d: number)." << std::endl;
        return false;
    }
    if (dist.find(schema::json_ted_max_degree) == dist.end()) {
        std::cerr << info << " (max_d)." << std::endl;
        return false;
    }
    if (!dist[schema::json_ted_max_degree].is_number()) {
        std::cerr << info << " (max_d: number)." << std::endl;
        return false;
    }
    if (dist_type == schema::json_dist_PowerLaw || dist_type == schema::json_dist_Exponential) {
        if (dist.find(schema::json_dist_lambda) == dist.end()) {
            std::cerr << info << " (power_law lambda)." << std::endl;
            return false;
        }
        if (!dist[schema::json_dist_lambda].is_number()) {
            std::cerr << info << " (power_law lambda: number)." << std::endl;
            return false;
        }
    } else if (dist_type == schema::json_dist_Normal || dist_type == schema::json_dist_LogNormal) {
        if (dist.find(schema::json_dist_mu) == dist.end()) {
            std::cerr << info << " ([log]normal mu)." << std::endl;
            return false;
        }
        if (!dist[schema::json_dist_mu].is_number()) {
            std::cerr << info << " ([log]normal mu: number)." << std::endl;
            return false;
        }
        if (dist.find(schema::json_dist_sigma) == dist.end()) {
            std::cerr << info << " ([log]normal sigma)." << std::endl;
            return false;
        }
        if (!dist[schema::json_dist_sigma].is_number()) {
            std::cerr << info << " ([log]normal sigma: number)." << std::endl;
            return false;
        }
    } else if (dist_type == schema::json_dist_Uniform) {
        //
    } else {
        std::cerr << info << " (unknown distribution or not implement: " << dist_type << ")." << std::endl;
        return false;
    }
    return true;
}



ProgressBar::ProgressBar() {
    mi_bar_width = 70;
    md_progress = 0.0;
}

ProgressBar::~ProgressBar() {

}

void ProgressBar::setProgress(double pg) {
    md_progress = pg;
    show();
}

void ProgressBar::show() {
    md_progress = std::min(md_progress, 1.0);
    std::cout << "[";
    int pos = mi_bar_width * md_progress;
    for (int i = 0; i < pos; ++i)
        std::cout << "=";
    std::cout << ">";
    for (int i = pos + 1; i < mi_bar_width; ++i)
        std::cout << " ";
    std::cout << "] " << (int)(md_progress * 100.0) << " %\r";
    std::cout.flush();
}

} //! namespace gtb
} //! namespace gl