/**
 * edit time : 2023-12-12
 * TEDDegree Implementation
 */
#include "TEDDegree.hpp"

namespace gl {
namespace gtb {

TEDDegree::TEDDegree() {
    min_degree = 0;
    max_degree = 0;
    num_nodes = 0;
    num_edges = 0;
    degree_range = 0;
    bucket = 1;
    bucket_step = 1.0 / num_nodes;
}

TEDDegree::TEDDegree(int_t mid, int_t mxd, int_t n, int_t m,
    std::unordered_map<std::string, double>& params) {
    min_degree = mid;
    max_degree = mxd;
    if (min_degree < 1) {
        min_degree = 1;
    }
    if (max_degree < 1) {
        max_degree = 1;
    }
    num_nodes = n;
    num_edges = m;
    degree_range = max_degree - min_degree + 1;
    // for power-law distribution
    if (params.count("lambda")) {
        params["lambda"] = -fabs(params["lambda"]);
    }
    theta = params;
    // bucket
    bucket = 1;
    bucket_step = 1.0 / (num_nodes * 1.0);
}

TEDDegree::~TEDDegree() {

}

// build auxiliary function
void TEDDegree::preProcess(bool out) {

    #ifdef DEBUG
    std::cout << "[TEDDegree::preProcess] Matching edges..." << std::endl;
    #endif

    if (out) {
        #ifdef DEBUG
        std::cout << "[TEDDegree::preProcess] Building for source nodes..." << std::endl;
        #endif
        matchEdges();
        buildForOd();
        #ifdef DEBUG
        std::cout << "[TEDDegree::preProcess] Building for source nodes: done." << std::endl;
        // printOd();
        #endif
    } else {
        #ifdef DEBUG
        std::cout << "[TEDDegree::preProcess] Building for target nodes..." << std::endl;
        #endif
        buildForId();
        #ifdef DEBUG
        std::cout << "[TEDDegree::preProcess] Building for target nodes: done." << std::endl;
        // printICdf();
        #endif
    }
}

void TEDDegree::setOutBucket(int x) {
    bucket = x;
    bucket_step = x * 1.0 / num_nodes;
}

// the number of nodes whose degree could be assigned to the min degree
int_t TEDDegree::numberOfMinDegree() {
    double min_pdf = pdf(min_degree);
    double sum_pdf = min_pdf;
    for (int_t i = min_degree + 1; i <= max_degree; ++i)
        sum_pdf += pdf(i);
    return Utility::mathRound(num_nodes * min_pdf / sum_pdf);
}

// the number of nodes whose degree could be assigned to the max degree
int_t TEDDegree::numberOfMaxDegree() {
    double max_pdf = pdf(max_degree);
    double sum_pdf = max_pdf;
    for (int_t i = min_degree; i < max_degree; ++i)
        sum_pdf += pdf(i);
    return Utility::mathRound(num_nodes * max_pdf / sum_pdf);
}

int_t TEDDegree::currentEdges() {
    degree_range = max_degree - min_degree + 1;
    double sum = 0.0;
    std::vector<double> p((int)(degree_range));
    for (int_t i = min_degree; i <= max_degree; ++i) {
        p[i - min_degree] = pdf(i);
        sum += p[i - min_degree];
    }
    double alpha = num_nodes * 1.0 / sum;
    int_t ans = 0;
    for (int i = 0; i < degree_range; ++i)
        ans += Utility::mathRound(alpha * p[i] * (i + min_degree));
    return ans;
}

void TEDDegree::matchEdges() {
    // adjust max_degree so that #nodes(degree = max_degree) >= 1
    max_degree = std::min(max_degree, num_nodes);
    int_t bin_l = max_degree + 1;
    int_t bin_r = max_degree;
    while (bin_r > 0 && numberOfMaxDegree() < 1) {
        bin_r = max_degree;
        max_degree /= 2;
        bin_l = max_degree;
    }
    while (bin_l < bin_r) {
        max_degree = bin_l + (bin_r - bin_l) / 2;
        if (numberOfMaxDegree() < 1)
            bin_r = max_degree;
        else
            bin_l = max_degree + 1;
    }
    max_degree = bin_l - 1;

    int_t actual_edges = currentEdges();

    // => expected #edges
    if (actual_edges < num_edges) {
        int_t temp = max_degree;
        bin_l = bin_r = max_degree + 1;
        while (max_degree < num_nodes && numberOfMaxDegree() > 0) {
            bin_l = max_degree;
            max_degree *= 2;
            bin_r = max_degree;
        }
        bin_r = std::min(bin_r, num_nodes);
        while (bin_l < bin_r) {
            max_degree = bin_l + (bin_r - bin_l) / 2;
            if (numberOfMaxDegree() < 1)
                bin_r = max_degree;
            else
                bin_l = max_degree + 1;
        }
        max_degree = bin_l - 1;
        actual_edges = currentEdges();
        if (actual_edges > num_edges) {
            bin_l = temp;
            bin_r = max_degree;
            while (bin_l < bin_r) {
                max_degree = bin_l + (bin_r - bin_l) / 2;
                actual_edges = currentEdges();
                if (actual_edges > num_edges)
                    bin_r = max_degree;
                else
                    bin_l = max_degree + 1;
            }
            max_degree = bin_l - 1;
            actual_edges = currentEdges();
        }
    } else if (actual_edges > num_edges) {
        bin_l = min_degree;
        bin_r = max_degree;
        while (bin_l < bin_r) {
            max_degree = bin_l + (bin_r - bin_l) / 2;
            actual_edges = currentEdges();
            if (actual_edges > num_edges)
                bin_r = max_degree;
            else
                bin_l = max_degree + 1;
        }
        max_degree = bin_l - 1;
        actual_edges = currentEdges();
    }

    degree_range = max_degree - min_degree + 1;
}

void TEDDegree::buildForOd() {
    // get CDF of out-degree distribution
    int size = (int)degree_range;
    std::vector<double> cdf(size);
    double sum = 0.0;
    double min_gap = 1.0;
    double pre = 0.0;
    for (int i = min_degree; i <= max_degree; ++i) {
        sum += pdf(i);
        cdf[i - min_degree] = sum;
    }
    double alpha = 1.0 / sum;
    for (int i = 0; i < size; ++i) {
        cdf[i] *= alpha;
        min_gap = std::min(min_gap, cdf[i] - pre);
        pre = cdf[i];
    }
    // std::cout << "[TEDDegree::buildForOd] min_gap: " << min_gap << std::endl;
    min_gap = std::max(min_gap, 0.00001);

    int length = (int)Utility::mathCeil(1.0 / min_gap) + 1;
    od_degree.resize(length);
    pre = 0.0;
    int j = 0;
    int_t d = min_degree;

    // build
    for (int i = 0; i < degree_range; ++i) {
        int steps = (int)((cdf[i] - pre) / min_gap);
        for (int p = 0; p < steps; ++p)
            od_degree[j ++] = d;
        pre += steps * min_gap;
        while (pre < cdf[i]) {
            od_degree[j ++] = d;
            pre += min_gap;
        }
        d = i + min_degree;
    }
    while (j < length)
        od_degree[j ++] = d;
    o_min_gap = min_gap;
}

void TEDDegree::buildForId() {
    int size = (int)degree_range;
    std::vector<double> temp(size);
    double cum = 0.0;
    for (int i = min_degree; i <= max_degree; ++i) {
        cum += pdf(i);
        temp[i - min_degree] = cum;
    }
    double alpha = num_nodes / cum;
    for (int i = 0; i < degree_range; ++i) {
        temp[i] *= alpha;
    }

    // build CDF corresponding to degree
    std::vector<double> i_CDF(size + 1);
    std::vector<int_t> i_index(size + 1);
    i_CDF[0] = 0.0;
    i_index[0] = -1;
    int_t pre_num = 0;
    double deg_sum = 0;
    for (int i = 1; i <= degree_range; ++i) {
        int_t num = (int_t)Utility::mathRound(temp[i - 1]);
        i_index[i] = num;
        deg_sum += (num - pre_num) * (i - 1 + min_degree);
        i_CDF[i] = deg_sum;
        pre_num = num;
    }

    i_min_rv = min_degree / i_CDF[size];
    i_CDF[0] = 0.0;
    i_index[0] = 0;
    i_index[size] = num_nodes - 1;

    i_min_gap = 1.0;
    for (int i = 1; i <= degree_range; ++i) {
        i_CDF[i] /= i_CDF[size];
        i_min_gap = std::min(i_min_gap, i_CDF[i] - i_CDF[i - 1]);
    }
    i_min_gap = std::max(i_min_gap, 0.001);
    
    // build hash tables
    int length = (int)Utility::mathCeil(1.0 / i_min_gap) + 1;
    id_hash_tid.resize(length);
    id_hash_cdf.resize(length);
    id_hash_ratio.resize(length);
    int_t i_tid_next;
    double i_cdf_next;
    int_t pre_v = 0;
    int pre_i = 0;
    for (int j = 1; j <= degree_range; ++j) {
        int i = (int)iHF(i_CDF[j]);
        i_tid_next = i_index[j];
        i_cdf_next = i_CDF[j];
        for (int k = pre_i; k < i; ++k) {
            id_hash_tid[k] = pre_v;
            id_hash_cdf[k] = i_CDF[j - 1];
            if ((i_cdf_next - id_hash_cdf[k]) < 1e-40)
                id_hash_ratio[k] = 0.0;
            else
                id_hash_ratio[k] = (i_tid_next - id_hash_tid[k]) / (i_cdf_next - id_hash_cdf[k]);
        }
        pre_v = i_index[j];
        pre_i = i;
    }
    i_tid_next = num_nodes - 1;
    i_cdf_next = 1.0;
    for (int k = pre_i; k < length; ++k) {
        id_hash_tid[k] = pre_v;
        id_hash_cdf[k] = i_CDF[size];
        if ((i_cdf_next - id_hash_cdf[k]) < 1e-40)
            id_hash_ratio[k] = 0.0;
        else
            id_hash_ratio[k] = (i_tid_next - id_hash_tid[k]) / (i_cdf_next - id_hash_cdf[k]);
    }
}

int_t TEDDegree::oHF(double x) {
    return std::min(Utility::mathRound(x / o_min_gap), (double)(od_degree.size() - 1));
}

int_t TEDDegree::iHF(double x) {
    return std::min(Utility::mathRound(x / i_min_gap), (double)(id_hash_tid.size() - 1));
}

int_t TEDDegree::genOutDegree(int_t n, int_t id) {
    double rv = (double)id / (double)n + rand.nextReal() * bucket_step;
    int index = (int)oHF(rv);
    return od_degree[index];
}

int_t TEDDegree::genTargetID_rv(double rv) {
    int index = (int)iHF(rv);
    int_t a;
    double c, r;
    if (rv >= id_hash_cdf[index]) {
        a = id_hash_tid[index];
        c = id_hash_cdf[index];
        r = id_hash_ratio[index];
    } else {
        a = id_hash_tid[index - 1];
        c = id_hash_cdf[index - 1];
        r = id_hash_ratio[index - 1];
    }
    int_t tid = a + Utility::mathRound((rv - c) * r);
    return tid;
}

int_t TEDDegree::genTargetID(int_t n) {
    // find corresponding position for n
    assert(n <= num_nodes);
    double rv = std::max(rand.nextReal(), i_min_rv);
    int_t tid = genTargetID_rv(rv);
    // transfer node ID in current range to another range
    tid = (double)tid * ((double)n / (double)num_nodes);
    return tid;
}

// =================== For Debug ====================

double TEDDegree::getOMinGap() {
    return o_min_gap;
}

double TEDDegree::getIMinGap() {
    return i_min_gap;
}

int_t TEDDegree::getOVecLen() {
    return od_degree.size();
}

int_t TEDDegree::getIVecLen() {
    return id_hash_tid.size();
}

void TEDDegree::printOd() {
    std::cout << "[TEDDegree::printOd] od_degree: ["; 
    for (auto n : od_degree)
        std::cout << n << ", ";
    std::cout << "]" << std::endl;
}

void TEDDegree::printICdf() {
    std::cout << "[TEDDegree::printICdf] len(id_hash_cdf)=" << id_hash_cdf.size() << ", num_nodes=" << num_nodes << std::endl;
    std::cout << "[TEDDegree::printICdf] id_hash_cdf: ["; 
    for (auto n : id_hash_cdf)
        std::cout << n << ", ";
    std::cout << "]" << std::endl;

    std::cout << "[TEDDegree::printICdf] id_hash_tid: ["; 
    for (auto n : id_hash_tid)
        std::cout << n << ", ";
    std::cout << "]" << std::endl;
}

} //! namespace gtb
} //! namespace gl