/**
 * edit time : 2023-12-12
 * TEDDegree hpp
 */
#pragma once

#include "headers.hpp"
#include "Random.hpp"
#include "Utility.hpp"

namespace gl {
namespace gtb {

class TEDDegree
{
private:
    int_t min_degree;
    int_t max_degree;
    int_t degree_range;
    int_t num_nodes;
    int_t num_edges;

    // arrays and variables for getting out-degree
    int bucket; // generate out-degree randomly
    double o_min_gap;
    double bucket_step; // [0,1]
    std::vector<int_t> od_degree;

    // arrays and variables for getting target vertex ID
    std::vector<int_t> id_hash_tid; // corresponding target vertex ID
    std::vector<double> id_hash_cdf; // cumulative distribution function for sum of degrees
    std::vector<double> id_hash_ratio;   // corresponding ratio, for getting a ID in [id_hash_tid[i], id_hash_tid[i+1]]
    double i_min_gap;
    double i_min_rv;

public:
    // parameters of a distribution
    std::unordered_map<std::string, double> theta;

    // random number generator
    Random rand;

    TEDDegree();

    TEDDegree(int_t mid, int_t mxd, int_t n, int_t m,
        std::unordered_map<std::string, double>& params);

    ~TEDDegree();

    // build auxiliary function
    void preProcess(bool out);

    void setOutBucket(int x);

    // override
    virtual double pdf(int_t x) = 0;

    int_t numberOfMinDegree();

    int_t numberOfMaxDegree();

    int_t currentEdges();

    void matchEdges();

    void buildForOd();

    void buildForId();

    int_t oHF(double x);

    int_t iHF(double x);

    virtual int_t genOutDegree(int_t n, int_t id);

    virtual int_t genTargetID(int_t n);

    virtual int_t genTargetID_rv(double rv);

// ==================== For Debug ====================

    double getOMinGap();

    double getIMinGap();

    int_t getOVecLen();

    int_t getIVecLen();

    void printICdf();

    void printOd();

}; //! class TEDDegree

class UniformDegree : public TEDDegree {
private:
    bool is_special;
    int_t num_nodes;
    int_t num_edges;
    int_t min_degree;
    int_t max_degree;
    int_t frequency;
    int_t cur_degree;
    int_t cum_degree;
    int_t cum_frequency;
    int_t cur_id;

    int_t base_id;

    Random rand;

public:
    UniformDegree() : TEDDegree() {}

    UniformDegree(int_t mid, int_t mxd, int_t n, int_t m, std::unordered_map<std::string, double>& params) {
        num_nodes = n;
        num_edges = m;

        int_t avg_degre = m / n + 1;
        if (avg_degre <= mid) { mid = 1; }
        mxd = 2 * avg_degre - mid;
        min_degree = mid;
        max_degree = mxd;
        is_special = (mid == mxd);
        frequency = n / (mxd - mid + 1);
        cur_degree = mid;
        cum_degree = 0;
        cum_frequency = 0;
        cur_id = 0;
        base_id = 0;
    }

    // not used
    double pdf(int_t x) {
        return 0.0;
    }

    virtual int_t genOutDegree(int_t id) {
        if (is_special)
            return min_degree;
        int_t bound = max_degree - min_degree;
        return min_degree + rand.nextInt(bound);
    }

    virtual int_t genTargetID() {
        int_t ans = cur_id;
        cur_id += frequency;
        if (cur_id >= num_nodes) {
            cur_id = ++base_id;
            if (base_id + 1 == num_nodes) {
                base_id = 0;
            }
        }
        return ans;
    }

}; //! class UniformDegree

class PowerLawDegree : public TEDDegree {
public:
    PowerLawDegree() : TEDDegree() {}
    PowerLawDegree(int_t mid, int_t mxd, int_t n, int_t m, std::unordered_map<std::string, double>& params) : TEDDegree(mid, mxd, n, m, params) {}

    double pdf(int_t x) {
        double lambda = theta["lambda"];
        return pow((double)(x), lambda);
    }
}; //! class PowerLawDegree 

class NormalDegree : public TEDDegree {
public:
    NormalDegree() : TEDDegree() {}
    NormalDegree(int_t mid, int_t mxd, int_t n, int_t m, std::unordered_map<std::string, double>& params) : TEDDegree(mid, mxd, n, m, params) {}

    double pdf(int_t x) {
        double mu = theta["mu"];
        double sigma = theta["sigma"];
        double a = (x + 0.1 - mu) / sigma;
        double b = (x - 0.1 - mu) / sigma;
        return Utility::normCdf(a) - Utility::normCdf(b);
    }
}; //! class NormalDegree 

class LogNormalDegree : public TEDDegree {
public:
    LogNormalDegree() : TEDDegree() {}
    LogNormalDegree(int_t mid, int_t mxd, int_t n, int_t m, std::unordered_map<std::string, double>& params) : TEDDegree(mid, mxd, n, m, params) {}

    double pdf(int_t x) {
        double mu = theta["mu"];
        double sigma = theta["sigma"];
        double a = log(x)/log(10)-mu;
        return 1 / (x * log(10) * sigma * sqrt(2 * PI)) * exp(-a*a/(2*sigma*sigma));
    }

}; //! class LogNormalDegree

class ExponentialDegree : public TEDDegree {
public:
    ExponentialDegree() : TEDDegree() {}
    ExponentialDegree(int_t mid, int_t mxd, int_t n, int_t m, std::unordered_map<std::string, double>& params) : TEDDegree(mid, mxd, n, m, params) {}

    double pdf(int_t x) {
        double lambda = - theta["lambda"];
        return lambda * exp(-lambda * x);
    }
}; //! class ExponentialDegree

} //! namespace gtb
} //! namespace gl