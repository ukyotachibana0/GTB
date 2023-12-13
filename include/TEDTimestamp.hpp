/**
 * edit time : 2023-12-12
 * TEDTimestamp hpp
 */
#pragma once

#include "headers.hpp"
#include "Random.hpp"
#include "Utility.hpp"

namespace gl {
namespace gtb {

class TEDTimestamp
{
private:
    const double MIN_STEP = 0.000001;
    // pre process data
    double step;
    std::vector<int_t> CDF_inv;  // CDF^{-1}

protected:
    // basic info data
    int_t ts_range; // translate domain defination of pdf from [ts_min, ts_max] to [1, ts_range]

public:
    // basic info data
    int_t min_timestamp;
    int_t max_timestamp;
    std::unordered_map<std::string, double> params;

    Random rand;
    
private:
    int_t mapRange(int_t mit, int_t mat, int_t ts); // map timestamp in current range to another range

public:
    TEDTimestamp () : min_timestamp(0), max_timestamp(0), ts_range(0), step(1.0) {}

    TEDTimestamp (int_t mit, int_t mat, const std::unordered_map<std::string, double>& ps) 
        : min_timestamp(mit), max_timestamp(mat), ts_range(mat - mit + 1), params(ps), step(1.0) {  
        if (params.count("lambda")) params["lambda"] = -fabs(params["lambda"]);
    }

    virtual double pdf(int_t x) = 0;

    void preProcess();

    virtual int_t genTimestamp(int_t mit, int_t mat);

    int_t genTimestamp(const std::vector<std::vector<int_t>>& window);  // for uniform distribution
}; //! class TEDTimestamp

class UniformTimestamp : public TEDTimestamp {
public:
    UniformTimestamp() : TEDTimestamp() {}
    UniformTimestamp(int_t mit, int_t mat, const std::unordered_map<std::string, double>& ps) : TEDTimestamp(mit, mat, ps) {}

    double pdf(int_t x) {
        return .0;
    }

    virtual int_t genTimestamp(int_t mit, int_t mat) {
        int_t bound = mat - mit;
        return mit + rand.nextInt(bound);
    }
}; //! class UniformTimestamp

class PowerLawTimestamp : public TEDTimestamp {
public:
    PowerLawTimestamp() : TEDTimestamp() {}
    PowerLawTimestamp(int_t mit, int_t mat, const std::unordered_map<std::string, double>& ps) : TEDTimestamp(mit, mat, ps) {}

    double pdf(int_t x) {
        double lambda = params["lambda"];
        return pow((double)(x), lambda);
    }
}; //! class PowerLawTimestamp

class NormalTimestamp : public TEDTimestamp {
public:
    NormalTimestamp() : TEDTimestamp() {}
    NormalTimestamp(int_t mit, int_t mat, const std::unordered_map<std::string, double>& ps) : TEDTimestamp(mit, mat, ps) {}

    double pdf(int_t x) {
        double mu = params.count("mu") ? params["mu"] : (ts_range >> 1);
        double sigma = params["sigma"];
        double a = (x + 0.1 - mu) / sigma;
        double b = (x - 0.1 - mu) / sigma;
        return Utility::normCdf(a) - Utility::normCdf(b);
    }
}; //! class NormalTimestamp

class LogNormalTimestamp : public TEDTimestamp {
public:
    LogNormalTimestamp() : TEDTimestamp() {}
    LogNormalTimestamp(int_t mit, int_t mat, const std::unordered_map<std::string, double>& ps) : TEDTimestamp(mit, mat, ps) {}

    double pdf(int_t x) {
        double mu = params.count("mu") ? params["mu"] : (ts_range >> 1);
        double sigma = params["sigma"];
        double a = log(x)/log(10)-mu;
        return 1 / (x * log(10) * sigma * sqrt(2 * PI)) * exp(-a*a/(2*sigma*sigma));
    }
}; //! class LogNormalTimestamp

class ExponentialTimestamp : public TEDTimestamp {
public:
    ExponentialTimestamp() : TEDTimestamp() {}
    ExponentialTimestamp(int_t mid, int_t mxd, std::unordered_map<std::string, double>& ps) : TEDTimestamp(mid, mxd, ps) {}

    double pdf(int_t x) {
        double lambda = -params["lambda"];
        return lambda * exp(-lambda * x);
    }
}; //! class ExponentialDegree

} //! namespace gtb
} //! namespace gl