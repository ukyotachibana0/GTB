/**
 * edit time : 2023-12-12
 * Sampler hpp
 */
#pragma once

#include "headers.hpp"
#include "Random.hpp"
#include "Utility.hpp"

namespace gl {
namespace gtb {

class Sampler
{
private:
    double epsilon;
    // pre process data
    double step = 1.0;
    std::vector<int_t> CDF_inv;  // CDF^{-1}

protected:
    // basic info data
    int_t val_range; // translate domain defination of pdf from [ts_min, ts_max] to [1, val_range]

public:
    // basic info data
    int_t val_min;
    int_t val_max;
    std::unordered_map<std::string, double> params;

    Random rand;

public:
    Sampler () : val_min(0), val_max(0), val_range(0), step(1.0) {}

    Sampler (int_t mi, int_t ma, const std::unordered_map<std::string, double>& ps)
        : val_min(mi), val_max(ma), val_range(ma - mi + 1), epsilon(1e-6), params(ps) {  
        if (params.count("lambda")) params["lambda"] = -fabs(params["lambda"]);
    }

    Sampler (int_t mi, int_t ma, double eps, const std::unordered_map<std::string, double>& ps)
        : val_min(mi), val_max(ma), val_range(ma - mi + 1), epsilon(eps), params(ps) {  
        if (params.count("lambda")) params["lambda"] = -fabs(params["lambda"]);
    }

    virtual double pdf(int_t x) = 0;

    void preProcess() {
        // compute PDF
        std::vector<double> PDF(val_range);
        double sum = .0;
        for (int i = 0; i < val_range; i++) { 
            PDF[i] = pdf(i + 1);    // [1, val_range]
            sum += PDF[i];
        }
        for (int i = 0; i < val_range; i++) {
            PDF[i] /= sum;
            step = (0 < PDF[i] && PDF[i] < step) ? PDF[i] : step;
        }
        step = (step < epsilon) ? epsilon : step;

        // compute CDF
        std::vector<double> CDF(val_range);
        std::partial_sum(PDF.begin(), PDF.end(), CDF.begin());
        
        // compute CDF_inv
        int steps = Utility::mathCeil(1.0 / step);
        CDF_inv.resize(steps);
        int_t j = 0;
        for (int i = 0; i < steps; i++) {
            double cur_step = i * step;
            while (j < val_range && CDF[j] < cur_step) { j++; }
            CDF_inv[i] = j + val_min;
        }
    }

    virtual int_t sample() {
        double y = rand.nextReal(1.0);
        return CDF_inv[y / step];
    }

    virtual int_t sample(double y) {
        assert(.0 <= y && y < 1.0);
        return CDF_inv[y / step];
    }

    int_t sample(const std::vector<std::vector<int_t>>& window) {
        int_t ans = 0;
        // gen mapped ans
        int_t bound = 0, n = window.size();
        for (auto& win : window) { bound += win[1] - win[0]; }
        int_t ans_map = rand.nextInt(bound);
        // gen original ans
        int_t last_sum = 0, sum = window[0][1] - window[0][0];
        for (int_t i = 1; i < n; i++) {
            if (ans_map <= sum) return (ans_map - last_sum) + window[i - 1][0];

            last_sum = sum;
            sum += window[i][1] - window[i][0];
        }
        return (ans_map - last_sum) + window[n - 1][0];
    }

    virtual int_t sampleLower(double p) {
        double y = rand.nextReal(p);
        std::cout << "[Sampler::sampleLower] y = " << y << std::endl;
        return CDF_inv[y / step];
    }

    virtual int_t sampleUpper(double p) {
        double y = rand.nextReal(1.0 - p) + p;
        std::cout << "[Sampler::sampleUpper] y = " << y << std::endl;
        return CDF_inv[y / step];
    }

    virtual std::vector<int_t> samplePool(int_t n) {
        std::vector<int_t> val_pool = std::vector<int_t>(n, -1);
        double step_n = 1.0 / n;
        double y = step_n / 2;
        for (int i = 0; i < n; ++i ) {
            val_pool[i] = CDF_inv[y / step];
            y += step_n;
        }
        return val_pool;
    }
}; //! class Sampler


class SamplerUniform : public Sampler {
public:
    SamplerUniform() : Sampler() {}
    SamplerUniform(int_t mi, int_t ma, const std::unordered_map<std::string, double>& ps)
        : Sampler(mi, ma, ps) {}

    double pdf(int_t x) {
        return .0;
    }

    virtual int_t sample() {
        int_t bound = val_max - val_min;
        return val_min + rand.nextInt(bound);
    }

    virtual int_t sample(double y) {
        int_t bound = y * (val_max - val_min);
        return val_min + bound;
    }

    virtual int_t sampleLower(double p) { // [mit, mit + bound]
        int_t bound = val_min + (val_max - val_min) * p;
        return val_min + rand.nextInt(bound);
    }

    virtual int_t sampleUpper(double p) { // [bound, mat]
        int_t bound = val_min + (val_max - val_min) * p;
        return val_max - rand.nextInt(val_max - bound);
    }

    virtual std::vector<int_t> samplePool(int_t n) {
        std::vector<int_t> ts_pool = std::vector<int_t>(n, -1);
        double step_n = (double)val_range / (double)n;
        double ts = val_min;
        for (int i = 0; i < n; ++i ) {
            ts_pool[i] = ts;
            ts += step_n;
        }
        return ts_pool;
    }
};


class SamplerPowerLaw : public Sampler {
public:
    SamplerPowerLaw() : Sampler() {}
    SamplerPowerLaw(int_t mi, int_t ma, const std::unordered_map<std::string, double>& ps)
        : Sampler(mi, ma, ps) {}
    SamplerPowerLaw(int_t mi, int_t ma, double eps, const std::unordered_map<std::string, double>& ps)
        : Sampler(mi, ma, eps, ps) {}

    double pdf(int_t x) {
        double lambda = params["lambda"];
        return pow((double)(x), lambda);
    }
};


class SamplerNormal : public Sampler {
public:
    SamplerNormal() : Sampler() {}
    SamplerNormal(int_t mi, int_t ma, const std::unordered_map<std::string, double>& ps)
        : Sampler(mi, ma, ps) {}
    SamplerNormal(int_t mi, int_t ma, double eps, const std::unordered_map<std::string, double>& ps)
        : Sampler(mi, ma, eps, ps) {}

    double pdf(int_t x) {
        double mu = params.count("mu") ? params["mu"] : (val_range >> 1);
        double sigma = params["sigma"];
        double a = (x + 0.1 - mu) / sigma;
        double b = (x - 0.1 - mu) / sigma;
        return Utility::normCdf(a) - Utility::normCdf(b);
    }
};


class SamplerLogNormal : public Sampler {
public:
    SamplerLogNormal() : Sampler() {}
    SamplerLogNormal(int_t mi, int_t ma, const std::unordered_map<std::string, double>& ps)
        : Sampler(mi, ma, ps) {}
    SamplerLogNormal(int_t mi, int_t ma, double eps, const std::unordered_map<std::string, double>& ps)
        : Sampler(mi, ma, eps, ps) {}

    double pdf(int_t x) {
        double mu = params.count("mu") ? params["mu"] : (val_range >> 1);
        double sigma = params["sigma"];
        double a = (log(x + 0.1) - mu) / sigma;
        double b = (log(x - 0.1) - mu) / sigma;
        return Utility::normCdf(a) - Utility::normCdf(b);
    }
};

} //! namespace gtb
} //! namespace gl