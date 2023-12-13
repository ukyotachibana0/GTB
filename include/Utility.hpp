/**
 * edit time : 2023-12-12
 * Utility functions
 */
#pragma once
#include "headers.hpp"
#include "Random.hpp"

namespace gl {
namespace gtb {

class Utility
{
public:
    static double mathCeil(double x);

    static double mathRound(double x);

    static double mathFloor(double x);

    static double normPdf(double x);

    static double normCdf(double x);

    static void show(const std::vector<int_t>& list);

    static void show(const std::vector<std::vector<int_t>>& list);

    static int_t min(const std::vector<int_t>& list);

    static int_t max(const std::vector<int_t>& list);

    static int_t max(const std::vector<std::vector<int_t>> list, int col);

    static int_t max2(const std::vector<int_t>& list);

    static std::vector<int_t> getColumn(const std::vector<std::vector<int_t>>& list, int col);

    static int_t max_diff(const std::vector<std::vector<int_t>> list);

    static int numOneBitInt(uint32_t x);

    static std::string strSeqNumber(int x, int n_bits=5);

    static void randomChoice(std::vector<int>& nums, int K);

    static void randomChoice(std::vector<int>& nums, std::vector<int>& accomp, int K);

}; //! class Utility

class ConfChecker
{
public:
    static bool checkJson(nlohmann::json& json_obj);
private:
    static bool checkJsonTemp(nlohmann::json& temp, std::string info);
    static bool checkJsonComm(nlohmann::json& comm, std::string info, bool required);
    static bool checkJsonDegr(nlohmann::json& degr, std::string info);
}; //! class ConfChecker

class ProgressBar
{
public:
    ProgressBar();
    ~ProgressBar();

    void setProgress(double pg);

private:
    int mi_bar_width;
    double md_progress;

    void show();

}; //! class ProgressBar

} //! namespace gtb
} //! namespace gl