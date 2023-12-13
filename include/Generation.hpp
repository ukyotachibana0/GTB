/**
 * edit time : 2023-12-12
 * Generation header
 */
#pragma once

#include "headers.hpp"
#include "Random.hpp"
#include "Utility.hpp"
#include "TEDDegree.hpp"
#include "TEDTimestamp.hpp"

namespace gl {
namespace gtb {

typedef struct _edge_ground_truth {
    int_t actual_edges;
    int_t between_edges;
}St_EdgeGroundTruth;

typedef struct _generate_edge_basic {
    std::string basedir;
    std::string e_source;
    std::string e_target;
    int_t s_nodes;
    int_t t_nodes;
    int_t n_edges;
    // degree distribution
    std::string outd_type;
    std::unordered_map<std::string, double> out_params;
    std::string ind_type;
    std::unordered_map<std::string, double> in_params;
    // community
    std::string comm_type;
    std::unordered_map<std::string, double> comm_params;
    std::vector<std::vector<int_t>> comm_split;
    // temporal
    std::string temp_type;
    std::unordered_map<std::string, double> temp_params;
    std::vector<std::vector<int_t>> wind_split;
    // overlap
    std::unordered_map<std::string, double> ol_params;
    std::vector<std::unordered_map<int_t, double>> ol_comm;
    // ground truth
    St_EdgeGroundTruth ground_truth;
} St_BasicEdgeGeneration;

typedef struct _generate_edge : public _generate_edge_basic {
    std::string e_label;
    bool b_temporal;
    bool b_social;
    bool b_overlap;
} St_EdgeGeneration; //! Topology Generation

class Generation
{
private:
    std::string store_dir;
    std::string json_fn;
    JSON::json json_obj;
    std::string g_format;

    Random rand;
    std::vector<St_EdgeGeneration> edge_gen_plan;           // Generation Plan
    
    int n_threads;
    const int thread_chunk_size = 16;

    double gp_progress;     // Generation Progress
    std::string gp_tag;     // Generation Tag (Node/Edge-{name/source_name-target_name})
    bool gb_gen_done;       // is Global Generation Progress Done
    bool gb_start_gen;      // has started generating

public:
    Generation();

    Generation(std::string& filename, std::string& dirname);
    
    ~Generation();

    void run();

    void generatePlan();

    void generateGraph(St_EdgeGeneration& st_edge);

    void GTB(St_EdgeGeneration& st_edge);   // the temporal graph Generation method featuring Time-Bound communities

    void GGG(St_EdgeGeneration& st_edge);   // the General Graph Generaiton method

    TEDDegree* getDistDeg(int_t mid, int_t mxd, int_t n, int_t m,
        std::unordered_map<std::string, double>& params,
        bool out_pp, std::string& d_type);

    TEDTimestamp* getDistTs(int_t mit, int_t mat, std::unordered_map<std::string,
        double>& params, const std::string& t_type);

    std::string currentGenerationTag();
    double currentGenerationProgress();
    bool isGenerationDone();
    bool hasGenerationStart();

    int_t betweenOutDegree(int_t max_degree, double c);

    bool pathExists(std::string& path);

    bool mkdir(std::string& path);

private:
    void outputGT(St_EdgeGeneration& st_edge);

    void outputGTOverall();
}; //! class Generation

} //! namespace gtb
} //! namespace gl