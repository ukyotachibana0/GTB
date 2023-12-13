/**
 * edit time : 2023-12-12
 * Generation implementation
 */
#include "Generation.hpp"

#include "headers.hpp"
#include "Utility.hpp"
#include "Grouping.hpp"
#include "Binding.hpp"
#include "Store.hpp"

#include <sys/types.h>
#include <sys/stat.h>
#include <cstdlib>
#include <omp.h>
#include <list>

namespace gl {
namespace gtb {

Generation::Generation() {
    gp_progress = 0.0;
    gp_tag = "";
    gb_gen_done = false;
    gb_start_gen = false;
}

Generation::Generation(std::string& filename, std::string& dirname) {
    #ifdef DEBUG
    std::cout << "[Generation::Generation] config filename: " << filename << " , store dirname: " << dirname << std::endl;
    std::cout << std::boolalpha << "[Generation::Generation] dirname exists: " << bool(pathExists(dirname)) << std::endl;
    std::cout << std::noboolalpha;
    #endif

    json_fn = filename;
    std::ifstream fin(filename.c_str());
    if (!fin.is_open()) {
        std::cerr << "[Generation::Generation] Cannot open " << filename << std::endl;
        return;
    }
    fin >> json_obj;

    #pragma omp parallel
    {
        n_threads = omp_get_num_threads();
    }
    #ifdef DEBUG
    std::cout << "#Threads = " << n_threads << std::endl;
    #endif
    
    while (dirname.back() == '/' || dirname.back() == '\\') { dirname.pop_back(); }
    store_dir = dirname;
    gp_progress = 0.0;
    gp_tag = "";
    gb_gen_done = false;
    gb_start_gen = false;
}

Generation::~Generation() {
    // pass
}

void Generation::generatePlan() {
    mkdir(store_dir);
    // graph name
    std::string graph_name = json_obj[schema::json_graph];
    store_dir += "/" + graph_name;
    mkdir(store_dir);
    // node schema
    auto& node_schema = json_obj[schema::json_node];
    std::unordered_map<std::string, int_t> node_label_amount;
    for (auto& one_node : node_schema) {
        std::string n_label = one_node[schema::json_node_label];
        int_t n_amount = one_node[schema::json_node_amount];
        node_label_amount[n_label] = n_amount;
    }
    // storage format
    g_format = schema::json_format_TSV;
    if (json_obj.find(schema::json_store_format) != json_obj.end()) {
        g_format = json_obj[schema::json_store_format];
    }
    if (g_format != schema::json_format_TSV && g_format != schema::json_format_ADJ && g_format != schema::json_format_CSR)
        g_format = schema::json_format_TSV;
    // edge schema
    auto& edge_schema = json_obj[schema::json_edge];
    for (auto& one_edge : edge_schema) {
        // construct a St_EdgeGeneration
        St_EdgeGeneration st_one_edge;

        // basic Information
        std::string e_label = one_edge[schema::json_edge_label];
        std::string e_source = one_edge[schema::json_edge_source];
        std::string e_target = one_edge[schema::json_edge_target];
        int_t s_nodes = node_label_amount[e_source];
        int_t t_nodes = node_label_amount[e_target];

        std::string sub_dir = store_dir + "/" + e_label;
        mkdir(sub_dir);

        int_t e_amount = one_edge[schema::json_edge_amount];
        // information
        st_one_edge.e_label = e_label;
        st_one_edge.e_source = e_source;
        st_one_edge.e_target = e_target;
        st_one_edge.n_edges = e_amount;
        st_one_edge.s_nodes = s_nodes;
        st_one_edge.t_nodes = t_nodes;
        st_one_edge.basedir = sub_dir + "/";

        auto& out_params = st_one_edge.out_params;
        auto& in_params = st_one_edge.in_params;
        auto& comm_params = st_one_edge.comm_params;
        auto& temp_params = st_one_edge.temp_params;
        auto& ol_params = st_one_edge.ol_params;

        // out-degree distribution configuration
        auto& out_dist = one_edge[schema::json_ted_out];
        st_one_edge.outd_type = out_dist[schema::json_dist_type];
        out_params[schema::json_ted_min_degree] = out_dist[schema::json_ted_min_degree];
        out_params[schema::json_ted_max_degree] = out_dist[schema::json_ted_max_degree];
        if (out_dist.find(schema::json_dist_lambda) != out_dist.end())
            out_params[schema::json_dist_lambda] = out_dist[schema::json_dist_lambda];
        if (out_dist.find(schema::json_dist_mu) != out_dist.end())
            out_params[schema::json_dist_mu] = out_dist[schema::json_dist_mu];
        if (out_dist.find(schema::json_dist_sigma) != out_dist.end())
            out_params[schema::json_dist_sigma] = out_dist[schema::json_dist_sigma];

        // in-degree distribution configuration
        auto& in_dist = one_edge[schema::json_ted_in];
        st_one_edge.ind_type = in_dist[schema::json_dist_type];
        in_params[schema::json_ted_min_degree] = in_dist[schema::json_ted_min_degree];
        in_params[schema::json_ted_max_degree] = in_dist[schema::json_ted_max_degree];
        if (in_dist.find(schema::json_dist_lambda) != in_dist.end())
            in_params[schema::json_dist_lambda] = in_dist[schema::json_dist_lambda];
        if (in_dist.find(schema::json_dist_mu) != in_dist.end())
            in_params[schema::json_dist_mu] = in_dist[schema::json_dist_mu];
        if (in_dist.find(schema::json_dist_sigma) != in_dist.end())
            in_params[schema::json_dist_sigma] = in_dist[schema::json_dist_sigma];

        // timestamp distribution configuration
        st_one_edge.b_temporal = false;
        if (one_edge.find(schema::json_ted_temp) != one_edge.end()) {
            auto& one_temporal = one_edge[schema::json_ted_temp];
            st_one_edge.b_temporal = true;
            st_one_edge.temp_type = one_temporal[schema::json_dist_type];
            temp_params[schema::json_ted_min_timestamp] = one_temporal[schema::json_ted_min_timestamp];
            temp_params[schema::json_ted_max_timestamp] = one_temporal[schema::json_ted_max_timestamp];
            if (one_temporal.find(schema::json_dist_lambda) != one_temporal.end())
                temp_params[schema::json_dist_lambda] = one_temporal[schema::json_dist_lambda];
            if (one_temporal.find(schema::json_dist_mu) != one_temporal.end())
                temp_params[schema::json_dist_mu] = one_temporal[schema::json_dist_mu];
            if (one_temporal.find(schema::json_dist_sigma) != one_temporal.end())
                temp_params[schema::json_dist_sigma] = one_temporal[schema::json_dist_sigma];
        }

        // community configuration (optional)
        st_one_edge.b_social = false;
        st_one_edge.b_overlap = false;
        if (one_edge.find(schema::json_comm) != one_edge.end()) {
            auto& comm = one_edge[schema::json_comm];
            comm_params[schema::json_comm_amount] = comm[schema::json_comm_amount];
            comm_params[schema::json_comm_lambda] = comm[schema::json_comm_lambda];
            comm_params[schema::json_comm_delta] = comm[schema::json_comm_delta];
            st_one_edge.b_social = true;
            // overlap
            if (comm.find(schema::json_comm_ol) != comm.end()) {
                st_one_edge.b_overlap = true;
                auto& comm_ol = comm[schema::json_comm_ol];
                ol_params[schema::json_comm_ol_m] = comm_ol[schema::json_comm_ol_m];
                ol_params[schema::json_comm_ol_min_omega] = comm_ol[schema::json_comm_ol_min_omega];
                ol_params[schema::json_comm_ol_max_omega] = comm_ol[schema::json_comm_ol_max_omega];
            } else {
                st_one_edge.b_overlap = false;
                ol_params[schema::json_comm_ol_m] = 0;
                ol_params[schema::json_comm_ol_min_omega] = 0.0;
                ol_params[schema::json_comm_ol_max_omega] = 0.0;
            }
        }

    ending:
        #ifdef DEBUG
        std::cout << std::boolalpha << "[Generation::generatePlan] " << st_one_edge.e_label << ".b_temporal: " << st_one_edge.b_temporal << std::endl;
        std::cout << std::boolalpha << "[Generation::generatePlan] " << st_one_edge.e_label << ".b_social: " << st_one_edge.b_social<< std::endl;
        std::cout << std::boolalpha << "[Generation::generatePlan] " << st_one_edge.e_label << ".b_overlap: " << st_one_edge.b_overlap << std::endl;
        std::cout << std::noboolalpha;
        #endif
        // add
        edge_gen_plan.push_back(st_one_edge);
    }
}

void Generation::run() {
    bool is_legal = ConfChecker::checkJson(json_obj);
    if (!is_legal) {
        std::cerr << "[Generation::run] JSON format is wrong." << std::endl;
        return;
    }
    gb_start_gen = true;

    generatePlan();
    
    for (auto& st_one_edge : edge_gen_plan) {
        generateGraph(st_one_edge);
        outputGT(st_one_edge);
    }
    outputGTOverall();
    gb_gen_done = true;
}

std::string Generation::currentGenerationTag() {
    return gp_tag;
}

double Generation::currentGenerationProgress() {
    return gp_progress;
}

bool Generation::isGenerationDone() {
    return gb_gen_done;
}

bool Generation::hasGenerationStart() {
    return gb_start_gen;
}

void Generation::generateGraph(St_EdgeGeneration& st_edge) {
    if (st_edge.b_temporal && st_edge.b_social) GTB(st_edge);
    else GGG(st_edge);
}

void Generation::GTB(St_EdgeGeneration& st_edge) {
    std::cout << "[Generation::GTB] " << st_edge.e_source << " -> " << st_edge.e_target << std::endl;
    // according to st_edge
    std::unordered_map<std::string, double>& out_params = st_edge.out_params;
    std::unordered_map<std::string, double>& in_params = st_edge.in_params;
    std::unordered_map<std::string, double>& temp_params = st_edge.temp_params;
    std::unordered_map<std::string, double>& comm_params = st_edge.comm_params;
    std::unordered_map<std::string, double>& ol_params = st_edge.ol_params;
    std::string& ind_type = st_edge.ind_type;
    std::string& outd_type = st_edge.outd_type;
    std::string& temp_type = st_edge.temp_type;
    int_t s_nodes = st_edge.s_nodes;
    int_t t_nodes = st_edge.t_nodes;
    int_t n_edges = st_edge.n_edges;
    std::string basename = st_edge.basedir + st_edge.e_source + "_" + st_edge.e_target;

    double comm_delta = comm_params[schema::json_comm_delta];
    
    /* Grouping */
    std::cout << "[Generation::GTB] Grouping..." << std::endl;
    Grouping grouper = Grouping();
    st_edge.comm_split = grouper.splitCommunity(s_nodes, t_nodes, 
        comm_params[schema::json_comm_amount], comm_params[schema::json_comm_lambda]);
    auto& comm_split = st_edge.comm_split;
    int_t n_comms = comm_split.size();
    st_edge.ol_comm = grouper.idenOlComm(n_comms, ol_params[schema::json_comm_ol_m], 
        ol_params[schema::json_comm_ol_min_omega], ol_params[schema::json_comm_ol_max_omega]);
    auto& ol_comm = st_edge.ol_comm;
    
    std::vector<int_t> sp_row_id(s_nodes, 0);   // which comm node i belongs to
    std::vector<int_t> cumu_row_id(s_nodes, 0); // relative pos of node i in its comm
    std::vector<int_t> cumu_col_psum(n_comms, 0);
    int_t gi = 0, sri = 0, cri = 0, psum = 0;
    for (int i = 0; i < n_comms; ++i) {
        cumu_col_psum[i] = psum;
        psum += comm_split[i][1];
        for (int_t j = 0; j < comm_split[i][0]; ++j) {
            sp_row_id[gi] = i;
            cumu_row_id[gi] = j;
            gi ++;
        }
    }
    #ifdef DEBUG
        std::cout << "[Generation::GTB] #Communities = " << n_comms << " , sizes: [";
        for (int i = 0; i < n_comms; ++i) std::cout << comm_split[i][0] << " ";
        std::cout << "]" << std::endl;
    #endif
    
    /* Binding */
    std::cout << "[Generation::GTB] Binding..." << std::endl;
    Binding binder = Binding(temp_params[schema::json_ted_min_timestamp], temp_params[schema::json_ted_max_timestamp]);
    auto comm_size_src = Utility::getColumn(comm_split, 0);
    st_edge.wind_split = binder.binding(comm_size_src);
    auto& wind_split = st_edge.wind_split;

    /* Building indexes in the TED model */
    std::cout << "[Generation::GTB] Building indexes..." << std::endl;
    // for source nodes (out-degree distribution)
    int_t od_min = out_params[schema::json_ted_min_degree];
    int_t od_max = out_params[schema::json_ted_max_degree];
    int_t sp_row_max = std::min(Utility::max(comm_split, 0), s_nodes);
    int_t edges_row_max = Utility::mathRound(sp_row_max * 1.0 / s_nodes * n_edges);
    TEDDegree* dist_row_max = getDistDeg(od_min, od_max, sp_row_max, edges_row_max, out_params, true, outd_type);
    // for target nodes (in-degree distribution)
    int_t id_min = in_params[schema::json_ted_min_degree];
    int_t id_max = in_params[schema::json_ted_max_degree];
    int_t sp_col_max = std::min(Utility::max(comm_split, 1), t_nodes);
    int_t edges_col_max = Utility::mathRound(sp_col_max * 1.0 / t_nodes * n_edges);
    TEDDegree* dist_col_max = getDistDeg(od_min, od_max, sp_col_max, edges_col_max, in_params, false, ind_type);
    // for timestamp (timestamp distribution)
    int_t longest_i = Utility::max_diff(wind_split);
    TEDTimestamp* timer_main = getDistTs(wind_split[longest_i][0], wind_split[longest_i][1], temp_params, temp_type);
    int_t ts_min = temp_params[schema::json_ted_min_timestamp];
    int_t ts_max = temp_params[schema::json_ted_max_timestamp];
    TEDTimestamp* timer_betw = getDistTs(ts_min, ts_max, temp_params, schema::json_dist_Uniform);

    /* Linking */
    std::cout << "[Generation::GTB] Linking..." << std::endl;
    gp_tag = st_edge.e_label;
    bool is_homo = (st_edge.e_source == st_edge.e_target);
    int_t actual_edges = 0;
    int_t between_edges = 0;
    #ifdef PARALLEL
    std::vector<Store*> store_list;
    for (int i = 0; i < n_threads; ++i) { store_list.push_back(new Store(basename + "_thd" + std::to_string(i), g_format)); }
    #else
    Store *store = new Store(basename, g_format);
    #endif
    // showing process
    double cur = 0.0;
    double progress = 0.0;
    ProgressBar progress_bar;
    // linking start
    #ifdef PARALLEL
    #pragma omp parallel for schedule (dynamic, thread_chunk_size)
    #endif
    for (int_t i = 0; i < s_nodes; ++i) {
        Store *store_ptr = nullptr;
        #ifdef PARALLEL
        int tid = omp_get_thread_num();
        store_ptr = store_list[tid];
        #else
        store_ptr = store;
        #endif
        // unique to each thread
        int_t cumu_row_i = cumu_row_id[i];
        int_t sp_row_i = sp_row_id[i];
        int_t size_src_i = comm_split[sp_row_i][0];
        int_t size_trg_i = comm_split[sp_row_i][1];
        int_t mit_i = wind_split[sp_row_i][0];
        int_t mat_i = wind_split[sp_row_i][1];
        int_t out_degree_main = dist_row_max->genOutDegree(size_src_i, cumu_row_i);
        int_t out_degree_betw = (rand.nextReal() < comm_delta) ? 
            betweenOutDegree(od_max - out_degree_main + 10, comm_delta + 1.0) : 0;

        // i's communities (i: the source node)
        std::set<int_t> comms_i({ sp_row_i });
        for (auto& olc : ol_comm[sp_row_i]) {
            double ol = olc.second;
            int_t thre_row_i = (int_t)(size_src_i * ol);
            if (cumu_row_i > thre_row_i) comms_i.insert(olc.first);
        }

        int_t cumu_col = 0;
        std::vector<std::pair<int_t, int_t>> nbrs;
        for (int_t sp_col_j = 0; sp_col_j < n_comms; sp_col_j++) {
            int_t size_src_j = comm_split[sp_col_j][0], size_trg_j = comm_split[sp_col_j][1];
            int_t mit_j = wind_split[sp_col_j][0], mat_j = wind_split[sp_col_j][1];
            int_t out_degree_betw_j = Utility::mathRound(out_degree_betw * 1.0 * size_trg_j / (1.0 * t_nodes));
            /* linking within */
            if (sp_col_j == sp_row_i) {
                for (int_t e = 0; e < out_degree_main; e++) {
                    int_t j = dist_col_max->genTargetID(size_trg_j);
                    // time window of community `sp_col_i`
                    int_t ts = timer_main->genTimestamp(mit_j, mat_j);
                    nbrs.push_back({ j + cumu_col, ts });
                }
            } else { // overlapping
                bool b_overlap_ij = ol_comm[sp_row_i].find(sp_col_j) != ol_comm[sp_row_i].end();
                if (b_overlap_ij && sp_row_i < sp_col_j) {
                    double ol = ol_comm[sp_row_i][sp_col_j];
                    int_t thre_row_i = size_src_i - (int_t)(size_src_i * ol);
                    if (cumu_row_i >= thre_row_i) {
                        // a->a
                        int_t ol_num = out_degree_main * ol;
                        int_t ol_size = (int_t)(size_trg_i * ol);
                        int_t sp_size = size_trg_i - ol_size;
                        for (int_t e = 0; e < ol_num; e++) {
                            int_t j = dist_col_max->genTargetID(ol_size);
                            // time window of community `sp_col_j`
                            int_t ts = timer_main->genTimestamp(mit_j, mat_j);
                            nbrs.push_back({ j + cumu_col_psum[sp_row_i] + sp_size, ts });
                        }

                        // a->b
                        ol_num = out_degree_main * ol; 
                        for (int_t e = 0; e < ol_num; e++) {
                            int_t j = dist_col_max->genTargetID(size_trg_j);
                            // time window of community `sp_col_j`
                            int_t ts = timer_main->genTimestamp(mit_j, mat_j);
                            nbrs.push_back({ j + cumu_col, ts });
                        }
                    }
                }
                if (b_overlap_ij && sp_row_i > sp_col_j) {
                    // b->a
                    double ol = ol_comm[sp_row_i][sp_col_j];
                    int_t ol_size = (int_t)(size_trg_j * ol);
                    int_t sp_size = size_trg_j - ol_size;
                    int_t ol_num = out_degree_main * ol;
                    for (int_t e = 0; e < ol_num; ++e) {
                        int_t j = dist_col_max->genTargetID(ol_size);
                        // time window of community `sp_col_j`
                        int_t ts = timer_main->genTimestamp(mit_j, mat_j);
                        nbrs.push_back({ j + cumu_col + sp_size, ts });
                    }
                }
            }
            int_t within_edges = nbrs.size();

            /* linking between */
            for (int_t e = 0; e < out_degree_betw_j; e++) {
                int_t j = rand.nextInt(size_trg_j - 1);
                if (sp_row_i == sp_col_j && is_homo && i == j + cumu_col) continue;
                // j's communities (j: the target node)
                std::set<int_t> comms_j({sp_col_j});
                for (auto& olc : ol_comm[sp_col_j]) {
                    double ol = olc.second;
                    int_t thre_col_j = (int_t)(size_trg_j * ol);
                    if (j > thre_col_j) comms_j.insert(olc.first);
                }
                std::set<int_t> comms_ij;   // intersection set
                set_intersection(comms_i.begin(), comms_i.end(), comms_j.begin(), comms_j.end(), 
                    inserter(comms_ij, comms_ij.begin()));
                if (comms_ij.empty()) {
                    // no common time window: assign any timestamp
                    int_t ts = timer_betw->genTimestamp(ts_min, ts_max);
                    nbrs.push_back({j + cumu_col, ts});
                } else {
                    // common time window: assign timestamp outside of window_ij
                    std::vector<std::vector<int_t>> window_ij;
                    for (auto comm : comms_ij) { window_ij.push_back(wind_split[comm]); }
                    auto window_ij_unn = Binding::unionWindow(window_ij);
                    auto window_ij_cpl = Binding::compleWindow(window_ij_unn,
                        temp_params[schema::json_ted_min_timestamp],
                        temp_params[schema::json_ted_max_timestamp]);
                    if (window_ij_cpl.empty()) continue;
                    int_t ts = timer_betw->genTimestamp(window_ij_cpl);
                    nbrs.push_back({ j + cumu_col, ts });
                }
            }

            int_t all_edges = nbrs.size();

            #ifdef PARALLEL
            #pragma omp atomic
            #endif
            actual_edges += all_edges;
            
            #ifdef PARALLEL
            #pragma omp atomic
            #endif
            between_edges += all_edges - within_edges;

            store_ptr->writeLine(i, nbrs);
            nbrs.clear();
            cumu_col += size_trg_j;
        }

        cur = (double)actual_edges / (double) n_edges;
        #ifndef PARALLEL
        if (cur - progress >= 0.01) {
            progress += 0.01;
            gp_progress = progress;
            progress_bar.setProgress(progress);
        }
        #endif
    }
    // linking end

    // deconstruct
    #ifdef PARALLEL
    for (int i = 0; i < n_threads; ++i) { store_list[i]->close(); }
    #else
    store->close();
    progress_bar.setProgress(1.0);
    std::cout << std::endl;
    #endif

    gp_progress = 1.0;

    auto& ground_truth = st_edge.ground_truth;
    ground_truth.actual_edges = actual_edges;
    ground_truth.between_edges = between_edges;

    std::cout << "[Generation::GTB] #Source Nodes = " << s_nodes << " , #Target Nodes = " << t_nodes << std::endl;
    std::cout << "[Generation::GTB] #Actual Edges = " << actual_edges << std::endl;
    std::cout << "[Generation::GTB] #Extra Edges = " << between_edges << std::endl;
    std::cout << "[Generation::GTB] #Expect Edges = " << n_edges << std::endl;
}

void Generation::GGG(St_EdgeGeneration& st_edge) {
    std::cout << "[Generation::GGG] " << st_edge.e_source << " -> " << st_edge.e_target << std::endl;
    // according to st_edge
    std::unordered_map<std::string, double>& out_params = st_edge.out_params;
    std::unordered_map<std::string, double>& in_params = st_edge.in_params;
    std::string& ind_type = st_edge.ind_type;
    std::string& outd_type = st_edge.outd_type;
    int_t s_nodes = st_edge.s_nodes;
    int_t t_nodes = st_edge.t_nodes;
    int_t n_edges = st_edge.n_edges;
    std::string basename = st_edge.basedir + st_edge.e_source + "_" + st_edge.e_target;
    
    /* Building indexes */
    std::cout << "[Generation::GGG] Building indexes..." << std::endl;
    int_t id_min = in_params[schema::json_ted_min_degree];
    int_t id_max = in_params[schema::json_ted_max_degree];
    int_t od_min = out_params[schema::json_ted_min_degree];
    int_t od_max = out_params[schema::json_ted_max_degree];
    TEDDegree *out_dist = getDistDeg(od_min, od_max, s_nodes, n_edges, out_params, true, outd_type);
    TEDDegree *in_dist = getDistDeg(id_min, id_max, t_nodes, n_edges, in_params, false, ind_type);

    /* Linking */
    std::cout << "[Generation::GGG] Linking..." << std::endl;
    gp_tag = st_edge.e_label;
    bool is_homo = (st_edge.e_source == st_edge.e_target);
    int_t actual_edges = 0;

    #ifdef PARALLEL
    std::vector<Store*> store_list;
    for (int i = 0; i < n_threads; ++i) { store_list.push_back(new Store(basename + "_thd" + std::to_string(i), g_format)); }
    #else
    Store *store = new Store(basename, g_format);
    #endif
    // show process
    double cur = 0.0;
    double progress = 0.0;
    ProgressBar progress_bar;
    // linking start
    #ifdef PARALLEL
    #pragma omp parallel for schedule (dynamic, thread_chunk_size)
    #endif
    for (int_t i = 0; i < s_nodes; ++i) {
        Store *store_ptr = nullptr;
        #ifdef PARALLEL
        int tid = omp_get_thread_num();
        store_ptr = store_list[tid];
        #else
        store_ptr = store;
        #endif
        // unique to each thread
        std::unordered_set<int_t> nbrs;
        int_t out_degree = out_dist->genOutDegree(s_nodes, i);
        for (int_t e = 0; e < out_degree; ++e) {
            int_t j = in_dist->genTargetID(t_nodes);
            while (is_homo && j == i) { j = in_dist->genTargetID(t_nodes); }
            nbrs.insert(j);
        }

        #ifdef PARALLEL
        #pragma omp atomic
        #endif
        actual_edges += nbrs.size();

        store_ptr->writeLine(i, nbrs);
        nbrs.clear();

        cur = (double)actual_edges / (double)n_edges;
        #ifndef PARALLEL
        if (cur - progress >= 0.01) {
            progress += 0.01;
            gp_progress = progress;
            progress_bar.setProgress(progress);
        }
        #endif
    }
    // linking end

    // deconstruct
    #ifdef PARALLEL
    for (int i = 0; i < n_threads; ++i) { store_list[i]->close(); }
    #else
    store->close();
    progress_bar.setProgress(1.0);
    std::cout << std::endl;
    #endif

    gp_progress = 1.0;

    auto& ground_truth = st_edge.ground_truth;
    ground_truth.actual_edges = actual_edges;

    std::cout << "[Generation::GGG] #Source Nodes = " << s_nodes << ", #Target Nodes = " << t_nodes << std::endl;
    std::cout << "[Generation::GGG] #Actual Edges = " << actual_edges << std::endl;
    std::cout << "[Generation::GGG] #Expect Edges = " << n_edges << std::endl;
}

TEDDegree* Generation::getDistDeg(int_t mid, int_t mxd, int_t n, int_t m,
        std::unordered_map<std::string, double>& params, bool out_pp, std::string& d_type) {
    TEDDegree* ans = nullptr;
    if (d_type == schema::json_dist_PowerLaw) {
        ans = new PowerLawDegree(mid, mxd, n, m, params);
        ans->preProcess(out_pp);
    } else if (d_type == schema::json_dist_Exponential) {
        ans = new ExponentialDegree(mid, mxd, n, m, params);
        ans->preProcess(out_pp);
    } else if (d_type == schema::json_dist_Normal) {
        ans = new NormalDegree(mid, mxd, n, m, params);
        ans->preProcess(out_pp);
    } else if (d_type == schema::json_dist_LogNormal) {
        ans = new LogNormalDegree(mid, mxd, n, m, params);
        ans->preProcess(out_pp);
    } else if (d_type == schema::json_dist_Uniform) {
        ans = new UniformDegree(mid, mxd, n, m, params);
    } else {
        std::cerr << "[Generation::getDistDeg] Unknown distribution: " << d_type << std::endl;
    }
    return ans;
}

TEDTimestamp* Generation::getDistTs(int_t mit, int_t mat, 
    std::unordered_map<std::string, double>& params, const std::string& t_type) {
    TEDTimestamp* ans = nullptr;
    if (t_type == schema::json_dist_PowerLaw) {
        ans = new PowerLawTimestamp(mit, mat, params);
        ans->preProcess();
    } else if (t_type == schema::json_dist_Exponential) {
        ans = new ExponentialTimestamp(mit, mat, params);
        ans->preProcess();
    } else if (t_type == schema::json_dist_Normal) {
        ans = new NormalTimestamp(mit, mat, params);
        ans->preProcess();
    } else if (t_type == schema::json_dist_LogNormal) {
        ans = new LogNormalTimestamp(mit, mat, params);
        ans->preProcess();
    } else if (t_type == schema::json_dist_Uniform) {
        ans = new UniformTimestamp(mit, mat, params);
    } else { std::cerr << "[Generation::getDistTs] Unknown distribution: " << t_type << std::endl; }
    return ans;
}

int_t Generation::betweenOutDegree(int_t max_degree, double c) {
    double a = exp(-1.0 / c);
    double b = a - exp(-max_degree * 1.0 / c);
    double y = rand.nextReal();
    int_t ans = (int_t)(-c * log(a - y * b));
    if (ans < 0) {
        ans = 0;
    }
    return ans;
}

bool Generation::pathExists(std::string& path) {
    struct stat info;
    if (stat(path.c_str(), &info))
        return false;
    else if (info.st_mode & S_IFDIR)
        return true;
    else
         return false;
}

bool Generation::mkdir(std::string& path) {
    if (pathExists(path)) {
        #ifdef DEBUG
        std::cout << "[Generation::mkdir] \"" << path << "\" already exists." << std::endl;
        #endif
        return true;
    }
    
    std::string pathname;
    std::string md_cmd = "mkdir ";
    #ifdef _WIN32
    for (char ch : path) {
        if (ch == '/')
            pathname.push_back('\\');
        else
            pathname.push_back(ch);
    }
    md_cmd += pathname;
    #else
    md_cmd += " -p " + path;
    #endif

    int res = system(md_cmd.c_str());

    #ifdef DEBUG
    if (res == 0) std::cout << "[Generation::mkdir] Create " << path << " successfully." << std::endl;
    else std::cout << "[Generation::mkdir] Create " << path << " failed." << std::endl;
    #endif

    return res == 0;
}

void Generation::outputGT(St_EdgeGeneration& st_edge) {
    JSON::json gt;

    bool b_homo = st_edge.e_source == st_edge.e_target;
    gt["is_homo"] = b_homo;
    gt["is_temporal"] = st_edge.b_temporal && st_edge.b_social;
    gt["is_social"] = st_edge.b_temporal && st_edge.b_social;
    gt["lbl_source"] = st_edge.e_source;
    gt["lbl_target"] = st_edge.e_target;
    gt["num_source"] = st_edge.s_nodes;
    gt["num_target"] = st_edge.t_nodes;
    gt["num_edge"] = st_edge.ground_truth.actual_edges;

    if (gt["is_temporal"])
        gt["window"] = {
            int(st_edge.temp_params[schema::json_ted_min_timestamp]), 
            int(st_edge.temp_params[schema::json_ted_max_timestamp])
        };

    auto& comm_split = st_edge.comm_split;
    auto& ol_comm = st_edge.ol_comm;
    int_t n_comms = comm_split.size(), psum_s = 0, psum_t = 0;
    std::vector<std::vector<int_t>> comm_split_psum(n_comms + 1, std::vector<int_t>(2));

    if (!gt["is_social"]) goto ending;

    for (int_t i = 0; i < n_comms; ++i) {
        comm_split_psum[i][0] = psum_s;
        comm_split_psum[i][1] = psum_t;
        psum_s += comm_split[i][0];
        psum_t += comm_split[i][1];
    }
    comm_split_psum[n_comms][0] = psum_s;
    comm_split_psum[n_comms][1] = psum_t;

    gt["num_edge_between"] = st_edge.ground_truth.between_edges;
    gt["num_comm"] = n_comms;
    gt["community"] = JSON::json::array();
    for (int_t i = 0; i < n_comms; ++i) {
        JSON::json comm;

        std::vector<std::vector<int_t>> node_homo;
        std::vector<std::vector<std::vector<int_t>>> node_hete;
        // temporal
        if (gt["is_temporal"]) comm["window"] = { st_edge.wind_split[i][0], st_edge.wind_split[i][1] };
        // overlap
        if (!ol_comm[i].empty()) {
            comm["overlap"] = ol_comm[i];
            if (b_homo) node_homo = Grouping::homoOlRange(comm_split_psum, ol_comm[i], i);
            else node_hete = Grouping::heteOlRange(comm_split_psum, ol_comm[i], i);
        }
        // node
        if (b_homo) {
            node_homo.push_back({ comm_split_psum[i][0], comm_split_psum[i + 1][0] - 1 });
            comm["node"] = node_homo;
        } else {
            node_hete.push_back({
                {comm_split_psum[i][0], comm_split_psum[i][1]},
                {comm_split_psum[i + 1][0] - 1, comm_split_psum[i + 1][1] - 1}
            });
            comm["node"] = node_hete;
        }
        
        gt["community"].push_back(comm);
    }

ending:
    std::string filename = st_edge.basedir + "ground_truth.json";
    std::ofstream of(filename);
    of << gt.dump(4);
    of.close();
    #ifdef DEBUG
    std::cout << "[Generation::outputGT] End." << std::endl;
    #endif
}

void Generation::outputGTOverall() {
    JSON::json gt;
    int_t num_node = 0, num_edge = 0, num_comm = 0;
    std::vector<std::string> lbl_node, lbl_edge, lbl_comm;
    // node
    auto& node_schema = json_obj[schema::json_node];
    for (auto& one_node : node_schema) {
        num_node += one_node[schema::json_node_amount].get<int>();
        lbl_node.push_back(one_node[schema::json_node_label]);
    }
    // edge
    for (auto& one_edge : edge_gen_plan) {
        num_edge += one_edge.ground_truth.actual_edges;
        lbl_edge.push_back(one_edge.e_label);
        // comm
        if (one_edge.b_social && one_edge.b_temporal) {
            num_comm += one_edge.comm_split.size();
            lbl_comm.push_back(one_edge.e_label);
        }
    }
    gt["num_node"] = num_node;
    gt["num_edge"] = num_edge;
    gt["lbl_node"] = lbl_node;
    gt["lbl_edge"] = lbl_edge;
    if (num_comm > 0 && lbl_comm.size() > 0) {
        gt["num_comm"] = num_comm;
        gt["lbl_comm"] = lbl_comm;
    }

    std::string filename = store_dir + "/ground_truth.json";
    std::ofstream of(filename);
    of << gt.dump(4);
    of.close();
    #ifdef DEBUG
    std::cout << "[Generation::outputGTOverall] End." << std::endl;
    #endif
}

} //! namespace gtb
} //! namespace gl
