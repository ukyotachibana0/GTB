/**
 * edit time : 2023-12-12
 * Schema for GTB
 */
#pragma once

namespace gl {
namespace gtb {

namespace schema {

const static std::string json_graph = "graph";
// node schema
const static std::string json_node = "node";
const static std::string json_node_label = "label";
const static std::string json_node_amount = "amount";
// edge schema
const static std::string json_edge = "edge";
const static std::string json_edge_label = "label";
const static std::string json_edge_source = "source";
const static std::string json_edge_target = "target";
const static std::string json_edge_amount = "amount";
// TED
const static std::string json_ted_in = "in";
const static std::string json_ted_out = "out";
const static std::string json_ted_temp = "ts";
const static std::string json_ted_min_timestamp = "min_ts";
const static std::string json_ted_max_timestamp = "max_ts";
const static std::string json_ted_min_degree = "min_d";
const static std::string json_ted_max_degree = "max_d";
const static std::string json_dist_type = "type";
const static std::string json_dist_PowerLaw = "power_law";
const static std::string json_dist_Exponential = "exponential";
const static std::string json_dist_Normal = "normal";   // Gaussian
const static std::string json_dist_LogNormal = "log_normal";
const static std::string json_dist_Uniform = "uniform";
const static std::string json_dist_lambda = "lambda";   // for power-law and exponential
const static std::string json_dist_mu = "mu";           // for normal and lognormal
const static std::string json_dist_sigma = "sigma";     // for normal and lognormal
// community
const static std::string json_comm = "community";
const static std::string json_comm_amount = "amount";
const static std::string json_comm_delta = "delta";
const static std::string json_comm_lambda = "lambda";
// overlap
const static std::string json_comm_ol = "overlap";
const static std::string json_comm_ol_m = "m";
const static std::string json_comm_ol_min_omega = "min_omega";
const static std::string json_comm_ol_max_omega = "max_omega";
// storage format
const static std::string json_store_format = "store_format";
const static std::string json_format_ADJ = "ADJ";
const static std::string json_format_TSV = "TSV";   // default format
const static std::string json_format_CSR = "CSR";
const static std::string json_format_CSROFF = "CSROFF"; // attached
// true or false
const static std::string json_true = "true";
const static std::string json_false = "false";
} //! namespace schema

} //! namespace gtb
} //! namespace gl
