/**
 * edit time : 2023-12-12
 * Store header
 */
#pragma once

#include <set>
#include "headers.hpp"
#include "BufferWriter.hpp"

namespace gl {
namespace gtb {

class Store
{
private:
    BufferWriter *bw;
    BufferWriter *bw_csroff;
    
    std::string filename;  // base filename
    std::string format;
    int buffer_size;

    int_t file_lines_upper;     // #Lines upper in a file
    int_t current_file_lines;   // #Lines in this file
    int_t current_file_no;      // #FileNo of current file
    int_t offset;               // CSR offset

public:
    Store() {
        bw = bw_csroff = nullptr;

        filename = "part_";
        format = schema::json_format_TSV;
        buffer_size = 1024*1024*32;

        file_lines_upper = 1000000;
        current_file_lines = 0;
        current_file_no = 0;
        offset = 0;

        nextFile();
    }

    Store(std::string _filename, std::string _format, int _buffer_size=1024*1024*32) {
        bw = bw_csroff = nullptr;

        filename = _filename;
        format = _format;
        buffer_size = _buffer_size;

        file_lines_upper = 1000000;
        current_file_lines = 0;
        current_file_no = 1;
        offset = 0;

        nextFile();
    }

    void nextFile() {
        close();
        std::string file_no = Utility::strSeqNumber(current_file_no);
        std::string fn = filename + "_" + file_no + "." + format;
        bw = new BufferWriter(fn.c_str(), buffer_size);
        if (format == schema::json_format_CSR) {
            offset = 0;
            std::string csr_off_filename = filename + "_" + file_no + "." + schema::json_format_CSROFF;
            bw_csroff = new BufferWriter(csr_off_filename.c_str(), buffer_size);
        }
        // update
        current_file_lines = 0;
        current_file_no ++;
    }

    void writeLine(int_t i, std::unordered_set<int_t>& adj) {
        if (adj.empty()) return;

        if (format == schema::json_format_TSV) { // TSV
            for (auto n : adj) {
                bw->write(i);
                bw->tab();
                bw->write(n);
                bw->newline();
                // update #Lines
                current_file_lines ++;
                if (current_file_lines == file_lines_upper) nextFile();
            }
        } else if (format == schema::json_format_ADJ) { // ADJ
            bw->write(i);
            for (auto n : adj) {
                bw->space();
                bw->write(n);
            }
            bw->newline();
            // update #Lines
            current_file_lines ++;
            if (current_file_lines == file_lines_upper) nextFile();
        } else { //  CSR
            bw->write(i);
            bw->tab();
            for (auto n : adj) {
                bw->write(n);
                bw->tab();
            }
            // offset
            bw_csroff->write(offset);
            bw_csroff->tab();
            offset += adj.size();
        }
    }

    void writeLine(int_t i, std::set<std::pair<int_t, int_t>>& adj) {
        if (adj.empty()) return;

        if (format == schema::json_format_TSV) {           // TSV
            for (auto n : adj) {
                bw->write(i);
                bw->tab();
                bw->write(n.first);
                bw->tab();
                bw->write(n.second);

                bw->newline();
                // update #Lines
                current_file_lines++;
                if (current_file_lines == file_lines_upper) {
                    nextFile();
                }
            }
        } else if (format == schema::json_format_ADJ) {    // ADJ
            bw->write(i);
            for (auto n : adj) {
                bw->space();
                bw->write(n.first);
                bw->space();
                bw->write(n.second);
            }
            bw->newline();
            // update #Lines
            current_file_lines++;
            if (current_file_lines == file_lines_upper) {
                nextFile();
            }
        } else {                                        //  CSR
            bw->write(i);
            bw->tab();
            for (auto n : adj) {
                bw->write(n.first);
                bw->tab();
                bw->write(n.second);
                bw->tab();
            }
            // offset
            bw_csroff->write(offset);
            bw_csroff->tab();
            offset += adj.size();
        }
    }

    void writeLine(int_t i, std::vector<std::pair<int_t, int_t>>& adj) {
        if (adj.empty()) return;

        if (format == schema::json_format_TSV) {           // TSV
            for (auto n : adj) {
                bw->write(i);
                bw->tab();
                bw->write(n.first);
                bw->tab();
                bw->write(n.second);

                bw->newline();
                // update #Lines
                current_file_lines++;
                if (current_file_lines == file_lines_upper) {
                    nextFile();
                }
            }
        } else if (format == schema::json_format_ADJ) {    // ADJ
            bw->write(i);
            for (auto n : adj) {
                bw->space();
                bw->write(n.first);
                bw->space();
                bw->write(n.second);
            }
            bw->newline();
            // update #Lines
            current_file_lines++;
            if (current_file_lines == file_lines_upper) {
                nextFile();
            }
        } else {                                          //  CSR
            bw->write(i);
            bw->tab();
            for (auto n : adj) {
                bw->write(n.first);
                bw->tab();
                bw->write(n.second);
                bw->tab();
            }
            // offset
            bw_csroff->write(offset);
            bw_csroff->tab();
            offset += adj.size();
        }
    }

    void flush() {
        if (bw) {
            bw->flush();
        }
        if (bw_csroff) {
            bw_csroff->flush();
        }
    }

    void close() {
        flush();
        if (bw) {
            bw->close();
            bw = nullptr;
        }
        if (bw_csroff) {
            bw_csroff->close();
            bw_csroff = nullptr;
        }
    }

    ~Store() {
        close();
    }
}; //! class Store

} //! namespace gtb
} //! namespace gl
