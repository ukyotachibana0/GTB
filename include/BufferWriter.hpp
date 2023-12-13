/**
 * edit time : 2023-12-12
 * Buffer Writer
 */

#pragma once

#include <cstdio>
#include <string>
#include "types.hpp"

namespace gl {
namespace gtb {

class BufferWriter
{
public:
    BufferWriter();
    
    BufferWriter(const char* filename, int buffer_size=32 * 1024 * 1024);

    BufferWriter(FILE* fout, int buffer_size=32 * 1024 * 1024);

    ~BufferWriter();

    void write(char ch);

    void write(int x);

    void write(double x);

    void write(int_t x);

    void write(const std::string& str);

    void write(const char* str);

    void space();

    void tab();

    void newline();

    void flush();

    void close();

private:
    FILE* fout;
    char* buffer;
    int chunk;
    int iptr;

}; //! class BufferWriter

} //! namespace gtb
} //! namespace gl