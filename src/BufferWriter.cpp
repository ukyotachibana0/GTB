/**
 * edit time : 2023-12-12
 * Buffer Writer
 */
#include "BufferWriter.hpp"
#include <sstream>
#include <iostream>

namespace gl {
namespace gtb {

BufferWriter::BufferWriter() {
    fout = nullptr;
    buffer = nullptr;
    chunk = iptr = 0;
}

BufferWriter::BufferWriter(const char* filename, int buffer_size) {
    fout = fopen(filename, "w");
    if (!fout) {
        std::cerr << "[BufferWriter Constructor] Cannot open " << filename << std::endl;
        return;
    }
    chunk = buffer_size;
    buffer = new char[buffer_size + 1];
    iptr = 0;
}

BufferWriter::BufferWriter(FILE* fout, int buffer_size) {
    this->fout = fout;
    chunk = buffer_size;
    buffer = new char[buffer_size + 1];
    iptr = 0;
}

BufferWriter::~BufferWriter() {
    close();
}

void BufferWriter::close() {
    if (fout) {
        flush();
        fclose(fout);
        fout = nullptr;
    }
    if (buffer) {
        delete [] buffer;
        buffer = nullptr;
    }
}

void BufferWriter::write(char ch) {
    buffer[iptr ++] = ch;
    if (iptr == chunk) {
        flush();
        iptr = 0;
    }
}

void BufferWriter::write(int x) {
    std::ostringstream oss;
    oss << x;
    write(oss.str());
}

void BufferWriter::write(double x) {
    std::ostringstream oss;
    oss << x;
    write(oss.str());
}

void BufferWriter::write(int_t x) {
    std::ostringstream oss;
    oss << x;
    write(oss.str());
}

void BufferWriter::write(const std::string& str) {
    for (char ch : str)
        write(ch);
}

void BufferWriter::write(const char* str) {
    const char* pch = str;
    while ((*pch) != 0) {
        write((*pch));
        pch ++;
    }
}

void BufferWriter::space() {
    write(' ');
}

void BufferWriter::tab() {
    write('\t');
}

void BufferWriter::newline() {
    write('\n');
}

void BufferWriter::flush() {
    if (!fout) {
        std::cerr << "[BufferWriter flush] fout is null" << std::endl;
        return;
    }
    fwrite(buffer, sizeof(char), iptr, fout);
    iptr = 0;
}

} //! namespace gtb
} //! namespace gl