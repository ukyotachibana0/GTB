/**
 * edit time : 2023-12-12
 * Main procedure
 */

#include "Generation.hpp"
#include <chrono>

using namespace std;
using namespace gl::gtb;

class TimeCounter {
private:
    std::chrono::high_resolution_clock::time_point tp_start;
    std::chrono::high_resolution_clock::time_point tp_stop;

public:
    void start() {
        tp_start = std::chrono::high_resolution_clock::now();
    }

    void stop() {
        tp_stop = std::chrono::high_resolution_clock::now();
    }

    double milliseconds() {
        auto ts = std::chrono::duration_cast<std::chrono::milliseconds>(tp_stop - tp_start);
        return ts.count();
    }

    double seconds() {
        auto ts = std::chrono::duration_cast<std::chrono::seconds>(tp_stop - tp_start);
        return ts.count();
    }
};


int main(int argc, char const *argv[])
{
    srand(time(NULL));
    TimeCounter tc;
    tc.start();
    const char* filename = argv[1];
    const char* dirname = argv[2];
    string str_filename(filename);
    string str_dirname(dirname);
    Generation gen(str_filename, str_dirname);
    gen.run();
    tc.stop();
    cout << "[Main] Elapsed time : " << tc.milliseconds() << " ms." << endl;

    return 0;
}
