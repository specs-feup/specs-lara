#include "aos2soa-config.h"

#include <vector>

struct TestStruct {
    long x;
    long y;
    long z;
    int padding[PADDING];
};

long test(std::vector<TestStruct> aos);