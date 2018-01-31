#ifndef AOS2SOA_H
#define AOS2SOA_H

#include "aos2soa-config.h"

#include <vector>

struct TestStruct {
    float x;
    float y;
    float z;
    int padding[PADDING];
};

float test(std::vector<TestStruct> aos);

#endif