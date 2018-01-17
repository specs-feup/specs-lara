#include "aos2soa.h"

#include <vector>
#include <iostream>
#include <chrono>


float test(std::vector<TestStruct> aos) {
	
	std::chrono::high_resolution_clock::time_point clava_timing_start_0 = std::chrono::high_resolution_clock::now();

	
    float acc = 0;
    for(int i=0; i<aos.size(); i++) {
        // Access the field directly from the array
        acc += aos[i].x;
        acc += aos[i].y;
        acc += aos[i].z;
    }

	std::chrono::high_resolution_clock::time_point clava_timing_end_0 = std::chrono::high_resolution_clock::now();
    auto clava_timing_duration_0 = std::chrono::duration_cast<std::chrono::milliseconds>(clava_timing_end_0 - clava_timing_start_0).count();
    std::cout << "kernel time: " << clava_timing_duration_0 << "ms" << std::endl;	
	
	return acc;
}