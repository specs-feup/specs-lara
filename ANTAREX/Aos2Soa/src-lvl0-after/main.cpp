#include "aos2soa.h"

#include <iostream>
#include <vector>
#include <chrono>


#define N 10000000

int main() {
	
	std::vector<TestStruct> aos;
	
	// Init
	long validation = 0;
	for(int i = 0; i< N; i++) {
		long x = i;
		long y = i+1;
		long z = i+2;
		
		aos.push_back({x, y, z, {}});
		validation += x + y + z;
		
		//std::cout << "x:" << x << ";y:" << y << ";z:" << z << "\n";
		//std::cout << "aos.x:" << aos[i].x << ";aos.y:" << aos[i].y << ";aos.z:" << aos[i].z << "\n";
	}
	
	
	std::chrono::high_resolution_clock::time_point clava_timing_start_0 = std::chrono::high_resolution_clock::now();
	long result = test(aos);
    std::chrono::high_resolution_clock::time_point clava_timing_end_0 = std::chrono::high_resolution_clock::now();
    auto clava_timing_duration_0 = std::chrono::duration_cast<std::chrono::milliseconds>(clava_timing_end_0 - clava_timing_start_0).count();
    std::cout << "total time: " << clava_timing_duration_0 << "ms" << std::endl;	
	
	std::cout << "result:" << result << "\n";

	if(result == validation) {
		std::cout << "PASSED\n";
	} else {
		std::cout << "FAILED, expected " << validation <<  "\n";
	}
	
	return 0;
}