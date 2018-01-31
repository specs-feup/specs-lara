#include "aos2soa.h"
#include "aos2soa-config.h"

#include <iostream>
#include <vector>
#include <chrono>
#include <cmath>

int main() {
	
	std::vector<TestStruct> aos(AOS_SIZE);
	
	// Init
	float validation = 0;
	for(int i = 0; i< AOS_SIZE; i++) {
		float x = i;
		float y = i+1;
		float z = i+2;
		
		aos.push_back({x, y, z, {}});
		validation += x + y + z;
		
		//std::cout << "x:" << x << ";y:" << y << ";z:" << z << "\n";
		//std::cout << "aos.x:" << aos[i].x << ";aos.y:" << aos[i].y << ";aos.z:" << aos[i].z << "\n";
	}
	
	
	std::chrono::high_resolution_clock::time_point clava_timing_start_0 = std::chrono::high_resolution_clock::now();
	float result = test(aos);
    std::chrono::high_resolution_clock::time_point clava_timing_end_0 = std::chrono::high_resolution_clock::now();
    auto clava_timing_duration_0 = std::chrono::duration_cast<std::chrono::milliseconds>(clava_timing_end_0 - clava_timing_start_0).count();
    std::cout << "total time: " << clava_timing_duration_0 << "ms" << std::endl;	
	
	std::cout << "result:" << result << "\n";

	if((std::abs(result-validation) / validation) < VALIDATION_ERROR) {
		std::cout << "PASSED\n";
	} else {
		std::cout << "FAILED, expected " << validation <<  "\n";
	}
	
	return 0;
}