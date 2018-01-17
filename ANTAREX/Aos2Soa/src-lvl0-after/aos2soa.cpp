#include "aos2soa.h"

#include <vector>
#include <iostream>
#include <chrono>


long test(std::vector<TestStruct> aos) {
	std::vector<long> aos_x(aos.size());
	std::vector<long> aos_y(aos.size());
	std::vector<long> aos_z(aos.size());

	for(std::vector<TestStruct>::iterator it = aos.begin(); it != aos.end(); ++it) {
		aos_x.push_back((*it).x);
		aos_y.push_back((*it).y);
		aos_z.push_back((*it).z);
		
		//std::cout << "aos.x:" << (*it).x << ";aos.y:" << (*it).y << ";aos.z:" << (*it).z << "\n";
	}
	
	std::chrono::high_resolution_clock::time_point clava_timing_start_0 = std::chrono::high_resolution_clock::now();
	
    long acc = 0;
	//for(int i=0; i<aos.size(); i++) {
    for(int i=0; i<aos_x.size(); i++) {
        // Access the field directly from the array
        //acc += aos[i].x;
		acc += aos_x[i];
        //acc += aos[i].y;
		acc += aos_y[i];
        //acc += aos[i].z;
        acc += aos_z[i];		
    }
	
	std::chrono::high_resolution_clock::time_point clava_timing_end_0 = std::chrono::high_resolution_clock::now();
    auto clava_timing_duration_0 = std::chrono::duration_cast<std::chrono::milliseconds>(clava_timing_end_0 - clava_timing_start_0).count();
    std::cout << "kernel time: " << clava_timing_duration_0 << "ms" << std::endl;		

	return acc;
}