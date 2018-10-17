#include "test_lib.h"

#include <iostream>
#include "CompType.h"
 
void foo(A& a) {
}

void foo(B& a) {
}

void foo(C& a) {
}

int main() {
	std::cout << "Hello\n";
	H5::CompType aCompType = hdf5type::AType::GetCompType();
}