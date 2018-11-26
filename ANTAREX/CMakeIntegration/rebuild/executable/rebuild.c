#include "rebuild.h"
#include "rebuild_lib.h"
#include "rebuild_lib_pre.h"

#include <math.h>
#include <stdio.h>

void rebuild_exe() {
	rebuild_lib();
	rebuild_lib_pre();
	
	printf("Pow(2,5): %f\n", pow(2.0, 5));
}

int main() {
	rebuild_exe();
}