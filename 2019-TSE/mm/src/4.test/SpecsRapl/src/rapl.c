#include "rapl.h"
#include "rapl_fs.h"

#include <stdio.h>

double rapl_monitor_report() {
	printf("[SPECS_RAPL] DEPRECATED: rapl_monitor_report() has been deprecated, use rapl_energy() instead\nIf this code was generated by Clava, please update Clava\n");
	return -1.0;
}

int rapl_monitor_start() {
	printf("[SPECS_RAPL] DEPRECATED: rapl_monitor_start() has been deprecated and is no longer needed, use rapl_energy() instead\nIf this code was generated by Clava, please update Clava\n");	
	return -1;
}


long long rapl_energy() {
	return rapl_fs_energy();
}


int hasLooped() {
	return hasLooped_fs();
}