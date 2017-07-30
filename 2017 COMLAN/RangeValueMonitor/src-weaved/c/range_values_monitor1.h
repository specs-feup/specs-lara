#ifndef RANGE_VALUES_MONITOR1_H_
#define RANGE_VALUES_MONITOR1_H_

extern double monitor1_range_min[2];
extern double monitor1_range_max[2];
void monitor1_range_update(unsigned int id, double value);
void monitor1_range_init();
void monitor1_print_ranges();

#endif
