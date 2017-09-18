#include <stdio.h>
#include <string>

#ifndef NUM
#define NUM 3.14
#endif

#ifndef NUM_FORMAT
#define NUM_FORMAT "d"
#endif

void kernel()
{
  float number = NUM;
  const std::string format_number = NUM_FORMAT;
  const std::string format_string = "Hello World. I print %" + format_number + "\n";

  printf(format_string.c_str(), number);
  return;
}
