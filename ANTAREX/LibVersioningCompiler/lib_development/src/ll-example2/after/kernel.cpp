#include <stdio.h>
#include "kernel.hpp"

extern "C" {
void kernel()
{
  printf("Hello World\n");
  return;
}
}
