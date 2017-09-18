#include <stdio.h>

extern "C" {
void kernel()
{
  printf("Hello World\n");
  return;
}
}
