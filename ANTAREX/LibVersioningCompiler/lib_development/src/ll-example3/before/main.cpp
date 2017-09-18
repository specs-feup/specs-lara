#include "main.hpp"

int main(int argc, char const *argv[]) {

  // some code here before function call

  // some DSL annotation here - VERSION THIS FUNCTION
  kernel();

  // some code here after function call

  return 0;
}

void kernel()
{
  printf("Hello world\n");
  return;
}
