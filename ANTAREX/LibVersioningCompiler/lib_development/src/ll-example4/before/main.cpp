#include "main.hpp"

int main(int argc, char const *argv[]) {

  // some code here before loop body

  for (int num = 0; num < 10; num ++) {
    const std::string format = (num % 3) ? "d" : "f";
    const float number_to_print = static_cast<float>(num) / 10;

    // some code here before function call

    // some DSL annotation here - VERSION THIS FUNCTION
    // add define NUM with value number_to_print
    // add define NUM_FORMAT with value format
    kernel();

    // some code here after function call
  }

  // some code here after loop body

  // this call has not to be versioned
  kernel();

  return 0;
}
