#include "main.hpp"
#ifndef KERNEL
#include "versioningCompiler/Utils.hpp"

void libVC_error_procedure(const std::string& msg);

int main(int argc, char const *argv[]) {

  vc::compiler_ptr_t my_compiler = vc::make_compiler<vc::SystemCompiler>();

  typedef void (*kernel_t)();

  // libVersioningCompiler setup
  vc::opt_list_t kernel_opt_list_1 = {
    vc::Option("","-O3"),
    vc::Option("","-std=c++11"),
  };

  vc::Version::Builder _my_builder;
  _my_builder._compiler = my_compiler;
  _my_builder._fileName_src = "../main.cpp";
  _my_builder._functionName = "kernel";
  _my_builder._optionList = kernel_opt_list_1;

  // automatically insert a macro with the name of the function
  // this is needed to isolate the versioned function from the rest of the code
  _my_builder.addFunctionFlag();

  auto v = _my_builder.build();

  if (!v) {
    libVC_error_procedure("Parameters error during version creation");
    return -1;
  }

  bool ok = v->compile();

  if (!ok) {
    libVC_error_procedure("Compilation Error");
    return -1;
  }

  kernel_t f = (kernel_t) v->getSymbol();

  // some code here before function call

  f();

  // some code here after function call

  return 0;
}

#endif
extern "C" {
void kernel()
{
  printf("Hello world\n");
  return;
}
}
#ifndef KERNEL

void libVC_error_procedure(const std::string& msg)
{
  std::cerr << "libVC has incourred in an error" << std::endl;
  std::cerr << msg << std::endl;
  return;
}

#endif
