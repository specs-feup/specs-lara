#include "main.hpp"
#include "versioningCompiler/Utils.hpp"

void libVC_error_procedure(const std::string& msg);

int main(int argc, char const *argv[]) {

  typedef void (*kernel_t)();

  vc::compiler_ptr_t my_compiler = vc::make_compiler<vc::SystemCompiler>();

  // libVersioningCompiler setup
  vc::opt_list_t kernel_opt_list_1 = {
    vc::Option("","-O3"),
    vc::Option("","-std=c++11"),
  };

  vc::Version::Builder _my_builder;
  _my_builder._compiler = my_compiler;
  _my_builder._fileName_src = "../kernel.cpp";
  _my_builder._functionName = "kernel";
  _my_builder._optionList = kernel_opt_list_1;

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

  kernel();

  return 0;
}

void libVC_error_procedure(const std::string& msg)
{
  std::cerr << "libVC has incourred in an error" << std::endl;
  std::cerr << msg << std::endl;
  return;
}
