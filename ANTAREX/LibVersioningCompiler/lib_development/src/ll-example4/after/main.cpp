#include "main.hpp"
#include "versioningCompiler/Utils.hpp"

void libVC_error_procedure(const std::string& msg);

typedef void (*_kernel_ptr_t)();

int main(int argc, char const *argv[]) {

  vc::compiler_ptr_t my_compiler = vc::make_compiler<vc::SystemCompiler>();

  // some code here before loop body

  for (int num = 0; num < 10; num ++) {
    const std::string format = (num % 3) ? "d" : "f";
    const float number_to_print = static_cast<float>(num) / 10;

    // some code here before function call

    const std::string _kernel_source_filename = "../kernel.cpp";
    const std::string _kernel_function_name = "kernel";
    const std::string _kernel_format_quoted_option = "\\\"" + format + "\\\"";
    const vc::opt_list_t _kernel_compiler_option_list = {
      vc::make_option("-std=c++11"), // some default option
      vc::make_option("-O3"),        // some other default option
      vc::Option("def_NUM", "-DNUM=", std::to_string(number_to_print)),
      vc::Option("def_NUM_FORMAT", "-DNUM_FORMAT=", _kernel_format_quoted_option),
    };

    vc::Version::Builder _my_builder;
    _my_builder._fileName_src = _kernel_source_filename;
    _my_builder._functionName =_kernel_function_name;
    _my_builder._compiler = my_compiler;
    _my_builder._optionList = _kernel_compiler_option_list;

    vc::version_ptr_t _kernel_v1 = _my_builder.build();

    bool ok = _kernel_v1->compile();
    if (!ok) {
      // handle possible errors
      // default message behaviour: error message + exit
      libVC_error_procedure("Error compiling kernel");
      return -1;
    }
    _kernel_ptr_t kernel_v1 = (_kernel_ptr_t) _kernel_v1->getSymbol();

    kernel_v1(); // actual dynamic function call

    // some code here after function call
  }

  // some code here after loop body

  // this call has not to be versioned
  kernel();

  return 0;
}

void libVC_error_procedure(const std::string& msg)
{
  std::cerr << "libVC has incourred in an error" << std::endl;
  std::cerr << msg << std::endl;
  return;
}
