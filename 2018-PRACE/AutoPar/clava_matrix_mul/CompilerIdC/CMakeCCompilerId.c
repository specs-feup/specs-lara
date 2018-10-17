/*cv-qualifiers did not exist in K&R C*/
/*Version number components: V=Version, R=Revision, P=Patch
Version date components:   YYYY=Year, MM=Month,   DD=Day*/
/*__INTEL_COMPILER = VRP*/
/*__INTEL_COMPILER_BUILD_DATE = YYYYMMDD*/
/*_MSC_VER = VVRR*/
/*__BORLANDC__ = 0xVRR*/
/*__WATCOMC__ = VVRR*/
/*__WATCOMC__ = VVRP + 1100*/
/*__SUNPRO_C = 0xVRRP*/
/*__SUNPRO_CC = 0xVRP*/
/*__HP_cc = VVRRPP*/
/*__DECC_VER = VVRRTPPPP*/
/*__IBMC__ = VRP*/
/*__IBMC__ = VRP*/
/*__IBMC__ = VRP*/
/*__TI_COMPILER_VERSION__ = VVVRRRPPP*/
/*_MSC_VER = VVRR*/
/*_MSC_VER = VVRR*/
/*_MSC_VER = VVRR*/
/*_MSC_FULL_VER = VVRRPPPPP*/
/*_MSC_FULL_VER = VVRRPPPP*/
/*__VISUALDSPVERSION__ = 0xVVRRPP00*/
/*__ARMCC_VERSION = VRRPPPP*/
/*__ARMCC_VERSION = VRPPPP*/
/*SDCC = VRP*/
/*_SGI_COMPILER_VERSION = VRP*/
/*_COMPILER_VERSION = VRP*/
/*These compilers are either not known or too old to define an
identification macro.  Try to identify the platform and guess that
it is the native compiler.*/
/*unknown compiler*/
/*Construct the string literal in pieces to prevent the source from
getting matched.  Store it in a pointer rather than an array
because some compilers will just produce instructions to fill the
array rather than assigning a pointer to a static array.*/
char const * info_compiler = "INFO" ":" "compiler[" COMPILER_ID "]";
/*Identify known platforms by name.*/
/*unknown platform*/
/*unknown platform*/
/*For windows compilers MSVC and Intel we can determine
the architecture of the compiler being used.  This is because
the compilers do not have flags that can change the architecture,
but rather depend on which compiler is being used
*/
/*unknown architecture*/
/*unknown architecture*/
/*unknown architecture*/
/*Convert integer to decimal digit literals.*/
/*Convert integer to hex digit literals.*/
/*Construct a string literal encoding the version number components.*/
char const info_version[50] = {'I', 'N', 'F', 'O', ':', 'c', 'o', 'm', 'p', 'i', 'l', 'e', 'r', '_', 'v', 'e', 'r', 's', 'i', 'o', 'n', '[', ('0' + (((3) / 10000000) % 10)), ('0' + (((3) / 1000000) % 10)), ('0' + (((3) / 100000) % 10)), ('0' + (((3) / 10000) % 10)), ('0' + (((3) / 1000) % 10)), ('0' + (((3) / 100) % 10)), ('0' + (((3) / 10) % 10)), ('0' + ((3) % 10)), '.', ('0' + (((8) / 10000000) % 10)), ('0' + (((8) / 1000000) % 10)), ('0' + (((8) / 100000) % 10)), ('0' + (((8) / 10000) % 10)), ('0' + (((8) / 1000) % 10)), ('0' + (((8) / 100) % 10)), ('0' + (((8) / 10) % 10)), ('0' + ((8) % 10)), '.', ('0' + (((0) / 10000000) % 10)), ('0' + (((0) / 1000000) % 10)), ('0' + (((0) / 100000) % 10)), ('0' + (((0) / 10000) % 10)), ('0' + (((0) / 1000) % 10)), ('0' + (((0) / 100) % 10)), ('0' + (((0) / 10) % 10)), ('0' + ((0) % 10)), ']', '\0'};
/*Construct a string literal encoding the internal version number.*/
/*Construct a string literal encoding the version number components.*/
/*Construct the string literal in pieces to prevent the source from
getting matched.  Store it in a pointer rather than an array
because some compilers will just produce instructions to fill the
array rather than assigning a pointer to a static array.*/
char const * info_platform = "INFO" ":" "platform[" PLATFORM_ID "]";
char const * info_arch = "INFO" ":" "arch[" ARCHITECTURE_ID "]";
char const * info_language_dialect_default = "INFO" ":" "dialect_default[" C_DIALECT "]";
/*--------------------------------------------------------------------------*/

int main(int argc, char * argv[]) {
   int require = 0;
   require += info_compiler[argc];
   require += info_platform[argc];
   require += info_arch[argc];
   require += info_version[argc];
   require += info_language_dialect_default[argc];
   (void) argv;
   
   return require;
}
