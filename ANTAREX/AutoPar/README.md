C source files of benchmarks NAS and Polybench, before applying AutoPar Clava, which automatically inserts OpenMP pragmas.

This folder contains CMake build files that perform the auto-parallelization. Instructions for replication (Linux-only):

- Download [clava-update](http://specs.fe.up.pt/tools/clava/clava-update) in the folder you want to install Clava and run it with `sudo` (necessary for installing the Clava CMake package)
- Clone the repository [specs-lara](https://github.com/specs-feup/specs-lara)
- Go to folder `ANTAREX/AutoPar`, create a build folder (e.g. `mkdir build`), enter the build folder
- Create the makefile with CMake (e.g. `cmake ..`). 

This will call Clava and perform the auto-parallelization for all examples. It can take a while, since there are several examples, and some are quite big (e.g. thousands of lines of code). For instance, consider that you want to just parallelize a single example, `Polybench/bicg`. To accomplish this, do the following:

- Comment lines 4 and 5 in this [specs-lara/ANTAREX/AutoPar/CMakeLists.txt](https://github.com/specs-feup/specs-lara/blob/master/ANTAREX/AutoPar/CMakeLists.txt) (e.g., `#add_subdirectory(NAS)`)
- Comment all lines from 4 to 31 in [specs-lara/ANTAREX/AutoPar/Polybench/CMakeLists.txt ](https://github.com/specs-feup/specs-lara/blob/master/ANTAREX/AutoPar/Polybench/CMakeLists.txt), except for lines 4 (common library) and 9 (`bicg`).

After running `cmake ..` (or simply `make`, after the first time you run CMake, in case your default is Unix Makefiles), the parallelized version of `bicg` should be in the build folder, in a path similar to `build/Polybench/bicg/bicg_clava_weave/woven/bicg.c`
