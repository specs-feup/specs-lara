To compile, create a build folder, enter it and use cmake:

 mkdir build
 cd build
 cmake -DLVL=0 -DAFTER=OFF ..
 
Options:

 - LVL: The version of the source to compile. Currently available "levels": 0
 
 - AFTER: Choose whether to compile the original version (OFF) or the version after AoS2SoA transformation (ON)
