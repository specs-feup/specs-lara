EXECUTABLE := syr2k.exe
CFILES := syr2k.c

#AMD
#OpenCL_SDK=/opt/AMDAPPSDK-3.0
#INCLUDE=-I${OpenCL_SDK}/include
#LIBPATH=-L${OpenCL_SDK}/lib/x86_64


#NVIDIA
OpenCL_SDK=/usr/local/cuda
INCLUDE=-I${OpenCL_SDK}/include -I"C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v8.0\include" 
#-I"C:\Users\JoaoBispo\Documents\MEGA\Work\Tasks\2017-12-04 Clava Half Experience\polybench-gpu-1.0"
LIBPATH=-L${OpenCL_SDK}/lib64 -L"C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v8.0\lib\x64"

LIB=-lOpenCL -lm

all:
	gcc -O3 ${INCLUDE} ${LIBPATH} ${CFILES} ${LIB} -o ${EXECUTABLE}

clean:
	rm -f *~ *.exe
