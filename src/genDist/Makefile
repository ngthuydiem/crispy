MPI_DIR = /usr/lib/openmpi
CUDA_DIR = /usr/lib/nvidia-cuda-toolkit

ARCH = sm_20	
NVCC = $(CUDA_DIR)/bin/nvcc
CXX = $(NVCC)

NVCCOPTIONS = -arch $(ARCH) --ptxas-options=-v --compiler-options=-Wall --compiler-options=-Wno-unused-result
NVCCFLAGS = -O3 $(NVCCOPTIONS) -Xcompiler -fopenmp -I$(CUDA_DIR)/include -L$(CUDA_DIR)/lib64 -maxrregcount 56 -lineinfo -DTHRUST_HOST_SYSTEM=THRUST_HOST_SYSTEM_OMP -DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_CUDA -DSET_AFFINITY #-DDEBUG

CXXFLAGS = $(NVCCFLAGS)
NVCCLIBS = -lm -lgomp -lcudart 

CXXOBJS = genMain.o genCPU.o 
NVCCOBJS = genKernel.o genCUDA.o 
EXEC = genDist

all: $(CXXOBJS) $(NVCCOBJS)
	$(NVCC) $(NVCCFLAGS) -o $(EXEC) $(CXXOBJS) $(NVCCOBJS) $(NVCCLIBS)
	strip $(EXEC)
	#cp $(EXEC) ../../executable/
genKernel.o:
	$(NVCC) $(NVCCFLAGS) -c genKernel.cu
genCUDA.o:
	$(NVCC) $(NVCCFLAGS) -c genCUDA.cu
clean: 
	rm -rf $(EXEC)  *.o 

