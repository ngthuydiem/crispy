COMPILER_DIR = /usr/bin/g++-4.5
MPI_DIR = /usr/lib/openmpi
CUDA_DIR = /usr/lib/nvidia-cuda-toolkit

ARCH = sm_20	#use sm_20 for Fermi architecture
NVCC = nvcc
CXX = $(NVCC)

NVCCOPTIONS = -arch $(ARCH) --ptxas-options=-v --compiler-options=-Wall --compiler-options=-Wno-unused-result
NVCCFLAGS = -O3 $(NVCCOPTIONS) -Xcompiler -fopenmp -I$(MPI_DIR)/include -I$(CUDA_DIR)/include -L$(MPI_DIR)/lib -L$(CUDA_DIR)/lib64 -lineinfo #-ccbin $(COMPILER_DIR)

CXXFLAGS = $(NVCCFLAGS)
NVCCLIBS = -lm -lgomp -lmpi -lcudart -lmpi_cxx

CXXOBJS = kmerMain.o kmerCPU.o 
NVCCOBJS =  kmerKernel.o kmerCUDA.o kmerMPI.o
EXEC = kmerDist

all: $(CXXOBJS) $(NVCCOBJS)
	$(NVCC) $(NVCCFLAGS) -o $(EXEC) $(CXXOBJS) $(NVCCOBJS) $(NVCCLIBS)
	strip $(EXEC)
kmerKernel.o:
	$(NVCC) $(NVCCFLAGS) -c kmerKernel.cu 
kmerCUDA.o: 
	$(NVCC) $(NVCCFLAGS) -c kmerCUDA.cu 
kmerMPI.o: 
	$(NVCC) $(NVCCFLAGS) -c kmerMPI.cu 
clean: 
	rm -rf $(EXEC)  *.o 


