MPI_DIR = /usr/lib/openmpi

ARCH = sm_50
NVCC = nvcc
CXX = $(NVCC)

NVCCOPTIONS = -arch $(ARCH) --ptxas-options=-v --compiler-options=-Wall --compiler-options=-Wno-unused-result
NVCCFLAGS = -O3 $(NVCCOPTIONS) -Xcompiler -fopenmp -I$(MPI_DIR)/include -L$(MPI_DIR)/lib -lineinfo

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


