ARCH = sm_50
NVCC = nvcc
CXX = $(NVCC)

NVCCOPTIONS = -arch $(ARCH) 
NVCCFLAGS = -O3 $(NVCCOPTIONS) -Xcompiler -fopenmp 

CXXFLAGS = $(NVCCFLAGS)
NVCCLIBS = -lm -lgomp -lcudart

CXXOBJS = kmerMain.o kmerCPU.o 
NVCCOBJS =  kmerKernel.o kmerCUDA.o
EXEC = kmerDist

all: $(CXXOBJS) $(NVCCOBJS)
	$(NVCC) $(NVCCFLAGS) -o $(EXEC) $(CXXOBJS) $(NVCCOBJS) $(NVCCLIBS)
	strip $(EXEC)
kmerKernel.o:
	$(NVCC) $(NVCCFLAGS) -c kmerKernel.cu 
kmerCUDA.o: 
	$(NVCC) $(NVCCFLAGS) -c kmerCUDA.cu 
clean: 
	rm -rf $(EXEC)  *.o 


