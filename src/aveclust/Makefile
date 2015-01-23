CXX = g++ -O3 -D_FILE_OFFSET_BITS=64 -std=gnu++0x

EXEC = aveclust

CXXOBJS = aveclust.o

all: $(CXXOBJS)
	$(CXX) -o $(EXEC) $(CXXOBJS)
	strip $(EXEC)
clean: 
	rm -rf $(EXEC)  *.o 


