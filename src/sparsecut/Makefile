CXX = g++ -O3 -Wall -Wextra -Wno-unused-result -D_FILE_OFFSET_BITS=64 -std=gnu++0x

EXEC = sparsecut

CXXOBJS = sparsecut.o

all: $(CXXOBJS)
	$(CXX) -o $(EXEC) $(CXXOBJS)
	strip $(EXEC)
clean: 
	rm -rf $(EXEC)  *.o 


