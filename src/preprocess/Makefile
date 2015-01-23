CXX = g++ -O3

EXEC = preprocess

CXXOBJS = preprocess.o

all: $(CXXOBJS)
	$(CXX) -o $(EXEC) $(CXXOBJS)
	strip $(EXEC)
clean: 
	rm -rf $(EXEC)  *.o 


