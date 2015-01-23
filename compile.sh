SM=$1
if [ -z "$SM" ]; then
	SM=50
fi

# remove existing executables
rm -r bin
mkdir bin

# compile CUDA codes
cd src/kmerDist
make clean
make ARCH=sm_$SM
cp kmerDist ../../bin
make clean
cd ../..

cd src/genDist
make clean
make ARCH=sm_$SM
cp genDist ../../bin
make clean
cd ../..

# compile C++ codes
cd src/preprocess
make clean
make
cp preprocess ../../bin
make clean
cd ../..

cd src/aveclust
make clean
make
cp aveclust ../../bin
make clean
cd ../..

cd src/sparsecut
make clean
make
cp sparsecut ../../bin
make clean
cd ../..
