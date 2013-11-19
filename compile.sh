rm -r bin
mkdir bin

cd src/preprocess
make clean
make
cp preprocess ../../bin
make clean
cd ../..

cd src/kmerDist
make clean
make
cp kmerDist ../../bin
make clean
cd ../..

cd src/genDist
make clean
make
cp genDist ../../bin
make clean
cd ../..

cd src/aveclust
make clean
make
cp aveclust ../../bin
make clean
cd ../..

cd src/hcluster
make clean
make
cp hcluster ../../bin
make clean
cd ../..
