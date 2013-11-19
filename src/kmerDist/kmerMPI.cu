/***********************************************
* # Copyright 2011. Thuy Diem Nguyen & Zejun Zheng
* # Contact: thuy1@e.ntu.edu.sg or zheng_zejun@sics.a-star.edu.sg
* #
* # GPL 3.0 applies.

* #
* ************************************************/

#include "kmerMain.h"
#include "kmerKernel.h"

// re-implement this
void computeKmerDist_MPI(READ* &readArray, FILE* pairFile, FILE * distFile, bool hasDistances, int numReads, int maxLen, float threshold, int arrayDim, int commRank, int commSize, int K) 
{	
	int devID = commRank%2;
	gpuDeviceInit(devID);
	printf("Rank: %d. Device %d.\n", commRank, devID);	 

	int i, j, stageX, stageY, row, offset, stageId, length;	
	unsigned long long totalNumPairs = 0;
	
	int maxNumTuples = maxLen - K + 1 + 1;
	int size = arrayDim * arrayDim;	
	int arraySize = size * NUM_STREAMS;
	int gridSize = (arrayDim + BLOCK_DIM - 1)/BLOCK_DIM;	
	int stageDim =  (numReads + arrayDim - 1)/arrayDim;
	
	// determine GRID_DIM and blockSize
	dim3 threadsPerBlock(BLOCK_DIM, BLOCK_DIM);	
	dim3 blocksPerGrid(gridSize, gridSize);

	// get number of SMs on this GPU
	printf("arraySize: %dx%d, stageDim: %dx%d\n", arrayDim, arrayDim, stageDim, stageDim);	
	printf("blockSize: %dx%d, gridSize: %dx%d\n", BLOCK_DIM, BLOCK_DIM, gridSize, gridSize); 	

	// declare host variables
	ushort *tupleSet;	
	// allocate host memory
	checkCudaErrors( cudaMallocHost((void**)&tupleSet, numReads * maxNumTuples * sizeof(ushort)) );
	
	for (i = 0; i < numReads; ++i)
	{
		row = i * maxNumTuples;
		length = readArray[i].length - K + 1;
		tupleSet[row + maxNumTuples - 1] = length;
		for (j = 0; j < length; ++j)
			tupleSet[row + j] = readArray[i].tuples[j];	
	}

	for (i = 0; i < numReads; ++i)
		readArray[i].finalize();
	free(readArray);
	
	// declare device variables
	float *d_distArray;
	float *h_distArray;	

	checkCudaErrors( cudaMalloc((void**)&d_distArray, arraySize * sizeof(float)) );	
	checkCudaErrors( cudaMallocHost((void**)&h_distArray, arraySize * sizeof(float)) );
	
	// use cudaArray to store tupleArraySet
	cudaChannelFormatDesc channelDesc=cudaCreateChannelDesc<ushort>();
	cudaArray *cuArray;
	size_t width, height;
	width = maxNumTuples*64;
	height = numReads/64;
	if ( (numReads&63) != 0) 
	 	++height;
	checkCudaErrors( cudaMallocArray(&cuArray, &channelDesc, width, height) );
	checkCudaErrors( cudaMemcpyToArray(cuArray, 0, 0, tupleSet, maxNumTuples * numReads * sizeof(ushort), cudaMemcpyHostToDevice) );
	checkCudaErrors( cudaBindTextureToArray(getTexRef(), cuArray, channelDesc) );
	checkCudaErrors( cudaFreeHost(tupleSet) );	
	
	cudaStream_t streams[NUM_STREAMS];

	for (i = 0; i < NUM_STREAMS; ++i) 
		checkCudaErrors( cudaStreamCreate(&streams[i]) );			
	
	int stageSize = stageDim * (stageDim + 1) / 2;			
	
	for (j = 0; j < stageSize; j += NUM_STREAMS)
	{		
		for (i = 0; i < NUM_STREAMS; ++i) {
			offset = i * size;		
			stageId = i + j;
			
			if (stageId < stageSize) {
				Trag_reverse_eq(stageId, stageDim, stageX, stageY);													
						        								
				launchKmerKernel(streams[i], blocksPerGrid, threadsPerBlock, d_distArray+offset, numReads, maxNumTuples, stageX, stageY, arrayDim);	
				
				checkCudaErrors( cudaMemcpyAsync(h_distArray+offset, d_distArray+offset, size * sizeof(float), cudaMemcpyDeviceToHost, streams[i]) );							
				 		
			}
		}	
					
		for (i = 0; i < NUM_STREAMS; ++i) {								
			offset = i * size;	
			stageId = i + j;				

			if (stageId < stageSize) {			
											
				checkCudaErrors( cudaStreamSynchronize(streams[i]) );	
								
				Trag_reverse_eq(stageId, stageDim, stageX, stageY);												
				
				writeToFile(pairFile, distFile, hasDistances, h_distArray+offset, stageX, stageY, arrayDim, threshold, totalNumPairs);						
			}
		}	
	}

	for (i = 0; i < NUM_STREAMS; ++i) 
		checkCudaErrors( cudaStreamDestroy(streams[i]) );
				
	// clean up host variables	
	checkCudaErrors( cudaFreeHost(h_distArray) );
	checkCudaErrors( cudaFree(d_distArray) );		
	
	// clean up device variables
	checkCudaErrors( cudaUnbindTexture(getTexRef()) );
	checkCudaErrors( cudaFreeArray(cuArray) );

	printf("totalNumPairs: %llu\n", totalNumPairs);	
}



