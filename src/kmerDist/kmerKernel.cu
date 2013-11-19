/***********************************************
* # Copyright 2011. Thuy Diem Nguyen & Zejun Zheng
* # Contact: thuy1@e.ntu.edu.sg or zheng_zejun@sics.a-star.edu.sg
* #
* # GPL 3.0 applies.
* #
* ************************************************/

// Note: don't use_fast_math option
#include "kmerKernel.h"

texture<ushort, 2, cudaReadModeElementType> texRef;

texture<ushort, 2, cudaReadModeElementType> &getTexRef(void)
{
        return texRef;
}

__global__ void kmerKernel(float *distArray, int numReads, int maxNumTuples, int stageX, int stageY, int arrayDim) 
{	
	int i = blockIdx.y * blockDim.y + threadIdx.y;	
	int j = blockIdx.x * blockDim.x + threadIdx.x;	

	int row = stageX * arrayDim + i;
	int col = stageY * arrayDim + j;
	int index = i * arrayDim + j;

	ushort matches = 0;
	ushort x, y, k, l;
	ushort tuple1Length, tuple2Length;

	if ( (row < col) && (col < numReads) )
	{			
		tuple1Length = tex2D(texRef, maxNumTuples * (row & 15) + (maxNumTuples-1), row/16); // tex2D( texRef, width, height )
		tuple2Length = tex2D(texRef, maxNumTuples * (col & 15) + (maxNumTuples-1), col/16); // tex2D( texRef, width, height )

		for (k = 0, l = 0; (k < tuple1Length) && (l < tuple2Length);)
		{
			x = tex2D(texRef, maxNumTuples * (row & 15) + k, row/16);
			y = tex2D(texRef, maxNumTuples * (col & 15) + l, col/16);

			matches = matches + (ushort)(x==y);
			k = k + (ushort)(x<=y);
			l = l + (ushort)(x>=y);	
		}

		distArray[index] = 1.0f - (float)matches / min(tuple1Length,tuple2Length);		
	}
	else
		distArray[index] = 1.1f;
}

void launchKmerKernel(cudaStream_t stream, dim3 blocksPerGrid, dim3 threadsPerBlock, float* d_distArray, int numReads, int maxNumTuples, int stageX, int stageY, int arrayDim) 
{	

	kmerKernel<<<blocksPerGrid, threadsPerBlock, 0, stream>>>(d_distArray, numReads, maxNumTuples, stageX, stageY, arrayDim);			
	
}


// General GPU Device CUDA Initialization
int gpuDeviceInit(int devID)
{
    int deviceCount;
    checkCudaErrors(cudaGetDeviceCount(&deviceCount));
    if (deviceCount == 0) {
        fprintf(stderr, "gpuDeviceInit() CUDA error: no devices supporting CUDA.\n");
        exit(-1);
    }
    if (devID < 0) 
        devID = 0;
    if (devID > deviceCount-1) {
        fprintf(stderr, "\n");
        fprintf(stderr, ">> %d CUDA capable GPU device(s) detected. <<\n", deviceCount);
        fprintf(stderr, ">> gpuDeviceInit (-device=%d) is not a valid GPU device. <<\n", devID);
        fprintf(stderr, "\n");
        return -devID;
    }
    
    cudaDeviceProp deviceProp;
    checkCudaErrors( cudaGetDeviceProperties(&deviceProp, devID) );
    if (deviceProp.major < 1) {
        fprintf(stderr, "gpuDeviceInit(): GPU device does not support CUDA.\n");
        exit(-1);                                                  \
    }
    
    checkCudaErrors( cudaSetDevice(devID) );
    printf("> gpuDeviceInit() CUDA device [%d]: %s\n", devID, deviceProp.name);
    return devID;
}
