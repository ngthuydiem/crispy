/***********************************************
* # Copyright 2011. Thuy Diem Nguyen
* # Contact: thuy1@e.ntu.edu.sg
* #
* # GPL 3.0 applies.
* #
* ************************************************/

#ifndef _GEN_KERNEL_H_
#define _GEN_KERNEL_H_

#include "genMain.h"
////////////////////////////////////////////////////////////////////////////////
// These are CUDA Helper functions

// This will output the proper CUDA error strings in the event that a CUDA host call returns an error
#define checkCudaErrors(err)           __checkCudaErrors (err, __FILE__, __LINE__)
inline void __checkCudaErrors( cudaError err, const char *file, const int line )
{
    if( cudaSuccess != err) {
        fprintf(stderr, "%s(%i) : CUDA Runtime API error %d: %s.\n",
                file, line, (int)err, cudaGetErrorString( err ) );
        exit(-1);
    }
}

// This will output the proper error string when calling cudaGetLastError
#define getLastCudaError(msg)      __getLastCudaError (msg, __FILE__, __LINE__)
inline void __getLastCudaError( const char *errorMessage, const char *file, const int line )
{
    cudaError_t err = cudaGetLastError();
    if( cudaSuccess != err) {
        fprintf(stderr, "%s(%i) : getLastCudaError() CUDA error : %s : (%d) %s.\n",
                file, line, errorMessage, (int)err, cudaGetErrorString( err ) );
        exit(-1);
    }
}

inline __device__ bool insideBand(int i, int j, int kBand)
{
	return ((i - j >= -kBand) && (i - j <= kBand));
}

texture<ushort, 2, cudaReadModeElementType> &getSeqTexRef(void);

void launchGenKernel_full(cudaStream_t stream, size_t blocksPerGrid, size_t threadsPerBlock, int *d_pairArray, float *d_distArray, int maxDigitalLen, int numPairs);

void launchGenKernel_band(cudaStream_t stream, size_t blocksPerGrid, size_t threadsPerBlock, int *d_pairArray, float *d_distArray, int maxDigitalLen, int numPairs, int band);

int gpuDeviceInit(int devID);

#endif
