/***********************************************
* # Copyright 2011. Thuy Diem Nguyen
* # Contact: thuy1@e.ntu.edu.sg
* #
* # GPL 3.0 applies.
* #
* ************************************************/
#ifndef _GEN_MAIN_H_
#define _GEN_MAIN_H_

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>
#include <sys/time.h>
#include <omp.h>

using namespace std;

#define BUF_SIZE 		4096
#define MAX_NUM_PAIRS_GPU 		1024*1024*32
#define MAX_NUM_PAIRS_CPU 		1024*1024*32

#define BLOCK_SIZE		576
#define GRID_SIZE		256
#define NUM_STREAMS		1
#define NUM_PAIRS		BLOCK_SIZE * GRID_SIZE * NUM_STREAMS

#define MAX_READ_LEN		1024
#define MAX_NUM_READS		1024*1024

#define GO			-10
#define GE			-1
#define MISMATCH	-4
#define MATCH		5
#define THRESHOLD	0.15
#define BAND 		2

#define SUBMAT_DIM	4
#define EPSILON 	0.00001

//#define max(a, b) ((a)>(b)?a:b)
//#define min(a, b) ((a)<(b)?a:b)

//-----------------------------------type defines
typedef struct read
{
	int id;
	int length;
	char * sequence;

	void initialize(int readId, const char* seq);
	void release();

} READ;

void* Malloc(size_t size);

void deleteSpaces(string &line);

void deleteSpacesInSeq(string &line);

int readFile(string inFileName, READ* &readArray, int &maxLen);

void convertToBinary(READ* &readArray, int numReads, unsigned short* &digitalReadArray, int maxDigitalLen);

int loadPairs(FILE* pairFile, int * &pairArray, int &EOFTag);

void help();

void getCommandOptions(int argc, char* argv[], string &inFileName, float &threshold, int &band, bool &useGPU, int& numThreads);

void computeGenDist_CPU_full(FILE* inPairFile, string pairFileName, string distFileName, READ* &readArray, int numReads, int maxLen, float threshold, int numThreads);

void computeGenDist_CPU_band(FILE* inPairFile, string pairFileName, string distFileName, READ* &readArray, int numReads, int maxLen, float threshold, int band, int numThreads);

void computeGenDist_CUDA(FILE* inPairFile, string pairFileName, string distFileName, READ* &readArray, int numReads, int maxLen, float threshold, int band);

#endif
