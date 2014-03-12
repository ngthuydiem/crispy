/***********************************************
* # Copyright 2011. Thuy Diem Nguyen
* # Contact: thuy1@e.ntu.edu.sg
* #
* # GPL 3.0 applies.
* #
* ************************************************/

#ifndef _KMER_MAIN_H_
#define _KMER_MAIN_H_

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <time.h>
#include <math.h>
#include <omp.h>
#include <string.h>
#include <mpi.h>
#include <sys/time.h>

using namespace std;

#define BUF_SIZE 			4096
#define BLOCK_DIM 			16
#define THRESHOLD 			0.5
#define EPSILON 			0.00001
#define DEFAULT_NUM_READS 	250000
#define MAX_READ_LEN		1024*2
#define MAX_NUM_READS		1024*1024

#define NUM_STREAMS			8

#define max(a, b) ((a)>(b)?a:b)
#define min(a, b) ((a)<(b)?a:b)

// define READ

typedef struct read {
	int id;
	int length;
	char* sequence;
	int numTuples;
	unsigned short* tuples;

	void initialize(int readId, const char* seq, int K);
	void finalize();
	void formTuples(int K);
	void sortTuples();

} READ;

// define function prototypes

void writeToFile(FILE *pairFile, FILE *distFile, bool hasDistances, float *h_distArray, int stageX, int stageY, int arrayDim, float threshold, unsigned long long & totalNumPairs);

void Trag_reverse_eq(int index, int N, int& row, int& col);

int compareTwoTuples(const void* t1, const void* t2);

void removeNewLine(string &line);

void removeInvalidChar(string &line);

int readFile(string inFileName, READ* &readArray, int &maxLen, int K);

void help();

void getCommandOptions(int argc, char* argv[], string &inFileName, float &threshold, bool &hasDistances, int &num_threads, bool &useGPU, bool &useMPI, int &K);

void computeKmerDist_CPU(READ* &readArray, FILE* pairFile, FILE* distFile, bool hasDistances, int numReads, float threshold, int arrayDim, int K);

void computeKmerDist_CUDA(READ* &readArray, FILE* pairFile, FILE* distFile, bool hasDistances, int numReads, int maxLen, float threshold, int arrayDim, int K);

void computeKmerDist_MPI(READ* &readArray, FILE* pairFile, FILE * distFile, bool hasDistances, int numReads, int maxLen, float threshold, int arrayDim, int commRank, int commSize, int K);

#endif
