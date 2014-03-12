/***********************************************
* # Copyright 2011. Thuy Diem Nguyen
* # Contact: thuy1@e.ntu.edu.sg
* #
* # GPL 3.0 applies.
* #
* ************************************************/

#include "genMain.h"

#include <thrust/host_vector.h>
#include <thrust/sort.h>
#include <thrust/copy.h>

void writeVectorToFile_CPU(thrust::host_vector< thrust::pair<unsigned int, unsigned int> > h_pairVector, thrust::host_vector< float > h_distVector, string pairFileName, string distFileName, unsigned long long count, int fileId) {
	FILE * pairFile, * distFile; 
	string tempStr;	
	char buf[10];
					
	sprintf(buf, "_%d", fileId);

	tempStr = pairFileName;
	tempStr.append(buf);
	pairFile = fopen(tempStr.c_str(), "wb");
	if (pairFile == NULL){
		printf("cannot open pairFile: %s\n", tempStr.c_str());
		exit(-1);
	}	
	tempStr = distFileName;
	tempStr.append(buf);
	distFile = fopen(tempStr.c_str(), "wb");
	if (distFile == NULL){
		printf("cannot open distFile: %s\n", tempStr.c_str());
		exit(-1);
	}
			
				
	thrust::sort_by_key(h_distVector.begin(), h_distVector.end(), h_pairVector.begin());				
								
	int pairArray[BUF_SIZE*2];
	float distArray[BUF_SIZE];	

	int h = 0;
	thrust::pair<unsigned int, unsigned int> aPair;						
	
	//cout << "write to : " << tempStr << " " << count << " pairs" << endl; 
				
	for (unsigned int i = 0; i < count; ++i)
	{					
		aPair = h_pairVector[i];	
		distArray[h] = h_distVector[i];
		pairArray[h*2] = aPair.first;
		pairArray[h*2+1] = aPair.second;		
		++h;		
	/*
		if (i <= 100)
			cout << aPair.first << "\t" << aPair.second << "\t" << distArray[i] << endl;	
	*/
		if (h == BUF_SIZE) {					
			fwrite(pairArray, sizeof(unsigned int), BUF_SIZE * 2, pairFile);		
			fwrite(distArray, sizeof(float), BUF_SIZE, distFile);		
			h = 0;
		}	
	}
	
	if (h > 0) {					
		fwrite(pairArray, sizeof(unsigned int), h * 2, pairFile);		
		fwrite(distArray, sizeof(float), h, distFile);
		h = 0;
	}	
		
	fclose(pairFile);
	fclose(distFile);
}

bool insideBand(int i, int j, int kBand)
{
	return ((i - j >= -kBand) && (i - j <= kBand));
}

void computeGenDist_CPU_full(FILE* inPairFile, string pairFileName, string distFileName, READ* &readArray, int numReads, int maxLen, float threshold, int numThreads)
{
	int pairIndex=0, i, j, numPairs;	
	unsigned long long totalNumPairs = 0, count = 0;
	int fileId=0;

	int **M, **ML, **AL;
	char **ULD;	
	
	// allocate memory
	M = (int**) Malloc(sizeof(int*) * numThreads * 2);
	ML = (int**) Malloc(sizeof(int*) * numThreads * 2);
	AL = (int**) Malloc(sizeof(int*) * numThreads * 2);
	ULD = (char**) Malloc(sizeof(char*) * numThreads * 2);

	for (i = 0; i < numThreads * 2; ++i) {
		M[i] = (int*) Malloc(sizeof(int) * (maxLen+1));
		ML[i] = (int*) Malloc(sizeof(int) * (maxLen+1));
		AL[i] = (int*) Malloc(sizeof(int) * (maxLen+1));
		ULD[i] = (char*) Malloc(sizeof(char) * (maxLen+1));	
	}

	int M_up, M_left, M_diag, alignedScore;
	int length1, length2, readIndex1, readIndex2;	

	bool isMatched;
	float dist;
	int EOFTag = 0;
	int n;
	int curRow, preRow;
	int curRowFlag, preRowFlag;
	int maxMLastRow, maxMLastCol, MLRow, MLCol, ALRow, ALCol;

	// Allocate space for the pair id array	
	int *kmerPairArray;
	kmerPairArray = (int*)malloc(NUM_PAIRS * 2 * sizeof(int));
	if (kmerPairArray == NULL) {
		cout << "Error: not enough memory! EXIT..." << endl;
		fclose(inPairFile);
		exit(-1);
	}		

	size_t maxNumPairs = BLOCK_SIZE*GRID_SIZE;
	thrust::host_vector< float > h_distVector (maxNumPairs * 2);	
	thrust::host_vector< thrust::pair<unsigned int, unsigned int> > h_pairVector (maxNumPairs * 2);	

	while (EOFTag != 1)
	{
		numPairs = loadPairs(inPairFile, kmerPairArray, EOFTag);				

		#pragma omp parallel for private(n, i, j, length1, length2, readIndex1, readIndex2, isMatched, alignedScore, M_up, M_left, M_diag, curRow, preRow, curRowFlag, preRowFlag, dist, maxMLastRow, maxMLastCol, MLRow, MLCol, ALRow, ALCol) schedule(static,1)
		for (pairIndex = 0; pairIndex < numPairs; ++pairIndex)
		{	
			maxMLastRow = maxMLastCol = -1000000;		
			MLRow = MLCol = ALRow = ALCol = 0;
			n = omp_get_thread_num();
			curRowFlag = 0, preRowFlag = 1;
			curRow = curRowFlag + 2 * n;
			preRow = preRowFlag + 2 * n;
			readIndex1 = kmerPairArray[pairIndex * 2];
			readIndex2 = kmerPairArray[pairIndex * 2 + 1];

			length1 = readArray[readIndex1].length;
			length2 = readArray[readIndex2].length;

			M[curRow][0] = 0;
			ULD[curRow][0] = 0;
			ML[curRow][0] = 0;
			AL[curRow][0] = 0;

			// Initialise the first row
			for (j = 1; j <= length1; ++j)
			{
				M[curRow][j] = 0;
				ULD[curRow][j] = 0;
				ML[curRow][j] = 0;
				AL[curRow][j] = j;			
			}
			
			// Recurrence relations
			// Process row by row
			for (i = 1; i <= length2; ++i)
			{
				curRowFlag = 1 - curRowFlag;
				preRowFlag = 1 - preRowFlag;
				curRow = curRowFlag + 2 * n;
				preRow = preRowFlag + 2 * n;

				M[curRow][0] = 0;
				ULD[curRow][0] = 0;
				ML[curRow][0] = 0;
				AL[curRow][0] = i;

				for (j = 1; j <= length1; ++j)
				{
					if (readArray[readIndex2].sequence[i - 1] == readArray[readIndex1].sequence[j - 1])
					{
						alignedScore = MATCH;
						isMatched = 1;
					}
					else
					{
						alignedScore = MISMATCH;
						isMatched = 0;
					}

					M_up = M[curRow][j-1] + (ULD[curRow][j-1] & 1) * GO + (ULD[curRow][j-1] >> 2 & 1) * GE;
					M_left = M[preRow][j] + (ULD[preRow][j] & 1) * GO + (ULD[preRow][j] >> 1 & 1) * GE;
					M_diag = M[preRow][j-1] + alignedScore;
					M[curRow][j] = max(max(M_diag, M_up), M_left);

					if (M[curRow][j] == M_diag) {
						ULD[curRow][j] = 1;
						ML[curRow][j] = (ML[preRow][j-1] - isMatched) + 1;
						AL[curRow][j] = AL[preRow][j-1] + 1;						
					}
					else if (M[curRow][j] == M_up) {
						ULD[curRow][j] = 4;
						ML[curRow][j] = ML[curRow][j-1] + 1;
						AL[curRow][j] = AL[curRow][j-1] + 1;
					}
					else {
						ULD[curRow][j] = 2;
						ML[curRow][j] = ML[preRow][j]  + 1;
						AL[curRow][j] = AL[preRow][j] + 1;
					}						

					if (i == length2 && M[curRow][j] >= maxMLastCol) {
						maxMLastCol = M[curRow][j];		
						MLCol = ML[curRow][j];
						ALCol = AL[curRow][j];
					}	

					if (j == length1 && M[curRow][j] >= maxMLastRow) {
						maxMLastRow = M[curRow][j];
						MLRow = ML[curRow][j];
						ALRow = AL[curRow][j];					
					}													
				}				
			}			
			
			if (maxMLastCol > maxMLastRow)
				dist = (float)MLCol/ALCol;
			else 
				dist = (float)MLRow/ALRow;

			#pragma omp critical 
			{
#if DEBUG			
			if (maxMLastCol > maxMLastRow)
				printf("%d: (%d-%d, %d-%d) %.3f (%d/%d)\n", pairIndex, readIndex1, length1, readIndex2, length2, dist, MLCol, ALCol);
			else
				printf("%d: (%d-%d, %d-%d) %.3f (%d/%d)\n", pairIndex, readIndex1, length1, readIndex2, length2, dist, MLRow, ALRow);
#endif			
			if (dist < threshold || fabs(dist-threshold) < EPSILON)
			{
				h_pairVector[count] = thrust::make_pair(readIndex1, readIndex2);
				h_distVector[count] = dist;
				++count;
			}	
			}
		}
		
		if (count >= maxNumPairs)
		{		
			h_pairVector.resize(count);
			h_distVector.resize(count);			
			
			writeVectorToFile_CPU(h_pairVector, h_distVector, pairFileName, distFileName, count, fileId);	
			
			++ fileId;
			totalNumPairs += count;
			count = 0;										
			
			h_pairVector.resize(maxNumPairs * 2);
			h_distVector.resize(maxNumPairs * 2);	
		}		
	}
	
	if (count > 0)
	{			
		h_pairVector.resize(count);
		h_distVector.resize(count);			

		writeVectorToFile_CPU(h_pairVector, h_distVector, pairFileName, distFileName, count, fileId);	
		totalNumPairs += count;
	}		
	
	//printf("totalNumPairs: %llu\n", totalNumPairs);
	printf("%llu\n", totalNumPairs);

	// clean up
	for (i = 0; i < numReads; ++i)
		readArray[i].release();
	free(readArray);
	for (i = 0; i < numThreads * 2; ++i) {
		free(M[i]);
		free(AL[i]);
		free(ML[i]);
		free(ULD[i]);		
	}
	free(M);
	free(AL);
	free(ML);
	free(ULD);
	free(kmerPairArray);
}


void computeGenDist_CPU_band(FILE* inPairFile, string pairFileName, string distFileName, READ* &readArray, int numReads, int maxLen, float threshold, int band, int numThreads)
{
	int pairIndex, i, j, numPairs;	
	unsigned long long totalNumPairs = 0, count = 0;
	int fileId = 0;

	int maxMLastRow, maxMLastCol, MLRow, MLCol, ALRow, ALCol;
	int M_up, M_left, M_diag, alignedScore;	
	int length1, length2, readIndex1, readIndex2;	

	bool isMatched;
	float dist;
	int n, kBand;
	int curRow, preRow;
	int EOFTag = 0, curRowFlag, preRowFlag;

	int **M, **ML, **AL;
	char **ULD;	
	
	// allocate memory
	M = (int**) Malloc(sizeof(int*) * numThreads * 2);
	ML = (int**) Malloc(sizeof(int*) * numThreads * 2);
	AL = (int**) Malloc(sizeof(int*) * numThreads * 2);
	ULD = (char**) Malloc(sizeof(char*) * numThreads * 2);
	if (M==NULL || ML==NULL || AL==NULL || ULD==NULL) {
		printf("error allocating memory");
		exit(-1);
	}

	for (i = 0; i < numThreads * 2; ++i) {
		M[i] = (int*) Malloc(sizeof(int) * (maxLen+1));
		ML[i] = (int*) Malloc(sizeof(int) * (maxLen+1));
		AL[i] = (int*) Malloc(sizeof(int) * (maxLen+1));
		ULD[i] = (char*) Malloc(sizeof(char) * (maxLen+1));	
		
		if (M[i]==NULL || ML[i]==NULL || AL[i]==NULL || ULD[i]==NULL) {
			printf("error allocating memory");
			exit(-1);
		}
	}
	
	// Allocate space for the pair id array	
	int *kmerPairArray;
	kmerPairArray = (int*)malloc(NUM_PAIRS * 2 * sizeof(int));	
	if (kmerPairArray == NULL) {
		cout << "Memory error! EXIT..." << endl;
		fclose(inPairFile);
		exit(-2);
	}	
	size_t maxNumPairs = BLOCK_SIZE*GRID_SIZE;
	thrust::host_vector< float > h_distVector (maxNumPairs * 2);	
	thrust::host_vector< thrust::pair<unsigned int, unsigned int> > h_pairVector (maxNumPairs * 2);	

	while (EOFTag != 1)
	{
		numPairs = loadPairs(inPairFile, kmerPairArray, EOFTag);		

		#pragma omp parallel for private(n, i, j, kBand, length1, length2, readIndex1, readIndex2, isMatched, alignedScore, M_up, M_left, M_diag, curRow, preRow, curRowFlag, preRowFlag, dist, maxMLastRow, maxMLastCol, MLRow, MLCol, ALRow, ALCol) schedule(static,1)
		for (pairIndex = 0; pairIndex < numPairs; ++pairIndex)
		{
			maxMLastRow = maxMLastCol = -1000000;
			MLRow = MLCol = ALRow = ALCol = 0;
			n = omp_get_thread_num();

			curRowFlag = 0, preRowFlag = 1;
			curRow = curRowFlag + 2 * n;
			preRow = preRowFlag + 2 * n;
			readIndex1 = kmerPairArray[pairIndex * 2];
			readIndex2 = kmerPairArray[pairIndex * 2 + 1];

			length1 = readArray[readIndex1].length;
			length2 = readArray[readIndex2].length;
			// input is sorted by sequence length
			// hence length1 is always smaller than length2
			//kBand = min (length1, length2) / band;
			kBand =  (length1 + length2) / (band * 2);			
						
			M[curRow][0] = 0;
			ULD[curRow][0] = 0;
			ML[curRow][0] = 0;
			AL[curRow][0] = 0;

			// Initialise the first row
			for (j = 1; j <= kBand; ++j)
			{
				M[curRow][j] = 0;
				ULD[curRow][j] = 0;
				ML[curRow][j] = 0;
				AL[curRow][j] = j;
			}

			// Recurrence relations
			// Process row by row
			for (i = 1; i <= length2; ++i)
			{
				curRowFlag = 1 - curRowFlag;
				preRowFlag = 1 - preRowFlag;
				curRow = curRowFlag + 2 * n;
				preRow = preRowFlag + 2 * n;

				if (i <= kBand)
				{
					M[curRow][0] = 0;
					ULD[curRow][0] = 0;
					ML[curRow][0] = 0;
					AL[curRow][0] = i;
				}

				for (j = 1; j <= length1; ++j)
				{
					if ( insideBand(i, j, kBand) )	
					{				
						if (readArray[readIndex2].sequence[i-1] == readArray[readIndex1].sequence[j-1])
						{
							alignedScore = MATCH;
							isMatched = 1;
						}
						else
						{
							alignedScore = MISMATCH;
							isMatched = 0;
						}
					
						if (insideBand(i,j-1,kBand))
							M_up = M[curRow][j-1] + (ULD[curRow][j-1] & 1) * GO + (ULD[curRow][j-1] >> 2 & 1) * GE;							
						else
							M_up = -1000000;
	
						if (insideBand(i-1,j,kBand))				
							M_left = M[preRow][j] + (ULD[preRow][j] & 1) * GO + (ULD[preRow][j] >> 1 & 1) * GE;
						else
							M_left = -1000000;
						
						M_diag = M[preRow][j-1] + alignedScore;
						M[curRow][j] = max(max(M_diag, M_up), M_left);

						if (M[curRow][j] == M_diag) {
							ULD[curRow][j] = 1;
							ML[curRow][j] = (ML[preRow][j-1] - isMatched) + 1;
							AL[curRow][j] = AL[preRow][j-1] + 1;						
						}
						else if (M[curRow][j] == M_up) {
							ULD[curRow][j] = 4;
							ML[curRow][j] = ML[curRow][j-1] + 1;
							AL[curRow][j] = AL[curRow][j-1] + 1;
						}
						else {
							ULD[curRow][j] = 2;
							ML[curRow][j] = ML[preRow][j]  + 1;
							AL[curRow][j] = AL[preRow][j] + 1;
						}		
	
						if (i == length2 && M[curRow][j] >= maxMLastCol) 							
						{							
							maxMLastCol = M[curRow][j];			
							MLCol = ML[curRow][j];
							ALCol = AL[curRow][j];
						}						
						
						if (j == length1 && M[curRow][j] >= maxMLastRow) 							
						{
							maxMLastRow = M[curRow][j];
							MLRow = ML[curRow][j];
							ALRow = AL[curRow][j];							
						}						
					}	
					else 
						M[curRow][j] = 0;
				}				
			}

			if (maxMLastCol > maxMLastRow)
				dist = (float)MLCol/ALCol;
			else
				dist = (float)MLRow/ALRow;	
				
							
			#pragma omp critical 
			{
#if DEBUG
			if (maxMLastCol > maxMLastRow)
				printf("thread %d, pair %d: (%d, %d) %.3f (%d/%d) score=%d\n", n, pairIndex, readIndex1, readIndex2, dist, MLCol, ALCol, maxMLastCol);
			else
				printf("thread %d, pair %d: (%d, %d) %.3f (%d/%d) score=%d\n", n, pairIndex, readIndex1, readIndex2, dist, MLRow, ALRow, maxMLastRow);
#endif			
			if (dist < threshold || fabs(dist-threshold) < EPSILON)
			{
				
				h_pairVector[count] = thrust::make_pair(readIndex1, readIndex2);
				h_distVector[count] = dist;
				++count;				
			}
			}	
		}
		if (count >= maxNumPairs)
		{		
			h_pairVector.resize(count);
			h_distVector.resize(count);			
			
			writeVectorToFile_CPU(h_pairVector, h_distVector, pairFileName, distFileName, count, fileId);	
			
			++ fileId;
			totalNumPairs += count;
			count = 0;										
			
			h_pairVector.resize(maxNumPairs * 2);
			h_distVector.resize(maxNumPairs * 2);	
		}				
	}
	
	if (count > 0)
	{			
		h_pairVector.resize(count);
		h_distVector.resize(count);			

		writeVectorToFile_CPU(h_pairVector, h_distVector, pairFileName, distFileName, count, fileId);	
		totalNumPairs += count;
	}		
	
	printf("%lld\n", totalNumPairs);

	// clean up
	for (i = 0; i < numReads; ++i)
		readArray[i].release();
	free(readArray);
	for (i = 0; i < numThreads * 2; ++i) {
		free(M[i]);
		free(AL[i]);
		free(ML[i]);
		free(ULD[i]);		
	}
	free(M);
	free(AL);
	free(ML);
	free(ULD);
	free(kmerPairArray);
}


