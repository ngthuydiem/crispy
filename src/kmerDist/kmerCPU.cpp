/***********************************************
* # Copyright 2011. Thuy Diem Nguyen
* # Contact: thuy1@e.ntu.edu.sg
* #
* # GPL 3.0 applies.
* #
* ************************************************/

#include "kmerMain.h"

void computeKmerDist_CPU(READ* &readArray, FILE* pairFile, FILE* distFile, bool hasDistances, int numReads, float threshold, int arrayDim, int K)
{
	int i, j, k, l;
	float dist;
	unsigned long long totalNumPairs = 0;
	unsigned short x, y, length1, length2, matches;

	int pairArray[BUF_SIZE*2];
	float distArray[BUF_SIZE];

	int h = 0;
	for (i = 0; i < numReads; ++i)
	{		
		#pragma omp parallel 
		{	
		#pragma omp for nowait  private(j, k, l, length1, length2, matches, x, y, dist) schedule(static,1)
		for (j = i + 1; j < numReads; ++j)
		{
			matches = 0;
			length1 = readArray[i].length - K + 1;
			length2 = readArray[j].length - K + 1;

			for (k = 0, l = 0; (k < length1) && (l < length2);)
			{
				x = readArray[i].tuples[k];
				y = readArray[j].tuples[l];

				matches = matches + (unsigned short)(x==y);
				k = k + (unsigned short)(x<=y);
				l = l + (unsigned short)(x>=y);	
			}

			dist = 1.0f - (float) matches / min(length1, length2);
			#pragma omp critical
			if (dist < threshold || fabs(dist-threshold) < EPSILON)
			{
				pairArray[h*2] = i;
				pairArray[h*2+1] = j;
				if (hasDistances)
					distArray[h] = dist;
				++h;
				if (h == BUF_SIZE) {
					fwrite(pairArray, sizeof(int),BUF_SIZE*2, pairFile);
					if (hasDistances)
						fwrite(distArray, sizeof(float),BUF_SIZE, distFile);
					totalNumPairs += BUF_SIZE;
					h = 0;
				}
			}
		}
		}
	}

	if (h > 0) {
		fwrite(pairArray, sizeof(int),h*2, pairFile);
		if (hasDistances)
			fwrite(distArray, sizeof(float),h, distFile);
		totalNumPairs += h;					
	}
	
	for (i = 0; i < numReads; ++i)
		readArray[i].finalize();
	free(readArray);	

	printf("totalNumPairs: %llu\n", totalNumPairs);	
}

