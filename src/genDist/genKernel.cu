/***********************************************
* # Copyright 2011. Thuy Diem Nguyen
* # Contact: thuy1@e.ntu.edu.sg
* #
* # GPL 3.0 applies.
* #
* ************************************************/

#include "genKernel.h"

texture<ushort, 2, cudaReadModeElementType> seqTexRef;

texture<ushort, 2, cudaReadModeElementType> &getSeqTexRef(void) {
	return seqTexRef;
}

__device__ void computeSubMatrix_full(ushort read1, ushort read2, 
 char *ULD, int *M,  uint* L, char4 &subULD, int4 &subM, uint4 &subL, 
 int numRows, int numCols, int i, int j, int endJ, int &maxMLastCol, uint &LCol)
{
	bool isMatched;
	ushort symbol1, symbol2;
	int k, rowIndex;
	int tempM, tempML, tempAL, M_up, M_left, M_diag, alignedScore;
	char subULD0;
	int subM0, diagonalM;
	uint subL0, diagonalL;		
		
	// Initialize only for the first row of the first rowBlock
	if (i == 0)
	{			
		subULD0 = 0;
		subM0 = 0;		
		subULD = make_char4(0,0,0,0);
		subM = make_int4(0,0,0,0);

		subL0 = j*SUBMAT_DIM;
		subL.x = subL0+1;
		subL.y = subL.x+1;
		subL.z = subL.y+1;
		subL.w = subL.z+1;	
	}

	for (k = 1; k <= numRows; ++k)
	{
		diagonalM = subM0;
		diagonalL = subL0;
		rowIndex = i * SUBMAT_DIM + k;

		if (j == 0)
		{
			subULD0 = 0;
			subM0 = 0;
			subL0 = rowIndex;
		}
		else
		{
			// copy from the whole matrix to the first column of the submatrix
			subULD0 = ULD[rowIndex];
			subM0 = M[rowIndex];
			subL0 = L[rowIndex];
		}

		// extract the symbol from the ushort
		// shift the 2 bits to the left end
		symbol1 = read1 << ((k - 1) * 2);
		// shift the 3 bits all the way to the right end
		symbol1 = (symbol1 >> 14);

		symbol2 = read2 >> 14;

		if (symbol1 == symbol2)
		{
			alignedScore = MATCH;
			isMatched = 1;
		}
		else
		{
			alignedScore = MISMATCH;
			isMatched = 0;
		}

		// vertical
		M_up = subM.x + (subULD.x & 1) * GO + (subULD.x >> 2 & 1) * GE;
		// horizontal
		M_left = subM0 + (subULD0 & 1) * GO + (subULD0 >> 1 & 1) * GE;
		// diagonal
		M_diag = diagonalM+alignedScore;
		
		tempM = max(max(M_diag, M_up), M_left);

		if (tempM == M_diag) {
			subULD.x = 1; // 001
			tempAL = ((diagonalL << 16) >> 16) + 1;
			tempML = ((diagonalL >> 16) - isMatched) + 1;
		}
		else if (tempM == M_up) {
			subULD.x = 4; // 100
			tempAL = ((subL.x << 16) >> 16) + 1;
			tempML = (subL.x >> 16) + 1;  
		}
		else {
			subULD.x = 2; // 010
			tempAL = ((subL0 << 16) >> 16) + 1;
			tempML = (subL0 >> 16) + 1;
		}
		
		diagonalM = subM.x;
		diagonalL = subL.x;

		subM.x = tempM;
		subL.x = ((tempML<<16) | (tempAL & 0xffff));

		if (j == endJ-1 && numCols == 1 && subM.x >= maxMLastCol) {
			maxMLastCol = subM.x;
			LCol = subL.x;			
		}	
			
		if (numCols > 1) {
			symbol2 = read2 << 2;
			symbol2 = (symbol2 >> 14);

			if (symbol1 == symbol2)
			{
				alignedScore = MATCH;
				isMatched = 1;
			}
			else
			{
				alignedScore = MISMATCH;
				isMatched = 0;
			}

			M_up = subM.y + (subULD.y & 1) * GO + (subULD.y >> 2 & 1) * GE;
			M_left = subM.x + (subULD.x & 1) * GO + (subULD.x >> 1 & 1) * GE;
			// diagonal
			M_diag = diagonalM+alignedScore;
		
			tempM = max(max(M_diag, M_up),M_left);

			if (tempM == M_diag) {
				subULD.y = 1; // 001
				tempAL = ((diagonalL << 16) >> 16) + 1;
				tempML = ((diagonalL >> 16) - isMatched) + 1;
			}
			else if (tempM == M_up) {
				subULD.y = 4; // 100
				tempAL = ((subL.y << 16) >> 16) + 1;
				tempML = (subL.y >> 16) + 1;  
			}
			else {
				subULD.y = 2; // 010
				tempAL = ((subL.x << 16) >> 16) + 1;
				tempML = (subL.x >> 16) + 1;
			}
		
		
			diagonalM = subM.y;
			diagonalL = subL.y;

			subM.y = tempM;
			subL.y = ((tempML<<16) | (tempAL & 0xffff));

			if (j == endJ-1 && numCols == 2 && subM.y >= maxMLastCol) {
				maxMLastCol = subM.y;
				LCol = subL.y;			
			}		
		} 
 
		if (numCols > 2) {
			symbol2 = read2 << 4;
			symbol2 = (symbol2 >> 14);

			if (symbol1 == symbol2)
			{
				alignedScore = MATCH;
				isMatched = 1;
			}
			else
			{
				alignedScore = MISMATCH;
				isMatched = 0;
			}

			M_up = subM.z + (subULD.z & 1) * GO + (subULD.z >> 2 & 1) * GE;
			M_left = subM.y + (subULD.y & 1) * GO + (subULD.y >> 1 & 1) * GE;
			// diagonal
			M_diag = diagonalM+alignedScore;
		
			tempM = max(max(M_diag, M_up),M_left);

			if (tempM == M_diag) {
				subULD.z = 1; // 001
				tempAL = ((diagonalL << 16) >> 16) + 1;
				tempML = ((diagonalL >> 16) - isMatched) + 1;
			}
			else if (tempM == M_up) {
				subULD.z = 4; // 100
				tempAL = ((subL.z << 16) >> 16) + 1;
				tempML = (subL.z >> 16) + 1;  
			}
			else {
				subULD.z = 2; // 010
				tempAL = ((subL.y << 16) >> 16) + 1;
				tempML = (subL.y >> 16) + 1;
			}
				
			diagonalM = subM.z;
			diagonalL = subL.z;

			subM.z = tempM;
			subL.z = ((tempML<<16) | (tempAL & 0xffff));
			if (j == endJ-1 && numCols == 3 && subM.z >= maxMLastCol) {
				maxMLastCol = subM.z;
				LCol = subL.z;				
			}		
		} 

		if (numCols > 3) {
			symbol2 = read2 << 6;
			symbol2 = (symbol2 >> 14);

			if (symbol1 == symbol2)
			{
				alignedScore = MATCH;
				isMatched = 1;
			}
			else
			{
				alignedScore = MISMATCH;
				isMatched = 0;
			}

			M_up = subM.w + (subULD.w & 1) * GO + (subULD.w >> 2 & 1) * GE;
			M_left = subM.z + (subULD.z & 1) * GO + (subULD.z >> 1 & 1) * GE;
			// diagonal
			M_diag = diagonalM+alignedScore;
		
			tempM = max(max(M_diag, M_up), M_left);

			if (tempM == M_diag) {
				subULD.w = 1; // 001
				tempAL = ((diagonalL << 16) >> 16) + 1;
				tempML = ((diagonalL >> 16) - isMatched) + 1;
			}
			else if (tempM == M_up) {
				subULD.w = 4; // 100
				tempAL = ((subL.w << 16) >> 16) + 1;
				tempML = (subL.w >> 16) + 1;  
			}
			else {
				subULD.w = 2; // 010
				tempAL = ((subL.z << 16) >> 16) + 1;
				tempML = (subL.z >> 16) + 1;
			}
		
			subM.w = tempM;
			subL.w = ((tempML << 16) | (tempAL & 0xffff));			

			if (j == endJ - 1 && subM.w >= maxMLastCol) {
				maxMLastCol = subM.w;
				LCol = subL.w;				
			}
			
			if (j < endJ - 1) // is not last col block
			{
				// copy the last row of the submatrix to the whole matrix
				ULD[rowIndex] = subULD.w;
				M[rowIndex] = subM.w;
				L[rowIndex] = subL.w;
			}				
		}			
	}
}

__global__ void genKernel_full(int* pairArray, float *distArray, int maxDigitalLen, int numPairs)
{
	int pairIndex = blockIdx.x * blockDim.x + threadIdx.x;

	char ULD[MAX_READ_LEN + 1]; // up left diagonal matrix
	int M[MAX_READ_LEN + 1];
	uint L[MAX_READ_LEN + 1]; // combination of ML (first 16 bits) and AL (last 16 bits))
	
	char4 subULD;
	int4 subM;
	uint4 subL;

	int maxMLastRow=-1000000, maxMLastCol=-1000000;
	uint LRow, LCol;
	
	ushort read1, read2, length1, length2;
	int i, j, numRows, numCols, readIndex1, readIndex2, endI, endJ;

	// For each ID pair
	if (pairIndex < numPairs)
	{
		readIndex1 = pairArray[pairIndex*2];
		readIndex2 = pairArray[pairIndex*2+1];

		length1 = tex2D(seqTexRef,  maxDigitalLen * (readIndex1 & 15) + (maxDigitalLen-1), readIndex1/16);
		length2 = tex2D(seqTexRef, maxDigitalLen * (readIndex2 & 15) + (maxDigitalLen-1), readIndex2/16);

		// Recurrence relations
		// Process row by row
		endJ = length2/SUBMAT_DIM + (int)(length2 % SUBMAT_DIM != 0);
		endI = length1/SUBMAT_DIM + (int)(length1 % SUBMAT_DIM != 0);

		for (j = 0; j < endJ; ++j)
		{		
			read2 = tex2D(seqTexRef, maxDigitalLen * (readIndex2 & 15) + j, readIndex2/16);	
			if (length2 % SUBMAT_DIM != 0 && j == endJ - 1) // last col block
				numCols = length2 % SUBMAT_DIM;
			else
				numCols = SUBMAT_DIM;		

			
			for (i = 0; i < endI; ++i)
			{							
				read1 = tex2D(seqTexRef, maxDigitalLen * (readIndex1 & 15) + i, readIndex1/16);
				if (length1 % SUBMAT_DIM != 0 && i == endI - 1) // last col block {
					numRows = length1 % SUBMAT_DIM;
				else
					numRows = SUBMAT_DIM;	
						
				// compute submatrices
				computeSubMatrix_full(read1, read2, ULD, M, L, subULD, subM, subL, numRows, numCols, i, j, endJ, maxMLastCol, LCol );					
				if (i == endI-1) {					
					if (subM.x >= maxMLastRow) {
						maxMLastRow = subM.x;
						LRow = subL.x;
					} 
					if (subM.y >= maxMLastRow) {
						maxMLastRow = subM.y;
						LRow = subL.y;
					} 
					if (subM.z >= maxMLastRow) {
						maxMLastRow = subM.z;
						LRow = subL.z;
					} 
					if (subM.w >= maxMLastRow) {
						maxMLastRow = subM.w;
						LRow = subL.w;
					}															
				}				
			}			
		}
	
		if (maxMLastCol > maxMLastRow) 
			distArray[pairIndex] = (float) (LCol >> 16) / ((LCol << 16) >> 16);	
		else	
			distArray[pairIndex] = (float) (LRow >> 16) / ((LRow << 16) >> 16);

#if DEBUG			
			if (maxMLastCol > maxMLastRow)
				printf("%d: (%d-%d, %d-%d) %.3f (%d/%d)\n", pairIndex, readIndex1, length1, readIndex2, length2, distArray[pairIndex], (LCol >> 16), ((LCol << 16) >> 16));
			else
				printf("%d: (%d-%d, %d-%d) %.3f (%d/%d)\n", pairIndex, readIndex1, length1, readIndex2, length2, distArray[pairIndex], (LRow >> 16), ((LRow << 16) >> 16));
#endif					
			
	}	
	else
		distArray[pairIndex] = 1.1;	
}

void launchGenKernel_full(cudaStream_t stream, size_t blocksPerGrid, size_t threadsPerBlock, int *d_pairArray, float *d_distArray, int maxDigitalLen, int numPairs) {
	
	genKernel_full<<<blocksPerGrid, threadsPerBlock, 0, stream>>>(d_pairArray, d_distArray, maxDigitalLen, numPairs);
		
}


__device__ void computeSubMatrix_band(unsigned short read1, unsigned short read2,
		char *ULD, int *M, unsigned int *L, char4 &subULD, int4 &subM, uint4 &subL,
		int numRows, int numCols, int i, int j, int kBand, int endI, int endJ, int &maxMLastCol, uint &LCol, int &maxMLastRow, uint &LRow)
{
	bool isMatched;
	unsigned short symbol1, symbol2;
	int k, rowIndex, colIndex;
	int tempM, tempML, tempAL, M_up, M_left, M_diag, alignedScore;
	char subULD0;
	int subM0, diagonalM;
	uint subL0, diagonalL;

	// Initialize only for the first row of the first rowBlock
	if (i == 0)
	{				
		subULD0 = 0;
		subM0 = 0;		
		subULD = make_char4(0,0,0,0);
		subM = make_int4(0,0,0,0);

		colIndex = j * SUBMAT_DIM;
		if (colIndex <= kBand)
			subL0 = colIndex;
		++colIndex;
		if (colIndex <= kBand)
			subL.x = colIndex;
		++colIndex;
		if (colIndex <= kBand)
			subL.y = colIndex;
		++colIndex;
		if (colIndex <= kBand)
			subL.z = colIndex;
		++colIndex;
		if (colIndex <= kBand)
			subL.w = colIndex;		
	}

	for (k = 1; k <= numRows; ++k)
	{
		diagonalM = subM0;
		diagonalL = subL0;
		rowIndex = i * SUBMAT_DIM + k;	
		if (j == 0)
		{
			if (rowIndex <= kBand)
			{
				subULD0 = 0;
				subM0 = 0;
				subL0 = rowIndex;
			}
		}
		else
		{
			// copy from the whole matrix to the first column of the submatrix
			colIndex = j * SUBMAT_DIM;
			if (insideBand(colIndex,rowIndex,kBand))
			{
				subULD0 = ULD[rowIndex];
				subM0 = M[rowIndex];
				subL0 = L[rowIndex];						
			}
		}

		// extract the symbol from the unsigned short
		// shift the 2 bits to the left end
		symbol1 = read1 << ((k - 1) * 2);
		// shift the 3 bits all the way to the right end
		symbol1 = (symbol1 >> 14);
		
		colIndex = j * SUBMAT_DIM + 1;

		if (insideBand(colIndex,rowIndex,kBand))
		{
			symbol2 = read2 >> 14;

			if (symbol1 == symbol2)
			{
				alignedScore = MATCH;
				isMatched = 1;
			}
			else
			{
				alignedScore = MISMATCH;
				isMatched = 0;
			}

				
			if (insideBand(colIndex,rowIndex-1,kBand))
				M_up = subM.x + (subULD.x & 1) * GO + (subULD.x >> 2 & 1) * GE;
			else
				M_up = -1000000;

			if (insideBand(colIndex-1,rowIndex,kBand))
				M_left = subM0 + (subULD0 & 1) * GO + (subULD0 >> 1 & 1) * GE;
			else
				M_left = -1000000;
				
			M_diag = diagonalM + alignedScore;
			tempM = max(max(M_diag, M_up), M_left);
			
			if (tempM == M_diag) {
				subULD.x = 1; // 001
				tempAL = ((diagonalL << 16) >> 16) + 1;
				tempML = ((diagonalL>> 16) - isMatched) + 1;				
			}
			else if (tempM == M_up) {
				subULD.x = 4; // 100
				tempAL = ((subL.x << 16) >> 16) + 1;
				tempML = (subL.x >> 16) + 1;  
			}
			else {
				subULD.x = 2; // 010
				tempAL = ((subL0 << 16) >> 16) + 1;
				tempML = (subL0 >> 16) + 1;
			}

			diagonalM = subM.x;
			diagonalL = subL.x;

			subM.x = tempM;
			subL.x = ((tempML<<16) | (tempAL & 0xffff));	

			if (j == endJ-1 && numCols == 1 && subM.x >= maxMLastCol) {
				maxMLastCol = subM.x;
				LCol = subL.x;
			}	
			
			if (i == endI-1 && k == numRows && subM.x >= maxMLastRow) {
				maxMLastRow = subM.x;
				LRow = subL.x;
			} 				
		} else {
			diagonalM = subM.x;
			diagonalL = subL.x;
		}
			
		++colIndex;

		if (insideBand(colIndex,rowIndex,kBand) && numCols > 1) {
			symbol2 = read2 << 2;
			symbol2 = (symbol2 >> 14);

			if (symbol1 == symbol2)
			{
				alignedScore = MATCH;
				isMatched = 1;
			}
			else
			{
				alignedScore = MISMATCH;
				isMatched = 0;
			}

			if (insideBand(colIndex,rowIndex-1,kBand))
				M_up = subM.y + (subULD.y & 1) * GO + (subULD.y >> 2 & 1) * GE;
			else
				M_up = -1000000;
				
			if (insideBand(colIndex-1,rowIndex,kBand))
				M_left = subM.x + (subULD.x & 1) * GO + (subULD.x >> 1 & 1) * GE;
			else
				M_left = -1000000;
				
			M_diag = diagonalM + alignedScore;
			tempM = max(max(M_diag, M_up), M_left);
			
			if (tempM == M_diag) {
				subULD.y = 1; // 001
				tempAL = ((diagonalL << 16) >> 16) + 1;
				tempML = ((diagonalL>> 16) - isMatched) + 1;				
			}
			else if (tempM == M_up) {
				subULD.y = 4; // 100
				tempAL = ((subL.y << 16) >> 16) + 1;
				tempML = (subL.y >> 16) + 1;  				
			}
			else {
				subULD.y = 2; // 010
				tempAL = ((subL.x << 16) >> 16) + 1;
				tempML = (subL.x >> 16) + 1;				
			}
		
			diagonalM = subM.y;
			diagonalL = subL.y;

			subM.y = tempM;
			subL.y = ((tempML<<16) | (tempAL & 0xffff));	

			if (j == endJ-1 && numCols == 2 && subM.y >= maxMLastCol) {
				maxMLastCol = subM.y;
				LCol = subL.y;
			}
			if (i == endI-1 && k == numRows && subM.y >= maxMLastRow) {
				maxMLastRow = subM.y;
				LRow = subL.y;
			} 					
		} else {
			diagonalM = subM.y;
			diagonalL = subL.y;
		}
			
		++colIndex;

		if (insideBand(colIndex,rowIndex,kBand) && numCols > 2) {
			symbol2 = read2 << 4;
			symbol2 = (symbol2 >> 14);

			if (symbol1 == symbol2)
			{
				alignedScore = MATCH;
				isMatched = 1;
			}
			else
			{
				alignedScore = MISMATCH;
				isMatched = 0;
			}

			if (insideBand(colIndex,rowIndex-1,kBand))
				M_up = subM.z + (subULD.z & 1) * GO + (subULD.z >> 2 & 1) * GE;
			else
				M_up = -1000000;	
				
			if (insideBand(colIndex-1,rowIndex,kBand))
				M_left = subM.y + (subULD.y & 1) * GO + (subULD.y >> 1 & 1) * GE;
			else
				M_left = -1000000;
				
			M_diag = diagonalM + alignedScore;
			tempM = max(max(M_diag, M_up), M_left);
			
			if (tempM == M_diag) {
				subULD.z = 1; // 001
				tempAL = ((diagonalL << 16) >> 16) + 1;
				tempML = ((diagonalL>> 16) - isMatched) + 1;
			}
			else if (tempM == M_up) {
				subULD.z = 4; // 100
				tempAL = ((subL.z << 16) >> 16) + 1;
				tempML = (subL.z >> 16) + 1;  				
			}
			else {
				subULD.z = 2; // 010
				tempAL = ((subL.y << 16) >> 16) + 1;
				tempML = (subL.y >> 16) + 1;
			}
		
			diagonalM = subM.z;
			diagonalL = subL.z;

			subM.z = tempM;
			subL.z = ((tempML<<16) | (tempAL & 0xffff));

			if (j == endJ-1 && numCols == 3 && subM.z >= maxMLastCol) {
				maxMLastCol = subM.z;
				LCol = subL.z;
			}	
			if (i == endI-1 && k == numRows && subM.z >= maxMLastRow) {
				maxMLastRow = subM.z;
				LRow = subL.z;
			} 		
		} else {
			diagonalM = subM.z;
			diagonalL = subL.z;
		}
			
		++colIndex;

		if (insideBand(colIndex,rowIndex,kBand) && numCols > 3) {
			symbol2 = read2 << 6;
			symbol2 = (symbol2 >> 14);

			if (symbol1 == symbol2)
			{
				alignedScore = MATCH;
				isMatched = 1;
			}
			else
			{
				alignedScore = MISMATCH;
				isMatched = 0;
			}

			if (insideBand(colIndex,rowIndex-1,kBand))
				M_up = subM.w + (subULD.w & 1) * GO + (subULD.w >> 2 & 1) * GE;
			else
				M_up = -1000000;
				
			if (insideBand(colIndex-1,rowIndex,kBand))
				M_left = subM.z + (subULD.z & 1) * GO + (subULD.z >> 1 & 1) * GE;
			else
				M_left = -1000000;
				
			M_diag = diagonalM + alignedScore;
			tempM = max(max(M_diag, M_up), M_left);
					
			if (tempM == M_diag) {
				subULD.w = 1; // 001
				tempAL = ((diagonalL << 16) >> 16) + 1;
				tempML = ((diagonalL>> 16) - isMatched) + 1;				
			}
			else if (tempM == M_up) {
				subULD.w = 4; // 100
				tempAL = ((subL.w << 16) >> 16) + 1;
				tempML = (subL.w >> 16) + 1;  				
			}
			else {
				subULD.w = 2; // 010
				tempAL = ((subL.z << 16) >> 16) + 1;
				tempML = (subL.z >> 16) + 1;
			}
					
			subM.w = tempM;
			subL.w = ((tempML<<16) | (tempAL & 0xffff));	
			if (j < endJ - 1) // is not last col block
			{
				// copy the last row of the submatrix to the whole matrix
				ULD[rowIndex] = subULD.w;
				M[rowIndex] = subM.w;
				L[rowIndex] = subL.w;
			}	

			if (j == endJ-1 && subM.w >= maxMLastCol) {
				maxMLastCol = subM.w;
				LCol = subL.w;
			}		
			if (i == endI-1 && k == numRows && subM.w >= maxMLastRow) {
				maxMLastRow = subM.w;
				LRow = subL.w;
			} 									
		} 
	}
}


__global__ void genKernel_band (int* pairArray, float *distArray, int maxDigitalLen, int numPairs, int band)
{
	int pairIndex = blockIdx.x * blockDim.x + threadIdx.x;

	char ULD[MAX_READ_LEN + 1]; // up left diagonal matrix
	int M[MAX_READ_LEN + 1];
	uint L[MAX_READ_LEN + 1]; // combination of ML (first 16 bits) and AL (last 16 bits))

	char4 subULD;
	int4 subM;
	uint4 subL;

	unsigned short read1, read2, readIndex1, readIndex2, length1, length2;
	int i, j, numRows, numCols, kBand, endI, endJ;
	int maxMLastRow=-1000000, maxMLastCol=-1000000;
	uint LRow, LCol;
	
	// For each ID pair
	if (pairIndex < numPairs)
	{
		readIndex1 = pairArray[pairIndex*2];
		readIndex2 = pairArray[pairIndex*2+1];

		length1 = tex2D(seqTexRef,  maxDigitalLen * (readIndex1 & 15) + (maxDigitalLen-1), readIndex1/16);
		length2 = tex2D(seqTexRef, maxDigitalLen * (readIndex2 & 15) + (maxDigitalLen-1), readIndex2/16);

		// Recurrence relations
		// Process row by row
		// input file is sorted by sequence length
		// hence length1 is always smaller than length2
		//kBand = min ((int)length1, (int)length2) / band;
		kBand = (length1 + length2) / (band * 2);

		endJ = (length2 + SUBMAT_DIM-1)/SUBMAT_DIM;
		endI = (length1 + SUBMAT_DIM-1)/SUBMAT_DIM;	
		
		for (j = 0; j < endJ; ++j)
		{
			if (length2 % SUBMAT_DIM != 0 && j == endJ - 1) // last col block
				numCols = length2 % SUBMAT_DIM;
			else
				numCols = SUBMAT_DIM;		

			for (i = 0; i < endI; ++i)
			{	
				if (length1 % SUBMAT_DIM != 0 && i == endI - 1) // last col block 
					numRows = length1 % SUBMAT_DIM;
				else
					numRows = SUBMAT_DIM;					
									
				read2 = tex2D(seqTexRef, maxDigitalLen * (readIndex2 & 15) + j, readIndex2/16);			
				read1 = tex2D(seqTexRef, maxDigitalLen * (readIndex1 & 15) + i, readIndex1/16);
				// compute submatrices
				computeSubMatrix_band(read1, read2, ULD, M, L, subULD, subM, subL, numRows, numCols,i, j, kBand, endI, endJ, maxMLastCol, LCol, maxMLastRow, LRow);
					
			}
		}	
		
		if (maxMLastCol > maxMLastRow)
			distArray[pairIndex] = (float) (LCol >> 16) / ((LCol << 16) >> 16);		
		else 	
			distArray[pairIndex] = (float) (LRow >> 16) / ((LRow << 16) >> 16);	

#if DEBUG			
			if (maxMLastCol > maxMLastRow)
				printf("%d: (%d-%d, %d-%d) %.3f (%d/%d)\n", pairIndex, readIndex1, length1, readIndex2, length2, distArray[pairIndex], (LCol >> 16), ((LCol << 16) >> 16));
			else
				printf("%d: (%d-%d, %d-%d) %.3f (%d/%d)\n", pairIndex, readIndex1, length1, readIndex2, length2, distArray[pairIndex], (LRow >> 16), ((LRow << 16) >> 16));
#endif	
		
	}	
	else
		distArray[pairIndex] = 1.1;			
}

void launchGenKernel_band(cudaStream_t stream, size_t blocksPerGrid, size_t threadsPerBlock, int * d_pairArray, float *d_distArray, int maxDigitalLen, int numPairs, int band) {
	
	genKernel_band<<<blocksPerGrid, threadsPerBlock, 0, stream>>>(d_pairArray,  d_distArray, maxDigitalLen, numPairs, band);
	
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



