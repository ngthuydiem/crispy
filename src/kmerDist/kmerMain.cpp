/***********************************************
* # Copyright 2011. Nguyen Thuy Diem
* # Contact: Kristin
* #          thuy1@e.ntu.edu.sg
* #
* # GPL 3.0 applies.
* #
* ************************************************/

#include "kmerMain.h"

int main(int argc, char* argv[]) {

	struct timeval startTime, endTime;
	
	gettimeofday(&startTime, NULL);
	
	READ *readArray = NULL;	

	string inFileName, pairFileName, distFileName;
	bool hasDistances = false, useGPU = true, useMPI = false;
	float threshold = -1;	

	int numThreads = 1, numReads, maxLen = 0, arrayDim, K = 6;
	// open the output file
	FILE* pairFile = NULL, *distFile = NULL;
		
	getCommandOptions(argc, argv, inFileName, threshold, hasDistances, numThreads, useGPU, useMPI, K);

	pairFileName = inFileName;
	pairFileName.append(".kpair");
	distFileName = inFileName;
	distFileName.append(".kdist");
	
	// read from input files
	numReads = readFile(inFileName, readArray, maxLen, K);		
	
	if (numReads > MAX_NUM_READS || maxLen > MAX_READ_LEN) {
		printf("Error: unsupported numReads: %d or maxLen: %d. Exit...", numReads, maxLen);
		exit(-1);
	}
	
	if (threshold < 0 || threshold > 1) {		
		if (numReads > 500)
			threshold = min(2/log2((double)numReads),1.0);	
		else
			threshold = min(2/log10((double)numReads),1.0);	
	}
		
	arrayDim = (int)(4 * pow(8.0, floor(log10((double)numReads))-1));	 

	if (arrayDim > 1536)
		arrayDim = 1536;		
	
	omp_set_num_threads(numThreads);

	// form the k-tuple array for each sequence from the file
#pragma omp parallel for schedule(static,1)
	for (int i = 0; i < numReads; ++i) {
		readArray[i].formTuples(K);
		readArray[i].sortTuples();
	}
		
	if (useGPU && useMPI) {		

		// Initialize MPI state
		MPI_Init(NULL, NULL);
		    		 
		// Get our MPI node number and node count
    	int commSize, commRank, len;
		char name[BUF_SIZE];
	    MPI_Comm_size(MPI_COMM_WORLD, &commSize);
	    MPI_Comm_rank(MPI_COMM_WORLD, &commRank);
		MPI_Get_processor_name(name, &len);
			
		char temp[5];
		sprintf(temp,"_%d",commRank);
		pairFileName.append(temp);
		pairFile = fopen(pairFileName.c_str(), "wb");
		if (hasDistances) {
			distFileName.append(temp);
			distFile = fopen(distFileName.c_str(), "wb");	
		}

		if (commRank == 0)
		{
			printf("\n----------------------------------------------------------------------\n");
			printf("                       COMPUTE KMER DISTANCES                           \n");
			printf("File name: %s. K = %d\n\n", inFileName.c_str(), K);
			printf("numReads: %d, maxLen: %d, threshold: %.2f\n  ", numReads, maxLen, threshold);
			printf("USE MPI. numCPUs: %d\n\n", commSize);
		}		
		computeKmerDist_MPI(readArray, pairFile, distFile, hasDistances, numReads, maxLen, threshold, arrayDim, commRank, commSize, K);			
		fclose(pairFile);
		if (hasDistances)
			fclose(distFile);

		if (commRank == 0) {
			gettimeofday(&endTime, NULL);	
			long elapsedTime = (endTime.tv_sec - startTime.tv_sec) * 1000u + (endTime.tv_usec - startTime.tv_usec) / 1.e3 + 0.5;
	
			printf("Time taken: %.3f s\n", elapsedTime/1.e3);
			printf("\n----------------------------------------------------------------------\n");
		}
	
		MPI_Finalize();

	}	
	else { 
		printf("\n----------------------------------------------------------------------\n");
		printf("                       COMPUTE KMER DISTANCES                           \n");
		printf("File name: %s. K = %d\n\n", inFileName.c_str(), K);
		printf("numReads: %d, maxLen: %d, threshold: %.2f\n  ", numReads, maxLen, threshold);
			
		pairFile = fopen(pairFileName.c_str(), "wb");
		if (hasDistances)
			distFile = fopen(distFileName.c_str(), "wb");	

		if (useGPU) {
			printf("USE GPU. numStreams: %d\n", NUM_STREAMS);
			computeKmerDist_CUDA(readArray, pairFile, distFile, hasDistances, numReads, maxLen, threshold, arrayDim, K);					
		}
		else {
			printf("USE CPU. numThreads: %d\n", numThreads);
			computeKmerDist_CPU(readArray, pairFile, distFile, hasDistances, numReads, threshold, arrayDim, K);
		}
		
		fclose(pairFile);
		if (hasDistances)
			fclose(distFile);
	
		gettimeofday(&endTime, NULL);	
		long elapsedTime = (endTime.tv_sec - startTime.tv_sec) * 1000u + (endTime.tv_usec - startTime.tv_usec) / 1.e3 + 0.5;
	
		printf("Time taken: %.3f s\n", elapsedTime/1.e3);
		printf("\n----------------------------------------------------------------------\n");
	}		
	
	return 0;
}

// implementation of READ

void READ::initialize(int readId, const char* seq, int K) {

	this->id = readId;
	this->length = strlen(seq);
	this->numTuples = this->length - K + 1;
	if (this->length < 2)
		cout << "Error: empty sequence!\n" << endl;
	this->sequence = (char *) malloc(sizeof(char) * (strlen(seq) + 1));
	strcpy(this->sequence, seq);	
}

void READ::finalize() {

	if (this->tuples!=NULL)
		free(this->tuples);
	this->tuples = NULL;
}

// digital forms of 4 symbols: A:00, G:01, T:10, C:11

void READ::formTuples(int K) {

	// calculate the number of tuples for each sequence
	int symbolIndex, tupleIndex;
	this->tuples = (unsigned short*) calloc(this->numTuples, sizeof(unsigned short));

	// for each symbol in the sequence
	for (symbolIndex = 0, tupleIndex = 0; symbolIndex < this->length; symbolIndex++) {
		if (symbolIndex == 0)
			tuples[tupleIndex] = 0;
		if (symbolIndex >= K) {
			++tupleIndex;
			this->tuples[tupleIndex] = (this->tuples[tupleIndex - 1] << (2 * (9
				- K)));
			this->tuples[tupleIndex] = (this->tuples[tupleIndex] >> (2
				* (8 - K)));
		} else {
			this->tuples[tupleIndex] = (this->tuples[tupleIndex] << 2);
		}

		switch (this->sequence[symbolIndex]) {
		case 'A': // 00
			break;
		case 'G': // 01
			this->tuples[tupleIndex] |= 1;
			break;
		case 'T': // 10
			this->tuples[tupleIndex] |= (1 << 1);
			break;
		case 'C': // 11
			this->tuples[tupleIndex] |= (1 << 1);
			this->tuples[tupleIndex] |= 1;
			break;
		default: 
			break;
		}
	}
}

void READ::sortTuples() {

	qsort(this->tuples, this->numTuples, sizeof(short int), compareTwoTuples);
}

void writeToFile(FILE *pairFile, FILE *distFile, bool hasDistances, float *h_distArray, int stageX, int stageY, int arrayDim, float threshold, unsigned long long & totalNumPairs) 
{

	int i, h=0, row, col, rowOffset, colOffset;	
	float dist;	
	int pairArray[BUF_SIZE*2];
	float distArray[BUF_SIZE];	
	
	int arraySize = arrayDim * arrayDim;
	
	rowOffset = stageX * arrayDim;
	colOffset = stageY * arrayDim;				

	// write result to output file
	for (i = 0; i < arraySize; ++i) 
	{
		row = rowOffset + (int)i / arrayDim;
		col = colOffset + (int)i % arrayDim;	
		dist = h_distArray[i];	
		if (dist < threshold || fabs(dist-threshold) < EPSILON)
		{
			pairArray[h*2] = row;
			pairArray[h*2+1] = col;
			if (hasDistances)
				distArray[h] = dist;	
			++h;
			if (h == BUF_SIZE) {
				fwrite(pairArray, sizeof(int),BUF_SIZE*2, pairFile);
				if (hasDistances)
					fwrite(distArray, sizeof(float),BUF_SIZE, distFile);
				h = 0;
				totalNumPairs += BUF_SIZE;
			}
		}							
	}
	
	if (h > 0) {
		fwrite(pairArray, sizeof(int),h*2, pairFile);
		if (hasDistances)
			fwrite(distArray, sizeof(float),h, distFile);
		totalNumPairs += h;
	}
}

void Trag_reverse_eq(int index, int N, int& row, int& col) 
{ 
   row = 0; 
   int keyafter; 
   do 
   { 
       row++; 
       keyafter = row * N - (row - 1) * row / 2; 
   } while (index >= keyafter); 
   row--; 
   col = N - keyafter + index; 
} 

void help()
{
	cout << "<-i inFileName>: FASTA file " << endl;	
	cout << "[-t threshold] : value range: 0.0 to 1.0, default: 0.5" << endl;
	cout << "[-n numThreads]: value range: 1 to maxNumProcs, default: 1" << endl;
	cout << "[-d]		    : output kmer distances " << endl;
	cout << "[-c]           : use CPU instead of GPU" << endl;
	cout << "[-m]           : use MPI version for GPU cluster" << endl;
	cout << "[-h]		    : help" << endl;
}

void getCommandOptions(int argc, char* argv[], string &inFileName, float &threshold, bool &hasDistances, int &numThreads, bool &useGPU, bool &useMPI, int &K)
{

	int numProcs = omp_get_num_procs();					
	if (numProcs > 2)
		numThreads = numProcs/2;
	else
		numThreads = 1;
		
	// process input arguments
	for (int i = 0; i < argc; ++i)
	{
		if (strcmp("-i", argv[i]) == 0) 
		{
			inFileName.assign(argv[i + 1]);		
			// check input file
			if (inFileName.length() < 2)
			{
				help();
				exit(-1);				
			}		
		}
		if (strcmp("-t", argv[i]) == 0) {
			threshold = (float)atof(argv[i + 1]);								
		}
		if (strcmp("-n", argv[i]) == 0) {
			numThreads = atoi(argv[i + 1]);				
			if (numThreads < 1 || numThreads > numProcs)
			{
				cout << "Warning: invalid number of threads (-n option). Set to " << numThreads << "." << endl;
			}
		}
		if (strcmp("-k", argv[i]) == 0)
			K = atoi(argv[i + 1]);		
		if (strcmp("-d", argv[i]) == 0)
			hasDistances = true;
		if (strcmp("-c", argv[i]) == 0)
			useGPU = false;
		if (strcmp("-m", argv[i]) == 0)
			useMPI = true;
		if (strcmp("-h", argv[i]) == 0) {
			help();
			exit(-1);
		}
	}	
}


int compareTwoTuples(const void* t1, const void* t2) {

	return (*(short int*) t1 - *(short int*) t2);

}

//--------------------- Trim the buff, de-whitespaces --------------------
void removeNewLine(string &line) {
	size_t found;
	found = line.find_last_not_of(" \n\r\t\v\f\b");
	if (found != string::npos)
		line = line.substr(0, found + 1);
}

void removeInvalidChar(string &line) {
	size_t found;
	found = line.find_last_of("ATGCatgc");
	if (found != string::npos)
		line = line.substr(0, found + 1);
}


//--------------------------------------------------------------------------
int readFile(string inFileName, READ* &readArray, int &maxLen, int K) {

	FILE * inFile;
	char buffer[BUF_SIZE];

	inFile = fopen(inFileName.c_str(), "rb");
	if (inFile == NULL) {
		cout << "Error: cannot open read file: " << inFileName << " Exit..." << endl;
		exit(-1);
	}

	readArray = (READ *) malloc(DEFAULT_NUM_READS * sizeof(READ));

	if (readArray == NULL) {
		cout << "Error: not enough memory. Require extra: " << DEFAULT_NUM_READS * sizeof(READ) / (1024 * 1024) << "MB. Exit..." << endl;		
		fclose(inFile);
		exit(-1);
	}

	string tempStr;
	int readIndex = 0;
	string seq = "", label;
	int stopFlag = 0;

	// initialise the reading process
	fgets(buffer, BUF_SIZE, inFile);
	tempStr.clear();
	tempStr.assign(buffer);
	if (tempStr[0] != '>') {
		cout << "Error: read file is not in FASTA Format. Exit..."<< endl;
		fflush(stdin);
		exit(-1);
	} else {
		removeNewLine(tempStr);
		label.clear();
		label = tempStr;
	}

	// read from input file
	while (!feof(inFile)) {
		buffer[0] = '$';
		fgets(buffer, BUF_SIZE, inFile);
		tempStr.clear();
		tempStr.assign(buffer);
		removeNewLine(tempStr);
		if (tempStr[0] != '>' && tempStr[0] != '$') {
			removeInvalidChar(tempStr);
			seq += tempStr;
			continue;
		} else if (seq.length() == 0) {
			if (buffer[0] == '$') {
				stopFlag = 1;
				break;
			}
			continue;
		} else {
			removeInvalidChar(seq);
			// if there are more sequences than the buffer length
			if (readIndex + 1 > DEFAULT_NUM_READS) 
				readArray = (READ *) realloc(readArray, (readIndex + 1) * sizeof(READ));
			if (readArray == NULL) {
				cout << "Error: not enough memory. Exit..."<< endl;
				fclose(inFile);
				exit(-1);
			}

			if (seq.length() > maxLen)
				maxLen = seq.length();
				
			readArray[readIndex].initialize(readIndex, seq.c_str(), K);
			readIndex++;
			seq.clear();

			if (buffer[0] == '$') {
				stopFlag = 1;
				break;
			}
			removeNewLine(tempStr);
			label.clear();
			label = tempStr;
			continue;
		}
	}

	if (stopFlag == 0 && seq.length() != 0) {
		if (readIndex + 1 > DEFAULT_NUM_READS)
			readArray = (READ *) realloc(readArray, (readIndex + 1) * sizeof(READ));
		if (readArray == NULL) {
			cout << "Error: not enough memory. Exit..."<< endl;
			fclose(inFile);
			exit(-1);
		}
		removeInvalidChar(seq);
		readArray[readIndex].initialize(readIndex, seq.c_str(), K);
		readIndex++;
	}
	
	fclose(inFile);
	return readIndex;
}

