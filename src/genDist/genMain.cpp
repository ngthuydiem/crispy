/***********************************************
* # Copyright 2011. Thuy Diem Nguyen
* # Contact: thuy1@e.ntu.edu.sg
* #
* # GPL 3.0 applies.
* #
* ************************************************/

#include "genMain.h"

int main(int argc, char* argv[])
{
	struct timeval startTime, endTime;
	
	gettimeofday(&startTime, NULL);

	READ *readArray;	
	string readFileName, inFileName, pairFileName, distFileName;
	int numReads, numThreads, band=BAND, maxLen = 0;
	FILE *inPairFile;
	float threshold = -1;
	bool useGPU = true;
	
	getCommandOptions(argc, argv, readFileName, threshold, band, useGPU, numThreads);
	inFileName = readFileName;	
	inFileName.append(".kpair");
	pairFileName = readFileName;
	pairFileName.append(".npair");
	distFileName = readFileName;
	distFileName.append(".ndist");

	// Loading a sequence from file after trimming
	numReads = readFile(readFileName, readArray, maxLen);		
	
	if (numReads > MAX_NUM_READS || maxLen > MAX_READ_LEN) {
		printf("Error: unsupported numReads or maxLen. Exit...");
		
		cout << numReads << "\t" << maxLen << endl;
		exit(-1);
	}
	
	if (threshold < 0)
		threshold = 1/log((double)numReads);

#if SET_AFFINITY	
  int num_cpus = 8;
	cpu_set_t *cpusetp;
	cpusetp = CPU_ALLOC(num_cpus);
	if (cpusetp == NULL) {
		printf("!!!!error cpu_alloc\n");	
		exit(-1);
	}
	
	size_t size = CPU_ALLOC_SIZE(num_cpus);
	CPU_ZERO_S(size, cpusetp);
		
	int cpu;
	for (cpu = 0; cpu < num_cpus; cpu ++)
		CPU_SET_S(cpu, size, cpusetp);
	printf("CPU_COUNT() of set:    %d\n", CPU_COUNT_S(size, cpusetp));
	
	if (sched_setaffinity(0, size, cpusetp) == - 1) {
		printf("!!!!error setting affinity\n");	
		exit(-1);
	}
	
	int afin;	
	afin = sched_getaffinity(0, size, cpusetp);
	printf("Affinity: %d, 0x%x\n", afin, cpusetp);
#else
	cpu_set_t mask;	
	sched_getaffinity(0, sizeof(mask), &mask);
#endif
	
		
	printf( "\n----------------------------------------------------------------------\n");
	printf("                    COMPUTE GENETIC DISTANCES                          \n");
	printf("File name: %s.\n", readFileName.c_str());
	printf("numReads: %d. band: %d. threshold: %f\n\n", numReads, band, threshold);

	inPairFile = fopen(inFileName.c_str(), "rb");
	if (inPairFile == NULL) {
		cout << "Error: cannot open pair file: " << inFileName.c_str() << ". Exit..." << endl;
		exit(-1);
	}
			
	if (useGPU) {
		printf("USE GPU. numStreams: %d\n", NUM_STREAMS);
		computeGenDist_CUDA(inPairFile, pairFileName, distFileName, readArray, numReads, maxLen, threshold, band);	
	}
	else {
		omp_set_num_threads(numThreads);
		
		printf("USE CPU. numThreads: %d\n", numThreads);

		if (band > 1)
			computeGenDist_CPU_band(inPairFile, pairFileName, distFileName, readArray, numReads, maxLen, threshold, band, numThreads);	
		else
			computeGenDist_CPU_full(inPairFile, pairFileName, distFileName, readArray, numReads, maxLen, threshold, numThreads);	
	}

	fclose(inPairFile);
		
	gettimeofday(&endTime, NULL);	
	long elapsedTime = (endTime.tv_sec - startTime.tv_sec) * 1000u + (endTime.tv_usec - startTime.tv_usec) / 1.e3 + 0.5;

	printf("Time taken: %.3f s\n", elapsedTime/1.e3);
	printf("\n----------------------------------------------------------------------\n");

	return 0;
}

void* Malloc(size_t size)
{
	void *ptr;
	if ((ptr=malloc(size))==NULL)
	{
		fprintf(stderr,"Out of Memory!\n");
		exit(-1);
	}
  	return ptr;
}


int loadPairs(FILE* pairFile, int * &pairArray, int &EOFTag)
{
	size_t readSize;

	readSize = fread(pairArray, sizeof(int), NUM_PAIRS*2, pairFile);	

	if (readSize < NUM_PAIRS*2)
		EOFTag = 1;

	return readSize/2;
}

//------------------ The above matrix is cited from ESPRIT -------------------

void READ::initialize(int readId, const char* seq)
{

	this->id = readId;
	this->length = strlen(seq);
	if (this->length < 2)
		cout << "Warning: empty sequence - ID " << readId << endl;
	this->sequence = (char *) malloc(sizeof(char) * (strlen(seq) + 1));
	if (this->sequence == NULL) {
		cout << "Error: Not enough memory. Exit..." << endl;
		exit (-1);
	}	
	strcpy(this->sequence, seq);

}

void READ::release()
{
	if (this->sequence!=NULL)
		free(this->sequence);
	this->sequence = NULL;
}

//---------------------------------------------------------------
void deleteSpaces(string &line)
{
	size_t found;
	found = line.find_last_not_of("\n\r \t\v\f\b");
	if (found != string::npos)
		line = line.substr(0, found + 1);
}

void deleteSpacesInSeq(string &line)
{
	size_t found;
	found = line.find_last_of("nNATGCatgc");
	if (found != string::npos)
		line = line.substr(0, found + 1);
}

//--------------------------------------------------------------------------
int readFile(string inFileName, READ* &readArray, int &maxLen) {

	FILE * inFile;
	char buffer[BUF_SIZE];

	inFile = fopen(inFileName.c_str(), "rb");
	if (inFile == NULL) {
		cout << "Error: cannot open read file! EXIT..." << endl;
		fflush(stdin);
		exit(-1);
	}

	readArray = (READ *) Malloc(MAX_NUM_READS * sizeof(READ));

	string tempStr;
	int readIndex = 0;
	string seq = "", label;
	int stopFlag = 0;

	// initialise the reading process
	fgets(buffer, BUF_SIZE, inFile);
	tempStr.clear();
	tempStr.assign(buffer);
	if (tempStr[0] != '>') {
		cout << "Error: the read file is not in FASTA format! EXIT..."<< endl;
		fflush(stdin);
		exit(0);
	} else {
		deleteSpaces(tempStr);
		label.clear();
		label = tempStr;
	}

	// read from input file
	while (!feof(inFile)) {
		buffer[0] = '$';
		fgets(buffer, BUF_SIZE, inFile);
		tempStr.clear();
		tempStr.assign(buffer);
		deleteSpaces(tempStr);
		if (tempStr[0] != '>' && tempStr[0] != '$') {
			deleteSpacesInSeq(tempStr);
			seq += tempStr;
			continue;
		} else if (seq.length() == 0) {
			if (buffer[0] == '$') { // skip empty sequences
				stopFlag = 1;
				break;
			}
			continue;
		} else {
			deleteSpacesInSeq(seq);
			// if there are more sequences than the buffer length
			if (readIndex + 1 > BUF_SIZE) 
				readArray = (READ *) realloc(readArray, (readIndex + 1) * sizeof(READ));

			if (readArray == NULL) {
				cout << "Error: not enought memory! EXIT..."<< endl;
				fclose(inFile);
				exit(-2);
			}

			if (seq.length() > maxLen)
				maxLen = seq.length();
			readArray[readIndex].initialize(readIndex, seq.c_str());
			++readIndex;
			seq.clear();

			if (buffer[0] == '$') {
				stopFlag = 1;
				break;
			}
			deleteSpaces(tempStr);
			label.clear();
			label = tempStr;
			continue;
		}
	}

	if (stopFlag == 0 && seq.length() != 0) {
		if (readIndex + 1 > BUF_SIZE)
			readArray = (READ *) realloc(readArray, (readIndex + 1) * sizeof(READ));
		if (readArray == NULL) {
			cout << "Error: not enough memory! EXIT..."<< endl;
			fclose(inFile);
			exit(-2);
		}
		deleteSpacesInSeq(seq);
		readArray[readIndex].initialize(readIndex, seq.c_str());
		++readIndex;
	}
	fclose(inFile);
	return readIndex;
}


void convertToBinary(READ* &readArray, int numReads,
	unsigned short* &binaryReadArray, int maxDigitalLen)
{
	int readIndex = 0, binaryReadIndex = 0;
	int i = 0, j = 0;

	//--------------------------------------------------
	//    A: 00, G: 01, T: 10, C: 11                      |
	//   SUBMAT_DIM characters were stored in the same UINT value |
	//--------------------------------------------------

#pragma omp parallel for private (i, j, binaryReadIndex) schedule(static,1)
	for (readIndex = 0; readIndex < numReads; ++readIndex) // for each read
	{		
		for (i = 0; i < readArray[readIndex].length; i += SUBMAT_DIM) // for each character of a sequence
		{
			binaryReadIndex = readIndex * maxDigitalLen + i / SUBMAT_DIM;

			for (j = 0; j < SUBMAT_DIM && i + j < readArray[readIndex].length; ++j)
			{
				switch (readArray[readIndex].sequence[i + j])
				{
				case 'A':
					break;
				case 'G':
					binaryReadArray[binaryReadIndex] += (1 /*01*/<< (14 - 2 * j));
					break;
				case 'T':
					binaryReadArray[binaryReadIndex] += (2 /*10*/<< (14 - 2 * j));
					break;
				case 'C':
					binaryReadArray[binaryReadIndex] += (3 /*11*/<< (14 - 2 * j));
					break;
				default:
					break;
				}
			}
		}
		// Store the length at the end of each digital sequence
		binaryReadArray[readIndex * maxDigitalLen + maxDigitalLen - 1] = abs(readArray[readIndex].length);
	}
}


void help()
{	
	cout << "<-i inFileName>: FASTA file. Require: output from kmerDist (inFileName.kpair) " << endl;	
	cout << "[-t threshold] : value range: 0.0 to 1.0, default: 0.2" << endl;
	cout << "[-n numThreads]: value range: 1 to maxNumProcs, default: 1" << endl;
	cout << "[-b]			: use banded alignment, default: 10 i.e. 1/10 band width" << endl;
	cout << "[-c]           : use CPU instead of GPU" << endl;
	cout << "[-h]		    : help" << endl;
}

void getCommandOptions(int argc, char* argv[], string &readFileName, float &threshold, int &band, bool &useGPU, int &numThreads)
{
	numThreads = omp_get_num_procs();					
	
	for (int i = 0; i < argc; ++i) {
		if (strcmp("-i", argv[i]) == 0) {
			readFileName.assign(argv[i + 1]);
			if (readFileName.length() < 2) {
				help();
				exit(0);
			}
		}
		if (strcmp("-t", argv[i]) == 0) {
			threshold = (float)atof(argv[i + 1]);			
		}
		if (strcmp("-n", argv[i]) == 0) {
			numThreads = atoi(argv[i + 1]);				
			if (numThreads < 1)
			{
				cout << "Warning: invalid number of threads (-n option). Set to " << numThreads << "." <<  endl;
			}
		}
		if (strcmp("-b", argv[i]) == 0) {
			band = atoi(argv[i + 1]);
			// check distance threshold 
			if (band <= 0)
			{
				cout << "Warning: invalid band (-b option). Set to " << BAND << "."<< endl;
				band = BAND;
			}
		}
		if (strcmp("-c", argv[i]) == 0) 
			useGPU = false;			
	
		if (strcmp("-h", argv[i]) == 0) {
			help();
			exit(0);
		}
	}		
}


