/***********************************************
* # Copyright 2011. Zejun Zheng
* # Contact: zheng_zejun@sics.a-star.edu.sg
* #
* # GPL 3.0 applies.
* #
* ************************************************/

#include <string>
#include <iostream>
#include <fstream>
#include <map>
#include <list>
#include <string.h>
#include <cstdlib>
#include <math.h>
#include <time.h>
#include <algorithm>
using namespace std;

#define BUF_SIZE 	4096

#define min(x,y) ((x) < (y) ? (x) : (y))
#define max(x,y) ((x) > (y) ? (x) : (y))

typedef struct read
{
	unsigned id;
	string seq;
	string label;
	int freq;

	void initialize(string aSeq, string aLabel, int aFreq,
			unsigned int aId)
	{
		aLabel.erase(remove_if(aLabel.begin(), aLabel.end(), (int(*)(int))isspace), aLabel.end());		
		this->seq = aSeq;
		this->label = aLabel;
		this->freq = aFreq;
		this->id = aId;
	}
	bool operator <(const struct read &aRead) const
	{
		if (aRead.seq.length() == this->seq.length())
		{
			return strcmp(this->seq.c_str(), aRead.seq.c_str()) < 0;
		}
		else
			return this->seq.length() < aRead.seq.length();
	}
	bool operator ==(const struct read &aRead) const
	{
		return (this->seq.length() == aRead.seq.length() && strcmp(
				this->seq.c_str(), aRead.seq.c_str()) == 0);
	}
} READ;

void deleteSpaces(string &aString)
{
	size_t found;
	found = aString.find_last_not_of("\n\r \t\v\f\b"); 
	if (found != string::npos)
		aString = aString.substr(0, found + 1);
}
void deleteSpacesInSeq(string &aSeq)
{
	size_t found = string::npos;
	found = aSeq.find_last_of("nNATGCatgc"); 
	if (found != string::npos)
		aSeq = aSeq.substr(0, found + 1);
}

// Load reads from data file, delete repeated reads
int readFile(list<READ> &readList, string inFileName, char fileFormat)
{

	// open input file
	FILE * inFile;
	char buf[BUF_SIZE];
	inFile = fopen(inFileName.c_str(), "rb");
	if (inFile == NULL)
	{
		cout << "Can not open read file. Press any key to exit..." << endl;
		fflush(stdin);
		getchar();
		exit(-1);
	}
	string tempString;
	unsigned int readIndex = 0;
	string seq, label;
	string nSymbol = "Nn";
	READ tempRead;

	seq = "";
	int endFlag = 0;

	// Process a FASTA file
	if (fileFormat == 'a')
	{
		cout << "Format: FASTA." << endl;
		fgets(buf, BUF_SIZE, inFile);
		tempString.clear();
		tempString.assign(buf);
		if (tempString[0] != '>')
		{
			printf("Error: the read file is not in FASTA format. EXIT...\n");
			exit(-1);
		}
		else
		{
			deleteSpaces(tempString);
			label.clear();
			label = tempString;
		}
		while (!feof(inFile))
		{
			buf[0] = '$';
			fgets(buf, BUF_SIZE, inFile);
			tempString.clear();
			tempString.assign(buf);

			deleteSpaces(tempString);
			if (tempString[0] != '>' && tempString[0] != '$')
			{
				deleteSpacesInSeq(tempString);
				seq += tempString;
				continue;
			}
			else if (seq.length() == 0)
			{
				++readIndex;
				if (buf[0] == '$')
				{
					endFlag = 1;
					break;
				}
				continue;
			}
			else
			{
				
				deleteSpacesInSeq(seq);
				tempRead.initialize(seq, label, 1, readIndex);
				readList.push_back(tempRead);
				++readIndex;
				seq.clear();

				if (buf[0] == '$')
				{
					endFlag = 1;
					break;
				}
				deleteSpaces(tempString);
				label.clear();
				label = tempString;
				continue;
			}

		}
	}

	// process a FASTQ file
	else
	{
		cout << "Format: FASTQ." << endl;
		while (!feof(inFile))
		{
			buf[0] = '$';
			fgets(buf, BUF_SIZE, inFile);
			tempString.clear();
			tempString.assign(buf);
			if (feof(inFile) || buf[0] == '$')
			{
				endFlag = 1;
				break;
			}
			deleteSpaces(tempString);
			label.clear();
			label = tempString;
			fgets(buf, BUF_SIZE, inFile);

			tempString.clear();
			tempString.assign(buf);
			if (feof(inFile) || buf[0] == '$')
			{
				endFlag = 1;
				break;
			}
			deleteSpaces(tempString);
			seq.clear();
			seq = tempString;

			if (strcmp(seq.c_str(), "") == 0)
			{
				tempRead.initialize(seq, label, 0, readIndex);
				readList.push_back(tempRead);
				++readIndex;
				fgets(buf, BUF_SIZE, inFile);
				if (feof(inFile) || buf[0] == '$')
				{
					endFlag = 1;
					break;
				}
				fgets(buf, BUF_SIZE, inFile);
				if (feof(inFile) || buf[0] == '$')
				{
					endFlag = 1;
					break;
				}
				continue;
			}
			else
			{				
				tempRead.initialize(seq, label, 1, readIndex);
				readList.push_back(tempRead);
				++readIndex;
			}
			fgets(buf, BUF_SIZE, inFile);
			if (feof(inFile) || buf[0] == '$')
			{
				endFlag = 1;
				break;
			}
			fgets(buf, BUF_SIZE, inFile);
			if (feof(inFile) || buf[0] == '$')
			{
				endFlag = 1;
				break;
			}
		}
	}

	if (endFlag == 0 && seq.length() != 0)
	{
		deleteSpacesInSeq(seq);
		tempRead.initialize(seq, label, 1, readIndex);
		readList.push_back(tempRead);
		++readIndex;
	}
	fclose(inFile);
	return readIndex;
}

int mergeReads(list<READ> &readList)
{
	list<READ>::iterator it1, it2;
	int numMerges = 0;
	string seq1, seq2;

	// Sort the whole list
	readList.sort();

	for (it1 = readList.begin(); it1 != readList.end();)
	{
		if (it1->freq == 0)
		{
			++it1;
			continue;
		}
		seq1 = it1->seq;
		it2 = it1;
		++it2;
		while (it2->freq == 0 && it2 != readList.end())
			++it2;

		if (it2 == readList.end())
			break;
		seq2 = it2->seq;

		if (seq2 == seq1)
		{
			it2->freq += it1->freq;
			it1->freq = 0;
			++numMerges;
		}
		it1 = it2;
	}

	return numMerges;
}

int mergeSubstrings(list<READ> &readList)
{
	list<READ>::iterator it1, it2;
	int numMerges = 0;
	string seq1, seq2;

	// Sort the whole list
	readList.sort();

	for (it1 = readList.begin(); it1 != readList.end();)
	{
		if (it1->freq == 0)
		{
			++it1;
			continue;
		}
		seq1 = it1->seq;
		it2 = it1;
		++it2;
		while (it2->freq == 0 && it2 != readList.end())
			++it2;

		if (it2 == readList.end())
			break;
		seq2 = it2->seq;

		if (seq2.find(seq1) != string::npos)
		{
			it2->freq += it1->freq;
			it1->freq = 0;
			++numMerges;
		}
		it1 = it2;
	}

	return numMerges;
}

// Delete seqs with atypical lengths and not within 1SD from the mean length and 
// Delete seq that contain ambiguous nucleotides (N) (low quality reads)
int processReadLength(list<READ> &readList, double numSDAllowed)
{
	int seqIndex;
	double sumLen, aveLen, SD, diff, sumDiff, numReads;
	int numDiscards = 0;
	int minLen = 65535, maxLen = 0;

	sumLen = aveLen = SD = 0.0;
	numReads = sumDiff = 0.0;
	list<READ>::iterator it;

	for (it = readList.begin(); it != readList.end(); ++it)
	{
		if (it->freq == 0)
			continue;
		sumLen += (it->seq.length() + 0.0);
		++numReads;

	}
	if (numReads < 1)
	{
		cout << "Error in processReadLength"
				<< endl;
		return 0.0;
	}
	aveLen = sumLen/numReads;

	// Calculate standard deviation of the seq lengths
	for (it = readList.begin(); it != readList.end(); ++it)
	{
		if (it->freq == 0)
			continue;
		diff = it->seq.length() - aveLen;
		sumDiff += diff * diff;
		minLen = min(minLen, it->seq.length());
		maxLen = max(maxLen, it->seq.length());
	}
	sumDiff = sumDiff/numReads;
	SD = sqrt(sumDiff);
	
	cout << "minLen: " << minLen << "\t\tmaxLen: " << maxLen << endl;
	
	minLen = max(minLen, (int) (aveLen - numSDAllowed * SD));
	maxLen = min(maxLen, (int) (aveLen + numSDAllowed * SD));
	for (it = readList.begin(); it != readList.end(); ++it)
	{
		if (it->freq == 0)
			continue;
		if (it->seq.length() < minLen || it->seq.length() > maxLen)
		{
			it->freq = 0;
			++numDiscards;
			continue;
		}
		
		for (seqIndex = it->seq.length() - 1; seqIndex >= 0; seqIndex--)
		{
			if (it->seq[seqIndex] == 'n' || it->seq[seqIndex] == 'N') // remove reads which contain ambiguous nucleotides
			{
				it->freq = 0;
				++numDiscards;
				break;
			}
		}
	}
		
	cout << "aveLen: " << aveLen << "\t\tSD: " << SD << endl;
	cout << "minLen: " << minLen << "\t\tmaxLen: " << maxLen << endl;

	return numDiscards;
}

void SWLocalAlignment(const char* seqA, const char* seqB,
		int maxMismatchesAllowed, int &startIndex, int &endIndex) //-----0 offset
{
	int lengthA, lengthB;
	int i, j; //i---primer
	int **scoreMatrix, **startIndexA, **mismatchMatrix;
	int tempScore[3];
	int alignedScore;

	lengthA = strlen(seqA);
	lengthB = strlen(seqB);
	scoreMatrix = (int**) malloc(sizeof(int*) * (lengthA + 1));
	startIndexA = (int**) malloc(sizeof(int*) * (lengthA
			+ 1));
	mismatchMatrix = (int**) malloc(sizeof(int*) * (lengthA + 1));
	if (scoreMatrix == NULL || startIndexA == NULL
			|| mismatchMatrix == NULL)
	{
		cout << "Memory malloc error: 1" << endl;
		if (scoreMatrix != NULL)
			free(scoreMatrix);
		if (startIndexA != NULL)
			free(startIndexA);
		if (mismatchMatrix != NULL)
			free(mismatchMatrix);
		return;
	}

	for (i = 0; i < lengthA + 1; ++i)
	{
		scoreMatrix[i] = (int*) calloc(lengthB + 1, sizeof(int));
		startIndexA[i] = (int*) calloc(lengthB + 1,
				sizeof(int));
		mismatchMatrix[i] = (int*) calloc(lengthB + 1, sizeof(int));
		if (scoreMatrix[i] == NULL || startIndexA[i] == NULL
				|| mismatchMatrix[i] == NULL)
		{
			cout << "Memory malloc error: 2" << endl;
			if (scoreMatrix[i] != NULL)
				free(scoreMatrix);
			if (startIndexA[i] != NULL)
				free(startIndexA);
			if (mismatchMatrix[i] != NULL)
				free(mismatchMatrix);
			return;
		}
	}

	startIndex = endIndex = 0;
	int maxScore, startA, numMismatches, endA;
	maxScore = -9999;
	numMismatches = 999;
	startA = 0;
	endA = 0;
	for (i = 0; i < lengthA + 1; ++i)
		startIndexA[i][0] = i;
	for (i = 1; i < lengthA + 1; ++i)
	{
		for (j = 1; j < lengthB + 1; ++j)
		{
			// compute the temporary alignment scores
			alignedScore = (seqA[i - 1] == seqB[j - 1]) ? 4 : -1; // match = +4, mismatch = -1
			tempScore[0] = scoreMatrix[i - 1][j - 1] + alignedScore; // match or mismatch
			tempScore[1] = scoreMatrix[i - 1][j] - 2; // gap, gap penalty = -2
			tempScore[2] = scoreMatrix[i][j - 1] - 2; // gap

			// get the maximum value from 4 scores: temp[0], temp[1], temp[2] and 0
			scoreMatrix[i][j] = max(tempScore[0], tempScore[1]);
			scoreMatrix[i][j] = max(scoreMatrix[i][j], tempScore[2]);
			scoreMatrix[i][j] = max(scoreMatrix[i][j], 0);

			if (scoreMatrix[i][j] == 0)
			{
				startIndexA[i][j] = i - 1;
				mismatchMatrix[i][j] = 0;
			}
			else if (scoreMatrix[i][j] == tempScore[0])
			{
				startIndexA[i][j]
						= startIndexA[i - 1][j - 1];
				mismatchMatrix[i][j] = mismatchMatrix[i - 1][j - 1];
				if (alignedScore != 4)
					++mismatchMatrix[i][j];
			}
			else if (scoreMatrix[i][j] == tempScore[1])
			{
				startIndexA[i][j]
						= startIndexA[i - 1][j];
				mismatchMatrix[i][j] = mismatchMatrix[i - 1][j] + 1;
			}
			else
			{
				startIndexA[i][j] = startIndexA[i][j
						- 1];
				mismatchMatrix[i][j] = mismatchMatrix[i][j - 1] + 1;
			}
			if ((maxScore < scoreMatrix[i][j]) || (maxScore
					== scoreMatrix[i][j] && numMismatches
					> mismatchMatrix[i][j]))
			{
				maxScore = scoreMatrix[i][j];
				startA = startIndexA[i][j];
				endA = i;
				numMismatches = mismatchMatrix[i][j];
			}
		}
	}
	numMismatches += startA;
	numMismatches += abs(lengthA - endA);
	if (numMismatches <= maxMismatchesAllowed)
	{
		startIndex = startA;
		endIndex = endA;
	}
	else
	{
		startIndex = 0;
		endIndex = 0;
	}

	for (i = 0; i < lengthA + 1; ++i)
	{
		if (scoreMatrix[i] != NULL)
			free(scoreMatrix[i]);
		if (startIndexA[i] != NULL)
			free(startIndexA[i]);
		if (mismatchMatrix[i] != NULL)
			free(mismatchMatrix[i]);
	}

	free(scoreMatrix);
	free(mismatchMatrix);
	free(startIndexA);
}

// delete reads with more than one mismatch with the PCR primer at the beginning of a read
int processPrimer(list<READ> &readList,
		string primerSequence, int maxMismatchesAllowed)
{
	string tempSequence;
	int startIndex, endIndex;
	int numDiscards = 0;
	unsigned int i;
	int search_len;

	i = 0;
	for (list<READ>::iterator it = readList.begin(); it
			!= readList.end(); ++it)
	{
		if (it->freq == 0)
			continue;
		tempSequence.clear();
		tempSequence = it->seq;
		search_len = min(primerSequence.length() * 2, tempSequence.length());
		tempSequence = tempSequence.substr(0, search_len);
		SWLocalAlignment(primerSequence.c_str(), tempSequence.c_str(),
				maxMismatchesAllowed, startIndex, endIndex);
		if (startIndex == endIndex)
		{
			it->freq = 0;
			++numDiscards;
		}
		++i;
	}

	return numDiscards;
}

int writeToFile(list<READ> &readList, string outFileName,
		char fileFormat)
{
	string freqFileName = outFileName;
	freqFileName.append(".frq");
	
	FILE *outFile, *freqFile;
	outFile = fopen(outFileName.c_str(), "wb");
	freqFile = fopen(freqFileName.c_str(), "wb");
	if (outFile == NULL || freqFile == NULL) {
		printf("Error: cannot open output files. EXIT...");
		exit(-1);
	}
	
	int numReads = 0;

	for (list<READ>::const_iterator it = readList.begin(); it != readList.end(); ++it)
	{
		if (it->freq == 0)
			continue;

		if (fileFormat == 'a') {
			fprintf(outFile, "> %s\n", it->label.substr(1, it->label.length() - 1).c_str());
			fprintf(freqFile, "%s ", it->label.substr(1, it->label.length() - 1).c_str());
		}
		else {
			fprintf(outFile, "> %s\n", it->label.c_str());
			fprintf(freqFile, "%s ", it->label.c_str());
		}
			
		fprintf(outFile, "%s\n", it->seq.c_str());
		fprintf(freqFile, "%d\n", it->freq);
		++numReads;
		//cout << it->seq << endl;
	}
	
	fclose(outFile);
	fclose(freqFile);
	return numReads;
}

void help()
{
	cout << "<-i inFileName>" << endl;
	cout << "[-f input format a(FASTA) q(FASTQ), default: a" << endl;
	cout << "[-p primerFileName]" << endl;
	cout << "[-m number of mismatches to primers allowed, default 2]" << endl;
	cout << "[-s standard deviation from average length allowed, default 1.0]" << endl;
	cout << "[-h  help]" << endl;
}

void getCommandOptions(int argc, char** argv, string &inFileName, string &outFileName, char &fileFormat, 
int &primerCheckFlag, string &primerSequence, int &maxMismatches, int &lengthCheckFlag, double &SDAllowed)
{
	string primerFileName;
	ifstream primerFile;
	primerCheckFlag = 0;

	for (int i = 0; i < argc; ++i)
	{
	
		if (strcmp(argv[i], "-i") == 0)
		{
			inFileName.assign(argv[i + 1]);
			if (inFileName.length() < 2) {
				help();
				exit(-1);
			}						
		}
		
		if (strcmp(argv[i], "-o") == 0)
		{
			outFileName.assign(argv[i + 1]);
		}

		if (strcmp(argv[i], "-f") == 0)
		{
			if (strcmp(argv[i + 1], "q") == 0)
				fileFormat = 'q';
		}

		if (strcmp(argv[i], "-p") == 0)
		{			
			primerFileName.assign(argv[i + 1]);
			primerCheckFlag = 1;
		}

		if (strcmp(argv[i], "-m") == 0)
		{
			maxMismatches = atoi(argv[i + 1]);				
		}		

		if (strcmp(argv[i], "-s") == 0)
		{
			SDAllowed = atof(argv[i + 1]);
			lengthCheckFlag = 1;
		}
					
		if (strcmp(argv[i], "-h") == 0)
		{
			help();
			exit(-1);
		}
	}
	
	if (outFileName.length() < 2) {
		outFileName = inFileName;
		outFileName.append("_Clean");
	}
			
	if (primerCheckFlag == 1)
	{
		primerFile.open(primerFileName.c_str());
		if (primerFile.fail())
		{
			printf("Error: cannot open primer file. EXIT...");
			exit(-1);
		}
		getline(primerFile, primerSequence);
		deleteSpaces(primerSequence);
		primerFile.close();
	}
}

int main(int argc, char* argv[])
{
	char fileFormat = 'a'; //-f  option: a  or q
	int primerCheckFlag = 0, lengthCheckFlag =0; //-p  option: filename
	string primerSequence = "$"; //seq of primer
	int maxMismatches = 2; //-m  option: mismatches allowed
	double numSDAllowed = 1.0; //-s  option: standard deviation allowed
	string inFileName, outFileName; //-i  option: name and path of the input seq file
	clock_t startTime, endTime;
	int numReads, numReadsDeleted=0;

	startTime = clock();
	getCommandOptions(argc, argv, inFileName, outFileName, fileFormat,
			primerCheckFlag, primerSequence, maxMismatches, lengthCheckFlag, numSDAllowed);

	list<READ> readList;
	printf("\n--------------------------------------------------------------------\n");
	printf("\n                            PREPROCESS	               			  \n");
	
	numReads = readFile(readList, inFileName,fileFormat);	
	printf("File name: %s.\n", inFileName.c_str()); 	

	if (lengthCheckFlag) 
		numReadsDeleted += processReadLength(readList, numSDAllowed);
	
	if (primerCheckFlag)
		numReadsDeleted += processPrimer(readList, primerSequence, maxMismatches);

	numReadsDeleted += mergeReads(readList);

	printf("Total numReads: %d. Preprocessed numReads: %d\n", numReads, numReads- numReadsDeleted);

	writeToFile(readList, outFileName, fileFormat);

	readList.clear();
	endTime = clock();
	
	printf("\nTime taken: %.3f secs\n", ((double) (endTime - startTime)) / CLOCKS_PER_SEC);
	printf("\n--------------------------------------------------------------------\n");
	return 1;
}
