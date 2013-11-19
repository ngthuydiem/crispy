#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <iostream>
#include <map>
#include <vector>
#include <list>
#include <cmath>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <sys/time.h>
#include <sys/resource.h>
using namespace std;

#define BUF_SIZE 				8192
#define EPSILON 				0.00001
#define UPPER_BOUND_NUM_EDGES 	240000000
#define INVALID_ID 				1000000000
	
struct DistPair
{   
	unsigned int idX;
	unsigned int idY;
	unsigned int fileId;
	
	DistPair(unsigned int x, unsigned int y, unsigned int id): idX(x),idY(y),fileId(id) 
	{
	}
}; // end of struct DistPair

struct Link {
	unsigned int numEdges;
	float sumEdges;	
	Link() {
		numEdges = 0;
		sumEdges = 0.0f;
	}	
};

struct TreeNode {
	unsigned int ID;
	TreeNode* topParent;
	TreeNode* parent;
	TreeNode* left;
	TreeNode* right;	
	float dist;
	unsigned int numMembers;
	unordered_map<unsigned int, Link> linkMap;	
	
	TreeNode () {
		ID=INVALID_ID;
		topParent=0;
		parent=0;
		left=0;
		right=0;
		dist=0;
		numMembers=0;		
	}

	TreeNode (unsigned int id)
	{   
		ID=id;
		topParent=this;
		parent=0;
		left=0;
		right=0;
		dist=0;
		numMembers=1;			
	}	
}; // end of struct TreeNode

typedef unordered_map<unsigned int, Link> LinkMap;
typedef LinkMap::iterator LinkMapIter;

typedef unordered_set<TreeNode*> NodeSet;
typedef NodeSet::iterator NodeSetIter;

typedef list<TreeNode*> NodeList;
typedef NodeList::iterator NodeListIter;

typedef list<unsigned int> IDList;
typedef IDList::iterator IDListIter;

// global variables
TreeNode** leaves;
vector<unsigned int> freq;
unsigned int newId;
unsigned long long totalNumEdges = 0, totalUnlinked = 0;

vector<TreeNode*> vActiveNodes;
vector<TreeNode*> vExactNodes;
vector<float> vExactDist;	
vector<float> vInexactDist;

bool allLoaded = false;
float lambda = 0.0f;

void help()
{
	cout<<"-i <inFileName>"<<endl;
	cout<<"-n <numFiles>"<<endl;
	cout<<"-h help description."<<endl;
}

void getOptions(int argc, char** argv, string &inFileName, unsigned int &numFiles, float &stepSize, float &endLevel, unsigned long long &matrixSize)
{   
	for(int i=0; i<argc; ++i)
	{   
		if(strcmp("-i",argv[i])==0)
			inFileName.assign(argv[i+1]);
		if(strcmp("-n",argv[i])==0) {		
			numFiles=atoi(argv[i+1]);
		}	
		if(strcmp("-s",argv[i])==0) {		
			stepSize=atof(argv[i+1]);
		}	
		if(strcmp("-e",argv[i])==0) {		
			endLevel=atof(argv[i+1]);
		}	
		if(strcmp("-m",argv[i])==0) {		
			matrixSize=atol(argv[i+1]);
		}	
		if(strcmp("-h",argv[i])==0)
		{   
			help();
			exit(0);
		}
	} 

	if (inFileName.length() < 2)
	{
		cout<<"Please specify the input file name"<<endl;
		help();
		exit(0);   
	}
}

void getDistNameList(string inFileName, vector<string> &pairNameVector, vector<string>& distNameVector, unsigned int numFiles)
{   
	string tempStr;
	char buf[1000];

	for (unsigned int fileId=0; fileId<numFiles; fileId++) {
		sprintf(buf, "_%d", fileId);

		tempStr = inFileName;
		tempStr.append(".npair");

		tempStr.append(string(buf));	
		pairNameVector.push_back(tempStr);

		tempStr = inFileName;
		tempStr.append(".ndist");
		
		tempStr.append(string(buf));			
		distNameVector.push_back(tempStr);
	}	
}

unsigned int loadFile(FILE* pairFile, FILE* distFile, unsigned int * pairArray, float * distArray, bool & EOFTag)
{
	size_t readSize;

	readSize = fread(pairArray, sizeof(unsigned int), BUF_SIZE * 2, pairFile);	
	readSize = fread(distArray, sizeof(float), BUF_SIZE , distFile);	

	if (readSize < BUF_SIZE)
		EOFTag = true;

	return readSize;
}

bool loadAPair(int & readSize, int & index, bool & EOFTag,  unsigned int & idX, unsigned int & idY, float & dist, FILE* pairFile, FILE* distFile, unsigned int * pairArray, float * distArray)	
{
	if (readSize == 0 && !EOFTag) {
		// load from file
		readSize = loadFile(pairFile, distFile, pairArray, distArray, EOFTag);
		index = 0;
		if (EOFTag)
		{   			
			fclose(pairFile);
			fclose(distFile);
		}
	}
	
	if (readSize > 0) {
		// return value in array
		idX = pairArray[index*2];
		idY = pairArray[index*2+1];
		dist = distArray[index];
		--readSize;
		++index;
		//cout << readSize << "\t" << index << "\t" << dist << endl;
		return true;
	}
	else
		return false;
}

void addFreq(unsigned int numReads)
{
	TreeNode* aNode;
	int tempFreq;

	for(unsigned int i=0;i<numReads;++i)
		if(leaves[i]->parent!=0)
		{
			tempFreq=freq[leaves[i]->ID];
			if(tempFreq < 2)
				continue;
			else
			{
				tempFreq--;
				aNode=leaves[i];
				while(aNode->parent!=0)
				{
					aNode->numMembers += tempFreq;
					aNode = aNode->parent;
				}   
			} 
		}
}

void emptyTree(NodeSet roots)
{   
	TreeNode *tempNode = 0, *parentNode = 0;
	NodeSetIter setIter;
	NodeList nodeList;
	NodeListIter listIter;

	for(setIter=roots.begin(); setIter!=roots.end(); ++setIter)
	{   
		tempNode=0;
		parentNode=0;
		if(*setIter!=0)
		{   
			nodeList.push_front(*setIter);
			while (nodeList.size()!=0)
			{   
				listIter=nodeList.begin();
				tempNode=(*listIter);
				nodeList.pop_front();
				
				if (tempNode->right==0 && tempNode->left==0)
				{   
					parentNode=tempNode->parent;					
					if (parentNode->right->ID==tempNode->ID)
						parentNode->right = 0;
					else
						parentNode->left=0;
					delete tempNode;
					tempNode=0;
				}
				else
				{   					
					if(tempNode->right!=0)
						nodeList.push_front(tempNode->right);
					if(tempNode->left!=0)
						nodeList.push_front(tempNode->left);
				}


			}
		}
		nodeList.clear();
	}
}

int printClusters(NodeSet roots, IDList orphanNodes,
	string clusterListName, string clusterName, 
	float stepSize, float endLevel)
{   
	TreeNode *tempNode = 0;
	NodeSetIter setIter;

	NodeList nodeList, tempList;
	NodeListIter nodeIt, tempIt;
	
	IDList OTU;
	IDListIter it;
	
	unsigned int size, numOTUs;
	FILE *clusterListFile, *clusterFile;

	clusterListFile = fopen(clusterListName.c_str(),"wb");
	clusterFile = fopen(clusterName.c_str(),"wb");
	if(clusterListFile == NULL|| clusterFile == NULL)
	{   
		cout << "Cannot open output files. Skipped" << endl;
		return 0;
	}
	printf("\n");

	// for each distance level
	for(float distLevel=stepSize; distLevel<endLevel || fabs(distLevel-endLevel) < EPSILON; distLevel+=stepSize)
	{   
		numOTUs = 0;
		nodeList.clear();
		
		// extract the valid nodes for each distance level
		for(setIter=roots.begin(); setIter!=roots.end(); ++setIter)
		{   
			tempNode=0;
			if(*setIter != 0)
			{   
				if((*setIter)->dist < distLevel || fabs((*setIter)->dist-distLevel) < EPSILON)
				{   
					nodeList.push_front(*setIter);
					continue;
				}

				tempList.push_front(*setIter);
				while (tempList.size()!=0)
				{   
					tempIt=tempList.begin();
					tempNode=(*tempIt);
					tempList.pop_front();

					if (tempNode->left->dist < distLevel || fabs(tempNode->left->dist-distLevel) < EPSILON)
						nodeList.push_front(tempNode->left);						
					else
						tempList.push_front(tempNode->left);

					if (tempNode->right->dist < distLevel || fabs(tempNode->right->dist-distLevel) < EPSILON)
						nodeList.push_front(tempNode->right);
					else
						tempList.push_front(tempNode->right);					
				}
			}
			tempList.clear();
		}

		fprintf(clusterListFile," %.3f ", distLevel);
		fprintf(clusterFile," %.3f ", distLevel);
		
		// write the nodeList to file
		tempList.clear();
		for(nodeIt=nodeList.begin(); nodeIt!=nodeList.end(); ++nodeIt)
		{   
			// clean up and initialize
			fprintf(clusterFile,"|");
			tempNode=0;			
			size=0;
			OTU.clear();
			
			tempList.push_front(*nodeIt);
			
			while(tempList.size()!=0)
			{   
				tempIt=tempList.begin();
				tempNode=(*tempIt);
				tempList.pop_front();
				
				if(tempNode->left==0 && tempNode->right==0)
				{   
					OTU.push_back(tempNode->ID);
					size+=tempNode->numMembers;
				}				
				if (tempNode->right!=0)
					tempList.push_front(tempNode->right);
				if(tempNode->left!=0 )
					tempList.push_front(tempNode->left);
				
			}
			tempList.clear();					
			// print to clusterFile
			it=OTU.begin();
			fprintf(clusterFile,"%u",(*it));
			++it;
			for(;it!=OTU.end(); ++it)
				fprintf(clusterFile," %u",(*it));
			
			fprintf(clusterListFile, "%d ", size);	
			++numOTUs;			
		}
		for (it = orphanNodes.begin(); it != orphanNodes.end(); ++it) {
			fprintf(clusterFile,"|%u",(*it));
			fprintf(clusterListFile, "1 ");
		}
		numOTUs += orphanNodes.size();
		fprintf(clusterFile,"|\n");
		fprintf(clusterListFile, "\n");		
		printf("Dist: %.3f. numOTUs: %u. numSingletons: %lu\n", distLevel, numOTUs, orphanNodes.size());
	}
	
	printf("\n");
	OTU.clear();    
	fclose(clusterListFile);
	fclose(clusterFile);
	return 1;
}

unsigned int loadDistFile(FILE* distFile, bool &EOFTag)
{
	size_t readSize;

	float distArray[BUF_SIZE];
	
	readSize = fread(distArray, sizeof(float), BUF_SIZE , distFile);	

	if (readSize < BUF_SIZE)
		EOFTag = true;

	return readSize;
}

unsigned int loadFreq(string freqFileName)
{
	unsigned int tempFreq;
	unsigned int count = 0;
	FILE *freqFile;
	char buf[BUF_SIZE];

	if ((freqFile=fopen(freqFileName.c_str(),"r"))==NULL)
	{
		fprintf(stderr,"Cannot find frequency file %s\n", freqFileName.c_str());
		exit(1);
	}
	while (fscanf(freqFile,"%s %u",buf,&tempFreq)!=EOF)	{
		freq.push_back(tempFreq);
		++count;
	}

	fclose(freqFile);
	return count;
}


void updateTopParent(TreeNode* &aNode, TreeNode* &parentNode) {
	if (aNode == NULL)
		return;
	aNode->topParent = parentNode;
	updateTopParent(aNode->left, parentNode);
	updateTopParent(aNode->right, parentNode);
}

// idX > idY
float computeDist(TreeNode* nodeX, TreeNode* nodeY, float endLevel) {

	float dist = 1.0f;
	unsigned int totalLinks, numUnlinked;
	float sumEdges;
	
	unsigned int idX = nodeX->ID;
	unsigned int idY = nodeY->ID;

	if (nodeX->linkMap.find(idY) != nodeX->linkMap.end()) {
		totalLinks = nodeX->numMembers * nodeY->numMembers;
		sumEdges = nodeX->linkMap[idY].sumEdges;
		numUnlinked = totalLinks - nodeX->linkMap[idY].numEdges;

		if (numUnlinked == 0 || allLoaded) {
			dist = (sumEdges + endLevel*numUnlinked) / totalLinks;
			if (dist < vExactDist[idX]) {
				vExactDist[idX] = dist;
				vExactNodes[idX] = nodeY;
			}
		}		
		else {
			dist = (sumEdges + lambda*numUnlinked) / totalLinks;
			if (dist < vInexactDist[idX]) 
				vInexactDist[idX] = dist;									
		}		
	} 				

	return dist;		
}

void updateMin(TreeNode* nodeX, float endLevel) {
	
	TreeNode *tempNode;	
	unsigned int tempId;
	vector<TreeNode*>::iterator iter;
	LinkMap::iterator mapIter;

	vExactDist[nodeX->ID] = 1.0f;
	
	if (nodeX != 0 && !nodeX->linkMap.empty()) 
	{
		for (mapIter = nodeX->linkMap.begin(); mapIter != nodeX->linkMap.end(); ++mapIter) 
		{
			tempId = mapIter->first;
			tempNode = vActiveNodes[tempId];
			
			if (nodeX->ID > tempId)	
				computeDist(nodeX, tempNode, endLevel);												
			else
				computeDist(tempNode, nodeX, endLevel);												
		}
	}
}


// idX > idY
void addNode(TreeNode* &nodeX, TreeNode* &nodeY, TreeNode* &parentNode, float newDist, float endLevel)
{   
	TreeNode *tempNode;
	unsigned int tempId, idX, idY;
	LinkMapIter mapIter;
	
	idX = nodeX->ID;
	idY = nodeY->ID;
	
	// set child nodes for the parent TreeNode
	parentNode->left=nodeX;
	parentNode->right=nodeY;

	// set parent TreeNode for the child nodes	
	nodeX->parent=parentNode;
	nodeY->parent=parentNode;
	
	updateTopParent(nodeX, parentNode);
	updateTopParent(nodeY, parentNode);

	// set dist between two child nodes
	parentNode->dist=newDist;
	parentNode->numMembers = nodeX->numMembers + nodeY->numMembers;
	
	// pass the link map of nodeX to the parentNode	
	nodeX->linkMap.erase(idY);	

	--totalNumEdges;

	parentNode->linkMap = nodeX->linkMap;	

	// for each item in the link map of nodeY
	for(mapIter = nodeY->linkMap.begin(); mapIter != nodeY->linkMap.end(); ++mapIter) {
		tempId = mapIter->first;
		parentNode->linkMap[tempId].numEdges += mapIter->second.numEdges;		
		parentNode->linkMap[tempId].sumEdges += mapIter->second.sumEdges;															
	}				
	
	nodeX->linkMap.clear();
	nodeY->linkMap.clear();				
	totalNumEdges = totalNumEdges - (nodeX->linkMap.size() + nodeY->linkMap.size());	


	// for nodes with larger ID than nodeY that hold the links to nodeY
	for(unsigned int i = idY + 1; i < idX; ++i)
	{
		tempNode = vActiveNodes[i];			
		if (tempNode && (tempNode->linkMap.find(idY) != tempNode->linkMap.end())) {
			parentNode->linkMap[i].sumEdges += tempNode->linkMap[idY].sumEdges;
			parentNode->linkMap[i].numEdges += tempNode->linkMap[idY].numEdges;
			
			tempNode->linkMap.erase(idY);	
			--totalNumEdges;										
		}
	}	
				
	// for nodes with larger ID than nodeX that hold the links to nodeX
	for(unsigned int i = idX + 1; i < newId; ++i)
	{
		tempNode = vActiveNodes[i];		
		if (tempNode) {
					
			if (tempNode->linkMap.find(idX) != tempNode->linkMap.end()) {
				parentNode->linkMap[i].sumEdges += tempNode->linkMap[idX].sumEdges;
				parentNode->linkMap[i].numEdges += tempNode->linkMap[idX].numEdges;
			
				tempNode->linkMap.erase(idX);	
				--totalNumEdges;																	
			}
				
			if (tempNode->linkMap.find(idY) != tempNode->linkMap.end()) {
				parentNode->linkMap[i].sumEdges += tempNode->linkMap[idY].sumEdges;
				parentNode->linkMap[i].numEdges += tempNode->linkMap[idY].numEdges;
			
				tempNode->linkMap.erase(idY);	
				--totalNumEdges;					
			}
		}
	}				
	
	totalNumEdges = totalNumEdges + parentNode->linkMap.size();
	
	vActiveNodes[idX] = vActiveNodes[idY] = 0;
	vExactNodes[idX] = vExactNodes[idY] = 0;
	vExactDist[idX] = vExactDist[idY] = 1.0f;
	vInexactDist[idX] = vInexactDist[idY] = 1.0f;
	
	for(mapIter = parentNode->linkMap.begin(); mapIter != parentNode->linkMap.end(); ++mapIter) {
		tempId = mapIter->first;	
		tempNode = vActiveNodes[tempId];		

		computeDist(parentNode, tempNode, endLevel);													
		
		if (vExactNodes[tempId] == nodeX || vExactNodes[tempId] == nodeY) {
			//cout << "update node " << tempNode->ID << endl;
			updateMin(tempNode, endLevel);									
		}							
	}		
}	

// merge two clusters
void merge(TreeNode * nodeX, TreeNode * nodeY, float dist, float endLevel)
{   
	TreeNode *parentNode=NULL;	
	LinkMapIter mapIter;
	
	parentNode = new TreeNode(newId);	
	vActiveNodes[newId] = parentNode;	
	
	if (nodeX->ID > nodeY->ID)	
		addNode(nodeX, nodeY, parentNode, dist, endLevel);	
	else
		addNode(nodeY, nodeX, parentNode, dist, endLevel);	
	
	// DEBUG
	//cout << totalNumEdges << "\t" << dist << "\t" << nodeX->ID << "\t" << nodeY->ID << "\t" << parentNode->ID << "\r";		
	++ newId;	
}

// absorb a TreeNode
void absorb(unsigned int anIdX, unsigned int anIdY, float dist)
{   
	TreeNode *nodeX, *nodeY;
	unsigned int idX, idY;

	nodeX = leaves[anIdX]->topParent;
	nodeY = leaves[anIdY]->topParent;

	idX = nodeX->ID;
	idY = nodeY->ID;		
	
	if (idX > idY) {		
		
		if (nodeX->linkMap.find(idY) == nodeX->linkMap.end()) 
			++totalNumEdges;
			
		nodeX->linkMap[idY].sumEdges += dist;
		++nodeX->linkMap[idY].numEdges;					
	} 
	
	else if (idY > idX) {
	
		if (nodeY->linkMap.find(idX) == nodeY->linkMap.end()) 
			++totalNumEdges;
			
		nodeY->linkMap[idX].sumEdges += dist;
		++nodeY->linkMap[idX].numEdges;							
	}
}

void lookAhead(unsigned int anIdX, unsigned int anIdY, float dist)
{   
	TreeNode *nodeX, *nodeY;
	unsigned int idX, idY;

	nodeX = leaves[anIdX]->topParent;
	nodeY = leaves[anIdY]->topParent;

	idX = nodeX->ID;
	idY = nodeY->ID;		

	if (idX > idY) {		
	
		if (nodeX->linkMap.find(idY) == nodeX->linkMap.end()) {
			++totalUnlinked;
		}
		else {
			nodeX->linkMap[idY].sumEdges += dist;
			++ nodeX->linkMap[idY].numEdges;				
		}
		
	} else if (idY > idX) {
	
		if (nodeY->linkMap.find(idX) == nodeY->linkMap.end()) {
			++totalUnlinked;
		}
		else {
			nodeY->linkMap[idX].sumEdges += dist;
			++ nodeY->linkMap[idX].numEdges;				
		}
	
	}
}

void updateAllMin(float endLevel) {
	
	TreeNode *nodeX, *nodeY;	
	unsigned int tempId;
	vector<TreeNode*>::iterator iter;
	LinkMap::iterator mapIter;

	//cout << "update all min..." << endl;
	fill(vExactDist.begin(), vExactDist.end(), 1.0f);
	fill(vInexactDist.begin(), vInexactDist.end(), 1.0f);
	
	for (iter = vActiveNodes.begin(); iter != vActiveNodes.end(); ++iter) {
		nodeX = *iter;
		if (nodeX != 0 && !nodeX->linkMap.empty()) {
			for (mapIter = nodeX->linkMap.begin(); mapIter != nodeX->linkMap.end(); ++mapIter) {
				tempId = mapIter->first;
				nodeY = vActiveNodes[tempId];
				
				if (nodeX->ID > nodeY->ID)
					computeDist(nodeX, nodeY, endLevel);										
				else
					computeDist(nodeY, nodeX, endLevel);										
			}
		}
	}
}

int main(int argc, char* argv[])
{   
	struct rlimit r;
	getrlimit(RLIMIT_NOFILE, &r);
	cout << "current rlimit: " << r.rlim_cur << endl;
	r.rlim_cur = 2048;
	setrlimit(RLIMIT_NOFILE, &r);
	cout << "change rlimit to: " << r.rlim_cur << endl;
	
	struct timeval startTime, endTime;
	
	gettimeofday(&startTime, NULL);
	
	unsigned int numFiles=1;
	vector<string> pairNameVector, distNameVector;

	FILE **distFileList, **pairFileList;
	unsigned int ** inPairArray;
	float ** inDistArray;

	string inFileName = "";
		
	printf("\n----------------------------------------------------------------------\n");
	printf("                  AVERAGE CLUSTERING GENETIC DISTANCES                  \n");

	float stepSize = 0.01f, endLevel = 0.20f;
	unsigned long long matrixSize=0, maxNumEdges=0;		

	getOptions(argc, argv, inFileName, numFiles, stepSize, endLevel, matrixSize);
	getDistNameList(inFileName, pairNameVector, distNameVector, numFiles);

	if(pairNameVector.size()==0)
	{   
		cout<<"Error: No distance file loaded."<<endl;
		exit(-1);
	}
	
	unsigned int fileId;
	
	pairFileList=(FILE**)malloc(sizeof(FILE*)*pairNameVector.size());
	distFileList=(FILE**)malloc(sizeof(FILE*)*distNameVector.size());
	
	if(distFileList==NULL || pairFileList==NULL)
	{   
		cout<<"Error: Not enough memory" << endl;
		exit(-1);
	}

	for(fileId = 0; fileId < numFiles; ++fileId)
	{
		pairFileList[fileId]=fopen(pairNameVector[fileId].c_str(),"rb");
		if(pairFileList[fileId]==NULL)
		{   
			cout<<"Error: Cannot open file" << pairNameVector[fileId].c_str() << endl;
			exit(-1);
		}
	}
	
	for(fileId = 0; fileId < numFiles; ++fileId)
	{
		distFileList[fileId]=fopen(distNameVector[fileId].c_str(),"rb");
		if(distFileList[fileId]==NULL)
		{   
			cout<<"Error: Cannot open file" << distNameVector[fileId].c_str() << endl;
			exit(-1);
		}
	}

	if(numFiles!=distNameVector.size() || numFiles <= 0)
	{   
		cout<<"Error: invalid number of files!EXIT..."<<endl;
		exit(-1);
	}
	cout<<"Use "<<numFiles<<" distance file(s)."<<endl;
	cout<<"endLevel: " << endLevel << " stepSize: " << stepSize << endl;
						
	unsigned long long totalNumPairs = 0;

	multimap<float, DistPair> nodeMap;
	multimap<float, DistPair>::iterator iter;
	unsigned int idX, idY;
	float dist;
	inPairArray = (unsigned int **) malloc(sizeof(unsigned int*) * numFiles);
	inDistArray = (float **) malloc(sizeof(float*) * numFiles);
	int * indices = (int*) malloc(sizeof(int) * numFiles);
	int * readSizes = (int*) malloc(sizeof(int) * numFiles);
	bool * EOFTags = (bool*) malloc(sizeof(bool) * numFiles);
	bool suc;

	for(fileId=0; fileId<numFiles; fileId++)
	{   
		// initialize
		inPairArray[fileId] = (unsigned int*) malloc(sizeof(unsigned int) * BUF_SIZE * 2);
		inDistArray[fileId] = (float*) malloc(sizeof(float) * BUF_SIZE);
		indices[fileId] = 0;
		readSizes[fileId] = 0;
		EOFTags[fileId] = false;				

		// add the first pair of each file to the nodeMap
		suc = loadAPair(readSizes[fileId], indices[fileId], EOFTags[fileId], idX, idY, dist, pairFileList[fileId], distFileList[fileId], inPairArray[fileId], inDistArray[fileId]);				
	
		if (suc)
			nodeMap.insert(pair<float, DistPair>(dist,DistPair(idX,idY,fileId)));
	}		

				
	unsigned int i, iteration = 0;
	unsigned int numReads=0;
		
	unsigned int minExactIndex;
	TreeNode *nodeX=0, *nodeY=0;
	vector<float>::iterator minExactIter;
	float minExactDist, minInexactDist;
	LinkMap::iterator mapIter;	
						
	string freqFileName = inFileName;
	freqFileName.append(".frq");
	numReads=loadFreq(freqFileName);

	vActiveNodes.resize(2*numReads+1);
	vExactNodes.resize(2*numReads+1);
	vExactDist.resize(2*numReads+1);	
	vInexactDist.resize(2*numReads+1);

	leaves=(TreeNode**)malloc(sizeof(TreeNode*)*numReads);

	if(leaves==0)
	{   
		cout<<"Error: Not enough memory" << endl;
		exit(-1);
	}
	
	for(i = 0; i < numReads; ++i) {
		vActiveNodes[i] = leaves[i] = new TreeNode(i);	
		vExactNodes[i] = 0;
	}
	
	for(i = numReads+1; i < 2*numReads+1; ++i) {
		vActiveNodes[i] = 0;
		vExactNodes[i] = 0;
	}
	
	newId = numReads;	

	fill(vExactDist.begin(), vExactDist.end(), 1.0f);
	fill(vInexactDist.begin(), vInexactDist.end(), 1.0f);
	
	float alpha = sqrt(numReads)/4;
	float beta = log10(alpha)/2;
	if (beta < 1.0)
		beta = 1.1f;
	maxNumEdges = min((float)UPPER_BOUND_NUM_EDGES, matrixSize/alpha);
	cout << "numReads: " << numReads << "\tmaxNumEdges: " << maxNumEdges << endl;
	cout << "alpha: " << alpha << "\tbeta: " << beta << endl;

	//while (!nodeMap.empty())
	while (!allLoaded)
	{				
		// the clustering can't continue and the edges exceed the capacity of the RAM
		if (totalNumEdges < maxNumEdges) 
		{
			while(totalNumEdges <= maxNumEdges && !nodeMap.empty())
			{   			
				// get the first item in the nodeMap
				iter = nodeMap.begin();
				fileId = iter->second.fileId;        

				// write to output
				idX = iter->second.idX;
				idY =  iter->second.idY;
				dist = iter->first;
						
				if (dist < endLevel || fabs(dist-endLevel) < EPSILON)				
					absorb(idX, idY, dist);	
			
				// remove the current item from the nodeMap
				nodeMap.erase(iter);
			
				suc = loadAPair(readSizes[fileId], indices[fileId], EOFTags[fileId], idX, idY, dist, pairFileList[fileId], distFileList[fileId], inPairArray[fileId], inDistArray[fileId]);				
				if (suc) 
					nodeMap.insert(pair<float, DistPair>(dist,DistPair(idX,idY,fileId)));	
			}
			
			if (dist > lambda) 
				lambda = dist;														
		}	
		else 
		{	
			maxNumEdges *= beta;
			cout << "new maxNumEdges: " << maxNumEdges << "\tnewId: " << newId << "\tlambda: " << lambda << endl;
			
			if (maxNumEdges > UPPER_BOUND_NUM_EDGES)
			{
				/*	
				cout << "LOOK AHEAD: current node: " << newId << "\tnum edges: " << totalNumEdges << endl;	
				
				while(!nodeMap.empty())
				{   
					// get the first item in the nodeMap
					iter = nodeMap.begin();
					fileId = iter->second.fileId;        

					// write to output
					idX = iter->second.idX;
					idY =  iter->second.idY;
					dist = iter->first;

					lookAhead(idX, idY, dist);	
					// remove the current item from the nodeMap
					nodeMap.erase(iter);
			
					suc = loadAPair(readSizes[fileId], indices[fileId], EOFTags[fileId], idX, idY, dist, pairFileList[fileId], distFileList[fileId], inPairArray[fileId], inDistArray[fileId]);
					if (suc) 
						nodeMap.insert(pair<float, DistPair>(dist,DistPair(idX,idY,fileId)));	
				}	
				*/
				endLevel = lambda;				
				allLoaded = true;
				cout << "new endLevel: " << endLevel << endl;
			}							
		}	
				
		if (nodeMap.empty())
			allLoaded = true;					
			
		updateAllMin(endLevel);		
				
		minInexactDist = *min_element(vInexactDist.begin(), vInexactDist.end());
		minExactIter = min_element(vExactDist.begin(),vExactDist.end());
		minExactDist = *minExactIter;
		
		//cout << "lambda_" << iteration << " = " << lambda << "\t" << lambda/2 << "\t" << minInexactDist << endl;
		
		while ((minExactDist < 1.0f) && (minExactDist < minInexactDist || fabs(minExactDist-minInexactDist) < EPSILON))
		{				
			minExactIndex = minExactIter-vExactDist.begin();
				
			nodeX = vActiveNodes[minExactIndex]->topParent;				
			nodeY = vExactNodes[minExactIndex]->topParent;
				
			merge(nodeX, nodeY, minExactDist, endLevel);	
			
			minInexactDist = *min_element(vInexactDist.begin(), vInexactDist.end());
			minExactIter = min_element(vExactDist.begin(),vExactDist.end());
			minExactDist = *minExactIter;
		}	
		
		//cout << "cannot progress: " << newId << "\t" << totalNumEdges << "\t" << minInexactDist << "\t" << minExactDist << "\t" << minExactDist - lambda/2 << "\n";	
		++iteration;							
	} 
		
	cout << "DONE!" << endl;
	cout << "current node: " << newId << "\tnum unlinked: " << totalUnlinked << endl;	

	// get root nodes and orphan nodes
	NodeSet roots;
	IDList orphanNodes;	
	TreeNode *aLeaf = 0;
     	
	for(i=0; i< numReads; ++i)
	{   
		aLeaf = leaves[i];
		if(aLeaf->parent==0)  // find nodes with no parent
			orphanNodes.push_back(i); 				
		else
			roots.insert(aLeaf->topParent); 										
	}
	addFreq(numReads);

 	// print output to files
 	
	string clusterListName, clusterName;

	clusterListName=inFileName;
	clusterListName.append(".Cluster_List");
	clusterName=inFileName;
	clusterName.append(".Cluster");

	printClusters(roots, orphanNodes, clusterListName, clusterName, stepSize, endLevel);			
	
	// clear memory
	emptyTree(roots);	
	freq.clear();
	roots.clear();	
	orphanNodes.clear();	
	free(leaves);
	vActiveNodes.clear();
	vExactNodes.clear();
	vExactDist.clear();
	
	// clean up
	for(fileId=0; fileId<numFiles; ++fileId) {
		free(inDistArray[fileId]);
		free(inPairArray[fileId]);
	}
	free(inDistArray);
	free(inPairArray);
	free(indices);
	free(readSizes);
	free(EOFTags);
	free(pairFileList);
	free(distFileList);	
	
	gettimeofday(&endTime, NULL);	
	long elapsedTime = (endTime.tv_sec - startTime.tv_sec) * 1000u + (endTime.tv_usec - startTime.tv_usec) / 1.e3 + 0.5;
	
	printf("totalNumPairs: %llu\n", totalNumPairs);
	printf("Time taken: %.3f s\n", elapsedTime/1.e3);
	printf("\n----------------------------------------------------------------------\n");
	return 0;
}

