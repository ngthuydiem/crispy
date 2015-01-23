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
#define UPPER_BOUND_NUM_EDGES 	200000000
#define INVALID_ID 				1000000000

struct TreeNode {
	unsigned int ID;
	TreeNode* topParent;
	TreeNode* parent;
	TreeNode* left;
	TreeNode* right;	
	float dist;
	unsigned int numMembers;
	
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


typedef unordered_set<TreeNode*> NodeSet;
typedef NodeSet::iterator NodeSetIter;

typedef list<TreeNode*> NodeList;
typedef NodeList::iterator NodeListIter;

typedef list<unsigned int> IDList;
typedef IDList::iterator IDListIter;

// global variables
TreeNode** leaves;
unsigned int newId;

void help()
{
	cout<<"-i <inFileName>"<<endl;
	cout<<"-c <cutoffFileName>"<<endl;	
	cout<<"-n <numReads>"<<endl;
	cout<<"-h help description."<<endl;
}

void getOptions(int argc, char** argv, string &inFileName, string &cutoffFileName, int &numReads)
{   
	numReads = 0;
	for(int i=0; i<argc; ++i)
	{   
		if(strcmp("-i",argv[i])==0)
			inFileName.assign(argv[i+1]);
		if(strcmp("-c",argv[i])==0)
			cutoffFileName.assign(argv[i+1]);		
		if(strcmp("-n",argv[i])==0) {		
			numReads=atoi(argv[i+1]);
		}	
		if(strcmp("-h",argv[i])==0)
		{   
			help();
			exit(0);
		}
	} 

	if (inFileName.length() < 2 || numReads < 1)
	{
		help();
		exit(0);   
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
	vector<float> cutoffs)
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

	vector<float>::iterator c;
	float distLevel;
	for(c = cutoffs.begin(); c != cutoffs.end(); c++)
	{   
		distLevel = *(c);
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

		fprintf(clusterListFile," %.6f ", distLevel);
		fprintf(clusterFile," %.6f ", distLevel);
		
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
		
		for (it=orphanNodes.begin(); it != orphanNodes.end(); ++it) {
			fprintf(clusterFile,"|%u",(*it));
			fprintf(clusterListFile, "1 ");
		}
		numOTUs += orphanNodes.size();
		
		fprintf(clusterFile,"|\n");
		fprintf(clusterListFile, "\n");		
		printf("Dist: %.6f. numOTUs: %u. numSingletons: %lu\n", distLevel, numOTUs, orphanNodes.size());
	}
	
	printf("\n");
	OTU.clear();    
	fclose(clusterListFile);
	fclose(clusterFile);
	return 1;
}

void updateTopParent(TreeNode* &aNode, TreeNode* &parentNode) {
	if (aNode == NULL)
		return;
	aNode->topParent = parentNode;
	updateTopParent(aNode->left, parentNode);
	updateTopParent(aNode->right, parentNode);
}

// idX > idY
void addNode(TreeNode* &nodeX, TreeNode* &nodeY, TreeNode* &parentNode, float newDist)
{   	
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
}	

// merge two clusters
void merge(TreeNode * nodeX, TreeNode * nodeY, float dist)
{   
	TreeNode *parentNode=NULL;	
	
	parentNode = new TreeNode(newId);	
	leaves[newId] = parentNode;	
	
	if (nodeX->ID > nodeY->ID)	
		addNode(nodeX, nodeY, parentNode, dist);	
	else
		addNode(nodeY, nodeX, parentNode, dist);	
	
	++ newId;	
}

int main(int argc, char* argv[])
{   	
	struct timeval startTime, endTime;
	
	gettimeofday(&startTime, NULL);
	
	string inFileName = "", cutoffFileName = "";
		
	printf("\n----------------------------------------------------------------------\n");
	printf("                  DENDROGRAM CUTTING GENETIC DISTANCES                  \n");

	int i, numReads;

	getOptions(argc, argv, inFileName, cutoffFileName, numReads);
	
	FILE * cutoffFile = NULL;
	vector<float> cutoffs;
	float cutoff;
	cutoffFile = fopen(cutoffFileName.c_str(),"r");
	while(!feof(cutoffFile)) {
		fscanf(cutoffFile, "%f\n", &cutoff);
		cutoffs.push_back(cutoff);
		cout << cutoff << endl;
	}
	fclose(cutoffFile);	

	leaves=(TreeNode**)malloc(sizeof(TreeNode*)*(2*numReads+1));

	if(leaves==0)
	{   
		cout<<"Error: Not enough memory" << endl;
		exit(-1);
	}
	
	for(i = 0; i < numReads; ++i) 
		leaves[i] = new TreeNode(i);	
	
	for(i = numReads+1; i < 2*numReads+1; ++i) 
		leaves[i] = 0;
	
	int idX, idY;
	float dist;		
	string mergeFileName=inFileName;
	mergeFileName.append("_Merge");	
	
	FILE * mergeFile;
	mergeFile = fopen(mergeFileName.c_str(), "r");
	
	newId = numReads;	
	while ( fscanf(mergeFile, "%d %d %f", &idX, &idY, &dist) == 3 ) 
		merge(leaves[idX-1], leaves[idY-1], dist);	
			
	fclose(mergeFile);
	cout << "DONE!" << endl;
	cout << "current node: " << newId << endl;	

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

 	// print output to files
 	
	string clusterListName, clusterName;

	clusterListName=inFileName;
	clusterListName.append(".Cluster_List");
	clusterName=inFileName;
	clusterName.append(".Cluster");

	printClusters(roots, orphanNodes, clusterListName, clusterName, cutoffs);		
	
	// clear memory
	emptyTree(roots);	
	roots.clear();	
	orphanNodes.clear();	
	free(leaves);
	
	gettimeofday(&endTime, NULL);	
	long elapsedTime = (endTime.tv_sec - startTime.tv_sec) * 1000u + (endTime.tv_usec - startTime.tv_usec) / 1.e3 + 0.5;

	printf("Time taken: %.3f s\n", elapsedTime/1.e3);
	printf("\n----------------------------------------------------------------------\n");
	return 0;
}

