#include<stdexcept>
#include"MPIStructure.h"
using namespace std;
void RootProcess::init(DataPartition* dg){
	 //only root may init
	if(dg->comRank!=rank) return;
	printer = new TerminalPrinter;
}


void RootProcess::allocate(DataPartition* dg){ //prepare for gather;
	if(dg->comRank!=rank) return; //only in root
	rootArrayBuffer = new double[rootNGlobal];
}


void RootProcess::clean(){
	delete rootElems;
	delete rootVerts;
	rootElems = NULL;
	rootVerts = NULL;
}

void RootProcess::read(DataPartition* dg,const string& title){
	if(dg->comRank!=rank) return; //only in root
	int itemp;
	string line;
	int numberOfNodesInElementTypeOf[8] = {
		0,2,3,4,4,8,6,5
	};
	printf("start reading in root\n");
	std::ifstream infile(title.c_str());
	
	//skip 6 lines
	if(!infile) throw logic_error("cant find grid file\n");
	infile>>line;
	infile>>line;
	infile>>line;
	infile>>line;
	infile>>line;
	infile>>line;

	infile>>rootNVert;
	if(rootNVert<=0) throw logic_error("vertex less than 0\n");

	/**************************************************
	 *	PHASE 1: 
	 **************************************************/
	//read verts
	rootVerts = new InputVert[rootNVert];
	for(int i=0;i!=rootNVert;++i)
		infile>>itemp
		      >>rootVerts[i].x
		      >>rootVerts[i].y
		      >>rootVerts[i].z;

	infile>>line;
	infile>>line;

	//read cells
	infile>>rootNElement;
	if(rootNElement<=0) throw logic_error("elements of the mesh file less than 0\n");
	rootElems = new InputElement[rootNElement];

	int ntag;	
	for(int i=0;i!=rootNElement;++i){
		infile>>itemp>>rootElems[i].type;
		infile>>ntag;
		for(int j=0;j!=ntag;++j)
			infile>>rootElems[i].tag[j];	
		infile>>ntag;
		for(int j=0 ; j != numberOfNodesInElementTypeOf[rootElems[i].type]; ++j)
			infile>>rootElems[i].vertex[j];	//read all the vertex of the msh
	}

	infile.close();
	printf("complete reading in root\n");
	printf("now the input array is \n");
}


void RootProcess::partition(int N,DataPartition* dg){
	if(dg->comRank!=rank) return;//only in root
	printf("start partitioning in root \n");
	/*****************DATA PARTITION*******************/
	rootgridList = new std::vector<int>(N,0);
	int n = rootNGlobal / N;
	int counter = 0;
	for(std::vector<int>::iterator it = rootgridList->begin(); it!=rootgridList->end()-1; ++it){
		*it =  n;
		counter+=n;
	}
	rootgridList->back() = rootNGlobal - counter;
	/**************************************************/
	printf("complete partitioning in root \n");
}


void RootProcess::write(DataPartition* dg){
	if(dg->comRank!=rank) return; //only in root
	printf("writing....");
	std::ofstream outfile("result.dat");
	char temp[256];
	for(int i=0;i!=rootNGlobal;++i){
		//sprintf(temp,"%15d\tx:%15e\tb:%15e\n",i+1,rootuBuffer[i],rootbuBuffer[i]);
		//outfile<<temp;
	}
	outfile.close();
}



