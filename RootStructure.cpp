#include<stdexcept>
#include<algorithm>
#include"metis.h"
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
	for(int i=0;i!=rootNVert;++i){
		infile>>itemp
		      >>rootVerts[i].x
		      >>rootVerts[i].y
		      >>rootVerts[i].z;
		rootVerts[i].tag = i;// tag is the original order of the vertex
	}

	infile>>line;
	infile>>line;

	//read cells
	infile>>rootNElement;
	if(rootNElement<=0) throw logic_error("elements of the mesh file less than 0\n");

	rootElems = new InputElement[rootNElement];

	int ntag=0;	
	for(int i=0;i!=rootNElement;++i){
		infile>>itemp>>rootElems[i].type;
		infile>>ntag;
		for(int j=0;j!=ntag;++j)
			infile>>rootElems[i].tag[j];	
		infile>>ntag;
		nVertSum += nv; //record the sum of vertex in all element
		for(int j=0 ; j != numberOfNodesInElementTypeOf[rootElems[i].type]; ++j)
			infile>>rootElems[i].vertex[j];	//read all the vertex of the msh
	}

	infile.close();
	printf("complete reading in root\n");
}

/***************functor object used by partition routine*****************/
struct _SortAccordingToPart{
	InputElement* head;
	idx_t* partid;
	bool operator()(InputElement& lhs,InputElement& rhs){
		int i=&lhs - head;
		int j=&rhs - head;
		return partid[i]<partid[j] ? true : false ;
	}
};
struct _SortVertAccordingToPart{
	InputVert* head;
	idx_t* partid;
	bool operator()(InputVert& lhs,InputVert& rhs){
		int i=&lhs - head;
		int j=&rhs - head;
		return partid[i]<partid[j] ? true : false ;
	}
};
struct _AssigntoOriginIdx{
	int*  /*stop here, restart tomorrow*/
	int* originVertOrder;
	void operator()(InputVert& vert){
		oriv[vert.tag] = vert.tag;
	}
};
void RootProcess::partition(int N,DataPartition* dg){
	if(dg->comRank!=rank) return;//only in root

	printf("start partitioning in root \n");

	/*****************DATA PARTITION*******************
	 * 	data partition with METIS
	***************************************************/
	int numberOfNodesInElementTypeOf[8] = {
		0,2,3,4,4,8,6,5
	};

	int nv=0;
	for(int i=0;i!=rootNElement;++i)
		nv += numberOfNodesInElementTypeOf[rootElems[i].type];	

	rootgridList = new std::vector<int>(N,0);
	idx_t* eptr = new idx_t[rootNElement+1];
	idx_t* eind = new idx_t[nv];
	idx_t* epart = new idx_t[rootNElement]; //return value of the metis routine
	idx_t* npart = new idx_t[rootNVert];
	idx_t edgecut;

	eptr[0]=0;
	for(int i=0;i!=rootNElement;++i){
		nv = numberOfNodesInElementTypeOf[rootElems[i].type];
		for(int j=0;j!=nv;++j)
			eind[eptr[i]+j] =  rootElems[i].vertex[j]
		eptr[i+1]=eptr[i]+nv;
	}
	
	//call metis
	int ret = METIS_PartMeshDual(rootNElement,rootNVert,eptr,eind,
			NULL,NULL,		/*no weight edge*/
			3/*ncommon*/,	N/*npart*/,
			NULL,NULL,		/*now wighit partition*/
			&edgecut, epart, npart);

	if(ret!=METIS_OK)throw runtime_error("METIS Partition error!\n");
	printf("METIS return successfuly, partitioned into%d part\n",N);
	/*****************DATA PARTITION*******************
	 * 	reorder elements according to epart;
	***************************************************/

	for(int i=0;i!=N;++i)
		rootgridList->at(i) = std::count(epart,epart + rootNElement,i);

	_SortAccordingToPart op;
	op.head = rootElems;
	op.partid = epart;
	std::sort(rootElems,rootElems + rootNElement, op); //sort according to npart;

	

	/*****************DATA PARTITION*******************
	 * 	reorder vertex according to npart, update index in rootElems !!
	***************************************************/
	_SortVertAccordingToPart op2;
	op2.head = rootVerts;
	op.partid = npart;
	std::sort(rootVerts,rootVerts + rootNElement, op2);

	int* originVertOrder = new int[rootNVert];
	_AssigntoOriginIdx op3;	
	op3.originVertOrder = originVertOrder;
	std::for_each(rootVerts,rootVerts+rootNVert,op3);	
	delete originVertOrder; originVertOrder = NULL;



	delete []etpr,eind,epart,npart;
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



