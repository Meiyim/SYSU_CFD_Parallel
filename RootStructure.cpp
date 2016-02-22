#include<stdexcept>
#include<algorithm>
#include"metis.h"
#include"MPIStructure.h"

#define NCOMMON 3
using namespace std;
void RootProcess::init(DataPartition* dg){
	 //only root may init
	if(dg->comRank!=rank) return;
	//printer = new TerminalPrinter;
	//other initiation stuff 
}


void RootProcess::allocate(DataPartition* dg){ //prepare for gather;
	if(dg->comRank!=rank) return; //only in root
	rootArrayBuffer = new double[rootNGlobal];
}


void RootProcess::clean(){
	delete rootElems;
	delete rootVerts;
	delete nodesPool;
	rootElems = NULL;
	rootVerts = NULL;
	nodesPool = NULL;
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
	}

	infile>>line;
	infile>>line;

	//read cells
	infile>>rootNElement;
	if(rootNElement<=0) throw logic_error("elements of the mesh file less than 0\n");

	rootElems = new InputElement[rootNElement];

	int ntag=0;	
	int _v=0;
	for(int i=0;i!=rootNElement;++i){
		infile>>itemp>>rootElems[i].type;
		infile>>ntag;
		for(int j=0;j!=ntag;++j)
			infile>>rootElems[i].tag[j];	
		for(int j=0 ; j != numberOfNodesInElementTypeOf[rootElems[i].type]; ++j){
			infile>>_v;
			rootElems[i].vertex[j] = _v-1;	//read all the vertex of the msh
		}
	}

	infile.close();
	printf("complete reading in root\n");
}

/***************functor object used by partition routine*****************/
struct _SortAccordingToPart{
	bool operator()(InputElement& lhs,InputElement& rhs){
		return lhs.pid<rhs.pid ? true : false;
	}
};
/*
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
	int*  //stop here, restart tomorrow
	int* originVertOrder;
	void operator()(InputVert& vert){
		originVertOrder[vert.tag] = vert.tag;
	}
};
*/
void RootProcess::partition(DataPartition* dg, int N){
	if(dg->comRank!=rank) return;//only in root

	printf("start partitioning in root \n");

	/*****************DATA PARTITION*******************
	 * 	phase2 :data partition with METIS
	***************************************************/
	int numberOfNodesInElementTypeOf[8] = {
		0,2,3,4,4,8,6,5
	};

	int nv=0;
	for(int i=0;i!=rootNElement;++i)
		nv += numberOfNodesInElementTypeOf[rootElems[i].type];	

	rootgridList = new std::vector<int>(N,0);
	idx_t* eptr  = new idx_t[rootNElement+1];
	idx_t* eind  = new idx_t[nv];
	idx_t* epart = new idx_t[rootNElement]; //return value of the metis routine
	idx_t* npart = new idx_t[rootNVert];
	idx_t edgecut;

	eptr[0]=0;
	for(int i=0;i!=rootNElement;++i){
		nv = numberOfNodesInElementTypeOf[rootElems[i].type];
		for(int j=0;j!=nv;++j)
			eind[eptr[i]+j] =  rootElems[i].vertex[j];
		eptr[i+1]=eptr[i]+nv;
	}
	
	//call metis
	idx_t _ncommon = NCOMMON;
	idx_t _ne = rootNElement;
	idx_t _nv = rootNVert;
	idx_t _np = N;
	int ret = METIS_PartMeshDual(&_ne,&_nv,eptr,eind,
			NULL,NULL,		/*no weight edge*/
			&_ncommon,		/*ncommon*/
			&_np,			/*npart*/
			NULL,NULL,		/*no weight partition*/
			&edgecut, epart, npart);

	if(ret!=METIS_OK)throw runtime_error("METIS Partition error!\n");
	printf("METIS return successfuly, partitioned into %d part, totally edgecut%lld \n",N,edgecut);
	//for(int i=0;i!=rootNVert;++i)
	//	PetscPrintf(MPI_COMM_WORLD,"!!!!!!!!!!!!!!!!npart: %d , %lld\n",i,npart[i]);
	delete []eptr;delete eind; eptr=NULL;eptr=NULL;

	for(int i=0;i!=rootNElement;++i)
		rootElems[i].pid = epart[i];

	/*****************DATA PARTITION*******************
	 * 	phase3 :reorder elements according to epart;
	***************************************************/

	for(int i=0;i!=N;++i)
		rootgridList->at(i) = std::count(epart,epart + rootNElement,i);
	delete []epart, epart = NULL;
	_SortAccordingToPart op;
	std::sort(rootElems,rootElems + rootNElement, op); //sort according to npart;


	//for(int i=0;i!=rootNElement;++i)
	//	PetscPrintf(MPI_COMM_WORLD,"rootElems: pid %d, type %d\n",rootElems[i].pid,rootElems[i].type);
	

	/*****************DATA PARTITION*******************
	 * 	phase4 :create nodes pool and interfaces info
	***************************************************/

	nodesPool = new map<int,unordered_set<int> >[N]; // one for each partition,<partID,nodesInthisParts>
	//interfacesInfo = new map<int,int>[N]; 			//one for each partition,<partID,width>
	int iE = 0;	
	for(int i=0;i!=N;++i){ //each part
		int nElementsPerPart = rootgridList->at(i);

		nodesPool[i][i].reserve(nElementsPerPart*4);
		for(int j=0;j!=nElementsPerPart;++j){ //each elements
			InputElement& _thisEle = rootElems[iE++];
			nv = numberOfNodesInElementTypeOf[_thisEle.type]; 
			for(int k=0;k!=nv;++k){	     //each vertex
				int _thisVert = _thisEle.vertex[k];
				int _partid = npart[_thisVert];

				nodesPool[i][_partid].insert(  _thisVert ); //if vert exists, nothing will do

				if(_partid != i){ //if is a boundary
					//interfacesInfo[i][_partid]++; //count the number of body & boundaries
					//interfacesInfo[_partid][i]++; // update the part info on the other side of the boundary	
					nodesPool[_partid][_partid].erase(_thisVert); //correct the nodesPool of the partition on the
					nodesPool[_partid][i].insert( _thisVert );    // other side of the boundary
				}
			}	
		}
	}
	
	delete []epart;delete npart;
	printf("complete partitioning in root \n");
	printf("report partition result:\n");
	for(int i=0;i!=N;++i){
		printf("part: %d, elements: %d, interfaces:%lu\n",i,rootgridList->at(i),nodesPool[i].size());
		for(auto it = nodesPool[i].begin();it!=nodesPool[i].end();++it){
			printf("\t pid: %d: %lu\n",it->first,it->second.size());	
		}
	}

}


void RootProcess::write(DataPartition* dg){
	if(dg->comRank!=rank) return; //only in root
	printf("writing....");
	std::ofstream outfile("result.dat");
	//char temp[256];
	for(int i=0;i!=rootNGlobal;++i){
		//sprintf(temp,"%15d\tx:%15e\tb:%15e\n",i+1,rootuBuffer[i],rootbuBuffer[i]);
		//outfile<<temp;
	}
	outfile.close();
}



