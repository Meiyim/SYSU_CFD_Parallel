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
	for(int i=0;i!=rootNElement;++i){
		delete rootElems[i];
	}
	delete []rootElems;
	delete []rootVerts;
	delete []nodesPool;
	delete []boundNodesPool;
	rootElems = NULL;
	rootVerts = NULL;
	nodesPool = NULL;
	boundNodesPool = NULL;
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

	rootElems = new InputElement* [rootNElement];

	int _type;
	int _ntag;
	int _v=0;
	rootNGlobal = 0;
	for(int i=0;i!=rootNElement;++i){
		infile>>itemp>>_type>>_ntag;
		rootElems[i] = new InputElement(_type,_ntag,numberOfNodesInElementTypeOf[_type]);

		for(int j=0;j!=_ntag;++j)
			infile>>rootElems[i]->tag[j];	
		for(int j=0 ; j != numberOfNodesInElementTypeOf[_type]; ++j){
			infile>>_v;
			rootElems[i]->vertex[j] = _v-1;	//read all the vertex of the msh
		}
		if(rootElems[i]->type > 3) ++rootNGlobal; 
	}

	infile.close();
	printf("complete reading in root\n");
}

/***************functor object used by partition routine*****************/
struct _SortAccordingToPart{
	bool operator()(InputElement* lhs,InputElement* rhs){
		return lhs->pid<rhs->pid ? true : false;
	}
};

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
		nv += numberOfNodesInElementTypeOf[rootElems[i]->type];	

	rootgridList = new std::vector<int>(N,0);
	rootNCells = new std::vector<int>(N,0);
	idx_t* eptr  = new idx_t[rootNElement+1];
	idx_t* eind  = new idx_t[nv];
	idx_t* epart = new idx_t[rootNElement]; //return value of the metis routine
	idx_t* npart = new idx_t[rootNVert];

	idx_t edgecut;

	eptr[0]=0;
	for(int i=0;i!=rootNElement;++i){
		nv = numberOfNodesInElementTypeOf[rootElems[i]->type];
		for(int j=0;j!=nv;++j)
			eind[eptr[i]+j] =  rootElems[i]->vertex[j];
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
	delete []npart;npart = NULL;
	for(int i=0;i!=rootNElement;++i)
		rootElems[i]->pid = epart[i];

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

	boundNodesPool = new map<int,unordered_set<int> > [N]; // one for each partition,<partID,nodesInthisParts>
	nodesPool = new map<int,int> [N];		      //one for each partition,<partID,width>
	//interfacesInfo = new map<int,int>[N]; 			
	int iE = 0;	
	for(int i=0;i!=N;++i){ //each part
		int nElementsPerPart = rootgridList->at(i);
	
		//nodesPool[i].reserve(nElementsPerPart*4);  // map does NOT have a reserve method!!!
		//boundNodesPool[i][i].reserve(nElementsPerPart);// a guess to the boundary nodes No. to be optimized

		for(int j=0;j!=nElementsPerPart;++j){ //each elements
			InputElement* _thisEle = rootElems[iE++];
			nv = numberOfNodesInElementTypeOf[_thisEle->type]; 
			if(_thisEle->type > 3){
				(rootNCells->at(i))++; 	//record the ncells for each part
			} 
			for(int k=0;k!=nv;++k){	     //each vertex
				int _thisVert = _thisEle->vertex[k];
				//int _partid = npart[_thisVert];

				nodesPool[i].insert(make_pair<int,int>(_thisVert,0) );//nodesPool[i][i] contains all verts of a partition, 
				
				/*
				if(_partid != i){ //if is a boundary, nodesPool[i][_partid] contains verts on each boundary

					boundNodesPool[i][_partid].insert(  _thisVert ); //if vert exists, nothing will do
					//interfacesInfo[i][_partid]++; //count the number of body & boundaries
					//interfacesInfo[_partid][i]++; // update the part info on the other side of the boundary	
					//boundNodesPool[_partid][_partid].erase(_thisVert); //correct the nodesPool of the partition on the
					boundNodesPool[_partid][i].insert( _thisVert );    // other side of the boundary
				}
				*/
			}	

		}
	}	


	//bug fixed! boundNodesPool should obtained by intersection
	for(int thisPart = 0;thisPart!=N;++thisPart){
		for(int thatPart = thisPart+1;thatPart!=N;++thatPart){
			vector<pair<int,int> > res;
			set_intersection(nodesPool[thisPart].begin(),nodesPool[thisPart].end(),
					 nodesPool[thatPart].begin(),nodesPool[thatPart].end(),
					 back_inserter(res));
			if(res.empty()){
				//no interface with this and that
			}else{
				for(vector<pair<int,int> >::iterator iter = res.begin();iter!=res.end();++iter){
					boundNodesPool[thisPart][thatPart].insert(iter->first);
					boundNodesPool[thatPart][thisPart].insert(iter->first);
				}
			}
		}
	}

	printf("complete partitioning in root \n");
	/*****************DATA PARTITION*******************
	 * 	phase5 : translate global idx of vertex list to local idx
	 * 		 translate global idx of interface pool to local idx
	 * 		 
	 * 	THIS STEP IS DEFERED TO SENDING ROUTINE
	***************************************************/
	
	printf("report partition result:\n");
	for(int i=0;i!=N;++i){
		printf("part: %d, elements: %d, interfaces:%lu\n",i,rootgridList->at(i),boundNodesPool[i].size());
		for(auto it = boundNodesPool[i].begin();it!=boundNodesPool[i].end();++it){
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


/*********************************************
 *	tecplot printer
 *	call after allocate
 *	call after data collection
 *********************************************/
void RootProcess::printTecplotHead(DataPartition* dg, ofstream& file){
	if(dg->comRank!=rank) return; //root only
	file<<"variables="<<"\"x\","<<"\"y\","<<"\"z\""
		<<"\"p\","<<"\"u\","<<"\"v\","<<"\"w\","<<"\"ro\","<<"\"T\","
		<<"\"Mach\","<<"\"mu\""<<endl;
}


/*********************************************
 * 	build MPI transfer buffer for vertex  element interfaceinfo
 * 	
 * 	WARNING: this 3 function MUST be called in specifig sequence!
 * 	calling sequence :
 * 		getvertex[pid]
 * 		getelement[pid]
 * 		getinterface[pid]
 *
 * 	WARNING: it is the user's responsibility to free buffer return bu these functions
 *
 *********************************************/


int RootProcess::getVertexSendBuffer(DataPartition*dg,int pid, double** buffer){
	if(dg->comRank!=rank) return 0; //root only

	/******************************************
	 *	get the local vertex id for part pid
	 ******************************************/



	// build buffer
	size_t nV = nodesPool[pid].size();
	*buffer = new double[nV*3];

	size_t counter = 0;
	size_t localID = 0;
	for(map<int,int>::iterator it = nodesPool[pid].begin(); it!=nodesPool[pid].end(); ++it){
		it->second = localID++; // the second is local id
		(*buffer)[counter++] = rootVerts[ it->first ].x;
		(*buffer)[counter++] = rootVerts[ it->first ].y;
		(*buffer)[counter++] = rootVerts[ it->first ].z;
	}

	//printf("vertex buffer for %d prepared\n",pid);
	return counter;
}


int RootProcess::getElementSendBuffer(DataPartition* dg,int pid ,int** buffer){
	if(dg->comRank!=rank) return 0; //root only

	int numberOfNodesInElementTypeOf[8] = {
		0,2,3,4,4,8,6,5
	};

	
	/**************************************************
	 *	transfer element Vertex list to local idx
	 **************************************************/
	size_t iEhead = 0;
	for(int i=0;i!=pid;++i){
		iEhead+=rootgridList->at(i);
	}
	size_t iE = iEhead;
	size_t counter = 0; //calculate the buffer length
	for(int j=0;j!=rootgridList->at(pid);++j){
		InputElement* _thisEle = rootElems[iE++];
		int nv = numberOfNodesInElementTypeOf[_thisEle->type];
		for(int k=0;k!=nv;++k){
			int globalID = _thisEle->vertex[k];
			int localID = nodesPool[pid][globalID];
			_thisEle->vertex[k] = localID;
		}
		counter+= (2 + _thisEle->ntag + nv);
	}

	// build buffer
	(*buffer) = new int [ counter ];//

	counter = 0;
	iE = iEhead;
	for(int j=0;j!=rootgridList->at(pid);++j){ //each elements
		InputElement* _thisEle = rootElems[iE++];
		(*buffer)[counter++] = _thisEle->type;
		(*buffer)[counter++] = _thisEle->ntag;
		for(int k=0;k!=_thisEle->ntag;++k){
			(*buffer)[counter++] = _thisEle->tag[k];
		}	
		for(int k=0;k!=numberOfNodesInElementTypeOf[_thisEle->type];++k){
			(*buffer)[counter++] = _thisEle->vertex[k];
		}
	}

	//printf("element buffer for %d prepared\n",pid);
	return counter;
}





int RootProcess::getInterfaceSendBuffer(DataPartition* dg,int pid ,int** buffer){
	if(dg->comRank!=rank) return 0; //root only


	// build buffer
	//
	size_t nI = 0;	

	nI++; // head contains the number of interfaces;
	//the data structure of this buffer:
	//
	//numberOfInterfaces, interfacePartID, interfaceWidth, vert1, ver2, ver3,... interfacePartID2, interfaceWidth2, ...
	
	for(map<int,unordered_set<int> >::iterator ittup = boundNodesPool[pid].begin(); ittup!=boundNodesPool[pid].end(); ++ittup){
		nI += 2;			 //each bounds has a head to record the interface info;
		nI += ittup->second.size();
	}

	*buffer = new int[nI];

	size_t counter = 0;	
	(*buffer)[counter++] = boundNodesPool[pid].size();

	for(map<int,unordered_set<int> >::iterator ittup = boundNodesPool[pid].begin(); ittup!=boundNodesPool[pid].end(); ++ittup){
		(*buffer)[counter++] = ittup->first;
		(*buffer)[counter++] = ittup->second.size();	 
		for(unordered_set<int>::iterator itnodes = ittup->second.begin(); itnodes!=ittup->second.end(); ++itnodes){
			int globalID = *itnodes;
			int localID = nodesPool[pid][globalID];
			(*buffer)[counter++] = localID;  //transfer boundNodesPool to local id
		}
	}

	//printf("interface buffer for %d prepared\n",pid);
	return nI;

}


/***************************************************
 * 	 print screen
 * 	 root only
 * *************************************************/
void RootProcess::printStarter(DataPartition* dg){
	if(dg->comRank != rank) return ;
	system("clear");
	printf("    1111111    111        111    1111111         111111         1111111          \n");
	printf("  11111111111  111        111  11111111111      11111111      11111111111        \n");
	printf(" 1111     1111 111        111 1111     1111    1111  1111    1111     1111       \n");
	printf(" 1111            111    111   1111            1111    1111   1111     1111       \n");
	printf(" 1111              111111     1111           1111      1111    111               \n");
	printf(" 1111               1111      1111           1111      1111      11111           \n");
	printf(" 1111               1111      1111           11111111111111         1111         \n");
	printf(" 1111               1111      1111           1111      1111  1111     1111       \n");
	printf(" 1111     1111      1111       1111    1111  1111      1111  1111     1111       \n");
	printf("  11111111111       1111        1111111111   1111      1111   11111111111        \n");
   	printf("    1111111         1111          111111     1111      1111     1111111          \n");
	printf("\n");
	printf("\n");
    	printf("000000000000000000000000000000000000000000000000000000000000000000000000000000\n");
	printf("0 0000 0 0    00     0      0     00 0000 00     0      0 0000 0     00      0\n");
	printf("0 0000 0  000 0 00000000  000 0000 0 0000 0 00000000  000 0000 0 0000 0 000000\n");
	printf("0 0000 0 0000 00    0000  000     00 0000 0 00000000  000 0000 0     00      0\n");
	printf("0 0000 0 0000 000000 000  000 000 00 0000 0 00000000  000 0000 0 000 00 000000\n");
	printf("00    00 0000 0     0000  000 0000 00    000     000  0000    00 0000 0      0\n");
    	printf("000000000000000000000000000000000000000000000000000000000000000000000000000000\n");


	printf("                   *******************************************\n");
	printf("                     COMPUTATIONAL CODES FOR FLUID DYNAMICS\n");
	printf("                   *******************************************\n");
	printf("                   VERSION 2.0                     27 NOV 2015 \n");
	printf("                             EXECUTABLE ATTRIBUTES \n\n");
	printf("   All rights reserved. Unauthorized use, distribution or duplication is \n");
	printf("   prohibited. This product is subject to China Nuclear Power Technology \n");
	printf("   Research Institute and State Nuclear Power Software Development Center.\n");
	printf("--------------------------------------------------------------------------------\n\n");
}

void RootProcess::printEnding(DataPartition* dg){
	if(dg->comRank!=rank) return;
	printf( " -------------------------------------------------------------------------------- \n\n");
	printf( "                      CPU REQUIREMENTS OF NUMERICAL SOLUTION\n\n");
	printf( " -------------------------------------------------------------------------------- \n");
	printf( "    The total problem used   : %10.1f seconds of CPU time\n", 0.0);
	printf( "    includes system time     : %10.1f seconds of CPU time\n", 0.0);
	printf(	"    for an average of        : %10.1f  MICROSECONDS/CELL/CYCLE.\n", 0.0);
	printf( "    Total wall clock time of : %10.1f  seconds \n", 0.0);
}

void RootProcess::printStepStatus(DataPartition*dg, int step,int piter ,double time,double dt,double res){
	if(dg->comRank!=rank) return;
	printf("%15f\t%10d\t%10d\t%13.5f\t%13.5f\n",time,step,piter,dt,res);
}
void RootProcess::printSteadyStatus(DataPartition*dg,int step,double res){
	if(dg->comRank!=rank) return;
	printf("%15s\t%10d\t%10s\t%13s\t%15f\n","---",step,"---","---",res);
}
void RootProcess::printSectionHead(DataPartition* dg,double timeElapse){;
	if(dg->comRank!=rank) return;
	printf("%15s\t%10s\t%10s\t%13s\t%15s\n","TIME","CAL STEPE","ITER","DELT","MAX RES");
}
