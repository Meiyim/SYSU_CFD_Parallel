#include<stdexcept>
#include<algorithm>
#include"metis.h"
#include"tools.h"
#include"MPIStructure.h"
extern "C"{
#include "cgnslib.h"
}

#define NCOMMON 3
#define CRITICAL 1.e-10

#define TOT_FILE_TITLE "totResult.tot"
#define MONITOR_FILE_TITLE "monitorResult.mon"

using namespace std;

int numberOfNodesInElementTypeOf[8] = {//command
	0,2,3,4,4,8,6,5
};

void RootProcess::init(DataPartition* dg){
	 //only root may init
	if(dg->comRank!=rank) return;
	totFile.open(TOT_FILE_TITLE);
	monitorFile.open(MONITOR_FILE_TITLE);
	//printer = new TerminalPrinter;
	//other initiation stuff 
}


void RootProcess::allocate(DataPartition* dg){ //prepare for gather;
	if(dg->comRank!=rank) return; //only in root
	rootArrayBuffer = new double[rootNGlobal];
}


void RootProcess::clean(){
	
	if(rootElems!=NULL)
		for(size_t i=0;i!=rootNElement;++i){
			delete rootElems[i];
	}
	
	delete []rootElems;
	delete []rootVerts;

	rootElems = NULL; //put NULL avoid dangling pointer. a good habit
	rootVerts = NULL;
}


/*******************************************************/
// USING GMSH Notation
// read a Compressed Gmesh file
// Output:
// rootNGlobal
// rootNElement
// rootElems
/*******************************************************/
void RootProcess::readBin(DataPartition* dg, const string& title){
	if(dg->comRank!=rank) return; //only in root
	int itemp;
	string line;

	printf("start reading in root...");
	fflush(stdout);
	std::ifstream infile(title.c_str(),std::ios::binary);
	
	//skip 6 lines
	if(!infile) throw logic_error("cant find grid file\n");

	infile.read((char*)&rootNVert,sizeof(size_t));
	if(rootNVert<=0) throw logic_error("vertex less than 0\n");

	/**************************************************
	 *	PHASE 1: 
	 **************************************************/
	//read verts
	rootVerts = new InputVert[rootNVert];

	for(size_t i=0;i!=rootNVert;++i){
		infile.read((char*)&itemp,sizeof(int));
		infile.read((char*)&(rootVerts[i].x),sizeof(double));
		infile.read((char*)&(rootVerts[i].y),sizeof(double));
		infile.read((char*)&(rootVerts[i].z),sizeof(double));
	}
	//read cells

	infile.read((char*)&rootNElement,sizeof(size_t));
	if(rootNElement<=0) throw logic_error("elements of the mesh file less than 0\n");
	rootElems = new InputElement* [rootNElement];

	int _type;
	int _ntag;
	int _v=0;
	rootNGlobal = 0;
	for(size_t i=0;i!=rootNElement;++i){
		infile.read((char*)&itemp,sizeof(int));
		infile.read((char*)&_type,sizeof(int));
		infile.read((char*)&_ntag,sizeof(int));

		rootElems[i] = new InputElement(_type,_ntag,numberOfNodesInElementTypeOf[_type]);

		for(int j=0;j!=_ntag;++j)
			infile.read((char*)&(rootElems[i]->tag[j]),sizeof(int));
		for(int j=0 ; j != numberOfNodesInElementTypeOf[_type]; ++j){
			infile.read((char*)&_v, sizeof(int));
			rootElems[i]->vertex[j] = _v-1;	//read all the vertex of the msh
		}
		if(rootElems[i]->type > 3) ++rootNGlobal; 
	}

	infile.close();
	printf("done\n");
}



/*******************************************************/
// USING CGNS Notation
// Output:
// rootNGlobal
// rootNElement, rootNVerts
// rootElems,rootVerts
/*******************************************************/
#define CGNS_VERBOSE
void RootProcess::readCGNS(DataPartition* dg, const string& title){
	if(dg->comRank != rank) return;	
	printf("start reading CGNS in root");
	fflush(stdout);
    int icg;
    // open files
	float version;
    if( cg_open(title.c_str(), CG_MODE_READ, &icg) )
    	errorHandler.fatalLogicError("CGNS error",cg_get_error());
    // check CGNS Library Version used for file creation:
    if( cg_version(icg,&version) )
    	errorHandler.fatalLogicError("CGNS error",cg_get_error());

#ifdef CGNS_VERBOSE
    cout<<"\nLibrary Version used for file creation:"<<version<<endl;
#endif
    
    int  i,j, ier;
    // read zone and grid types
    int        ibase = 1,
               zone  = 1;
    char       zonename[200];
    ZoneType_t zonetype;
    cgsize_t        size;
    if(  cg_zone_read(icg, ibase, zone, zonename, &size) )//this size is node size
    	errorHandler.fatalLogicError("CGNS error",cg_get_error());

#ifdef CGNS_VERBOSE
    cout<<"Name of Zone"<<zone<<" is"<<zonename<<endl;
    printf("zone %d name: %s, size: %d\n",zone,zonename,size);
#endif
    zone=1;
    if(cg_zone_type(icg, ibase, zone, &zonetype))
    	errorHandler.fatalLogicError("CGNS error",cg_get_error());
#ifdef CGNS_VERBOSE
    cout<<"  Zone type is "<< cg_ZoneTypeName(zonetype)<<endl;
#endif
    if(      zonetype == Structured )
    {
    	errorHandler.fatalLogicError("structured grids. stop");
    }
    else if( zonetype == Unstructured ){
    	//pass
    }
    else
    {
    	errorHandler.fatalLogicError("unknown zonetype.");
    }
    
    // grid coordinates
    cgsize_t range_min=1,range_max=size;
    int ncoords, NDim,
        nnode     = size;
    DataType_t  datatype;
    char        coordname[200];
    if( cg_ncoords(icg,ibase,zone, &ncoords) )
    	errorHandler.fatalLogicError("CGNS error",cg_get_error());

#ifdef CGNS_VERBOSE
    cout<< "ncoords: "<< ncoords<<endl;
#endif
    NDim   = ncoords;
    assert(NDim==3);
   	rootNVert = size;
   	rootVerts = new InputVert[rootNVert];
    for( int idim=1; idim<=NDim; idim++ )
    {
        if(  cg_coord_info(icg,ibase,zone, idim,  &datatype, coordname) )
	    	errorHandler.fatalLogicError("CGNS error",cg_get_error());
        if( std::string(DataTypeName[datatype]) != "RealDouble" ){
        	errorHandler.fatalLogicError("error, change coordinate datatype to single/Double");
        }
        double* temp = new double[size];
        if(  cg_coord_read( icg,ibase,zone, coordname, datatype, 
                             &range_min,&range_max, temp ) )
	    	errorHandler.fatalLogicError("CGNS error",cg_get_error());
        if(idim==1){
	        for(int i=0;i!=size;++i)
	            rootVerts[i].x = temp[i];
        }if(idim==2){
	        for(int i=0;i!=size;++i)
	            rootVerts[i].y = temp[i];
        }if(idim==3){
	        for(int i=0;i!=size;++i)
	            rootVerts[i].z = temp[i];
        }
        delete [] temp;
#ifdef CGNS_VERBOSE
        printf("%s: type: %s %d --> %d\n",coordname,DataTypeName[datatype],range_min,range_max);
#endif
    }

    int ncell =0;
    // cell connectivity and boundary
    int       nelem =0, nelem0,count, nsection, parent_flag;
    cgsize_t  ElementDataSize, *elements;
    char      elementSectionName[100];
    cgsize_t  start,end;
    int nbndry,npe;
    ElementType_t   type;

    if(cg_nsections(icg,ibase,zone,&nsection) )
    	errorHandler.fatalLogicError("CGNS error",cg_get_error());
#ifdef CGNS_VERBOSE
    printf("section number %d\n",nsection);
    cout<<"read Element Connectivity:\n";    
#endif
    rootNElement = 0;
    rootNGlobal = 0;
    for(int section = 1;section<=nsection;++section){
        ier = cg_section_read(icg,ibase,zone,section, elementSectionName, &type,
                 &start,&end,&nbndry,&parent_flag);
        rootNElement+=(end-start)+1;
        if(type>9) //volumn section
        	rootNGlobal+=(end-start)+1;
    }
    rootElems = new  InputElement* [rootNElement];

//		rootElems[i] = new InputElement(_type,_ntag,numberOfNodesInElementTypeOf[_type]);

    size_t iE = 0;
    for( int section=1; section<=nsection; section++ )
    {
        ier = cg_section_read(icg,ibase,zone,section, elementSectionName, &type,
                 &start,&end,&nbndry,&parent_flag);
 #ifdef CGNS_VERBOSE
        cout<< elementSectionName << ", type:" << cg_ElementTypeName(type) << endl;
 #endif
        ier = cg_ElementDataSize(icg,ibase,zone,section, &ElementDataSize);
        elements = new cgsize_t[ElementDataSize];
        ier = cg_elements_read(icg,ibase,zone,section, elements, NULL );
        int nelem_sec= end-start+1 ;
        // single-type grids
        if( type != MIXED )
        {
            // pure elements
            ier = cg_npe(type, &npe);
            // store in cell elements/boundary
            assert(NDim==3);
            if( NDim==3 && type>= 5)
            {
                for( i=0; i< nelem_sec; i++ )
				{
					int _type = -1;
					switch (type){
					case TRI_3:
						_type = 2;
						break;
					case QUAD_4:
						_type = 3;
						break;
                    case  TETRA_4:
						_type = 4;
                        break;
                    case  HEXA_8:
						_type = 5;
                        break;
                    default:
                    	errorHandler.fatalLogicError("unknown element type",cg_ElementTypeName(type));
                    }

                   	rootElems[i+iE] = new InputElement(_type,2,npe);
                    for( j=0;j!=npe; j++ )
                    	rootElems[i+iE]->vertex[j] = elements[i*npe+j] - 1;//why -1 ???
                    for(map<int,BdRegion>::const_iterator iter = regionMap->begin();iter!=regionMap->end();++iter){
                    	if(iter->second.name == string(elementSectionName) ){
		                    rootElems[i+iE]->tag[0] = iter->first;
		                    break;
                    	}else if(iter==regionMap->end()){
                    		assert(false);
                    	}
                    }
 #ifdef CGNS_VERBOSE 
                    //printf("found rid %d\n",rootElems[i+iE]->tag[0]);
 #endif                    
				}
                iE += nelem_sec;
            }else
                cout<< "useless data?";
        }
        // mixed-type grids
        else
        {
            // mixed element. Assume it only belongs to cells? not boundary
            //to be implement
            assert(false);
            /*
            nelem0= nelem ;
            nelem = nelem + nelem_sec ;
            count = 0 ;
            for( i=1; i<=nelem_sec; i++ )
            {
                count= count + 1 ;
                type = (ElementType_t)elements[count] ;
                ier  = cg_npe( type, &npe );
                //KNode[i]= npe ;
				for( j=1;j<=npe;j++ )
					ENode[i][j] = elements[count+j];  // (count+1:count+npe);
                count = count+npe ;
            }
            */
        }
        delete [] elements;
    }
#ifdef CGNS_VERBOSE
    printf("nelement %d, nvert %d, nglobal %d\n",rootNElement,rootNVert,rootNGlobal);
#endif
    cg_close(icg);
    cout<<"done"<<endl;
}




/*******************************************************/
// USING GMSH Notation
// Output:
// rootNGlobal
// rootNElement
// rootElems
/*******************************************************/
void RootProcess::read(DataPartition* dg,const string& title){
	if(dg->comRank!=rank) return; //only in root
	int itemp;
	string line;

	printf("start reading in root...");
	fflush(stdout);
	std::ifstream infile(title.c_str());
	
	//skip 6 lines
	if(!infile) throw logic_error("cant find grid file\n");

	do{
		infile>>line;
	}while(line!="$Nodes");

	infile>>rootNVert;
	if(rootNVert<=0) throw logic_error("vertex less than 0\n");

	/**************************************************
	 *	PHASE 1: 
	 **************************************************/
	//read verts
	rootVerts = new InputVert[rootNVert];
	for(size_t i=0;i!=rootNVert;++i){
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
	for(size_t i=0;i!=rootNElement;++i){
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
	printf("done\n");
}


//temporary class only used in this file;
struct PeriodicBcElement
{
	int originalIdx;
	double n[3];
	double x[3];
};
void getNormalVector(double** verts, double* normalVec){
	double dx1[3],dx2[3];
	vec_minus( dx1, verts[0],verts[2], 3 ); // dx1= vert[n1] - vert[n3]
    vec_minus( dx2, verts[1],verts[3], 3 ); // dx2= vert[n2] - vert[n4]
    vec_cross( normalVec, dx1, dx2 );     // face.n is normal area
}
int getAdjacent(idx_t* xadj, idx_t* adjncy,int ib){
	int head = xadj[ib];
	int end = xadj[ib+1];
	assert(end == head+1);
	return adjncy[head];
}
void RootProcess::buildPeriodicInterface(InputElement** elements,idx_t* xadj, idx_t* adjncy,int thisbid, int thatbid){
	//check period bc pair:
	size_t nThisBnd = 0;
	size_t nThatBnd = 0;
	for(int i=0;i!=rootNElement;++i){
		int bid = elements[i]->tag[0];
		if(bid == thisbid) nThisBnd++;
		if(bid == thatbid) nThatBnd++;
	}
	if(nThatBnd==0 || nThisBnd == 0){
		errorHandler.fatalLogicError("cant find periodic bc, bid",thisbid);
	}
	if(nThisBnd != nThatBnd){
		errorHandler.fatalLogicError("periodic bc has differend boundary cells,bid",thisbid);
	}

	PeriodicBcElement* thisBnds = new PeriodicBcElement[nThisBnd];
	PeriodicBcElement* thatBnds = new PeriodicBcElement[nThatBnd];

	size_t iThisbnd = 0;
	size_t iThatbnd = 0;
	double thisCenter[3] = {0.,0.,0.};
	double thatCenter[3] = {0.,0.,0.};

	for(int i=0;i!=rootNElement;++i){
		int bid = elements[i]->tag[0];
		if(bid == thisbid || bid == thatbid){
			int elemType = elements[i]->type;
			int nv = numberOfNodesInElementTypeOf[elemType];

			double** verts = new_Array2D<double>(4,3);

			for(int j = 0;j!=nv;++j){
				verts[j][0] = rootVerts[ elements[i]->vertex[j] ].x;
				verts[j][1] = rootVerts[ elements[i]->vertex[j] ].y;
				verts[j][2] = rootVerts[ elements[i]->vertex[j] ].z;
			}

			if(nv==3){
				verts[3][0] = verts[2][0];
				verts[3][1] = verts[2][1];
				verts[3][2] = verts[2][2];
			}

			PeriodicBcElement toInsert;
			toInsert.originalIdx = i;
			toInsert.x[0] = toInsert.x[1] = toInsert.x[2] = 0.0;
			for(int j=0;j!=nv;++j){
				toInsert.x[0] += verts[j][0];
				toInsert.x[1] += verts[j][1];
				toInsert.x[2] += verts[j][2];
			}

			getNormalVector(verts,toInsert.n);
			toInsert.x[0] /= (double)nv;
			toInsert.x[1] /= (double)nv;
			toInsert.x[2] /= (double)nv;
			delete_Array2D(verts,4,3);

			if(bid == thisbid) {
				thisBnds[iThisbnd++] = toInsert;//default copy function
				thisCenter[0] += toInsert.x[0];
				thisCenter[1] += toInsert.x[1];
				thisCenter[2] += toInsert.x[2];
			}
			if(bid == thatbid){
				thatBnds[iThatbnd++] = toInsert;
				thatCenter[0] += toInsert.x[0];
				thatCenter[1] += toInsert.x[1];
				thatCenter[2] += toInsert.x[2];
			}
		}
	}
	assert(iThisbnd==nThisBnd);
	assert(iThatbnd==nThatBnd);

	thisCenter[0] /= nThisBnd;
	thisCenter[1] /= nThisBnd;
	thisCenter[2] /= nThisBnd;

	thatCenter[0] /= nThatBnd;
	thatCenter[1] /= nThatBnd;
	thatCenter[2] /= nThatBnd;

	double connectionVec[3];
	vec_minus(connectionVec,thatCenter,thisCenter,3);
	assert( (*regionMap)[thisbid].type1==6 );
	(*regionMap)[thisbid].initvalues[0] = connectionVec[0];
	(*regionMap)[thisbid].initvalues[1] = connectionVec[1];
	(*regionMap)[thisbid].initvalues[2] = connectionVec[2];


	for(int i=0;i!=nThisBnd;++i){ //move to origin point
		for(int j=0;j!=3;++j){
			thisBnds[i].x[j] -= thisCenter[j];
			thatBnds[i].x[j] -= thatCenter[j];
		}
	}



	double* thisNormal = thisBnds[0].n;// bnd should in a plain surface !
	double* thatNormal = thatBnds[0].n;
	double cross[3];
	vec_cross(cross,thisNormal,thatNormal);
	if(vec_len(cross,3) < CRITICAL){
		//translation periodic
		set<int> foundSet;
		for(int i=0;i!=nThisBnd;++i){
			for(int j=0;j!=nThatBnd;++j){
				if(foundSet.find(j)!=foundSet.end()){
					continue;
				}
				double distance[3];
				vec_minus(distance,thatBnds[j].x,thisBnds[i].x,3);
				if( vec_len(distance,3) < CRITICAL){
					int iThisbnd = thisBnds[i].originalIdx;
					int iThatbnd = thatBnds[j].originalIdx;
					int thiscell = getAdjacent(xadj,adjncy,iThisbnd);
					int thatcell = getAdjacent(xadj,adjncy,iThatbnd);
					int thatPart = elements[ thatcell ]->pid;

					vector<int> _res(1,thatPart);
					_res.push_back(elements[ thatcell ]->idx);

					set<int> _thisSet;
					set<int> _thatSet;
					for(int k=0;k!=numberOfNodesInElementTypeOf[elements[thiscell]->type];++k){
						_thisSet.insert(elements[thiscell]->vertex[k]);
					}
					for(int k=0;k!=numberOfNodesInElementTypeOf[elements[iThisbnd]->type];++k){//that set is this bnd...
						_thatSet.insert(elements[iThisbnd]->vertex[k]);
					}
					set_intersection(_thisSet.begin(),_thisSet.end(),_thatSet.begin(),_thatSet.end(),back_inserter(_res));
					assert(_res.size()>=5 && _res.size() <= 6);

					elements[ thiscell ]->interfaceInfo.push_back(_res);
					foundSet.insert(j);
					break;
				}
			}
		}
		assert(foundSet.size()==nThatBnd);
	}else{
		//rotation periodic
		errorHandler.fatalLogicError("cannot match periodic boundary",thisbid);
	}

	delete [] thisBnds;
	delete [] thatBnds;

}


/***************functor object used by partition routine*****************/
struct _SortAccordingToPart{
	bool operator()(InputElement* lhs,InputElement* rhs){
		return lhs->pid<rhs->pid ? true : false;
	}
};

void RootProcess::partition(DataPartition* dg, int N){
	if(dg->comRank!=rank) return;//only in root

	printf("start partitioning in root...");
	fflush(stdout);

	/*****************DATA PARTITION*******************
	 * 	phase2 :data partition with METIS
	***************************************************/


	int nv=0;
	for(size_t i=0;i!=rootNElement;++i)
		nv += numberOfNodesInElementTypeOf[rootElems[i]->type];	

	rootgridList.resize(N);
	rootNCells.resize(N); 
	idx_t* eptr  = new idx_t[rootNElement+1];
	idx_t* eind  = new idx_t[nv];
	idx_t* epart = new idx_t[rootNElement]; //return value of the metis routine

	idx_t edgecut;

	eptr[0]=0;
	for(size_t i=0;i!=rootNElement;++i){
		nv = numberOfNodesInElementTypeOf[rootElems[i]->type];
		for(int j=0;j!=nv;++j)
			eind[eptr[i]+j] =  rootElems[i]->vertex[j];
		eptr[i+1]=eptr[i]+nv;
	}
	
	//call metis
	idx_t _ncommon = NCOMMON;
	idx_t _numflg = 0;
	idx_t _ne = rootNElement;
	idx_t _nv = rootNVert;
	idx_t _np = N;
	idx_t *xadj = NULL;
	idx_t *adjncy = NULL;

	int ret = METIS_MeshToDual(&_ne,&_nv,eptr,eind,
			&_ncommon,
			&_numflg,
			&xadj,&adjncy);		//METIS Mesh to Graph

	if(ret!=METIS_OK) errorHandler.fatalRuntimeError("METIS Partition error!\n");

	delete []eptr;delete [] eind; eptr=NULL;eptr=NULL;

	idx_t _ncon = 1; // number of balancing constraints
	//idx_t options[METIS_NOPTIONS];
	//options[METIS_OPTION_NUMBERING] = 0;//C-stype numbering

	ret = METIS_PartGraphKway(&_ne,
			&_ncon,
			xadj,adjncy,
			NULL,//vertex weight
			NULL,//vertex size
			NULL,//edgh weight
			&_np,//nparts
			NULL,//weight for each partition/ constraints
			NULL,//weigh tolerece
			NULL,//other options
			&edgecut,
			epart);


	if(ret!=METIS_OK) errorHandler.fatalRuntimeError("METIS Partition error!\n");
	printf("METIS partition successfuly, partitioned into %d part, totally edgecut%ld \n",N,edgecut);
	
		
	for(size_t i=0;i!=rootNElement;++i){
		rootElems[i]->pid = epart[i];
	}
	for(int i=0;i!=N;++i)
		rootgridList.at(i) = std::count(epart,epart + rootNElement,i);
	delete []epart, epart = NULL;


	/*****************DATA PARTITION*******************
	 * 	phase3 :reorder elements according to epart;
	***************************************************/

	InputElement** elementsInOriginalOrdering = new InputElement* [rootNElement];
	for(size_t i=0;i!=rootNElement;++i){
		elementsInOriginalOrdering[i] = rootElems[i];
	}

	_SortAccordingToPart op;
	std::sort(rootElems,rootElems + rootNElement, op); //sort according to npart;

	int iE = 0;	
	int iEhead = 0;
	for(int i=0;i!=N;++i){ //each part
		int localCellIdxCounter = 0;
		int nElementsPerPart = rootgridList.at(i);
		for(int j=0;j!=nElementsPerPart;++j){ 	//each elements
			if(rootElems[iE]->type > 3){    //is a volume
				rootElems[iE]->idx = localCellIdxCounter + iEhead;// local idx of matrx + partition Segmentor(remove after sort in _connectionMap)
				localCellIdxCounter++;
			}else{
				rootElems[iE]->idx = -1; //boundary does not have a idx
			}
			iE++;
		}
		rootNCells[i] = localCellIdxCounter; //record Ncel;
		iEhead += rootgridList[i];  
	}	


	/*****************DATA PARTITION*******************
	 * 	phase4 :add interface info to interfaceCell;
	***************************************************/
	map<int,int>* _boundCells = new map<int,int> [N]; // one for each partition,<partID,interfaceWidth>

	for(size_t i=0;i!=rootNElement;++i){//original order
		InputElement* _thisEle = elementsInOriginalOrdering[i];
		int head = xadj[i];
		int end = xadj[i+1];
		int thispart = _thisEle->pid;
		for(int j=head;j!=end;++j){ //for each neighboring
			InputElement* _thatEle = elementsInOriginalOrdering[adjncy[j]];
			int thatpart = _thatEle->pid;
		
			if(thatpart != thispart){
				set<int> _thisSet,_thatSet;
				int thatIdx = _thatEle->idx;
				vector<int> _res(1,thatpart);
				_res.push_back(thatIdx);
				
				for(int k=0;k!=numberOfNodesInElementTypeOf[_thisEle->type];++k){
					_thisSet.insert(_thisEle->vertex[k]);
				}
				for(int k=0;k!=numberOfNodesInElementTypeOf[_thatEle->type];++k){
					_thatSet.insert(_thatEle->vertex[k]);
				}
				set_intersection(_thisSet.begin(),_thisSet.end(),_thatSet.begin(),_thatSet.end(),back_inserter(_res));
				assert(_res.size()>=5 && _res.size() <= 6);
					
				elementsInOriginalOrdering[i]->interfaceInfo.push_back(_res);
				
				_boundCells[thispart][thatpart]++;
			}
		}
	}


	//build periodic boundary connectivity;
	for(map<int,BdRegion>::iterator iter = regionMap->begin();iter!=regionMap->end();++iter){
		if(iter->second.type1 == 6 ){
			int thatbid = iter->second.type2;
			if(regionMap->find(thatbid)==regionMap->end()){
				errorHandler.fatalLogicError("cant find corresponding boundary for period bc ",iter->first);
			}
			if((*regionMap)[thatbid].type1 !=6 || (*regionMap)[thatbid].type2 != iter->first){
				errorHandler.fatalLogicError("periodic bc not correspond, bid ",iter->first);
			}
			
			buildPeriodicInterface(elementsInOriginalOrdering,xadj,adjncy,iter->first,thatbid);
		}
	}

	delete []elementsInOriginalOrdering;elementsInOriginalOrdering=NULL;

	METIS_Free(xadj);	xadj = NULL;
	METIS_Free(adjncy);	adjncy = NULL;

	
	printf("done\n");
	/*****************DATA PARTITION*******************
	 * 	phase5 : translate global idx of vertex list to local idx
	 * 		 
	 * 	THIS STEP IS DEFERED TO SENDING ROUTINE
	***************************************************/
	
	printf("report partition result:\n");
	for(int i=0;i!=N;++i){
		printf("part: %d, elements: %d, interfaces:%lu\n",i,rootgridList.at(i),_boundCells[i].size());
		for(map<int,int>::iterator it = _boundCells[i].begin();it!=_boundCells[i].end();++it){
			printf("\t pid: %d: %d\n",it->first,it->second);	
		}
	}
	delete []_boundCells;
 

}



void RootProcess::writeMonitorFile(DataPartition* dg,const char* chs){
	if(dg->comRank!=rank) return; //only in root
	monitorFile<<chs;
	monitorFile.flush();
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
int RootProcess::getBuffer(DataPartition* dg, int pid,int* sendSizes, double** vbuffer, int** ebuffer, int** ibuffer){
	if(dg->comRank!=rank) return 0; //root only
	map<int,int> nodesPool;
	sendSizes[0] = getVertexSendBuffer(   pid, vbuffer,&nodesPool);
	sendSizes[1] = getElementSendBuffer(  pid, ebuffer,&nodesPool);
	sendSizes[2] = getInterfaceSendBuffer(pid, ibuffer,&nodesPool);	
	return 0;
}


int RootProcess::getVertexSendBuffer(int pid, double** buffer,map<int,int>* nodesPool){



	size_t iEhead = 0;
	for(int i=0;i!=pid;++i){
		iEhead+=rootgridList.at(i);
	}
	size_t iE = iEhead;
	for(int j=0;j!=rootgridList.at(pid);++j){
		InputElement* _thisEle = rootElems[iE++];
		int nv = numberOfNodesInElementTypeOf[_thisEle->type];
		for(int k=0;k!=nv;++k){	     //vertex collection
			int _thisVert = _thisEle->vertex[k];
			nodesPool->insert(make_pair<int,int>(_thisVert,0) );
		}	
	}	


	

	// build buffer
	size_t nV = nodesPool->size();
	*buffer = new double[nV*3 + 1];//nodesPoolSize, nodes1.x, nodes1.y, nodes1.z, nodes2.x, nodes2.y, ...

	size_t counter = 0;
	(*buffer)[counter++] = nV;

	size_t localID = 0;
	for(map<int,int>::iterator it = nodesPool->begin(); it!=nodesPool->end(); ++it){
		it->second = localID++; // the second is local id globalId -> localId
		(*buffer)[counter++] = rootVerts[ it->first ].x;
		(*buffer)[counter++] = rootVerts[ it->first ].y;
		(*buffer)[counter++] = rootVerts[ it->first ].z;
	}

	//printf("vertex buffer for %d prepared\n",pid);
	return counter;
}
/*
bool elementIsOnPeriodicBoundary(InputElement* thisEle,const map<int,BdRegion>& thisRegionMap){
	if(thisEle->type==2||thisEle->type==3){
		int bid = thisEle->tag[0];
		map<int,BdRegion>::const_iterator iter = thisRegionMap.find(bid);
		assert(iter!=thisRegionMap.end() );
		if(iter->second.type1==6){
			return true;
		}
	}
	return false;
}
*/
int RootProcess::getElementSendBuffer(int pid ,int** buffer,map<int,int>* nodesPool){


	
	/**************************************************
	 *	transfer element Vertex list to local idx
	 **************************************************/
	size_t iEhead = 0;
	for(int i=0;i!=pid;++i){
		iEhead+=rootgridList.at(i);
	}
	size_t iE = iEhead;
	size_t counter = 0; //calculate the buffer length
	for(int j=0;j!=rootgridList.at(pid);++j){
		InputElement* _thisEle = rootElems[iE++];
		/*
		if(elementIsOnPeriodicBoundary(_thisEle,regionMap)){
			continue;	//should not send periodic boundary, which will not appear in the final mesh;
		}	
		*/
		int nv = numberOfNodesInElementTypeOf[_thisEle->type];
		for(int k=0;k!=nv;++k){ //debug
			assert(nodesPool->find(_thisEle->vertex[k])!=nodesPool->end());//debug
		}
		counter+= (2 + _thisEle->ntag + nv);
	}

	// build buffer
	(*buffer) = new int [ counter ];//

	counter = 0;
	iE = iEhead;
	for(int j=0;j!=rootgridList.at(pid);++j){ //each elements
		InputElement* _thisEle = rootElems[iE++];
		/*
		if(elementIsOnPeriodicBoundary(_thisEle,regionMap)){
			continue;	//should not send periodic boundary, which will not appear in the final mesh;
		}	
		*/
		(*buffer)[counter++] = _thisEle->type;
		(*buffer)[counter++] = _thisEle->ntag;
		for(int k=0;k!=_thisEle->ntag;++k){
			(*buffer)[counter++] = _thisEle->tag[k];
		}	
		for(int k=0;k!=numberOfNodesInElementTypeOf[_thisEle->type];++k){
			int globalID = _thisEle->vertex[k];
			int localID = (*nodesPool)[globalID];//transfer vertex list in Cell to local ordering
			(*buffer)[counter++] = localID;
		}
	}

	//printf("element buffer for %d prepared\n",pid);
	return counter;
}





int RootProcess::getInterfaceSendBuffer(int pid ,int** buffer,map<int,int>* nodesPool){

	// build buffer

	typedef map<int,map<double,vector<int> > > ComplicatedType;

	ComplicatedType _connectionMap;//<interfaceID <thisID+thatID,<thisID,thatID> > > 
							       //must ensure that the order in  connection map of both side of interface are the SAME!
	
	size_t iEhead = 0;
	for(int i=0;i!=pid;++i)
		iEhead+=rootgridList.at(i);
	size_t iE = iEhead;
	for(int j=0;j!=rootgridList.at(pid);++j){
		InputElement* _thisEle = rootElems[iE++];
		for(vector<vector<int> >::iterator iter = _thisEle->interfaceInfo.begin();iter!=_thisEle->interfaceInfo.end();++iter){
			int interfaceID = iter->at(0);
			int thisID = _thisEle->idx;
			int thatID = iter->at(1);
			if(interfaceID != pid){
				double key = (double)(thisID+thatID)+(double)fabs(thisID-thatID)/(rootNGlobal);
				vector<int> vec(1,thisID - iEhead);// thisID - iE is the real local Index
				for(size_t i=2;i!=iter->size();++i){
					vec.push_back (iter->at(i) );
				}
				_connectionMap[interfaceID].insert(make_pair(key,vec));
			}else{ // update to fix periodic bc;
				//thisID < thatID
				double key = (double)CYCASMIN(thisID,thatID);
				vector<int> vec(1,thisID-iEhead);
				for(size_t i=2;i!=iter->size();++i){
					vec.push_back(iter->at(i));
				}
				assert(thisID!=thatID);
				if(thisID<thatID){
					_connectionMap[-1].insert(make_pair(key,vec));
				}else{
					_connectionMap[-2].insert(make_pair(key,vec));	
				}

			}

		}
	}	
	//if(pid==11) cout<<_connectionMap[10].size()<<endl;

	size_t nI = 0;
	nI++;	// number of interfaces
	//the data structure of this buffer:
	//printf("partition  %d:\n",pid);
	for(ComplicatedType::iterator iter = _connectionMap.begin(); iter!=_connectionMap.end(); ++iter){
		nI+=2;	//interfaceID, interfaceWidth;
		for(map<double,vector<int> >::iterator iterin = iter->second.begin();iterin!=iter->second.end();++iterin){
			nI += (iterin->second.size() + 1); //elemnInfoLeng, thisID, vert1, vert2,...
		}
	//	printf("--->%d : %d\n",iter->first,iter->second.size());
	}


	//debug
	/*
	iEhead = 0;
	InputElement** cells = new InputElement* [rootNCells[pid]];
	for(int i=0;i!=pid;++i){
		iEhead+=rootgridList[i];
	}
	int counterD = 0;
	for(int i=iEhead;i!=rootgridList[pid];++i){
		if(rootElems[i]->idx>=0){
			cells[counterD++] = rootElems[i];
		}
	}

	for(auto& iter : _connectionMap){
		int interfaceID = iter.first;
		printf("pid %d -> %d width: %d\n",pid,interfaceID,iter.second.size());
	}
	*/
	
	/*
	map<double,vector<int> >::iterator iter1 = _connectionMap[-1].begin();
	map<double,vector<int> >::iterator iter2 = _connectionMap[-2].begin();

	for(;iter2!=_connectionMap[-2].end();iter2++,iter1++){
		cout<<iter1->second[0]<<"--->"<<iter1->second[0]<<endl;
	}
	*/	


	*buffer = new int[nI];

	size_t counter = 0;	
	(*buffer)[counter++] = _connectionMap.size();


	for(ComplicatedType::iterator iterOut = _connectionMap.begin(); iterOut!=_connectionMap.end(); ++iterOut){
		int interfaceID = iterOut->first;
		(*buffer)[counter++] = interfaceID; 		//interfaceID
		(*buffer)[counter++] = iterOut->second.size();	//interfaceWidth 
		

		for(map<double,vector<int> >::iterator iterIn = iterOut->second.begin(); iterIn!=iterOut->second.end(); ++iterIn){
			 //iterIn->first is meaningless 
			(*buffer)[counter++] = iterIn->second.size();//
			(*buffer)[counter++] = iterIn->second[0];//in_Process ID
			for(size_t i=1;i!=iterIn->second.size();++i){
				int globalID = iterIn->second[i];
				int localID = (*nodesPool)[globalID];
				assert(nodesPool->find(globalID)!=nodesPool->end());
				(*buffer)[counter++] = localID; 
			}

		}
	}

	assert(counter==nI);
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

void RootProcess::printEnding(DataPartition* dg,double sec){
	if(dg->comRank!=rank) return;

	printf( " -------------------------------------------------------------------------------- \n\n");
	printf( "                      CPU REQUIREMENTS OF NUMERICAL SOLUTION\n\n");
	printf( " -------------------------------------------------------------------------------- \n");
	printf( "    The NSSolve used: %10e s of CPU time\n", sec);
}
void RootProcess::printSolutionNotGood(DataPartition* dg){
	if(dg->comRank!=rank) return;
	printf("**************************************************\n");	
	printf("********   CYCAS didnt get a good result    ******\n");	
	printf("********   Iteration exceed step limitaion  ******\n");	
	printf("**************************************************\n");	
}
void RootProcess::printStepStatus(DataPartition*dg, int step,int piter ,double time,double dt,double res){
	if(dg->comRank!=rank) return;
	printf("%15f\t%10d\t%10d\t%13.4f\t%15.10f\n",time,step,piter,dt,res);
}
void RootProcess::printSteadyStatus(DataPartition*dg,int outiter,double res){
	if(dg->comRank!=rank) return;
	printf("%15s\t%10s\t%10d\t%13s\t%15.10f\n","---","---",outiter,"---",res);
}
void RootProcess::printSectionHead(DataPartition* dg){
	if(dg->comRank!=rank) return;
		printf("%15s\t%10s\t%10s\t%13s\t%15s\n","TIME","CAL STEP","ITER","DELT","MAX RES");
}
