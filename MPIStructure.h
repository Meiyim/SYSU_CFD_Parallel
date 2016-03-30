#include <fstream>
#include <stdio.h>
#include <vector>
#include <map>
#include <set>
#include <string>
#include <stdexcept>
#include <petscksp.h>
#include "NS/BasicType.h"

#ifndef CYCAS_DEBUG_MODE
#define NDEBUG
#endif
#include <assert.h>


#ifndef _DATA_PROCESS_H_
#define _DATA_PROCESS_H_


#define MAX_ROW 5
#define MAX_LOCAL_PREALLOCATION 7
#define ELEMENT_TAG_LEN 2

#define CHECK(err) 

#define PRINT_LOG(VAR) printlog(VAR,#VAR)

class RootProcess;
class Interface;
class DataPartition;
class CellData;

using std::vector;
using std::map;
using std::set;
using std::string;

/******************************************************
 *  	each interface correspond to a edge of the connectivity graph
 ******************************************************/
class Interface{
public:
	// Constructor
	Interface():
		selfRank(-1),
		otherRank(-1),
		sendBuffer(NULL),
		sendBufferGradient(NULL),
		sendBufferCell(NULL)
	{}

	Interface(int s,int o,MPI_Comm c): //recvlist & sendlist must be setted by Caller dataGroup
		selfRank(s),	
		otherRank(o),
		sendBuffer(NULL),
		sendBufferGradient(NULL),
		sendBufferCell(NULL),
		comm(c)
	{
		CellData cellSample;
		int lena[2];
		MPI_Aint loca[2];
		MPI_Datatype typa[2];
		MPI_Aint baseAddress;

		MPI_Get_address(&cellSample,&baseAddress);

		lena[0] = 22; //22 Ints
		MPI_Get_address(&cellSample.nface,&loca[0]);
		loca[0] -= baseAddress;
		typa[0] = MPI_INT;

		lena[1] = 4; //4  Double
		MPI_Get_address(&cellSample.vol,&loca[1]);
		loca[1] -= baseAddress;
		typa[1] = MPI_DOUBLE;
		
		MPI_Type_create_struct(2,lena,loca,typa,&MPI_CellData);
		MPI_Type_commit(&MPI_CellData);

	}	

	~Interface(){
		if(sendBuffer!=NULL) delete [] sendBuffer;
		if(sendBufferGradient!=NULL) delete [] sendBufferGradient;
		if(sendBufferCell!=NULL) delete []sendBufferCell;
	}

	// MPI non blocking communication method !
	int send(MPI_Request* ,CellData* phi,int tag); 		 // non-blocking method!! return immediatly
	int recv(MPI_Request* ,CellData* phi,int tag); 		 // non-blocking method!! return immediatly
	int send(MPI_Request* ,double* phi, int tag); 		 // non-blocking method!! return immediatly
	int recv(MPI_Request* ,double* phi, int tag); 		 // non-blocking method!! return immediatly
	int send(MPI_Request* ,double* phi[3], int tag); 		 // non-blocking method!! return immediatly
	int recv(MPI_Request* ,double* phi[3], int tag); 		 // non-blocking method!! return immediatly

	void getData(CellData* phi);
	void getData(double* phi);
	void getData(double* phi[3]);

	size_t getWidth(){return sendposis.size();}
public:
	/****************************************/
	vector<set<int> > boundNodes;
	vector<int> sendposis;//indexes of cells needs to communicate
	int recvposi;         //index head of position > Ncel

private:	
	int selfRank;
	int otherRank;
	double *sendBuffer;		//copy to this buffer and send;
	double *sendBufferGradient;
	CellData *sendBufferCell;
	MPI_Comm comm;
	MPI_Datatype MPI_CellData; //a variable to hold the customized MPI TYPE, contains 22 Ints, 4 Doubles
};



/******************************************************
 *  	this data group lies each on a processor
 *  	a simple pack of u,v,w, etc.
 ******************************************************/
class DataPartition{ 
public:
	/***********PETSC_HANDLE***************/
	Vec bu;
	Vec bv;
	Vec bw;
	Vec bp;
	Vec bs;   //universal scarlar solver


	Vec xsol; //determined when solve

	Mat Au;
	Mat Ap;
	Mat As;

	PetscErrorCode ierr;
	MPI_Comm comm;

	/***********KSP CONTEXT***************/
	KSP ksp;	
	PC pc;

	/***********MPI_PARAMETER*************/
	int mpiErr;
	int comRank;
	int comSize;
	int nLocal; 	// number of local cell, same as Ncel
	int nVirtualCell;
	int nGlobal;    // number of global cell
	int nProcess;   // number of partitions
	int* gridList;  //size of nProcess , gridList[comRank] == nLocal;

	vector<MPI_Request> requests;
	int tagCounter; //used in interface communication
	/**********INTERFACE INFO**************/
	map<int,Interface> interfaces; //<partID,Interface>

	DataPartition():
		comm(MPI_COMM_WORLD),
		nLocal(0),
		nVirtualCell(0),
		nGlobal(0),
		nProcess(0),
		gridList(NULL),
		tagCounter(0)
		
	{
		mpiErr = MPI_Comm_rank(comm,&comRank);CHECK(mpiErr)
		mpiErr = MPI_Comm_size(comm,&comSize);CHECK(mpiErr)
		char temp[256];
		sprintf(temp,"log/log%d",comRank);
		logfile.open(temp);
	}
	
	~DataPartition(){
		/*******DESTROY PETSC OBJECTS**********/
		logfile.close();
		delete []gridList;	
	}

	int initPetsc(); //build petsc vectors, collective call

	int deinit(); // a normal deconstructor seems not working in MPI

	int buildInterfaceFromBuffer(int* buffer);

	int fetchDataFrom(RootProcess& root);  //MPI_ScatterV

	int pushDataTo(RootProcess& root);

	int solveVelocity_GMRES(double tol,int maxIter,double const* xu,double const* xv,double const* xw); //return 0 if good solve, retrun 1 if not converge

	int solvePressureCorrection(double tol, int maxIter,double const* xp,bool isSymmetric);

	
	/*************MPI INTERFACE COMMUNICATION***********************
	 *	collective
	 ***************************************************************/
	template<typename T>
	int interfaceCommunicationBegin(T var){
		MPI_Barrier(comm);
		size_t nInter = interfaces.size();
		MPI_Request* sendReq = new MPI_Request[nInter];	
		MPI_Request* recvReq = new MPI_Request[nInter];	
		size_t reqCounter = 0;
		for(map<int,Interface>::iterator iter = interfaces.begin(); iter!=interfaces.end(); ++iter){
			Interface& _inter = iter->second;
			_inter.send(sendReq+reqCounter,var,tagCounter);//non-blocking...
			_inter.recv(recvReq+reqCounter,var,tagCounter);//non-blocking
			reqCounter++;
		}	
		tagCounter++; //same tag in each interface communication

		for(int i=0;i!=nInter;++i){
			requests.push_back(sendReq[i]);
			requests.push_back(recvReq[i]);
		}
		delete []sendReq;
		delete []recvReq;
		
		return 0;
	}
	int interfaceCommunicationEnd(){
		size_t _requestsSize = requests.size();
		MPI_Request* _requests = &(requests[0]);

		int ierr = MPI_Waitall(_requestsSize,_requests,MPI_STATUSES_IGNORE);
		if(ierr!=MPI_SUCCESS){
			errorHandler.fatalRuntimeError("error occured when recving interface value\n");
		}
		requests.clear();
		tagCounter=0;//clear all communication

		return 0;

	}

	template<typename T>
	void printlog(const T& var,const char* varname){
		char temp[256];
		sprintf(temp,"RANK: %d, %s :",comRank,varname);
		logfile<<temp;
		logfile<<var<<std::endl;

	}

private:
	std::ofstream logfile; //for test purpose
	int pushVectorToRoot(const Vec& petscVec, double* rootbuffer,int rootRank);

};


/******************************************************
 *  	this class should work only in main processor
 *  	FUNCITON: read, partition, reorder the mesh in root process
 *  	FUNCTION: gather data and write output file
 *  	controls all the I/O
 ******************************************************/
class InputElement{// this small structure is only used in this file
public:
	int type;
	int pid; 	//use to sort, no need to send
	int idx;
	vector< vector<int> > interfaceInfo; //<partID,<globalIdx, intersectionNode1,2,3... > >
	int ntag; 
	int* tag;	//the length is fixed in order to send through MPI
	int* vertex;	
	InputElement(int ty, int nt,int nv):type(ty),ntag(nt){
		tag = new int[nt];
		vertex = new int[nv];
	}
	~InputElement(){
		delete []tag;
		delete []vertex;
	}

};


struct InputVert{
	double x,y,z;	
};

class RootProcess{
public:
	int rank; //rank of the root;
	/*******************owned by ROOT*****/ 
	size_t rootNGlobal;		//size of a global vector
	size_t rootNElement;  		//number of cells and bounds element
	size_t rootNVert;		// upperbound of int on a 32 bit machine is 2 billion, which is enough
	double* rootArrayBuffer; 	//used to collect data

	/***************used when partitioning*************************/
	InputElement** rootElems;   	 //a pointer array //for faster sorting
	InputVert* rootVerts; 		

	std::vector<int> rootgridList;	// the element number of each parition
	std::vector<int> rootNCells;	// the cell number of each partition
	/*********************************************************/
	RootProcess(int r):
		rank(r),
		rootNGlobal(-1),
		rootNVert(-1),
		rootArrayBuffer(NULL), //NULL if not root
		rootElems(NULL),
		rootVerts(NULL)
	{}
	~RootProcess(){
		clean();
		delete []rootArrayBuffer;
	}

	void init(DataPartition* dg); 		//init for patitioning

	void allocate(DataPartition* dg); 	//prepare for gathering

	void clean(); 			  	//clean when partition is done;
	/***************************************************
	 * 	 root Only
	 * *************************************************/
	void read(DataPartition* dg,const string& title);

	/***************************************************
	 * 	 root Only
	 **************************************************/
	void partition(DataPartition* dg, int N);



	int getBuffer(DataPartition* dg, int pid,int* sendsizes, double** vbuffer, int** ebuffer, int** ibuffer);


	/***************************************************
	 * 	 print screen
	 * *************************************************/
	void printStarter(DataPartition*);
	void printEnding(DataPartition*);
	void printStepStatus(DataPartition* , int,int,double,double,double);
	void printSteadyStatus(DataPartition*,int,double);
	void printSectionHead(DataPartition*, double);

	/***************************************************
	 * 	 writing result to root
	 * *************************************************/
	void write(DataPartition* dg);


private:
	/*********************************************
	 * 	build MPI transfer buffer for vertex  element interfaceinfo
	 * 	
	 * 	WARNING: this 3 function MUST be called in specifig sequence!
	 * 	calling sequence :
	 * 		getvertex[pid]
	 * 		getelement[pid]
	 * 		getinterface[pid]
	 *
	 * 	WARNING: it is the user's responsibility to free buffer return by these functions
	 *********************************************/
	int getElementSendBuffer(int pid, int** buffer,map<int,int>* nodesPool); 	 
	int getVertexSendBuffer(int pid, double** buffer,map<int,int>* nodesPool);	
	int getInterfaceSendBuffer(int pid, int** buffer,map<int,int>* nodesPool);   

};
#endif
