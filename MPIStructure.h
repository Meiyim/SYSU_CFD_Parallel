#include <fstream>
#include <stdio.h>
#include <vector>
#include <map>
#include <unordered_set>
#include <string>
#include <stdexcept>
#include <petscksp.h>
#include "NS/TerminalPrinter.h"
#ifndef _DATA_PROCESS_H_
#define _DATA_PROCESS_H_

#define MAX_ROW 5
#define ELEMENT_TAG_LEN 2

#ifdef SHOULD_CHECK_MPI_ERROR 
#define CHECK(err) throwError("MPI_ERROR");
#define PCHECK(err) throwError("PETSC_ERROR");
#else
#define CHECK(err)
#define PCHECK(err)
#endif

#define PRINT_LOG(VAR) printlog(VAR,#VAR)

class RootProcess;
class Interface;
class DataPartition;


/******************************************************
 *  	each interface correspond to a edge of the connectivity graph
 ******************************************************/
class Interface{
public:
	// Constructor
	Interface(int selfRank,int rRank,int w): //recvlist & sendlist must be setted by Caller dataGroup
		width(w),
		sendRank(selfRank),	
		recvRank(rRank),
		sendposis(new int[w]),
		sendBuffer(new double[w])
	{}	
	~Interface(){
		delete sendposis;
		delete sendBuffer;
		sendposis = NULL;
		sendBuffer = NULL;
	}
	// MPI non blocking communication method !
	int push(); 		 // non-blocking method!! return immediatly
	int update(); 		 // non-blocking method!! return immediatly
public:
	int width;
	int sendRank;
	int recvRank;

	/****************************************/
	int recvposi; 		 //index of the head of off_process cells
	int const * sendposis;   //indexes of on_process cells
private:	
	double const *sendBuffer;//copy to this buffer and send;
};


/******************************************************
 *  	this data group lies each on a processor
 *  	a simple pack of u,v,w, etc.
 ******************************************************/
class DataPartition{ 
public:
	/***********PETSC*********************/
	Vec u;
	Vec bu;
	Mat Au;
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
	int nGlobal;    // number of global cell
	int nProcess;   // number of partitions
	int* gridList;  //size of nProcess , gridList[comRank] == nLocal;

	/**********INTERFACE INFO**************/
	vector<Interface*> interfaces; //remember to free pointer

	DataPartition():
		comm(MPI_COMM_WORLD),
		//nLocal(0),
		nGlobal(0),
		nProcess(0),
		gridList(NULL),
		errorCounter(0)
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
		gridList = NULL;
		for(vector<Interface*>::iterator it = interfaces.begin();it!=interfaces.end();++it){
			if((*it)!=NULL)
				delete (*it);
		}
	}

	int initPetsc(); //build petsc vectors, collective call

	int deinit(); // a normal deconstructor seems not working in MPI

	int fetchDataFrom(RootProcess& root);  //MPI_ScatterV

	int pushDataTo(RootProcess& root);

	int buildMatrix();

	int solveGMRES(double tol,int maxIter); //return 0 if good solve, retrun 1 if not converge


	std::ofstream logfile; //for test purpose
	template<typename T>
	void printlog(const T& var,const char* varname){
		char temp[256];
		sprintf(temp,"RANK: %d, %s :",comRank,varname);
		logfile<<temp;
		logfile<<var<<std::endl;

	}

private:
	int errorCounter;
	void throwError(const std::string& msg){ // only check error when SHOULD_CHECK_MPI_ERROR is defined
		char temp[256];
		printf("*****************************ON ERROR***************************************");
		sprintf(temp,"MPI_Failure in process %d\n message: %s\n",comRank,msg.c_str()); 
		printf("error msg: %s",temp);
		printf("error in call: %d",errorCounter++);
		printf("*****************************END ERROR***************************************");
		throw std::runtime_error(temp);
	}
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
	std::map<int, unordered_set<int> >* boundNodesPool ; //one for each partition, <partID,nodesPool>
	std::map<int,int>* nodesPool;			     //one for each partition, <globalID, localID>;

	std::vector<int>* rootgridList;	// the element number of each parition
	std::vector<int>* rootNCells;	// the cell number of each partition
	TerminalPrinter* printer;
	/*********************************************************/
	RootProcess(int r):
		rank(r),
		rootNGlobal(-1),
		rootNVert(-1),
		rootArrayBuffer(NULL), //NULL if not root
		rootElems(NULL),
		rootVerts(NULL),
		boundNodesPool(NULL),
		nodesPool(NULL),
		rootgridList(NULL),
		printer(NULL)
	{}
	~RootProcess(){
		clean();
		delete printer;
		printer=NULL;
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
	 *********************************************/
	
	int getElementSendBuffer(DataPartition* dg,int pid, int** buffer); 	 
	int getVertexSendBuffer(DataPartition* dg,int pid, double** buffer);	
	int getInterfaceSendBuffer(DataPartition* dg,int pid, int** buffer);   



	/***************************************************
	 * 	 writing result to root
	 * *************************************************/
	void write(DataPartition* dg);

	void printTecplotHead(DataPartition* dg, std::ofstream& file);

private:

};
#endif
