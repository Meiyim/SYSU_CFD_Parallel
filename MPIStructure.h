#include <fstream>
#include <stdio.h>
#include <vector>
#include <string>
#include <stdexcept>
#include <petscksp.h>
#include "NS/TerminalPrinter.h"
#ifndef _DATA_PROCESS_H_
#define _DATA_PROCESS_H_

#define MAX_ROW 5

#ifdef SHOULD_CHECK_MPI_ERROR 

#define CHECK(err) throwError("MPI_ERROR");
#define PCHECK(err) throwError("PETSC_ERROR");
#else
#define CHECK(err)
#define PCHECK(err)

#endif
/******************************************************
 *  	this data group lies each on a processor
 *  	a simple pack of u,v,w, etc.
 ******************************************************/
class RootProcess;
struct DataPartition{ 
	Vec u;
	Vec bu;
	Mat Au;
	PetscErrorCode ierr;
	MPI_Comm comm;
	/***********KSP CONTEXT***************/
	KSP ksp;	
	PC pc;
	/*************************************/
	int mpiErr;
	int comRank;
	int comSize;
	int nLocal;
	int nGlobal;
	int nProcess;
	int* gridList;  //size of nGlobal , gridList[comRank] == nLocal;
	double** Avals; //temp for file reading test;
	int** Aposi; 
	DataPartition():
		comm(MPI_COMM_WORLD),
		nLocal(0),
		nGlobal(0),
		nProcess(0),
		gridList(NULL),
		Avals(NULL), //preassume the size
		Aposi(NULL),
		errorCounter(0)
	{
		mpiErr = MPI_Comm_rank(comm,&comRank);CHECK(mpiErr)
		mpiErr = MPI_Comm_size(comm,&comSize);CHECK(mpiErr)
	}
	~DataPartition(){}
	int init(RootProcess& root); // communicate to get local size on each processes , collective call

	int deinit(); // a normal deconstructor seems not working in MPI

	int fetchDataFrom(RootProcess& root);

	int pushDataTo(RootProcess& root);

	int buildMatrix();

	int solveGMRES(double tol,int maxIter); //return 0 if good solve, retrun 1 if not converge

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
private:
	int pushVectorToRoot(const Vec& petscVec, double* rootbuffer,int rootRank);

};



/******************************************************
 *  	each interface correspond to a edge of the connectivity graph
 ******************************************************/
class Interface{
public:
	int width;
	int sendRank;
	int recvRank;
	/****************FOR PROPERTY U************************/
	int* recvlist; //point to off-process cells
	int* sendlist; //point to on-process cells
	
	// MPI non blocking communication method !
	int push(); // non-blocking method!! return immediatly
	int update(); // non-blocking method!! return immediatly
	// Constructor
	Interface(int selfRank,int rRank,int w): //recvlist & sendlist must be setted by Caller dataGroup
		width(w),
		sendRank(selfRank),	
		recvRank(rRank)
	{}	
	~Interface(){
		delete recvlist;
		delete sendlist;
		recvlist = sendlist = NULL;
	}
	
};


/******************************************************
 *  	this class should work only in main processor
 *  	FUNCITON: read, partition, reorder the mesh in root process
 *  	FUNCTION: gather data and write output file
 *  	controls all the I/O
 ******************************************************/
struct InputElement{// this small structure is only used in this file
	int tag; //tag auxilary to sorting, reordering ,making boundary, etc;
	int type;
	int tag[2];	//the length is fixed in order to send through MPI
	int vertex[8];	
};

struct InputVert{
	double x,y,z;	
};

class RootProcess{
public:
	int rank; //rank of the root;
	/*******************this part should only own by ROOT*****/ 
	int rootNGlobal;	//size of a global vector
	int rootNElement;  	//number of cells and bounds element
	int rootNVert;		// upperbound of int on a 32 bit machine is 2 billion, which is enough
	double* rootArrayBuffer; //used to collect data
	InputElement* rootElems;   //used when partitioning
	InputVert* rootVerts;
	std::vector<int>* rootgridList;

	TerminalPrinter* printer;
	/*********************************************************/
	RootProcess(int r):
		rank(r),
		rootNGlobal(-1),
		rootNVert(-1),
		rootArrayBuffer(NULL), //NULL if not root
		rootElems(NULL),
		rootVerts(NULL),
		rootgridList(NULL),
		printer(NULL)
	{}
	~RootProcess(){
		clean();
		delete printer;
		printer=NULL;
	}
	void init(DataPartition* dg); //init for patitioning
	void allocate(DataPartition* dg); //prepare for gathering
	void clean(); 			  //clean when partition is done;
	/***************************************************
	 * 	 root Only
	 * *************************************************/
	void read(DataPartition* dg,const string& title);
	/***************************************************
	 * 	 root Only
	 * *************************************************/
	void partition(int N,DataPartition* dg);

	/***************************************************
	 * 	 writing result to root
	 * *************************************************/
	void write(DataPartition* dg);
};
#endif
