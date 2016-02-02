#include <fstream>
#include <stdio.h>
#include <vector>
#include <string>
#include <stdexcept>
#include <petscksp.h>
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
struct DataGroup{ 
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
	DataGroup():
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
	int init(RootProcess& root); // communicate to get local size on each processes , collective call
	int deinit(){ // a normal deconstructor seems not working in MPI
		printf("datagroup NO. %d died\n",comRank);
		ierr = VecDestroy(&u);CHKERRQ(ierr);
		ierr = VecDestroy(&bu);CHKERRQ(ierr);
		ierr = MatDestroy(&Au);CHKERRQ(ierr);
		ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
		delete []Avals;
		delete []Aposi;
		delete []gridList;
		Avals=NULL;
		gridList=NULL;
		Aposi=NULL;
		return 0;
	}
	~DataGroup(){
	}
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
 ******************************************************/
class RootProcess{
public:
	int rank; //rank of the root;
	/*******************this part should only own by ROOT*****/ 
	int rootNGlobal;
	double* rootuBuffer;
	double* rootbuBuffer;
	double* rootABuffer;
	int* rootAPosiBuffer;
	std::vector<int>* rootgridList;
	/*********************************************************/
	RootProcess(int r):
		rank(r),
		rootNGlobal(-1),
		rootuBuffer(NULL), //NULL if not root
		rootbuBuffer(NULL),
		rootABuffer(NULL),
		rootAPosiBuffer(NULL),
		rootgridList(NULL)
	{}
	~RootProcess(){
		clean();
	}
	void allocate(DataGroup* dg){ //prepare for gather;
		if(dg->comRank!=rank) return; //only in root
		rootuBuffer = new double[rootNGlobal];
		rootbuBuffer = new double[rootNGlobal];
	}
	void clean(){
		delete rootuBuffer;
		delete rootbuBuffer;
		delete []rootABuffer;
		delete []rootAPosiBuffer;
		rootuBuffer=rootABuffer=NULL;
		rootAPosiBuffer = NULL;
	}
	/***************************************************
	 * 	 this funciton should involk only in root process
	 * *************************************************/
	void read(DataGroup* dg){ 
		if(dg->comRank!=rank) return; //only in root
		printf("start reading in root\n");
		char ctemp[256];
		int itemp=0;
		int maxRow = 0;
		std::ifstream infile("Au.dat");
		infile>>ctemp>>rootNGlobal;
		infile>>ctemp>>itemp;
		infile>>ctemp>>maxRow;

		allocate(dg);
		rootABuffer = new double[rootNGlobal*MAX_ROW];
		rootAPosiBuffer = new int[rootNGlobal*MAX_ROW];

		for(int i=0;i!=rootNGlobal*MAX_ROW;++i){
			rootABuffer[i] = 0.0;
			rootAPosiBuffer[i] = -1; //-1 means no zeros in mat;
		}
		for(int i=0;i!=rootNGlobal;++i){
			size_t nCol = 0;
			double dump;
			double posi;
			double vals;
			double bVal;
			infile>>nCol>>ctemp;
			for(int j=0;j!=nCol;++j){
				infile>>posi>>vals;
				posi--;
				rootABuffer[i*MAX_ROW+j] = vals;
				rootAPosiBuffer[i*MAX_ROW+j] = posi;
			}
			infile>>ctemp>>bVal>>dump;
			rootuBuffer[i] = bVal;
		}
		infile.close();
		printf("complete reading in root\n");
		printf("now the input array is \n");
	}
	/***************************************************
	 * 	 this funciton should involk only in root process
	 * *************************************************/
	void partition(int N,DataGroup* dg){
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

	/***************************************************
	 * 	 writing result to root
	 * *************************************************/
	void write(DataGroup* dg){
		if(dg->comRank!=rank) return; //only in root
		printf("writing....");
		std::ofstream outfile("result.dat");
		char temp[256];
		for(int i=0;i!=rootNGlobal;++i){
			sprintf(temp,"%15d\tx:%15e\tb:%15e\n",i+1,rootuBuffer[i],rootbuBuffer[i]);
			outfile<<temp;
		}
		outfile.close();
	}
};
#endif
