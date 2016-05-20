
#include <math.h>
#include <stdlib.h>
#include <sstream>
#include <limits.h>
#include "navier.h"
#include "tools.h"
#include "terminalPrinter.h"


//length of options pool


using namespace std;

void NavierStokesSolver::NSSolve( )
{
	int iter;
	double ResMax,ResMax0;
	timespec tstart ={0,0};
	timespec tend ={0,0};
	
	//profiling
	PetscLogStage buildVelocityStage = 1;
	PetscLogStageRegister("BuildVelocityStage",&buildVelocityStage);
	PetscLogStage solveVelocityStage = 2;
	PetscLogStageRegister("SolveVelocityStage",&solveVelocityStage);
	PetscLogStage buildPressureStage = 3;
	PetscLogStageRegister("BuildPressureStage",&buildPressureStage);
	PetscLogStage solvePressureStage = 4;
	PetscLogStageRegister("SolvePressureStage",&solvePressureStage);
	PetscLogStage otherStage = 5;
	PetscLogStageRegister("otherStage",&otherStage);

	PetscLogStagePush(otherStage);//Profilling


	MPI_Barrier(dataPartition->comm);
	CYCAS_GET_TIME(tstart);

	cur_time = 0.0;
	if(IfSteady)
		dt = CYCASHUGE_D;

	cur_time += dt;
	if(IfSteady){
		MaxOuterStep = MaxStep;
	}
	//MaxOuterStep=2;//test
    for( step=1; step<=MaxStep; step++ ) //step : total time step
    {
		if( !IfSteady ){
			SaveTransientOldData( );
		}

		// outer iteration, drive residual to zero for every physical time step
		for( iter=1; iter<MaxOuterStep; iter++ ){

			//------------   Record   ------------//
			if( shouldBackup(step,iter,cur_time) )
				WriteBackupFile();
			if( shouldPostProcess(step,iter,cur_time) ){
				Output2Tecplot();
				root.printSectionHead(dataPartition);
			}

			//----------- record,tot file, restart file, etc.----//
			writeTotFile();
				
			//------------   SIMPLE_C ------------//
			MPI_Barrier(dataPartition->comm);
				
			std::fill(localRes,localRes+RESIDUAL_LEN,0.0);
			std::fill(Residual,Residual+RESIDUAL_LEN,0.0);
			
			PetscLogStagePop();//profilling
			CalculateVelocity ( ); //interface communication U, V, W, Apr
			
			CalculatePressure ( ); //calculate deltaP and correct P,[R], U,V,W
					       //interface communication U,V,W,P,R

			PetscLogStagePush(otherStage);//Profilling
			
			//------------   Physics Models   ------------//
			// scalar transportation
			//1. turbulence model
			if(TurModel==1) {
				UpdateTurKEpsilon( );
			}
			//CHECK_ARRAY(TE,Ncel);
			//CHECK_ARRAY(ED,Ncel);

			//2. energy couple
			if( SolveEnergy  ) {
				UpdateEnergy ( );
			}

			//3. species transport
			if( SolveSpecies ) UpdateSpecies( );//to be implemented
			//4. other physical models
			//    e.g., condensation, combustion, wall heat conduction

			/*-----------check if should break----------*/
			MPI_Allreduce(localRes,Residual,RESIDUAL_LEN,MPI_DOUBLE,MPI_SUM,dataPartition->comm);

			ResMax = vec_max( Residual,RESIDUAL_LEN );
			
			if( IfSteady ){
				//steady
				root.printSteadyStatus(dataPartition,iter,ResMax);
				if( ResMax<ResidualSteady )break;
			}else{
				//unsteady
				root.printStepStatus(dataPartition,step,iter,cur_time,dt,ResMax);
				if( iter == 1 ) ResMax0 = ResMax;
				if( ResMax<1.e-4 || ResMax0/(ResMax+1.e-16)>1000. ){
					break; // more reasonal to break : order drop 2
				}
			}

		}
		if(iter==MaxOuterStep){
			PetscLogStagePop();//profilling
			CYCAS_GET_TIME(tend);
			root.printSolutionNotGood(dataPartition);
			break;
		}		


		if(cur_time >= total_time){
			PetscLogStagePop();//profilling
			CYCAS_GET_TIME(tend);
			break;
		}else{//time advance
			cur_time+=dt;
		}	
	}

	//extra work before solve is complete 
	Output2Tecplot();

	root.printEnding(dataPartition,
			tend.tv_sec-tstart.tv_sec,
			tend.tv_nsec-tstart.tv_nsec);

	return;
}


/*******************************************************/
//	determint if should output to tecplot, vtk...
/*******************************************************/
bool NavierStokesSolver::shouldBackup(int timestep,int outiter,double now){
	if(IfSteady){
		return (outiter-1)%noutput == 0;
	}else{
		return false;
	}

}


/*******************************************************/
//	determint if should backUp
//	currently the same frequency as post process
/*******************************************************/
bool NavierStokesSolver::shouldPostProcess(int timestep,int outiter, double now){
	if(IfSteady){
		return (outiter-1)%noutput == 0;
	}else{
		return (timestep-1)%noutput==0 && outiter == 1;
	}
}


/*******************************************************/
//	fetch param and mesh Data From root
//	MPI_Broadcast routine
//	collective
/*******************************************************/
void NavierStokesSolver::broadcastSolverParam(){
	PetscPrintf(dataPartition->comm,"Broadcasting basic parameter...");
	fflush(stdout);
	int sourceRank = root.rank;

	MPI_Barrier(dataPartition->comm);

	//------bool option pool
	MPI_Bcast(iOptions,INT_OPTION_NO,MPI_INT,sourceRank,MPI_COMM_WORLD);

	//------double option pool
	MPI_Bcast(dbOptions,DB_OPTION_NO,MPI_DOUBLE,sourceRank,MPI_COMM_WORLD);

	//------region map;
	int nRegion = regionMap.size();
	int bufferSize = 0;
	int mapIndex = 0;
	double* regionBuffer = new double[1024];//preassume buffer lenth

	MPI_Bcast(&nRegion,1,MPI_INT,sourceRank,MPI_COMM_WORLD);

	std::map<int,BdRegion>::iterator iter = regionMap.begin();
	for(int i=0;i!=nRegion;++i){
		if(dataPartition->comRank == root.rank){
			bufferSize = iter->second.getSendBuffer(regionBuffer);
			mapIndex = iter->first;
			iter++;
		}
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Bcast(&mapIndex,1,MPI_INT,sourceRank,MPI_COMM_WORLD);
		MPI_Bcast(&bufferSize,1,MPI_INT,sourceRank,MPI_COMM_WORLD);
		MPI_Bcast(regionBuffer,bufferSize,MPI_DOUBLE,sourceRank,MPI_COMM_WORLD);

		//printf("initing region(%d) of size %d \n",mapIndex,bufferSize );
		if(dataPartition->comRank != root.rank){
			regionMap[mapIndex] = BdRegion(regionBuffer,bufferSize);
		}


	}

	delete []regionBuffer;
	MPI_Barrier(MPI_COMM_WORLD);
	PetscPrintf(MPI_COMM_WORLD,"done\n");

}


/*******************************************************/
//	fetch mesh Partition info From root
//	nProcess, gridlist, Ncel, Nbnd...
//	MPI_Broadcast routine
//	collective
/*******************************************************/
void NavierStokesSolver::broadcastPartitionInfo(){
	PetscPrintf(dataPartition->comm,"Broadcasting partition info...");
	fflush(stdout);
	int sourceRank = root.rank;
	int* _sendbuf = NULL;

	if(dataPartition->comRank == root.rank){//root only
		dataPartition->nProcess = root.rootgridList.size();
		dataPartition->nGlobal = root.rootNGlobal;
	}

	MPI_Bcast(&(dataPartition->nProcess),1,MPI_INT,sourceRank,dataPartition->comm);//nProcess

	MPI_Bcast(&(dataPartition->nGlobal),1,MPI_INT,sourceRank,dataPartition->comm);//nGlobal

	dataPartition->gridList = new int[dataPartition->nProcess];

	if(dataPartition->comRank == root.rank){//root only
		for(int i=0;i!=dataPartition->nProcess;++i)
			dataPartition->gridList[i] = root.rootgridList[i];

		_sendbuf = new int[dataPartition->nProcess]; //only significant in root

		for(int i=0;i!=dataPartition->nProcess;++i)
			_sendbuf[i] = root.rootNCells[i];

	}
	
	MPI_Bcast(dataPartition->gridList, dataPartition->nProcess, MPI_INT, sourceRank, dataPartition->comm);//gridlist

	//scatter: sendbuf, sendcount(number of element to EACH process), sendtype, recvbuf, recvcount, recvtype, rootrank, communicator
	MPI_Scatter(_sendbuf,1, MPI_INT, &(this->Ncel), 1, MPI_INT,root.rank,dataPartition->comm  );//nCell


	//some configuration...
	int _nlocalElement = dataPartition->gridList[dataPartition->comRank];

	this->Nbnd = _nlocalElement - Ncel;
	dataPartition->nLocal = Ncel;


	/*
	printf("check parameter\n");
	printf("vert------>%d\n",Nvrt);
	*/
	/*
	printf("rank %d, nProcess %d, nLocal %d, nGlobal %d, nCell %d, nBnd %d\n",
			dataPartition->comRank, dataPartition->nProcess, dataPartition->nLocal, dataPartition->nGlobal, Ncel, Nbnd);
	*/
	delete []_sendbuf;
	MPI_Barrier(dataPartition->comm);
	PetscPrintf(dataPartition->comm,"done\n");
}


/*******************************************************/
//	fetch mesh Data From root
//	MPI_Issend routine
//	collective
//
//	parameter is for input, pass NULL to this routine
//	this routien will malloc buffer
//	it is the USERs responsibility to free this buffer!
/*******************************************************/
#define SEND_TAG_ELEMENT 0
#define SEND_TAG_VERTEX 100
#define SEND_TAG_INTERFACEINFO 200
void NavierStokesSolver::scatterGridFile(int** elemBuffer, double** vertexBuffer, int** interfaceBuffer){
	//scatterv and multiple send operation are both good for this subroutine
	PetscPrintf(dataPartition->comm,"scattering geometry...");
	fflush(stdout);

	MPI_Request* sendRequests = NULL;	
	MPI_Request* recvRequtests = NULL;
	double** buffer1       = NULL; //vertex send buffer
	int** buffer2          = NULL; //elements send buffer
	int** buffer3          = NULL; //interfaces send buffer
	
	int ierr = 0;
	
	int _size = 0;
	int _sendSize[3] = {0,0,0};
	int _recvSize[3] = {0,0,0};

	MPI_Barrier(dataPartition->comm);

	if(dataPartition->comRank==root.rank){//root //sending side...
		_size = dataPartition->nProcess;
		sendRequests = new MPI_Request[ 3 * _size ];
		buffer1 = new double* [_size];
		buffer2 = new int* [_size];
		buffer3 = new int* [_size];

		for(int pid = 0;pid!=_size;++pid){
			//printf("root: sending geometry to grid: %d\n",pid);
			root.getBuffer(dataPartition,pid,_sendSize,buffer1+pid, buffer2+pid, buffer3+pid);
			//printf("prepared buffer for pid %d, %d,%d,%d\n",pid,_sendSize[0],_sendSize[1],_sendSize[2]);

			//non-blocking buffered send
			MPI_Bsend(_sendSize,3,MPI_INT,pid,244,dataPartition->comm);

			//send vertex	
			MPI_Issend(buffer1[pid], _sendSize[0], MPI_DOUBLE,
					pid,		     	     //dest
					SEND_TAG_VERTEX,	     //send tag
					dataPartition->comm,
					sendRequests+3*pid+0); 	     //return request
			//send elements
			MPI_Issend(buffer2[pid], _sendSize[1], MPI_INT,
					pid,		      	     //dest
					SEND_TAG_ELEMENT,	     //send tag
					dataPartition->comm,
					sendRequests+3*pid+1); 	     //return request

			//send interfaceinfo
			MPI_Issend(buffer3[pid], _sendSize[2], MPI_INT,
					pid,			     //dest
					SEND_TAG_INTERFACEINFO,  //send tag
					dataPartition->comm,
					sendRequests+3*pid+2); 	     //return request
		}			




	} 

	//receiveing side ...
	//root is sending to itself...
	//printf("rank: %d is receiving geometry...\n",dataPartition->comRank);
	recvRequtests = new MPI_Request[ 3 ];
	//blocking receive!
	MPI_Recv(_recvSize,3,MPI_INT,root.rank,244,dataPartition->comm,MPI_STATUS_IGNORE);

	
	//printf("rank:%d recving buffersize: %d, %d, %d\n",dataPartition->comRank,_recvSize[0],_recvSize[1],_recvSize[2]);

	*vertexBuffer = new double[_recvSize[0]];
	*elemBuffer = new int[_recvSize[1]];
	*interfaceBuffer = new int[_recvSize[2]];

	//recv elements
	MPI_Irecv(*vertexBuffer,_recvSize[0], MPI_DOUBLE,
			root.rank,			//source
			SEND_TAG_VERTEX,
			dataPartition->comm,
			recvRequtests+0); 

	MPI_Irecv(*elemBuffer,_recvSize[1], MPI_INT,
			root.rank,			//source
			SEND_TAG_ELEMENT,		//tag
			dataPartition->comm,
			recvRequtests+1); 			//return value
		//return value
	MPI_Irecv(*interfaceBuffer,_recvSize[2], MPI_INT,
			root.rank,			//source
			SEND_TAG_INTERFACEINFO,
			dataPartition->comm,
			recvRequtests+2); 			//return value

	MPI_Status status[3];
	ierr = MPI_Waitall(3,recvRequtests,status);
	if(ierr!=MPI_SUCCESS){
		char temp[256];
		sprintf(temp,"MPI receive failure in geometry transfering\n");
		errorHandler.fatalRuntimeError(temp);
	}
	
	//printf("rank: %d received geometry buffer\n",dataPartition->comRank);

	writeGeometryBackup(_recvSize[0],*vertexBuffer,_recvSize[1],*elemBuffer,_recvSize[2],*interfaceBuffer);

	delete []recvRequtests;


	if(dataPartition->comRank == root.rank){//check completion in root
		ierr = MPI_Waitall(3*_size,sendRequests,MPI_STATUS_IGNORE);
		if(ierr!=MPI_SUCCESS){
			errorHandler.fatalRuntimeError("MPI send failure in geometry transfering!\n ");
		}

		for(int i=0;i!=_size;++i){
			delete []buffer1[i];
			delete []buffer2[i];
			delete []buffer3[i];
		}
		delete []buffer1;
		delete []buffer2;
		delete []buffer3;
		delete []sendRequests;
		root.clean();//OK to free root data
	}
	MPI_Barrier(dataPartition->comm);
	PetscPrintf(dataPartition->comm,"done\n");
}



/*******************************************************/
//	inition of MPI related variables
//	the inition of Params can be overwritten by Command line option and input parameters
//	this function is root ONLY
/*******************************************************/
void NavierStokesSolver::initSolverParam() 
{
	if(dataPartition->comRank != root.rank) return; //root ONLY

	root.printStarter(dataPartition);
	int i;
	// default values, can be reset in initflow
	MaxStep      = 10000 ;
	MaxOuterStep = 50  ;
	IfReadBackup = false;
	IfSteady     = true ;  ResidualSteady= 1.e-6;
	TurModel     = 0;
	DensityModel = 0;
	SolveEnergy  = false;
	SolveSpecies = false;
	Nspecies     = 0;

	PressureReference= 1.01325e5; cellPressureRef=0;
	gama			 = 1.4;
	ga1				 = 0.4;
	cp				 = 1006.;
	cv				 = cp/gama;
	prl				 = 0.72;
	prte			 = 0.9;
	Rcpcv			 = cp-cv;
	TempRef			 = 273.15;
	for( i=0; i<3; i++ )
		gravity[i] = 0.;

	//-- numerical scheme
	// relaxation factor
	URF[0]= 0.6;  // u
	URF[1]= 0.6;  // v
	URF[2]= 0.6;  // w
	URF[3]= 0.5;  // p
	URF[4]= 0.8;  // T
	URF[5]= 0.8;  // k
	URF[6]= 0.8;  // e
	URF[7]= 0.8;  // scalar

	limiter = 0;

	total_time   = 0.    ;
	dt           = 1.0   ;    // no meaning in steady case
	TimeScheme   = 1     ;

	noutput  = 1;
	outputFormat = 0;  // 0 for tecplot; 1 for vtk

	// boundary

	/*
	// specific problem parameters
	// cylinder
	IfSteady     = true ;  dt=0.1;  TimeScheme=1;
	SolveEnergy  = false ;
	TurModel     = 1     ;
	DensityModel = 0     ;
	*/

	/*
	// heat exchanger
	IfSteady     = true ;
	MaxStep      = 500  ;
	SolveEnergy  = true ;
	TurModel     = 1    ;
	DensityModel = 1    ; */


	// init parameters from param.in
	
	ReadParamFile( );
	// some parameters check work, e.g. 
	
}


void NavierStokesSolver::InitFlowField( ){
	int i;

	if( IfReadBackup ) 
		ReadBackupFile( );
	else{
		if(regionMap.find(0)!=regionMap.end()){
			for( i=0; i<Ncel; i++ )
			{
				double* initvalues = regionMap[0].initvalues;
				Un[i] = initvalues[0];
				Vn[i] = initvalues[1];
				Wn[i] = initvalues[2];
				Rn[i] = initvalues[3];
				Tn[i] = initvalues[4];
				Pn[i] = 0.;
				if( DensityModel== 1 ) Rn[i]= (Pn[i]+PressureReference)/(Rcpcv*Tn[i]);
				VisLam[i]= initvalues[5]; // 0.6666667e-2;  // 1.458e-6 * pow(Tn[i],1.5) /(Tn[i]+110.5) ;
				VisTur[i]= 0.;
				if( TurModel==1 )
				{
					TE[i]    = initvalues[6];  // 1.e-4*(Un[i]*Un[i]+Vn[i]*Vn[i]+Wn[i]*Wn[i]);
					ED[i]    = initvalues[7];    // TurKEpsilonVar::Cmu * pow(TE[i],1.5) / 1.;
					VisTur[i]= Rn[i]*TurKEpsilonVar::Cmu * TE[i]*TE[i]/(ED[i]+SMALL);
				}
			}
		}

		for( i=0; i<Ncel; i++ )
		{
			int rid = Cell[i].rid;
			assert(regionMap.find(rid)!=regionMap.end());
			double* initvalues = regionMap[rid].initvalues;
			Un[i] = initvalues[0];
			Vn[i] = initvalues[1];
			Wn[i] = initvalues[2];
			Rn[i] = initvalues[3];
			Tn[i] = initvalues[4];
			Pn[i] = 0.;
		    if( DensityModel== 1 ) Rn[i]= (Pn[i]+PressureReference)/(Rcpcv*Tn[i]);
			VisLam[i]= initvalues[5]; // 0.6666667e-2;  // 1.458e-6 * pow(Tn[i],1.5) /(Tn[i]+110.5) ;
			VisTur[i]= 0.;
			if( TurModel==1 )
			{
				TE[i]    = initvalues[6];  // 1.e-4*(Un[i]*Un[i]+Vn[i]*Vn[i]+Wn[i]*Wn[i]);
				ED[i]    = initvalues[7];    // TurKEpsilonVar::Cmu * pow(TE[i],1.5) / 1.;
				VisTur[i]= Rn[i]*TurKEpsilonVar::Cmu * TE[i]*TE[i]/(ED[i]+SMALL);
			}
		}

	}
	std::fill(deltaP,deltaP+Ncel,0.0);

	int ridCount[7] = {0,0,0,0,0,0,0};

	for(i=0;i!=Nbnd;++i){
		assert(regionMap.find(Bnd[i].rid)!=regionMap.end());
		BdRegion& reg = regionMap[Bnd[i].rid];
		ridCount[reg.type1]++;
		if(reg.type1==1){
			if(reg.type2==0){
				BTem[i] = reg.fixedValue;
			}else if(reg.type2 == 1){
				Bnd[i].q = reg.fixedValue;
			}else if(reg.type2 == 2){ //coupled bounday
				Bnd[i].q = 0.0;
			}
		}else if(reg.type1==2){//inlet, u,v,w,p,ro,t,te,ed
			BU[i] = reg.initvalues[0];
			BV[i] = reg.initvalues[1];
			BW[i] = reg.initvalues[2];
			//pressure is ignored  [3]
			BRo[i] = reg.initvalues[4];
			BTem[i]= reg.initvalues[5];
			BTE[i] = reg.initvalues[6];
			BED[i] = reg.initvalues[7];
	 	}else if(reg.type1 == 3){
	 		//pass
		}else if(reg.type1==4){//sym
	 		Bnd[i].q =0.0;//not implement yet
	 	}else if(reg.type1 == 6){//periodic
	 		//pass 
	 	}else {
	 		errorHandler.fatalLogicError("bounday didnt match, bid",Bnd[i].rid);
	 	}
	}
	printf("bounday summary: wall:%d, inlet:%d, outlet:%d, sym:%d\n",ridCount[1],ridCount[2],ridCount[3],ridCount[4]);

	// change grid boundary to actual boundary
	for( i=0;i<Nfac;i++ )
		RUFace[i] = 0.;

	//---------init interface-------
	dataPartition->interfaceCommunicationBegin(Un);
	dataPartition->interfaceCommunicationBegin(Vn);
	dataPartition->interfaceCommunicationBegin(Wn);
	dataPartition->interfaceCommunicationBegin(Pn);
	dataPartition->interfaceCommunicationBegin(Rn);
	dataPartition->interfaceCommunicationBegin(Tn);
	if(TurModel==1){
		dataPartition->interfaceCommunicationBegin(TE);
		dataPartition->interfaceCommunicationBegin(ED);
	}
	
	dataPartition->interfaceCommunicationBegin(deltaP);
	dataPartition->interfaceCommunicationBegin(VisLam);
	dataPartition->interfaceCommunicationBegin(VisTur);

	dataPartition->interfaceCommunicationEnd();
}


void NavierStokesSolver::SaveTransientOldData( )
{
	int i,j;
	if(      TimeScheme==1 ){  // Euler forwards
		for(i=0;i<Ncel;i++)
		{
			Rnp [i]= Rn[i];
			Unp [i]= Un[i];
			Vnp [i]= Vn[i];
			Wnp [i]= Wn[i];
			Tnp [i]= Tn[i];
			for( j=0; j<Nspecies; j++ )
				RSnp[i][j]= RSn[i][j];
			if( TurModel==1 )
			{
				TEp[i]= TE[i];
				EDp[i]= ED[i];
			}
		}
	}
	else if( TimeScheme==2 ){  // 2nd order BDF
		for(i=0;i<Ncel;i++)
		{
			//?? this may go wrong if compiler does not execute from right to left ??
			Rnp2[i]= Rnp[i]= Rn[i];
			Unp2[i]= Unp[i]= Un[i];
			Vnp2[i]= Vnp[i]= Vn[i];
			Wnp2[i]= Wnp[i]= Wn[i];
			Tnp2[i]= Tnp[i]= Tn[i];
			for( j=0; j<Nspecies; j++ )
				RSnp2[i][j]= RSnp[i][j]= RSn[i][j];
			if( TurModel==1 )
			{
				TEp2[i]= TEp[i]= TE[i];
				EDp2[i]= EDp[i]= ED[i];
			}
		}
	}
	else
	{
		errorHandler.fatalLogicError("No unsteady time advance scheme? Are you kidding?\n");
	}
}



/***********************************************/
// 	CONSTRUCTOR !!!	
/***********************************************/
NavierStokesSolver::NavierStokesSolver():
	outputCounter(0), 
	dataPartition(new DataPartition),
	root(0),			// the root rank is 0
	field(NULL),
	oldField(NULL),
	oldField2(NULL),
	iOptions(new int[INT_OPTION_NO]),
	dbOptions(new double[DB_OPTION_NO]),
	//option sets
	//bool
	IfReadBackup		(iOptions[0]),
	IfSteady			(iOptions[1]),
	SolveEnergy			(iOptions[2]),
	SolveSpecies		(iOptions[3]),
	//int
	MaxOuterStep		(iOptions[4]),
	TurModel			(iOptions[5]),
	DensityModel		(iOptions[6]),
	limiter				(iOptions[7]),
	TimeScheme			(iOptions[8]),
	noutput				(iOptions[9]),
	outputFormat		(iOptions[10]),
	Nspecies			(iOptions[11]),
	cellPressureRef		(iOptions[12]),
	MaxStep				(iOptions[13]),
	//double
	PressureReference	(dbOptions[0]),
	gama				(dbOptions[1]),
	ga1					(dbOptions[2]),
	cp					(dbOptions[3]),
	cv					(dbOptions[4]),
	prl					(dbOptions[5]),
	prte				(dbOptions[6]),
	Rcpcv				(dbOptions[7]),
	TempRef				(dbOptions[8]),
	total_time			(dbOptions[9]),
	dt					(dbOptions[10]),
	gravity				(dbOptions+11), //gravity components: 11,12,13
	URF					(dbOptions+14),  //URF 	14~21 //length 8
	ResidualSteady		(dbOptions[22]),

	//all put NULL to avoid wild pointer
	Vert(NULL),Face(NULL),Cell(NULL),Bnd(NULL),
	Rn(NULL),Un(NULL),Vn(NULL),Wn(NULL),Tn(NULL),TE(NULL),ED(NULL),
	RSn(NULL),

	Pn(NULL),
	
	Rnp(NULL),Unp(NULL),Vnp(NULL),Wnp(NULL),Tnp(NULL),TEp(NULL),EDp(NULL),RSnp(NULL),
	Rnp2(NULL),Unp2(NULL),Vnp2(NULL),Wnp2(NULL),Tnp2(NULL),TEp2(NULL),EDp2(NULL),RSnp2(NULL),

	deltaP(NULL),
	VisLam(NULL),VisTur(NULL),
	dPdX(NULL),dUdX(NULL),dVdX(NULL),dWdX(NULL),Apr(NULL),dPhidX(NULL),

	RUFace(NULL),

	BRo(NULL),BU(NULL),BV(NULL),BW(NULL),BPre(NULL),BTem(NULL),BRS(NULL),BTE(NULL),BED(NULL)
{
}


/***********************************************/
// 	DECONSTRUCTOR
/***********************************************/
NavierStokesSolver::~NavierStokesSolver()
{
	// output the result before error
	// Output2Tecplot ();
	
	delete field;
	delete oldField;
	delete oldField2;
	

	delete []iOptions;
	delete []dbOptions;
	// delete primitive variable
	delete [] Vert;
	delete [] Face;
	delete [] Cell;	
	delete [] Bnd;

		
	//these are just pointers.
	// delete variables
	/*
   	delete [] Rn;
	delete [] Un;
	delete [] Vn;
	delete [] Wn;
	delete [] Pn;
	delete [] Tn;
	delete [] TE;
	delete [] ED;
	delete_Array2D<double>(RSn,Nspecies,Ncel+dataPartition->nVirtualCell);
	

	delete [] Rnp;
	delete [] Unp;
	delete [] Vnp;
	delete [] Wnp;
	delete [] Tnp;
	delete [] TEp;
	delete [] EDp;
	delete_Array2D(RSnp,Nspecies,Ncel);
	delete [] Rnp2;
	delete [] Unp2;
	delete [] Vnp2;
	delete [] Wnp2;
	delete [] Tnp2;
	delete [] TEp2;
	delete [] EDp2;
	delete_Array2D(RSnp2,Nspecies,Ncel);
	*/

	delete [] deltaP;
	delete [] VisLam;
	delete [] VisTur;
	delete_Array2D( dPdX,Ncel,3 );
	delete_Array2D( dUdX,Ncel,3 );
	delete_Array2D( dVdX,Ncel,3 );
	delete_Array2D( dWdX,Ncel,3 );
	delete [] Apr;
	delete_Array2D( dPhidX,Ncel,3 );

	
	delete [] BRo ;
	delete [] BU  ;
	delete [] BV  ;
	delete [] BW  ;
	delete [] BTem;
	delete [] BPre;


	delete_Array2D(BRS,Nspecies,Nbnd);
	delete [] BTE;
	delete [] BED;
	delete [] RUFace;
	/*****************THIS PART IS ADDED BY CHENXUYI*******************/

	
	delete dataPartition;

	/*
        V_Destr ( &bs );
	V_Destr ( &bu );
	V_Destr ( &bv );
	V_Destr ( &bw );
	V_Destr ( &bp );
	V_Destr ( &xsol );
	*/
}

