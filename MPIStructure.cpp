
#include "MPIStructure.h"

using namespace std;


int DataPartition::initPetsc(){ //collcetive

	MPI_Barrier(comm);

	PetscPrintf(comm,"PETSC initing\n");
		
	//init PETSC vec
	ierr = VecCreateMPI(comm,nLocal,nGlobal,&bu);CHKERRQ(ierr); 
	ierr = VecDuplicate(bu,&bv);CHKERRQ(ierr);
	ierr = VecDuplicate(bu,&bw);CHKERRQ(ierr);
	ierr = VecDuplicate(bu,&bp);CHKERRQ(ierr);
	ierr = VecDuplicate(bu,&bs);CHKERRQ(ierr);
	ierr = VecDuplicate(bu,&xdp);CHKERRQ(ierr);

	ierr = VecSet(bu,0.0);CHKERRQ(ierr);
	ierr = VecSet(bv,0.0);CHKERRQ(ierr);
	ierr = VecSet(bw,0.0);CHKERRQ(ierr);
	ierr = VecSet(bp,0.0);CHKERRQ(ierr);
	ierr = VecSet(bs,0.0);CHKERRQ(ierr);
	ierr = VecSet(xdp,0.0);CHKERRQ(ierr);


	//init PETSC mat
	ierr = MatCreate(comm,&Au);CHKERRQ(ierr);     
	ierr = MatSetSizes(Au,nLocal,nLocal,nGlobal,nGlobal);CHKERRQ(ierr);
	ierr = MatSetType(Au,MATAIJ);CHKERRQ(ierr);
	ierr = MatMPIAIJSetPreallocation(Au,MAX_LOCAL_PREALLOCATION,NULL,MAX_LOCAL_PREALLOCATION,NULL);CHKERRQ(ierr);	


	ierr = MatCreate(comm,&Ap);CHKERRQ(ierr); //symmetric or not, determined when solve
	ierr = MatSetSizes(Ap,nLocal,nLocal,nGlobal,nGlobal);CHKERRQ(ierr);
	ierr = MatSetType(Ap,MATAIJ);CHKERRQ(ierr);
	ierr = MatMPIAIJSetPreallocation(Ap,MAX_LOCAL_PREALLOCATION,NULL,MAX_LOCAL_PREALLOCATION,NULL);CHKERRQ(ierr);	


	ierr = MatCreate(comm,&As);CHKERRQ(ierr);     
	ierr = MatSetSizes(As,nLocal,nLocal,nGlobal,nGlobal);CHKERRQ(ierr);
	ierr = MatSetType(As,MATAIJ);CHKERRQ(ierr);
	ierr = MatMPIAIJSetPreallocation(As,MAX_LOCAL_PREALLOCATION,NULL,MAX_LOCAL_PREALLOCATION,NULL);CHKERRQ(ierr);	

	//init KSP context
	ierr = KSPCreate(comm,&ksp);

	printf("dataPartition PETSC NO. %d init complete, dimension %d x %d = %d\n",comRank,nLocal,nProcess,nGlobal);

	return 0;
}


int DataPartition::deinit(){
	printf("datagroup NO. %d died\n",comRank);
	ierr = VecDestroy(&bu);CHKERRQ(ierr);
	ierr = VecDestroy(&bv);CHKERRQ(ierr);
	ierr = VecDestroy(&bw);CHKERRQ(ierr);
	ierr = VecDestroy(&bp);CHKERRQ(ierr);
	ierr = VecDestroy(&bs);CHKERRQ(ierr);
	ierr = MatDestroy(&Au);CHKERRQ(ierr);
	ierr = MatDestroy(&Ap);CHKERRQ(ierr);
	ierr = MatDestroy(&As);CHKERRQ(ierr);
	ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
	delete []gridList;
	gridList=NULL;
	return 0;
};


/*****************************************************
 *	MPI ScatterV / Send routine
 *	scatter geometry data to each partition
 *****************************************************/
int DataPartition::fetchDataFrom(RootProcess& root){ //collective 
	int sourceRank=root.rank;

	int destCount=0;
	int* sourceCount = new int[nProcess];
	int* offsets = new int[nProcess];
	double* localArray = NULL; 

	MPI_Barrier(comm);

	printf("start fetching data from root\n");
	offsets[0] = 0;
	sourceCount[0] = gridList[0];
 	//*************************fetching vecotr U***************************	
	for(int i=1;i!=nProcess;++i){
		sourceCount[i] = gridList[i];
		offsets[i] = offsets[i-1] + sourceCount[i];
	}
	destCount = nLocal;
	ierr = VecGetArray(bu,&localArray);CHKERRQ(mpiErr); //fetch raw pointer of PetscVector;

	//mpiErr = MPI_Scatterv(root.rootuBuffer,sourceCount,offsets,MPI_DOUBLE,
	//		localArray,destCount,MPI_DOUBLE,sourceRank,comm); CHECK(mpiErr)
	/*
	if(comRank==root.rank)
	for(int i=0;i!=nLocal;++i)	
		localArray[i] = root.rootuBuffer[i];
	*/	

	ierr = VecRestoreArray(bu,&localArray);CHKERRQ(mpiErr);

 	//*************************fetching Matrix A***************************	


	delete offsets;
	delete sourceCount;
	printf("complete fetching data from root\n");
	return 0;
}

int DataPartition::pushDataTo(RootProcess& root){//collective reverse progress of the fetch function
	MPI_Barrier(comm);
	printf("begin pushing data to root\n");

	root.allocate(this);
	pushVectorToRoot(bu,root.rootArrayBuffer,root.rank);

	printf("complete pushing data to root\n");
	return 0;
}

//********last 2 parameter only significant at root, sending PETSC vector only!***********
//			collective call
//************************************************************************
int DataPartition::pushVectorToRoot(const Vec& petscVec,double* rootBuffer,int rootRank){
	double* sendbuf = NULL;
	int sendcount = 0;

	double* recvbuf = NULL;
	int* recvcount = NULL; // significant only at root
	int* offsets = NULL;    // significant only at root
	
	recvcount = new int[nProcess];
	offsets = new int[nProcess];
	offsets[0] = 0;
 	//*************************pushing vecotr U***************************	
	recvcount[0] = gridList[0];
	for(int i=1;i!=nProcess;++i){
		recvcount[i] = gridList[i];
		offsets[i] = offsets[i-1] + recvcount[i];
	}
	sendcount = nLocal;
	recvbuf = rootBuffer;
	ierr = VecGetArray(petscVec,&sendbuf);CHKERRQ(ierr);
	mpiErr = MPI_Gatherv( sendbuf,sendcount,MPI_DOUBLE,
			      recvbuf,recvcount,offsets,MPI_DOUBLE,
			      rootRank,comm
			    );CHECK(mpiErr)

	ierr = VecRestoreArray(petscVec,&sendbuf);CHKERRQ(ierr);


	delete recvcount;
	delete offsets;

	return 0;

}



/*
int DataPartition::buildMatrix(){ //local but should involked in each processes
	PetscInt linecounter = 0;
	int ibegin=0, iend=0;
	PetscInt iInsert=0;
	PetscInt* jInsert = new PetscInt[MAX_ROW]; // it is necesasry to prescribe the max index of a row
	PetscScalar* vInsert = new PetscScalar[MAX_ROW];

	ierr = MatGetOwnershipRange(Au,&ibegin,&iend);CHKERRQ(ierr);//get range in global index
	for(int i=0;i!=nLocal;++i){ //this loop should be optimized with local parallel, tbb , openMP, etc.
		linecounter = 0;
		for(int j=0;j!=MAX_ROW;++j){
			//matrix building!
			linecounter++;
		}
		iInsert = ibegin + i;
		ierr = MatSetValues(Au,1,&iInsert,linecounter,jInsert,vInsert,INSERT_VALUES);CHKERRQ(ierr); //build matrix line by line
	}

	ierr = MatAssemblyBegin(Au,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	//!!!!!!!!!!!!!!!!!!!!!!END ASSEMBLY IS AT SOLVE GMRES!!!!!!!!!!!!!//

	printf("process%d complete buildMatrix\n",comRank);

	delete jInsert;
	delete vInsert;
	return 0;
}
*/


int DataPartition::solveVelocity_GMRES(double tol, int maxIter,double const *xu,double const *xv,double const* xw){
	KSPConvergedReason reason;
	int iters;
	double residule;

	MPI_Barrier(comm);
	PetscPrintf(comm,"begin GMRES velocity solve\n");

	ierr = MatAssemblyEnd(Au,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	//doubtful: Au Structurally Symmetric ???
	//ierr = MatSetOption(Au,MAT_STRUCTURALLY_SYMMETRIC,PETSC_TRUE);CHKERRQ(ierr);
	
	KSPSetOperators(ksp,Au,Au);
	KSPSetType(ksp,KSPGMRES);
	KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);

	/***************************************
	 *      SET  TOLERENCE
	 ***************************************/
	KSPSetTolerances(ksp,tol,PETSC_DEFAULT,PETSC_DEFAULT,maxIter);	

	/***************************************
	 * 	ILU preconditioner:
	 ***************************************/
	//KSPGetPC(ksp,&pc);
	KSPSetFromOptions(ksp);//can override settings from command line
	KSPSetUp(ksp); //the precondition is done at this step


	/***************************************
	 * 	SOLVE U!
	 ***************************************/
	ierr = VecCreateMPIWithArray(comm,1,nLocal,nGlobal,xu,&xsol);CHKERRQ(ierr);
	ierr = VecAssemblyBegin(xsol);CHKERRQ(ierr);

	ierr = VecAssemblyEnd(bu);		      CHKERRQ(ierr);
	ierr = VecAssemblyEnd(xsol);		      CHKERRQ(ierr);
	ierr = KSPSolve(ksp,bu,xsol);		      CHKERRQ(ierr);
	
	//ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

	KSPGetConvergedReason(ksp,&reason);

	if(reason<0){
		KSPGetIterationNumber(ksp,&iters);
		KSPGetResidualNorm(ksp,&residule);
		throw ConvergeError(iters,residule,"Un");
	}else if(reason ==0){
		PetscPrintf(comm,"why is this program still running?\n");
	}else{
		KSPGetIterationNumber(ksp,&iters);
		PetscPrintf(comm,"KSP GMRES - U converged in %d step! :)\n",iters);
	}
	//VecView(bu,PETSC_VIEWER_STDOUT_WORLD);//test
	double unorm;
	double bnorm;
	VecNorm(xsol,NORM_2,&unorm);
	VecNorm(bu,NORM_2,&bnorm);

	PetscPrintf(comm,"unorm:  %e\n",unorm);
	PetscPrintf(comm,"bnorm:  %e\n",bnorm);
	
	ierr = VecDestroy(&xsol);CHKERRQ(ierr);
	return 0;
	

	/***************************************
	 * 	SOLVE V!
	 * 	no need to reset KSP context and Mat
	 ***************************************/
	ierr = VecCreateMPIWithArray(comm,1,nLocal,nGlobal,xv,&xsol);CHKERRQ(ierr); 
	ierr = VecAssemblyBegin(xsol);			CHKERRQ(ierr);

	ierr = VecAssemblyEnd(bv);			CHKERRQ(ierr);
	ierr = VecAssemblyEnd(xsol);			CHKERRQ(ierr);
	ierr = KSPSolve(ksp,bv,xsol);			CHKERRQ(ierr);
	//ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);	CHKERRQ(ierr);
	ierr = KSPGetConvergedReason(ksp,&reason);	CHKERRQ(ierr);

	if(reason<0){
		KSPGetIterationNumber(ksp,&iters);
		KSPGetResidualNorm(ksp,&residule);
		throw ConvergeError(iters,residule,"Vn");
	}else if(reason ==0){
		PetscPrintf(comm,"why is this program still running?\n");
	}else{
		KSPGetIterationNumber(ksp,&iters);
		PetscPrintf(comm,"KSP GMRES - V converged in %d step! :)\n",iters);
	}
	ierr = VecDestroy(&xsol);CHKERRQ(ierr);


	/***************************************
	 * 	SOLVE W!
	 * 	no need to reset KSP context and Mat
	 ***************************************/
	ierr = VecCreateMPIWithArray(comm,1,nLocal,nGlobal,xw,&xsol);CHKERRQ(ierr); 
	ierr = VecAssemblyBegin(xsol);			CHKERRQ(ierr);

	ierr = VecAssemblyEnd(bw);			CHKERRQ(ierr);
	ierr = VecAssemblyEnd(xsol);			CHKERRQ(ierr);
	ierr = KSPSolve(ksp,bw,xsol);			CHKERRQ(ierr);
	//ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);	CHKERRQ(ierr);
	ierr = KSPGetConvergedReason(ksp,&reason);	CHKERRQ(ierr);

	if(reason<0){
		KSPGetIterationNumber(ksp,&iters);
		KSPGetResidualNorm(ksp,&residule);
		throw ConvergeError(iters,residule,"Wn");
	}else if(reason ==0){
		PetscPrintf(comm,"why is this program still running?\n");
	}else{
		KSPGetIterationNumber(ksp,&iters);
		PetscPrintf(comm,"KSP GMRES - W converged in %d step! :)\n",iters);
	}
	ierr = VecDestroy(&xsol);CHKERRQ(ierr);

	return 0;	
}


int DataPartition::solvePressureCorrection(double tol, int maxIter,bool isSymmetric){
	KSPConvergedReason reason;
	int iters;
	double residule;

	MPI_Barrier(comm);
	PetscPrintf(comm,"begin Pressure Correction solve\n");

	ierr = MatAssemblyEnd(Ap,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

	KSPSetOperators(ksp,Ap,Ap);

	if(isSymmetric){
		ierr = MatSetOption(Ap,MAT_SYMMETRIC,PETSC_TRUE);CHKERRQ(ierr);
		ierr = KSPSetType(ksp,KSPCG);CHKERRQ(ierr);
	}else{
		ierr = KSPSetType(ksp,KSPGMRES);CHKERRQ(ierr);
	}

	KSPSetInitialGuessNonzero(ksp,PETSC_FALSE);//xdp(:)=0

	/***************************************
	 *      SET  TOLERENCE
	 ***************************************/
	KSPSetTolerances(ksp,tol,PETSC_DEFAULT,PETSC_DEFAULT,maxIter);	//absolute residule

	/***************************************
	 * 	ILU preconditioner:
	 ***************************************/
	//KSPGetPC(ksp,&pc);
	KSPSetFromOptions(ksp);//can override settings from command line
	KSPSetUp(ksp); //the precondition is done at this step


	/***************************************
	 * 	SOLVE Pressure Correction!
	 ***************************************/
	/*
	ierr = VecAssemblyBegin(xdp);  CHKERRQ(ierr);// no need to assembly xdp
	ierr = VecAssemblyEnd(xdp); CHKERRQ(ierr);
	*/
	
	ierr = VecAssemblyEnd(bp);  CHKERRQ(ierr);

	//ierr = VecView(bp,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

	ierr = MatView(Ap,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

	ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
	
	ierr = KSPSolve(ksp,bp,xdp);CHKERRQ(ierr);

	KSPGetConvergedReason(ksp,&reason);

	if(reason<0){
		KSPGetIterationNumber(ksp,&iters);
		KSPGetResidualNorm(ksp,&residule);
		throw ConvergeError(iters,residule,"Pn");
	}else if(reason ==0){
		PetscPrintf(comm,"why is this program still running?\n");
	}else{
		KSPGetIterationNumber(ksp,&iters);
		PetscPrintf(comm,"KSP GMRES - Pressure Correction converged in %d step! :)\n",iters);
	}
	return 0;

}



/*****************************************
 *	implement of Interface
****************************************/
int Interface::send(MPI_Request* req, double* phi){
	size_t width = getWidth();
	if(sendBuffer==NULL)
		sendBuffer = new double[width];

	for(int i=0;i!=width;++i) //openMP optimizeable
		sendBuffer[i] = phi[sendposis[i]];
	MPI_Issend(sendBuffer,	width, MPI_DOUBLE, 
			otherRank,
			selfRank, //TAG==self Rank
			comm,
			req);
	return 0;

}


int Interface::recv(MPI_Request* req, double* phi){
	size_t width = getWidth();
	if(recvBuffer==NULL)
		recvBuffer = new double[width];

	MPI_Irecv(recvBuffer,width,MPI_DOUBLE,
			otherRank,
			otherRank,//TAG = source
			comm,
			req);

	return 0;
}

int Interface::send(MPI_Request* req, CellData* phi){
	size_t width = getWidth();
	if(sendBufferCell==NULL)
		sendBufferCell = new CellData[width];

	for(int i=0;i!=width;++i) //openMP optimizeable
		sendBufferCell[i] = phi[sendposis[i]];

	MPI_Issend(sendBufferCell,width, MPI_CellData,
			otherRank,
			selfRank, //TAG==self Rank
			comm,
			req);
	return 0;

}


int Interface::recv(MPI_Request* req, CellData* phi){
	size_t width = getWidth();
	if(recvBufferCell==NULL)
		recvBufferCell = new CellData[width];

	MPI_Irecv(recvBufferCell,width,MPI_CellData,
			otherRank,
			otherRank,//TAG = source
			comm,
			req);

	return 0;
}

int Interface::send(MPI_Request* req, double* phi[3]){
	size_t width = getWidth();
	if(sendBufferGradient==NULL)
		sendBufferGradient = new double[3*width];

	for(int i=0;i!=width;++i) //openMP optimizeable
		for(int j=0;j!=3;++j)
			sendBufferGradient[i*3+j] = phi[sendposis[i]][j];


	MPI_Issend(sendBufferGradient,	3*width, MPI_DOUBLE, 
			otherRank,
			selfRank, //TAG==self Rank
			comm,
			req);
	return 0;

}


int Interface::recv(MPI_Request* req, double* phi[3]){
	size_t width = getWidth();
	if(recvBufferGradient==NULL)
		recvBufferGradient = new double[3*width];

	MPI_Irecv(recvBufferGradient,3*width,MPI_DOUBLE,
			otherRank,
			otherRank,//TAG = source
			comm,
			req);

	return 0;
}


void Interface::getData(CellData* phi){
	size_t width = getWidth();
	for(int i=0;i!=width;++i){
		phi[recvposis[i]] = recvBufferCell[i];
	}

}

void Interface::getData(double* phi){
	size_t width = getWidth();
	for(int i=0;i!=width;++i)
		phi[recvposis[i]] = recvBuffer[i];
}


void Interface::getData(double* phi[3]){
	size_t width = getWidth();
	for(int i=0;i!=width;++i)
		for(int j=0;j!=3;++j)
			phi[recvposis[i]][j] = recvBufferGradient[i*3+j];
}




