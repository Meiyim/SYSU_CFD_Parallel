
#include "MPIStructure.h"

using namespace std;


int DataPartition::initPetsc(){ //collcetive

	MPI_Barrier(comm);

		
	//init PETSC vec
	ierr = VecCreateMPI(comm,nLocal,nGlobal,&bu);CHKERRQ(ierr); 
	ierr = VecDuplicate(bu,&bv);CHKERRQ(ierr);
	ierr = VecDuplicate(bu,&bw);CHKERRQ(ierr);
	ierr = VecDuplicate(bu,&bp);CHKERRQ(ierr);
	ierr = VecDuplicate(bu,&bs);CHKERRQ(ierr);

	ierr = VecSet(bu,0.0);CHKERRQ(ierr);
	ierr = VecSet(bv,0.0);CHKERRQ(ierr);
	ierr = VecSet(bw,0.0);CHKERRQ(ierr);
	ierr = VecSet(bp,0.0);CHKERRQ(ierr);
	ierr = VecSet(bs,0.0);CHKERRQ(ierr);


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

//	printf("PETSC NO. %d init complete, dimension %d x %d = %d\n",comRank,nLocal,nProcess,nGlobal);
	MPI_Barrier(comm);

	return 0;
}


int DataPartition::deinit(){
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


int DataPartition::buildInterfaceFromBuffer(int* buffer){
	size_t counter = 0;	
	size_t ninterfaces = buffer[counter++];
	for(size_t i=0;i!=ninterfaces;++i){
		int interfaceID  = buffer[counter++];
		size_t interfaceWidth = buffer[counter++];
		
		if(interfaces.find(interfaceID) == interfaces.end()){
			interfaces.insert(make_pair(interfaceID,Interface(comRank,interfaceID,comm)));
		}
		Interface* theInterface = &interfaces[interfaceID];
		for(size_t j=0;j!=interfaceWidth;++j){ //for each cell
			int recordLength = buffer[counter++];
			theInterface->sendposis.push_back(buffer[counter++]); //in_processor position //the sequence of senposis is really important
			set<int> cellIntersection;
			for(int k=0;k!=recordLength-1;++k){		      //record Length - 1 is the length of nodes list
				cellIntersection.insert( buffer[counter++]  );
			}
			theInterface->boundNodes.push_back(cellIntersection);
		}
	}


	
	printf("partition %d: ninterface: %lu\n",comRank,ninterfaces);
	for( auto& iter : interfaces ){
		printf("-->%d : width %d\n",iter.first,iter.second.getWidth());
		//printf("\n");
	}
	
	
	return 0;
}


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


	delete []offsets;
	delete []sourceCount;
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
#ifdef PETSC_SOLVE_VERBOSE
	PetscPrintf(comm,"begin GMRES velocity solve\n");
#endif

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
#ifdef PETSC_SOLVE_VERBOSE
		KSPGetIterationNumber(ksp,&iters);
		PetscPrintf(comm,"KSP GMRES - U converged in %d step! :)\n",iters);
#endif
	}
	//VecView(bu,PETSC_VIEWER_STDOUT_WORLD);//debug

#ifdef PETSC_SOLVE_VERBOSE
	double unorm;
	VecNorm(xsol,NORM_2,&unorm);
#endif

	
	ierr = VecDestroy(&xsol);CHKERRQ(ierr);
	

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
#ifdef PETSC_SOLVE_VERBOSE
		KSPGetIterationNumber(ksp,&iters);
		PetscPrintf(comm,"KSP GMRES - V converged in %d step! :)\n",iters);
#endif
	}

#ifdef PETSC_SOLVE_VERBOSE
	double vnorm;
	VecNorm(xsol,NORM_2,&vnorm);
#endif

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
#ifdef PETSC_SOLVE_VERBOSE
		KSPGetIterationNumber(ksp,&iters);
		PetscPrintf(comm,"KSP GMRES - W converged in %d step! :)\n",iters);
#endif
	}
#ifdef PETSC_SOLVE_VERBOSE
	double wnorm;
	VecNorm(xsol,NORM_2,&wnorm);
	PetscPrintf(comm,"Unorm %f\n",unorm);
	PetscPrintf(comm,"Vnorm %f\n",vnorm);
	PetscPrintf(comm,"Wnorm %f\n",wnorm);
#endif	
	ierr = VecDestroy(&xsol);CHKERRQ(ierr);

	ierr = MatZeroEntries(Au);CHKERRQ(ierr);
	


	return 0;	
}


int DataPartition::solveScarlar_GMRES(double tol,int maxIter,double const* xs){
	KSPConvergedReason reason;
	int iters;
	double residule;

	MatAssemblyBegin(As,MAT_FINAL_ASSEMBLY);
	VecAssemblyBegin(bs);
	MatAssemblyEnd(As,MAT_FINAL_ASSEMBLY);
	VecAssemblyEnd(bs);	

	MPI_Barrier(comm);
#ifdef PETSC_SOLVE_VERBOSE
	PetscPrintf(comm,"begin Scarlar Correction solve\n");
#endif


	KSPSetOperators(ksp,As,As);

	ierr = KSPSetType(ksp,KSPGMRES);CHKERRQ(ierr);

	KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);

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
	 * 	SOLVE Scarlar
	 ***************************************/
	ierr = VecCreateMPIWithArray(comm,1,nLocal,nGlobal,xs,&xsol);CHKERRQ(ierr); 
	ierr = VecAssemblyBegin(xsol);			CHKERRQ(ierr);
	
#ifdef PETSC_SOLVE_VERBOSE
	//double snorm;
	//VecNorm(bs,NORM_2,&snorm);
	//PetscPrintf(comm,"snorm %e\n",spnorm);
#endif
	ierr = KSPSolve(ksp,bs,xsol);CHKERRQ(ierr);

	KSPGetConvergedReason(ksp,&reason);

	if(reason<0){
		KSPGetIterationNumber(ksp,&iters);
		KSPGetResidualNorm(ksp,&residule);
		throw ConvergeError(iters,residule,"Tn");
	}else if(reason ==0){
		PetscPrintf(comm,"why is this program still running?\n");
	}else{
#ifdef PETSC_SOLVE_VERBOSE
		KSPGetIterationNumber(ksp,&iters);
		PetscPrintf(comm,"KSP GMRES - Pressure Correction converged in %d step! :)\n",iters);
#endif
	}

	ierr = VecDestroy(&xsol);CHKERRQ(ierr);
	ierr = MatZeroEntries(As);CHKERRQ(ierr);
	
	
	return 0;
}


int DataPartition::solvePressureCorrection(double tol, int maxIter,double const* xp,bool isSymmetric){
	KSPConvergedReason reason;
	int iters;
	double residule;

	MPI_Barrier(comm);
#ifdef PETSC_SOLVE_VERBOSE
	PetscPrintf(comm,"begin Pressure Correction solve\n");
#endif


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
	ierr = VecCreateMPIWithArray(comm,1,nLocal,nGlobal,xp,&xsol);CHKERRQ(ierr); 
	ierr = VecAssemblyBegin(xsol);			CHKERRQ(ierr);
	
#ifdef PETSC_SOLVE_VERBOSE
	double bpnorm;
	VecNorm(bp,NORM_2,&bpnorm);
	PetscPrintf(comm,"bpnorm %e\n",bpnorm);
#endif

	//ierr = VecView(bp,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

	//ierr = MatView(Ap,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

	//ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
	
	ierr = KSPSolve(ksp,bp,xsol);CHKERRQ(ierr);

	KSPGetConvergedReason(ksp,&reason);

	if(reason<0){
		KSPGetIterationNumber(ksp,&iters);
		KSPGetResidualNorm(ksp,&residule);
		throw ConvergeError(iters,residule,"Pn");
	}else if(reason ==0){
		PetscPrintf(comm,"why is this program still running?\n");
	}else{
#ifdef PETSC_SOLVE_VERBOSE
		KSPGetIterationNumber(ksp,&iters);
		PetscPrintf(comm,"KSP GMRES - Pressure Correction converged in %d step! :)\n",iters);
#endif
	}

#ifdef PETSC_SOLVE_VERBOSE
	double pnorm;
	VecNorm(xsol,NORM_2,&pnorm);
	PetscPrintf(comm,"pnorm %e\n",pnorm);
#endif
	
	ierr = VecDestroy(&xsol);CHKERRQ(ierr);
	ierr = MatZeroEntries(Ap);CHKERRQ(ierr);
	
	
	return 0;

}



/*****************************************
 *	implement of Interface
****************************************/
int Interface::send(MPI_Request* req, double* phi,int tag, const map<int,BdRegion>* rm = NULL){
	size_t width = getWidth();
	if(sendBuffer==NULL)
		sendBuffer = new double[width];

	tag+=sendTagOffset;

	for(size_t i=0;i!=width;++i){ //openMP optimizeable
		sendBuffer[i] = phi[sendposis[i]];
	}

	MPI_Issend(sendBuffer,	width, MPI_DOUBLE, 
			otherRank,
			tag, //TAG
			comm,
			req);
	return 0;

}


int Interface::recv(MPI_Request* req, double* phi,int tag){
	size_t width = getWidth();
	double* recvBuffer = phi+recvposi;

	tag+=recvTagOffset;

	MPI_Irecv(recvBuffer,width,MPI_DOUBLE,
			otherRank,
			tag,//TAG
			comm,
			req);

	return 0;
}

int Interface::send(MPI_Request* req, CellData* phi,int tag, const map<int,BdRegion>* rm = NULL){
	size_t width = getWidth();
	if(sendBufferCell==NULL)
		sendBufferCell = new CellData[width];

	tag+=sendTagOffset;

	for(size_t i=0;i!=width;++i){ //openMP optimizeable
		sendBufferCell[i] = phi[sendposis[i]];
		if(needsTranslate.find(sendposis[i])!=needsTranslate.end()  ){
			int bid = needsTranslate[sendposis[i]]; //perform translation !
			map<int,BdRegion>::const_iterator iter = rm->find(bid);
			assert(iter->second.type1 == 6);
			const BdRegion& reg = iter->second;
			sendBufferCell[i].x[0] += reg.initvalues[0];
			sendBufferCell[i].x[1] += reg.initvalues[1];
			sendBufferCell[i].x[2] += reg.initvalues[2];
		}
	}

	MPI_Issend(sendBufferCell,width, MPI_CellData,
			otherRank,
			tag, //TAG
			comm,
			req);
	return 0;

}


int Interface::recv(MPI_Request* req, CellData* phi,int tag){
	size_t width = getWidth();

	CellData* recvBufferCell = phi+recvposi;

	tag+=recvTagOffset;

	MPI_Irecv(recvBufferCell,width,MPI_CellData,
			otherRank,
			tag,//TAG
			comm,
			req);

	return 0;
}

int Interface::send(MPI_Request* req, double* phi[3],int tag, const map<int,BdRegion>* rm = NULL){
	size_t width = getWidth();
	if(sendBufferGradient==NULL)
		sendBufferGradient = new double[3*width];

	tag+=sendTagOffset;

	for(size_t i=0;i!=width;++i) //openMP optimizeable
		for(int j=0;j!=3;++j)
			sendBufferGradient[i*3+j] = phi[sendposis[i]][j];


	MPI_Issend(sendBufferGradient,	3*width, MPI_DOUBLE, 
			otherRank,
			tag, //TAG
			comm,
			req);
	return 0;

}


int Interface::recv(MPI_Request* req, double* phi[3],int tag){
	size_t width = getWidth();

	double* recvBufferGradient = &phi[recvposi][0];
	
	tag+=recvTagOffset;

	MPI_Irecv(recvBufferGradient,3*width,MPI_DOUBLE,
			otherRank,
			tag,//TAG 
			comm,
			req);

	return 0;
}


void Interface::getData(CellData* phi){
	/*
	size_t width = getWidth();
	for(int i=0;i!=width;++i){
		phi[recvposis[i]] = recvBufferCell[i];
	}
	*/

}

void Interface::getData(double* phi){
	/*
	size_t width = getWidth();
	for(int i=0;i!=width;++i)
		phi[recvposis[i]] = recvBuffer[i];
	*/
}



void Interface::getData(double* phi[3]){
	/*
	size_t width = getWidth();
	for(int i=0;i!=width;++i)
		for(int j=0;j!=3;++j)
			phi[recvposis[i]][j] = recvBufferGradient[i*3+j];
	*/
}





