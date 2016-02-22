#include "MPIStructure.h"
using namespace std;


int DataPartition::init(RootProcess& root){ //collcetive

	MPI_Barrier(comm);

	printf("process %d initing\n",comRank);
	nProcess = comSize;
	if(comRank==root.rank) // in root processes
		nProcess = root.rootgridList->size();

	mpiErr = MPI_Bcast(&nProcess, 1, MPI_INT,
			   root.rank, comm);CHECK(mpiErr);
	gridList = new int[nProcess];

	if(comRank==root.rank) // in root process;
 	for(int i=0;i!=nProcess;++i)
		gridList[i] = root.rootgridList->at(i);

	mpiErr = MPI_Bcast(gridList, nProcess, MPI_INT,
			   root.rank, comm);CHECK(mpiErr)
		
	//printf("process: %d gridlist: ",comRank);
	nGlobal = 0;
	for(int i=0;i!=nProcess;++i){
		nGlobal+=gridList[i];
	//	printf("%d",gridList[i]);
	}
	nLocal = gridList[comRank];
	//printf("\n");
	//init PETSC vec
	ierr = VecCreateMPI(comm,nLocal,nGlobal,&u);CHKERRQ(ierr); 
	ierr = VecDuplicate(u,&bu);CHKERRQ(ierr);

	ierr = VecSet(bu,0.0);CHKERRQ(ierr);
	ierr = VecSet(u,1.0);CHKERRQ(ierr);
	ierr = VecAssemblyBegin(u);CHKERRQ(ierr);
	ierr = VecAssemblyEnd(u);CHKERRQ(ierr);
	ierr = VecAssemblyBegin(bu);CHKERRQ(ierr);
	ierr = VecAssemblyEnd(bu);CHKERRQ(ierr);

	//init PETSC mat
	ierr = MatCreate(comm,&Au);PCHECK(ierr);CHKERRQ(ierr);     
	ierr = MatSetSizes(Au,nLocal,nLocal,nGlobal,nGlobal);CHKERRQ(ierr);
	ierr = MatSetType(Au,MATAIJ);CHKERRQ(ierr);
	ierr = MatMPIAIJSetPreallocation(Au,MAX_ROW,NULL,MAX_ROW,NULL);CHKERRQ(ierr);	

	//init KSP context
	ierr = KSPCreate(comm,&ksp);

	printf("datagrounp NO. %d init complete, dimension %d x %d = %d\n",comRank,nLocal,nProcess,nGlobal);

	return 0;
}


int DataPartition::deinit(){
	printf("datagroup NO. %d died\n",comRank);
	ierr = VecDestroy(&u);CHKERRQ(ierr);
	ierr = VecDestroy(&bu);CHKERRQ(ierr);
	ierr = MatDestroy(&Au);CHKERRQ(ierr);
	ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
	delete []gridList;
	gridList=NULL;
	return 0;
};


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
	pushVectorToRoot(u,root.rootArrayBuffer,root.rank);
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


int DataPartition::solveGMRES(double tol, int maxIter){
	KSPConvergedReason reason;
	int iters;
	ierr = MatAssemblyEnd(Au,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	MPI_Barrier(comm);
	
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
	KSPSetFromOptions(ksp);
	KSPSetUp(ksp); //the precondition is done at this step


	/***************************************
	 * 	SOLVE!
	 ***************************************/
	KSPSolve(ksp,bu,u);
	//KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);
	KSPGetConvergedReason(ksp,&reason);
	if(reason<0){
		printf("seems the KSP didnt converge :(\n");
		return 1;
	}else if(reason ==0){
		printf("why is this program still running?\n");
	}else{
		KSPGetIterationNumber(ksp,&iters);
		printf("KSP converged in %d step! :)\n",iters);
		return 0;
	}
	



	return 0;	
}


