/***********************************************************
 *
 *	CYCAS 2 Parallel Power by PETSC
 *         CHEN XU YI-------2016/02/01	
 *
 * 
***********************************************************/

#undef __FUNCT__
#define __FUNCT__
#include<stdio.h>
#include<iostream>
#include<stdexcept>
#include<petscksp.h>
#include"MPIStructure.h"
#include"NS/navier.h"
int main(int argc, char* argv[])try{
	PetscErrorCode ierr;
	PetscBool shouldReadLocal = PETSC_FALSE;

	ierr = PetscInitialize(&argc,&argv,NULL,NULL); CHKERRQ(ierr);
	PetscOptionsGetBool(NULL,"-readLocally",&shouldReadLocal,NULL);

	/******************************************
	 * START UP
	 ******************************************/

	NavierStokesSolver* nsSolver = new NavierStokesSolver;
	nsSolver->initSolverParam(); 	//root only : read param.in and check
	nsSolver->readAndPartition();	//root only : read msh and partition
	nsSolver->broadcastSolverParam(); 	//collective fetch param from root

	int* elementBuffer 	= NULL;
	double* vertexBuffer 	= NULL;
	int*  interfaceBuffer 	= NULL;

	if(shouldReadLocal == PETSC_FALSE){ //transfer geometry through MPI
		nsSolver->scatterGridFile(&elementBuffer,&vertexBuffer,&interfaceBuffer);//collective
	}else{
		//read geometry locally
	}

	map<int,set<int> >* interfaceNodes= new map<int,set<int> >;

	//parse the gridfile as original, buffer freeed, boundInfo got;
	nsSolver->ReadGridFile(elementBuffer,vertexBuffer,interfaceBuffer,interfaceNodes);

	nsSolver->dataPartition->initPetsc();

	//build faces. same sequence as Original;
	nsSolver->CreateFaces();
	nsSolver->CellFaceInfo(interfaceNodes);
	delete interfaceNodes;
	nsSolver->CheckAndAllocate();
	nsSolver->InitFlowField();
	
	/******************************************
	 * NS_solve
	 * MAIN CFD 
	 ******************************************/
	nsSolver->NSSolve();

	/******************************************
	 * Post Process
	 * VTK, Tecplot, etc.
	 ******************************************/

	MPI_Barrier(MPI_COMM_WORLD);
	PetscPrintf(MPI_COMM_WORLD,"done\n");
	getchar();
	ierr = PetscFinalize(); CHKERRQ(ierr);
	return 0;


}catch(std::logic_error& err){
	int r;
	MPI_Comm_rank(MPI_COMM_WORLD,&r);
	printf("!!!!!!!!!!!!!!!!Logic error occured in rank: %d!!!!!!!!!!!!!!!!!\n",r);
	std::cout<<err.what();
	printf("!!!!!!!!!!!!!!!!System will Abort!!!!!!!!!!!!!!!!!\n");
	getchar();
	MPI_Abort(MPI_COMM_WORLD,0);
	PetscFinalize(); //is this necessary?
}catch(std::runtime_error& err){
	int r;
	MPI_Comm_rank(MPI_COMM_WORLD,&r);
	printf("!!!!!!!!!!!!!!!!run time error occured in rank: %d!!!!!!!!!!!!!!!!!\n",r);
	std::cout<<err.what();
	printf("!!!!!!!!!!!!!!!!System will Abort!!!!!!!!!!!!!!!!!\n");
	getchar();
	MPI_Abort(MPI_COMM_WORLD,0);
	PetscFinalize();
}

