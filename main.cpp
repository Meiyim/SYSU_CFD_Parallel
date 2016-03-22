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

ErrorHandler errorHandler;//abort when fatal error

int main(int argc, char* argv[]){
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
		nsSolver->broadcastPartitionInfo();
		nsSolver->scatterGridFile(&elementBuffer,&vertexBuffer,&interfaceBuffer);//collective
	}else{
		//read geometry locally
	}


	//parse the gridfile as original, buffer freeed, boundInfo got;
	nsSolver->ReadGridFile(elementBuffer,vertexBuffer,interfaceBuffer);

	nsSolver->dataPartition->initPetsc();

	//build faces. same sequence as Original;
	//
	nsSolver->CreateFaces();
	nsSolver->CellFaceInfo();
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
	nsSolver->dataPartition->deinit();
	delete nsSolver;
	PetscPrintf(MPI_COMM_WORLD,"done\n");
	getchar();
	ierr = PetscFinalize(); CHKERRQ(ierr);
	return 0;


}

