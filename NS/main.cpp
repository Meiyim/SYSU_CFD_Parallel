/***********************************************************
 *
 *	CYCAS 2 Parallel Power by PETSC
 *         CHEN XU YI-------2016/02/01	
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

	ierr = PetscInitialize(&argc,&argv,NULL,NULL); CHKERRQ(ierr);
	std::map<std::string,bool> commandline;
	observeCommand(commandline,"-readLocally");
	observeCommand(commandline,"-mshBinary");
	observeCommand(commandline,"-outBinary");
	observeCommand(commandline,"-cgns");
	if(!parseCommand(commandline)){
		ierr = PetscFinalize(); CHKERRQ(ierr);
		return 0;
	}

	/******************************************
	 * START UP
	 ******************************************/

	NavierStokesSolver* nsSolver = new NavierStokesSolver;

	nsSolver->initSolverParam(); 	//root only : read param.in and check
	nsSolver->readCommand(commandline);
	nsSolver->broadcastSolverParam(); 	//collective fetch param from root

	int* elementBuffer 	= NULL;
	double* vertexBuffer 	= NULL;
	int*  interfaceBuffer 	= NULL;
	if(!commandline["-readLocally"]){ //transfer geometry through MPI
		nsSolver->readAndPartition(getCommand(commandline,"-mshBinary"), getCommand(commandline,"-cgns"));	//root only : read msh and partition
		nsSolver->broadcastPartitionInfo();
		nsSolver->scatterGridFile(&elementBuffer,&vertexBuffer,&interfaceBuffer);//collective
	}else{
		//read geometry locally
	}


	//parse the gridfile as original, buffer freeed, boundInfo got;
	nsSolver->ReadGridFile(elementBuffer,vertexBuffer,interfaceBuffer);
	nsSolver->dataPartition->initPetsc();

	//reorder cells of different body rid (FLUID & SOLID) within each partition...
	nsSolver->reorderCell();

	//build faces. same sequence as Original;
	nsSolver->CreateFaces();
	nsSolver->CellFaceInfo();
	nsSolver->CheckAndAllocate();

	//Reorder Fluid & Solid
	nsSolver->buildCoupledFace();

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
	//PetscPrintf(MPI_COMM_WORLD,"done\n");
	PetscPrintf(MPI_COMM_WORLD,"done\n");
	MPI_Barrier(MPI_COMM_WORLD);
	ierr = PetscFinalize(); CHKERRQ(ierr);
	return 0;

}

