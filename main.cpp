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
#include"NS/navier.h"
#include"MPIStructure.h"
int main(int argc, char* argv[])try{
	PetscErrorCode ierr;

	ierr = PetscInitialize(&argc,&argv,NULL,NULL); CHKERRQ(ierr);

	NavierStokesSolver* nsSolver = new NavierStokesSolver;
	//PetscPrintf(MPI_COMM_WORLD,"init\n");	
	nsSolver->Init();	     	//root only : read param.in and check

	//PetscPrintf(MPI_COMM_WORLD,"init complete\n");	
	nsSolver->readAndPartition();	//root only : read msh and partition

	//PetscPrintf(MPI_COMM_WORLD,"read\n");
	//nsSolver->InitSolverParam(); 	//collective fetch param from root
	//nsSolver->ReadGridFile();    	//collective fetch geometry from root, 
					//and buid CellData , BoundaryData like it is in the serial version
	//PetscPrintf(MPI_COMM_WORLD,"done\n");
	MPI_Barrier(MPI_COMM_WORLD);
	getchar();
	ierr = PetscFinalize(); CHKERRQ(ierr);
	return 0;


}catch(std::logic_error& err){
	printf("!!!!!!!!!!!!!!!!Logic error occured!!!!!!!!!!!!!!!!!\n");
	std::cout<<err.what();
	printf("!!!!!!!!!!!!!!!!System will Abort!!!!!!!!!!!!!!!!!\n");
	getchar();
	MPI_Abort(MPI_COMM_WORLD,0);
	PetscFinalize();
}catch(std::runtime_error& err){
	printf("!!!!!!!!!!!!!!!!run time error occured!!!!!!!!!!!!!!!!!\n");
	std::cout<<err.what();
	printf("!!!!!!!!!!!!!!!!System will Abort!!!!!!!!!!!!!!!!!\n");
	getchar();
	MPI_Abort(MPI_COMM_WORLD,0);
	PetscFinalize();
}

