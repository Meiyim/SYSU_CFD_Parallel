/***********************************************************
 *
 *	CYCAS 2 Parallel Power by PETSC
 *         CHEN XU YI-------2016/02/01	
 *
 * 
***********************************************************/

#undef __FUNCT__
#define __FUNCT__
static char help[256] = "CYCAS 2 Parallel Power by PETSC";
#include<stdio.h>
#include<iostream>
#include<stdexcept>
#include<petscksp.h>
#include"NS/navier.h"
#include"MPIStructure.h"
int main(int argc, char* argv[])try{
	PetscErrorCode ierr;

	ierr = PetscInitialize(&argc,&argv,(char*)0,help); CHKERRQ(ierr);

	NavierStokesSolver* nsSolver = new NavierStokesSolver;
	
	nsSolver->Init();	     	//root only : read param.in and check
	nsSolver->readAndPartition();	//root only : read msh and partition

	nsSolver->InitSolverParam(); 	//collective fetch param from root
	nsSolver->ReadGridFile();    	//read& buid CellData , BoundaryData

	ierr = PetscFinalize(); CHKERRQ(ierr);

	delete nsSolver;
	return 0;


}catch(std::logic_error& err){
	printf("!!!!!!!!!!!!!!!!Logic error occured!!!!!!!!!!!!!!!!!\n");
	std::cout<<err.what();
	printf("!!!!!!!!!!!!!!!!System will Abort!!!!!!!!!!!!!!!!!\n");
	getchar();
	PetscFinalize();
}catch(std::runtime_error& err){
	printf("!!!!!!!!!!!!!!!!run time error occured!!!!!!!!!!!!!!!!!\n");
	std::cout<<err.what();
	printf("!!!!!!!!!!!!!!!!System will Abort!!!!!!!!!!!!!!!!!\n");
	getchar();
	PetscFinalize();
}

