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
#include<petscksp.h>
#include"NS/navier.h"
#include"dataProcess.h"
int main(int argc, char* argv[]){
	PetscErrorCode ierr;

	ierr = PetscInitialize(&argc,&argv,(char*)0,help); CHKERRQ(ierr);
	NavierStokesSolver* nsSolver = new NavierStokesSolver;
	
	nsSolver->InitSolverParam();

	ierr = PetscFinalize(); CHKERRQ(ierr);
	delete nsSolver;
	return 0;
}
