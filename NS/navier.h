#pragma once
#include <stdlib.h>
#include <iostream>
#include <string>
#include <fstream>
#include <algorithm>
#include "BasicType.h"
#ifndef CYCAS_DEBUG_MODE
#define NDEBUG
#endif
#include <assert.h>
#include "../MPIStructure.h"
#include "tools.h"




using std::cout;
using std::endl;

//Error class used at Krylov iteration



// geometry, face & cell data
class NavierStokesSolver
{
public:
	NavierStokesSolver();
	~NavierStokesSolver();
    
	// this zone is added by CHENXUYI
	
	//main data sets;
	size_t outputCounter;
	DataPartition* dataPartition;	
	RootProcess root;

	//option Sets
	int* iOptions; 	  //len 13       //handles for integer & bool option pool, for bool, 0 = false, 1 = true
	double* dbOptions;//len 33	//handles for double option pool

	//original bools
	/******Length - 4***********/
	int& IfReadBackup;
	int& IfSteady;
	int& SolveEnergy;
	int& SolveSpecies;

	/******Length - 10***********/
	int& MaxOuterStep;
	int& TurModel;		// TurModel=0: laminar flow; =1: k-epsilon; =2: wait to implement
	int& DensityModel;       // DensityModel=0: constant density; =1: perfect gas ro=p/(RT); =2: wait to implement
	int& limiter;
	int& TimeScheme;   	// TimeScheme=1, Euler (default) ;  =2, BDF 2nd
	int& noutput;		// output step
	int& outputFormat;	// 0 for tec; 1 for vtk
	int& Nspecies;
	int& cellPressureRef;
	int& MaxStep;      
	int     *rtable;//[100];				// preset rtables only used with debugging

	/******Length - 32***********/
	double& PressureReference; 
	double& gama;
	double& ga1;
	double& cp;
	double& cv;
	double& prl;
	double& prte;
	double& Rcpcv;
	double& TempRef;
	double& total_time;
	double& dt;
	double& uin,&vin,&win,&roin,&Tin,&tein,&edin, &Twall, &pin,&pout; //input parameters ,simple implementation
	double *gravity;    // length of 3
	double *URF; 	    //numerical scheme relaxation factor , currently length 8
	double *initvalues; // len:100
	


	// above is added by CHENXUYI


	char   GridFileName[100];

    	// --- time evolution
	int    step;
    	double cur_time, Residual[10],ResidualSteady;


    	// geometry, local, build from mesh files
    	
	int          Nvrt, Ncel, Nfac, Nbnd; // the local size of the mesh ,build from mesh files
	double       **Vert; // coordinate x,y,z
    	FaceData     *Face;
   	CellData     *Cell;
    	BoundaryData *Bnd;

    	// physical variables at cell center, all local
	double  *Rn, *Un, *Vn, *Wn,  *Pn, *Tn, *TE, *ED;     // primitive vars CXY: this vars is to be replaced by DataPartition
	double  **RSn;                             // species density
	double  *VisLam, *VisTur;
	double  **dPdX, **dUdX,**dVdX,**dWdX, *Apr, **dPhidX;
	// variables at face center
	double  *RUFace;   // RUFace[Nfac]
	// Boundary faces values
	double  *BRo,*BU,*BV,*BW,*BPre, *BTem,**BRS, *BTE,*BED;
   	
        //CXY: should upgrade to PETSC	
	// laspack library work array 
	//QMatrix As,Ap;                // As for non-sysmmetric, Ap for sysmmetric matrix
	//Vector  bs,bu,bv,bw,bp, xsol; // right-hand-side vector

	// dual time unsteady simulation backup data, p = previous
	double  *Rnp, *Unp,  *Vnp,  *Wnp,  *Tnp,  *TEp,  *EDp, **RSnp,    // Euler
		    *Rnp2,*Unp2, *Vnp2, *Wnp2, *Tnp2, *TEp2, *EDp2,**RSnp2;   // BDF 2nd

/********************************************	
 *	   METHODS
*******************************************/
//CXY:
	void readAndPartition();


 	//first 3 parameter is input buffer:
	//last 1 parameter is output boundaryinfo;
	////read geometry from a buffer,
	int  ReadGridFile    (int*,double*,int*,map<int,set<int> >*);

// Geometry

	void OutputGrid      (map<int,set<int> >*); //modified by CXY
	int  CreateFaces     ( );
	void FindFace( int, int,int,int,int, int&, int*,int** );
	int  CellFaceInfo    (map<int,set<int> >* );//modified by CXY //the return value indicates the number of virtual cell beyond Ncel
	int  CheckAndAllocate();		    //modified by CXY: input the number of virtual cells

// Init flow field
    // read solver param, material, post, everything except 
	void initSolverParam(); 		//CXY: root ONLY	
	void broadcastSolverParam();	//CXY: now a MPI BROADCAST routine
	void scatterGridFile(int** elemBuffer,double** vertexBuffer,int** interfaceBuffer);		//CXY: a MPI ScatterV routine

	void InitFlowField  ( );
	


// Fluid calculation
	void NSSolve ( );

    // velocity
	int  CalculateVelocity( );
	void BuildVelocityMatrix( Mat&, Vec&, Vec&, Vec&);
	void CalRUFace ( );
	void CalRUFace2( );
	// pressure
	int  CalculatePressure( );
	void BuildPressureMatrix( Mat&, Vec&);
	void CorrectRUFace2(double*);
	// scalar. temperature, other passive variables
	// void ScalarTranport   ( double *Phi, double *BPhi, double *DiffCoef, double *source );
	void BuildScalarMatrix(int iSca,double *Phi, double *BPhi, double *DiffCoef, double *source, double *App);
	     // iSca= 1:temperature; =2: ke; =3: ed;  = 4-larger: species concentration , CXY: the difference is at bound?
	void UpdateEnergy     ( );
	void UpdateSpecies    ( );
	void UpdateTurKEpsilon( );
	void SaveTransientOldData( );

    // gradients and divergence
	int  Gradient    ( double*, double*, double** );
    // int  Divergence  ( double*, double*, double*, double*, double*, double*, double* );
	int  Limiter_MLP  (double[],double **);
	int  Limiter_Barth(double[],double **);
	int  Limiter_WENO (double[],double **);

	// Boundary condition. one defact is arrays cost too much memory, especially 2D
	void SetBCVelocity( double*br,double*bu,double*bv,double*bw );
	void SetBCPressure( double*bp );
	void SetBCDeltaP  ( double*bp, double *dp );
	void SetBCTemperature( double *bt );
	void SetBCSpecies ( double **brs );
	void SetBCKEpsilon( double *TESource,double *EDSource,double *ApTE,double*ApED,double *Prod);

// Post process
	void OutputMoniter  ( );
	void Output2Tecplot ( );
	void Output2VTK     ( );
	void WriteBackupFile( );
	void ReadBackupFile ( );

private:
	//backup
	void writeGeometryBackup(int* ebuffer, double* vbuffer,int* ibuffer); //local binary backup of the grid
	//initiation
	void ReadParamFile   ( );
	//post process
	bool shouldPostProcess(int step);		//should output
	void writeTecZoneParallel(const string& title); //collective on MPI
};

namespace TurKEpsilonVar
{
	const double 
	kappa    = 0.419 ,
	Cmu      = 0.09  ,
	Ceps1    = 1.44  ,
	Ceps2    = 1.92  ,
	LenSc    = 0.1   ,
	Sigma_k  = 1.0   ,    //! turbulence diff. coef. factors
	Sigma_e  = 1.3   ,    //! ie. turbulent Prandtl numbers
	Sigma_s  = 0.9   ;
}

extern class NavierStokesSolver NSSolver;
