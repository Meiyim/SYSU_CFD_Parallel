#include <iostream>
#include <string>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include "navier.h"
#include "tools.h"
#include "terminalPrinter.h"

using namespace std;

void NavierStokesSolver::NSSolve( )
{
	int iter;
	double ResMax,ResMax0, tstart,tend;
	tstart = ttime( );
    for( step=1; step<=MaxStep; step++ )
    {
		if( (step-1)%10==0 )printer->printSectionHead(cur_time);
		if( !IfSteady ){
			cur_time += dt;
			SaveTransientOldData( );
		}

		// outer iteration, drive residual to zero for every physical time step
		for( iter=1; iter<MaxOuterStep; iter++ ){
			CalculateVelocity ( );
				// Output2Tecplot ( );
			CalculatePressure ( );

			// scalar transportation
			//1. turbulence model
			if( TurModel==1  ) UpdateTurKEpsilon( );
			//2. energy couple
			if( SolveEnergy  ) UpdateEnergy ( );
			//3. species transport
			if( SolveSpecies ) UpdateSpecies( );
			//4. other physical models
			//    e.g., condensation, combustion, wall heat conduction

			ResMax = vec_max( Residual,10 );
			if( IfSteady ) 
				break;
			else  // unsteady
			{
				if( iter==1    ) ResMax0 = ResMax;
				if( ResMax<1.e-4 || ResMax0/(ResMax+1.e-16)>1000. ){
					printer->printStepStatus(step,iter,cur_time,dt,ResMax);
					break; // more reasonal to break : order drop 2
				}
			}
		}
		if( step%noutput==0 ){
			tend   = ttime( );
			//cout<<"  spend time "<<tend-tstart<<" seconds"<<endl;
			tstart = tend;
			//cout<<"  output to tecplot ......."<<endl;
			if( outputFormat==0 ) Output2Tecplot ( );  // exit(0);
			if( outputFormat==1 ) Output2VTK     ( );
			WriteBackupFile( );
		}
		if( IfSteady ){
			printer->printSteadyStatus(step,ResMax);
			if( ResMax<ResidualSteady )break;
		}
		else if( cur_time >= total_time )break;
	}
	printer->printEnding();
	//if( outputFormat==0 ) Output2Tecplot ( );
	//if( outputFormat==1 ) Output2VTK     ( );
	//WriteBackupFile( );
}

void NavierStokesSolver::InitSolverParam( )
{
	int i;
	// default values, can be reset in initflow
	MaxStep      = 10000 ;
	MaxOuterStep = 50  ;
	IfReadBackup = false;
	IfSteady     = true ;  ResidualSteady= 1.e-6;
	TurModel     = 0;
	DensityModel = 0;
	SolveEnergy  = false;
	SolveSpecies = false;
	Nspecies     = 0;

	PressureReference= 1.01325e5; cellPressureRef=0;
	gama			 = 1.4;
	ga1				 = 0.4;
	cp				 = 1006.;
	cv				 = cp/gama;
	prl				 = 0.72;
	prte			 = 0.9;
	Rcpcv			 = cp-cv;
	TempRef			 = 273.15;
	for( i=0; i<3; i++ )
		gravity[i] = 0.;

	//-- numerical scheme
	// relaxation factor
	URF[0]= 0.6;  // u
	URF[1]= 0.6;  // v
	URF[2]= 0.6;  // w
	URF[3]= 0.5;  // p
	URF[4]= 0.8;  // T
	URF[5]= 0.8;  // k
	URF[6]= 0.8;  // e
	URF[7]= 0.8;  // scalar
	limiter = 0;

	total_time   = 0.    ;
	dt           = 1.0   ;    // no meaning in steady case
	TimeScheme   = 1     ;

	noutput  = 1;
	outputFormat = 0;  // 0 for tecplot; 1 for vtk

	// boundary
	pin = 0.; 
	pout= 0.;

	/*
// specific problem parameters
	// cylinder
	IfSteady     = true ;  dt=0.1;  TimeScheme=1;
	SolveEnergy  = false ;
	TurModel     = 1     ;
	DensityModel = 0     ;
	*/

	/*
	// heat exchanger
	IfSteady     = true ;
	MaxStep      = 500  ;
	SolveEnergy  = true ;
	TurModel     = 1    ;
	DensityModel = 1    ; */


	// init parameters from param.in
	ReadParamFile( );

	// some parameters check work, e.g. 
}


void NavierStokesSolver::InitFlowField( )
{
	int i;

	if( IfReadBackup ) 
		ReadBackupFile( );
	else{
	for( i=0; i<Ncel; i++ )
	{
		Un[i] = initvalues[0];
		Vn[i] = initvalues[1];
		Wn[i] = initvalues[2];
		Rn[i] = initvalues[3];
		Tn[i] = initvalues[4];
		Pn[i] = 0.;
		if( DensityModel== 1 ) Rn[i]= (Pn[i]+PressureReference)/(Rcpcv*Tn[i]);

		VisLam[i]= initvalues[5]; // 0.6666667e-2;  // 1.458e-6 * pow(Tn[i],1.5) /(Tn[i]+110.5) ;
		VisTur[i]= 0.;
		if( TurModel==1 )
		{
			TE[i]    = initvalues[6];  // 1.e-4*(Un[i]*Un[i]+Vn[i]*Vn[i]+Wn[i]*Wn[i]);
			ED[i]    = initvalues[7];    // TurKEpsilonVar::Cmu * pow(TE[i],1.5) / 1.;
			VisTur[i]= Rn[i]*TurKEpsilonVar::Cmu * TE[i]*TE[i]/(ED[i]+SMALL);
		}
	}
	}
	// change grid boundary to actual boundary
	int rid;  // rtable[10]={0,1,2,3,4,0},
	for( i=0; i<Nbnd; i++ )
	{
		rid = Bnd[i].rid;
		Bnd[i].rid = rtable[ rid ];
	}
	

	for( i=0;i<Nfac;i++ )
		RUFace[i] = 0.;
}

void NavierStokesSolver::SaveTransientOldData( )
{
	int i,j;
	if(      TimeScheme==1 ){  // Euler forwards
		for(i=0;i<Ncel;i++)
		{
			Rnp [i]= Rn[i];
			Unp [i]= Un[i];
			Vnp [i]= Vn[i];
			Wnp [i]= Wn[i];
			Tnp [i]= Tn[i];
			TEp [i]= TE[i];
			EDp [i]= ED[i];
			for( j=0; j<Nspecies; j++ )
				RSnp[i][j]= RSn[i][j];
			if( TurModel==1 )
			{
				TEp[i]= TE[i];
				EDp[i]= ED[i];
			}
		}
	}
	else if( TimeScheme==2 ){  // 2nd order BDF
		for(i=0;i<Ncel;i++)
		{
			//?? this may go wrong if compiler does not execute from right to left ??
			Rnp2[i]= Rnp[i]= Rn[i];
			Unp2[i]= Unp[i]= Un[i];
			Vnp2[i]= Vnp[i]= Vn[i];
			Wnp2[i]= Wnp[i]= Wn[i];
			Tnp2[i]= Tnp[i]= Tn[i];
			TEp2[i]= TEp[i]= TE[i];
			EDp2[i]= EDp[i]= ED[i];
			for( j=0; j<Nspecies; j++ )
				RSnp2[i][j]= RSnp[i][j]= RSn[i][j];
			if( TurModel==1 )
			{
				TEp2[i]= TEp[i]= TE[i];
				EDp2[i]= EDp[i]= ED[i];
			}
		}
	}
	else
	{
		cout<<"No unsteady time advanced scheme? Are you kidding?"<<endl;
		exit(0);
	}
}

void OutArray2File(double arr[],int N, ofstream &of)
{
	for(int i=0; i<N; i++ ){
		of<<arr[i]<<"  ";
		if( i%5==0 ) of<<endl;
	}
}
void NavierStokesSolver::Output2Tecplot()
{
	int i,j;
	double *tmp,nvar;
	ofstream of;
	char tecTitle[256];
	sprintf(tecTitle,"res%04d.dat",int(this->outputCounter++));
	of.open(tecTitle);
	of<<"variables="<<"\"x\","<<"\"y\","<<"\"z\""
		<<"\"p\","<<"\"u\","<<"\"v\","<<"\"w\","<<"\"ro\","<<"\"T\","
		<<"\"Mach\","<<"\"mu\"";
	nvar = 11;
	if( TurModel==1 ){
		of<<"\"te\""<<"\"ed\"";
		nvar = 13;
	}
	of<<endl;
	of<<"zone n="<<Nvrt<<", e="<<Ncel<<", VARLOCATION=([1-3]=NODAL,[4-"<<nvar<<"]=CELLCENTERED)"
		<<"DATAPACKING=BLOCK, ZONETYPE=FEBRICK"
		<<endl;
	for( j=0; j<3; j++ ){
		for( i=0; i<Nvrt; i++ ){
			of<<Vert[i][j]<<" ";
			if( i%5==0 ) of<<endl;
		}
	}
	OutArray2File( Pn,Ncel,of );
	OutArray2File( Un,Ncel,of );
	OutArray2File( Vn,Ncel,of );
	OutArray2File( Wn,Ncel,of );
	OutArray2File( Rn,Ncel,of );
	OutArray2File( Tn,Ncel,of );
	tmp = new double[Ncel];
	for( i=0; i<Ncel; i++ ){  // Mach / velocity magnitude
		if( DensityModel==1 )
			tmp[i]= sqrt( (Un[i]*Un[i]+Vn[i]*Vn[i]+Wn[i]*Wn[i])/(gama*(Pn[i]+PressureReference)/Rn[i]) );
		else
			tmp[i]= sqrt( (Un[i]*Un[i]+Vn[i]*Vn[i]+Wn[i]*Wn[i]) ); 
	}
	OutArray2File( tmp,Ncel,of );
	for( i=0; i<Ncel; i++ ){  // mu
		tmp[i]= VisLam[i] + VisTur[i];
	}
	OutArray2File( tmp,Ncel,of );
	if( TurModel==1 ){
		OutArray2File( TE,Ncel,of );
		OutArray2File( ED,Ncel,of );
	}

	for( i=0; i<Ncel; i++ ){
		for( j=0;j<8;j++ )
			of<<Cell[i].vertices[j]+1<<" ";
		of<<endl;
	}

	of.close();
	delete []tmp;
}

void outputVTKScalar( const char name[], double arr[],int N, ofstream &of)
{
	int i;
	of<<"SCALARS "<<name<<" FLOAT"<<endl;
	of<<"LOOKUP_TABLE default"<<endl;
	for( i=0; i<N; i++ )
		of<<float(arr[i])<<endl;
}
void NavierStokesSolver::Output2VTK( )
{
	int i,j;
	double *tmp;
	ofstream of;
	of.open("res.vtk");
	of<<"# vtk DataFile Version 2.0"<<endl;
	of<<"Navier-Stokes solver result"<<endl;
	of<<"ASCII"<<endl;
	of<<"DATASET UNSTRUCTURED_GRID"<<endl;
	of<<"POINTS "<<Nvrt<<" FLOAT "<<endl;
	for( i=0; i<Nvrt; i++ ){
		for(j=0;j<3;j++) of<<Vert[i][j]<<" ";
		of<<endl;
	}
	of<<endl;

	of<<"CELLS "<<Ncel<<" "<<Ncel*9<<endl;
	for( i=0; i<Ncel; i++ ){
		of<<8<<" ";
		for(j=0;j<8;j++) of<<Cell[i].vertices[j]<<" ";
		of<<endl;
	}
	of<<"CELL_TYPES "<<Ncel<<endl;
	for( i=0; i<Ncel; i++ ){
		of<<12<<" ";
		if( (i+1)%8==0 ) of<<endl;
	}
	of<<endl;

	of<<"CELL_DATA "<<Ncel<<endl;
	
	of<<"VECTORS velocity FLOAT"<<endl;
	for( i=0; i<Ncel; i++ )
		of<<Un[i]<<" "<<Vn[i]<<" "<<Wn[i]<<endl;
	of<<endl;
	
	outputVTKScalar("P",Pn,Ncel,of);
	of<<endl;
	outputVTKScalar("ro"  ,Rn, Ncel,of);
	of<<endl;
	outputVTKScalar("T"   ,Tn, Ncel,of);
	of<<endl;

	tmp = new double[Ncel];
	for( i=0; i<Ncel; i++ ){  // Mach / velocity magnitude
		if( DensityModel==1 )
			tmp[i]= sqrt( (Un[i]*Un[i]+Vn[i]*Vn[i]+Wn[i]*Wn[i])/(gama*(Pn[i]+PressureReference)/Rn[i]) );
		else
			tmp[i]= sqrt( (Un[i]*Un[i]+Vn[i]*Vn[i]+Wn[i]*Wn[i]) ); 
	}
	outputVTKScalar("mach",tmp,Ncel,of);
	of<<endl;

	for( i=0; i<Ncel; i++ ){  // mu
		tmp[i]= VisLam[i] + VisTur[i];
	}
	outputVTKScalar("mu"  ,tmp, Ncel,of);
	of<<endl;

	if( TurModel==1 ){
	outputVTKScalar("te",TE,Ncel,of);
	of<<endl;
	outputVTKScalar("ed",ED,Ncel,of);
	of<<endl;
	}

	for( i=0; i<Ncel; i++ )
		tmp[i]= 2.0;
	outputVTKScalar("Elemtype" ,tmp, Ncel,of);
	of<<endl;

	of.close();

	delete []tmp;
}

void NavierStokesSolver::WriteBackupFile( )
{
	int i;
	ofstream of;
	cout<<"backup data..."<<endl;
	of.open("res.sav");
	for( i=0; i<Ncel; i++ )
	{
		of<<Rn[i]<<" "<<Un[i]<<" "<<Vn[i]<<" "<<Wn[i]
		  <<" "<<Pn[i]<<" "<<Tn[i]<<" "<<Rn[i]<<" "<<VisLam[i]<<endl;
		if( TurModel==1 )
			of<<VisTur[i]<<" "<<TE[i]<<" "<<ED[i]<<endl;
	}
	of.close( );
}
void NavierStokesSolver::ReadBackupFile( )
{
	int i;
	ifstream inf;
	cout<<"read data from backup..."<<endl;
	inf.open("res.sav");
	for( i=0; i<Ncel; i++ )
	{
		inf>>Rn[i]>>Un[i]>>Vn[i]>>Wn[i]
		   >>Pn[i]>>Tn[i]>>Rn[i]>>VisLam[i];
		if( TurModel==1 )
			inf>>VisTur[i]>>TE[i]>>ED[i];
	}
	inf.close( );
}

void NavierStokesSolver::OutputMoniter( )
{
    
}

NavierStokesSolver::~NavierStokesSolver()
{
	// output the result before error
	// Output2Tecplot ( );

	cout<<"desturct object and free space"<<endl;
	// delete variables
   	delete [] Rn;
	delete [] Un;
	delete [] Vn;
	delete [] Wn;
	delete [] Pn;
	delete [] Tn;
	delete [] VisLam;
	delete [] VisTur;
	delete_Array2D( dPdX,Ncel,3 );
	delete_Array2D( dUdX,Ncel,3 );
	delete_Array2D( dVdX,Ncel,3 );
	delete_Array2D( dWdX,Ncel,3 );
	delete [] Apr;
	delete_Array2D( dPhidX,Ncel,3 );

	delete [] BRo ;
	delete [] BU  ;
	delete [] BV  ;
	delete [] BW  ;
	delete [] BTem;
	delete [] BPre;

	delete [] RUFace;

	/*
        V_Destr ( &bs );
	V_Destr ( &bu );
	V_Destr ( &bv );
	V_Destr ( &bw );
	V_Destr ( &bp );
	V_Destr ( &xsol );
	*/
}

