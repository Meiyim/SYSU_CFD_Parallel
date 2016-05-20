#include <math.h>
#include <iostream>
#include <fstream>
#include "navier.h"

using namespace std;

void NavierStokesSolver::SetBCVelocity( double *br, double *bu,double *bv,double *bw )
{
	int    i,rid,iface,ic;
	double unormal,sav1n,sav2n,sav3n;

	for( i=0; i<Nbnd; i++ )
	{
		rid   = Bnd[i].rid;
		iface = Bnd[i].face;
		ic    = Face[iface].cell1;
		switch( regionMap[rid].type1 ){
		case(1):  // wall
			br[i]= Rn[ic];
			bu[i]= 0.;
			bv[i]= 0.;
			bw[i]= 0.;
			break;
		case(2):  // inlet
			//remain initialvalue;
			/*
			double* initvalues = regionMap[rid].initvalues;
			br[i]= initvalues[4];
			bu[i]= initvalues[0];
			bv[i]= initvalues[1];
			bw[i]= initvalues[2];
			*/

			// if( DensityModel==1 ) br[i]= PressureReference/(Rcpcv*298.);
			break;
		case(3):  // outlet
			br[i]= Rn[ic];
			bu[i]= Un[ic];
			bv[i]= Vn[ic];
			bw[i]= Wn[ic];
			break;
		case(4):  // symmetric
			sav1n= Face[iface].n[0]/Face[iface].area;
			sav2n= Face[iface].n[1]/Face[iface].area;
			sav3n= Face[iface].n[2]/Face[iface].area;
			unormal = Un[ic]*sav1n + Vn[ic]*sav2n + Wn[ic]*sav3n;
			bu[i]= Un[ic] - unormal * sav1n;
			bv[i]= Vn[ic] - unormal * sav2n;
			bw[i]= Wn[ic] - unormal * sav3n;
			br[i]= Rn[ic];
			break;
		case(6)://periodic
			//pass
			break;

		default:
			errorHandler.fatalLogicError("no such boundary type in velocity bc:",regionMap[rid].type1);
		}
	}

	// correct outlet boundary condition
	if( DensityModel==0 ){
		double localSum[3] = {0.0, 0.0, 0.0};
		double globalSum[3] = {0.0, 0.0, 0.0};
		double &massflowin=localSum[0], &massflowout=localSum[1], &areaout=localSum[2];
		double &massflowinGlobal = globalSum[0], &massflowoutGlobal=globalSum[1],&areaoutGlobal=globalSum[2];
		double rate;
		for( i=0; i<Nbnd; i++ )
		{
			rid   = Bnd[i].rid;
			iface = Bnd[i].face;
			if( regionMap[rid].type1==2 )
				massflowin  += br[i]*( 	bu[i]*Face[iface].n[0] +
						       bv[i]*Face[iface].n[1] +
							   bw[i]*Face[iface].n[2] );
			else if( regionMap[rid].type1==3 ){
				massflowout += br[i]*( bu[i]*Face[iface].n[0] +
						       bv[i]*Face[iface].n[1] +
			                   bw[i]*Face[iface].n[2] );
				areaout += br[i]*Face[i].area;
			}
		}
		//MPI communication
		MPI_Allreduce(localSum,globalSum,3,MPI_DOUBLE,MPI_SUM,dataPartition->comm);	
		//
		if( fabs(massflowoutGlobal)>SMALL ){
			rate = - massflowinGlobal / massflowoutGlobal;
			for( i=0; i<Nbnd; i++ )
			{
				rid   = Bnd[i].rid;
				iface = Bnd[i].face;
				ic    = Face[iface].cell1;
				if( regionMap[rid].type1==3 )
				{
					BU[iface] *= rate ;
					BV[iface] *= rate ;
					BW[iface] *= rate ;
				}
			}
		}
		else
		{
			rate = (-massflowinGlobal - massflowoutGlobal)/areaoutGlobal;
			for( i=0; i<Nbnd; i++ )
			{
				rid   = Bnd[i].rid;
				iface = Bnd[i].face;
				ic    = Face[iface].cell1;
				if( regionMap[rid].type1==3 )
				{
					BU[iface] += rate * Face[iface].n[0]/Face[iface].area;
					BV[iface] += rate * Face[iface].n[1]/Face[iface].area;
					BW[iface] += rate * Face[iface].n[2]/Face[iface].area;
				}
			}
		}
	}
}

void NavierStokesSolver::SetBCPressure(double*bp)
{
	int    i,rid,iface,ic;
	for( i=0; i<Nbnd; i++ )
	{
		rid   = Bnd[i].rid;
		iface = Bnd[i].face;
		ic    = Face[iface].cell1;
		switch( regionMap[rid].type1 ){
		case(1):  // wall
			bp[i] = Pn[ic];
			break;
		case(2):  // inlet
			/*
			bp[i] = pin;
			break;
			*/
		case(3):  // outlet, back step
			/*
			bp[i] = pout;
			break;
			*/
		case(4):
			bp[i] = Pn[ic];
			break;
		case(6)://periodic
			//pass
			break;
		default:
			char temp[256];
			sprintf(temp,"no such boundary type %d \n ",rid);
			errorHandler.fatalLogicError(temp);
		}
	}
}

void NavierStokesSolver::SetBCDeltaP(double*bp, double *dp)
{
	int    i,rid,iface,ic;
	for( i=0; i<Nbnd; i++ )
	{
		rid   = Bnd[i].rid;
		iface = Bnd[i].face;
		ic    = Face[iface].cell1;
		switch( regionMap[rid].type1 ){
		case(1):  // wall
		case(2):  // inlet
		case(3):
		case(4):
			bp[i] = dp[ic];
			break;
		case(6)://periodic
			//pass
			break;
		default:
			char temp[256];
			sprintf(temp,"no such boundary type %d \n ",rid);
			errorHandler.fatalLogicError(temp);
		}
	}
}

void NavierStokesSolver::SetBCTemperature( double *bt )
{
	int    i,rid,iface,ic;
	for( i=0; i<Nbnd; i++ )
	{
		rid   = Bnd[i].rid;
		iface = Bnd[i].face;
		ic    = Face[iface].cell1;
		switch( regionMap[rid].type1 ){
		case(1):  // wall
			//remain initial
			break;
		case(2):  // inlet
			//remain initial
			//printf("inlet bt is %e\n",bt[i]);
			break;
		case(3):
			bt[i]= Tn[ic];
			break;
		case(4):
			bt[i]= Tn[ic];
			break;
		case(6)://periodic
			//pass
			break;
		default:
			char temp[256];
			sprintf(temp,"no such boundary type %d \n ",rid);
			errorHandler.fatalLogicError(temp);
		}
	}
}

void NavierStokesSolver::SetBCSpecies ( double **brs )
{
}

void NavierStokesSolver::SetBCKEpsilon(double *TESource,double *EDSource,double *ApTE,double *ApED, double *Prod )
{
	int    i,rid,iface,ic;
	double Cmu25,vnor,vt1,vt2,vt3,vel,utau,yplus,tauw, eps,sav1n,sav2n,sav3n,vol;
	using namespace TurKEpsilonVar;

	/*ofstream of;
	of.open("wallcell.dat");
	for( i=0; i<Nbnd; i++ )
	{
		rid   = Bnd[i].rid;
		iface = Bnd[i].face;
		ic    = Face[iface].cell1;
		if( regionMap[rid].type1==1 ) of<<ic<<endl;
	}
	of.close();*/

	/*
	double* debugArray = new double[Nbnd];

	CHECK_ARRAY(ED,Ncel);
	for(i=0;i<Nbnd;++i){
		debugArray[i] = 0.0;
	}
	*/

	

	for( i=0; i<Nbnd; i++ )
	{
		rid   = Bnd[i].rid;
		iface = Bnd[i].face;
		ic    = Face[iface].cell1;

		if(regionMap[rid].type1==1){
		 	 // wall
			BTE[i]= TE[ic];
			BED[i]= ED[ic];

			vol = Cell[ic].vol;
			// calculate the skin friction velocity
			sav1n = Face[iface].n[0]/(Face[iface].area+1.e-16);
			sav2n = Face[iface].n[1]/(Face[iface].area+1.e-16);
			sav3n = Face[iface].n[2]/(Face[iface].area+1.e-16);
			vnor= Un[ic]*sav1n + Vn[ic]*sav2n + Wn[ic]*sav3n;
			vt1 = Un[ic] - vnor*sav1n;
			vt2 = Vn[ic] - vnor*sav2n;
			vt3 = Wn[ic] - vnor*sav3n;
			vel = sqrt( vt1*vt1 + vt2*vt2 + vt3*vt3 );

			Cmu25 = pow(Cmu,0.25);
			utau  = Cmu25*sqrt( TE[ic] );
			yplus = Bnd[i].distance * utau * Rn[ic] / (VisLam[ic]+SMALL);
			tauw  = VisLam[ic] * vel/Bnd[i].distance;  // ?? which viscosity ??
			if( yplus<11. )  // inner
				utau = sqrt(tauw/Rn[ic]);
			else             // outer
				utau = Cmu25 * sqrt(TE[ic]);
			Prod[ic]= tauw * utau / (kappa*Bnd[i].distance);
			eps     = pow(Cmu,0.75)*pow(TE[ic],1.5)/(kappa*Bnd[i].distance);
			TESource[ic]  = Prod[ic] * vol ;
			ApTE    [ic] -= Rn[ic] * ED[ic] /TE[ic] * vol;
			ApTE    [ic] += Rn[ic] * eps    /TE[ic] * vol;

			ED      [ic] = eps;
			BED     [i ] = eps;
		}
	}

	for( i=0; i<Nbnd; i++ )
	{
		rid   = Bnd[i].rid;
		iface = Bnd[i].face;
		ic    = Face[iface].cell1;

		//remain initial value
		if(regionMap[rid].type1==2){
			/*
			BTE[i]= tein ;  // turbulence intensity
			BED[i]= edin ;
			*/
		}
	}
	
	for( i=0; i<Nbnd; i++ )
	{
		rid   = Bnd[i].rid;
		iface = Bnd[i].face;
		ic    = Face[iface].cell1;

		if(regionMap[rid].type1==3){
			BTE[i]= TE[ic];
			BED[i]= ED[ic];
		}
	}

	for( i=0; i<Nbnd; i++ )
	{
		rid   = Bnd[i].rid;
		iface = Bnd[i].face;
		ic    = Face[iface].cell1;

		if(regionMap[rid].type1==4){
			BTE[i]= TE[ic];
			BED[i]= ED[ic];
		}
	}

	/*
	for(int i=0;i!=Nbnd;++i){
		debugArray[i] = B[i];
	}
	*/

	/*
	CHECK_ARRAY(debugArray,Nbnd);
	CHECK_ARRAY(ED,Ncel);
	CHECK_ARRAY(BED,Nbnd);
	delete [] debugArray;
	*/

	/******************  INTERFACE COMMUNICATION  ********************/
	dataPartition->interfaceCommunicationBegin(ED);
	dataPartition->interfaceCommunicationEnd();

}
