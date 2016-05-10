#include <iostream>
#include <fstream>
#include <math.h>
#include "navier.h"
#include "tools.h"

using namespace std;

int NavierStokesSolver::CalculatePressure( )
{

	PetscLogStagePush(3);

	SetBCVelocity( BRo,BU,BV,BW );
	CalRUFace2   ( );    	 // it doesn't need re-calculation ???
                        	 //calculated with u*   this is the m*
			
	// Build matrix
	BuildPressureMatrix(dataPartition->Ap,dataPartition->bp );

	PetscLogStagePop();
	PetscLogStagePush(4);

	
	//Solve Matrix
	// Solve pressure Poisson equation. Symmetric sparse linear equations solver
	try{
		dataPartition->solvePressureCorrection(1.e-6,1000,deltaP,DensityModel==0);
	}catch(ConvergeError& err){
		char temp[256];	
		sprintf(temp,"%s not converge in iter: %d, res %f\n",err.varname.c_str(),err.iter,err.residual);
		errorHandler.fatalRuntimeError(temp);//perhaps calculation might continue?
	}

	PetscLogStagePop();
	PetscLogStagePush(3);
   	// update other variables
	
	dataPartition->interfaceCommunicationBegin(deltaP);
	dataPartition->interfaceCommunicationEnd();

   	SetBCDeltaP( BPre,deltaP);
        Gradient   ( deltaP, BPre, dPdX );  // note the first parameter

	double deltaP_reference = 0.0;
	if(dataPartition->comRank == root.rank){//in root:
		deltaP_reference = deltaP[cellPressureRef];
	}

	MPI_Bcast(&deltaP_reference,1,MPI_DOUBLE,root.rank,dataPartition->comm);//communicate to get reference point data;

	for( int i=0; i<Ncel; i++ )//optimizeable
	{
		localRes[3] += fabs( deltaP[i] - deltaP_reference )*Cell[i].vol;
	
		// pressure (and density) correction
		if(      DensityModel==0 )
			Pn[i] +=  URF[3] * (deltaP[i]-deltaP_reference);  //under-relaxation factor: 
		else if( DensityModel==1 ){ // perfect gas
			Pn[i] +=  URF[3] *  deltaP[i];
			Rn[i] += deltaP[i]/(Rcpcv*Tn[i]);
			// Rn[i] = (Pn[i]+PressureReference)/(Rcpcv*Tn[i]);
		}

		// velocity correction
		double coef  = Cell[i].vol*Apr[i];
	        Un[i] -= coef*dPdX[i][0];
       	 	Vn[i] -= coef*dPdX[i][1];
        	Wn[i] -= coef*dPdX[i][2];
	}

	dataPartition->interfaceCommunicationBegin(Un);
	dataPartition->interfaceCommunicationBegin(Vn);
	dataPartition->interfaceCommunicationBegin(Wn);
	dataPartition->interfaceCommunicationBegin(Pn);
	dataPartition->interfaceCommunicationBegin(Rn);

	dataPartition->interfaceCommunicationEnd();

	/*
	printf("deltaP ref  %e\n",deltaP_reference);
	checkArray(deltaP,Ncel,"deltaP");
	checkArray(Rn,Ncel,"R'");
	checkArray(Pn,Ncel,"P'");
	checkArray(Un,Ncel,"U'");	
	checkArray(Vn,Ncel,"V'");	
	checkArray(Wn,Ncel,"W'");	
	*/

	SetBCVelocity(BRo, BU,BV,BW);

	// correct face normal velocity to satify the mass balance equ
	CorrectRUFace2( deltaP );
	
	PetscLogStagePop();
	return 0;
}

void NavierStokesSolver::BuildPressureMatrix(Mat& Ap, Vec& bp) //no second pressure correctio is used in this funtion .CXY e.g. equation 8.62
{
	int i,j, nj,in, cn[7], iface,bnd,rid;
	double Acn[7], roapf,lambda,lambda2,valcen,bpv,tmp,tmp2,vol, rof,Tf,RUnormal;



	for( i=0; i<Ncel; i++ )
	{
        	valcen = 0.;
		bpv    = 0.;
        	nj     = 0 ;
		// compressible gas, perfect gas
		// ?? d_ro/d_t = d_p/d_t /(RT) ?? source term and diagonal terms ,how to change?
		if( !IfSteady && DensityModel==1 )
		{
			// only Euler, or BDF 2nd can also be used ? Both, I guess
			valcen += -Cell[i].vol/( dt*Rcpcv*Tn[i] );
		}
		for( j=0; j<Cell[i].nface; j++ )
        	{
            		iface= Cell[i].face[j];
			// right hand side
			if( i==Face[iface].cell1 )
				bpv += RUFace[iface];
			else
				bpv -= RUFace[iface];

            		in   = Cell[i].cell[j]; //j-face neighbor cell
            		if( in<0  ){  		// ???????????????????
				if( DensityModel==0 )continue;
				else if( DensityModel==1){
					bnd= Face[iface].bnd;
					rid= Bnd[bnd].rid;
					RUnormal = RUFace[iface];
					if(      regionMap[rid].type1==2 ) // inlet
					{
						rof = BRo[bnd];
						Tf  = Tn[i];
						valcen  +=  CYCASMIN(RUnormal,0.) / (rof*Rcpcv*Tf);
					}
					else if( regionMap[rid].type1==3 ) // outlet
					{
						rof = Rn[i];
						Tf  = Tn[i];
						valcen  +=  CYCASMAX(RUnormal,0.) / (rof*Rcpcv*Tf);
					}
				}else{
					char temp[256];
					sprintf(temp,"unknown DensityModel found: %d\n",DensityModel);
					errorHandler.fatalLogicError(temp);
				}
			}
			else
			{
				if( i == Face[iface].cell1 ){
					lambda  = Face[iface].lambda;
					RUnormal=  RUFace[iface];
				}
				else{
					lambda  = 1.-Face[iface].lambda ;
					RUnormal= -RUFace[iface];
				}
				lambda2= 1.-lambda;
				roapf  = Rn[i]*Apr[i]*lambda + Rn[in]*Apr[in]*lambda2;
				vol    = Cell[i].vol*lambda + Cell[in].vol*lambda2;
				//-- centeral cell
				tmp    = roapf*vol*Face[iface].rlencos;
				valcen +=  -tmp;

				// incompressible flow
				//-- neighbor cell. Since it's symmetric, only consider upper triangle
				if( DensityModel==0 ){
					Acn[nj] = tmp;
					cn [nj] = Cell[in].globalIdx;
					nj ++ ;
				}
				// compressible correction
				// perfect gas, density change --> pressure change
				else if( DensityModel==1 ){
					rof = lambda*Rn[i] + lambda2*Rn[in];
					if( RUnormal>0. ){
						Tf  = Tn [i];
						tmp2     =  -RUnormal/(rof*Rcpcv*Tf); //CXY: addition to diagnal A due to compressible correction
						valcen  += tmp2;
					}

					Acn[nj] = tmp ;
					if( RUnormal<0. ){
						Tf  = Tn [in];
						Acn[nj] += -RUnormal/(rof*Rcpcv*Tf); //CXY: addition to neighbooring Coef due to compressible correction
					}
					cn [nj] = Cell[in].globalIdx;
					nj ++ ;
				}else{
					char temp[256];
					sprintf(temp,"unknown DensityModel found: %d\n",DensityModel);
					errorHandler.fatalLogicError(temp);
				}

			}
        	}

		// (i+1)th row, nj+1 non-zero values
		
		PetscInt row = Cell[i].globalIdx;
		assert(row>=0&&row<dataPartition->nGlobal);//check
		for(int ii=0; ii!=nj; ii++ )
			assert(cn[ii]>=0&&cn[ii]<dataPartition->nGlobal);//check

		MatSetValue(Ap,row,row,valcen,INSERT_VALUES);// diagonal
		MatSetValues(Ap,1,&row,nj,cn,Acn,INSERT_VALUES);// off-diagonal

		// right hand side
		VecSetValue(bp,row,bpv,INSERT_VALUES);

		// check
		 /*
		dataPartition->PRINT_LOG(row);
		for(int ii=0;ii!=nj;++ii){
			dataPartition->PRINT_LOG(cn[ii]);
			dataPartition->PRINT_LOG(Acn[ii]);
		}
		dataPartition->PRINT_LOG(bpv);
		*/
    	}

	//begin assembly
	MatAssemblyBegin(Ap,MAT_FINAL_ASSEMBLY);
	VecAssemblyBegin(bp);
	MatAssemblyEnd(Ap,MAT_FINAL_ASSEMBLY);
	VecAssemblyEnd(bp);

}


void NavierStokesSolver::CorrectRUFace2( double *dp )
{
	int i, c1,c2;
	int bnd;
	double dx1[3],dx2[3],lambda,lambda2,vol,roapf,coef,dp1,dp2,rop,rof;
	double ruf,rvf,rwf;

	// only for inner faces
	for( i=0; i<Nfac; i++ )
	{
		c1= Face[i].cell1;
		c2= Face[i].cell2;
		if( c2<0 ){
			bnd = Face[i].bnd;
			ruf = BU[bnd]*BRo[bnd];
			rvf = BV[bnd]*BRo[bnd];
			rwf = BW[bnd]*BRo[bnd];
			RUFace[i]=	 ruf * Face[i].n[0]
					+rvf * Face[i].n[1]
					+rwf * Face[i].n[2];
			continue;
		}


		lambda = Face[i].lambda;
		lambda2= 1.-lambda;
		vec_minus( dx1, Face[i].Xpac, Cell[c1].x, 3 );
		vec_minus( dx2, Face[i].Xnac, Cell[c2].x, 3 );

		// (ro+ro')*(u+u')= ro*u + ro*u' + ro'*u + ro'*u' (last term always omitted)
		// ro'*u (only compressible)
		if( DensityModel==1 ){
			rop = dp[c1]/(Rcpcv*Tn[c1])*lambda + dp[c2]/(Rcpcv*Tn[c2])*lambda2;
			rof = Rn[c1]*lambda                + Rn[c2]*lambda2;
			RUFace[i] += RUFace[i]*rop/rof;
		}
		// ro*u' (incoompressible or compressible)
		roapf = Rn[c1]*Apr[c1]*lambda + Rn[c2]*Apr[c2]*lambda2;
		vol   = Cell[c1].vol  *lambda + Cell[c2].vol  *lambda2;
		coef  = roapf*vol*Face[i].rlencos;
		dp1   = dp[c1] + URF[3]* vec_dot(dPdX[c1],dx1,3);  // must note this URF[3]
		dp2   = dp[c2] + URF[3]* vec_dot(dPdX[c2],dx2,3);
		RUFace[i] += -coef* ( dp2 - dp1 );
	}
	// check if sum of RUFace[i] in one cell is ZERO
	/*
	double *sum = new double[Ncel+dataPartition->nVirtualCell];
	for(int i=0;i!=Ncel+dataPartition->nVirtualCell;++i){
		sum[i] = 0.0;
	}
	
	for(int i=0;i!=Nfac;++i){
		c1 = Face[i].cell1;
		c2 = Face[i].cell2;
		if(c2<0) continue;
		sum[c1] -=  RUFace[i];
		sum[c2] += RUFace[i];
	}
	for(int i=0;i!=Ncel;++i){
		dataPartition->PRINT_LOG(sum[i]);
	}
	checkArray(sum,Ncel,"sumRUFace");
	delete []sum;
	*/
}
