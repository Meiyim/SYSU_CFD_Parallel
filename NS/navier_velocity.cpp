#include <iostream>
#include <fstream>
#include <math.h>
#include "navier.h"
#include "tools.h"

using namespace std;


int NavierStokesSolver::CalculateVelocity( )
{
//	ofstream of;
//	Q_Constr(&As,   "matrixU",   Ncel, False, Rowws, Normal, True);
	PetscLogStagePush(1);
	BuildVelocityMatrix(dataPartition->Au,dataPartition->bu,dataPartition->bv,dataPartition->bw );
	PetscLogStagePop();
	PetscLogStagePush(2);
	
	double *const _array1 = new double[Nfluid];
	double *const _array2 = new double[Nfluid];
	double *const _array3 = new double[Nfluid];
	std::copy(Un,Un+Nfluid,_array1);
	std::copy(Vn,Vn+Nfluid,_array2);
	std::copy(Wn,Wn+Nfluid,_array3);

	// solve U,V,W
	try{
		dataPartition->solveVelocity_GMRES(1.e-8,1000,Un,Vn,Wn);
	}catch(ConvergeError& err){
		char temp[256];	
		sprintf(temp,"%s not converge in iter: %d, res %f\n",err.varname.c_str(),err.iter,err.residual);
		errorHandler.fatalRuntimeError(temp);//perhaps calculation might continue?
	}
	
	PetscLogStagePop();
	PetscLogStagePush(1);

	//check	
	
	dataPartition->interfaceCommunicationBegin(Un);//update boundary
	dataPartition->interfaceCommunicationBegin(Vn);
	dataPartition->interfaceCommunicationBegin(Wn);
	dataPartition->interfaceCommunicationBegin(Apr);

	dataPartition->interfaceCommunicationEnd();
	

	for(int i=0;i!=Nfluid;++i){//optimizeable
		localRes[0] += fabs(_array1[i]-Un[i])*Cell[i].vol;
		localRes[1] += fabs(_array2[i]-Vn[i])*Cell[i].vol;
		localRes[2] += fabs(_array3[i]-Wn[i])*Cell[i].vol;
	}

	delete [] _array1;
	delete [] _array2;
	delete [] _array3;
	PetscLogStagePop();

	return 0;
}


void NavierStokesSolver::BuildVelocityMatrix(Mat& Au, Vec&bu, Vec& bv, Vec& bw)
{
	int i,j,iface, ip,in,nj,rid,bnd;
	PetscInt ani[6];
	double app,apn[6],lambda,lambda2, Visc,dxc[3], f,
		dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,
		sav1,sav2,sav3,RUnormal,ViscAreaLen, 
		su,sv,sw, ufi,vfi,wfi,ufh,vfh,wfh,fcs[3]={0.}, fde[3],fdi[3], dx[3],
		coef;


	SetBCVelocity( BRo,BU,BV,BW );
	SetBCPressure( BPre );
	Gradient ( Un, BU,   dUdX );
	Gradient ( Vn, BV,   dVdX );
	Gradient ( Wn, BW,   dWdX );
	Gradient ( Pn, BPre, dPdX );
	//debug check


	//PetscPrintf(dataPartition->comm,"begin build velocity matrix\n");
	
	/* This will lead to turbulence problem, wierd !!
	   limiter on velocity seems to affect turbulence generation */
	/*if( limiter==1 ){
		Limiter_MLP( Un, dUdX );
		Limiter_MLP( Vn, dVdX );
		Limiter_MLP( Wn, dWdX );
		Limiter_MLP( Pn, dPdX );
	}
	else if( limiter==2 ){
		Limiter_WENO( Un, dUdX );
		Limiter_WENO( Vn, dVdX );
		Limiter_WENO( Wn, dWdX );
		Limiter_WENO( Pn, dPdX );
	}*/
	Checker ck("ck in velocity");

	for( i=0; i<Nfluid; i++ )
	{
		app = 0.;
		nj  = 0 ;
		su=0.; sv=0.; sw=0.; //s1=0.;s2=0.;
		for( j=0;j<6;j++ ) apn[j] = 0.;

		// unsteady terms, here it can be optimzed by stored old data
		if( !IfSteady ){
			if(      TimeScheme==1 ){  // Euler forwards
				coef  = Rn[i]/dt*Cell[i].vol;
				app  += coef;
				su   += coef * Unp[i];
				sv   += coef * Vnp[i];
				sw   += coef * Wnp[i];
			}
			else if( TimeScheme==2 ){  // 2nd order BDF
				coef  = Rn[i]/dt*Cell[i].vol;
				app  += 1.5*coef;
				su   += coef * (2*Unp[i]-0.5*Unp2[i]);
				sv   += coef * (2*Vnp[i]-0.5*Vnp2[i]);
				sw   += coef * (2*Wnp[i]-0.5*Wnp2[i]);
			}
		}

		for( j=0;j<Cell[i].nface;j++ )
		{
			iface  = Cell[i].face[j];
			ip     = i;
			in     = Cell[i].cell[j];

			if( in<0 ) // boundary, i=ip naturally
			{
				bnd= Face[iface].bnd;
				assert(bnd>=0);
				rid= Bnd[bnd].rid;
				sav1   = Face[iface].n[0];
				sav2   = Face[iface].n[1];
				sav3   = Face[iface].n[2];
				RUnormal   = RUFace[iface];

				// diffusion boundary
				Visc   = VisLam[i]+VisTur[i]; // This Should be changed using boundary condition. e.g., BVisTur[bnd]
				dudx   = dUdX[i][0];
				dudy   = dUdX[i][1];
				dudz   = dUdX[i][2];
				dvdx   = dVdX[i][0];
				dvdy   = dVdX[i][1];
				dvdz   = dVdX[i][2];
				dwdx   = dWdX[i][0];
				dwdy   = dWdX[i][1];
				dwdz   = dWdX[i][2];
				vec_minus( dxc, Face[iface].x, Cell[i].x, 3 );
				switch( regionMap[rid].type1 ){
				case(1):
					ViscAreaLen = Visc*Face[iface].rlencos;
					app   += ViscAreaLen;

					Bnd[bnd].shear[0] = (dudx+dudx)*sav1 + (dudy+dvdx)*sav2 + (dudz+dwdx)*sav3;//modified by CXY:	
					Bnd[bnd].shear[1] = (dvdx+dudy)*sav1 + (dvdy+dvdy)*sav2 + (dvdz+dwdy)*sav3;
					Bnd[bnd].shear[2] = (dwdx+dudz)*sav1 + (dwdy+dvdz)*sav2 + (dwdz+dwdz)*sav3;
					Bnd[bnd].shear[0] *= Visc;
					Bnd[bnd].shear[1] *= Visc;
					Bnd[bnd].shear[2] *= Visc;
					Bnd[bnd].shear[0] += Pn[i] * sav1;
					Bnd[bnd].shear[1] += Pn[i] * sav2;
					Bnd[bnd].shear[2] += Pn[i] * sav3;

					fde[0]= ViscAreaLen*Un[i];
					fde[1]= ViscAreaLen*Vn[i];
					fde[2]= ViscAreaLen*Wn[i];
					fdi[0]= ViscAreaLen*Un[i];  // why is this part is the same? CXY should refering 8.79
					fdi[1]= ViscAreaLen*Vn[i];
					fdi[2]= ViscAreaLen*Wn[i];
					break;

				case(2):// ??????????? This still has problem ????????????????
				case(3):
				case(4):
					// diffusion to implicit, only to central cell
					ViscAreaLen = Visc*Face[iface].rlencos;
					app   += ViscAreaLen;

					fde[0] = (dudx+dudx)*sav1 + (dudy+dvdx)*sav2 + (dudz+dwdx)*sav3 ;
					fde[1] = (dvdx+dudy)*sav1 + (dvdy+dvdy)*sav2 + (dvdz+dwdy)*sav3;
					fde[2] = (dwdx+dudz)*sav1 + (dwdy+dvdz)*sav2 + (dwdz+dwdz)*sav3;
					fde[0] *= Visc;
					fde[1] *= Visc;
					fde[2] *= Visc;
					fdi[0] = ViscAreaLen*( dudx*dxc[0] + dudy*dxc[1] + dudz*dxc[2] - BU[bnd] );
					fdi[1] = ViscAreaLen*( dvdx*dxc[0] + dvdy*dxc[1] + dvdz*dxc[2] - BV[bnd] );
					fdi[2] = ViscAreaLen*( dwdx*dxc[0] + dwdy*dxc[1] + dwdz*dxc[2] - BW[bnd] );
					break;
				default:
					ErrorStop("no such rid");
				}
				

				// convection boundary
				switch( regionMap[rid].type1 ){
				case(1):     //---- Wall ----
					// convection to implicit, nothing
					fcs[0] = 0.;
					fcs[1] = 0.;
					fcs[2] = 0.;
					break;
				case(2):     //---- Inlet ----
					// convection to implicit
					if( RUnormal>0 ){
						cout<<"warning!! reverse flow get out of inlet. It may stop!"<<endl;
						// exit(0);
					}
					f      = CYCASMIN( RUnormal , 0.0 );
					app   -= f;
					fcs[0] = f*BU[bnd];
					fcs[1] = f*BV[bnd];
					fcs[2] = f*BW[bnd];
					break;
				case(3):     //---- Outlet ----??????????? boundary equal inner cell, so ???
					if( RUnormal<0 ){
						// cout<<"warning reverse flow get in of outlet. It may stop!"<<endl;
						// exit(0);
					}
					f      = CYCASMIN( RUnormal , 0.0 );
					app   -= f;
					fcs[0] = f*BU[bnd];
					fcs[1] = f*BV[bnd];
					fcs[2] = f*BW[bnd];
					break;
				case(4):     //---- Symmetric ----
					// RUnormal = 0.
					fcs[0] = 0.;
					fcs[1] = 0.;
					fcs[2] = 0.;
					break;
				default:
					char temp[256];
					sprintf(temp,"unknown rid: %d\n in BuildVelocityMatrix",rid);
					errorHandler.fatalRuntimeError(temp);
				}

				
			}
			else // inner cell
			{
				// force i as face-left-cell, in as face-right-cell
				if( i != ip ){
					in  = ip;
					sav1    = -Face[iface].n[0];
					sav2    = -Face[iface].n[1];
					sav3    = -Face[iface].n[2];
					RUnormal= -RUFace[iface];
					lambda  = 1.- Face[iface].lambda;
				}
				else{
					sav1    = Face[iface].n[0];
					sav2    = Face[iface].n[1];
					sav3    = Face[iface].n[2];
					RUnormal= RUFace[iface];
					lambda  = Face[iface].lambda;
				}
				
				lambda2= 1.- lambda;
				Visc= lambda*(VisLam[i]+VisTur[i]) + lambda2*(VisLam[in]+VisTur[in]);
				dudx= lambda*dUdX[i][0] + lambda2*dUdX[in][0];
				dudy= lambda*dUdX[i][1] + lambda2*dUdX[in][1];
				dudz= lambda*dUdX[i][2] + lambda2*dUdX[in][2];
				dvdx= lambda*dVdX[i][0] + lambda2*dVdX[in][0];
				dvdy= lambda*dVdX[i][1] + lambda2*dVdX[in][1];
				dvdz= lambda*dVdX[i][2] + lambda2*dVdX[in][2];
				dwdx= lambda*dWdX[i][0] + lambda2*dWdX[in][0];
				dwdy= lambda*dWdX[i][1] + lambda2*dWdX[in][1];
				dwdz= lambda*dWdX[i][2] + lambda2*dWdX[in][2];

				// convection to implicit
				if( RUnormal<0. ){
					apn[nj] += RUnormal;
					app     -= RUnormal;
				}

				// diffusion to implicit
				ViscAreaLen =  Visc*Face[iface].rlencos;
				app        +=  ViscAreaLen;
				apn[nj]    += -ViscAreaLen;
				ani[nj]     =  Cell[in].globalIdx;// in is the idx of virtual cell, cellglobalidx[in] is the global idx for matrix build
				nj ++ ;
				ck.check(Face[iface].rlencos);

				// convection to source term. (high order schemes)
				if( RUnormal>0. )
				{
					// low order interpolation, 1st upwind
					ufi= Un[i];
					vfi= Vn[i];
					wfi= Wn[i];
					// high order interpolation, e.g., 2nd upwind
					vec_minus( dx, Face[iface].x, Cell[i].x, 3);
					ufh= Un[i] + vec_dot( dUdX[i], dx, 3 );
					vfh= Vn[i] + vec_dot( dVdX[i], dx, 3 );
					wfh= Wn[i] + vec_dot( dWdX[i], dx, 3 );
				}
				else
				{
					// low order interpolation, 1st upwind
					ufi= Un[in];
					vfi= Vn[in];
					wfi= Wn[in];
					// high order interpolation, e.g., 2nd upwind
					vec_minus( dx, Face[iface].x, Cell[in].x, 3);
					ufh= Un[in] + vec_dot( dUdX[in], dx, 3 );
					vfh= Vn[in] + vec_dot( dVdX[in], dx, 3 );
					wfh= Wn[in] + vec_dot( dWdX[in], dx, 3 );
				}
				fcs[0]= RUnormal*(ufh-ufi);
				fcs[1]= RUnormal*(vfh-vfi);
				fcs[2]= RUnormal*(wfh-wfi);
				// first order upwind
				/* fcs[0] = 0.;
				fcs[1] = 0.;
				fcs[2] = 0.; */
				
				// diffusion to source term. ( compressible or incompressible )
				fde[0] = (dudx+dudx)*sav1 + (dudy+dvdx)*sav2 + (dudz+dwdx)*sav3 ; 
				// dudx, dudy, dudz is the diffusion term in the generic equation,
                                //dudx,dvdx,dwdx will vanish if Visc & Rou is constant. 
			 	//Is this part considered to be the Compressible part?
                                // why didnt see the divergence of U part.
				fde[1] = (dvdx+dudy)*sav1 + (dvdy+dvdy)*sav2 + (dvdz+dwdy)*sav3;
				fde[2] = (dwdx+dudz)*sav1 + (dwdy+dvdz)*sav2 + (dwdz+dwdz)*sav3;
				fde[0] *= Visc;
				fde[1] *= Visc;
				fde[2] *= Visc;
				// fdi[0] = ViscAreaLen*( Un[in] - Un[i] );
				// fdi[1] = ViscAreaLen*( Vn[in] - Vn[i] );
				// fdi[2] = ViscAreaLen*( Wn[in] - Wn[i] );
				vec_minus( dxc, Cell[in].x, Cell[i].x, 3 );
				fdi[0] = ViscAreaLen*( dudx*dxc[0]+dudy*dxc[1]+dudz*dxc[2] );
				fdi[1] = ViscAreaLen*( dvdx*dxc[0]+dvdy*dxc[1]+dvdz*dxc[2] );
				fdi[2] = ViscAreaLen*( dwdx*dxc[0]+dwdy*dxc[1]+dwdz*dxc[2] );

				//dataPartition->PRINT_LOG(Face[iface].rlencos);
				//dataPartition->PRINT_LOG(fde[0]);
				//dataPartition->PRINT_LOG(fdi[0]);
				//dataPartition->PRINT_LOG(fcs[0]);

			}

			su += fde[0]-fdi[0] - fcs[0]; // convective and diffusion part to source term aka right hand side of the LA
			sv += fde[1]-fdi[1] - fcs[1];
			sw += fde[2]-fdi[2] - fcs[2];
		}

		// central cell coef is stored for later use
		// app   += Rn[i]/dt*Cell[i].vol;

		app   /= URF[0];  // relaxation

		// SIMPLE method
		Apr[i] = 1./ (app);  // *URF[0],  ???? relaxation factor put with (URF[0]*app)  to recover ???      URF = under relaxation factor
		// SIMPLEC method
		/*  App[i] = app;
		for( j=0; j<nj; j++ )
			App[i] += apn[j];
		App[i] = 1./App[i];  */



		PetscInt row = Cell[i].globalIdx;
		assert(row>=0&&row<dataPartition->nGlobal);
		for(int ii=0;ii!=nj;++ii){
			assert(ani[ii]>=0&&ani[ii]<dataPartition->nGlobal);
		}
		MatSetValue(Au,row,row,app,INSERT_VALUES);	  //center cell
		MatSetValues(Au,1,&row,nj,ani,apn,INSERT_VALUES); //off-diagonal


		// right hand side, including
		//   pressure, gravity and part of convection&diffusion terms (explicit - implicit), 
		//   relaxation terms

		VecSetValue(bu,row,su + (1.-URF[0])*app*Un[i] + ( -dPdX[i][0] + Rn[i]*gravity[0] )*Cell[i].vol , INSERT_VALUES);
		VecSetValue(bv,row,sv + (1.-URF[0])*app*Vn[i] + ( -dPdX[i][1] + Rn[i]*gravity[1] )*Cell[i].vol , INSERT_VALUES);
		VecSetValue(bw,row,sw + (1.-URF[0])*app*Wn[i] + ( -dPdX[i][2] + Rn[i]*gravity[2] )*Cell[i].vol , INSERT_VALUES);


	}

	ck.report();
	MatAssemblyBegin(Au,MAT_FINAL_ASSEMBLY);
	VecAssemblyBegin(bu);
	VecAssemblyBegin(bv);
	VecAssemblyBegin(bw);

	return;
	

}

void NavierStokesSolver::CalRUFace( )
{
	int i, bnd, c1,c2;
	double ruf,rvf,rwf,lambda,lambda2;
	for( i=0; i<NfluidFac; i++ )
	{
		c1 = Face[i].cell1;
		c2 = Face[i].cell2;
		assert(Face[i].bnd!=INNER_FACE_BOUNDARY_SOLID);
		if( Face[i].bnd>=0 )
		{
			bnd= Face[i].bnd;
			ruf= BU[bnd]*BRo[bnd];
			rvf= BV[bnd]*BRo[bnd];
			rwf= BW[bnd]*BRo[bnd];
		}
		else
		{
			lambda = Face[i].lambda;
			lambda2= 1.-lambda;
			ruf= lambda*Rn[c1]*Un[c1] + lambda2*Rn[c2]*Un[c2];
			rvf= lambda*Rn[c1]*Vn[c1] + lambda2*Rn[c2]*Vn[c2];
			rwf= lambda*Rn[c1]*Wn[c1] + lambda2*Rn[c2]*Wn[c2];
		}
		RUFace[i] = ruf*Face[i].n[0] + rvf*Face[i].n[1] + rwf*Face[i].n[2];
	}
}

//--------------------------------------------
// Rhie-Chow momentum interpolation
//-------------------------------------
void NavierStokesSolver::CalRUFace2( )
{
	int i, bnd, c1,c2;
	double rf,ruf,rvf,rwf,lambda,lambda2, /*sav1n,sav2n,sav3n*/dpx,dpy,dpz,dpn,
		dx[3],P1,P2,d12,aprf,vol;
	for( i=0; i<NfluidFac; i++ )
	{
		c1 = Face[i].cell1;
		c2 = Face[i].cell2;
		assert(Face[i].bnd!=INNER_FACE_BOUNDARY_SOLID);
		if((bnd=Face[i].bnd)>=0)
		{
			ruf= BU[bnd]*BRo[bnd];
			rvf= BV[bnd]*BRo[bnd];
			rwf= BW[bnd]*BRo[bnd];
			RUFace[i] = ruf*Face[i].n[0] + rvf*Face[i].n[1] + rwf*Face[i].n[2];
		}
		else
		{
			lambda = Face[i].lambda;
			lambda2= 1.-lambda;
			rf = lambda*Rn[c1]        + lambda2*Rn[c2];
			ruf= lambda*Rn[c1]*Un[c1] + lambda2*Rn[c2]*Un[c2];
			rvf= lambda*Rn[c1]*Vn[c1] + lambda2*Rn[c2]*Vn[c2];
			rwf= lambda*Rn[c1]*Wn[c1] + lambda2*Rn[c2]*Wn[c2];
			RUFace[i] = ruf*Face[i].n[0] + rvf*Face[i].n[1] + rwf*Face[i].n[2];
			
			// edge pressure gradient
			vec_minus( dx, Face[i].Xpac, Cell[c1].x, 3 );
			P1  = Pn[c1] + vec_dot( dPdX[c1], dx, 3);
			vec_minus( dx, Face[i].Xnac, Cell[c2].x, 3 );
			P2  = Pn[c2] + vec_dot( dPdX[c2], dx, 3);
			vec_minus( dx, Face[i].Xnac, Face[i].Xpac, 3);
			d12  = vec_len( dx, 3 );
			aprf = lambda*Apr [c1]     + lambda2* Apr[c2];
			vol  = lambda*Cell[c1].vol + lambda2*Cell[c2].vol;

			// interpolated pressure gradient
			//sav1n = Face[i].n[0]/Face[i].area;
			//sav2n = Face[i].n[1]/Face[i].area;
			//sav3n = Face[i].n[2]/Face[i].area;
			dpx= lambda*dPdX[c1][0] + lambda2*dPdX[c2][0];
			dpy= lambda*dPdX[c1][1] + lambda2*dPdX[c2][1];
			dpz= lambda*dPdX[c1][2] + lambda2*dPdX[c2][2];
			// dpn= dpx*sav1n + dpy*sav2n + dpz*sav3n;
			dpn= dpx*dx[0] + dpy*dx[1] + dpz*dx[2];

			//correction
			// RUFace[i] -= rf*aprf*Face[i].area*vol *( (P2-P1)/d12 - dpn );  // check the unit on right hand side
			RUFace[i] -= rf*Face[i].area *aprf/d12*vol *( (P2-P1) - dpn );
		}
	}

/*
	Checker ck("chkecker in cal ruface");
	for(i=0;i!=Ncel;++i){
		for(int j=0;j!=Cell[i].nface;++j){
			int iface = Cell[i].face[j];
			if(Face[iface].cell2 < 0)continue;
			if(Face[iface].cell1 == i ) 
				ck.check(RUFace[iface]);
			if(Face[iface].cell2 == i ) 
				ck.check(-RUFace[iface]);

		}
	}
	ck.report();
	*/
}
