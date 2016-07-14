#include <iostream>
#include <fstream>
#include <math.h>
#include "navier.h"

using namespace std;

//-----------------------------------
// k-epsilon turbulence model
//------------------------------
void NavierStokesSolver::UpdateTurKEpsilon( )
{
	int i,boundaryType,ic,iface;
	double *dudx=NULL,*dvdx=NULL,*dwdx=NULL,
		   s1,s2,s3,Dis=0, vol=0,coef=0,
		   fact=0, 
		   *Prod=NULL,*TESource=NULL,*EDSource=NULL,*VisTE=NULL,*VisED=NULL,*ApTE=NULL,*ApED=NULL;
	ofstream of;
	using namespace TurKEpsilonVar;

	Prod     = new double[Nfluid];
	TESource = new double[Nfluid];
	EDSource = new double[Nfluid];
	VisTE    = new double[Ncel+dataPartition->nVirtualCell];//note the Size
	VisED    = new double[Ncel+dataPartition->nVirtualCell];
	ApTE     = new double[Nfluid];
	ApED     = new double[Nfluid];
	vec_init( TESource,Nfluid,0. );
	vec_init( EDSource,Nfluid,0. );
	vec_init( ApTE,    Nfluid,0. );
	vec_init( ApED,    Nfluid,0. );
	// viscosity and source term
	for( i=0; i<Nfluid; i++ )
	{
		vol = Cell[i].vol;
		// turbulent kinetic viscosity
    	VisTE[i] = VisLam[i] + VisTur[i] / Sigma_k ;
    	// turbulent dissipation viscosity
    	VisED[i] = VisLam[i] + VisTur[i] / Sigma_e ;

		// source terms
		dudx = dUdX[i];
		dvdx = dVdX[i];
		dwdx = dWdX[i];
    	// rate of production of turbulent energy (eq. 9.40)
    	s1 = (dudx[0]+dudx[0])*dudx[0] + (dudx[1]+dvdx[0])*dudx[1] + (dudx[2]+dwdx[0])*dudx[2];
        s2 = (dvdx[0]+dudx[1])*dvdx[0] + (dvdx[1]+dvdx[1])*dvdx[1] + (dvdx[2]+dwdx[1])*dvdx[2];
   		s3 = (dwdx[0]+dudx[2])*dwdx[0] + (dwdx[1]+dvdx[2])*dwdx[1] + (dwdx[2]+dwdx[2])*dwdx[2];

    	Prod[i] = VisTur[i] * ( s1 + s2 + s3 ) ;
    	// dissipation
    	Dis     = Rn[i] * ED[i];
		// bouyancy production term: - Gi/(sigma_h,t rho) drho/dx
        	// Pbouy =
		TESource[i] = Prod[i] * vol ;   //   - Dis  ;
		ApTE    [i] = Dis  /TE[i] * vol;        // add to diagonal element of linear system As

		//// ED production & dissipation. similarity to TE
		// fact = ED[i]/(TE[i]+SMALL);
		// Pked = Ceps1 * fact * Prod[i] * vol ;
		// Dised= Ceps2 * fact * Dis     * vol ;
		// EDSource[i] = Pked * vol;  //   - Dised;
		// ApED    [i] = Dised/ED[i] * vol;
	}
	// boundary for k, epsilon
	//
	SetBCKEpsilon( TESource,EDSource,ApTE,ApED,Prod );//CXY: setting BC would change the ED variables, interfaceCommunication needed!

	// other source term. do not put the before the boundary, for boundary change Pk,epsilon mandatorily

	if( !IfSteady ){
		if(      TimeScheme==1 ){  // Euler forwards
			for( i=0; i<Nfluid; i++ ){
				coef        = Rn[i]/dt*Cell[i].vol;
				ApTE[i]    += coef;
				TESource[i]+= coef * TEp[i];
				ApED[i]    += coef;
				EDSource[i]+= coef * EDp[i];
			}
		}
		else if( TimeScheme==2 ){  // 2nd order BDF
			for( i=0; i<Nfluid; i++ ){
				coef        = Rn[i]/dt*Cell[i].vol;
				ApTE[i]    += 1.5*coef;
				TESource[i]+= coef * (2*TEp[i]-0.5*TEp2[i]);
				ApED[i]    += 1.5*coef;
				EDSource[i]+= coef * (2*EDp[i]-0.5*EDp2[i]);
			}
		}
	}


	//Solve turbulence kinetic energy
	// build matrix
	BuildScalarMatrix( 2,TE,BTE,VisTE,TESource, ApTE,dataPartition->As,dataPartition->bs);
	// Solve equations
	try{
		dataPartition->solveScarlar_GMRES(1.e-8,500,TE);
	}catch(ConvergeError& err){
		char temp[256];	
		sprintf(temp,"TE not converge in iter: %d, res %f\n",err.iter,err.residual);
		errorHandler.fatalRuntimeError(temp);
	}

	for(i=0;i<Nfluid;++i){
		if(TE[i]<1.e-9){
			TE[i] = 1.e-9;
		}
	}

	/*****************************
	 *	solve terbulence dissipation
	 *****************************/
	
	for( i=0; i<Nfluid; i++ )
	{
		// ED production & dissipation. similarity to TE
		fact = ED[i]/(TE[i]+SMALL)*Cell[i].vol;
		EDSource[i]  = Ceps1* fact * Prod[i];
		ApED    [i] += Ceps2* fact * Rn[i]  ;
	}


	// Solve turbulence dissipation rate
	// build matrix
	BuildScalarMatrix( 3,ED,BED,VisED,EDSource,ApED,dataPartition->As,dataPartition->bs);

	// ?? specieal treatment, force ED[i] on boundary cells to be BED[ib]  ??
	for( i=0; i<Nbnd; i++ ){
		boundaryType   = regionMap[ Bnd[i].rid ].type1 ;
		if( boundaryType==1 ){
			iface = Bnd[i].face;
			ic    = Face[iface].cell1;
			int ani[7];
			double apn[7];
			int nj = 0;
			int row = Cell[ic].globalIdx;

			ani[nj] = row;
			apn[nj] = 1.0;
			++nj;
			for(int j=0;j!=Cell[ic].nface;++j){
				int icn = Cell[ic].cell[j];
				if(icn==VOID_CELL_ON_BOUNDARY){
					continue;	
				}
				ani[nj] = Cell[ icn ].globalIdx;
				apn[nj] = 0.0;
				++nj;
			}
			MatSetValues(dataPartition->As,1,&row,nj,ani,apn,INSERT_VALUES);
			VecSetValue(dataPartition->bs,row,ED[ic],INSERT_VALUES);
		}
	}

	//---solve matrix
	try{
		dataPartition->solveScarlar_GMRES(1.e-8,500,ED);
	}catch(ConvergeError& err){
		char temp[256];	
		sprintf(temp,"ED not converge in iter: %d, res %f\n",err.iter,err.residual);
		errorHandler.fatalRuntimeError(temp);
	}

	for(i=0;i<Nfluid;++i){
		if(ED[i]<1.e-12){
			ED[i] = 1.e-12;
		}
	}

	// Calculate the turbulent viscosity
	for( i=0; i<Nfluid; i++ )
    	{
		if( ED[i]>1.e-12 )
			VisTur[i]= Cmu * Rn[i] * TE[i]*TE[i]/ (ED[i]+SMALL);
		else
			VisTur[i]= 0.;//CXY: seems meaningless 
    	}

	//*******************MPI INTERFACE COMMUNICATION************************//
	//CAUTION: no gradient/interpolation or other function that involve vertual cell before this point
	
	dataPartition->interfaceCommunicationBegin(TE);
	dataPartition->interfaceCommunicationBegin(ED);
	dataPartition->interfaceCommunicationBegin(VisLam);
	dataPartition->interfaceCommunicationBegin(VisTur);
	dataPartition->interfaceCommunicationEnd();

	delete []Prod;
	delete []TESource;
	delete []EDSource;
	delete []VisTE;
	delete []VisED;
	delete []ApTE;
	delete []ApED;

}

//------------------------------------------------------------
// Energy transport equation : 
// p_(ro*T)/p_t + div(ro*V*T) = div(romda/cp*grad(T)) + QT/cp
//  note that cp is removed from time derivative
//-------------------------------------------------------
void NavierStokesSolver::UpdateEnergy( )
{
	int i;
	double *kcond=NULL,*ESource=NULL,*ApE=NULL,coef;

//	Q_Constr(&As,   "matrixU",   Ncel, False, Rowws, Normal, True);
	kcond  = new double[Ncel + dataPartition->nVirtualCell]; //note the size;
	ESource= new double[Nfluid];
	ApE    = new double[Nfluid];
	vec_init( ESource, Nfluid, 0. );
	vec_init( ApE,     Nfluid, 0. );
	vec_init( kcond, Nfluid+dataPartition->nVirtualCell,CYCASHUGE_D);
	// prepare the diffusion coefficient and source terms
	for( i=0; i<Nfluid; i++ )
	{
		kcond[i] = cp*(VisLam[i]/prl + VisTur[i]/prte)/cp;//only for gas i surpose
	}

	// part of viscous terms
	if( DensityModel==1 ){
	/*
	int c1,c2;
	double VisL,VisT,vmul,dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,
		vlac,div,txx,tyy,tzz,txy,txz,tyz, vxg,vyg,vzg,btx,bty,btz,fvis,lambda,lambda2;
	for( i=0; i<Nfac; i++ )
	{
		c1 = Face[i].cell1;
		c2 = Face[i].cell2;
		if(c1>=Nfluid) continue;

		if( Face[i].bnd>=0 )
		{
			VisL = VisLam[c1];
			VisT = VisTur[c1];
			dudx = dUdX[c1][0];
			dudy = dUdX[c1][1];
			dudz = dUdX[c1][2];
			dvdx = dVdX[c1][0];
			dvdy = dVdX[c1][1];
			dvdz = dVdX[c1][2];
			dwdx = dWdX[c1][0];
			dwdy = dWdX[c1][1];
			dwdz = dWdX[c1][2];
			vxg  = Un[c1];
			vyg  = Vn[c1];
			vzg  = Wn[c1];
		}
		else
		{
			lambda = Face[i].lambda;
			lambda2= 1.-lambda;
			VisL = lambda*VisLam[c1]  + lambda2*VisLam[c2]  ;
			VisT = lambda*VisTur[c1]  + lambda2*VisTur[c2]  ;
			dudx = lambda*dUdX[c1][0] + lambda2*dUdX[c2][0] ;
			dudy = lambda*dUdX[c1][1] + lambda2*dUdX[c2][1] ;
			dudz = lambda*dUdX[c1][2] + lambda2*dUdX[c2][2] ;
			dvdx = lambda*dVdX[c1][0] + lambda2*dVdX[c2][0] ;
			dvdy = lambda*dVdX[c1][1] + lambda2*dVdX[c2][1] ;
			dvdz = lambda*dVdX[c1][2] + lambda2*dVdX[c2][2] ;
			dwdx = lambda*dWdX[c1][0] + lambda2*dWdX[c2][0] ;
			dwdy = lambda*dWdX[c1][1] + lambda2*dWdX[c2][1] ;
			dwdz = lambda*dWdX[c1][2] + lambda2*dWdX[c2][2] ;
			vxg  = lambda*Un[c1]      + lambda2*Un[c2];
			vyg  = lambda*Vn[c1]      + lambda2*Vn[c2];
			vzg  = lambda*Wn[c1]      + lambda2*Wn[c2];
		}

		vmul= VisL + VisT;
		vlac=  -2./3.*vmul;
		div= dudx+ dvdy+ dwdz;
		txx= 2.*vmul*dudx+ vlac*div ; //CXY: pressure should goes here ?
		tyy= 2.*vmul*dvdy+ vlac*div ;
		tzz= 2.*vmul*dwdz+ vlac*div ;
		txy= vmul*(dudy+dvdx);
		txz= vmul*(dudz+dwdx);
		tyz= vmul*(dvdz+dwdy);
		btx= vxg*txx+ vyg*txy+ vzg*txz ;
		bty= vxg*txy+ vyg*tyy+ vzg*tyz ;
		btz= vxg*txz+ vyg*tyz+ vzg*tzz ;
		fvis= ( btx*Face[i].n[0] + bty*Face[i].n[1] + btz*Face[i].n[2] ) / cp;//CXY: where is work done by pressure?
		ESource[c1] += fvis;
		if( c2>=0 )
		ESource[c2] -= fvis; //CXY: dissipation source term due to compressibility
                             //CX: why is this integral done on face? i supposed it should be done in volumn
	}
		*/
	}//end if : densityModel == 1

	// boundary
	SetBCTemperature( BTem);
	// source terms, e.g., energy release, condensation/vaporization
	if( !IfSteady ){
	if(      TimeScheme==1 ){  // Euler forwards
		for( i=0; i<Nfluid; i++ ){
			coef       = Rn[i]/dt * Cell[i].vol;
			ApE[i]    += coef;
			ESource[i]+= coef * Tnp[i];
		}
	}
	else if( TimeScheme==2 ){  // 2nd order BDF
		for( i=0; i<Nfluid; i++ ){
			coef       = Rn[i]/dt * Cell[i].vol;
			ApE[i]    += 1.5*coef;
			ESource[i]+= coef * (2*Tnp[i]-0.5*Tnp2[i]);
		}
	}
	}
// build matrix
	BuildScalarMatrix( 1, Tn,BTem,kcond,ESource,ApE,dataPartition->As,dataPartition->bs );
	// record original solve
	double* const _array = new double[Nfluid];
	std::copy(Tn,Tn+Nfluid,_array);
	// Solve equations
	try{
		dataPartition->solveScarlar_GMRES(1.e-8,500,Tn);
	}catch(ConvergeError& err){
		char temp[256];	
		sprintf(temp,"Energy not converge in iter: %d, res %f\n",err.iter,err.residual);
		errorHandler.fatalRuntimeError(temp);
	}

	//********************MPI INTERFACE COMMUNICATION*****************************//
	dataPartition->interfaceCommunicationBegin(Tn);
	dataPartition->interfaceCommunicationEnd();

	for(int i=0;i!=Nfluid;++i){
		localRes[4] += fabs(_array[i] - Tn[i])*Cell[i].vol;
	}
	delete [] _array;

	if(SolveConjungateHeat){
		SetBCSolidTemperature(dPhidX,kcond);
	}

	
	// clipping work
	delete [] kcond;
	delete [] ESource;
	delete [] ApE;
}

//------------------------------------------------------------
// species transport equation : 
// p_(ro*c)/p_t + div(ro*V*c) = div(D*grad(c)) + Qm
//    c is mass fraction, Qm could be injection,phase transition etc
//--------------------------------------------
void NavierStokesSolver::UpdateSpecies( )
{
	/*
	int i,is;
	double **DiffC, **ScSource,*ApS, coef;

	ApS     = new double[Ncel]; 
	DiffC   = new_Array2D<double>(Nspecies,Ncel);
	ScSource= new_Array2D<double>(Nspecies,Ncel);
	init_Array2D( ScSource,Nspecies,Ncel, 0. );
	// prepare the diffusion coefficient and source terms
	for( is=0; is<Nspecies; is++){
		for( i=0; i<Ncel; i++ )
		{
			DiffC[is][i] = 1.e-5; // this should be calculated in Material class
			if( ! IfSteady )
				ScSource[is][i] += 0.;
		}
	}
	// boundary
	SetBCSpecies( BRS );
	// source terms, e.g., energy release

	for( is=0; is<Nspecies; is++ )
	{
		vec_init( ApS, Ncel, 0. );

		// transient time to source term and diagonal element
		if( !IfSteady ){
			if(      TimeScheme==1 ){  // Euler forwards
				for( i=0; i<Ncel; i++ ){
					coef            = Rn[i]/dt*Cell[i].vol;
					ApS[i]         += coef;
					ScSource[is][i]+= coef * RSnp[is][i];
				}
			}
			else if( TimeScheme==2 ){  // 2nd order BDF
				for( i=0; i<Ncel; i++ ){
					coef            = Rn[i]/dt*Cell[i].vol;
					ApS[i]         += 1.5*coef;
					ScSource[is][i]+= coef * (2*RSnp[is][i]-0.5*RSnp2[is][i]);
				}
			}
		}

		//Q_Constr(&As,   "matrixU",   Ncel, False, Rowws, Normal, True);

		// build matrix
		
		//BuildScalarMatrix( 4+is, RSn[is],BRS[is],DiffC[is],ScSource[is],ApS );

		// Solve equations
		for( i=0; i<Ncel; i++ ) ;
//			xsol.Cmp[i+1]= RSn[is][i];

		for( i=0; i<Ncel; i++ );
//			RSn[is][i] = xsol.Cmp[i+1];

//		Q_Destr ( &As );
	}
		*/
}


/***************************************************
 *	A universal scarlar equation solve routine
 *	note that Diffcoef is cross-process and the size needs to be Ncel+NVirtual
 *	source and App is local. the size is Ncel;
 ***************************************************/
void NavierStokesSolver::BuildScalarMatrix( int iSca, double *Phi,double *BPhi,double *DiffCoef, double *source, double *App,Mat& As,Vec& bs )
{
	int i,j,iface, ip,in,ani[6],nj,boundaryType,bnd;
	double app,apn[6],lambda,lambda2, Visc,dxc[3],
		dphidx,dphidy,dphidz, f,
		sav1,sav2,sav3,RUnormal,ViscAreaLen, 
		fcs=0., fde,fdi, dx[3],  sphi, pfi,pfh;

	Gradient ( Phi, BPhi,  dPhidX );

	dataPartition->interfaceCommunicationBegin(DiffCoef);
	dataPartition->interfaceCommunicationEnd();

	for( i=0; i<Nfluid; i++ )
	{

		app = App[i];
		nj  = 0 ;
		sphi= source[i];
		for( j=0;j<6;j++ ) apn[j] = 0.;

		for( j=0;j<Cell[i].nface;j++ )
		{
			iface  = Cell[i].face[j];
			ip     = i;
			in     = Cell[i].cell[j];
			
			if( in<0 ) // boundary, i=ip naturally
			{
				bnd = Face[iface].bnd;
				boundaryType = regionMap[Bnd[bnd].rid].type1;
				sav1    = Face[iface].n[0];
				sav2    = Face[iface].n[1];
				sav3    = Face[iface].n[2];
				RUnormal= RUFace[iface];

				// diffusion boundary
				Visc   = DiffCoef[i]; // This Should be changed using boundary condition. e.g., BVisTur[bnd]
				dphidx = dPhidX[i][0];
				dphidy = dPhidX[i][1];
				dphidz = dPhidX[i][2];
				vec_minus( dxc, Face[iface].x, Cell[i].x, 3 );
				switch( boundaryType ){
				case(1):
				case(2):
				case(3):
				case(4):
					/*fde = 0.;
					fdi = 0.;
					break;*/
					// diffusion to implicit, only to central cell
					ViscAreaLen = Visc*Face[iface].rlencos;
					app += ViscAreaLen;
					fde  = Visc*( dphidx*sav1 + dphidy*sav2 + dphidz*sav3 );
					fdi  = ViscAreaLen*( dphidx*dxc[0] + dphidy*dxc[1] + dphidz*dxc[2] - BPhi[bnd] );
					break;
				default:
					ErrorStop("no such bid");
				}
				
				// convection boundary
				switch( boundaryType ){
				case(1):     //---- Wall ----
					// convection to implicit, nothing
					fcs = 0.;
					break;
				case(2):     //---- Inlet ----
					// convection to implicit
					if( RUnormal>0 ){
						cout<<"reverse flow get out of inlet. stop!"<<endl;
						exit(0);
					}
					f    = CYCASMIN( RUnormal , 0.0 );
					app -= f;
					fcs  = f*BPhi[bnd];
					break;
				case(3):     //---- Outlet ----??????????? boundary equal inner cell, so ???
					if( RUnormal<0 ){
						cout<<"reverse flow get in of outlet. stop!"<<endl;
						// exit(0);
					}
					f    = CYCASMIN( RUnormal , 0.0 );
					app -= f;
					fcs  = f*BPhi[bnd];
					break;
				case(4):     //---- Symmetric ----
					// RUnormal = 0.
					fcs = 0.;
					break;
				default:
					cout<<"no this type of boundary! bid="<<boundaryType<<endl;
					exit(0);
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
				Visc   = lambda*DiffCoef[i]  + lambda2*DiffCoef[in];
				dphidx = lambda*dPhidX[i][0] + lambda2*dPhidX[in][0];
				dphidy = lambda*dPhidX[i][1] + lambda2*dPhidX[in][1];
				dphidz = lambda*dPhidX[i][2] + lambda2*dPhidX[in][2];

				// convection to implicit
				if( RUnormal<0. ){
					apn[nj] += RUnormal;
					app     -= RUnormal;
				}

				// diffusion to implicit
				ViscAreaLen =  Visc*Face[iface].rlencos;
				app        +=  ViscAreaLen;
				apn[nj]    += -ViscAreaLen;
				ani[nj]     =  Cell[in].globalIdx;
				nj ++ ;

				// convection to source term. (high order schemes)
				if( RUnormal>0. )
				{
					// low order interpolation, 1st upwind
					pfi= Phi[i];
					// high order interpolation, e.g., 2nd upwind
					vec_minus( dx, Face[iface].x, Cell[i].x, 3);
					pfh= Phi[i] + vec_dot( dPhidX[i], dx, 3 );
				}
				else
				{
					// low order interpolation, 1st upwind
					pfi= Phi[in];
					// high order interpolation, e.g., 2nd upwind
					vec_minus( dx, Face[iface].x, Cell[in].x, 3);
					pfh= Phi[in] + vec_dot( dPhidX[in], dx, 3 );
				}
				fcs = RUnormal*(pfh-pfi);
				
				// diffusion to source term. ( compressible or incompressible )
				fde = Visc*( dphidx*sav1 + dphidy*sav2 + dphidz*sav3 );
				vec_minus( dxc, Cell[in].x, Cell[i].x, 3 );
				fdi = ViscAreaLen*( dphidx*dxc[0]+dphidy*dxc[1]+dphidz*dxc[2] );
			}
			sphi += fde - fdi - fcs;
		}
		
		// right hand side, including
		//   pressure, gravity and part of diffusion terms (explicit - implicit), 
		//   relaxation terms
		
		// central cell coef is stored for later use
		
		app   /= URF[5];  // relaxation


		PetscInt row = Cell[i].globalIdx;
		assert(row>=0&&row<dataPartition->nGlobal);//check
		for(int ii=0; ii!=nj; ii++ )
			assert(ani[ii]>=0&&ani[ii]<dataPartition->nGlobal);//check

		MatSetValue(As,row,row,app,INSERT_VALUES);// diagonal
		MatSetValues(As,1,&row,nj,ani,apn,INSERT_VALUES);// off-diagonal

		// right hand side
		double bsv = sphi + (1.-URF[5])*app*Phi[i];
		VecSetValue(bs,row,bsv,INSERT_VALUES);

	}


}
