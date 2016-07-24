#include <BasicType.h>
#include <navier.h>

/***************************************************
	reorder cells of different body rid (FLUID & SOLID) within each partition...
	reorder is done before face is constructed
***************************************************/
//cell
class ReorderFluidSolidOperator{
public:
	map<int,BdRegion>* rm;
	bool operator()(CellData& cell){
		map<int,BdRegion>::const_iterator iter = rm->find(cell.rid);
		return (iter->second).type2==0;//return true for fluid
	}
};

/**********************************************************
 *	reorder operator for bnd and face
 * ********************************************************/
//bnd
class ReorderBndOperator{
public:
	FaceData* faces;
	CellData* cells;
	map<int,BdRegion>* rm;
	bool operator()(BoundaryData& bnd){
		int iface = bnd.face;
		int ic = faces[iface].cell1;
		map<int,BdRegion>::const_iterator iter = rm->find(cells[ic].rid);
		return (iter->second).type2 == 0;
	}
};

//face
class ReorderFaceOperator{
public:
	CellData* cells;
	map<int,BdRegion>* rm;
	bool operator()(FaceData& face){
		int ic = face.cell1;
		map<int,BdRegion>::const_iterator iter = rm->find(cells[ic].rid);
		return (iter->second).type2 == 0;
	}
};


int NavierStokesSolver::reorderCell(){
	if(!SolveConjungateHeat){
		Nfluid = Ncel;
		Nsolid = 0;
		return 0;
	}
	Nfluid = dataPartition->nLocal;
	Nsolid = dataPartition->nLocalSolid;

	map<int,int> _reorder;
	for(int i=0;i!=Ncel;++i){
		Cell[i].globalIdx = i;// globalIdx is used to mark original position...
	}

	//partition...
	ReorderFluidSolidOperator op;
	op.rm = &regionMap;
	for(int i=0;i!=Ncel;++i){
		assert(Cell[i].nface!=0);
	}
	std::partition(Cell,Cell+Ncel,op);
	for(int i=0;i!=Ncel;++i){
		assert(Cell[i].nface!=0);
	}

	for(int i=0;i!=Ncel;++i){
		BdRegion& reg = regionMap[Cell[i].rid];
		assert(reg.type1==5);
		int originalIdx = Cell[i].globalIdx;
		_reorder.insert(make_pair(originalIdx,i));
	}
	//assert(_nf==dataPartition->nLocal);
	//assert(_ns==dataPartition->nLocal+dataPartition->nLocalSolid);


	//update cell index in interfaceInfo
	for(map<int,Interface>::iterator iter = dataPartition->interfaces.begin();iter!=dataPartition-> interfaces.end();++iter){
		for(vector<int>::iterator iterIn = iter->second.sendposis.begin();iterIn!=iter->second.sendposis.end();++iterIn){
			map<int,int>::const_iterator _reorderIter = _reorder.find(*iterIn);
			if(_reorderIter!=_reorder.end()){
				(*iterIn) = (_reorderIter->second);
			}else{
				assert(false);
			}
		}
	}
	return 0;

}



/***************************************************
	find coupled boundary between FLUID & SOLID
***************************************************/
int NavierStokesSolver::buildCoupledFace(){
	if(!SolveConjungateHeat){
		NfluidFac = Nfac;	
		NfluidBnd = Nbnd;
		return 0;
	}
	for(int i=0;i!=Ncel;++i){
		assert(Cell[i].nface != 0);
	}

	set<int> boundFaceIDSet;
	NcoupledBnd = 0;
	for(int i=0;i!=Nfac;++i){
		int c1 = Face[i].cell1;
		int c2 = Face[i].cell2;
		int rid1 = regionMap[Cell[c1].rid].type2;
		int rid2 = regionMap[Cell[c2].rid].type2;
		if( (rid1==0 && rid2==1 )||
			(rid1==1 && rid2==0 ) ){
			boundFaceIDSet.insert(i);
			NcoupledBnd++;
			int* thisCell = std::find(Cell[c1].cell, Cell[c1].cell+Cell[c1].nface,c2);
			if(thisCell!= Cell[c1].cell+Cell[c1].nface){
				*thisCell = VOID_CELL_ON_BOUNDARY;
			}else{
				assert(false);
			}
			thisCell = std::find(Cell[c2].cell, Cell[c2].cell+Cell[c2].nface,c1);
			if(thisCell!= Cell[c2].cell+Cell[c2].nface){
				*thisCell = VOID_CELL_ON_BOUNDARY;
			}else{
				cout<<c2<<endl;
				assert(c2>=Ncel);//interface might not point to this face
			}
		}else if(rid1==1&&rid2==1){//face inside solid body
			Face[i].bnd = INNER_FACE_BOUNDARY_SOLID;
		}

		/*
		if(rid1==1&&rid2==0){//coupled face must point from FLUID --> SOLID
			Face[i].reverse();
		}
		*/
	}

	Nbnd += 2*NcoupledBnd;
	Nfac += NcoupledBnd;
	Bnd = (BoundaryData*) realloc(Bnd,Nbnd*sizeof(BoundaryData)); //realloc will not trigger Constructor
	Face= (FaceData*)     realloc(Face,Nfac*sizeof(FaceData));
	if(Face==NULL){
		errorHandler.fatalRuntimeError("run out of memory when creating coupeld boundary");
	}
	if(Bnd==NULL){
		errorHandler.fatalRuntimeError("run out of memory when creating coupeld boundary");
	}

	//build coupeld face...
	int icpFace = Nfac - NcoupledBnd;
	set<int>::const_iterator cpFaceIter;
	for(cpFaceIter = boundFaceIDSet.begin();cpFaceIter!=boundFaceIDSet.end();cpFaceIter++){
		int c1 = Face[(*cpFaceIter)].cell1;
		int c2 = Face[(*cpFaceIter)].cell2;
		//int rid1 = regionMap[Cell[c1].rid].type2;
		//int rid2 = regionMap[Cell[c2].rid].type2;
		
		//using default copy constructor
		Face[icpFace] = Face[(*cpFaceIter)];
		Face[icpFace].reverse();
		Face[icpFace].cell1 = c2;
		Face[(*cpFaceIter)].cell2 = VOID_CELL_ON_BOUNDARY-icpFace;
		Face[icpFace].cell2 = VOID_CELL_ON_BOUNDARY-(*cpFaceIter);
		Face[icpFace].lambda = 1.0;
		Face[(*cpFaceIter)].lambda = 1.0;
		
		double dx[3];
		vec_minus( dx,Face[(*cpFaceIter)].x,Cell[c1].x,3 );
		Face[(*cpFaceIter)].rlencos = Face[(*cpFaceIter)].area / (vec_dot( dx,Face[(*cpFaceIter)].n,3 )/Face[(*cpFaceIter)].area);
		vec_minus( dx,Face[(icpFace)].x,Cell[c2].x,3 );
		Face[(icpFace)].rlencos = Face[(icpFace)].area / (vec_dot( dx,Face[(icpFace)].n,3 )/Face[(icpFace)].area);

		//Xpac , Xnac on boundary face is meaning less ...
		Face[icpFace].Xpac[0] = Face[(*cpFaceIter)].Xpac[0] = 0.0;
		Face[icpFace].Xpac[1] = Face[(*cpFaceIter)].Xpac[1] = 0.0;
		Face[icpFace].Xpac[2] = Face[(*cpFaceIter)].Xpac[2] = 0.0;

		//congfig pointer from cell
		int* originalFaceID = std::find(Cell[c2].face,Cell[c2].face+Cell[c2].nface, (*cpFaceIter));
		if(originalFaceID!=Cell[c2].face+Cell[c2].nface){
			(*originalFaceID) = icpFace;
		}else{
			assert(c2>=Ncel);
		}
		icpFace++;
	}
	assert(icpFace==Nfac); 


	//insert Face[Nfac-NcoupledBnd : Nfac] ... the newly created coupled face
	for(int i=Nfac - NcoupledBnd;i!=Nfac;++i){
		boundFaceIDSet.insert(i);
	}
	assert(boundFaceIDSet.size()==(size_t) 2*NcoupledBnd);


	//Structure of Bnd: 0 ----> (Nbnd-NcoupledBnd) ---> Nbnd;
	int ridForCoupledBnd = regionMap.size();
	regionMap[ridForCoupledBnd].type1 = 1;// insert a new BdRegion;
	regionMap[ridForCoupledBnd].type2 = 2;
	cpFaceIter = boundFaceIDSet.begin();
	for(int i=Nbnd-2*NcoupledBnd;i!=Nbnd;++i){
		Bnd[i].rid = ridForCoupledBnd;
		int iface = *(cpFaceIter++);
		assert(iface<Nfac);
		Bnd[i].face = iface;
		Face[iface].bnd = i;
		assert(Face[iface].cell1>=0);
		for(int j=0;j!=4;++j){
			Bnd[i].vertices[j] = Face[iface].vertices[j];
		}
		for(int j=0;j!=3;++j)	
    		Bnd[i].shear[j] = 0.;
	}

	//reorder bnd
	for(long i=0;i!=Nbnd;++i){
		Bnd[i].info = (void*) i;
	}
	ReorderBndOperator op;
	op.rm = &regionMap;
	op.cells = Cell;
	op.faces = Face;
	NfluidBnd = std::count_if(Bnd,Bnd+Nbnd,op);
	std::partition(Bnd,Bnd+Nbnd,op);

	map<int,int> _reorder;
	for(int i=0;i!=Nbnd;++i){
		int originalIdx = (long) Bnd[i].info;
		_reorder.insert(make_pair(originalIdx,i));
	}
	//update index in face...
	for(int i=0;i!=Nfac;++i){
		if(Face[i].bnd<0) continue;
		map<int,int>::const_iterator iter = _reorder.find(Face[i].bnd);
		assert(iter!=_reorder.end());
		Face[i].bnd = iter->second;
	}

	//reorder face
	for(long i=0;i!=Nfac;++i){
		Face[i].info = (void*) i;
	}
	ReorderFaceOperator op2;
	op2.rm = &regionMap;
	op2.cells = Cell;
	NfluidFac = std::count_if(Face,Face+Nfac,op2);
	std::partition(Face,Face+Nfac,op2);
	_reorder.clear();
	for(int i=0;i!=Nfac;++i){
		int originalIdx = (long) Face[i].info;
		_reorder.insert(make_pair(originalIdx,i));
	}
	//update index in cell and bnd,and Face !...
	for(int i=0;i!=Ncel;++i){
		for(int j=0;j!=Cell[i].nface;++j){
			map<int,int>::const_iterator iter = _reorder.find( Cell[i].face[j] );
			assert(iter!=_reorder.end());
			Cell[i].face[j] = iter->second;
		}
	}
	for(int i=0;i!=Nbnd;++i){
		map<int,int>::const_iterator iter = _reorder.find(Bnd[i].face);
		assert(iter!=_reorder.end());
		Bnd[i].face = iter->second;
	}
	for(int i=0;i!=Nfac;++i){
		if(Face[i].cell2<VOID_CELL_ON_BOUNDARY){
			int originalIdx = COUPLED_FACE_ID(Face[i].cell2);
			map<int,int>::const_iterator iter = _reorder.find(originalIdx);
			assert(iter!=_reorder.end());
			Face[i].cell2 = VOID_CELL_ON_BOUNDARY - iter->second;
		}
	}

	//assertion check
	for(int i=0;i!=NfluidFac;++i){
		assert(regionMap[Cell[Face[i].cell1].rid].type2==0);
		assert(Face[i].cell1<Nfluid || Face[i].cell1>=Ncel);
	}
	for(int i=NfluidFac;i!=Nfac;++i){
		assert(regionMap[Cell[Face[i].cell1].rid].type2==1);
		assert(Face[i].cell1>=Nfluid);
	}
	for(int i=0;i!=NfluidBnd;++i){
		int iface = Bnd[i].face;
		assert(iface<NfluidFac);
		BdRegion& reg = regionMap[Bnd[i].rid];
		if(reg.type1==1&&reg.type2==2){
			int cpface = COUPLED_FACE_ID( Face[iface].cell2 );
			assert(cpface>=NfluidFac);
		}
	}
	for(int i=NfluidBnd;i!=Nbnd;++i){
		int iface = Bnd[i].face;
		assert(iface>=NfluidFac);
		assert(regionMap[Cell[Face[iface].cell1].rid].type2 == 1);
		BdRegion& reg = regionMap[Bnd[i].rid];
		if(reg.type1==1&&reg.type2==2){
			int cpface = COUPLED_FACE_ID( Face[iface].cell2 );
			assert(regionMap[Cell[Face[cpface].cell1].rid].type2 == 0);
			assert(cpface<NfluidFac);
		}
	}
	printf("Partition %d Coupled Boundary: %d\n",dataPartition->comRank,NcoupledBnd);
	
	return 0;
}


/***************************************************
		MAIN Solid Temperature Field Solver
***************************************************/
void NavierStokesSolver::UpdateSolidTemperature(){
	int i;
	double *ESource=NULL,*ApE=NULL,coef;


	ESource= new double[Nsolid];
	ApE    = new double[Nsolid];
	vec_init( ESource, Nsolid, 0. );
	vec_init( ApE,     Nsolid, 0. );
	ESource -= Nfluid;	//special treatment to make this array accessible by Nfluid-->Ncel index
	ApE     -= Nfluid;
	// prepare the diffusion coefficient and source terms

	// source terms, e.g., energy release, condensation/vaporization
	if( !IfSteady ){
		if(      TimeScheme==1 ){  // Euler forwards
			for( i=Nfluid; i!=Ncel; i++ ){
				coef       = Rn[i]/dt * Cell[i].vol;
				ApE[i]    += coef;
				ESource[i]+= coef * Tnp[i];
			}
		}
		else if( TimeScheme==2 ){  // 2nd order BDF
			for( i=Nfluid; i!=Ncel; i++ ){
				coef       = Rn[i]/dt * Cell[i].vol;
				ApE[i]    += 1.5*coef;
				ESource[i]+= coef * (2*Tnp[i]-0.5*Tnp2[i]);
			}
		}
	}
	// build matrix
	BuildSolidMatrix(Tn,BTem,VisLam,ESource,ApE,dataPartition->ASolid,dataPartition->bSolid );
	// Solve equations
	try{
		dataPartition->solveSolidTemp_GMRES(1.e-8,500,Tn+Nfluid);
	}catch(ConvergeError& err){
		char temp[256];	
		sprintf(temp,"solid temperature not converge in iter: %d, res %f\n",err.iter,err.residual);
		errorHandler.fatalRuntimeError(temp);
	}

	//********************MPI INTERFACE COMMUNICATION*****************************//
	dataPartition->interfaceCommunicationBegin(Tn);
	dataPartition->interfaceCommunicationEnd();

	
	// clipping work
	ESource+=Nfluid;
	ApE += Nfluid;
	delete [] ESource;
	delete [] ApE;
}


void NavierStokesSolver::BuildSolidMatrix(double *Phi, double *BPhi, double *DiffCoef, double *source, double *App,Mat& A,Vec& b){
	int i,j,iface, ip,in,ani[6],nj,bnd;
	double app,apn[6],lambda,lambda2, Visc,dxc[3],
		dphidx,dphidy,dphidz, 
		sav1,sav2,sav3,ViscAreaLen, 
		fde,fdi,   sphi;

	solidGradient ( Phi, BPhi,  dPhidX );

	dataPartition->interfaceCommunicationBegin(DiffCoef);
	dataPartition->interfaceCommunicationEnd();

	for( i=Nfluid; i<Ncel; i++ )
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
				BdRegion& reg = regionMap[Bnd[bnd].rid];
				//The face point from outside to inside on coupled bound!
				sav1    = Face[iface].n[0];
				sav2    = Face[iface].n[1];
				sav3    = Face[iface].n[2];

				// diffusion boundary
				Visc   = DiffCoef[i]; // This Should be changed using boundary condition. e.g., BVisTur[bnd]
				dphidx = dPhidX[i][0];
				dphidy = dPhidX[i][1];
				dphidz = dPhidX[i][2];
				vec_minus( dxc, Face[iface].x, Cell[i].x, 3 );
				if(reg.type1 == 1 && reg.type2 == 0){ //given T
  					ViscAreaLen = Visc*Face[iface].rlencos;
					app += ViscAreaLen;
					fde  = Visc*( dphidx*sav1 + dphidy*sav2 + dphidz*sav3 );
					fdi  = ViscAreaLen*( dphidx*dxc[0] + dphidy*dxc[1] + dphidz*dxc[2] - BPhi[bnd] );
				}else if(reg.type1 ==1 && reg.type2 == 2){//coupled boundary
					//Visc = Bnd[bnd].h / (Face[iface].rlencos / Face[iface].area);
					//ViscAreaLen = Face[iface].area * Bnd[bnd].h;
					fde = 0. - Bnd[bnd].q * Face[iface].area;// coupled boundary  
					fdi = 0.;	
				}else if(reg.type1 ==4 || reg.type1 ==1){ //symmetric and given flux
					fde = 0. - Bnd[bnd].q * Face[iface].area;// flux boundary to source
					fdi = 0.;	
				}else{
					errorHandler.fatalRuntimeError("such boundary should not in SOLID BODY!");
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
					lambda  = 1.- Face[iface].lambda;
				}
				else{
					sav1    = Face[iface].n[0];
					sav2    = Face[iface].n[1];
					sav3    = Face[iface].n[2];
					lambda  = Face[iface].lambda;
				}
				
				lambda2= 1.- lambda;
				Visc   = lambda*DiffCoef[i]  + lambda2*DiffCoef[in];
				dphidx = lambda*dPhidX[i][0] + lambda2*dPhidX[in][0];
				dphidy = lambda*dPhidX[i][1] + lambda2*dPhidX[in][1];
				dphidz = lambda*dPhidX[i][2] + lambda2*dPhidX[in][2];


				// diffusion to implicit
				ViscAreaLen =  Visc*Face[iface].rlencos;
				app        +=  ViscAreaLen;
				apn[nj]    += -ViscAreaLen;
				ani[nj]     =  Cell[in].globalIdx;
				nj ++ ;

				// diffusion to source term. ( compressible or incompressible )
				fde = Visc*( dphidx*sav1 + dphidy*sav2 + dphidz*sav3 );
				vec_minus( dxc, Cell[in].x, Cell[i].x, 3 );
				fdi = ViscAreaLen*( dphidx*dxc[0]+dphidy*dxc[1]+dphidz*dxc[2] );
			}
			sphi += fde - fdi;
		}
		
		// right hand side, including
		//   pressure, gravity and part of diffusion terms (explicit - implicit), 
		//   relaxation terms
		
		// central cell coef is stored for later use
		
		app   /= URF[7];  // relaxation

		PetscInt row = Cell[i].globalIdx;
		assert(row>=0&&row<dataPartition->nGlobalSolid);//check
		for(int ii=0; ii!=nj; ii++ )
			assert(ani[ii]>=0&&ani[ii]<dataPartition->nGlobalSolid);//check

		MatSetValue(A,row,row,app,INSERT_VALUES);// diagonal
		MatSetValues(A,1,&row,nj,ani,apn,INSERT_VALUES);// off-diagonal

		// right hand side
		double bsv = sphi + (1.-URF[7])*app*Phi[i];
		VecSetValue(b,row,bsv,INSERT_VALUES);

	}
}


void NavierStokesSolver::SetBCSolidTemperature(double** gradientT,double* fluidDiffCoef ){
	int    i,rid,iface;
	double fluxrecord = 0.0;
	for( i=NfluidBnd; i<Nbnd; i++ )
	{
		rid   = Bnd[i].rid;
		iface = Bnd[i].face;
		//int ic    = Face[iface].cell1; // THE ONLY PART DIFFER FROM FLUID
		switch( regionMap[rid].type1 ){
		case(1):  // wall
			//remain initial
			if(regionMap[rid].type2==2){//coupled boundary
				int cpface = COUPLED_FACE_ID( Face[iface].cell2 );
				assert(cpface<NfluidFac && cpface >= 0);
				int ifluidcell = Face[cpface].cell1;
				assert(ifluidcell<Nfluid||ifluidcell>=Ncel);
				double dtdn = (
					gradientT[ifluidcell][0] * Face[cpface].n[0] +
					gradientT[ifluidcell][1] * Face[cpface].n[1] +
					gradientT[ifluidcell][2] * Face[cpface].n[2]
				) / Face[cpface].area;
				dataPartition->PRINT_LOG(dtdn);
				//flux is deducted from fluid side;
				int cpbnd = Face[cpface].bnd;
				assert(cpbnd>=0 && cpbnd<NfluidBnd);
				Bnd[cpbnd].q = (0. - dtdn * fluidDiffCoef[ifluidcell])*Face[cpface].area;//flux from fluid to solid;
				Bnd[i].q = -Bnd[cpbnd].q;
				fluxrecord += Bnd[i].q;
				//temperature is deducted from solid side;		
			}
			break;
		case(2):  // inlet
		case(3):
		case(4):
		case(6)://periodic
			//pass
			break;
		default:
			char temp[256];
			sprintf(temp,"no such boundary type %d \n ",rid);
			errorHandler.fatalLogicError(temp);
		}
	}
	double fluxSum;
	MPI_Reduce(&fluxrecord,&fluxSum,1,MPI_DOUBLE,MPI_SUM,root.rank,dataPartition->comm);
	PetscPrintf(dataPartition->comm,"flux to solid %e\n",fluxSum);
}


#define PETSC_SOLVE_VERBOSE
int DataPartition::solveSolidTemp_GMRES(double tol,int maxIter,double const* xs){
	KSPConvergedReason reason;
	int iters;
	double residule;

	MatAssemblyBegin(ASolid,MAT_FINAL_ASSEMBLY);
	VecAssemblyBegin(bSolid);
	MatAssemblyEnd(ASolid,MAT_FINAL_ASSEMBLY);
	VecAssemblyEnd(bSolid);	

	MPI_Barrier(comm);
#ifdef PETSC_SOLVE_VERBOSE
	PetscPrintf(comm,"begin solid temperature solve\n");
#endif


	KSPSetOperators(kspSolid,ASolid,ASolid);

	ierr = KSPSetType(kspSolid,KSPGMRES);CHKERRQ(ierr);

	KSPSetInitialGuessNonzero(kspSolid,PETSC_TRUE);

	/***************************************
	 *      SET  TOLERENCE
	 ***************************************/
	KSPSetTolerances(kspSolid,tol,PETSC_DEFAULT,PETSC_DEFAULT,maxIter);	//absolute residule

	/***************************************
	 * 	ILU preconditioner:
	 ***************************************/
	//KSPGetPC(kspSolid,&pcSolid);
	KSPSetFromOptions(kspSolid);//can override settings from command line
	KSPSetUp(kspSolid); //the precondition is done at this step


	/***************************************
	 * 	SOLVE Scarlar
	 ***************************************/
	ierr = VecCreateMPIWithArray(comm,1,nLocalSolid,nGlobalSolid,xs,&xsolSolid);CHKERRQ(ierr); 
	ierr = VecAssemblyBegin(xsolSolid);			CHKERRQ(ierr);
	ierr = VecAssemblyEnd(xsolSolid);			CHKERRQ(ierr);
	int vn,bn,mn,mm;
	VecGetSize(xsolSolid,&vn);
	VecGetSize(bSolid,&bn);
	MatGetSize(ASolid,&mn,&mm);
	assert(vn==nGlobalSolid);
	assert(bn==nGlobalSolid);
	assert(mn==nGlobalSolid &&mm==nGlobalSolid);
	
	ierr = KSPSolve(kspSolid,bSolid,xsolSolid);CHKERRQ(ierr);

	KSPGetConvergedReason(kspSolid,&reason);

	if(reason<0){
		KSPGetIterationNumber(kspSolid,&iters);
		KSPGetResidualNorm(kspSolid,&residule);
		throw ConvergeError(iters,residule,"Tn");
	}else if(reason ==0){
		PetscPrintf(comm,"why is this program still running?\n");
	}else{
#ifdef PETSC_SOLVE_VERBOSE
		KSPGetIterationNumber(kspSolid,&iters);
		PetscPrintf(comm,"KSP GMRES - SOLID Temperature converged in %d step! :)\n",iters);
#endif
	}

#ifdef PETSC_SOLVE_VERBOSE
	double sbnorm,snorm;
	VecNorm(xsolSolid,NORM_2,&snorm);
	VecNorm(bSolid,NORM_2,&sbnorm);
	PetscPrintf(comm,"solid temp norm %e, solid temp bnorm %e\n",snorm,sbnorm);
#endif
	ierr = VecDestroy(&xsolSolid);CHKERRQ(ierr);
	ierr = MatZeroEntries(ASolid);CHKERRQ(ierr);
	
	
	return 0;
}


int NavierStokesSolver::solidGradient(double *phi, double *Bphif, double **phigd){
	int    i,g, c1,c2;
	double lambda,pf;
	// using Gauss theorem
	for( i=Nfluid; i<Ncel; i++ )
    	for(g=0;g<3;g++)
        		phigd[i][g]= 0.;


	for( i=NfluidFac;i<Nfac; i++ )
	{
       	lambda = Face[i].lambda;
       	c1     = Face[i].cell1;
       	c2     = Face[i].cell2;
        assert(Face[i].bnd!=INNER_FACE_BOUNDARY) ;//ONLY place differ from fluid part
            //face in fluid part
        if( Face[i].bnd<0){
		pf = lambda*phi[c1] + (1.-lambda)*phi[c2];
		for( g=0;g<3;g++ ){
			phigd[c1][g] += pf * Face[i].n[g];
			phigd[c2][g] -= pf * Face[i].n[g];//might add to interface cell, meaningless but no harmful
		}
        }else{
	 	pf = Bphif[Face[i].bnd]; // how to add boundary condition ?
		for( g=0;g<3;g++ )
			phigd[c1][g] += pf * Face[i].n[g];
            		//assure face is on solid-fluid boundary
	}
	}

	//CHECK_ARRAY(phi,Ncel);
	//CHECK_ARRAY(ED,Ncel);

   	for( i=Nfluid; i<Ncel; i++ ){
       	for( g=0; g<3; g++ ){
       	    phigd[i][g] /= Cell[i].vol;
		}
    }

    switch(limiter){
    case 0:	
  		//pass 
	break;
    case 1:
    case 2:
    case 3:
	Limiter_Barth( phi, phigd,Nfluid,Ncel);
	break;
    default:
    	errorHandler.fatalLogicError("no such limiter",limiter);
    }

	//CHECK_ARRAY(phigd[0],3*Ncel);
	dataPartition->interfaceCommunicationBegin(phigd);			
	dataPartition->interfaceCommunicationEnd();
	return 0;
}

