
#include"navier.h"
void OutArray2File(double arr[],int N, string& str,int step){
	char* tmp = new char[20*N + (N/5)+1] ;
	int posi = 0;
	for(int i=0; i< N*step ; i+=step ){
		sprintf(tmp+posi,"%19.5f ",arr[i]);
		posi+=20;
		if( (i/step)%5==4 ){
			sprintf(tmp+posi,"\n");
			posi++;
		}
	}
	str.append(tmp);
	str.append("\n");
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


void NavierStokesSolver::OutputGridBinary(){
	size_t numberofInt = 8*Ncel+5; //partition info 3
	size_t numberofDouble = 3*Nvrt;
	size_t bufferSize = sizeof(int)*numberofInt + sizeof(double)*numberofDouble;
	char* buffer = new char[bufferSize];
	int* _ptr = (int*)buffer;
	_ptr[0] = dataPartition->comRank;
	_ptr[1] = dataPartition->nProcess;
	_ptr[2] = Ncel;
	_ptr[3] = Nfluid;
	_ptr[4] = Nvrt;

	double* _dptr = (double*)(buffer+4*sizeof(int));
	for(int i=0;i!=Nvrt;++i)
		for(int j=0;j!=3;++j)
			_dptr[i*3+j] = Vert[i][j];

	_ptr = (int*)(buffer + 4*sizeof(int)+ + 3*Nvrt*sizeof(double));
	for(int i=0;i!=Ncel;++i)
		for(int j=0;j!=8;++j)
			_ptr[i*8+j] = Cell[i].vertices[j] + 1;

	ofstream of("tec/geo",ios::binary);
	parallelWriteBuffer("tec/geo",buffer,bufferSize,dataPartition,0);
	delete [] buffer;
	of.close( );
}


// output grid to tecplot for validation
void NavierStokesSolver::OutputGrid(const string& title,int type,int N)
{
	int i;
	ofstream of;
	int head =0;
	
	of.open(title.c_str());
	string sbuffer;
	char cbuffer[256];

	sprintf(cbuffer,"title = \"cycas partition grid\"\nvariables=\"x\",\"y\",\"z\",\"p\"\n");
	sbuffer = cbuffer;
	head = sbuffer.size();
	if(dataPartition->comRank == root.rank){//print head in root
		of<<sbuffer;
	}

	sbuffer="";
	//sprintf(cbuffer,"ZONE T=part%d N=%d, E=%d, VARLOCATION=([1-3]=NODAL,[4-%d]=CELLCENTERED) DATAPACKING=BLOCK, ZONETYPE=FEBRICK\n",dataPartition->comRank,Nvrt,Ncel,4);
	sprintf(cbuffer,"zone n=%d, e=%d, f=fepoint, et=BRICK\n",Nvrt,N);
	sbuffer=cbuffer;

	map<int, set<int> > nodesSets;
	for(map<int,Interface>::iterator it = dataPartition->interfaces.begin(); it!=dataPartition->interfaces.end(); ++it){
		set<int> nodesSet;
		for(vector<set<int> >::iterator iter = it->second.boundNodes.begin(); iter!=it->second.boundNodes.end(); ++iter){
			for(set<int>::iterator iter2 = iter->begin();iter2!=iter->end();++iter2){
				nodesSet.insert(*iter2);
			}
		}
		nodesSets[it->first] = nodesSet;
	}
	
	for(i =0;i!=Nvrt;++i){
		double _iinfo = 0.0;
		for(map<int,set<int> >::iterator iter = nodesSets.begin();iter!=nodesSets.end();++iter){
			if( iter->second.find(i)!=iter->second.end()){
				_iinfo = (double) iter->first + 1; 
				break;
			}
		}
		sprintf(cbuffer,"%f %f %f %f\n",Vert[i][0],Vert[i][1],Vert[i][2],_iinfo);
		sbuffer.append(cbuffer);

	}

	for( i=0; i<Ncel; i++ ){
		if(regionMap[Cell[i].rid].type2 != type) continue;
		sprintf(cbuffer,"%8d %8d %8d %8d %8d %8d %8d %8d\n",
				Cell[i].vertices[0]+1,
				Cell[i].vertices[1]+1,
				Cell[i].vertices[2]+1,
				Cell[i].vertices[3]+1,
				Cell[i].vertices[4]+1,
				Cell[i].vertices[5]+1,
				Cell[i].vertices[6]+1,
				Cell[i].vertices[7]+1
				);
		sbuffer.append(cbuffer);
	}

	//printf("rank:%d printing %lu to buffer\n",dataPartition->comRank,sbuffer.size());

	parallelWriteBuffer(title,sbuffer.c_str(),sbuffer.size(),dataPartition,head);
	of.close( );
	//PetscPrintf(dataPartition->comm,"complete writing grid.dat\n");
}

void OutArray2Buffer(const double* val, double* buffer,int& iptr,int first,int N){
	for(int i=first;i!=N;++i)
		buffer[iptr+i] = val[i];
	iptr+=N;
}

void NavierStokesSolver::Output2TecplotBinary(){
	PetscPrintf(dataPartition->comm,"*****************************    OUTPUT TECPLOT   *********************************\n");
	size_t solidVar = 0;
	size_t bufferSize = 3 * sizeof(int) + field->nVar*Nfluid*sizeof(double);
	size_t solidBufferSize = 0;
	if(SolveConjungateHeat){
		solidVar =  1;
		solidBufferSize = (Ncel-Nfluid)*sizeof(double);
	}
	char* buffer = new char[bufferSize + solidBufferSize];
	int* ptr = (int*)buffer;
	ptr[0] = dataPartition->comRank;
	ptr[1] = field->nVar;
	ptr[2] = solidVar;
	double* dptr = (double*)(buffer+3*sizeof(int));
	int iptr = 0;
	OutArray2Buffer(Pn,dptr,iptr,0,Nfluid);	
	OutArray2Buffer(Un,dptr,iptr,0,Nfluid);	
	OutArray2Buffer(Vn,dptr,iptr,0,Nfluid);	
	OutArray2Buffer(Wn,dptr,iptr,0,Nfluid);	
	OutArray2Buffer(Rn,dptr,iptr,0,Nfluid);	
	OutArray2Buffer(Tn,dptr,iptr,0,Nfluid);	
	if(TurModel==1){
		OutArray2Buffer(TE,dptr,iptr,0,Nfluid);
	       	OutArray2Buffer(ED,dptr,iptr,0,Nfluid);
	}
	assert(iptr==Nfluid*field->nVar);

	if(SolveConjungateHeat){//only 1 solid var
	       	OutArray2Buffer(Tn,dptr,iptr,Nfluid,Ncel);
	}
	ofstream of;
	of.open("tec/data",ios::binary);
	/*****************MPI PARALLEL I/O APIs*************************/
	parallelWriteBuffer("tec/data",buffer,bufferSize,dataPartition,0);//collective
	of.close();

	delete [] buffer;
}

void NavierStokesSolver::Output2Tecplot()

{
	PetscPrintf(dataPartition->comm,"*****************************    OUTPUT TECPLOT   *********************************\n");

	ofstream of;
	char tecTitle[256];
	sprintf(tecTitle,"tec/res%04d.dat",int(this->outputCounter));
	of.open(tecTitle);
	
	string head;
	char tmpchar[256];
	sprintf(tmpchar,"variables= \"x\", \"y\", \"z\", \"p\", \"u\", \"v\", \"w\", \"ro\", \"T\",\"Ma\",\"mu\" ");
	head = tmpchar;
	if(TurModel==1){
		sprintf(tmpchar,"\"te\", \"ed\"");
		head.append(tmpchar);
	}
	head.append("\n");

	if(dataPartition->comRank==root.rank){
		of<<head;
	}

	MPI_Barrier(dataPartition->comm);

	/*****PRINTING OTHER PARTS...**********/

	string tecBuffer;

	sprintf(tmpchar,"ZONE N=%d, E=%d, VARLOCATION=([1-3]=NODAL,[4-%d]=CELLCENTERED) DATAPACKING=BLOCK, ZONETYPE=FEBRICK\n",Nvrt,Ncel,3+field->nVar+2);
	tecBuffer = tmpchar;



	OutArray2File(Vert[0]+0,Nvrt,tecBuffer,3);//x
	OutArray2File(Vert[0]+1,Nvrt,tecBuffer,3);//y
	OutArray2File(Vert[0]+2,Nvrt,tecBuffer,3);//z

	OutArray2File( Pn,Ncel, tecBuffer,1);
	OutArray2File( Un,Ncel, tecBuffer,1);
	OutArray2File( Vn,Ncel, tecBuffer,1);
	OutArray2File( Wn,Ncel, tecBuffer,1);
	OutArray2File( Rn,Ncel, tecBuffer,1);
	OutArray2File( Tn,Ncel, tecBuffer,1);
	
	double* tmp = new double[Ncel];
	for( int i=0; i<Ncel; i++ ){  // Mach / velocity magnitude
		if( DensityModel==1 )
			tmp[i]= sqrt( (Un[i]*Un[i]+Vn[i]*Vn[i]+Wn[i]*Wn[i])/(gama*(Pn[i]+PressureReference)/Rn[i]) );
		else
			tmp[i]= sqrt( (Un[i]*Un[i]+Vn[i]*Vn[i]+Wn[i]*Wn[i]) ); 
	}

	OutArray2File( tmp,Ncel,tecBuffer,1);
	for(int i=0; i<Ncel; i++ ){  // mu
		tmp[i]= VisLam[i] + VisTur[i];
	}
	OutArray2File( tmp,Ncel,tecBuffer,1);
	delete []tmp;
	tmp=NULL;

	if( TurModel==1 ){
		OutArray2File( TE,Ncel,tecBuffer,1 );
		OutArray2File( ED,Ncel,tecBuffer,1 );
	}

	for(int i=0; i<Ncel; i++ ){ //connectivity
			sprintf(tmpchar,"%d %d %d %d %d %d %d %d\n ",
					Cell[i].vertices[0]+1,
					Cell[i].vertices[1]+1,
					Cell[i].vertices[2]+1,
					Cell[i].vertices[3]+1,
					Cell[i].vertices[4]+1,
					Cell[i].vertices[5]+1,
					Cell[i].vertices[6]+1,
					Cell[i].vertices[7]+1
					);
			tecBuffer.append(tmpchar);
	}

	/*****************MPI PARALLEL I/O APIs*************************/
	parallelWriteBuffer(tecTitle,tecBuffer.c_str(),tecBuffer.size(),dataPartition, head.size());//collective
	of.close();
	
	
}



/***************************************************
 *	write binary geometry backup to local storage
 *	this routine is local
 ***************************************************/
void NavierStokesSolver::writeGeometryBackup(int vsize,double* vbuffer,int esize, int* ebuffer, int isize, int* ibuffer){
	char title[256];	
	sprintf(title,"localGeometryBackup/grid_backup_rank%04d.dat",dataPartition->comRank);
	ofstream of(title);

	//partition info	
	of.write((char*)&(dataPartition->nProcess),sizeof(dataPartition->nProcess));
	of.write((char*)&(dataPartition->nLocal),sizeof(dataPartition->nLocal));	
	for(int i=0;i!=dataPartition->nProcess;++i){
		of.write((char*)&(dataPartition->gridList[i]),sizeof(dataPartition->gridList[i]));	
	}

	//geometry buffer
	for(int i=0;i!=vsize;++i)
		of.write((char*)&(vbuffer[i]),sizeof(vbuffer[i]));

	for(int i=0;i!=esize;++i)
		of.write((char*)&(ebuffer[i]),sizeof(ebuffer[i]));


	for(int i=0;i!=isize;++i)
		of.write((char*)&(ibuffer[i]),sizeof(ibuffer[i]));


	of.close();
	//PetscPrintf(dataPartition->comm,"backup file written: %s\n",title);
	//printf("backup file written: %s\n",title);
}


/*******************************************************
 *	original cycas member
 *******************************************************/
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


void NavierStokesSolver::writeTotFile(){
	return; //to be implemented
}


void NavierStokesSolver::WriteBackupFile( )
{
	int i;
	ofstream of;
	char title[256];
	sprintf(title,"localGeometryBackup_%d/res.sav",dataPartition->comRank);
	of.open(title);
	
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
	map<int,vector<double> > regionRecord;
	//wall monitor
	for(map<int,BdRegion>::iterator iter = regionMap.begin();iter!=regionMap.end();++iter){
		for(int i=0;i!=3;++i)
			iter->second.vectorRecord[i] = 0.0;
		iter->second.scarlarRecord = 0.0;	
	}
	for(int i=0;i!=Nbnd;++i){
		map<int,BdRegion>::iterator iter = regionMap.find( Bnd[i].rid );
		BdRegion& reg = iter->second;
		if(reg.type1 == 1){ 	//record wall  drag / lift force
			reg.vectorRecord[0]  += Bnd[i].shear[0];
			reg.vectorRecord[1]  += Bnd[i].shear[1];
			reg.vectorRecord[2]  += Bnd[i].shear[2];
			reg.scarlarRecord += Bnd[i].q;
		}
	}

	for(int i=0;i!=Ncel;++i){
		map<int,BdRegion>::iterator iter = regionMap.find( Cell[i].rid );
		BdRegion& reg = iter->second;
		if(reg.type1==5 && reg.type2==0){
			reg.scarlarRecord += pow(Wn[i] ,2) * Cell[i].vol;
		}
	}

	for(map<int,BdRegion>::iterator iter = regionMap.begin();iter!=regionMap.end();++iter){
		int rid  = iter->second.type1;
		if(rid != 1 && rid !=5 ) continue;
		if(iter->first == 0 ) continue;

		double sumScar = 0.0;
		double sumVec[3] = {0.,0.,0.};

		MPI_Reduce(&( iter->second.scarlarRecord ),&sumScar,1,MPI_DOUBLE,MPI_SUM,root.rank,MPI_COMM_WORLD);	
		MPI_Reduce( iter->second.vectorRecord , sumVec,3,MPI_DOUBLE,MPI_SUM,root.rank,MPI_COMM_WORLD);	
		vector<double> toInsert(sumVec,sumVec+3);
		toInsert.push_back(sumScar);
		regionRecord.insert( make_pair( iter->first, toInsert  ) );
		MPI_Barrier(dataPartition->comm);
	}

	//Init output files
	if(dataPartition->comRank == root.rank){
		char temp[4096];	
		sprintf(temp,"%20e || ",cur_time);
		root.writeMonitorFile(dataPartition,temp);
		for(map<int,vector<double> >::const_iterator iter = regionRecord.begin();iter!=regionRecord.end();++iter){

			sprintf(temp,"%20e, %20e, %20e | %20e || ",
					iter->second[0],
					iter->second[1],
					iter->second[2],
					iter->second[3]
				);
			root.writeMonitorFile(dataPartition,temp);
		}
		root.writeMonitorFile(dataPartition,"\n");
	}
	MPI_Barrier(dataPartition->comm);

}
