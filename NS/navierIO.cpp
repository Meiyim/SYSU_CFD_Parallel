
#include"navier.h"
void OutArray2File(double arr[],int N, stringstream& of){
	for(int i=0; i<N; i++ ){
		of<<arr[i]<<"  ";
		if( i%5==0 ) of<<endl;
	}
}

void outputVTKScalar( const char name[], double arr[],int N, ofstream &of)
{
	int i;
	of<<"SCALARS "<<name<<" FLOAT"<<endl;
	of<<"LOOKUP_TABLE default"<<endl;
	for( i=0; i<N; i++ )
		of<<float(arr[i])<<endl;
}

 /*******************************************************
  *		parallel I/O
  *		collective on MPI
  *******************************************************/
void NavierStokesSolver::writeTecZoneParallel(const string& title){
	size_t charLen = TECPLOT_NVAR * Ncel + 256; //256: zone info length
	char* tecBuffer = NULL;
	stringstream _strm;

	_strm<<"zone n="<<Nvrt<<", e="<<Ncel<<", VARLOCATION=([1-3]=NODAL,[4-"<<TECPLOT_NVAR<<"]=CELLCENTERED)"
		<<"DATAPACKING=BLOCK, ZONETYPE=FEBRICK"
		<<endl;
	for( int j=0; j<3; j++ ){
		for( int i=0; i<Nvrt; i++ ){
			_strm<<Vert[i][j]<<" ";
			if( i%5==0 ) _strm<<endl;
		}
	}	

	OutArray2File( Pn,Ncel, _strm);
	OutArray2File( Un,Ncel, _strm);
	OutArray2File( Vn,Ncel, _strm);
	OutArray2File( Wn,Ncel, _strm);
	OutArray2File( Rn,Ncel, _strm);
	OutArray2File( Tn,Ncel, _strm);
	
	double* tmp = new double[Ncel];
	for( int i=0; i<Ncel; i++ ){  // Mach / velocity magnitude
		if( DensityModel==1 )
			tmp[i]= sqrt( (Un[i]*Un[i]+Vn[i]*Vn[i]+Wn[i]*Wn[i])/(gama*(Pn[i]+PressureReference)/Rn[i]) );
		else
			tmp[i]= sqrt( (Un[i]*Un[i]+Vn[i]*Vn[i]+Wn[i]*Wn[i]) ); 
	}

	OutArray2File( tmp,Ncel,_strm );
	for(int i=0; i<Ncel; i++ ){  // mu
		tmp[i]= VisLam[i] + VisTur[i];
	}
	OutArray2File( tmp,Ncel,_strm );
	if( TurModel==1 ){
		OutArray2File( TE,Ncel,_strm );
		OutArray2File( ED,Ncel,_strm );
	}

	for(int i=0; i<Ncel; i++ ){
		for( int j=0;j<8;j++ )
			_strm<<Cell[i].vertices[j]+1<<" ";
		_strm<<endl;
	}
	delete []tmp;

	tecBuffer = new char[charLen];
	_strm >> tecBuffer;	
	/*****************MPI PARALLEL I/O APIs*************************/

	delete tecBuffer;
}


void NavierStokesSolver::Output2Tecplot()
{
	ofstream of;
	char tecTitle[256];
	sprintf(tecTitle,"res%04d.dat",int(this->outputCounter++));
	of.open(tecTitle);

	root.printTecplotHead(dataPartition,of);//root only
	
	writeTecZoneParallel(tecTitle);//collective call

	/*****PRINTING OTHER PARTS...**********/


}


/***************************************************
 *	write binary geometry backup to local storage
 *	this routine is local
 ***************************************************/
void NavierStokesSolver::writeGeometryBackup(int* ebuffer, double* vbuffer, int* ibuffer){
	char title[256];	
	sprintf(title,"localGeometryBackup/grid_backup_rank%04d.dat",dataPartition->comRank);
	ofstream of(title);
	of.write((char*)&Ncel,sizeof(Ncel));
	of.write((char*)&Nbnd,sizeof(Nbnd));
	of.write((char*)&(dataPartition->nProcess),sizeof(dataPartition->nProcess));
	for(int i=0;i!=dataPartition->nProcess;++i){
		of.write((char*)&(dataPartition->gridList[i]),sizeof(dataPartition->gridList[i]));
	}
	for(int i=0;i!=Ncel+Nbnd;++i)
		of.write((char*)&(ebuffer[i]),sizeof(ebuffer[i]));

	of.write((char*)&Nvrt,sizeof(Nvrt));
	for(int i=0;i!=Nvrt;++i)
		of.write((char*)&(vbuffer[i]),sizeof(vbuffer[i]));

	int ninterface = ibuffer[0];
	of.write((char*)&(ninterface),sizeof(ninterface));

	size_t counter = 1;
	for(int i=0;i!=ninterface;++i){
		int pid = ibuffer[counter++];
		int width = ibuffer[counter++];
		of.write((char*)&(pid),sizeof(pid));
		of.write((char*)&(width),sizeof(width));

		for(int j=0;j!=width;++j)
			of.write((char*)&(ibuffer[counter++]),sizeof(int));
	}

	of.close();
	printf("backup file written: %s\n",title);
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
