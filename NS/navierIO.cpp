
#include"navier.h"
void OutArray2File(double arr[],int N, stringstream& of){
	for(int i=0; i<N; i++ ){
		of<<arr[i]<<"  ";
		if( i%5==0 ) of<<endl;
	}
	of<<endl;
}

void outputVTKScalar( const char name[], double arr[],int N, ofstream &of)
{
	int i;
	of<<"SCALARS "<<name<<" FLOAT"<<endl;
	of<<"LOOKUP_TABLE default"<<endl;
	for( i=0; i<N; i++ )
		of<<float(arr[i])<<endl;
}


// output grid to tecplot for validation
void NavierStokesSolver::OutputGrid()
{
	int i;
	string title("grid.dat");
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
	sprintf(cbuffer,"zone n=%d, e=%d, f=fepoint, et=BRICK\n",Nvrt,Ncel);
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

	parallelWriteBuffer(title,sbuffer,dataPartition,head);


	of.close( );
	//PetscPrintf(dataPartition->comm,"complete writing grid.dat\n");
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
