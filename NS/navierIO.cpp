
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




void NavierStokesSolver::Output2Tecplot()
{
	PetscPrintf(dataPartition->comm,"*****************************    OUTPUT TECPLOT   *********************************\n");

	ofstream of;
	char tecTitle[256];
	sprintf(tecTitle,"tec/res%04d.dat",int(this->outputCounter++));
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
	parallelWriteBuffer(tecTitle,tecBuffer,dataPartition, head.size());//collective
	
	
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
	map<int,vector<double> > wallForces;
	//wall monitor
	for(map<int,BdRegion>::const_iterator iter = regionMap.begin();iter!=regionMap.end();++iter){
		MPI_Barrier(dataPartition->comm);
		int bid  = iter->first;
		if (iter->second.type1 == 1) {//for walls
			double f[3] = {0.,0.,0.};
			double F[3] = {0.,0.,0.};
			for(int i=0;i!=Nbnd;++i){
				if(Bnd[i].rid == bid){
					f[0] += Bnd[i].shear[0];		
					f[1] += Bnd[i].shear[1];		
					f[2] += Bnd[i].shear[2];		
				}
			}
			MPI_Reduce(f,F,3,MPI_DOUBLE,MPI_SUM,root.rank,MPI_COMM_WORLD);	
			if(dataPartition->comRank == root.rank){
				wallForces.insert(make_pair(bid,vector<double> (F,F+3)));
			}
		}
	}
	if(dataPartition->comRank == root.rank){
		char temp[4096];	
		for(map<int,vector<double> >::const_iterator iter = wallForces.begin();iter!=wallForces.end();++iter){
			sprintf(temp,"wall: %d force: %e, %e, %e\n",
					iter->first,
					iter->second[0],
					iter->second[1],
					iter->second[2]
				);
			root.monitorFile<<temp;
		}
	}
	MPI_Barrier(dataPartition->comm);

}
