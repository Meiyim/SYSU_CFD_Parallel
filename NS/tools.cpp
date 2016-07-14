#include <iostream>
#include <math.h>
#include <string>
#include <ctime>
#include <mpi.h>
#include "tools.h"

using namespace std;
// vector manipulation. should be defined as inline function
void   vec_init  ( double a[], int n, double val )  // a[0:n-1] = val
{
    for( int i=0; i<n; i++ )
        a[i] = val;
}

void   vec_minus (double *x1, double *x2, double *x3, int n) // x1= x2 - x3
{
    for( int i=0; i<n; i++ )
        x1[i] = x2[i] - x3[i];
}

double vec_dot   (double *a, double *b, int n)  // Return = a . b
{
    double s=0.;
    for( int i=0; i<n; i++ )
        s += a[i]*b[i];
    return s;
}

double vec_len   (double *a, int n)
{
    double s=0.;
    for( int i=0; i<n; i++ )
        s += a[i]*a[i];
    return sqrt(s);
}

void   vec_cross (double a[], double b[], double c[])  // only for a[3]= b[3] x c[3];
{
    a[0]= b[1]*c[2] - b[2]*c[1];
    a[1]=-b[0]*c[2] + b[2]*c[0];
    a[2]= b[0]*c[1] - b[1]*c[0];
}

double vec_max   (double *a, int n)
{
    double s=0.;
    for( int i=0; i<n; i++ )
        s = CYCASMAX( s, a[i] );
    return s;
}
size_t vec_max_id   (double *a, int n){
     size_t s=0;
     double max = 0.;
    for( int i=0; i<n; i++ )
	    if(a[i]>max){
		    max = a[i];
		    s = i;
	    }
    return s;
}

// solution procedure : JacobiIter, SORForwIter, SORBackwIter, and SSORIter
//               ChebyshevIter, CGIter, CGNIter, GMRESIter, BiCGIter, QMRIter, CGSIter, BiCGSTABIter
// precondition : JacobiPrecond, SSORPrecond, ILUPrecond, NULL. the ILUPrecond is not what I expected

/*
void SolveLinearEqu(Vector* Func(QMatrix*, Vector*, Vector*, int,PrecondProcType, double),
					QMatrix *qa, Vector *x, Vector *b, int MaxIter, PrecondProcType PreCond, double omega,
					   double epsilon, int *Iter, double *IterRes)
{
	SetRTCAccuracy( epsilon );
	Func(qa, x, b, MaxIter, PreCond, omega);
	*Iter    = GetLastNoIter  ();
    *IterRes = GetLastAccuracy();
}

*/
void ErrorStop( string str )
{
	errorHandler.fatalRuntimeError(str);
}

char *trimwhitespace(char *str)
{
  char *end;

  // Trim leading space
  while(isspace(*str)) str++;

  if(*str == 0)  // All spaces?
    return str;

  // Trim trailing space
  end = str + strlen(str) - 1;
  while(end > str && isspace(*end)) end--;

  // Write new null terminator
  *(end+1) = 0;

  return str;
}

double ttime (void)
{
    double sec;
	sec = clock()/double(CLOCKS_PER_SEC);
    return (sec);
}

//WARNING: if string hash is used, Gmesh rid in param.in should not greeter than 9 !!!
int stringHash(const std::string& str){
	int ret = 1;
	for(size_t i=0;i!=str.size();++i)
		ret += (int)str[i] - (int)'0';
	return ret;
}

/********************************************
 * below is implement by CHENXUYI
 ********************************************/

/********************************************
 * 	DEBUG TOOL
 *******************************************/
std::ostream& operator<<(std::ostream& out,const CellData& cel){
	out<<"nface:"<<cel.nface<<std::endl;
	out<<"globalIDx "<<cel.globalIdx<<std::endl;
	out<<"face\tvertices:\tcell\n";
	for(int i=0;i!=6;++i)
		out<<cel.face[i]<<'\t'<<cel.vertices[i]<<'\t'<<cel.cell[i]<<std::endl;
	out<<"vol: "<<cel.vol<<std::endl;
	out<<"x:\t";
	for(int i=0;i!=3;++i)
		out<<cel.x[i]<<'\t';
	out<<std::endl;
	return out;

}
std::ostream& operator<<(std::ostream& out,const FaceData& fac){
	out<<"cel1: "<<	fac.cell1<<" cel2: "<<fac.cell2<<endl;
	out<<"area: "<<fac.area<<endl;
	out<<"lambda: "<<fac.lambda<<endl;
	out<<"rlencos: "<<fac.rlencos<<endl;
	out<<"x\tn"<<endl;
	for(int i=0;i!=3;++i){
		out<<fac.x[i]<<"\t"<<fac.n[i]<<endl;
	}

	out<<"xpac\txnac"<<endl;
	for(int i=0;i!=3;++i){
		out<<fac.Xpac[i]<<"\t"<<fac.Xnac[i]<<endl;
	}

	return out;
}


/********************************************
 * MPI Parallel I/O
 * collective
 * must ensure passing the buffer of same length to the function
 ********************************************/
int parallelWriteBuffer(const string& title,const void* buffer,int myBufSize,DataPartition* dg, int head){//collective
	MPI_File thefile;
	MPI_Barrier(MPI_COMM_WORLD);

	int* bufSizes = new int[dg->comSize];
	int disp = 0;
	
	MPI_Allgather(&myBufSize,1,MPI_INT,bufSizes,1,MPI_INT,dg->comm);

	for(int i=0;i!=dg->comRank;++i){
		disp+=bufSizes[i];
	}

	int ret = MPI_File_open(MPI_COMM_WORLD,(char*)title.c_str(), MPI_MODE_RDWR, MPI_INFO_NULL, &thefile);
	if(ret!=MPI_SUCCESS){ 
		throw runtime_error("MPI Parallel I/O fail: cant open file \n");
	}

	ret = MPI_File_set_view(thefile, (MPI_Offset)head*sizeof(char) + (MPI_Offset)disp * sizeof(char), MPI_CHAR,MPI_CHAR,"native",MPI_INFO_NULL);
	if(ret!=MPI_SUCCESS){ 
		throw runtime_error("MPI Parallel I/O fail: cant set view \n");
	}

	ret = MPI_File_write(thefile,(char*)buffer,myBufSize,MPI_CHAR,MPI_STATUS_IGNORE);
	if(ret!=MPI_SUCCESS){ 
		throw runtime_error("MPI Parallel I/O fail: cant write buffer \n");
	}

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_File_close(&thefile);
	delete []bufSizes;
	return 0;
}


//command line parsor
void observeCommand(map<string,bool>& cl, const string& command){
	cl.insert(make_pair(command,false));
}

bool parseCommand(map<string,bool>& cl){
	PetscBool temp = PETSC_FALSE;
	PetscBool temp2 = PETSC_FALSE;
	PetscOptionsGetBool(NULL,NULL,"-help",&temp,NULL);
	PetscOptionsGetBool(NULL,NULL,"-helpCycas",&temp2,NULL);
	if(temp==PETSC_TRUE || temp2==PETSC_TRUE){
		for(map<string,bool>::iterator iter = cl.begin();iter!=cl.end();++iter){
			PetscPrintf(MPI_COMM_WORLD,"%s\n",iter->first.c_str());
		}	
		return false;
	}
	for(map<string,bool>::iterator iter = cl.begin();iter!=cl.end();++iter){
		PetscOptionsGetBool(NULL,NULL,iter->first.c_str(),&temp,NULL);
		iter->second = temp==PETSC_TRUE;
		temp = PETSC_FALSE;
	}
	return true;
}

bool getCommand(const map<string,bool>& cl, const string& command){
	map<string,bool>::const_iterator iter = cl.find(command);	
	assert(iter!=cl.end());
	return iter->second;
}



