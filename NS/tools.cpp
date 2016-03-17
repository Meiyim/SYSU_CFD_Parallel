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
int parallelWriteBuffer(const string& title,const string& buffer,DataPartition* dg, int head){//collective
	MPI_File thefile;
	MPI_Barrier(MPI_COMM_WORLD);

	int myBufSize = buffer.size();
	int* bufSizes = new int[dg->comSize];
	int disp = 0;
	
	MPI_Allgather(&myBufSize,1,MPI_INT,bufSizes,1,MPI_INT,dg->comm);

	for(int i=0;i!=dg->comRank;++i){
		disp+=bufSizes[i];
	}

	int ret = MPI_File_open(MPI_COMM_WORLD,title.c_str(), MPI_MODE_RDWR, MPI_INFO_NULL, &thefile);
	if(ret!=MPI_SUCCESS){ 
		throw runtime_error("MPI Parallel I/O fail: cant open file \n");
	}

	ret = MPI_File_set_view(thefile, (MPI_Offset)head*sizeof(char) + (MPI_Offset)disp * sizeof(char), MPI_CHAR,MPI_CHAR,"native",MPI_INFO_NULL);
	if(ret!=MPI_SUCCESS){ 
		throw runtime_error("MPI Parallel I/O fail: cant set view \n");
	}

	ret = MPI_File_write(thefile,buffer.c_str(),myBufSize,MPI_CHAR,MPI_STATUS_IGNORE);
	if(ret!=MPI_SUCCESS){ 
		throw runtime_error("MPI Parallel I/O fail: cant write buffer \n");
	}

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_File_close(&thefile);
	delete []bufSizes;
	return 0;
}



