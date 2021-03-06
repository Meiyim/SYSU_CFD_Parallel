#ifndef  TOOLS_H
#define  TOOLS_H

#include <time.h>
#include <string>
#include <string.h>
#include <stdlib.h>
#include "MPIStructure.h"
using namespace std;

#define CYCASMAX(x,y)  ((x)>(y)?(x):(y))
#define CYCASMIN(x,y)  ((x)<(y)?(x):(y))
#define CYCASSIGN(x)   ((x)>0?1:(-1))
#define CYCAS_GET_TIME(ts) gettimeofday(&ts,NULL);
#define CYCAS_GET_ELAPSE_TIME(t1,t2) ((double)(t1.tv_sec - t2.tv_sec) - (t1.tv_usec-t2.tv_usec)/1.e6)

// maybe SIGN still has problem, not for 0
// vector manipulation. should be defined as inline function for higher efficiency
void   vec_init  (double[], int, double );
void   vec_minus (double *x1, double *x2, double *x3, int n); // x1= x2 - x3
double vec_dot   (double[], double[], int);
double vec_len   (double[], int );
void   vec_cross (double[], double[], double[]); // only for C[3]= A[3] x B[3];
double vec_max   (double[], int );


/*
void SolveLinearEqu( Vector* Func(QMatrix*, Vector*, Vector*, int,PrecondProcType, double),
			QMatrix *qa, Vector *x, Vector *b, int MaxIter, PrecondProcType PreCond, double omega,
					   double epsilon, int *Iter, double *IterRes);
*/
void ErrorStop( string str );
char *trimwhitespace(char *str);
double ttime (void);
int stringHash(const std::string& str);


template <typename T>
T** new_Array2D(int row, int col)
{
    int size = sizeof(T);
    int point_size = sizeof(T*);
    T **arr = (T **) malloc(point_size * row + size * row * col);
    if (arr != NULL)
    {   
        T *head = (T*)(arr + row); //made some change
        for (int i = 0; i < row; ++i)
        {
            arr[i] = &head[i*col]; //made some change
            for (int j = 0; j < col; ++j)
                new (&arr[i][j]) T;
        }
    }
    return (T**)arr;
}
template <typename T>
void delete_Array2D(T **arr, int row, int col)
{
	if(arr==NULL) return;
    for (int i = 0; i < row; ++i)
        for (int j = 0; j < col; ++j)
            arr[i][j].~T();
    free((void**)arr);
}
template <typename T>
void init_Array2D(T **arr, int row, int col, T val )
{
    for (int i = 0; i < row; ++i)
        for (int j = 0; j < col; ++j)
            arr[i][j] = val;
}


/********************************************
 * below is implement by CHENXUYI
 ********************************************/

//------for PRINT_LOG----------------------
//
//-----------------------------------------
std::ostream& operator<<(std::ostream&,const CellData& cel);
std::ostream& operator<<(std::ostream&,const FaceData& fac);

class Checker{
	string varname;
	vector<double> pool;
public:
	Checker(string n):varname(n){}
	void check(double v){
		pool.push_back(v);
	}
	void report(){
		double ret = 0.0;
		double reduceResult = 0.0;
		int count = pool.size();
		int countRes = 0;
		for(vector<double>::iterator it=pool.begin();it!=pool.end();++it)
			ret += (*it);
		MPI_Reduce(&ret,&reduceResult,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);	
		MPI_Reduce(&count,&countRes,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);	
		PetscPrintf(MPI_COMM_WORLD,"%s:\t%e\tcount:%d\n",varname.c_str(),reduceResult,countRes);

	}
};


#define CHECK_ARRAY(arr,len) checkArray(arr,len,#arr)
#define CHECK_MEMBER_ARRAY(arr,mem,len) checkMemberArray(arr,len,mem,#arr #mem)
template<typename T>
void checkArray(T* arr, size_t len,char* name){
	double ret = 0.0;
	double reduceResult = 0.0;
	for(int i=0;i!=len;++i){
		ret += (arr[i])*(arr[i]);
	}

	MPI_Reduce(&ret,&reduceResult,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);	
	
	char temp[256];
	sprintf(temp,"%s: norm %30.9e\n",name,sqrt(reduceResult));
	PetscPrintf(MPI_COMM_WORLD,"%s",temp);
}

template<typename T>
void checkMemberArray(T* arr, size_t len, double T::* m_ptr,char* name){
	double ret = 0.0;
	double reduceResult = 0.0;
	for(int i=0;i!=len;++i){
		ret += (arr[i].*m_ptr) * (arr[i].*m_ptr);

	}

	MPI_Reduce(&ret,&reduceResult,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);	
	
	char temp[256];
	sprintf(temp,"%s: norm %30.9e\n",name,sqrt(reduceResult));
	PetscPrintf(MPI_COMM_WORLD,"%s",temp);
}


/********************************************
 * MPI Parallel I/O
 * collective
 * must ensure passing the buffer of same length to the function
 ********************************************/
int parallelWriteBuffer(const string& title,const void* buffer,int myBufSize,DataPartition* dg, int head);


//Commandline Parsor
//accepte bool command line only
void observeCommand(map<string,bool>& cl, const string& command);
bool parseCommand(map<string,bool>& cl);
bool getCommand(const map<string,bool>& cl, const string& command);


#endif
