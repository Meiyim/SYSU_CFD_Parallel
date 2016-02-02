#ifndef  TOOLS_H
#define  TOOLS_H

#include <string>
#include <stdlib.h>
using namespace std;

#define CYCASMAX(x,y)  ((x)>(y)?(x):(y))
#define CYCASMIN(x,y)  ((x)<(y)?(x):(y))
#define CYCASSIGN(x)   ((x)>0?1:(-1))
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
    for (int i = 0; i < row; ++i)
        for (int j = 0; j < col; ++j)
            arr[i][j].~T();
    if (arr != NULL)
        free((void**)arr);
}
template <typename T>
void init_Array2D(T **arr, int row, int col, T val )
{
    for (int i = 0; i < row; ++i)
        for (int j = 0; j < col; ++j)
            arr[i][j] = val;
}


#endif
