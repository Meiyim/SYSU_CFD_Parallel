#include<iostream>
#include<petscksp.h>
#include<limits>
#include<numeric>
#include<string>
#include<algorithm>
/********************************************************
 * define the basic type in NavierStorkesSolver
 ********************************************************/

#define CYCAS_DEBUG_MODE

#define SMALL 1.e-16
#define CYCASHUGE_D std::numeric_limits<double>::max()
#define CYCASHUGE_I std::numeric_limits<int>::max()
#define INT_OPTION_NO 16 
#define DB_OPTION_NO  23
#define TECPLOT_NVAR  13


#ifndef BASIC_TYPE_H
#define BASIC_TYPE_H

class ErrorHandler{
public:
	void fatalRuntimeError(std::string msg){
		int r;
		MPI_Comm_rank(MPI_COMM_WORLD,&r);
		printf("!!!!!!!!!!!!!!!!run time error occured in rank: %d!!!!!!!!!!!!!!!!!\n",r);
		std::cout<<msg;
		printf("!!!!!!!!!!!!!!!!System will Abort!!!!!!!!!!!!!!!!!\n");
		getchar();
		MPI_Abort(MPI_COMM_WORLD,0);
		PetscFinalize();
	}
	void fatalLogicError(std::string msg, int id){
		char temp[1024];
		sprintf(temp,"%s: %d\n",msg.c_str(),id);
		fatalLogicError(temp);
	}
	void fatalLogicError(std::string msg, const char* detal){
		char temp[1024];
		sprintf(temp,"%s: %s\n",msg.c_str(),detal);
		fatalLogicError(temp);
	}
	void fatalLogicError(std::string msg){
		int r;
		MPI_Comm_rank(MPI_COMM_WORLD,&r);
		printf("!!!!!!!!!!!!!!!!run time error occured in rank: %d!!!!!!!!!!!!!!!!!\n",r);
		std::cout<<msg<<std::endl;
		printf("!!!!!!!!!!!!!!!!System will Abort!!!!!!!!!!!!!!!!!\n");
		getchar();
		MPI_Abort(MPI_COMM_WORLD,0);
		PetscFinalize();
	}
};
extern ErrorHandler errorHandler;


/******************************************
 *	Throw Error when not Converge
 ******************************************/
class ConvergeError{
public:
	int iter;
	double residual;
	std::string varname;
	ConvergeError(int i,double r,std::string v):iter(i),residual(r),varname(v){}
};


/******************************************
 *	Navier Stokes Data Structure
 ******************************************/
class FlowField{
public:
	int nVar;
private:
	int basicSize;
	int turbSize;
	int speciesSize;
	double* const basicflowField;
	double* turbField;
	double* speciesField;
public:
	double* const R;	
	double* const U;	
	double* const V;	
	double* const W;	
	double* const P;	
	double* const T;	
	double* TE;	
	double* ED;	
	double** RS;	
	FlowField(size_t size1,bool useturb, bool useMultySpecies,size_t size2)://size2 is nSpecies
		nVar(6),
		basicSize(nVar*size1),
		basicflowField(new double[basicSize]),
		turbField(NULL),
		speciesField(NULL),
		R(basicflowField + 0*size1),
		U(basicflowField + 1*size1),
		V(basicflowField + 2*size1),
		W(basicflowField + 3*size1),
		P(basicflowField + 4*size1),
		T(basicflowField + 5*size1),
		TE(NULL),
		ED(NULL),
		RS(NULL)
	{
		if(useturb){
			nVar+=2;
			turbSize =2;
			turbField = new double[size1*turbSize];
			TE = turbField;
			ED = turbField + size1;
		}
		if(useMultySpecies){
			printf("hals \n");
			speciesSize = size2;
			speciesField = new double[speciesSize*size1];
			RS = new double* [size2];
			for(size_t i=0;i!=size2;++i){
				RS[i] = speciesField + i*size1;
			}
		}
	}
	~FlowField(){
		delete []basicflowField;
		delete []turbField;
		delete []speciesField;
		delete []RS;
	}
	void attachSpecices(double***rs){
		*rs = RS;	
	}
	void attachTurb(double** te, double**ed){
		*te = TE;
		*ed = ED;
	}
	void attachBasic(double** r,double** u,double** v,double** w,double** p,double** t){
		*r = R; 
		*u = U;
		*v = V;
		*w = W;
		*p = P;
		*t = T;
	}
};


class FaceData
{
public:
    int    bnd;
    int    vertices[4];
    int    cell1,cell2;
    double x[3],n[3],area;   // face center and normal vector
    double lambda;   // lambda for left, (1.-lambda) for right
    double rlencos;     // area/(|Xpn|*vect_cosangle(n,Xpn))
    double Xpac[3],Xnac[3]; // auxiliary points
    void* info; //some auxiliary data...
    FaceData():
	    bnd(CYCASHUGE_I),
	    cell1(CYCASHUGE_I),
	    cell2(CYCASHUGE_I),
	    area(CYCASHUGE_D),
	    lambda(CYCASHUGE_D),
	    rlencos(CYCASHUGE_D),
	    info(NULL)
	{
		for(int i=0;i!=4;++i)
			vertices[i] = CYCASHUGE_I;
		for(int i=0;i!=3;++i)
			x[i] = n[i] = Xpac[i] = Xnac[i] = CYCASHUGE_D;
	}
	void reverse(){
		std::swap(cell1,cell2);
		lambda = 1. - lambda;
		for(int i=0;i!=3;++i){
			n[i] = 0. - n[i];
			std::swap(Xpac[i],Xnac[i]);
		}
	}
};


class CellData{
public:
    int nface;
    int rid;					//point to BdRegion
    int face[6], cell[6], vertices[8]; // maybe wasterful a bit. Make it dynamics to save memory
    				       //  face[nface]: for the faces index
    				       //  all elements are treated as 8 nodes hexahdron,
    				       //  but with different number of faces, so judges will be used for avoid faces with one vertex
    int  globalIdx; 		       //added by CXY: len: Ncel, global idx for each cell
    double vol;
    double x[3];
    CellData():
	    nface(CYCASHUGE_I),
	    rid(CYCASHUGE_I),
	    globalIdx(CYCASHUGE_I),
	    vol(CYCASHUGE_D)
	    {
		    for(int i=0;i!=6;++i)
			    face[i] = cell[i]  = CYCASHUGE_I;
		    for(int i=0;i!=8;++i)
			    vertices[i] = CYCASHUGE_I;
		    for(int i=0;i!=3;++i)
			    x[i] = CYCASHUGE_D;
	    }

    CellData& operator=(CellData& rhs){
	    this->nface = rhs.nface;
	    this->rid = rhs.rid;
	    this->globalIdx = rhs.globalIdx;
	    for(int i=0;i!=6;++i){
		   this->face[i] = rhs.face[i];
		   this->cell[i] = rhs.cell[i];
	    }
	    for(int i=0;i!=8;++i)
		    this->vertices[i] = rhs.vertices[i];

	    this->vol = rhs.vol;
	    for(int i=0;i!=3;++i)
		    this->x[i] = rhs.x[i];


	    return *this;
    }

};



// connect boundary to faces, boundary to BdRegion
struct BoundaryData
{
    int face;                // belongs to face...
    int vertices[4];         // the 4 vertices, to be done allocatable
    int rid;                 // region id as set in rtable
    double distance;         // normal distance from cell face
                             // center to cell center
    double  yplus;           // y+
    double  uplus;           // u+
    double  shear[3];        // shear force component, the meaning of this attribute is modified by CXY:
    double  h;               // local heattransfer coef.
    double  q;               // local heat flux (in W/m2)
    double  T;               // local wall temperature
    void* info;	//some auxilary data...
    BoundaryData():
    yplus(0.0),uplus(0.0),h(0.0),q(0.0),T(0.0),
    info(NULL)
    {
    		for(int i=0;i!=3;++i)	
    		shear[i] = 0.;
    }
};


// store the data set in one type of boundary
class BdRegion{
public:
    std::string name;
    int  type1;  //=1: wall
                 //=2: inlet
                 //=3: outlet
                 //=4: symmetric
                 //=5: body
    			 //=6: period
    // temperature
    int  type2; // for wall:(type1==1)
                // =0 : given temperature; 
                // =1 : heat flux (0 for default adiabatic, others )
                // =2 : coupled

                // for body:(type1==5)
                // =0 : fluid field: u, v, w, ro, p, t ...
                // =1 : solid: ro, diffCoef

    			// for period:(type1 == 6)
    			// type2 = corresponding bid;


    double fixedValue;// init or bnd values
    double initvalues[10];//u,v,w,p,ro,t,te,ed

    double  scarlarRecord; // for monitor output
    double  vectorRecord[3];// 
    BdRegion():
        type1(-1),
        type2(-1),
        fixedValue(CYCASHUGE_I),
        scarlarRecord(0.0)
    {
    	vectorRecord[0] = vectorRecord[1] = vectorRecord[2] = 0.0;
    }

    BdRegion(double* buffer,int len):	 //construct from buffer
    	scarlarRecord(0.0)
    {
    	int icounter = 0;
   		type1 = (int)buffer[icounter++];
   		type2 = (int)buffer[icounter++];
   		fixedValue = buffer[icounter++];
    	for(int i=0;i!=len-3;++i){
    		initvalues[i] = buffer[icounter++];
    	}
    	vectorRecord[0] = vectorRecord[1] = vectorRecord[2] = 0.0;
    }

    int getSendBuffer(double* buffer){		//user's responsibility to alloc&free memory
    	//name field is ignored ... who cares;
    	int len = 13;
    	buffer[0] = (double) type1;
    	buffer[1] = (double) type2;
    	buffer[2] = fixedValue;
    	for(int i=0;i!=10;++i){
    		buffer[3+i] = initvalues[i];

    	}
    	return len;
    }

};


#endif


