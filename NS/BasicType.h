#include<petscksp.h>
/********************************************************
 * define the basic type in NavierStorkesSolver
 ********************************************************/

#ifndef BASIC_TYPE_H
#define BASIC_TYPE_H

struct FaceData
{
    int    bnd;
    int    vertices[4];
    int    cell1,cell2;
    double x[3],n[3],area;   // face center and normal vector
    double lambda,   // lambda for left, (1.-lambda) for right
		rlencos;     // area/(|Xpn|*vect_cosangle(n,Xpn))
	double Xpac[3],Xnac[3]; // auxiliary points
};


class CellData{
public:
    int nface;
    int face[6], cell[6], vertices[8]; // maybe wasterful a bit. Make it dynamics to save memory
    //face[nface]: for the faces index
    //  all elements are treated as 8 nodes hexahdron,
    //  but with different number of faces, so judges will be used for avoid faces with one vertex
    double vol, x[3];
    PetscInt  globalIdx; //added by CXY: len: Ncel, global idx for each cell
    CellData():
	    nface(-1),
	    vol(-1.0),
	    globalIdx(-1)
	    {
		    for(int i=0;i!=6;++i)
			    face[i] = cell[i]  = -1;
		    for(int i=0;i!=8;++i)
			    vertices[i] = -1;
		    for(int i=0;i!=3;++i)
			    x[i] = -1;
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
    double  shear[3];        // shearstress components
    double  h;               // local heattransfer coef.
    double  q;               // local heat flux (in W/m2)
    double  T;               // local wall temperature
};

// store the data set in one type of boundary
struct BdRegion
{
	char name[50];
	int  itype;
	// momentum and pressure
	double vel[3],massflow[3],veldirection[3], pressure;
	// temperature
	int    temType; // =0 : given temperature; 
	                // =1 : heat flux (0 for default adiabatic, others )
	double tem,qflux;
	// turbulence
	double te,ed, turbIntensity;
	// species
	int    speType; // = 0 : specific mass fraction (e.g., saturated concentration);
	                // = 1 : specific mass flux (0 for default,means no diffusion flux in Fluent)
	double *speCon; // mass fraction
};

#endif


