#include <iostream>
#include <string.h>
#include <stdlib.h>
#include <fstream>
#include "material.h"
#include "navier.h"
using namespace std;

class Material *MatSp;

int MaterialInit( )
{
    int    i,j,k,m, ntotal, order[100];
    char   DB_Names [100][50];
    double DB_CvCoef[100][5], DB_DijCoef[100][100][5], DB_VisCoef[100][5], DB_ThermalCoef[100][5];

    // read in all database file, linear one
	ifstream inf;
	inf.open("glibc.dat");
	inf>>ntotal;
	for( i=0; i<ntotal; i++ )
	{
		inf>> order[i];
		for(j=0;j<order[i];j++) inf>>DB_CvCoef     [i][j];
		for(j=0;j<order[i];j++) inf>>DB_VisCoef    [i][j];
		for(j=0;j<order[i];j++) inf>>DB_ThermalCoef[i][j];
		
		for( j=0;j<ntotal-1;j++ )
		{
			for( k=0; k<order[i]; k++ )
				inf>> DB_DijCoef[i][j][k];
		}
	}
	inf.close( );


    // search for every species and set their properties
	for( i=0; i<NSSolver.Nspecies; i++ )
    {
        for( j=0; j<ntotal; j++ )
            if( strcmp(DB_Names[j],MatSp[i].name)==1 ) break;

        if(j>=ntotal){
            cout<<"cannot find "<<MatSp[i].name<<" in database"<<endl;
            exit(0);
        }

		MatSp[i].order = order[i];
		for( k=0; k<MatSp[i].order; k++ ){
			MatSp[i].CvCoef     [k]= DB_CvCoef     [j][k];
			MatSp[i].VisCoef    [k]= DB_VisCoef    [j][k];
			MatSp[i].ThermalCoef[k]= DB_ThermalCoef[j][k];

			for( m=0; m<ntotal-1; j++ )
				MatSp[i].DijCoef[m][k]= DB_DijCoef  [j][m][k];
		}
    }
	return 0;
}


double Material::CalCv(double Tem)
{
	double ff;
	switch( order )
	{
	case(0):
		ff= CvCoef[0];
	case(1):
		ff= CvCoef[0] + Tem* CvCoef[1];
	case(2):
		ff= CvCoef[0] + Tem*(CvCoef[1]+Tem* CvCoef[2]);
	case(3):
		ff= CvCoef[0] + Tem*(CvCoef[1]+Tem*(CvCoef[2]+Tem*CvCoef[3]));
	case(4):
		ff= CvCoef[0] + Tem*(CvCoef[1]+Tem*(CvCoef[2]+Tem*(CvCoef[3]+Tem*CvCoef[4])));
	default:
		ErrorStop("no such order.");
	}
	return ff;
}

double Material::CalVis(double Tem)
{
	double ff;
	switch( order )
	{
	case(0):
		ff= VisCoef[0];
	case(1):
		ff= VisCoef[0] + Tem* VisCoef[1];
	case(2):
		ff= VisCoef[0] + Tem*(VisCoef[1]+Tem* VisCoef[2]);
	case(3):
		ff= VisCoef[0] + Tem*(VisCoef[1]+Tem*(VisCoef[2]+Tem* VisCoef[3]));
	case(4):
		ff= VisCoef[0] + Tem*(VisCoef[1]+Tem*(VisCoef[2]+Tem*(VisCoef[3]+Tem*VisCoef[4])));
	default:
		ErrorStop("no such order.");
	}
	return ff;
}

double Material::CalTherCoef(double Tem)
{
	double ff;
	switch( order )
	{
	case(0):
		ff= ThermalCoef[0];
	case(1):
		ff= ThermalCoef[0] + Tem* ThermalCoef[1];
	case(2):
		ff= ThermalCoef[0] + Tem*(ThermalCoef[1]+Tem* ThermalCoef[2]);
	case(3):
		ff= ThermalCoef[0] + Tem*(ThermalCoef[1]+Tem*(ThermalCoef[2]+Tem* ThermalCoef[3]));
	case(4):
		ff= ThermalCoef[0] + Tem*(ThermalCoef[1]+Tem*(ThermalCoef[2]+Tem*(ThermalCoef[3]+Tem*ThermalCoef[4])));
	default:
		ErrorStop("no such order.");
	}
	return ff;
}

// diffusion coefficients
double Material::CalBinDiff (int type, double Tem)
{
	double ff,*v;
	v = DijCoef[type];
	switch( order )
	{
	case(0):
		ff= v[0];
	case(1):
		ff= v[0] + Tem* v[1];
	case(2):
		ff= v[0] + Tem*(v[1]+Tem* v[2]);
	case(3):
		ff= v[0] + Tem*(v[1]+Tem*(v[2]+Tem* v[3]));
	case(4):
		ff= v[0] + Tem*(v[1]+Tem*(v[2]+Tem*(v[3]+Tem*v[4])));
	default:
		ErrorStop("no such order.");
	}
	return ff;
}