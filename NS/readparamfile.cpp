#include <iostream>
#include <string>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include "navier.h"
#include "petscsys.h"

void NavierStokesSolver::ReadParamFile( )
{
	int i,j, ikey;
	ifstream fin;
	char 	  _commandlinetitle[128];
	std::string*	  title = NULL;
	PetscBool flg;
	char      *keyw[20];
	std::string str;

	PetscOptionsGetString(NULL,"-param",_commandlinetitle,sizeof(_commandlinetitle),&flg);
	if(flg){
		title = new string(_commandlinetitle);
	}else{
		title = new string("param.in");
	}	

	for( i=0; i<20; i++ )
		keyw[i] = new char[30];
	
	fin.open(title->c_str());
	if( ! fin.is_open( ) ){
		cout<<"param file does not exist!"<<endl;
		exit(0);
	}
	int _linecounter=0;
	while( getline( fin,str ) )
	{
		_linecounter++;
		
		// skip comment line
		if( str[0]=='/' && str[1]=='/' ) continue;
		if( str.empty() )continue;
		// clear keyw[] before putting value
		for( i=0; i<20; i++ )
			keyw[i][0] = '\0';
		// put commands+parameters into string array
		j=0; ikey=0;
		for( i=0; str[i]!='\0'; i++ )
		{
			if( str[i]!=',' )
			{
				keyw[ikey][j]= str[i];
				j++;
			}
			else
			{
				keyw[ikey][j]='\0';
				ikey++;
				j=0;
			}
		}
		keyw[ikey][j]='\0';
		
		// trim spaces in keywords
		for( i=0; i<=ikey; i++ )
		{
			keyw[i] = trimwhitespace( keyw[i] );
		}


		// analyze command and parameters
		//---- solver parameters
		if(      strcmp(keyw[0],"gridfile")==0 )
		{
			strcpy( GridFileName, keyw[1] );
		}
		else if( strcmp(keyw[0],"steady")==0 )
		{
			IfSteady = true;
			MaxStep  = atoi(keyw[1]);
			ResidualSteady= atof(keyw[2]);
		}
		else if( strcmp(keyw[0],"transient")==0 )
		{
			IfSteady = false;
			dt = atof(keyw[1]);
			total_time= atof(keyw[2]);
			if(      strcmp(keyw[3],"Euler")==0 )
				TimeScheme = 1;
			else if( strcmp(keyw[3],"Dual")==0 )
				TimeScheme = 2;
			else
			{
				char temp[256];
				sprintf(temp,"Error in raading param.in\tline: %d\n unknown time scheme: %s\n",_linecounter,keyw[3]);
				errorHandler.fatalLogicError(temp);
			}
		}
		else if( strcmp(keyw[0],"energy")==0 )
		{
			if(      strcmp(keyw[1],"on")==0 )
				SolveEnergy= true;
			else if( strcmp(keyw[1],"off")==0 )
				SolveEnergy= false;
			else{
				char temp[256];
				sprintf(temp,"Error in raading param.in\tline: %d\n unknown energy model\n",_linecounter);
				errorHandler.fatalLogicError(temp);
			}
		}
		else if( strcmp(keyw[0],"density")==0 )
		{
			DensityModel = atoi(keyw[1]);
		}
		else if( strcmp(keyw[0],"turbulence")==0 )
		{
			if(      strcmp(keyw[1],"no")==0 )
				TurModel= 0;
			else if( strcmp(keyw[1],"ke")==0 )
				TurModel= 1;
			else{
				char temp[256];
				sprintf(temp,"Error in raading param.in\tline: %d\n unknown turbulence model\n",_linecounter);
				errorHandler.fatalLogicError(temp);
			}
		}
		else if( strcmp(keyw[0],"gravity")==0 )
		{
			gravity[0] = atof( keyw[1] );
			gravity[1] = atof( keyw[2] );
			gravity[2] = atof( keyw[3] );
		}
		else if( strcmp(keyw[0],"PressureRef")==0 )
		{
			PressureReference = atof(keyw[1]);
			if(strlen(keyw[2])>0) cellPressureRef   = atoi(keyw[2])-1;
		}

		//--- numerical scheme
		else if( strcmp(keyw[0],"relaxation")==0 )
		{
			URF[0] = URF[1] = URF[2] = atof( keyw[1] );
			URF[3] = atof( keyw[2] );
			URF[4] = atof( keyw[3] );
		}
		else if( strcmp(keyw[0],"limiter")==0 )
		{
			if(      strcmp(keyw[1],"no")==0 )
				limiter = 0;
			else if( strcmp(keyw[1],"Barth" )==0 )
				limiter = 1;
			else if( strcmp(keyw[1],"MLP" )==0 )
				limiter = 2;
			else if( strcmp(keyw[1],"WENO" )==0 )
				limiter = 3;
			else{
				char temp[256];
				sprintf(temp,"Error in raading param.in\tline: %d\n unknown limiter definition\n",_linecounter);
				errorHandler.fatalLogicError(temp);
			}
		}


		//--- initial condition
		else if( strcmp(keyw[0],"restart")==0 )
		{
			if(      strcmp(keyw[1],"yes")==0 )
				IfReadBackup = true;
			else if( strcmp(keyw[1],"no" )==0 )
				IfReadBackup = false;
			else{
				char temp[256];
				sprintf(temp,"Error in raading param.in\tline: %d\n unknown restart command\n",_linecounter);
				errorHandler.fatalLogicError(temp);
			}
		}
		else if( strcmp(keyw[0],"initflow")==0 )
		{
			for( i=1; i<=ikey; i++ )
				initvalues[i-1]= atof( keyw[i] );
		}

		//--- boundary condition
		else if( strcmp(keyw[0],"bound")==0 )
		{
			int bid= atoi( keyw[1] );
			if(      strcmp(keyw[2],"wall")==0 )
			{
				rtable[bid] = 1;
				if( SolveEnergy ) Twall = atof( keyw[3] );
			}
			else if( strcmp(keyw[2],"inlet")==0 ){
				rtable[bid] = 2;
				uin = atof(keyw[3]);
				vin = atof(keyw[4]);
				win = atof(keyw[5]);
				pin = atof(keyw[6]);
				roin= atof(keyw[7]);
				Tin = atof(keyw[8]);
				if( TurModel==1 ){
					tein = atof(keyw[9]);
					edin = atof(keyw[10]);
				}
			}
			else if( strcmp(keyw[2],"outlet")==0 ){
				rtable[bid] = 3;
				pout = atof(keyw[3]);
			}
			else if( strcmp(keyw[2],"sym")==0 )
				rtable[bid] = 4;
			else
			{
				char temp[256];
				sprintf(temp,"Error in reading param.in\t line: %d \nunknown boundary type:%s\n",_linecounter,keyw[2]);
				errorHandler.fatalLogicError(temp);
			}
		}

		//--- post process 
		else if( strcmp(keyw[0],"output")==0 )
		{
			noutput = atoi( keyw[1] );
			if(      strcmp(keyw[2],"vtk")==0 )
				outputFormat = 1;
			else if( strcmp(keyw[2],"tecplot")==0 )
				outputFormat = 0;
			else
				outputFormat = 0;
		}
		else
		{
			char temp[256];
			sprintf(temp,"Error in reading param.in\t line: %d \nno such command: %s\n",_linecounter,keyw[0]);
			errorHandler.fatalLogicError(temp);
		}

	}

	if( DensityModel == 1 )
		SolveEnergy = true;

	fin.close( );
	delete title;
}
