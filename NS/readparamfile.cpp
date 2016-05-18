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
		errorHandler.fatalLogicError("param file not found\n check your input\n");
	}
	
	int _linecounter=0;
	while( getline( fin,str ) )
	{
		_linecounter++;
		
		// skip comment line
		if( str[0]=='/' && str[1]=='/' ) continue;
		if( str.size()==1 )continue;
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

			BdRegion& reg = regionMap[0];//inseret;
			reg.name= "default fluid";
			reg.type1 = 5;
			reg.type2 = 0;//fluid
			for( i=1; i<=ikey; i++ )
				reg.initvalues[i-1]= atof( keyw[i] );
		}

		//--- boundary condition
		else if( strcmp(keyw[0],"bound")==0 )
		{
			int bid= atoi( keyw[1] );
			if(bid == 0){
				errorHandler.fatalLogicError("default rid 0 is occupied");
			}
			if(regionMap.find(bid) != regionMap.end()){
					errorHandler.fatalLogicError("repeating rid definition",bid);
			}
			if( strcmp(keyw[2],"Twall")==0 )
			{
				if(ikey!=3){
					errorHandler.fatalLogicError("temperature wall field requied 1 parameters");
				}
				regionMap[bid].name="Twall";
				regionMap[bid].type1 = 1;
				regionMap[bid].type2 = 0;//given T
				regionMap[bid].fixedValue = atof(keyw[3]);
		
			}
			else if( strcmp(keyw[2],"Hwall")==0 ){
				if(ikey!=3){
					errorHandler.fatalLogicError("heat flux wall field requied 1 parameters");
				}
				regionMap[bid].name="Hwall";
				regionMap[bid].type1 = 1;
				regionMap[bid].type2 = 1;//given flux
				regionMap[bid].fixedValue = atof(keyw[3]);
			}
			else if( strcmp(keyw[2],"sym")==0 ){
				if(ikey!=2){
					errorHandler.fatalLogicError("symmetric field requied no parameters");
				}
				regionMap[bid].name="sym";
				regionMap[bid].type1 = 4;
				regionMap[bid].type2 = 0;
				regionMap[bid].fixedValue = 0.0;
			}
			else if( strcmp(keyw[2],"inlet")==0 ){
				if(ikey!=10){
					errorHandler.fatalLogicError("inlet field requied 8 parameters");
				}
				regionMap[bid].name="inlet";
				regionMap[bid].type1 = 2;
				regionMap[bid].type2 = 0;
				double* initvalues = regionMap[bid].initvalues;
				initvalues[0] = atof(keyw[3]);//u
				initvalues[1] = atof(keyw[4]);//v
				initvalues[2] = atof(keyw[5]);//w
				initvalues[3] = atof(keyw[6]);//p
				initvalues[4] = atof(keyw[7]);//r
				initvalues[5] = atof(keyw[8]);//t
				if( TurModel==1 ){
					initvalues[6] = atof(keyw[9]);//te
					initvalues[7] = atof(keyw[10]);//ed
				}
			}	
			else if( strcmp(keyw[2],"outlet")==0 ){
				if(ikey!=3){
					errorHandler.fatalLogicError("outlet field requied 1 parameter");
				}
				regionMap[bid].name="outlet";
				regionMap[bid].type1 = 3;
				regionMap[bid].type2 = 0;
				regionMap[bid].fixedValue = atof(keyw[3]); //pout
			}
			else if( strcmp(keyw[2],"period")==0) {
				if(ikey != 3){
					errorHandler.fatalLogicError("periodic bc requied 2 parameters");
				}
				regionMap[bid].name="period";
				regionMap[bid].type1 = 6;
				regionMap[bid].type2 = atoi(keyw[3]); //corresponding boundary bid;
				regionMap[bid].fixedValue = 0.0;
			}
			else
			{
				errorHandler.fatalLogicError("unknown boundry type in bound",keyw[2]);
			}
		}

		//volumn definition
		else if( strcmp(keyw[0],"volumn")==0 )
		{
			int bid= atoi( keyw[1] );
			if(bid == 0){
				errorHandler.fatalLogicError("default rid 0 is occupied");
			}
			if(regionMap.find(bid) != regionMap.end()){
					errorHandler.fatalLogicError("repeating rid definition",bid);
			}

			if(      strcmp(keyw[2],"fluid")==0 )//fluid
			{
				regionMap[bid].name= keyw[3];
				regionMap[bid].type1 = 5;
				regionMap[bid].type2 = 0;//fluid
				if(ikey!=10){
					errorHandler.fatalLogicError("fluid field required 8 parameters");
				}
				for(int i=4;i<=ikey;++i){
					regionMap[bid].initvalues[i-4]= atof( keyw[i] );
				}
			}
			else if( strcmp(keyw[2],"solid")==0 ){// solid
				regionMap[bid].name=keyw[3];
				regionMap[bid].type1 = 5;
				regionMap[bid].type2 = 1;//solid
				if(ikey!=5){
					errorHandler.fatalLogicError("solid field required 3 parameters");
				}
				for(int i=4;i<=ikey;++i){
					regionMap[bid].initvalues[i-4]= atof( keyw[i] );
				}
			}

			else
			{
				errorHandler.fatalLogicError("unknown boundry type in bound",keyw[2]);
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
	if(dataPartition->comRank == root.rank){
		root.regionMap = this->regionMap;
	}
}
