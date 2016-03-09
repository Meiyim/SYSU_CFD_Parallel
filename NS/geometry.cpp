#include <iostream>
#include <string>
#include <math.h>
#include <stdexcept>
#include <cstdlib>
#include <fstream>
#include <cassert>
#include <iterator>
#include <algorithm>
#include "navier.h"
#include "tools.h"

using namespace std;

/*******************************************************/
//	read msh file, call Metis to partition, find boundary and connectivity of partitions
//	1. read msh file
//	2. partition with Metis
//	3. reorder cell
//	4. create nodes pool and interface info
//	
//	root ONLY
/*******************************************************/
void NavierStokesSolver::readAndPartition(){
	root.init(dataPartition);
	root.read(dataPartition,string(GridFileName));
	root.partition(dataPartition, dataPartition->comSize);
	//root.reorder
}



/*******************************************************/
//	original msh reading function
//	read geometry from MPI_Buffer
//	local
//	last parameter is ouput value, containing the boundaryInfo
/*******************************************************/
int NavierStokesSolver::ReadGridFile(int* elementBuffer,double* vertexBuffer,int* interfaceBuffer,map<int,set<int> >* interfaceNodes)// free buffer after use
{
	string   line;
	int      NumNodeInCell[8]={0,2,3,4,4,8,6,5},vertices[8],
		i, elem_type,ntags,p,tag[10];
		// bndType[10]={4,2,3,1,4,4,0}; // bndType change the tag in gmsh file to navier_bc types ( 4 types currently )

	MPI_Barrier(dataPartition->comm);
	PetscPrintf(dataPartition->comm,"begin parsing grid file buffer\n");

	writeGeometryBackup(elementBuffer,vertexBuffer,interfaceBuffer);//back up binary

	/******************************************
	 *	Read Vertices
	 ******************************************/
    	assert (Nvrt > 0);
	Vert = new_Array2D<double>(Nvrt,3);

	dataPartition->PRINT_LOG(Nvrt);
	/*
	for(int i=0;i!=3*Nvrt;++i){
		dataPartition->PRINT_LOG(vertexBuffer[i]);
	}
	*/

	size_t counter = 0;
	for(i=0; i<Nvrt; i++){
             Vert[i][0] = vertexBuffer[counter++];
             Vert[i][1] = vertexBuffer[counter++];
             Vert[i][2] = vertexBuffer[counter++];
    	}
	delete []vertexBuffer;

	/******************************************
	 *	Read Cells
	 ******************************************/
	Bnd = new BoundaryData[Nbnd];	//pre_allocation
	Cell = new CellData[Ncel];
	//printf("read cell: %d and %d\n",Nbnd,Ncel);

   	counter = 0; 
	int icel = 0;
	int ibnd = 0;
	for(i=0; i!=dataPartition->gridList[dataPartition->comRank] ; ++i)
	{

		elem_type = elementBuffer[counter++];
		ntags = elementBuffer[counter++]; //temperary now fixed
		for( p=0; p<ntags; p++ )    // Dummy tags, maybe useful later
			tag[p] = elementBuffer[counter++];
		for( p=0; p<NumNodeInCell[elem_type];++p){
			vertices[p] = elementBuffer[counter++];
			if(p>=Nvrt) throw logic_error("vertex list out of range! bug may be at global -> local transfer part\n");
		}
		/*
		if(dataPartition->comRank == 3){
			printf("element read: type: %d, ntag %d, vertices:",elem_type,ntags);
			for( p=0; p<NumNodeInCell[elem_type];++p){
				printf("%d, ",vertices[p]);
			}
			printf("\n");
		}
		*/

	//-----------//

        // Some gmsh files have 2 and some have 3 tags
        	assert( ntags==2 || ntags==3 );
  
	        if( elem_type==2 || elem_type==3 ){ // Boundary, Tri/Quad face
			//Bnd = (BoundaryData*) realloc(Bnd,(Nbnd+1)*sizeof(BoundaryData)); //deprecated : pool performace!
			Bnd[ibnd].rid = tag[0];  // bndType[tag[0]];      // First tag is face type

			if(elem_type==2 ){
       		         	Bnd[ibnd].vertices[0]= vertices[0];
		                Bnd[ibnd].vertices[1]= vertices[1];
       			        Bnd[ibnd].vertices[2]= vertices[2];
		                Bnd[ibnd].vertices[3]= vertices[2];
			}else if(elem_type==3 ){
       	        		Bnd[ibnd].vertices[0]= vertices[0];
		                Bnd[ibnd].vertices[1]= vertices[1];
       		        	Bnd[ibnd].vertices[2]= vertices[2];
       		         	Bnd[ibnd].vertices[3]= vertices[3];
	       		 }else{
				throw logic_error("unknown boundary type in local grid file\n");
			}
			ibnd ++ ;
       		 }else if( elem_type==4 || elem_type==5 ||  // Tetrahedral/Hexa/Prism/Pyramid cell
       		          elem_type==6 || elem_type==7 )
       		 {
            		//Cell = (CellData*) realloc(Cell,(Ncel+1)*sizeof(CellData));//deprecated: pool performance!
            		// set all the cell types to 8-vertex hexa
            		if(     elem_type==4 )   // tetra
            		{
				// 0-3
				for( p=0; p<3; p++ )
					Cell[icel].vertices[p]= vertices[p];
				Cell[icel].vertices[3]= vertices[2];
				// 4-7
                		for( p=4; p<8; p++ )
					Cell[icel].vertices[p]= vertices[3];

	            	}else if( elem_type==5 )  // hexa
			{
				for( p=0; p<8; p++ )
				Cell[icel].vertices[p]= vertices[p];
			}else if( elem_type==6 )  // prism
       		     	{
				// 0-3
				for( p=0; p<3; p++ )
       	        			 Cell[icel].vertices[p]= vertices[p];
				Cell[icel].vertices[3]= vertices[2];
				// 4-7
       		         	for( p=4; p<7; p++ )
       		         		Cell[icel].vertices[p]= vertices[p-1];
	       		         Cell[icel].vertices[7]= vertices[5];
            		}
            		else if( elem_type==7 )  // pyramid
            		{
				// 0-3
				for( p=0; p<4; p++ )
					Cell[icel].vertices[p]= vertices[p];
				// 4-7
               			for( p=4; p<8; p++ )
               		 		Cell[icel].vertices[p]= vertices[4];
	       		 }
		   	icel++ ;
       		}else{
			throw logic_error("unknown element type in local grid file\n");
        	}

	}
	delete[] elementBuffer;

	/******************************************
	 *	Read Interface Info
	 ******************************************/
	counter = 0;
	int _nInterfaces = interfaceBuffer[counter++];	
	for(int i=0;i!=_nInterfaces;++i){
		int _interfaceID = interfaceBuffer[counter++];
		int _width = interfaceBuffer[counter++];
		set<int> _nodesSet;
		for(int _inode = 0; _inode!=_width; ++_inode){
			_nodesSet.insert(interfaceBuffer[counter++]);
		}
		interfaceNodes->insert(make_pair<int,set<int> >(_interfaceID,_nodesSet));
	}		

	//printf("rank: %d complete reading geo: icell %d, ncel %d\tibnd %d,nbnd %d\t,nvert%d \n",dataPartition->comRank,icel,Ncel,ibnd,Nbnd,Nvrt);


	OutputGrid(interfaceNodes); //output tecplot: grid.dat
	PetscPrintf(dataPartition->comm,"complete parsing grid file buffer\n");
	MPI_Barrier(dataPartition->comm);
    	return 0;
}

// output grid to tecplot for validation
void NavierStokesSolver::OutputGrid(map<int,set<int> >* boundinfo)
{
	int i;
	string title("grid.dat");
	ofstream of;
	int head =0;
	
	of.open(title.c_str());
	string sbuffer;
	char cbuffer[256];

	sprintf(cbuffer,"variables=\"x\",\"y\",\"z\",\"i\"\n");
	sbuffer = cbuffer;
	head = sbuffer.size();
	if(dataPartition->comRank == root.rank){//print head in root
		of<<sbuffer;
	}

	sbuffer="";
	sprintf(cbuffer,"zone n=%d, e=%d, f=fepoint, et=BRICK\n",Nvrt,Ncel);
	sbuffer=cbuffer;

	for( i=0; i<Nvrt; i++ ){
		double _iinfo = 0.0;
		for(map<int,set<int> >::iterator it=boundinfo->begin(); it!=boundinfo->end(); ++it){
			set<int>& _bset = it->second;
			if( _bset.find(i) != _bset.end() ){
				_iinfo = it->first + 1;
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

	/*
	char temp[256];
	sprintf(temp,"gridtemp%0d.dat",dataPartition->comRank);
	ofstream _of(temp);
	_of<<sbuffer;
	_of.close();

	// output cells
	of<<"zone n="<<Nvrt<<", e="<<Ncel<<", f=fepoint, et=BRICK"<<endl;
	for( i=0; i<Nvrt; i++ )
		of<<Vert[i][0]<<" "<<Vert[i][1]<<" "<<Vert[i][2]<<endl;
	for( i=0; i<Ncel; i++ ){
		for(j=0;j<8;j++)
			of<<Cell[i].vertices[j]+1<<" ";
		of<<endl;
	}
	*/

	/* // output boundary
	of<<"zone n="<<Nvrt<<", e="<<Ncel<<", f=fepoint, et=Quadrilateral"<<
		"VARSHARELIST = ([1-3]=1)"<<endl;
	for( i=0; i<Nbnd; i++ )
	{
	} 
	*/

	of.close( );
	PetscPrintf(dataPartition->comm,"complete writing grid.dat\n");
}
#define VOID_CELL_ON_BOUNDARY -10
#define VOID_CELL_ON_INTERFACE -20

int NavierStokesSolver::CreateFaces( )
{
    int i,j,id,n1,n2,n3,n4,n5,n6,n7,n8;
    int *NumNodeFace, **NodeFace;
    
    NumNodeFace = new int[Nvrt];
    NodeFace    = new_Array2D<int>( Nvrt,500 ); //CXY: not a good idea
	for( i=0; i<Nvrt; i++ )
		NumNodeFace[i] = 0;
	for( i=0; i<Ncel; i++ )
		Cell[i].nface = 0;
	/* //CXY: seems not useful
	for( i=0; i<Nfac; i++ ){
		Face[i].bnd = -10;
		Face[i].lambda = 1.;
	}
	*/
    
	PetscPrintf(dataPartition->comm,"begin faces construction\n");
    	//Nfac = 0;
   	// first boundary faces
   	 
	Nfac = 0;
	Face = new FaceData[Nbnd];	// modified by CXY

   	 for( i=0; i<Nbnd; i++ )
    	{
		//Face = (FaceData*)realloc( Face, (Nfac+1)*sizeof(FaceData) ); //deleted! bad performance
      	 	Bnd[i].face = Nfac;
      	 	for( j=0;j<4;j++ )
       			Face[Nfac].vertices[j] = Bnd[i].vertices[j];
	        Face[Nfac].bnd     = i   ;  // bnd(i)%rid CXY: Face.bnd is the idx of Bnd
       		Face[Nfac].cell2   = VOID_CELL_ON_BOUNDARY;  // cell1 is not known, set in FindFace
        	for( j=0;j<4;j++ )
        	{
            		id = Face[Nfac].vertices[j];
            		NodeFace[id][NumNodeFace[id]++]= Nfac;
        	}
		Nfac++;  // face accumulation is at last
    	}

	if(Nfac!=Nbnd) throw logic_error("error occur finding faces\n");
    	// interior faces
    	for( i=0;i<Ncel;i++ )
    	{
        	n1= Cell[i].vertices[0];
        	n2= Cell[i].vertices[1];
        	n3= Cell[i].vertices[2];
        	n4= Cell[i].vertices[3];
        	n5= Cell[i].vertices[4];
        	n6= Cell[i].vertices[5];
        	n7= Cell[i].vertices[6];
        	n8= Cell[i].vertices[7];
        	FindFace( i,n1,n2,n3,n4, Nfac, NumNodeFace,NodeFace );//Nface++
        	FindFace( i,n5,n6,n7,n8, Nfac, NumNodeFace,NodeFace );//Nface++
        	FindFace( i,n1,n2,n6,n5, Nfac, NumNodeFace,NodeFace );//Nface++
        	FindFace( i,n2,n3,n7,n6, Nfac, NumNodeFace,NodeFace );//Nface++
        	FindFace( i,n4,n3,n7,n8, Nfac, NumNodeFace,NodeFace );//Nface++
        	FindFace( i,n1,n4,n8,n5, Nfac, NumNodeFace,NodeFace );//Nface++
    	}

	printf("Rank: %d Number of faces: %d\n",dataPartition->comRank,Nfac);

	delete [] NumNodeFace;
	delete_Array2D<int>( NodeFace, Nvrt, 500 );
    	PetscPrintf(dataPartition->comm,"finish face construction\n");

	return 0;
}

void NavierStokesSolver::FindFace( int ic, int n1, int n2, int n3, int n4, int &nf, 
        int *NumNodeface, int **NodeFace )
{
    int fid,k,irepeat,iface,ifn1,ifn2,ifn3,ifn4,id;
    irepeat = 0;
    if( n1==n2 || n1==n3 || n1==n4 ) irepeat++;
    if( n2==n3 || n2==n4           ) irepeat++;
    if( n3==n4                     ) irepeat++;
    if( irepeat>1 ) return;
    // find 
    fid= -100;
    for( k=0; k<NumNodeface[n1]; k++ )
    {
        iface= NodeFace[n1][k];
        ifn1 = Face[iface].vertices[0];
        ifn2 = Face[iface].vertices[1];
        ifn3 = Face[iface].vertices[2];
        ifn4 = Face[iface].vertices[3];
        if( (ifn1==n1 || ifn1==n2 || ifn1==n3 || ifn1==n4 ) &&
            (ifn2==n1 || ifn2==n2 || ifn2==n3 || ifn2==n4 ) &&
            (ifn3==n1 || ifn3==n2 || ifn3==n3 || ifn3==n4 ) &&
            (ifn4==n1 || ifn4==n2 || ifn4==n3 || ifn4==n4 ) )
        {
            fid= iface;
            break;
        }
    }

    if( fid<0 )
    {
	Face = (FaceData*)realloc( Face, (Nfac+1)*sizeof(FaceData) );//CXY: not a good idea
        Face[nf].vertices[0]= n1;
        Face[nf].vertices[1]= n2;
        Face[nf].vertices[2]= n3;
        Face[nf].vertices[3]= n4;
        Face[nf].cell1      = ic;   
        Face[nf].cell2      = VOID_CELL_ON_INTERFACE; //added by CXY: cell2 might not be set in the following FindFace,
       						      //thus results in a interface cell
        Face[nf].bnd        = -10;  // inner face
        Cell[ic].face[ Cell[ic].nface++ ] = nf;
        // add face to NodeFace
        for( k=0; k<4; k++ )// CXY: is this loop really necessary?
        {
            id = Face[nf].vertices[k];
            NodeFace[id][ NumNodeface[id]++ ] = nf;
        }
	nf++; // face accumulation is at last
    }
    else
    {
        if( Face[fid].bnd >= 0 )    // boundary face
            Face[fid].cell1  = ic;  // cell2 == -10
        else                       // inner face
            Face[fid].cell2  = ic;
        
        Cell[ic].face[ Cell[ic].nface++ ] = fid;
    }
}


/////////////////////////////////////////////////////
/// Other cell and face information
///--------------------------------------------------

double Tri3DArea(double x1[], double x2[], double x3[] )
{
    double a1,a2,a3, v1[3],v2[3];
    vec_minus(v1, x2,x1, 3 );
    vec_minus(v2, x3,x1, 3 );
    a1=  v1[1]*v2[2] - v1[2]*v2[1];
    a2= -v1[0]*v2[2] + v1[2]*v2[0];
    a3=  v1[0]*v2[1] - v1[1]*v2[0];
    return 0.5 * sqrt( a1*a1 + a2*a2 + a3*a3 );
}
double TetVolum(double x1[], double x2[], double x3[], double x4[])
{
    double v1[3], dx1[3],dx2[3],dx3[3];
    vec_minus( dx1, x2, x1, 3);
    vec_minus( dx2, x3, x1, 3);
    vec_minus( dx3, x4, x1, 3);
    vec_cross( v1,  dx1,dx2);
    return 1./6.*fabs( vec_dot( v1, dx3, 3 ) );
}

/************************************************
 *	modified by CXY:
 *	input : interfaceNodes, the nodesPool of interfaces
 *	return value: indicates the number of virtual cell;
 ************************************************/
int NavierStokesSolver::CellFaceInfo(map<int,set<int> >* interfaceNodes)
{
    int i,j, n1,n2,n3,n4, ic1,ic2,ic, iface,nc,jf,jc;
    double area1,area2, r1,r2, 
       xc1[3],xc2[3],xc3[3],xc4[3],xc5[3],xc6[3], dx[3],dx1[3],dx2[3], 
       xv1[3],xv2[3],xv3[3],xv4[3],xv5[3],xv6[3],xv7[3],xv8[3], 
       V1,V2,V3,V4,V5,V6;

    // cell ceneter estimation, not final ones
    for( i=0; i<Ncel; i++ )
    {
        for( j=0; j<3; j++ )
        Cell[i].x[j]=1./8.*( Vert[ Cell[i].vertices[0] ][j] + 
                             Vert[ Cell[i].vertices[1] ][j] + 
                             Vert[ Cell[i].vertices[2] ][j] + 
                             Vert[ Cell[i].vertices[3] ][j] + 
                             Vert[ Cell[i].vertices[4] ][j] + 
                             Vert[ Cell[i].vertices[5] ][j] + 
                             Vert[ Cell[i].vertices[6] ][j] +
                             Vert[ Cell[i].vertices[7] ][j] );
    }

    //-- face
    for( i=0; i<Nfac; i++ )
    {
        n1= Face[i].vertices[0];
        n2= Face[i].vertices[1];
        n3= Face[i].vertices[2];
        n4= Face[i].vertices[3];

        // x, area, n, lambda
        // face center, (gravity center???) I doubt that ???
        area1 = Tri3DArea( Vert[n1], Vert[n2], Vert[n4] );
        area2 = Tri3DArea( Vert[n2], Vert[n3], Vert[n4] );

		// face center is bary-center
        for(j=0;j<3;j++){
        xc1[j]= 1./3*( Vert[n1][j]+Vert[n2][j]+Vert[n4][j] );
        xc2[j]= 1./3*( Vert[n2][j]+Vert[n3][j]+Vert[n4][j] );
        }
        for(j=0;j<3;j++)
            Face[i].x[j] = ( area1*xc1[j] + area2*xc2[j] ) / (area1+area2);
	
		// face center is bary-center
		/*  for(j=0;j<3;j++)
		{
			if( n3==n4 ) // triangle
				Face[i].x[j]=1./3.*(Vert[n1][j]+Vert[n2][j]+Vert[n3][j] );
			else         // quadrilateral
				Face[i].x[j]=1./4.*(Vert[n1][j]+Vert[n2][j]+Vert[n3][j]+Vert[n4][j]);
		}  */

        // face normal vector and area
        vec_minus( dx1, Vert[n1], Vert[n3], 3 ); // dx1= vert[n1] - vert[n3]
        vec_minus( dx2, Vert[n2], Vert[n4], 3 ); // dx2= vert[n2] - vert[n4]
        vec_cross( Face[i].n, dx1, dx2 );     // face.n is normal area
        for( j=0;j<3;j++ )
            Face[i].n[j] = 0.5 * Face[i].n[j];
        ic = Face[i].cell1;
        vec_minus( dx, Face[i].x, Cell[ic].x, 3);
        if( vec_dot(Face[i].n, dx,3)<0. ){   // normal points from cell1 to cell2
			for( j=0;j<3;j++ )
				Face[i].n[j] =  -Face[i].n[j] ;
        }
        Face[i].area= vec_len( Face[i].n,3 );
        if( fabs( (area1+area2-Face[i].area)/Face[i].area )>1.e-4 ){
            cout<<i<<" area is not correct "<<area1+area2<<" "<<Face[i].area<<endl;
			cout<<Face[i].x[0]<<" "<<Face[i].x[1]<<" "<<Face[i].x[2]<<endl;
			cout<<Face[i].vertices[0]<<" "<<Face[i].vertices[1]<<" "
			    <<Face[i].vertices[2]<<" "<<Face[i].vertices[3]<<endl;
            exit(0);
        }

        // face interpolation
        ic1= Face[i].cell1;
        ic2= Face[i].cell2;
		if( ic2>=0 ){
			vec_minus( dx1, Face[i].x, Cell[ic1].x, 3);
			vec_minus( dx2, Face[i].x, Cell[ic2].x, 3);
			r1 = vec_len( dx1,3 );
			r2 = vec_len( dx2,3 );
			Face[i].lambda = r2 / (r1+r2); // inverse-distance = (1/r1) / ( 1/r1 + 1/r2 )
		}
		else
			Face[i].lambda = 1.;
    }
	
	
    //-- cell
    int interfaceVoidCellCounter = Ncel;//out of Ncel range!
    int nVirtualCell = 0;
    for( i=0; i<Ncel; i++ )
    {
        for( j=0; j<3; j++ )
        {
            xv1[j] = Vert[ Cell[i].vertices[0] ] [j];
            xv2[j] = Vert[ Cell[i].vertices[1] ] [j];
            xv3[j] = Vert[ Cell[i].vertices[2] ] [j];
            xv4[j] = Vert[ Cell[i].vertices[3] ] [j];
            xv5[j] = Vert[ Cell[i].vertices[4] ] [j];
            xv6[j] = Vert[ Cell[i].vertices[5] ] [j];
            xv7[j] = Vert[ Cell[i].vertices[6] ] [j];
            xv8[j] = Vert[ Cell[i].vertices[7] ] [j];
        }

	// cell volume
        // way-1: |V|=Integ_V( div(x,0,0) ) = Integ_S( (x,0,0).n ).
        /* Cell[i].vol = 0.;
        for( j=0; j<Cell[i].nface; j++ )
        {
            int iface= Cell[i].face[j];
			if( Face[iface].cell1==i )
				Cell[i].vol += Face[iface].x[0] * Face[iface].n[0];
			else
				Cell[i].vol -= Face[iface].x[0] * Face[iface].n[0];
        } */
		// way-2:
		/* double ddp[3],dde[3],ddz[3];
		for( j=0;j<3;j++ ){
			ddp[j]=0.125*( -xv1[j]+xv2[j]+xv3[j]-xv4[j]-xv5[j]+xv6[j]+xv7[j]-xv8[j] );
			dde[j]=0.125*( -xv1[j]-xv2[j]+xv3[j]+xv4[j]-xv5[j]-xv6[j]+xv7[j]+xv8[j] );
			ddz[j]=0.125*( -xv1[j]-xv2[j]-xv3[j]-xv4[j]+xv5[j]+xv6[j]+xv7[j]+xv8[j] );
		}
		Cell[i].vol = ddp[0]*( dde[1]*ddz[2] - dde[2]*ddz[1] )
					 -ddp[1]*( dde[0]*ddz[2] - dde[2]*ddz[0] )
					 +ddp[2]*( dde[0]*ddz[1] - dde[1]*ddz[0] );
		Cell[i].vol *= 8.; */

		// cell-gravity
        // Hexa is divided into 6 Tets
        V1= TetVolum(xv1,xv2,xv4,xv6);
        V2= TetVolum(xv1,xv5,xv4,xv6);
        V3= TetVolum(xv4,xv5,xv8,xv6);
        V4= TetVolum(xv2,xv3,xv4,xv7);
        V5= TetVolum(xv2,xv6,xv4,xv7);
        V6= TetVolum(xv8,xv4,xv6,xv7);
	Cell[i].vol= V1 + V2 + V3 + V4 + V5 + V6;

        for( j=0; j<3; j++ ){
            xc1[j]= 0.25*( xv1[j]+xv2[j]+xv4[j]+xv6[j] );
            xc2[j]= 0.25*( xv1[j]+xv5[j]+xv4[j]+xv6[j] );
            xc3[j]= 0.25*( xv4[j]+xv5[j]+xv8[j]+xv6[j] );
            xc4[j]= 0.25*( xv2[j]+xv3[j]+xv4[j]+xv7[j] );
            xc5[j]= 0.25*( xv2[j]+xv6[j]+xv4[j]+xv7[j] );
            xc6[j]= 0.25*( xv8[j]+xv4[j]+xv6[j]+xv7[j] );
        }
	double Vsum = V1+V2+V3+V4+V5+V6;
        for( j=0; j<3; j++ )
            Cell[i].x[j]=1./Vsum*( V1*xc1[j] + V2*xc2[j] +
                                   V3*xc3[j] + V4*xc4[j] +
                                   V5*xc5[j] + V6*xc6[j] );
        if( fabs((Cell[i].vol - (V1+V2+V3+V4+V5+V6))/Cell[i].vol) > 1.e-3 ){
		char temp[256];
		sprintf(temp,"error in cell volume calculation%d vol:%f %f\n",i,Cell[i].vol,V1+V2+V3+V4+V5+V6 );
		throw logic_error(temp);
        }

        // neighbored cells
        for( j=0; j<Cell[i].nface; j++ )
        {
        	iface = Cell[i].face[j];
		nc    = Face[iface].cell1;
        	if( nc==i )
        		nc = Face[iface].cell2;
		if(nc==VOID_CELL_ON_INTERFACE){ //insert virtual cell for the interfaces
			set<int> verts;
			for(int i=0;i!=4;++i)
				verts.insert(Face[iface].vertices[i]);
			for(map<int,set<int> >::iterator it = interfaceNodes->begin();it!=interfaceNodes->end();++it){//for each interfaces
				vector<int> _intersectionVector;//hold result of intersection
				set<int>& _st = it->second;
				int partID    = it->first;
				set_difference(verts.begin(),verts.end(),_st.begin(),_st.end(),back_inserter(_intersectionVector));
				
				if(_intersectionVector.empty()){//verts belong to _st
					map<int,Interface>& _interfaces = dataPartition->interfaces;
					if(_interfaces.find(partID) == _interfaces.end() ){
						//create a new interface
						_interfaces.insert(make_pair<int,Interface>(partID,Interface(dataPartition->comRank,partID)));
					}
					_interfaces[partID].sendposis.push_back(i);	//i --> icell
					_interfaces[partID].recvposis.push_back(interfaceVoidCellCounter);
					assert(Face[iface].cell2==VOID_CELL_ON_INTERFACE);

					Cell[i].cell[j] = interfaceVoidCellCounter;	//add void cell
					Face[iface].cell2 = interfaceVoidCellCounter;

					//dataPartition->PRINT_LOG(partID);
					//dataPartition->PRINT_LOG(Face[iface].cell1);
					//dataPartition->PRINT_LOG(Face[iface].cell2);
					
					interfaceVoidCellCounter++;
					nVirtualCell++;
					break;	
				}
				
			}		

		}else{
        		Cell[i].cell[j]= nc; //nc = -10 when boundary
		}
        }
    }


	/* // output 
	ofstream of;
	of.open("face.dat");
	of<<Nfac<<endl;
	for( i=0; i<Nfac; i++ )
	{
		for( j=0; j<4; j++ )
			of<<Face[i].vertices[j]+1<<" ";
		of<<endl;
		of<<Face[i].n[0]<<" "<<Face[i].n[1]<<" "<<Face[i].n[2]<<" ";
		of<<Face[i].lambda<<" "<<Face[i].cell1+1<<" "<<Face[i].cell2+1;
		of<<endl;
	}
	of.close( );
	of.open("cell.dat");
	of<<Ncel<<endl;
	for( i=0; i<Ncel; i++ ){
		of<<Cell[i].vol<<" "<<Cell[i].x[0]<<" "<<Cell[i].x[1]<<" "<<Cell[i].x[2];
		of<<endl;
	}
	of.close(); */

    
    //-- boundary
    	for( i=0; i<Nbnd; i++ )
    	{
        	jf  = Bnd [i ].face  ;
        	jc  = Face[jf].cell1 ;
        	vec_minus( dx, Face[jf].x, Cell[jc].x,3 );
       		Bnd[i].distance =  vec_dot(dx, Face[jf].n, 3) / Face[jf].area;
	        if( Bnd[i].distance<0. )
        	{
			char temp[256];
			sprintf(temp,"wall distance is negatifve: %d\n",i);
			throw logic_error(temp);
       		}
	}

	//-- Face[].rlencos
	for( i=0; i<Nfac; i++ )
	{
		ic1= Face[i].cell1;
		ic2= Face[i].cell2;
		if( ic2<0 )
		{
			vec_minus( dx,Face[i].x,Cell[ic1].x,3 );
			// Face[i].rlencos= Face[i].area/vec_len(dx,3);
			Face[i].rlencos = Face[i].area / (vec_dot( dx,Face[i].n,3 )/Face[i].area);
		}
		else
		{
			vec_minus( dx,Cell[ic2].x,Cell[ic1].x,3 );
			// Face[i].rlencos= Face[i].area/vec_len(dx,3);
			Face[i].rlencos = Face[i].area / (vec_dot( dx,Face[i].n,3 )/Face[i].area);
		}

		// left and right auxiliary points Xpac, Xnac
		if( ic2>=0 ){
		double dx[3], ss;
		vec_minus( dx, Face[i].x, Cell[ic1].x, 3);
		ss = vec_dot(dx,Face[i].n,3)/ Face[i].area;
		for( j=0; j<3; j++ )
			dx[j] = ss * Face[i].n[j]/Face[i].area;
		vec_minus( Face[i].Xpac, Face[i].x, dx, 3);

		vec_minus( dx, Face[i].x, Cell[ic2].x, 3);
		ss = vec_dot(dx,Face[i].n,3)/ Face[i].area;
		for( j=0; j<3; j++ )
			dx[j] = ss * Face[i].n[j]/Face[i].area;
		vec_minus( Face[i].Xnac, Face[i].x, dx, 3);
		}
	}

    return nVirtualCell;
}

int NavierStokesSolver::CheckAndAllocate(int nVirtualCell)
{
	int i,c1,c2;
	// check if all boundaries are marked
	for( i=0; i<Nfac; i++ )
	{
		c1= Face[i].cell1;
		c2= Face[i].cell2;
		if( c2==VOID_CELL_ON_BOUNDARY || c2>=0 ) continue;
		char temp[256];
		sprintf(temp,"error in face right hand side\nfaceID: %d, cell1:%d cell2 %d\n",i,c1,c2);
		throw logic_error(temp);
	}

	
	//below is moved into NavierStokerSolver::Init
	// allocate variables
   	Rn = new double[Ncel+nVirtualCell];
	Un = new double[Ncel+nVirtualCell];
	Vn = new double[Ncel+nVirtualCell];
	Wn = new double[Ncel+nVirtualCell];
	Pn = new double[Ncel+nVirtualCell];
	Tn = new double[Ncel+nVirtualCell];
	TE = new double[Ncel+nVirtualCell];
	ED = new double[Ncel+nVirtualCell];
	RSn= new_Array2D<double>(Nspecies,Ncel+nVirtualCell);

	if( !IfSteady ){
		if( TimeScheme>=1 ){
			Rnp = new double[Ncel+nVirtualCell];
			Unp = new double[Ncel+nVirtualCell];
			Vnp = new double[Ncel+nVirtualCell];
			Wnp = new double[Ncel+nVirtualCell];
			Tnp = new double[Ncel+nVirtualCell];
			TEp = new double[Ncel+nVirtualCell];
			EDp = new double[Ncel+nVirtualCell];
			RSnp= new_Array2D<double>(Nspecies,Ncel+nVirtualCell);
		}
		if( TimeScheme>=2 ){
			Rnp2 = new double[Ncel+nVirtualCell];
			Unp2 = new double[Ncel+nVirtualCell];
			Vnp2 = new double[Ncel+nVirtualCell];
			Wnp2 = new double[Ncel+nVirtualCell];
			Tnp2 = new double[Ncel+nVirtualCell];
			TEp2 = new double[Ncel+nVirtualCell];
			EDp2 = new double[Ncel+nVirtualCell];
			RSnp2= new_Array2D<double>(Nspecies,Ncel+nVirtualCell);
		}
	}


	VisLam = new double[Ncel+nVirtualCell];
	VisTur = new double[Ncel+nVirtualCell];
	dPdX   = new_Array2D<double>(Ncel+nVirtualCell,3);
	dUdX   = new_Array2D<double>(Ncel+nVirtualCell,3);
	dVdX   = new_Array2D<double>(Ncel+nVirtualCell,3);
	dWdX   = new_Array2D<double>(Ncel+nVirtualCell,3);
	Apr    = new double[Ncel+nVirtualCell];
	dPhidX = new_Array2D<double>(Ncel+nVirtualCell,3);

	BRo = new double[Nbnd];
	BU  = new double[Nbnd];
	BV  = new double[Nbnd];
	BW  = new double[Nbnd];
	BTem= new double[Nbnd];
	BPre= new double[Nbnd];
	BRS = new_Array2D<double>(Nspecies,Nbnd);
	BTE = new double[Nbnd];
	BED = new double[Nbnd];

	RUFace = new double[Nfac];
	
	// laspack working array
	//V_Constr(&bs,   "rightU",    Ncel, Normal, True);
	//V_Constr(&bu,   "rightU",    Ncel, Normal, True);
	//V_Constr(&bv,   "rightU",    Ncel, Normal, True);
	//V_Constr(&bw,   "rightU",    Ncel, Normal, True);
	//V_Constr(&bp,   "rightU",    Ncel, Normal, True);
	//V_Constr(&xsol, "rightU",    Ncel, Normal, True);
	dataPartition->initPetsc();

	cur_time = 0.;
		

	return 0;
}
