import sys
import struct
import numpy as np

def compressIn(filename):
    print 'compressing...'
    f = open(filename,"r")
    foundNode = False
    foundElement = False;
    of = open('_'+filename,'wb')
    line = f.readline()
    while len(line )>0:
        #raw_input()
        if line.startswith("$Nodes"):
            foundNode = True
        elif line.startswith("$Elements"):
            foundElement = True;
        elif line.startswith("$EndNodes") or line.startswith("$EndElements"):
            pass
        else:
            numbers = line.split()
            if not foundNode:
                pass
            elif len(numbers)==1:
                of.write(struct.pack('q',int(numbers[0])))
            elif foundElement:
                for number in numbers:
                    of.write(struct.pack('i',int(number)))
            elif foundNode:
                of.write(struct.pack('i',int(numbers[0])))
                of.write(struct.pack('d',float(numbers[1])))
                of.write(struct.pack('d',float(numbers[2])))
                of.write(struct.pack('d',float(numbers[3])))
        line = f.readline()
    print 'done'

def readBuf(directory):
    geofile = open(directory+'geo','rb')
    rank,nProcess = 0, 99999;
    grid = [];
    connectivity = [];
    gridDim = [];
    while True:
        data = geofile.read(struct.calcsize('5i'));
        rank, nProcess,nc,nf,nv = struct.unpack('5i',data)
        print 'processing rank: %d/%d nc: %d,nf: %d, nv: %d' %(rank,nProcess,nc,nf,nv)
        gridDim.append((nc,nf,nv));
        _grid = np.zeros((nv,3));
        for i in range(0,nv):
            data = geofile.read(struct.calcsize('3d'));
            point = struct.unpack('3d',data);
            _grid[i,:]  = point;
            #print point
            #raw_input();
        #print _grid;

        _connectivity = np.zeros((nc,8),int);
        for i in range(0,nc):
            data = geofile.read(struct.calcsize('8i'));
            _connectivity[i,:] = struct.unpack('8i',data)

        print _connectivity;
        grid.append(_grid);
        connectivity.append(_connectivity);
        if rank == nProcess - 1:
            break;
    geofile.close();
    #process data file
    datafile = open(directory+'data','rb')
    fluidData = [] 
    solidData = [] 
    rank = 0;
    while True:
        data = datafile.read(struct.calcsize('3i'));
        rank,nvarF,nvarS = struct.unpack('3i',data);
        print 'processing rank: %d nvarfluid:%d, nvarrsolid %d' % (rank,nvarF,nvarS)
        nc = gridDim[rank][0];
        nf = gridDim[rank][1];
        _fluidData = np.zeros((nvarF,nf));
        for i in range(0,nvarF):
            for j in range(0,nf):
                data = datafile.read(struct.calcsize('d'));
                _fluidData[i,j], = struct.unpack('d',data);
        fluidData.append(_fluidData)
        _solidData = np.zeros((nvarS,nc-nf));
        for i in range(0,nvarS):
            for j in range(nf,nc):
                data = datafile.read(struct.calcsize('d'));
                _solidData[i,j-nf], = struct.unpack('d',data);
        solidData.append(_solidData)
        #print _fluidData;
        if rank == nProcess-1:
            break;
    datafile.close();
    return gridDim, grid,connectivity,fluidData,solidData

def writeArray2File(f,array):
    for i,number in enumerate(array):
        f.write(str(number)+' ')
        if i % 5 == 4:
            f.write('\n')
    f.write('\n')

def calcMach(u,v,w):
    for i,j,k in zip(u,v,w):
        yield i**2 + j**2 + k**2

def buildTec(directory):
    gridDim, grid, connectivity,fluidData,solidData = readBuf(directory)
    of = open('tec.dat','w')
    of.write('variables= "x","y","z","p","u","v","w","ro","T"')
    if fluidData[0].shape[0] >6 :
        of.write('",te","ed"')
    of.write('",Mach"\n')
    #print fluid part...
    for rank in range(0,len(gridDim)):
        nc,nf, nv = gridDim[rank];
        nvar = fluidData[rank].shape[0];

        nodePool = {};
        icounter = 0;
        for i in range(0,nf):
            for connidx in connectivity[rank][i]:
                nodePool[connidx] = icounter++
        subgridx = filter(lambda (i,j): i in nodePool, enumerate(grid[rank][:,0]))
        subgridy = filter(lambda (i,j): i in nodePool, enumerate(grid[rank][:,1]))
        subgridz = filter(lambda (i,j): i in nodePool, enumerate(grid[rank][:,2]))

        of.write("ZONE N=%d, E=%d, VARLOCATION=([1-3]=NODAL,[4-%d]=CELLCENTERED) DATAPACKING=BLOCK, ZONETYPE=FEBRICK\n"%(len(subgridx),nf,3+nvar+1) );
        writeArray2File(of,subgridx); # x
        writeArray2File(of,subgridy); # y      
        writeArray2File(of,subgridz); # z       
        for i in range(0,nvar):
            writeArray2File(of,fluidData[rank][i,:]);
        #calc param goes here...
        writeArray2File(of,calcMach(fluidData[rank][0,:],\
                                    fluidData[rank][1,:],\
                                    fluidData[rank][2,:])\
                       )

        for i in range(0,nf):
            for connidx in map(lambda i : nodePool[i], connectivity[rank][i]): # map to local index
                of.write(str(connidx)+' ');
            of.write('\n')
        of.write('\n')

    #print solid part...
    for rank in range(0,len(gridDim)):
        nc,nf, nv = gridDim[rank];
        nvar = solidData[rank].shape[0];
        nodePool = {};
        icounter = 0;
        for i in range(nf,nc):
            for connidx in connectivity[rank][i]:
                nodePool[connidx] = icounter++
        subgridx = filter(lambda (i,j): i in nodePool, enumerate(grid[rank][:,0]))
        subgridy = filter(lambda (i,j): i in nodePool, enumerate(grid[rank][:,1]))
        subgridz = filter(lambda (i,j): i in nodePool, enumerate(grid[rank][:,2]))

        of.write("ZONE N=%d, E=%d, VARLOCATION=([1-3]=NODAL,[4-%d]=CELLCENTERED) DATAPACKING=BLOCK, ZONETYPE=FEBRICK\n"%(len(subgridx),nc-nf,3+1) );
        writeArray2File(of,subgridx); # x
        writeArray2File(of,subgridy); # y      
        writeArray2File(of,subgridz); # z       
        #for i in range(0,nvar):
        writeArray2File(of,solidData[rank][0,:]);
                       )
        for i in range(nf,nc):
            for connidx in map(lambda i : nodePool[i], connectivity[rank][i]): # map to local index
                of.write(str(connidx)+' ');
            of.write('\n')
        of.write('\n')
 
    
    of.close()

if __name__ == "__main__":
    if len(sys.argv)<= 1 or sys.argv[1] == "-help":
        print '-in: filename'
        print '-out: dir'
    elif  sys.argv[1] == "-in":
        compressIn(sys.argv[2])
    elif  sys.argv[1] == '-out':
        buildTec(sys.argv[2]);


