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
        data = geofile.read(struct.calcsize('4i'));
        rank, nProcess,nc,nv = struct.unpack('4i',data)
        print 'processing rank: %d/%d nc: %d, nv: %d' %(rank,nProcess,nc,nv)
        gridDim.append((nc,nv));
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
    rank = 0;
    while True:
        data = datafile.read(struct.calcsize('2i'));
        rank,nvar = struct.unpack('2i',data);
        print 'processing rank: %d nvar:%d' % (rank,nvar)
        nc = gridDim[rank][0];
        _fluidData = np.zeros((nvar,nc));
        for i in range(0,nvar):
            for j in range(0,nc):
                data = datafile.read(struct.calcsize('d'));
                _fluidData[i,j], = struct.unpack('d',data);
        #print _fluidData;
        fluidData.append(_fluidData)
        if rank == nProcess-1:
            break;
    datafile.close();
    return gridDim, grid,connectivity,fluidData

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
    gridDim, grid, connectivity,fluidData = readBuf(directory)
    of = open('tec.dat','w')
    of.write('variables= "x","y","z","p","u","v","w","ro","T"')
    if fluidData[0].shape[0] >6 :
        of.write('",te","ed"')
    of.write('",Mach"\n')
    for rank in range(0,len(gridDim)):
        nc, nv = gridDim[rank];
        nvar = fluidData[rank].shape[0];
        of.write("ZONE N=%d, E=%d, VARLOCATION=([1-3]=NODAL,[4-%d]=CELLCENTERED) DATAPACKING=BLOCK, ZONETYPE=FEBRICK\n"%(nv,nc,3+nvar+1) );
        writeArray2File(of,grid[rank][:,0]); # x
        writeArray2File(of,grid[rank][:,1]); # y      
        writeArray2File(of,grid[rank][:,2]); # z       
        for i in range(0,nvar):
            writeArray2File(of,fluidData[rank][i,:]);
        #calc param goes here...
        writeArray2File(of,calcMach(fluidData[rank][0,:],\
                                    fluidData[rank][1,:],\
                                    fluidData[rank][2,:])\
                       )

        for i in range(0,nc):
            for connidx in connectivity[rank][i]:
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


