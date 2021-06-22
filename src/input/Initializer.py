import numpy as np
import h5py

class Initializer:
    
    
    def __init__(self):


        self.inputFile = ".../plasmoids/pl012/tilt_laser3D_6.h5"
        self.f = h5py.File(self.inputFile, 'r')
        
        self.Bx = np.array(self.f['vars/8_0'])
        self.By = np.array(self.f['vars/8_1'])
        self.Bz = np.array(self.f['vars/8_2'])
        
        self.Ex = np.array(self.f['vars/0_0'])
        self.Ey = np.array(self.f['vars/0_1'])
        self.Ez = np.array(self.f['vars/0_2'])
        
        self.f.close()
        
        # BOX
        Lx = 140.0
        Ly = 80.0
        Lz = 70.0
        
        self.flowDirection = 2 # 0 - X, 1 - Y, 2 - Z
        
        res = 0.4
        
        self.boxSize    = [Lx, Ly, Lz] # in d0 - ion inertia length
        self.boxSizePxl = [Lx/res, Ly/res, Lz/res]

        core = 1
        self.mpiCores  = [core,core,core]
        
        # time
        self.ts = 0.00005
        self.maxtsnum = 501
        self.outputStride = 250
        
        # output
        self.outputDir = "test_out/"
        self.fileTemplate = "proto_"
        
        # Particles
        self.ppc = 20
        self.mass  = 1
        self.charge = 1
        
        self.vel = [0.0, 0.0, 10000.0]

        self.magAmpl = 4
        self.eleAmpl = 1
        
        self.detectorPosition = 80*Lx
        self.sourcePosition = 10*Lx
        
        self.detectorSize = [2*Lx, 1*Ly]
    
    
    def getFlowDirection(self):
        return self.flowDirection
    
    #   spatial: left lower corner is (0,0,0)
    #   spatial: total box length in normalized units
    def getXright(self):
        return self.boxSize[0]
    
    def getYright(self):
        return self.boxSize[1]
    
    def getZright(self):
        return self.boxSize[2]
    
    # if number of pixels is less than 2, there is no direvative
    #   total box length in pixels Lx
    def getXresolution(self):
        return self.boxSizePxl[0]
    
    #   total box length in pixels Ly
    def getYresolution(self):
        return self.boxSizePxl[1]
    
    #   total box length in pixels Lz
    def getZresolution(self):
        return self.boxSizePxl[2]
    
    def getXmpiDomainNum(self):
        return self.mpiCores[0]
    
    def getYmpiDomainNum(self):
        return self.mpiCores[1]
    
    def getZmpiDomainNum(self):
        return self.mpiCores[2]
    
    #   time
    def getTimestep(self):
        return self.ts
    
    def getMaxTimestepsNum(self):
        return self.maxtsnum

    
    #   output:
    #   output: folder must be created manually
    def getOutputDir(self):
        return self.outputDir
    
    def getOutputFilenameTemplate(self):
        return self.fileTemplate
    
    def getOutputTimestep(self):
        return self.outputStride
    
    
    #   physics: particles
    def getParticlesPerCellNumber(self):
        return self.ppc
    
    #           species 1
    def getMass(self):
        return self.mass
    
    def getCharge(self):
        return self.charge
    
    
    #        modulus of intial velocity for protons
    def getVelocityX(self, x, y, z):
        return self.vel[0]
    
    def getVelocityY(self, x, y, z):
        return self.vel[1]
    
    def getVelocityZ(self, x, y, z):
        return self.vel[2]
    
    
    def loadVarCompAtPoint(self, fld, i, j, k):
        
        [Lx, Ly, Lz] = fld.shape
        
        if i >= Lx:
            i = Lx-1
        
        if j >= Ly:
            j = Ly-1
        
        if k >= Lz:
            k = Lz-1
        
        data = fld[int(i),int(j),int(k)]
        return data

    #        b field
    def getBfieldX(self, i, j, k):
        return self.magAmpl*self.loadVarCompAtPoint(self.Bx, i, j, k)
    
    def getBfieldY(self, i, j, k):
        return self.magAmpl*self.loadVarCompAtPoint(self.By, i, j, k)
    
    
    def getBfieldZ(self, i, j, k):
        return self.magAmpl*self.loadVarCompAtPoint(self.Bz, i, j, k)
    
    #        e field
    def getEfieldX(self, i, j, k):
        return self.eleAmpl*self.loadVarCompAtPoint(self.Ex, i, j, k)
    
    def getEfieldY(self, i, j, k):
        return self.eleAmpl*self.loadVarCompAtPoint(self.Ey, i, j, k)
    
    def getEfieldZ(self, i, j, k):
        return self.eleAmpl*self.loadVarCompAtPoint(self.Ez, i, j, k)
    
    #       detector
    def getDetectorPosition(self):
        return self.detectorPosition

    def getSourcePosition(self):
        return self.sourcePosition

    def getDetectorSizeX(self):
        return self.detectorSize[0]

    def getDetectorSizeY(self):
        return self.detectorSize[1]







