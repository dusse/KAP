#include "ModelInitializer.hpp"



using namespace std;
using namespace chrono;


ModelInitializer::ModelInitializer(shared_ptr<Loader> load,
                                   shared_ptr<GridManager> grid,
                                   shared_ptr<Pusher> push):loader(move(load)), gridMng(move(grid)), pusher(move(push)){
    logger.reset(new Logger());
    
    initialize();
    logger->writeMsg("[ModelInitializer] create...OK", DEBUG);
}

void ModelInitializer::initialize()
{
    logger->writeMsg("[ModelInitializer] initialize...OK", DEBUG);
    
    initElectroMagneticField();
    initParticles();
    
}


void ModelInitializer::initElectroMagneticField(){
    logger->writeMsg("[ModelInitializer] init Electro Magnetic field...", DEBUG);
    
    int xSize = loader->resolution[0];
    int ySize = loader->resolution[1];
    int zSize = loader->resolution[2];
    
    vector<double> bField = {0,0,0};
    vector<double> eField = {0,0,0};
    
    int i,j,k, idx, idxG1, idxG2;
    
     for ( i=0; i<xSize+1; i++ ){
         for ( j=0; j<ySize+1; j++ ){
             for ( k=0; k<zSize+1; k++ ){
                 
                 idxG1 = IDX(i  ,j  ,k  ,xSize+1,ySize+1,zSize+1);
                 idxG2 = IDX(i+1,j+1,k+1,xSize+2,ySize+2,zSize+2);
                 
                 bField = loader->getBfield(i,j,k);
                 gridMng->setVectorVariableForNodeG1(idxG1, VectorVar(MAGNETIC, bField));
                 
                 eField = loader->getEfield(i,j,k);
                 gridMng->setVectorVariableForNodeG2(idxG2, VectorVar(ELECTRIC, eField));
                 
             }
         }
     }

    gridMng->smoothG1(MAGNETIC);
    gridMng->smoothG1(MAGNETIC);
    gridMng->smooth(ELECTRIC);
    
    logger->writeMsg("[ModelInitializer] init Electro Magnetic field...OK", DEBUG);
    
}
                                         



void ModelInitializer::initParticles(){
    logger->writeMsg("[ModelInitializer] init particles...", DEBUG);
    double domainShiftX = loader->boxCoordinates[0][0];
    double domainShiftY = loader->boxCoordinates[1][0];
    double domainShiftZ = loader->boxCoordinates[2][0];
    
    int xRes = loader->resolution[0];
    int yRes = loader->resolution[1];
    int zRes = loader->resolution[2];
    
    double cellSizeX = loader->spatialSteps[0];
    double cellSizeY = loader->spatialSteps[1];
    double cellSizeZ = loader->spatialSteps[2];
    
    double detSizeA = loader->detectorSize[0];
    double detSizeB = loader->detectorSize[1];
    
    double detCellSizeA, detCellSizeB;
    
    switch (loader->flowDirection) {
        case 0:
            detCellSizeA = detSizeA/yRes;
            detCellSizeB = detSizeB/zRes;
            break;
        case 1:
            detCellSizeA = detSizeA/xRes;
            detCellSizeB = detSizeB/zRes;
            break;
        case 2:
            detCellSizeA = detSizeA/xRes;
            detCellSizeB = detSizeB/yRes;
            break;
            
        default:
            throw runtime_error("no such direction!");
    }
    
    
    double halfDetSizeA = 0.5*detSizeA;
    double halfDetSizeB = 0.5*detSizeB;
    
    int numOfSpecies = 1;
    
    int particle_idx=0, idxOnG2, idx, ptclIDX, spn;
    double pos[3], vpb[3];
    vector<double> vel;
    double r1, r2;
    
    int ppc = loader->getParticlesPerCellNumber();
    int totalCellsInPlaneXY = xRes*yRes;
    totParticleNumberInDomain = totalCellsInPlaneXY*ppc;
    
    pusher->initParticles(totParticleNumberInDomain);
    
    double totalPrtclNumber = pusher->getTotalParticleNumber();
    
    pusher->setParticleMass(loader->getMass());
    pusher->setParticleCharge(loader->getCharge());
    pusher->setParticleWeight(1);

    double requiredPrtclNum;
    
    int i, j, k;
    
    switch (loader->flowDirection) {
            
        case 0:
            for( j = 0; j < yRes; j++){
                for( k = 0; k < zRes; k++) {
                    
                    for (ptclIDX=0; ptclIDX < ppc; ptclIDX++){
                        r1 = RNM;
                        r2 = RNM;
                        pos[0] = cellSizeX + EPS8;
                        pos[1] = (j + r1) * cellSizeY + domainShiftY;
                        pos[2] = (k + r2) * cellSizeZ + domainShiftZ;
                        
                        double pos2Save[6] = {pos[0], pos[1], pos[2],
                            pos[0], pos[1], pos[2]};
                        pusher->setParticlePosition(particle_idx, pos2Save);
                        
                        pusher->setParticleType(particle_idx,  NOTFINISHED);
                        
                        vel = loader->getVelocity(pos[0], pos[1], pos[2]);
                        
                        double detY = (j + r1) * detCellSizeA - halfDetSizeA;
                        double detZ = (k + r2) * detCellSizeB - halfDetSizeB;
                        
                        double dir[3] = {loader->sourcePosition
                            + loader->detectorPosition, detY, detZ};
                        
                        double mod = sqrt(dir[0]*dir[0]
                                          +dir[1]*dir[1]
                                          +dir[2]*dir[2]);
                        
                        double velMod = vel[0];
                        
                        vpb[0] = velMod*dir[0]/mod;
                        vpb[1] = velMod*dir[1]/mod;
                        vpb[2] = velMod*dir[2]/mod;
                        
                        double vel2Save[6] = {vpb[0], vpb[1], vpb[2],
                            vpb[0], vpb[1], vpb[2]};
                        
                        pusher->setParticleVelocity(particle_idx, vel2Save);
                        
                        particle_idx++;
                        
                    }
                }
            }

            break;
            
        case 1:
            for( i = 0; i < xRes; i++){
                for( k = 0; k < zRes; k++) {
                    
                    for (ptclIDX=0; ptclIDX < ppc; ptclIDX++){
                        r1 = RNM;
                        r2 = RNM;
                        pos[0] = (i + r1) * cellSizeX + domainShiftX;
                        pos[1] = cellSizeY + EPS8;
                        pos[2] = (k + r2) * cellSizeZ + domainShiftZ;
                        
                        double pos2Save[6] = {pos[0], pos[1], pos[2],
                            pos[0], pos[1], pos[2]};
                        pusher->setParticlePosition(particle_idx, pos2Save);
                        
                        pusher->setParticleType(particle_idx,  NOTFINISHED);
                        
                        vel = loader->getVelocity(pos[0], pos[1], pos[2]);
                        
                        double detX = (i + r1) * detCellSizeA - halfDetSizeA;
                        double detZ = (k + r2) * detCellSizeB - halfDetSizeB;
                        
                        double dir[3] = {detX, loader->sourcePosition
                            + loader->detectorPosition, detZ};
                        
                        double mod = sqrt(dir[0]*dir[0]
                                          +dir[1]*dir[1]
                                          +dir[2]*dir[2]);
                        
                        double velMod = vel[1];
                        
                        vpb[0] = velMod*dir[0]/mod;
                        vpb[1] = velMod*dir[1]/mod;
                        vpb[2] = velMod*dir[2]/mod;
                        
                        double vel2Save[6] = {vpb[0], vpb[1], vpb[2],
                            vpb[0], vpb[1], vpb[2]};
                        
                        pusher->setParticleVelocity(particle_idx, vel2Save);
                        
                        particle_idx++;
                        
                    }
                }
            }
            
            break;

        case 2:
            for( i = 0; i < xRes; i++){
                for( j = 0; j < yRes; j++) {
                    
                    for (ptclIDX=0; ptclIDX < ppc; ptclIDX++){
                        r1 = RNM;
                        r2 = RNM;
                        pos[0] = (i + r1) * cellSizeX + domainShiftX;
                        pos[1] = (j + r2) * cellSizeY + domainShiftY;
                        pos[2] = cellSizeZ + EPS8;
                        
                        double pos2Save[6] = {pos[0], pos[1], pos[2],
                            pos[0], pos[1], pos[2]};
                        pusher->setParticlePosition(particle_idx, pos2Save);
                        
                        pusher->setParticleType(particle_idx,  NOTFINISHED);
                        
                        vel = loader->getVelocity(pos[0], pos[1], pos[2]);
                        
                        // origin by default is left lower coner
                        double detX = (i + r1) * detCellSizeA - halfDetSizeA;
                        double detY = (j + r2) * detCellSizeB - halfDetSizeB;
                        
                        
                        double dir[3] = {detX, detY,
                            loader->sourcePosition
                            + loader->detectorPosition};
                        
                        double mod = sqrt(dir[0]*dir[0]
                                          +dir[1]*dir[1]
                                          +dir[2]*dir[2]);
                        
                        double velMod = vel[2];
                        
                        vpb[0] = velMod*dir[0]/mod;
                        vpb[1] = velMod*dir[1]/mod;
                        vpb[2] = velMod*dir[2]/mod;
                        
                        double vel2Save[6] = {vpb[0], vpb[1], vpb[2],
                            vpb[0], vpb[1], vpb[2]};
                        
                        pusher->setParticleVelocity(particle_idx, vel2Save);
                        
                        particle_idx++;
                        
                    }
                }
            }
            break;
            
        default:
            throw runtime_error("no such direction!");
    }
    
    
    pusher->initTrackedParticles();
    
    logger->writeMsg("[ModelInitializer] init particles...OK", DEBUG);
}


ModelInitializer::~ModelInitializer(){
    finilize();
    logger->writeMsg("[ModelInitializer] delete...OK", DEBUG);
}

void ModelInitializer::finilize(){
    
}
