#include "BoundaryManager.hpp"

using namespace std;
using namespace chrono;

BoundaryManager::BoundaryManager(shared_ptr<Loader> ldr):loader(move(ldr)){
    logger.reset(new Logger());
    initialize();
    string msg ="[BoundaryManager] init...OK";
    logger->writeMsg(msg.c_str(), DEBUG);
}



void BoundaryManager::initialize(){
    logger->writeMsg("[BoundaryManager] initialize() ...", DEBUG);
    leavingParticles.reserve(NUM_OF_LEAVING_PACTICLES);
    for(int t = 0; t < 27; t++) {
        domain2send[t] = 0;
    }
    
    logger->writeMsg("[BoundaryManager] initialize() ...OK", DEBUG);
}



int BoundaryManager::isPtclOutOfDomain(double pos[3]){
    
    double domainShiftX = loader->boxCoordinates[0][0];
    double domainShiftY = loader->boxCoordinates[1][0];
    double domainShiftZ = loader->boxCoordinates[2][0];
    double domainXSize = loader->resolution[0];
    double domainYSize = loader->resolution[1];
    double domainZSize = loader->resolution[2];
    double dx = loader->spatialSteps[0];
    double dy = loader->spatialSteps[1];
    double dz = loader->spatialSteps[2];
    
    double x, y, z;
    int a, b, c, t;
    
    x = (pos[0] - domainShiftX)/domainXSize/dx;
    y = (pos[1] - domainShiftY)/domainYSize/dy;
    z = (pos[2] - domainShiftZ)/domainZSize/dz;
    
    a = (int)floor(x);
    b = (int)floor(y);
    c = (int)floor(z);
    
    t = (1+c)+3*((1+b)+3*(1+a));
    
    return (t == 13) ? IN : t;
}


int BoundaryManager::getNumberOfFinished(){
    return finished;
}

void BoundaryManager::applyBC(Particle** particles,  vector<shared_ptr<Particle>> &particles2add){
    logger->writeMsg(("[BoundaryManager] applyBC .. leavingParticles = "
                     +to_string(leavingParticles.size())).c_str(), DEBUG);
    
    double domainShiftX = loader->boxCoordinates[0][0];
    double domainShiftY = loader->boxCoordinates[1][0];
    double domainShiftZ = loader->boxCoordinates[2][0];
    double domainXSize = loader->resolution[0];
    double domainYSize = loader->resolution[1];
    double domainZSize = loader->resolution[2];
    double dx = loader->spatialSteps[0];
    double dy = loader->spatialSteps[1];
    double dz = loader->spatialSteps[2];
    
    double Lx = loader->boxSizes[0], Ly = loader->boxSizes[1], Lz = loader->boxSizes[2];
    
    
    int idx, t, a, b, c;
    
    double xloc, yloc, zloc;
    double xnew, ynew, znew;
    double xnewCor, ynewCor, znewCor;
    double x0, y0, z0;
    double x0Cor, y0Cor, z0Cor;
    double x0Pre, y0Pre, z0Pre;
    
    double *sendBuf[27];
    double *recvBuf[27];
    int partcls2send[27];
    int partcls2recv[27];
    
    for(t = 0; t < 27; t++) {
        sendBuf[t] = new double[domain2send[t]*PARTICLES_SIZE*sizeof(double)];
        recvBuf[t] = new double[EXPECTED_NUM_OF_PARTICLES*PARTICLES_SIZE*sizeof(double)];
        partcls2send[t] = 0;
        partcls2recv[t] = EXPECTED_NUM_OF_PARTICLES;
    }

    double* prtclPos;
    
    vector<int> removeFromLeaving;
    removeFromLeaving.reserve(leavingParticles.size());
    
    
    for(int ptclNum = 0; ptclNum < leavingParticles.size(); ptclNum++){
        idx = leavingParticles[ptclNum];
        
        prtclPos = particles[idx]->getPosition();
        
        x0 = prtclPos[0];
        y0 = prtclPos[1];
        z0 = prtclPos[2];
        
        //check whom to send
        xloc = (x0 - domainShiftX)/domainXSize/dx;
        yloc = (y0 - domainShiftY)/domainYSize/dy;
        zloc = (z0 - domainShiftZ)/domainZSize/dz;
        
        a = (int) floor(xloc);
        b = (int) floor(yloc);
        c = (int) floor(zloc);
        
        t = (1+c)+3*((1+b)+3*(1+a));
        
        if(applyOutflowBC(particles[idx]) == 1){
            //no need to send particle but need to remove
            continue;
        }
        
        
        if(t != 13){
            particles[idx]->serialize(sendBuf[t], PARTICLES_SIZE*partcls2send[t]);
            partcls2send[t] += 1;
        }else{
            removeFromLeaving.push_back(ptclNum);
        }
    }
    
    string msd ="[BoundaryManager] total leaving ="+to_string(leavingParticles.size());
    logger->writeMsg(msd.c_str(), DEBUG);
    
    for(int ptclNum = removeFromLeaving.size()-1; ptclNum >= 0 ; ptclNum--){
        leavingParticles.erase(leavingParticles.begin()+ptclNum);
    }

    
    MPI_Status st;
    int receivedTot;
    for(t=0; t < 27; t++){
            if(t != 13){
                
                MPI_Sendrecv(sendBuf[t], PARTICLES_SIZE*sizeof(double)*partcls2send[t],
                             MPI_DOUBLE, loader->neighbors2Send[t], t,
                             recvBuf[t], PARTICLES_SIZE*sizeof(double)*partcls2recv[t],
                             MPI_DOUBLE, loader->neighbors2Recv[t], t,
                             MPI_COMM_WORLD, &st);
                
                MPI_Get_count(&st, MPI_DOUBLE, &receivedTot);

                receivedTot = receivedTot/PARTICLES_SIZE/sizeof(double);
    
                for (int ptclNum = 0; ptclNum < receivedTot; ptclNum++){
                    particles2add.push_back(shared_ptr<Particle>(new Particle));
                    int idxOfAdded = particles2add.size()-1;
                    particles2add[idxOfAdded]->deserialize(recvBuf[t], PARTICLES_SIZE*ptclNum);
                }
            }
    }
    
    for(t = 0; t < 27; t++) {
                delete [] sendBuf[t];
                delete [] recvBuf[t];
    }
}


int BoundaryManager::applyOutflowBC(Particle* particle){
    int remove = 0;
    double* prtclPos = particle->getPosition();
    double L, r0;
    
    for(int i=0;i<3;i++){
        L  = loader->boxSizes[i];
        r0 = prtclPos[i];
        
        if (r0 < 0.0){
            particle->setType(FINISHED);
            finished++;
            remove = 1;
            break;
        }
        if (r0 > L){
            particle->setType(FINISHED);
            finished++;
            remove = 1;
            break;
        }
    }
    return remove;
}


vector<double> BoundaryManager::getLeavingParticlesIdxs(){
    return leavingParticles;
}

void BoundaryManager::reset(){
    logger->writeMsg("[BoundaryManager] reset ..", DEBUG);
    leavingParticles.clear();
    for(int t = 0; t < 27; t++) {
        domain2send[t] = 0;
    }
}

void BoundaryManager::storeParticle(int idx, double pos[3]){
    double domainShiftX = loader->boxCoordinates[0][0];
    double domainShiftY = loader->boxCoordinates[1][0];
    double domainShiftZ = loader->boxCoordinates[2][0];
    double domainXSize = loader->resolution[0];
    double domainYSize = loader->resolution[1];
    double domainZSize = loader->resolution[2];
    double dx = loader->spatialSteps[0];
    double dy = loader->spatialSteps[1];
    double dz = loader->spatialSteps[2];
    
    double x, y, z;
    int a, b, c, domain;
    
    x = (pos[0] - domainShiftX)/domainXSize/dx;
    y = (pos[1] - domainShiftY)/domainYSize/dy;
    z = (pos[2] - domainShiftZ)/domainZSize/dz;
    
    a = (int)floor(x);
    b = (int)floor(y);
    c = (int)floor(z);
    
    domain = (1+c)+3*((1+b)+3*(1+a));

    leavingParticles.push_back(idx);
    domain2send[domain] += 1;// need to know for memory preallocation
}












