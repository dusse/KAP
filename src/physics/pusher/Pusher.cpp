#include "Pusher.hpp"

using namespace std;
using namespace chrono;


Pusher::Pusher(shared_ptr<Loader> ldr,
               shared_ptr<GridManager> gridMnr,
               shared_ptr<BoundaryManager> boundMnr):loader(move(ldr)),
                gridMgr(move(gridMnr)), boundaryMgr(move(boundMnr)){
    
    logger.reset(new Logger());
    initialize();
    logger->writeMsg("[Pusher] create...OK", DEBUG);
}

Pusher::~Pusher(){
    delete[] particles;
    delete[] trackedParticles;
}
               
void Pusher::initialize(){
    
}



void Pusher::setParticlePosition(int idx, double input[6] ){
    particles[idx]->setPosition(input);
}

void Pusher::setParticleVelocity(int idx, double input[6]){
    particles[idx]->setVelocity(input);
}

void Pusher::setParticleType(int idx, int input){
    particles[idx]->setType(input);
}


void Pusher::initParticles(int num){
    
    double pos[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double vel[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    
    currentPartclNumOnDomain = num;
    int const ALLOCATION_FACTOR = 2;
    totalNum  = ALLOCATION_FACTOR*num;
    particles = new Particle*[totalNum];
    for(int i=0; i<totalNum; i++){
        Particle* pa = new Particle();
        particles[i] = pa;
        particles[i]->setPosition(pos);
        particles[i]->setVelocity(vel);
        particles[i]->setType(NOTFINISHED);
    }
    int TOT_IN_BOX = 0;
    MPI_Allreduce(&currentPartclNumOnDomain, &TOT_IN_BOX, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    logger->writeMsg(("[Pusher] total particles initialized in the box = "
                                                +to_string(TOT_IN_BOX)).c_str(),  INFO);
    totinBoxInit = TOT_IN_BOX;
    
    
}

void Pusher::initTrackedParticles(){
    
    trackedParticles = new double[NUMBER_TRACKED_PARTICLES*3];
    
    for( int idx=0; idx < NUMBER_TRACKED_PARTICLES; idx++){
        trackedParticles[3*idx + 0] = 0.0;
        trackedParticles[3*idx + 1] = 0.0;
        trackedParticles[3*idx + 2] = 0.0;
    }
    
    int countAdded = 0;
    double* prtclPos;
    double x,y, r;
    
    double Lx = loader->boxSizes[0],  Ly = loader->boxSizes[1];
    
    for( int idx=0; idx < currentPartclNumOnDomain; idx++){
        
        prtclPos = particles[idx]->getPosition();
        
        x = prtclPos[0]-0.5*Lx;
        y = prtclPos[1]-0.5*Ly;
        
        r = sqrt(x*x + y*y);
        
       
        
        if( r < 2 && countAdded < NUMBER_TRACKED_PARTICLES){
            trackList.insert(idx);
            countAdded++;
        }
        
    }

    
    
    logger->writeMsg(("[Pusher] tracked Particles number = "+to_string(countAdded)).c_str(),  DEBUG);
}



void Pusher::setParticleWeight(double weight){
    this->weight = weight;
    logger->writeMsg(("[Pusher] weight = "+to_string(weight)).c_str(),  DEBUG);
}

void Pusher::setParticleCharge(double charge){
    this->charge = charge;
    logger->writeMsg(("[Pusher] charge = "+to_string(charge)).c_str(),  DEBUG);
}

void Pusher::setParticleMass(double mass){
    this->mass = mass;
    logger->writeMsg(("[Pusher] mass = "+to_string(mass)).c_str(),  DEBUG);
}

int Pusher::getTotalParticleNumber(){
    return currentPartclNumOnDomain;
}



/*                        E(i,j+1)           |              E(i+1,j+1)
 *                           x-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-x
 *                           .               |               .
 *                           |               |   CELL i,j    |
 *                           .               |               .
 *                         dy/2              |   O(x,y)      |
 *                           .               |     ->        .
 *                           |               |i=int(x/dx+0.5)|
 *                           .               |j=int(y/dy+0.5).
 *                           |               |               |
 *            .B(i-1,j)______________________.B(i,j)____________________.B(i+1,j)
 *                           |               |               |
 *                           .               |               .
 *                           |               |               |
 *                         dy/2              |               .
 *                           |               |               |
 *                           .               |               .
 *                           |               |               |
 *                           . E(i,j)        |               .
 *                           x-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-x E(i+1,j)
 *                                  dx/2     |      dx/2
 */


Particle** Pusher::getParticles(){
    return particles;
}


void Pusher::applyFreeSpace(){
    
    double dx = loader->spatialSteps[0];
    double dy = loader->spatialSteps[1];
    double dz = loader->spatialSteps[2];
    
    int xSize = loader->resolution[0];
    int ySize = loader->resolution[1];
    int zSize = loader->resolution[2];
    
    double domainShiftX = loader->boxCoordinates[0][0];
    double domainShiftY = loader->boxCoordinates[1][0];
    double domainShiftZ = loader->boxCoordinates[2][0];

    double Lx = loader->boxSizes[0],
           Ly = loader->boxSizes[1],
           Lz = loader->boxSizes[2];
    
    double x, y, z;
    int idxG1, idxG2;
    
    double* prtclPos;
    double* prtclVel;
    
    double freeTime, freeSpaceX, freeSpaceY, freeSpaceZ;
    
    int finished = boundaryMgr->getNumberOfFinished();
    
    for( int idx=0; idx < currentPartclNumOnDomain; idx++){
        
        prtclPos = particles[idx]->getPosition();
        prtclVel = particles[idx]->getVelocity();
        
        switch (loader->flowDirection) {
            case 0:
                
                y = prtclPos[1];
                z = prtclPos[2];
                
                freeTime = (loader->detectorPosition - Lx)/ prtclVel[0];
                
                freeSpaceY = prtclVel[1]*freeTime;
                freeSpaceZ = prtclVel[2]*freeTime;
                
                y += freeSpaceY;
                z += freeSpaceZ;
                
                particles[idx]->setPosition(4, y);
                particles[idx]->setPosition(5, z);
                
                break;
                
            case 1:
                
                x = prtclPos[0];
                z = prtclPos[2];
                
                freeTime = (loader->detectorPosition - Ly)/ prtclVel[1];
                
                freeSpaceX = prtclVel[0]*freeTime;
                freeSpaceZ = prtclVel[2]*freeTime;
                
                x += freeSpaceX;
                z += freeSpaceZ;
                
                particles[idx]->setPosition(3, x);
                particles[idx]->setPosition(5, z);
                
                break;
                
            case 2:
                
                x = prtclPos[0];
                y = prtclPos[1];
                
                freeTime = (loader->detectorPosition - Lz)/ prtclVel[2];
                
                freeSpaceX = prtclVel[0]*freeTime;
                freeSpaceY = prtclVel[1]*freeTime;
                
                x += freeSpaceX;
                y += freeSpaceY;
                
                particles[idx]->setPosition(3, x);
                particles[idx]->setPosition(4, y);
                
                break;
                
            default:
                throw runtime_error("no such direction!");
        }
        
        

        
//        if(trackList.count(idx) && particles[idx]->getType() == FINISHED){
//            set<int>::iterator it = trackList.find(idx);
//            int trackedIdx = distance(trackList.begin(), it);;
//            trackedParticles[3*trackedIdx + 0] = x;
//            trackedParticles[3*trackedIdx + 1] = y;
//            trackedParticles[3*trackedIdx + 2] = loader->detectorPosition;
//        }
        
    }
    
    double detSizeA = loader->detectorSize[0];
    double detSizeB = loader->detectorSize[1];
    
    vector<int> detRes = gridMgr->getDetectorResolution();
    
    double detCellSizeA = detSizeA/detRes[0];
    double detCellSizeB = detSizeB/detRes[1];
    
    double halfDetSizeA = 0.5*detSizeA;
    double halfDetSizeB = 0.5*detSizeB;
    
    double halfLxBox = 0.5*Lx;
    double halfLyBox = 0.5*Ly;
    double halfLzBox = 0.5*Lz;
    
    double shift2DetectorA;
    double shift2DetectorB;
    
    switch (loader->flowDirection) {
        case 0:
            shift2DetectorA = halfDetSizeA - halfLyBox;
            shift2DetectorB = halfDetSizeB - halfLzBox;
            break;
        case 1:
            shift2DetectorA = halfDetSizeA - halfLxBox;
            shift2DetectorB = halfDetSizeB - halfLzBox;
            break;
        case 2:
            shift2DetectorA = halfDetSizeA - halfLxBox;
            shift2DetectorB = halfDetSizeB - halfLyBox;
            break;
            
        default:
            throw runtime_error("no such direction!");
    }

    int pastParticles = 0;
    int i, j, k;
    for( int idx=0; idx < currentPartclNumOnDomain; idx++){
        
        prtclPos = particles[idx]->getPosition();
        
        
        switch (loader->flowDirection) {
                
            case 0:
                
                
                y  = prtclPos[4] + shift2DetectorA;
                z  = prtclPos[5] + shift2DetectorB;
                
                y = y > detSizeA ? detSizeA : y;
                y = y < 0 ? 0 : y;
                
                z = z > detSizeB ? detSizeB : z;
                z = z < 0 ? 0 : z;
                
                if(y == 0 || y == detSizeA || z == 0 || z == detSizeB){
                    pastParticles++;
                    continue;
                }
                
                
                j  = int(y / detCellSizeA );
                k  = int(z / detCellSizeB );
                
                idxG2 = IDX(j, k, 0, detRes[0], detRes[1], 1);
                
                break;
                
            case 1:
                x  = prtclPos[3] + shift2DetectorA;
                z  = prtclPos[5] + shift2DetectorB;
                
                x = x > detSizeA ? detSizeA : x;
                x = x < 0 ? 0 : x;
                
                z = z > detSizeB ? detSizeB : z;
                z = z < 0 ? 0 : z;
                
                if(x == 0 || x == detSizeA || z == 0 || z == detSizeB){
                    pastParticles++;
                    continue;
                }
                
                
                i  = int(x / detCellSizeA );
                k  = int(z / detCellSizeB );
                
                idxG2 = IDX(i, k, 0, detRes[0], detRes[1], 1);
                
                break;

                
            case 2:
                
                x  = prtclPos[3] + shift2DetectorA;
                y  = prtclPos[4] + shift2DetectorB;
                
                x = x > detSizeA ? detSizeA : x;
                x = x < 0 ? 0 : x;
                
                y = y > detSizeB ? detSizeB : y;
                y = y < 0 ? 0 : y;
                
                if(x == 0 || x == detSizeA || y == 0 || y == detSizeB){
                    pastParticles++;
                    continue;
                }
                
                i  = int(x / detCellSizeA );
                j  = int(y / detCellSizeB );
                
                idxG2 = IDX(i ,j ,0, detRes[0], detRes[1], 1);
                
                break;
                
            default:
                throw runtime_error("no such direction!");
        }
        
        gridMgr->addVectorVariableForDetector(idxG2, DOSE_FAR, 0, 1);
        
    }
    
    logger->writeMsg(("[Pusher] finished = "+to_string(finished)
                      +"\n pastParticles = "+to_string(pastParticles)
                      +"\n currentPartclNumOnDomain = "+to_string(currentPartclNumOnDomain)).c_str(), DEBUG);
    

}


void Pusher::push(){
    
   logger->writeMsg(("[Pusher] start pushing "+to_string(currentPartclNumOnDomain)
                                    +" particles...").c_str(), DEBUG);
    
    auto start_time = high_resolution_clock::now();
    
    double E[3], B[3];
    
    double dx = loader->spatialSteps[0];
    double dy = loader->spatialSteps[1];
    double dz = loader->spatialSteps[2];
    
    int xSize = loader->resolution[0];
    int ySize = loader->resolution[1];
    int zSize = loader->resolution[2];
    
    double domainShiftX = loader->boxCoordinates[0][0];
    double domainShiftY = loader->boxCoordinates[1][0];
    double domainShiftZ = loader->boxCoordinates[2][0];
    
    double F;
    double Fsquare;
    double Bsquare, G;
    
    double curV[3] = {0.0, 0.0, 0.0};
    double curU[3] = {0.0, 0.0, 0.0};
    double VBprod[3], UBprod[3];
    double new_position[3] = {0.0, 0.0, 0.0};
    double new_velocity[3] = {0.0, 0.0, 0.0};
    double ts = loader->getTimeStep();
    double qm;
    int coord, idx;
    
    double x, y, z;
    double lx4B, ly4B, lz4B, lx4E, ly4E, lz4E;
    int i4B, j4B, k4B, i4E, j4E, k4E;

    
    int idxG1, idxG2;
    int idx_x, idx_y, idx_z;
    double alpha, betta, gamma, weightE, weightB;
    double neighbourhood[8][3] = {{0,0,0}, {1,0,0}, {0,1,0}, {1,1,0},
                                  {0,0,1}, {1,0,1}, {0,1,1}, {1,1,1}};

    int type;
    
    VectorVar** Efield = gridMgr->getVectorVariableOnG2(ELECTRIC);
    VectorVar** Bfield = gridMgr->getVectorVariableOnG1(MAGNETIC);
    
    int G2nodesNumber = (xSize+2)*(ySize+2)*(zSize+2);
    int G1nodesNumber = (xSize+1)*(ySize+1)*(zSize+1);
    
    double* prtclPos;
    double* prtclVel;
    
    logger->writeMsg(("[Pusher] pushing "+to_string(currentPartclNumOnDomain)
                      +" particles...").c_str(), DEBUG);
    
    for( int idx=0; idx < currentPartclNumOnDomain; idx++){
        
        prtclPos = particles[idx]->getPosition();
        prtclVel = particles[idx]->getVelocity();
        type     = particles[idx]->getType();
        
        if(type == FINISHED){
            continue;
        }
        
//        if(trackList.count(idx)){
//            set<int>::iterator it = trackList.find(idx);
//            int trackedIdx = distance(trackList.begin(), it);
//            trackedParticles[3*trackedIdx + 0] = prtclPos[0];
//            trackedParticles[3*trackedIdx + 1] = prtclPos[1];
//            trackedParticles[3*trackedIdx + 2] = prtclPos[2];
//        }
//        
        for(int coord=0; coord < 3; coord++){
            E[coord] = 0.0;
            B[coord] = 0.0;
        }
        
        x  = prtclPos[0];
        y  = prtclPos[1];
        z  = prtclPos[2];
        
        if(isnan(x) || isnan(y)|| isnan(z)){
            logger->writeMsg(("[Pusher] problem  x = "+to_string(x)
                              +" y = "+to_string(y)+" z = "+to_string(z)).c_str(), DEBUG);
        }

        
        x = (x - domainShiftX)/dx;
        y = (y - domainShiftY)/dy;
        z = (z - domainShiftZ)/dz;
        
        i4B  = int(x);
        j4B  = int(y);
        k4B  = int(z);
        
        if(i4B<0 || j4B<0 || k4B<0){
            logger->writeMsg(("[Pusher] problem  i4B = "+to_string(i4B)
                              +" j4B = "+to_string(j4B)+" k4B = "+to_string(k4B)).c_str(), DEBUG);
            
        }
        
        lx4B = x-i4B;
        ly4B = y-j4B;
        lz4B = z-k4B;
        
        double alphasB[8] = {1.0-lx4B, lx4B, 1.0-lx4B, lx4B,
                             1.0-lx4B, lx4B, 1.0-lx4B, lx4B};
        
        double bettasB[8] = {1.0-ly4B, 1.0-ly4B, ly4B, ly4B,
                             1.0-ly4B, 1.0-ly4B, ly4B, ly4B};
        
        double gammasB[8] = {1.0-lz4B, 1.0-lz4B, 1.0-lz4B, 1.0-lz4B,
                                 lz4B,     lz4B,     lz4B,     lz4B};
        
        x  = prtclPos[0];
        y  = prtclPos[1];
        z  = prtclPos[2];
        
        x = (x - domainShiftX)/dx+0.5;
        y = (y - domainShiftY)/dy+0.5;
        z = (z - domainShiftZ)/dz+0.5;
        
        i4E  = int(x);
        j4E  = int(y);
        k4E  = int(z);
        
        if(i4E<0 || j4E<0 || k4E<0){
            logger->writeMsg(("[Pusher] problem  i4E = "+to_string(i4E)
                              +" j4E = "+to_string(j4E)).c_str(), DEBUG);
        
        }
        idxG2 = IDX(i4E ,j4E ,k4E, xSize+2, ySize+2, zSize+2);
        gridMgr->addVectorVariableForNodeG2(idxG2, DOSE, 0, 1);
        
        lx4E = x-i4E;
        ly4E = y-j4E;
        lz4E = z-k4E;
        
        double alphasE[8] = {1.0-lx4E, lx4E, 1.0-lx4E, lx4E,
                             1.0-lx4E, lx4E, 1.0-lx4E, lx4E};
        
        double bettasE[8] = {1.0-ly4E, 1.0-ly4E, ly4E, ly4E,
                             1.0-ly4E, 1.0-ly4E, ly4E, ly4E};
        
        double gammasE[8] = {1.0-lz4E, 1.0-lz4E, 1.0-lz4E, 1.0-lz4E,
                                 lz4E,     lz4E,     lz4E,     lz4E};

        
        for (int neigh_num=0; neigh_num < 8; neigh_num++){
        
            // - set b field on g1
            idx_x = i4B + neighbourhood[neigh_num][0];
            idx_y = j4B + neighbourhood[neigh_num][1];
            idx_z = k4B + neighbourhood[neigh_num][2];
            
            idxG1 = IDX(idx_x ,idx_y ,idx_z, xSize+1, ySize+1, zSize+1);
            
            alpha = lx4B;
            betta = ly4B;
            gamma = lz4B;
            
            alpha = alphasB[neigh_num];
            betta = bettasB[neigh_num];
            gamma = gammasB[neigh_num];
            weightB = alpha*betta*gamma;
            
            // - set e field on g2
            idx_x = i4E + neighbourhood[neigh_num][0];
            idx_y = j4E + neighbourhood[neigh_num][1];
            idx_z = k4E + neighbourhood[neigh_num][2];
            
            idxG2 = IDX(idx_x ,idx_y ,idx_z, xSize+2, ySize+2, zSize+2);
            
            alpha = alphasE[neigh_num];
            betta = bettasE[neigh_num];
            gamma = gammasE[neigh_num];
            weightE = alpha*betta*gamma;
            if(idxG2>=(xSize+2)*(ySize+2)*(zSize+2)){
                logger->writeMsg(("[Pusher] problem  weightE = "+to_string(weightE)
                                  +" idxG2 = "+to_string(idxG2)).c_str(), DEBUG);
            }
            

 
            const double* ef = Efield[idxG2]->getValue();
            const double* bf = Bfield[idxG1]->getValue();
            for(int coord=0; coord < 3; coord++){
                E[coord] += weightE*ef[coord];
                B[coord] += weightB*bf[coord];
            }
        }
        
        qm = charge/mass;
        F = 0.5*qm*ts;
        Fsquare = F*F;
        Bsquare = B[0]*B[0]+B[1]*B[1]+B[2]*B[2];
        G = 2.0/(1.0+Bsquare*Fsquare);
        
        /* __ half acceleration in e field __ */
        curV[0] = prtclVel[0]+F*E[0];
        curV[1] = prtclVel[1]+F*E[1];
        curV[2] = prtclVel[2]+F*E[2];
        
        /* __ half rotation in b field __ */
        VBprod[0] = curV[1] * B[2] - curV[2] * B[1];
        VBprod[1] = curV[2] * B[0] - curV[0] * B[2];
        VBprod[2] = curV[0] * B[1] - curV[1] * B[0];
        
        curU[0] = curV[0] + F*VBprod[0];
        curU[1] = curV[1] + F*VBprod[1];
        curU[2] = curV[2] + F*VBprod[2];
        
        UBprod[0] = curU[1] * B[2] - curU[2] * B[1];
        UBprod[1] = curU[2] * B[0] - curU[0] * B[2];
        UBprod[2] = curU[0] * B[1] - curU[1] * B[0];
        
        for( coord=0; coord < 3; coord++){
            new_velocity[coord] = curV[coord]+(G*UBprod[coord]+E[coord])*F;
            new_position[coord] = prtclPos[coord] + new_velocity[coord]*ts;
            if(isnan(UBprod[coord])){
                logger->writeMsg(("[Pusher] Nan velo "+to_string(new_velocity[coord])
                                  +" coord = "+to_string(coord)
                                   +" UBprod[coord] = "+to_string(UBprod[coord])
                                   +" E[coord] = "+to_string(E[coord])
                                   +" G = "+to_string(G)
                                   +" F = "+to_string(F)
                                   +" type = "+to_string(type)).c_str(), DEBUG);
                throw runtime_error("nan velo");
            }
            particles[idx]->setPosition(coord, new_position[coord]);
            particles[idx]->setVelocity(coord, new_velocity[coord]);
        }
//         particles[idx]->setPosition(2, new_position[2]);
        
        int domainNum = boundaryMgr->isPtclOutOfDomain(new_position);
        if( domainNum != IN){
            prtclPos = particles[idx]->getPosition();
            double _pos[3] = {prtclPos[0],prtclPos[1],prtclPos[2]};
//            boundaryMgr->storeParticle(idx, _pos);
            particles[idx]->setType(FINISHED);
        }
        
    }
    
    
//    gridMgr->sendBoundary2Neighbor(DOSE);
    
    auto end_time = high_resolution_clock::now();
//    string msgs ="[Pusher] solve(): stage 0 duration = "
//                    +to_string(duration_cast<milliseconds>(end_time - start_time).count())+" ms";
//    logger->writeMsg(msgs.c_str(),  DEBUG);
//    
//    vector<shared_ptr<Particle>> particles2add;
//    particles2add.reserve(EXPECTED_NUM_OF_PARTICLES);
//        
//    boundaryMgr->applyBC(particles, particles2add);
//    vector<double> leavingParticles = boundaryMgr->getLeavingParticlesIdxs();
//    int tot2add = particles2add.size();
//    int tot2remove = leavingParticles.size();
//    int prevNumOfPartcl = currentPartclNumOnDomain;
//        
//        for(int i=0; i<tot2add; i++){
//            if(i<tot2remove){
//                particles[int(leavingParticles[i])]->reinitializeUsingParticle(particles2add[i]);
//            }else{
//                particles[currentPartclNumOnDomain]->reinitializeUsingParticle(particles2add[i]);
//                currentPartclNumOnDomain++;
//            }
//           
//        }
//    
//    //need to start from the end otherwise there is a risk
//    //to damage particle inside with a particle that has to be removed from the end
//    for(int i=tot2remove-1; i>=tot2add; i--){
//            int idxLeave = int(leavingParticles[i]);
//            int idxToUse = currentPartclNumOnDomain-1;// last particle
//            if(idxLeave == idxToUse){
//                currentPartclNumOnDomain--;
//            }else{
//                particles[idxLeave]->reinitializeUsingParticle(particles[idxToUse]);
//                prtclPos = particles[idxLeave]->getPosition();
//                currentPartclNumOnDomain--;
//            }
//                
//     }
//    
//    boundaryMgr->reset();
//    particles2add.clear();
//    leavingParticles.clear();
//    int finished = boundaryMgr->getNumberOfFinished();
//    string msg003 ="[Pusher] tot2add = "+to_string(tot2add)
//                        +"; tot2remove = "+to_string(tot2remove)+"; currentPartclNumOnDomain = "
//                            +to_string(currentPartclNumOnDomain)+"; finished = "
//                            +to_string(finished);
//    logger->writeMsg(msg003.c_str(),  DEBUG);
//    
//    logger->writeMsg(("[Pusher] On Domain  "+to_string(currentPartclNumOnDomain)
//                                    +" particles...").c_str(), DEBUG);
//    int TOT_IN_BOX = 0;
//    MPI_Allreduce(&currentPartclNumOnDomain, &TOT_IN_BOX, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
//    logger->writeMsg(("[Pusher] total particles in the box = "+to_string(TOT_IN_BOX)).c_str(),  DEBUG);

    
    end_time = high_resolution_clock::now();
    string msg ="[Pusher] pushing...OK: solve() duration = "
                +to_string(duration_cast<milliseconds>(end_time - start_time).count())+" ms";
    logger->writeMsg(msg.c_str(),  DEBUG);
        
}
    
    
    
