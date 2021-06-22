#include "GridManager.hpp"

using namespace std;
using namespace chrono;

GridManager::GridManager(shared_ptr<Loader> loader){
    this->loader=loader;
    logger.reset(new Logger());
    initialize();
    string msg ="[GridManager] init...OK {Node number:"+to_string(totalNodeNumber)+"}";
    logger->writeMsg(msg.c_str(), DEBUG);
}


GridManager::~GridManager(){
    delete nodesG2vars;
    delete nodesG1vars;
    delete nodesDetvars;
    
    for (int t = 0; t < 27; t++) {
        delete [] sendIdx[t];
        delete [] recvIdx[t];
    }
}


void GridManager::initialize(){
    logger->writeMsg("[GridManager] initialize() ...", DEBUG);
    
    int xRes = loader->resolution[0];
    int yRes = loader->resolution[1];
    int zRes = loader->resolution[2];
    
    this->totalNodeNumber = xRes*yRes*zRes;
    
    this->G1nodesNumber = (xRes+1)*(yRes+1)*(zRes+1);
    this->G2nodesNumber = (xRes+2)*(yRes+2)*(zRes+2);
    
    double Lx = loader->boxSizes[0], Ly = loader->boxSizes[1];
    
    detLxRes = int(xRes*loader->detectorSize[0]/Lx);
    detLyRes = int(yRes*loader->detectorSize[1]/Ly);

//    detLxRes = int(5000);
//    detLyRes = int(5000);
    
    this->detectorNodesNumber = (detLxRes)*(detLyRes);

    
    nodesG2vars = new VectorVar*[G2nodesNumber*SIZEG2];
    nodesG1vars = new VectorVar*[G1nodesNumber*SIZEG1];
    
    nodesDetvars = new VectorVar*[detectorNodesNumber*SIZEDET];
    
    initG1Nodes();
    initG2Nodes();
    
    initDetectorNodes();

    sendRecvIndecis4MPI();
    sendRecvIndecis4MPIonG4();
    
    getVectorVariablesForAllNodes();
    logger->writeMsg("[GridManager] initialize() ...OK", DEBUG);
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

//G1 grid describes vecticies of the cells
void GridManager::initG1Nodes(){
    auto start_time = high_resolution_clock::now();
    logger->writeMsg("[GridManager] init Node on G1 start ", DEBUG);
    int xRes = loader->resolution[0];
    int yRes = loader->resolution[1];
    int zRes = loader->resolution[2];
    int xResG1 = xRes+1, yResG1 = yRes+1, zResG1 = zRes+1;
    
    int i,j,k,idx;
    int idx_count = 0;
    for ( i=0; i<xResG1; i++){
        for ( j=0; j<yResG1; j++){
            for ( k=0; k<zResG1; k++){
                idx   = IDX(i,j,k,xResG1,yResG1,zResG1);
                VectorVar* mag = new VectorVar(MAGNETIC, {0.0, 0.0, 0.0});
                nodesG1vars[G1nodesNumber*MAGNETIC+idx] = mag;
                
            }
        }
    }
    
    auto end_time = high_resolution_clock::now();
    string msg ="[GridManager] init Node G1 OK..duration = "
                +to_string(duration_cast<milliseconds>(end_time - start_time).count())+" ms";
    logger->writeMsg(msg.c_str(), DEBUG);
}


//detector grid describes
void GridManager::initDetectorNodes(){
    auto start_time = high_resolution_clock::now();
    logger->writeMsg("[GridManager] init Node on detector start ", DEBUG);
    
    for ( int idx=0; idx<detectorNodesNumber; idx++){
        VectorVar* dose_far = new VectorVar(DOSE_FAR, {0.0});
        nodesDetvars[idx] = dose_far;
    }
    
    auto end_time = high_resolution_clock::now();
    string msg ="[GridManager] init Node on detector OK..duration = "
    +to_string(duration_cast<milliseconds>(end_time - start_time).count())+" ms";
    logger->writeMsg(msg.c_str(), DEBUG);
}


void GridManager::addVectorVariableForDetector(int idx, int name, int dim, double value){
    nodesDetvars[detectorNodesNumber*name+idx]->addValue(dim, value);
}


vector<VectorVar> GridManager::getDetectorDoseVar(){
    
    vector<VectorVar> result;
    result.reserve(detectorNodesNumber);
    
    for ( int idx=0; idx<detectorNodesNumber; idx++){
        result.push_back(*nodesDetvars[idx]);
    }
    
    return result;
}

vector<int> GridManager::getDetectorResolution(){
    return {detLxRes, detLyRes};
}

//G2 grid describes extended grid G1
void GridManager::initG2Nodes(){
    auto start_time = high_resolution_clock::now();
    logger->writeMsg("[GridManager] init Node start ", DEBUG);
    int xRes = loader->resolution[0], yRes = loader->resolution[1], zRes = loader->resolution[2];
    int xResG2 = xRes+2, yResG2 = yRes+2, zResG2 = zRes+2;
    
    int i,j,k,idx, idxG1;
    for ( i=0; i<xResG2; i++){
        for ( j=0; j<yResG2; j++){
            for ( k=0; k<zResG2; k++){
                idx = IDX(i,j,k,xResG2,yResG2,zResG2);
                VectorVar* ele = new VectorVar(ELECTRIC, {0.0, 0.0, 0.0});
                VectorVar* dose = new VectorVar(DOSE, {0.0, 0.0});
                nodesG2vars[G2nodesNumber*ELECTRIC+idx] = ele;
                nodesG2vars[G2nodesNumber*DOSE+idx] = dose;
                
            }
        }
    }
    auto end_time = high_resolution_clock::now();
    string msg ="[GridManager] init Node G2 duration = "
                +to_string(duration_cast<milliseconds>(end_time - start_time).count())+" ms";
    logger->writeMsg(msg.c_str(), DEBUG);
}


void GridManager::sendRecvIndecis4MPI(){
    int xRes = loader->resolution[0], yRes = loader->resolution[1], zRes = loader->resolution[2];
    int xResG2 = xRes+2, yResG2 = yRes+2, zResG2 = zRes+2;
    int idx, i, j, k, a, b, c, t;
    vector<int> idxX, idxY, idxZ;
    
    for (a = -1; a <= 1; a++) {
        for (b = -1; b <= 1; b++) {
            for (c = -1; c <= 1; c++) {
                t = (1 + c) + 3 * ((1 + b) + 3 * (1 + a));

                switch (a) {
                    case -1: idxX = {1   , 1   , xRes + 1, xRes + 1}; break;
                    case  0: idxX = {1   , xRes, 1       , xRes    }; break;
                    case  1: idxX = {xRes, xRes, 0       , 0       }; break;
                }
                switch (b) {
                    case -1: idxY = {1   , 1   , yRes + 1, yRes + 1}; break;
                    case  0: idxY = {1   , yRes, 1       , yRes    }; break;
                    case  1: idxY = {yRes, yRes, 0       , 0       }; break;
                }
                switch (c) {
                    case -1: idxZ = {1   , 1   , zRes + 1, zRes + 1}; break;
                    case  0: idxZ = {1   , zRes, 1       , zRes    }; break;
                    case  1: idxZ = {zRes, zRes, 0       , 0       }; break;
                }
                
                counter[t] = (idxX[1] - idxX[0] + 1)*
                             (idxY[1] - idxY[0] + 1)*
                             (idxZ[1] - idxZ[0] + 1);
                sendIdx[t] = new int[counter[t]];
                recvIdx[t] = new int[counter[t]];
                
                int lineIDX = 0;
                for (i = idxX[0]; i <= idxX[1]; i++) {
                    for (j = idxY[0]; j <= idxY[1]; j++) {
                        for (k = idxZ[0]; k <= idxZ[1]; k++) {
                            idx = IDX(i,j,k, xResG2, yResG2, zResG2);
                            sendIdx[t][lineIDX] = idx;
                            lineIDX++;
                        }
                    }
                }

                lineIDX = 0;
                for (i = idxX[2]; i <= idxX[3]; i++) {
                    for (j = idxY[2]; j <= idxY[3]; j++) {
                        for (k = idxZ[2]; k <= idxZ[3]; k++) {
                            idx = IDX(i,j,k, xResG2, yResG2, zResG2);
                            recvIdx[t][lineIDX] = idx;
                            lineIDX++;
                        }
                    }
                }
            }
        }
    }
}


void GridManager::sendBoundary2Neighbor(int varName){
    
    auto start_time = high_resolution_clock::now();
    int varDim = nodesG2vars[G2nodesNumber*varName]->getSize();
    int t, i;
    double *sendBuf[27];
    double *recvBuf[27];
    const double* vectorVar;
    
    for (t = 0; t < 27; t++) {
        sendBuf[t] = new double[counter[t]*varDim*sizeof(double)];
        recvBuf[t] = new double[counter[t]*varDim*sizeof(double)];
        for(i = 0; i < counter[t]; i++){
            vectorVar = nodesG2vars[G2nodesNumber*varName+sendIdx[t][i]]->getValue();
            for (int dim=0;dim<varDim;dim++) {
                sendBuf[t][varDim*i+dim] = vectorVar[dim];
            }
        }
    }
    
    MPI_Status st;
    for (t = 0; t < 27; t++) {
        if (t != 13){
            MPI_Sendrecv(sendBuf[t], counter[t]*varDim, MPI_DOUBLE, loader->neighbors2Send[t], t,
                         recvBuf[t], counter[t]*varDim, MPI_DOUBLE, loader->neighbors2Recv[t], t,
                         MPI_COMM_WORLD, &st);
        }
    }
    
    for (t = 0; t < 27; t++) {
        if (t != 13 && loader->neighbors2Recv[t] != MPI_PROC_NULL) {
            for(i = 0; i < counter[t]; i++){
                for (int dim=0;dim<varDim;dim++) {
                    nodesG2vars[G2nodesNumber*varName+recvIdx[t][i]]
                    ->setValue( dim, recvBuf[t][varDim*i+dim]);
                }
            }
        }
    }
    
    
    for (t = 0; t < 27; t++) {
        delete [] sendBuf[t];
        delete [] recvBuf[t];
    }
    
    auto end_time = high_resolution_clock::now();
    string msg ="[GridManager] send boundary for var = "+to_string(varName)
                +" duration = "
                +to_string(duration_cast<milliseconds>(end_time - start_time).count())+" ms";
    logger->writeMsg(msg.c_str(), DEBUG);
}




void GridManager::sendRecvIndecis4MPIonG4(){
    int xRes = loader->resolution[0],
    yRes = loader->resolution[1],
    zRes = loader->resolution[2];
    int xResG2 = xRes+2, yResG2 = yRes+2, zResG2 = zRes+2;
    int xResG4 = xRes+4, yResG4 = yRes+4, zResG4 = zRes+4;
    int idx, i, j, k, a, b, c, t;
    vector<int> idxX, idxY, idxZ;
    
    for (a = -1; a <= 1; a++) {
        for (b = -1; b <= 1; b++) {
            for (c = -1; c <= 1; c++) {
                
                switch (a) {
                    case -1: idxX = {3   , 3       , xRes + 3, xRes + 3}; break;
                    case  0: idxX = {1   , xRes + 2, 1       , xRes + 2}; break;
                    case  1: idxX = {xRes, xRes    , 0       , 0       }; break;
                }
                switch (b) {
                    case -1: idxY = {3   , 3       , yRes + 3, yRes + 3}; break;
                    case  0: idxY = {1   , yRes + 2, 1       , yRes + 2}; break;
                    case  1: idxY = {yRes, yRes    , 0       , 0       }; break;
                }
                switch (c) {
                    case -1: idxZ = {3   , 3       , zRes + 3, zRes + 3}; break;
                    case  0: idxZ = {1   , zRes + 2, 1       , zRes + 2}; break;
                    case  1: idxZ = {zRes, zRes    , 0       , 0       }; break;
                }
                
                t = (1 + c) + 3 * ((1 + b) + 3 * (1 + a));
                
                counterOnG4[t] = (idxX[1] - idxX[0] + 1)*
                (idxY[1] - idxY[0] + 1)*
                (idxZ[1] - idxZ[0] + 1);
                sendIdxOnG4[t] = new int[counterOnG4[t]];
                recvIdxOnG4[t] = new int[counterOnG4[t]];
                
                int lineIDX = 0;
                for (i = idxX[0]; i <= idxX[1]; i++) {
                    for (j = idxY[0]; j <= idxY[1]; j++) {
                        for (k = idxZ[0]; k <= idxZ[1]; k++) {
                            idx = IDX(i,j,k, xResG4, yResG4, zResG4);
                            sendIdxOnG4[t][lineIDX] = idx;
                            lineIDX++;
                        }
                    }
                }
                lineIDX = 0;
                for (i = idxX[2]; i <= idxX[3]; i++) {
                    for (j = idxY[2]; j <= idxY[3]; j++) {
                        for (k = idxZ[2]; k <= idxZ[3]; k++) {
                            idx = IDX(i,j,k, xResG4, yResG4, zResG4);
                            recvIdxOnG4[t][lineIDX] = idx;
                            lineIDX++;
                        }
                    }
                }
            }
        }
    }
}



void GridManager::smooth(int varName){
    
    int i, j, k, idx, idxG4;
    int xRes = loader->resolution[0],
    yRes = loader->resolution[1],
    zRes = loader->resolution[2];
    int xResG4 = xRes+4, yResG4 = yRes+4, zResG4 = zRes+4;
    int xResG2 = xRes+2, yResG2 = yRes+2, zResG2 = zRes+2;
    int totG4 = xResG4*yResG4*zResG4;
    auto start_time = high_resolution_clock::now();
    int varDim = nodesG2vars[G2nodesNumber*varName]->getSize();
    int t;
    double *sendBuf[27];
    double *recvBuf[27];
    const double* vectorVar;
    
    double* varValues = new double[totG4*varDim*sizeof(double)];
    
    for ( i=0; i<xResG2; i++){
        for ( j=0; j<yResG2; j++){
            for ( k=0; k<zResG2; k++){
                idx   = IDX(i  ,j  ,k  ,xResG2,yResG2,zResG2);
                idxG4 = IDX(i+1,j+1,k+1,xResG4,yResG4,zResG4);
                vectorVar = nodesG2vars[G2nodesNumber*varName+idx]->getValue();
                for (int dim=0;dim<varDim;dim++) {
                    varValues[varDim*idxG4+dim] = vectorVar[dim];
                }
            }
        }
    }
    
    for (t = 0; t < 27; t++) {
        sendBuf[t] = new double[counterOnG4[t]*varDim*sizeof(double)];
        recvBuf[t] = new double[counterOnG4[t]*varDim*sizeof(double)];
        for(i = 0; i < counterOnG4[t]; i++){
            int si = sendIdxOnG4[t][i];
            for (int dim=0;dim<varDim;dim++) {
                sendBuf[t][varDim*i+dim] = varValues[varDim*si+dim];
                
            }
        }
    }
    
    MPI_Status st;
    for (t = 0; t < 27; t++) {
        if (t != 13){
            MPI_Sendrecv(sendBuf[t], counterOnG4[t]*varDim, MPI_DOUBLE, loader->neighbors2Send[t], t,
                         recvBuf[t], counterOnG4[t]*varDim, MPI_DOUBLE, loader->neighbors2Recv[t], t,
                         MPI_COMM_WORLD, &st);
        }
    }
    for (t = 0; t < 27; t++) {
        if (t != 13 && loader->neighbors2Recv[t] != MPI_PROC_NULL) {
            for(i = 0; i < counterOnG4[t]; i++){
                int ri = recvIdxOnG4[t][i];
                for (int dim=0;dim<varDim;dim++) {
                    varValues[varDim*ri+dim] = recvBuf[t][varDim*i+dim];
                }
            }
        }
    }
    
    
    for (t = 0; t < 27; t++) {
        delete [] sendBuf[t];
        delete [] recvBuf[t];
    }
    
    int zeroOrderNeighb[6][3]  =
    {{-1,0 ,0 }, {+1,0 ,0 },
        {0 ,-1,0 }, {0 ,+1,0 },
        {0 ,0 ,-1}, {0 ,0 ,+1}};
    
    int firstOrderNeighb[12][3] =
    {{-1,-1, 0}, {-1,+1, 0},
        {-1, 0,-1}, {-1, 0,+1},
        {0 ,-1,-1}, { 0,-1,+1},
        {0 ,+1,-1}, { 0,+1,+1},
        {+1, 0,-1}, {+1, 0,+1},
        {+1,-1, 0}, {+1,+1, 0}};
    
    int secndOrderNeighb[8][3] =
    {{-1,-1,-1}, {-1,-1,+1},
        {-1,+1,-1}, {-1,+1,+1},
        {+1,-1,-1}, {+1,-1,+1},
        {+1,+1,-1}, {+1,+1,+1}};
    
    const double k2 = 0.125;// 1/8
    const double k3 = 0.0625;// 1/16
    const double k4 = 0.03125;// 1/32
    const double k5 = 0.015625;// 1/64
    //                    +--------
    //  0.125*1+0.0625*6+0.03125*12+0.015625*8 = 1
    
    int foi, soi, toi;
    int idxFo, idxSo, idxTo;
    
    for ( i=1; i<xResG2+1; i++){
        for ( j=1; j<yResG2+1; j++){
            for ( k=1; k<zResG2+1; k++){
                idxG4 = IDX(i,j,k,xResG4,yResG4,zResG4);
                
                for (int dim=0;dim<varDim;dim++) {
                    
                    double smoothedVal = k2*varValues[varDim*idxG4+dim];
                    
                    for(foi = 0; foi<6; foi++){
                        idxFo = IDX(i+zeroOrderNeighb[foi][0],
                                    j+zeroOrderNeighb[foi][1],
                                    k+zeroOrderNeighb[foi][2],
                                    xResG4,yResG4,zResG4);
                        smoothedVal += k3*varValues[varDim*idxFo+dim];
                    }
                    
                    for(soi = 0; soi<12; soi++){
                        idxSo = IDX(i+firstOrderNeighb[soi][0],
                                    j+firstOrderNeighb[soi][1],
                                    k+firstOrderNeighb[soi][2],
                                    xResG4,yResG4,zResG4);
                        smoothedVal += k4*varValues[varDim*idxSo+dim];
                    }
                    
                    for(toi = 0; toi<8; toi++){
                        idxTo = IDX(i+secndOrderNeighb[toi][0],
                                    j+secndOrderNeighb[toi][1],
                                    k+secndOrderNeighb[toi][2],
                                    xResG4,yResG4,zResG4);
                        smoothedVal += k5*varValues[varDim*idxTo+dim];
                    }
                    idx = IDX(i-1,j-1,k-1,xResG2,yResG2,zResG2);
                    nodesG2vars[G2nodesNumber*varName+idx]->setValue(dim, smoothedVal);
                }
            }
        }
    }
    
    
    delete varValues;
    auto end_time = high_resolution_clock::now();
    string msg ="[GridManager] smooth "+to_string(varName)
    +" duration = "+to_string(duration_cast<milliseconds>(end_time - start_time).count())+" ms";
    logger->writeMsg(msg.c_str(), DEBUG);
    
}

void GridManager::smoothG1(int varName){
    
    int i, j, k, idx, idxG4;
    int xRes = loader->resolution[0],
    yRes = loader->resolution[1],
    zRes = loader->resolution[2];
    int xResG4 = xRes+4, yResG4 = yRes+4, zResG4 = zRes+4;
    int xResG2 = xRes+2, yResG2 = yRes+2, zResG2 = zRes+2;
    int xResG1 = xRes+1, yResG1 = yRes+1, zResG1 = zRes+1;
    int totG4 = xResG4*yResG4*zResG4;
    auto start_time = high_resolution_clock::now();
    int varDim = nodesG1vars[G1nodesNumber*varName]->getSize();
    int t;
    double *sendBuf[27];
    double *recvBuf[27];
    const double* vectorVar;
    
    double* varValues = new double[totG4*varDim*sizeof(double)];
    
    for ( i=1; i<xResG2; i++){
        for ( j=1; j<yResG2; j++){
            for ( k=1; k<zResG2; k++){
                idx   = IDX(i-1  ,j-1  ,k-1  ,xResG1,yResG1,zResG1);
                idxG4 = IDX(i+1,j+1,k+1,xResG4,yResG4,zResG4);
                vectorVar = nodesG1vars[G2nodesNumber*varName+idx]->getValue();
                for (int dim=0;dim<varDim;dim++) {
                    varValues[varDim*idxG4+dim] = vectorVar[dim];
                }
            }
        }
    }
    
    for ( i=0; i<1; i++){
        for ( j=0; j<1; j++){
            for ( k=0; k<1; k++){
                idxG4 = IDX(i+1,j+1,k+1,xResG4,yResG4,zResG4);
                for (int dim=0;dim<varDim;dim++) {
                    varValues[varDim*idxG4+dim] = 0.0;
                }
            }
        }
    }
    
    for (t = 0; t < 27; t++) {
        sendBuf[t] = new double[counterOnG4[t]*varDim*sizeof(double)];
        recvBuf[t] = new double[counterOnG4[t]*varDim*sizeof(double)];
        for(i = 0; i < counterOnG4[t]; i++){
            int si = sendIdxOnG4[t][i];
            for (int dim=0;dim<varDim;dim++) {
                sendBuf[t][varDim*i+dim] = varValues[varDim*si+dim];
                
            }
        }
    }
    
    MPI_Status st;
    for (t = 0; t < 27; t++) {
        if (t != 13){
            MPI_Sendrecv(sendBuf[t], counterOnG4[t]*varDim, MPI_DOUBLE, loader->neighbors2Send[t], t,
                         recvBuf[t], counterOnG4[t]*varDim, MPI_DOUBLE, loader->neighbors2Recv[t], t,
                         MPI_COMM_WORLD, &st);
        }
    }
    for (t = 0; t < 27; t++) {
        if (t != 13 && loader->neighbors2Recv[t] != MPI_PROC_NULL) {
            for(i = 0; i < counterOnG4[t]; i++){
                int ri = recvIdxOnG4[t][i];
                for (int dim=0;dim<varDim;dim++) {
                    varValues[varDim*ri+dim] = recvBuf[t][varDim*i+dim];
                }
            }
        }
    }
    
    
    for (t = 0; t < 27; t++) {
        delete [] sendBuf[t];
        delete [] recvBuf[t];
    }
    
    int zeroOrderNeighb[6][3]  =
    {{-1,0 ,0 }, {+1,0 ,0 },
        {0 ,-1,0 }, {0 ,+1,0 },
        {0 ,0 ,-1}, {0 ,0 ,+1}};
    
    int firstOrderNeighb[12][3] =
    {{-1,-1, 0}, {-1,+1, 0},
        {-1, 0,-1}, {-1, 0,+1},
        {0 ,-1,-1}, { 0,-1,+1},
        {0 ,+1,-1}, { 0,+1,+1},
        {+1, 0,-1}, {+1, 0,+1},
        {+1,-1, 0}, {+1,+1, 0}};
    
    int secndOrderNeighb[8][3] =
    {{-1,-1,-1}, {-1,-1,+1},
        {-1,+1,-1}, {-1,+1,+1},
        {+1,-1,-1}, {+1,-1,+1},
        {+1,+1,-1}, {+1,+1,+1}};
    
    const double k2 = 0.125;// 1/8
    const double k3 = 0.0625;// 1/16
    const double k4 = 0.03125;// 1/32
    const double k5 = 0.015625;// 1/64
    //                    +--------
    //  0.125*1+0.0625*6+0.03125*12+0.015625*8 = 1
    
    int foi, soi, toi;
    int idxFo, idxSo, idxTo;
    
    for ( i=2; i<xResG2+1; i++){
        for ( j=2; j<yResG2+1; j++){
            for ( k=2; k<zResG2+1; k++){
                idxG4 = IDX(i,j,k,xResG4,yResG4,zResG4);
                
                for (int dim=0;dim<varDim;dim++) {
                    
                    double smoothedVal = k2*varValues[varDim*idxG4+dim];
                    
                    for(foi = 0; foi<6; foi++){
                        idxFo = IDX(i+zeroOrderNeighb[foi][0],
                                    j+zeroOrderNeighb[foi][1],
                                    k+zeroOrderNeighb[foi][2],
                                    xResG4,yResG4,zResG4);
                        smoothedVal += k3*varValues[varDim*idxFo+dim];
                    }
                    
                    for(soi = 0; soi<12; soi++){
                        idxSo = IDX(i+firstOrderNeighb[soi][0],
                                    j+firstOrderNeighb[soi][1],
                                    k+firstOrderNeighb[soi][2],
                                    xResG4,yResG4,zResG4);
                        smoothedVal += k4*varValues[varDim*idxSo+dim];
                    }
                    
                    for(toi = 0; toi<8; toi++){
                        idxTo = IDX(i+secndOrderNeighb[toi][0],
                                    j+secndOrderNeighb[toi][1],
                                    k+secndOrderNeighb[toi][2],
                                    xResG4,yResG4,zResG4);
                        smoothedVal += k5*varValues[varDim*idxTo+dim];
                    }
                    idx = IDX(i-2,j-2,k-2,xResG1,yResG1,zResG1);
                    nodesG1vars[G1nodesNumber*varName+idx]->setValue(dim, smoothedVal);
                }
            }
        }
    }
    
    
    delete varValues;
    auto end_time = high_resolution_clock::now();
    string msg ="[GridManager] smooth G1 "+to_string(varName)
    +" duration = "+to_string(duration_cast<milliseconds>(end_time - start_time).count())+" ms";
    logger->writeMsg(msg.c_str(), DEBUG);
    
}


vector<vector<VectorVar>> GridManager::getVectorVariablesForAllNodes(){
    
    int idxG2, idxG1;
    int xRes = loader->resolution[0], yRes = loader->resolution[1], zRes = loader->resolution[2];
    int xResG1 = xRes+1, yResG1 = yRes+1, zResG1 = zRes+1;
    int xResG2 = xRes+2, yResG2 = yRes+2, zResG2 = zRes+2;
    vector<vector<VectorVar>> result;
    result.reserve(xRes*yRes*zRes);
    int i,j,k;
    smooth(DOSE);
    
    for ( i=0; i<xRes; i++){
        for ( j=0; j<yRes; j++){
            for ( k=0; k<zRes; k++){
                idxG2 = IDX(i+1,j+1,k+1,xResG2,yResG2,zResG2);
                vector<VectorVar> allVars;
                allVars.push_back(*nodesG2vars[G2nodesNumber*DOSE+idxG2]);
                allVars.push_back(*nodesG2vars[G2nodesNumber*ELECTRIC+idxG2]);
                
                idxG1 = IDX(i,j,k,xResG1,yResG1,zResG1);
                allVars.push_back(*nodesG1vars[G1nodesNumber*MAGNETIC+idxG1]);

                result.push_back(allVars);
            }
        }
    }
    
    return result;
}



void GridManager::setVectorVariableForNodeG1(int idx, VectorVar variable){
    nodesG1vars[G1nodesNumber*variable.getName()+idx]->setValue(variable.getValue());
}


void GridManager::setVectorVariableForNodeG2(int idx, VectorVar variable){
    nodesG2vars[G2nodesNumber*variable.getName()+idx]->setValue(variable.getValue());
}

void GridManager::setVectorVariableForNodeG2(int idx, int name, int dim, double value){
    nodesG2vars[G2nodesNumber*name+idx]->setValue(dim, value);
}

void GridManager::setVectorVariableForNodeG1(int idx, int name, int dim, double value){
    nodesG1vars[G1nodesNumber*name+idx]->setValue(dim, value);
}

void GridManager::addVectorVariableForNodeG2(int idx, int name, int dim, double value){
    nodesG2vars[G2nodesNumber*name+idx]->addValue(dim, value);
}


VectorVar** GridManager::getVectorVariableOnG1(int varName){
    return &nodesG1vars[G1nodesNumber*varName];
}


VectorVar** GridManager::getVectorVariableOnG2(int varName){
     return &nodesG2vars[G2nodesNumber*varName];
}

