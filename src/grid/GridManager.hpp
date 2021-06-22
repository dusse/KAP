#ifndef GridManager_hpp
#define GridManager_hpp
#include <stdio.h>
#include <vector>
#include <map>
#include <set>
#include <cmath>
#include <memory>
#include <chrono>
#include <mpi.h>
#include "../misc/Logger.hpp"
#include "../misc/Misc.hpp"
#include "../input/Loader.hpp"

#include "../common/variables/VectorVar.hpp"


enum G1VAR{
    MAGNETIC,
    SIZEG1
};

enum G2VAR{
    ELECTRIC,
    DOSE,
    SIZEG2
};

enum GDETVAR{
    DOSE_FAR,
    SIZEDET
};


class GridManager{
    
    
private:
    
    std::unique_ptr<Logger> logger;
    std::shared_ptr<Loader> loader;
    
    
    int G1nodesNumber;
    int G2nodesNumber;
    int totalNodeNumber;
    int detectorNodesNumber;
    int detLxRes, detLyRes;
    
    VectorVar** nodesG1vars;
    VectorVar** nodesG2vars;
    VectorVar** nodesDetvars;
    
    int counter[27], counterOnG4[27];
    int* sendIdx[27];
    int* recvIdx[27];
    int* sendIdxOnG4[27];
    int* recvIdxOnG4[27];
    
    void initialize();
    
    void initG1Nodes();
    void initG2Nodes();
    void initDetectorNodes();
    
    void sendRecvIndecis4MPI();
    void sendRecvIndecis4MPIonG4();
    
public:
    
    GridManager(std::shared_ptr<Loader>);
    ~GridManager();
    
    VectorVar** getVectorVariableOnG1(int);
    VectorVar** getVectorVariableOnG2(int);
    
    std::vector<VectorVar> getDetectorDoseVar();
    std::vector<int> getDetectorResolution();
    
    std::vector<std::vector<VectorVar>> getVectorVariablesForAllNodes();
    
    void setVectorVariableForNodeG1(int, VectorVar);
    void setVectorVariableForNodeG2(int, VectorVar);
    
    void setVectorVariableForNode(int, VectorVar);
    
    void setVectorVariableForNodeG1(int , int , int, double);
    void setVectorVariableForNodeG2(int , int , int, double);
    void addVectorVariableForNodeG2(int , int , int, double);
    void addVectorVariableForDetector(int , int , int, double);
    
    void sendBoundary2Neighbor(int);
    void smooth(int);
    void smoothG1(int);
    
    
};
#endif /* GridManager_hpp */
