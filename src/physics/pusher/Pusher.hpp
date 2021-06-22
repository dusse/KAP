#ifndef Pusher_hpp
#define Pusher_hpp

#include <stdio.h>

#include <chrono>
#include <stdio.h>
#include <iostream>
#include <cmath>
#include <string>
#include <memory>
#include <set>


#include "../../grid/GridManager.hpp"
#include "../../grid/boundary/BoundaryManager.hpp"

#include "../../particles/Particle.hpp"

#include "../../input/Loader.hpp"
#include "../../misc/Misc.hpp"


const int NUMBER_TRACKED_PARTICLES = 5;

class Pusher{
    
private:
    std::unique_ptr<Logger> logger;
    std::shared_ptr<Loader> loader;
    std::shared_ptr<GridManager> gridMgr;
    std::shared_ptr<BoundaryManager> boundaryMgr;
    
    double weight;
    double charge;
    double mass;
    
    int totinBoxInit = 0;
    int totalNum;
    
    std::set<int> trackList;
    
    
    
    int currentPartclNumOnDomain = 0;
    Particle** particles;
    
    void initialize();
    
    
public:
    
    double* trackedParticles;
    
    
    Pusher(std::shared_ptr<Loader>, std::shared_ptr<GridManager>, std::shared_ptr<BoundaryManager>);
    ~Pusher();
    void push();
    
    int getTotalParticleNumber();
    void initTrackedParticles();
    void setParticleWeight(double);
    void setParticleCharge(double);
    void setParticleMass(double);
    void initParticles(int);
    void addParticles(std::vector<std::shared_ptr<Particle>>);
    
    void setParticlePosition(int, double[6]);
    void setParticleVelocity(int, double[6]);
    void setParticleType(int, int);
    
    void applyFreeSpace();
    Particle** getParticles();
    

};


#endif
