//
//  ModelInitializer.hpp

#ifndef ModelInitializer_hpp
#define ModelInitializer_hpp

#include <chrono>
#include <stdio.h>
#include <iostream>
#include <string>
#include <hdf5.h>
#include <memory>

#include "../grid/GridManager.hpp"
#include "../input/Loader.hpp"
#include "../misc/Misc.hpp"
#include "../physics/pusher/Pusher.hpp"


class ModelInitializer
{
    
private:
    
    std::unique_ptr<Logger> logger;
    
    std::shared_ptr<Loader> loader;
    std::shared_ptr<GridManager> gridMng;
    std::shared_ptr<Pusher> pusher;
    
    int totParticleNumberInDomain;
    
    void initElectroMagneticField();
    void initParticles();
    
public:
    ModelInitializer(std::shared_ptr<Loader>,
                     std::shared_ptr<GridManager>,
                     std::shared_ptr<Pusher>);
    
    void initialize();
    
    void finilize();
    ~ModelInitializer();
};
#endif
