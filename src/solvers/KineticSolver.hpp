//
//  KineticSolver.hpp

#ifndef KineticSolver_hpp
#define KineticSolver_hpp

#include <chrono>
#include <stdio.h>
#include <iostream>
#include <string>
#include <memory>

#include "../grid/GridManager.hpp"
#include "../input/Loader.hpp"
#include "../misc/Misc.hpp"

#include "../physics/pusher/Pusher.hpp"


const static int  SOLVE_OK   = 0;
const static int  SOLVE_FAIL = 1;

class KineticSolver{
    
private:
    
    std::unique_ptr<Logger> logger;
    
    std::shared_ptr<Loader> loader;
    std::shared_ptr<GridManager> gridMng;
    std::shared_ptr<Pusher> pusher;

    
public:
    KineticSolver(std::shared_ptr<Loader>,
                  std::shared_ptr<GridManager>,
                  std::shared_ptr<Pusher>);
        
    void initialize();
    int solve(double);
    void finilize();
    ~KineticSolver();
};
#endif
