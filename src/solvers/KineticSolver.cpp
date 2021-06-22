#include "KineticSolver.hpp"



using namespace std;
using namespace chrono;




KineticSolver::KineticSolver(shared_ptr<Loader> load,
                           shared_ptr<GridManager> grid,
                           shared_ptr<Pusher> pshr):loader(move(load)),
                            gridMng(move(grid)), pusher(move(pshr))
{
    logger.reset(new Logger());
    
    initialize();
    logger->writeMsg("[KineticSolver] create...OK", DEBUG);
}

void KineticSolver::initialize()
{
    logger->writeMsg("[KineticSolver] initialize...OK", DEBUG);
}


int KineticSolver::solve(double time)
{
    int rank ;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    try{
        pusher->push();
        
    }catch(...){
        return SOLVE_FAIL;
    }

    if (rank == 0){
        double ts = loader->getTimeStep();
        string msg ="[KineticSolver] SOLVER step ="
                        +to_string(time)+"; time = "+to_string(time*ts);
        logger->writeMsg(msg.c_str(), DEBUG);
    }
    
    return SOLVE_OK;
}



KineticSolver::~KineticSolver(){
    finilize();
    logger->writeMsg("[KineticSolver] delete...OK", DEBUG);
}

void KineticSolver::finilize(){
    
}
