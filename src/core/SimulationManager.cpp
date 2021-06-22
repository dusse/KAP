#include "SimulationManager.hpp"
#include <iostream>
#include <string>
#include <thread>
using namespace std;
using namespace chrono;


SimulationManager::SimulationManager(int ac, char **av) {
    
    logger.reset(new Logger());
    
    int rank ;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    //init random
    srand(time(NULL) + rank);
    
    if(rank == 0){
        char buff[20];
        time_t now = time(NULL);
        strftime(buff, 20, "%Y-%m-%d %H:%M:%S", localtime(&now));
        
        logger->writeMsg("*****************************************************", INFO);
        logger->writeMsg("****                                             ****", INFO);
        logger->writeMsg(string("****     START  "+string(buff)+"              ****").c_str(), INFO);
        logger->writeMsg("****                                             ****", INFO);
        logger->writeMsg("****       0  A       A       PPPPPP             ****", INFO);
        logger->writeMsg("****       O  A      A A      P    P             ****", INFO);
        logger->writeMsg("****       O A      A   A     P    P             ****", INFO);
        logger->writeMsg("****       0A      AAAAAAA    PPPPPP             ****", INFO);
        logger->writeMsg("****       O A    A       A   P                  ****", INFO);
        logger->writeMsg("****       O  A  A         A  P                  ****", INFO);
        logger->writeMsg("****    Kinetic Algorithm for Proton-radiography ****", INFO);
        logger->writeMsg("*****************************************************", INFO);
        
    }
    
    
    
    if (ac > 1){
        setenv("PYTHONPATH", av[1] , 1);
    }else{
        if(rank == 0){
            string msg ="[SimulationManager] ATTENTION!!! you use default Initializer.py from src/input/ !!!";
            logger->writeMsg(msg.c_str(), CRITICAL);
        }
        const char  *PRJ_PATH = "src/input/";
        setenv("PYTHONPATH", PRJ_PATH , 1);
    }
    initialize();
}


void SimulationManager::initialize() {
    
    loader.reset(new Loader());
    loader->load();
    gridMng.reset(new GridManager(loader));
    boundMng.reset(new BoundaryManager(loader));
    
    pusher.reset(new Pusher(loader, gridMng, boundMng));
    initMng.reset(new ModelInitializer(loader, gridMng, pusher));
    
    solver.reset(new KineticSolver(loader, gridMng, pusher));
    writer.reset(new Writer(loader, gridMng, pusher));
    
    logger->writeMsg("[SimulationManager] init...OK", DEBUG);
}

void SimulationManager::runSimulation(int ac, char **av) {
    logger->writeMsg("[SimulationManager] launch simulation...OK", DEBUG);
    auto start_time_tot = high_resolution_clock::now();
    
    int maxTimeStep = loader->getMaxTimestepsNum();
    int maxTimeStep2Write = loader->getTimestepsNum2Write();
    int fileNumCount = 0;
    int i_time;
    int STOP_SIMULATION = 1;

    int rank ;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
  
    for ( i_time=0; i_time < maxTimeStep; i_time++ ){
        auto start_time = high_resolution_clock::now();
        if(rank == 0){
            logger->writeMsg("*****************************************************", INFO);
            logger->writeMsg("****                                             ****", INFO);
            logger->writeMsg(string("[SimulationManager] step = "+to_string(i_time)).c_str(), INFO);
            logger->writeMsg(string("[SimulationManager] time = "
                                        +to_string(i_time*loader->getTimeStep())).c_str(), INFO);
            logger->writeMsg("****                                             ****", INFO);
            logger->writeMsg("*****************************************************", INFO);
        }
        if( i_time % maxTimeStep2Write == 0 && i_time != 0){
            
            pusher->applyFreeSpace();
            writer->write(fileNumCount);
            fileNumCount++;
        }
        

        if(STOP_SIMULATION == 0){
             break;
        }
        
        
        if( solver->solve(i_time) == SOLVE_FAIL ){
            STOP_SIMULATION = 0;
            
            if(rank == 0){
                logger->writeMsg("*****************************************************", CRITICAL);
                logger->writeMsg("****                                             ****", CRITICAL);
                logger->writeMsg("****  STOP SIMULATION!!!    UNEXPECTED PROBLEM   ****", CRITICAL);
                logger->writeMsg("****                                             ****", CRITICAL);
                logger->writeMsg("****   try debug mode  'make -DLOG'              ****", CRITICAL);
                logger->writeMsg("*****************************************************", CRITICAL);
            }
        }
        
        if(rank == 0){
            auto end_time = high_resolution_clock::now();
            string msg = "[SimulationManager] Step duration = "
                +to_string(duration_cast<milliseconds>(end_time - start_time).count())+" ms";
            logger->writeMsg(msg.c_str(), INFO);
            logger->writeMsg("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx", INFO);
        }
    }
    
    
    
    if(rank == 0){
        auto end_time_tot = high_resolution_clock::now();
        logger->writeMsg("*****************************************************", INFO);
        logger->writeMsg("****                                             ****", INFO);
        logger->writeMsg("****          FINISH SIMULATION                  ****", INFO);
        string msg = "**** Total duration for "+to_string(i_time)+" steps : "
            +to_string(duration_cast<minutes>(end_time_tot - start_time_tot).count())+" min";
        logger->writeMsg(msg.c_str(), INFO);
        logger->writeMsg("****                                             ****", INFO);
        logger->writeMsg("*****************************************************", INFO);
    }

}



void SimulationManager::finilize() {
    logger->writeMsg("[SimulationManager] finalize...OK", DEBUG);
}

SimulationManager::~SimulationManager() {
    finilize();
}

