#include "Logger.hpp"

#include <iostream>
#include <string>
using namespace std;



Logger::Logger()
{
}

void Logger::writeMsg(const char* input, int level)
{
    if (level <= MINIMAL_LEVEL){
        cout << "[" << level << "] " << input <<  "." << endl;
    }
}

void Logger::writeMsg(const char* input)
{
    int rank ;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0){
        cout << "[INFO] " << input <<  "." << endl;
    }
}



Logger::~Logger()
{
}


