#ifndef Loader_hpp
#define Loader_hpp

#include <Python.h>
#include <mpi.h>
#include <stdio.h>
//#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <math.h>
#include "../misc/Logger.hpp"
#include "../misc/Misc.hpp"

class Loader{
    
private:
    Logger logger;
    double timeStep;
    
    int maxTimestepsNum;
    int timestepsNum2Write;

    int ppc;
    
    std::string outputDir;
    std::string fileNameTemplate;

    PyObject *pInstance;
    
    double callPyFloatFunction( PyObject*, const std::string, const std::string);
    
    double callPyFloatFunctionWith3args( PyObject*, const std::string, const std::string,
                                                double, double, double);
    
    long callPyLongFunction( PyObject*, const std::string, const std::string);
    
    std::string callPyStringFunction( PyObject*, const std::string, const std::string);
    
    PyObject* getPyMethod(PyObject* , const std::string , const std::string);
    
public:
    
    int flowDirection;
    
    int resolution[3];
    int totPixelsPerBoxSide[3];
    int offsetInPixels[3];
    int mpiDomains[3];
    int mpiCoords[3];
    const int SIMULATION_SIZE = 3;
    int dim;
    double boxCoordinates[3][2];
    double boxSizes[3];
    double spatialSteps[3];
    double detectorSize[2];
    double sourcePosition;
    double detectorPosition;
    
       //MPI staff
    std::vector<int> neighbors2Send;//27
    std::vector<int> neighbors2Recv;//27
    
    
    Loader();
    ~Loader();
    void load();
    std::vector<double> getVelocity(double,double,double);
    double getTimeStep();
    int getMaxTimestepsNum();
    int getTimestepsNum2Write();
    double getMass();
    double getCharge();
    double getParticlesPerCellNumber();
    std::string getOutputDir() const;
    std::string getFilenameTemplate();
    PyObject * getPythonClassInstance(std::string className);
    
    std::vector<double> getBfield(double,double,double);
    std::vector<double> getEfield(double,double,double);
};
#endif /* Loader_hpp */
