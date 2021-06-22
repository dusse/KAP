#include "Loader.hpp"

using namespace std;

Loader::Loader(){
    logger.writeMsg("[Loader] Loader start!\n", DEBUG);
}

const char  *INIT_CLASS_NAME = "Initializer";
const string  BRACKETS ="()";
const string  BRACKETS_3DOUBLE = "(ddd)";
const string  GET = "get";
const string  GET_TIMESTEP = "getTimestep";
const string  GET_MAX_TIMESTEPS_NUM = "getMaxTimestepsNum";

const string  GET_TIMESTEP_WRITE = "getOutputTimestep";
const string  GET_OUTPUT_DIR = "getOutputDir";

const string  GET_FILENAME_TEMPLATE = "getOutputFilenameTemplate";
const string  GET_DENSITY = "getDensity";
const string  GET_NUM_OF_PARTICLES = "getParticlesPerCellNumber";

const string   GET_DETECTOR_POSITION  = "getDetectorPosition";
const string   GET_SOURCE_POSITION    = "getSourcePosition";

const string   GET_DETECTOR_SIZE  = "getDetectorSize";

const string  GET_VELOCITY = "getVelocity";
const string  GET_MASS   = "getMass";
const string  GET_CHARGE = "getCharge";

const string  GET_BFIELD = "getBfield";
const string  GET_EFIELD = "getEfield";

const string  GET_FLOW_DIRECTION = "getFlowDirection";

const string  GET_MPI_DOMAIN_NUM = "mpiDomainNum";



const string  dirs[] = {"X", "Y", "Z"};

void Loader::load()
{
        MPI_Comm com;
        int periods[3]; //1 - periodic for MPI
        periods[0] = 1;
        periods[1] = 1;
        periods[2] = 1;
    
        int rank ;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        int cs[3];
        int ss[3];


        string msg = "[Loader] initialization...for proc = " + to_string(rank);
        logger.writeMsg(msg.c_str(), DEBUG);
        Py_Initialize();
        this->pInstance = getPythonClassInstance(INIT_CLASS_NAME);


        for (int n=0; n<3; n++){
             string domainNum =GET+dirs[n]+GET_MPI_DOMAIN_NUM;
             this->mpiDomains[n] = (int) callPyLongFunction( pInstance, domainNum, BRACKETS );
        }
    
        MPI_Cart_create(MPI_COMM_WORLD, 3, mpiDomains, periods, 1, &com);
        MPI_Cart_coords(com, rank, 3, cs);
        for(int i = 0; i < 3; i++){
            this->mpiCoords[i] = cs[i];
        }
        int wr;

        for(int i = 0; i < 27; i++){
            this->neighbors2Send.push_back(MPI_PROC_NULL);
            this->neighbors2Recv.push_back(MPI_PROC_NULL);
        }
        for (int a = -1; a <= +1; a++){
            for (int b = -1; b <= +1; b++){
                for (int c = -1; c <= +1; c++){
                        /* __ guess the coordinate of the neighbor __ */
                        ss[0] = cs[0]+a;
                        ss[1] = cs[1]+b;
                        ss[2] = cs[2]+c;

                    // TODO for perf cond add check
                    /* __ if coordinate out of range : no neighbor __ */
                    //  if ((ss[0] >= 0 && ss[0] < mpiDomains[0]) && (ss[1] >= 0 && ss[1] < mpiDomains[1]) && (ss[2] >= 0 && ss[2] < mpiDomains[2])){
                    //  For a process group with cartesian structure, the function MPI_CART_RANK translates the logical process coordinates to process ranks as they are used by the point-to-point routines.
                    //  For dimension i with periods(i) = true, if the coordinate, coords(i), is out of range, that is, coords(i) < 0 or coords(i)  dims(i), it is shifted back to the interval 0 <= coords(i) < dims(i) automatically.
                    //
                    //  Out-of-range coordinates are erroneous for non-periodic dimensions. Versions of MPICH before 1.2.2 returned MPI_PROC_NULL for the rank in this case.
                    
                    MPI_Cart_rank(com, ss, &wr);
                    //-1 -> left, bottom, back
                    // 0 -> same position
                    //+1 -> right, top, front
                    this->neighbors2Send[(1+c)+3*((1+b)+3*(1+a))] = wr;
                    this->neighbors2Recv[(1-c)+3*((1-b)+3*(1-a))] = wr;
                    
                    }
                }
            }


        double boxSizePerDomain[3];
        for (int n=0; n<3; n++){
        
            string res = GET + dirs[n] + "resolution";
            string right_str= GET + dirs[n] + "right";
            this->boxCoordinates[n][0] = 0.0; //origin point (0,0,0) is left lower coner i=0 j=0 k=0
            this->boxCoordinates[n][1] = callPyFloatFunction( pInstance, right_str, BRACKETS );
            
            
            this->boxSizes[n] = boxCoordinates[n][1];
            //  length of the box in normalized units
            
            //  length of the domain in normalized units
            boxSizePerDomain[n] = boxSizes[n]/mpiDomains[n];
            
            this->boxCoordinates[n][0] = boxCoordinates[n][0] + cs[n]*boxSizePerDomain[n];
            this->boxCoordinates[n][1] = boxCoordinates[n][0] + boxSizePerDomain[n];

            //  number of pixel per box
            this->totPixelsPerBoxSide[n] = (int) callPyLongFunction( pInstance, res, BRACKETS);
            
            if( boxSizes[n] <= 0.0 ){
                throw runtime_error("box has zero length!");
            }

            if (totPixelsPerBoxSide[n] == 1){
                // length per pixel equals to total size
                this->spatialSteps[n]=boxSizes[n];
            }else{
                this->spatialSteps[n]=boxSizes[n]/(double)(totPixelsPerBoxSide[n]);
            }
           
            this->resolution[n] = round(totPixelsPerBoxSide[n]/mpiDomains[n]);
            // number of pixel per domain
            // todo for problem rations
            this->offsetInPixels[n] = resolution[n]*cs[n];
            
        
    }
    
    this->dim                = SIMULATION_SIZE;
    this->timeStep           = callPyFloatFunction( pInstance, GET_TIMESTEP, BRACKETS );
    
    this->ppc                = (int) callPyLongFunction( pInstance, GET_NUM_OF_PARTICLES, BRACKETS);
    
    this->maxTimestepsNum    = callPyFloatFunction( pInstance, GET_MAX_TIMESTEPS_NUM, BRACKETS );
    
    this->timestepsNum2Write = callPyFloatFunction( pInstance, GET_TIMESTEP_WRITE, BRACKETS );
    
    this->outputDir          = callPyStringFunction( pInstance, GET_OUTPUT_DIR, BRACKETS );
    
    this->fileNameTemplate   = callPyStringFunction( pInstance, GET_FILENAME_TEMPLATE, BRACKETS );
    
    this->sourcePosition     = callPyFloatFunction( pInstance, GET_SOURCE_POSITION, BRACKETS );
    
    this->detectorPosition   = callPyFloatFunction( pInstance, GET_DETECTOR_POSITION, BRACKETS );
    
    this->flowDirection      = (int) callPyLongFunction( pInstance, GET_FLOW_DIRECTION, BRACKETS);
    
    
    
    for (int n=0; n<2; n++){
        
        string detsize = GET_DETECTOR_SIZE + dirs[n];
        this->detectorSize[n] =  callPyFloatFunction( pInstance, detsize, BRACKETS );
    
    }
    
    
    if(rank == 0){


        printf("[Loader]  Box size: \n");

        for (int n=0; n<3; n++){
            printf("[Loader]  [%s] [%1.5f, %1.5f] l = %1.5f res = %3i => step  = %1.6f\n",
                   dirs[n].c_str(), boxCoordinates[n][0],boxCoordinates[n][1],
                   boxSizes[n],resolution[n],spatialSteps[n]);
        }
        
        string flowDir = (flowDirection == 0 ? "X" : (flowDirection == 1 ? "Y" : "Z")  );
        msg = "[Loader] tested direction = "+flowDir;
        logger.writeMsg(msg.c_str(), INFO);
        
        msg = "[Loader] timeStep = "+to_string(timeStep);
        logger.writeMsg(msg.c_str(), INFO);
        msg = "[Loader] Max timeStep number = "+to_string(maxTimestepsNum);
        logger.writeMsg(msg.c_str(), INFO);
        msg = "[Loader]  timeStep number to write to file = "+to_string(timestepsNum2Write);
        logger.writeMsg(msg.c_str(), INFO);
        msg = "[Loader]  outputDir = "+outputDir;
        logger.writeMsg(msg.c_str(), INFO);
        msg = "[Loader]  filename template = "+fileNameTemplate;
        logger.writeMsg(msg.c_str(), INFO);
        msg = "[Loader]  ppc = "+to_string(ppc);
        logger.writeMsg(msg.c_str(), INFO);
        logger.writeMsg("[Loader] initialization...OK!", DEBUG);
    }

}

double Loader::getTimeStep(){
    return timeStep;
}

int Loader::getMaxTimestepsNum(){
    return maxTimestepsNum;
}

int Loader::getTimestepsNum2Write(){
    return timestepsNum2Write;
}

string Loader::getOutputDir() const {
    return outputDir;
}

string Loader::getFilenameTemplate(){
    return fileNameTemplate;
}

vector<double> Loader::getVelocity(double x, double y, double z){
    PyObject *pValue;
    vector<double> velocity;
    
    for (int n=0; n<3; n++){
        string varName = GET_VELOCITY+dirs[n];
        double val = callPyFloatFunctionWith3args( pInstance, varName, BRACKETS_3DOUBLE, x, y, z );
        velocity.push_back(val);
    }
    return velocity;
}


vector<double> Loader::getBfield(double x, double y, double z){
    PyObject *pValue;
    vector<double> bfield;
    
    for (int n=0; n<3; n++){
        string varName = GET_BFIELD+dirs[n];
        double val = callPyFloatFunctionWith3args( pInstance, varName, BRACKETS_3DOUBLE, x, y, z );
        bfield.push_back(val);
    }
    return bfield;
}

vector<double> Loader::getEfield(double x, double y, double z){
    PyObject *pValue;
    vector<double> efield;
    
    for (int n=0; n<3; n++){
        string varName = GET_EFIELD+dirs[n];
        double val = callPyFloatFunctionWith3args( pInstance, varName, BRACKETS_3DOUBLE, x, y, z );
        efield.push_back(val);
    }
    return efield;
}



double Loader::getMass(){
    return callPyFloatFunction( pInstance, GET_MASS, BRACKETS );
}

double Loader::getCharge(){
    return callPyFloatFunction( pInstance, GET_CHARGE, BRACKETS );
}


double Loader::getParticlesPerCellNumber(){
    return ppc;
}


Loader::~Loader(){
    Py_Finalize();
    logger.writeMsg("[Loader] FINALIZE...OK!", DEBUG);
}


/*#################################################################################################*/


PyObject * Loader::getPythonClassInstance(string className){
    PyObject  *pName, *pClass, *pModule, *pDict;
    string msg = "[Loader] Start to instantiate the Python class " + className;
    logger.writeMsg(msg.c_str(), DEBUG);
    
    #if PY_MAJOR_VERSION >= 3
    pName = PyUnicode_FromString(className.c_str());
    #else
    pName   = PyString_FromString(className.c_str());
    #endif
    
    pModule = PyImport_Import(pName);
    
    if( pModule == NULL ){
        logger.writeMsg("*****************************************************", CRITICAL);
        logger.writeMsg("****                                             ****", CRITICAL);
        logger.writeMsg("****  STOP SIMULATION!!!    INPUT FILE PROBLEM   ****", CRITICAL);
        logger.writeMsg("****                                             ****", CRITICAL);
        logger.writeMsg("****   try debug mode 'make -DLOG'               ****", CRITICAL);
        logger.writeMsg("****   check indents in python input file        ****", CRITICAL);
        logger.writeMsg("****   check python enviroment and imports       ****", CRITICAL);
        logger.writeMsg("*****************************************************", CRITICAL);
        exit(-1);
    }
    
    pDict = PyModule_GetDict(pModule);
    pClass = PyDict_GetItemString(pDict, className.c_str());
    
    if( PyCallable_Check(pClass) ){
        pInstance = PyObject_CallObject(pClass, NULL);
    }else{
        logger.writeMsg("[Loader] Cannot instantiate the Python class", CRITICAL);
        pInstance = nullptr;
    }
    
    logger.writeMsg("[Loader] finish to instantiate the Python class ", DEBUG);
    
    return pInstance;
}
    


double Loader::callPyFloatFunction( PyObject* instance,
                                   const string funcName,
                                   const string brackets){
    return PyFloat_AsDouble(getPyMethod(instance,funcName,brackets));
}

double Loader::callPyFloatFunctionWith3args( PyObject* instance,
                                            const string funcName,
                                            const string brackets,
                                            double x,double y,double z){
    
    PyObject *methodReturn = PyObject_CallMethod(instance, strdup(funcName.c_str()),
                                 strdup(BRACKETS_3DOUBLE.c_str()),x,y,z);
    
    if( methodReturn ){
        return PyFloat_AsDouble(methodReturn);
    }else{
        PyErr_Print();
        throw runtime_error("Problem to call python function: "+funcName);
    }
    
}

long Loader::callPyLongFunction( PyObject* instance,
                                const string funcName,
                                const string brackets){
    
    return PyLong_AsLong(getPyMethod(instance,funcName,brackets));
}

string Loader::callPyStringFunction( PyObject* instance,
                                    const string funcName,
                                    const string brackets){
    
    #if PY_MAJOR_VERSION >= 3
    PyObject* str = PyUnicode_AsUTF8String(getPyMethod(instance,funcName,brackets));
    return PyBytes_AsString(str);
    #else
    return PyString_AS_STRING(getPyMethod(instance,funcName,brackets));
    #endif

}

PyObject* Loader::getPyMethod(PyObject* instance,
                              const string funcName,
                              const string brackets){
    
    return PyObject_CallMethod(instance, strdup(funcName.c_str()),
                               strdup(brackets.c_str()));
}

/*#################################################################################################*/





