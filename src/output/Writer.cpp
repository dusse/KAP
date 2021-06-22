#include "Writer.hpp"

using namespace std;
using namespace chrono;

Writer::Writer(shared_ptr<Loader> load, shared_ptr<GridManager> gridMnr, shared_ptr<Pusher> pshr):
                loader(move(load)), gridMgr(move(gridMnr)), pusher(move(pshr)){
    logger.reset(new Logger());

    this->outputDir       = loader->getOutputDir();
    this->fileNamePattern = loader->getFilenameTemplate();
    
    logger->writeMsg("[Writer] initialize...OK", DEBUG);
}

void Writer::write(int fileNum){
    auto start_time = high_resolution_clock::now();
   
    logger->writeMsg("[Writer] writing...", DEBUG);
    int subdomainXSize = loader->resolution[0];
    int subdomainYSize = loader->resolution[1];
    int subdomainZSize = loader->resolution[2];
    
    int size[3] = {subdomainXSize, subdomainYSize, subdomainZSize};
    
    int totalNodeNum = subdomainXSize*subdomainYSize*subdomainZSize;
    int ijNode;
    vector<vector<VectorVar>> vectorVars = gridMgr->getVectorVariablesForAllNodes();
    
    vector<VectorVar> doseFar = gridMgr->getDetectorDoseVar();
    vector<int> detRes        = gridMgr->getDetectorResolution();
    
    MPI_Info info = MPI_INFO_NULL;
    
    hid_t access = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(access, MPI_COMM_WORLD, info);
    
    string fileName = outputDir + fileNamePattern + to_string(fileNum) + ".h5";
    hid_t fileID = H5Fcreate(fileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, access);
    
    const string groupname = "/vars";
    hid_t group   = H5Gcreate(fileID, groupname.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    
    hid_t dxpl_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(dxpl_id, H5FD_MPIO_COLLECTIVE);
    
    int detnodes = detRes[0]*detRes[1];
    int varDim = 1;//doseFar[0].getSize();
    int sizeDet[3] = {detRes[0], detRes[1], 1};
    for (int dir = 0; dir < varDim; dir++) {
        
        string varName = "dosefar_"+to_string(dir);
        double* var = new double[detnodes*sizeof(double)];
        
        for (ijNode = 0; ijNode < detnodes; ijNode++) {
            var[ijNode] = doseFar[ijNode].getValue()[dir];
        }
        
        writeDetectorVar(fileID, group, dxpl_id, varName, var, sizeDet);
        
        delete [] var;
    }
    

    
    H5Gclose(group);
    H5Fflush(fileID, H5F_SCOPE_GLOBAL);
    H5Pclose(dxpl_id);
    H5Fclose(fileID);
    
    // TODO check for X and Y directions
//    writeTrackedParticles(fileNum);
    
    auto end_time = high_resolution_clock::now();
    string msg ="[Writer] writing duration = "
                +to_string(duration_cast<milliseconds>(end_time - start_time).count())+" ms";
    logger->writeMsg(msg.c_str(), DEBUG);
}


void Writer::writeTrackedParticles(int fileNum){
    auto start_time = high_resolution_clock::now();
    
    logger->writeMsg("[Writer] writing tracked...", DEBUG);
    
    MPI_Info info = MPI_INFO_NULL;
    
    hid_t access = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(access, MPI_COMM_WORLD, info);
    
    string fileName = outputDir + fileNamePattern + "tracked_" + to_string(fileNum) + ".h5";
    hid_t fileID = H5Fcreate(fileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, access);
    
    const string groupname = "/vars";
    hid_t group   = H5Gcreate(fileID, groupname.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    
    hid_t dxpl_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(dxpl_id, H5FD_MPIO_COLLECTIVE);
    
    hsize_t offset[3];
    hsize_t nLoc[3];
    hsize_t nGlob[3];
    hid_t dset;
    
    int sizepartcl[3] = {NUMBER_TRACKED_PARTICLES, 3, 1};
    
    string varName = "tracked_protons";
    
    double* var = new double[NUMBER_TRACKED_PARTICLES*3*sizeof(double)];
    
    for (int dir = 0; dir < 3; dir++) {
        for (int i = 0; i < NUMBER_TRACKED_PARTICLES; i++) {
            var[3*i+dir] = pusher->trackedParticles[3*i+dir];
        }
     }
        
    for (int dir = 0; dir < 3; dir++) {
        offset[dir] = 0.0;
        nLoc[dir]   = sizepartcl[dir];
        nGlob[dir]  = sizepartcl[dir];
    }
        
    hid_t memspace  = H5Screate_simple(3, nLoc , NULL);
    hid_t filespace = H5Screate_simple(3, nGlob, NULL);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, nLoc, NULL);
    
    
    dset  = H5Dcreate(group, varName.c_str(), H5T_NATIVE_DOUBLE, filespace,
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        
    H5Dwrite(dset, H5T_NATIVE_DOUBLE, memspace, filespace, dxpl_id, var);
    H5Dclose(dset);
        
    H5Sclose(memspace);
    H5Sclose(filespace);
    
    H5Gclose(group);
    H5Fflush(fileID, H5F_SCOPE_GLOBAL);
    H5Pclose(dxpl_id);
    H5Fclose(fileID);
    
    delete [] var;
    
    auto end_time = high_resolution_clock::now();
    string msg ="[Writer] writing tracked duration = "
    +to_string(duration_cast<milliseconds>(end_time - start_time).count())+" ms";
    logger->writeMsg(msg.c_str(), DEBUG);


}


void Writer::writeParallel(hid_t fileID, hid_t group, hid_t dxpl_id,
                           string dsetName , const double* data, int sizes[3]){

    int rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank);
    hsize_t offset[3];
    hsize_t nLoc[3];
    hsize_t nGlob[3];

    for (int dir = 0; dir < 3; dir++) {
        offset[dir] = loader->offsetInPixels[dir];
        nLoc[dir]   = sizes[dir];
        nGlob[dir]  = loader->totPixelsPerBoxSide[dir];
    }

    hid_t memspace  = H5Screate_simple(3, nLoc , NULL);
    hid_t filespace = H5Screate_simple(3, nGlob, NULL);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, nLoc, NULL);

    hid_t dset;
    dset  = H5Dcreate(group, dsetName.c_str(), H5T_NATIVE_DOUBLE, filespace,
                      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    
    H5Dwrite(dset, H5T_NATIVE_DOUBLE, memspace, filespace, dxpl_id, data);
    H5Dclose(dset);

    H5Sclose(memspace);
    H5Sclose(filespace);
}

void Writer::writeDetectorVar(hid_t fileID, hid_t group, hid_t dxpl_id,
                           string dsetName , const double* data, int sizes[3]){
    
    int rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank);
    hsize_t offset[3];
    hsize_t nLoc[3];
    hsize_t nGlob[3];
    
    for (int dir = 0; dir < 3; dir++) {
        offset[dir] = 0.0;
        nLoc[dir]   = sizes[dir];
        nGlob[dir]  = sizes[dir];
    }
    
    hid_t memspace  = H5Screate_simple(3, nLoc , NULL);
    hid_t filespace = H5Screate_simple(3, nGlob, NULL);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, nLoc, NULL);
    
    hid_t dset;
    dset  = H5Dcreate(group, dsetName.c_str(), H5T_NATIVE_DOUBLE, filespace,
                      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    
    H5Dwrite(dset, H5T_NATIVE_DOUBLE, memspace, filespace, dxpl_id, data);
    H5Dclose(dset);
    
    H5Sclose(memspace);
    H5Sclose(filespace);
}







