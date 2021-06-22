
*****************************************************
****---------------------------------------------****
****-------0--A-------A-------PPPPPP-------------****
****-------O--A------A A------P----P-------------****
****-------O-A------A   A-----P----P-------------****
****-------0A------AAAAAAA----PPPPPP-------------****
****-------O-A----A-------A---P------------------****
****-------O--A--A---------A--P------------------****
****----Kinetic Algorithm for Proton-radiography—****
*****************************************************

# KAP
 Kinetic Algorithm for Proton-radiography. \

* electro-magnetic fields are set in the input file\

* protons are described in a kinetic way, \
  dynamics is solved using first order interpolation of em fields\

*code pretends to be parallel.\

stack: C++11, MPI, HDF5, python2-3.\
_______________________
#       HOWTO
_______________________
1. before 'make' need to set in makefile\
    HDF5_PATH= path to hdf5 lib\
    MPI_PATH= path to mpi lib\
    PYTHON_INC= path to python2.7-3 include\
    PYTHON_LIB= path to python2.7-3 lib\

2. for running default example from src/input/Initializer.py\
    mpirun -n 1 aka.exe \

3. normally need to set input python file\
    mpirun -n 1 kap.exe PATH/TO/PYTHON/INPUT/FILE\

4. before running need to create output folder and set in the python file\

5. for visualization use python notebook KAP\


___________________________________
#     MPI and HDF5 installation
___________________________________

download openmpi v 4.0.5 from download.open-mpi.org/release/open-mpi/v4.0/openmpi-4.0.5.tar.gz

extract in folder /TEMP/FLD4INSTALLATION/, go to the folder

run in terminal:
1. ./configure --prefix=/FLD2LIBS/openmpi CFLAGS=-m64 CXXFLAGS=-m64 FFLAGS=-m64 FCFLAGS=-m64 \
// need to specify 64-bit compilation via additional flag that should be passed to ALL compilers
2. make
3. make install

download hdf5 v1.10.5 from http://hdfgroup.org/package/hdf5-1-10-5-tar-gz/?wpdmdl=13571

run in terminal:
1. ./configure --prefix=/FLD2LIBS/hdf5 --enable-parallel CC=/FLD2LIBS/openmpi/bin/mpicc
2. make
3. make install







