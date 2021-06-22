
#include "Misc.hpp"
using namespace std;

bool areSame(double a, double b)
{
    return fabs(a - b) < EPSILON_COMPARATOR;
}



void crossProd(double A[3], double B[3], double cross[3]){
    cross[0] = A[1] * B[2] - A[2] * B[1];
    cross[1] = A[2] * B[0] - A[0] * B[2];
    cross[2] = A[0] * B[1] - A[1] * B[0];
}




void rotate(int dir, double phi, double* vec){
    
    vector<vector<double>> rotMatrx;

    switch (dir) {
        case 0:
            rotMatrx = {{1 ,0       ,0        },
                        {0 ,cos(phi),-sin(phi)},
                        {0 ,sin(phi), cos(phi)}};break;
        case 1:
            rotMatrx = {{cos(phi), 0, sin(phi)},
                        {0       , 1, 0        },
                        {-sin(phi), 0, cos(phi)}};break;
        case 2:
            rotMatrx = {{cos(phi),-sin(phi), 0},
                        {sin(phi), cos(phi), 0},
                        { 0      , 0       , 1}};break;
        default: return;
    }
    
    double res[3];
    for(int i=0; i < 3; i++){
        res[i] = 0.0;
        for(int j=0; j < 3; j++){
            res[i] += rotMatrx[i][j]*vec[i];
        }
    }
    
    for(int i=0; i < 3; i++){
        vec[i] = res[i];
    }
}


