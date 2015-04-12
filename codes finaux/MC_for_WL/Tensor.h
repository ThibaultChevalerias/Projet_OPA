/***
This file is the header file associated with 'Tensor.cpp'. See the .cpp file for information about the functions.
***/
#ifndef TENSOR_H_INCLUDED
#define TENSOR_H_INCLUDED

#include <vector>

using namespace std;

class Tensor
{
    /* All the functions in this class are in the 'public:' section */
    public: 
        Tensor();
        Tensor(int n1, int n2, int n3);
        void init();
        void read_config(double T, double J2);
        void write_config(double T, double J2);
        void switchValue(int x, int y, int z);
        double getEnergy(double J0, double J1, double J2);
        int getValue(int x,int y,int z);
        double getDeltaE(int x, int y, int z, double J0, double J1, double J2);
        int getMagnetization();
        int getMagnetizationPlane(int z);
        int getAFMCriterium();
        int getSumMagnetizationPlanes();
        ~Tensor();
    
    /* All the attributes in this class are in the 'private:' section */
    private:
        int nx; // The number of spins along x
        int ny; // The number of spins along y
        int nz; // The number of spins along z
        vector <int> spins; // The table 'spins' contains the state of every spin in the system (1 for 'up' and -1 for 'down').
        /*** The spins are ordered as follows: first we put nx spins along x (for z = y = 0), then we go to y = 1 and continue ; when y = ny-1, we go to z = 1 and continue, etc. Basically, they are ordered first along x, then along y, then along z.
        In more detail:
        * spins[0], spins[1], ..., spins[nx-1] are the spins along x (from 0 to nx-1) for y = z = 0,
        * spins[nx], spins[nx+1], ..., spins[2*nx-1] are the spins along x for z = 0, y = 1,
        etc.,
        * spins[(ny-1)*nx], spins[(ny-1)*nx+1], ..., spins[ny*nx-1] are the spins along x for z = 0, y = ny-1,
        * spins[ny*nx], spins[ny*nx+1], ..., spins[(ny+1)*nx-1] are the spins along x for z = 1, y = 0,
        etc.,
        * spins[(nz*ny-1)*nx], spins[(nz*ny-1)*nx+1], ..., spins[nz*ny*nx-1] are the spins along x for z = nz-1, y = ny-1.
        ***/
};

#endif // TENSOR_H_INCLUDED
