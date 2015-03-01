#ifndef TENSOR_H_INCLUDED
#define TENSOR_H_INCLUDED

#include <vector>

using namespace std;

class Tensor
{
public:
    Tensor();
    Tensor(int n1, int n2, int n3);
    void init();
    void read_config(double T, double J2);
    void write_config(double T, double J2);
    double getEnergy(double J0, double J1, double J2);
    int getValue(int x,int y,int z);
    double getDeltaE(int x, int y, int z, double J0, double J1, double J2);
    int getMagnetization();
    int getMagnetizationPlane(int z);
    int getAFMCriterium();
    int getSumMagnetizationPlanes();
    void switchValue(int x, int y, int z);
    ~Tensor();

private:
    int nx;
    int ny;
    int nz;
    vector <int> spins;
};

#endif // TENSOR_H_INCLUDED
