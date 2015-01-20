#include "Tensor.h"
#include <cstdlib>
#include <cmath>
#include <string>
#include <fstream>
#include "fonctionsAnnexes.h"

using namespace std;

// The object Tensor makes more easy the storage and use of the spins on the 3D periodic lattice.
// In particular, the spins are stored in the Tensor object in a 1D array, but Tensor allows the user to call a 3D array with the coordinates x, y and z of the selected spin, which is more userfriendly.

Tensor::Tensor() : nx(10), ny(10), nz(10)
{
    
}

Tensor::Tensor(int n1, int n2, int n3) : nx(n1), ny(n2), nz(n3)
{
    
}

void Tensor::init() // Random initialization
{
    /* Random generation of the initial state : +1 for a spin up, -1 for a spin down */

    spins.clear(); // We remove all previous elements

    for(int i = 0; i < nx; i++)
    {
        for(int j = 0; j < ny; j++)
        {
            for(int k = 0; k < nz; k++)
            {
                if(rand() % 2 < 1)
                {
                    spins.push_back(-1);
                }
                else
                {
                    spins.push_back(1);
                }
            }
        }
    }
}

void Tensor::read_config(double T)
{
    string temperature;
    temperature = to_string(T);
    string const read_file("results/configs/config_T=" + temperature + ".dat");
    
    ifstream read_stream(read_file.c_str());
    
    char state;
    
    for(int i = 0; i < nx; i++)
    {
        for(int j = 0; j < ny; j++)
        {
            for(int k = 0; k < nz; k++)
            {
                read_stream.get(state);
                if(state == '1')
                {
                    spins.push_back(1);
                }
                else
                {
                    spins.push_back(-1);
                }
            }
        }
    }
}

void Tensor::write_config(double T)
{
    string temperature;
    temperature = to_string(T);
    string const write_file("results/configs/config_T=" + temperature + ".dat");
    
    ofstream write_stream(write_file.c_str());
    
    for(int x = 0; x < nx; x++)
    {
        for(int y = 0; y < ny; y++)
        {
            for(int z = 0; z < nz; z++)
            {
                if(spins[x + nx * y + nx * ny * z] == 1)
                {
                    write_stream << 1;
                }
                else
                {
                    write_stream << 0;
                }
            }
        }
    }
}

double Tensor::getEnergy(double J0, double J1, double J2) // Computes and returns the energy of the system
{
    double E = 0;

    for(int x = 0; x < nx; x++)
    {
        for(int y = 0; y < ny; y++)
        {
            for(int z = 0; z < nz; z++)
            {
                /* There is a -= because there is a minus sign in the exchange Hamiltonian */
                E -= spins[x + nx * y + nx * ny * z] *
                    (J0 * (spins[mod(x + 1, nx) + nx * y + nx * ny * z] + spins[mod(x - 1, nx) + nx * y + nx * ny * z]
                    + spins[x + nx * mod(y + 1, ny) + nx * ny * z] + spins[x + nx * mod(y - 1, ny) + nx * ny * z])
                    + J1 * (spins[x + nx * y + nx * ny * mod(z + 1, nz)] + spins[x + nx * y + nx * ny * mod(z - 1, nz)])
                    + J2 * (spins[x + nx * y + nx * ny * mod(z + 2, nz)] + spins[x + nx * y + nx * ny * mod(z - 2, nz)]));
            }
        }
    }

    return E;
}

int Tensor::getValue(int x,int y,int z) // return the value of the spin at the spatial coordinates x, y and z
{
    return spins[x + nx * y + nx * ny * z];
}

double Tensor::getDeltaE(int x, int y, int z, double J0, double J1, double J2) // Computes and returns the energy variation of the system when a spin is switch (from up to down or down to up)
{
    /* Computation of the generated DeltaE
    There is a *2 because Einitial = -Efinal, thus DeltaE = 2 * Efinal
    There is an other *2 because it is requiered to count the interaction of i with j, then of j with i (one counts 2 times each interaction)
    Thus in total there is a *4 in the formula */

    return 4 * spins[x + nx * y + nx * ny * z] *
        (J0 * (spins[mod(x + 1, nx) + nx * y + nx * ny * z] + spins[mod(x - 1, nx) + nx * y + nx * ny * z]
        + spins[x + nx * mod(y + 1, ny) + nx * ny * z] + spins[x + nx * mod(y - 1, ny) + nx * ny * z])
        + J1 * (spins[x + nx * y + nx * ny * mod(z + 1, nz)] + spins[x + nx * y + nx * ny * mod(z - 1, nz)])
        + J2 * (spins[x + nx * y + nx * ny * mod(z + 2, nz)] + spins[x + nx * y + nx * ny * mod(z - 2, nz)]));
}

int Tensor::getMagnetization() // Computes and returns the total magnetization of the system
{
    int magnetization = 0;

    for(int i = 0; i < nx; i++)
    {
        for(int j = 0; j < ny; j++)
        {
            for(int k = 0; k < nz; k++)
            {
                magnetization += spins[i + nx * j + nx * ny * k];
            }
        }
    }

    return magnetization;
}

int Tensor::getMagnetizationPlane(int z) // Computes and returns the value of the magnetisation of a plane of the system
{
    /* Used in our energy minimization study */
    
    int magnetizationPlane = 0;

    for(int i = 0; i < nx; i++)
    {
        for(int j = 0; j < ny; j++)
        {
            magnetizationPlane += spins[i + nx * j + nx * ny * z];
        }
    }

    return magnetizationPlane;
}

int Tensor::getSumMagnetizationPlanes() // Computes and returns the sum of the absolute values of the magnetizations of each plane of the system
{
    /* We use the sum of the plane's magnetizations as a criterium in the energy minimization study */

    int magnetization = 0;

    for(int k = 0; k < nz; k++)
    {
        int magnetizationPlane = 0;

        for(int i = 0; i < nx; i++)
        {
            for(int j = 0; j < ny; j++)
            {
                magnetizationPlane += spins[i + nx * j + nx * ny * k];
            }
        }

        magnetization += abs(magnetizationPlane); //Pour voir le caractère FM de chaque plan
    }

    return magnetization;
}

int Tensor::getAFMCriterium() // The AFM criterium consists in computing a sum where, when looking at a given spin, if the neighbouring spins (z axis) are antiparallel (AFM), then we count +1, and if they are parallel (FM) we count -1
{
    /* We use this criterium in the energy minimization study */
    
    int criterium = 0;

    for(int i = 0; i < nx; i++)
    {
        for(int j = 0; j < ny; j++)
        {
            for(int k = 0; k < nz; k++)
            {
                // If the z+1 and z-1 spins are opposed, the criterium gets +1, otherwise -1
                criterium += - spins[i + nx * j + nx * ny * mod(k - 1, nz)] * spins[i + nx * j + nx * ny * mod(k + 1, nz)];
            }
        }
    }
    return criterium;
}

void Tensor::switchValue(int x, int y, int z)
{
    spins[x + nx * y + nx * ny * z] *= -1; // We switch the sign of the spin
}

Tensor::~Tensor()
{
    
}
