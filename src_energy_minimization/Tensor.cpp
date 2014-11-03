#include "Tensor.h"
#include <cstdlib>
#include <cmath>
#include "fonctionsAnnexes.h"

using namespace std;

Tensor::Tensor() : nx(10), ny(10), nz(10)
{
    
}

Tensor::Tensor(int n1, int n2, int n3) : nx(n1), ny(n2), nz(n3)
{
    
}

void Tensor::init() // Random initialization
{
    /* Génération aléatoire de l'état inital : +1 pour un spin up, -1 pour un spin down*/

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

double Tensor::getEnergy(double J0, double J1, double J2)
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

int Tensor::getValue(int x,int y,int z)
{
    return spins[x + nx * y + nx * ny * z];
}

double Tensor::getDeltaE(int x, int y, int z, double J0, double J1, double J2)
{
    /* Calcul du DeltaE engendré
    Il y a un *2 car Einitiale = -Efinale, donc DeltaE = 2 * Efinale
    Il y a un autre *2 car il faut compter l'interaction de i avec j, puis de j avec i (on compte 2 fois chaque interaction)
    Donc au total il y a un *4 */

    return 4 * spins[x + nx * y + nx * ny * z] *
        (J0 * (spins[mod(x + 1, nx) + nx * y + nx * ny * z] + spins[mod(x - 1, nx) + nx * y + nx * ny * z]
        + spins[x + nx * mod(y + 1, ny) + nx * ny * z] + spins[x + nx * mod(y - 1, ny) + nx * ny * z])
        + J1 * (spins[x + nx * y + nx * ny * mod(z + 1, nz)] + spins[x + nx * y + nx * ny * mod(z - 1, nz)])
        + J2 * (spins[x + nx * y + nx * ny * mod(z + 2, nz)] + spins[x + nx * y + nx * ny * mod(z - 2, nz)]));
}

int Tensor::getMagnetization()
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

int Tensor::getMagnetizationPlane(int z)
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

int Tensor::getSumMagnetizationPlanes()
{
    /* We use the sum of the plane's magnetizations as a criterium in the energy minimization study */

    int magnetization = 0;

    for(int k = 0; k < nz; k++)
    {
        magnetization += abs(spins.getMagnetizationPlane(k)); //Pour voir le caractère FM de chaque plan
    }

    return magnetization;
}

int Tensor::getAFMCriterium()
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
