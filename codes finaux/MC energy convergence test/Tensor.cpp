/***
This file is a C++ object file. It implements an object called Tensor which models the system studied.
This file contains several functions associated with this object, which are detailed below.
***/
#include "Tensor.h"
#include <cstdlib> // Functions srand, rand
#include <cmath> // Function abs
#include <ctime> // Function time
#include <string>
#include <fstream>
#include <iostream>
#include "fonctionsAnnexes.h"

using namespace std;

Tensor::Tensor() : nx(12), ny(12), nz(12) // Basic constructor, makes automatically a 12x12x12 system.
{
    
}

Tensor::Tensor(int n1, int n2, int n3) : nx(n1), ny(n2), nz(n3) // With this constructor, you can choose the size of your system.
{
    
}

void Tensor::init()
{
    /***
    This function randomly initialises the system : +1 for a spin up, -1 for a spin down.
    ***/
    
    srand(time(NULL)); // Initialisation of rand.
    
    spins.clear(); // We remove all previous elements.

    for(int k = 0; k < nz; k++)
    {
        for(int j = 0; j < ny; j++)
        {
            for(int i = 0; i < nx; i++)
            {
                if(rand() % 2 < 1) // We have randomly 1 or 0, giving a spin up or down.
                /* Remark: '''if(rand() % 2){spins.push_back(1);}else{spins.push_back(-1);}''' should work exactly the same (because 0 is false and 1 is true), but has not been tested */
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

void Tensor::read_config(double T, double J2)
{
    /***
    This function reads the configuration in a configuration file and applies it to the system.
    ***/
    
    spins.clear(); // We remove all previous elements.
    
    /* We put the temperature and the value of J2 entered in the function parameters in strings. */
    string temperature;
    temperature = to_string(T);
    string J;
    J = to_string(J2);
    
    /* These strings allow us to find the right configuration file, according to the naming convention. */
    string const read_file("results/configs_" + to_string(nx) + "x" + to_string(ny) + "x" +to_string(nz) + "/config_J2=" + J + "_T=" + temperature + ".dat");
    
    ifstream read_stream(read_file.c_str()); // We open the file in reading mode.
    
    if(read_stream)
    {
        char state;
        
        for(int k = 0; k < nz; k++)
        {
            for(int j = 0; j < ny; j++)
            {
                for(int i = 0; i < nx; i++)
                {
                    read_stream.get(state); // We read each number (as a character, not an int) from the file and put it in a char variable called 'state'.
                    
                    if(state == '1') // Then we put new spins in the system according to the configuration file.
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
    else // If there is a problem with the opening of the file (naming convention, etc.)
    {
        cout << "Failed to open file " << read_file << endl;
    }
}

void Tensor::write_config(double T, double J2)
{
    /***
    This function creates a configuration file in which the system configuration is stored.
    ***/
    
    /* We put the temperature and the value of J2 entered in the function parameters in strings. */
    string temperature;
    temperature = to_string(T);
    string J;
    J = to_string(J2);
    
    /* These strings allow us to create the right configuration file name, according to the naming convention. */
    string const write_file("results/configs_" + to_string(nx) + "x" + to_string(ny) + "x" +to_string(nz) + "/config_J2=" + J + "_T=" + temperature + ".dat");
    
    ofstream write_stream(write_file.c_str());
    
    if(write_stream)
    {
        for(int z = 0; z < nz; z++)
        {
            for(int y = 0; y < ny; y++)
            {
                for(int x = 0; x < nx; x++)
                {
                    if(spins[x + nx * y + nx * ny * z] == 1) // We then write a 1 for a spin up and a 0 for a spin down.
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
    else // If there is a problem with the opening of the file (for instance, the 'results' directory has not been created).
    {
        cout << "Failed to open file " << write_file << endl;
    }
}

void Tensor::switchValue(int x, int y, int z)
{
    /***
    This function switches the spin at the coordinates (x,y,z).
    ***/
    spins[x + nx * y + nx * ny * z] *= - 1; // We simply switch the sign of the spin (up (+1) becomes down (-1), and vice versa).
}

double Tensor::getEnergy(double J0, double J1, double J2)
{
    /***
    This function computes and returns the energy of the system
    ***/
    
    double E = 0;

    for(int z = 0; z < nz; z++)
    {
        for(int y = 0; y < ny; y++)
        {
            for(int x = 0; x < nx; x++) // For each spin, we sum its energy.
            {
                /***
                There is a -= because there is a minus sign in the exchange Hamiltonian
                (see the exchange hamiltonian to understand the expression of the energy).
                -----
                The 'mod' function allow us to get the right number while passing through the borders of the system
                For example, if we want the neighbours along x of a spin located at x = nx - 1 (border because the numbering begins at 0): 
                We will use the 'x - 1' spin, but the 'x + 1' doesn't exist (it would be in 'nx' position).
                So, we use the 'mod(x + 1, nx)' spin, which is simply always 'x + 1' except in our case, where it is 'x + 1 - nx' (see 'fonctionsAnnexes.cpp').
                Thus, it is in this case the spin '0', the first spin.
                This implements our choice of the periodic boundary conditions : the spins at the boundaries interact with those at the other side of the system.
                ***/
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
    /***
    This function returns the value of the spin at the spatial coordinates x, y and z.
    ***/
    return spins[x + nx * y + nx * ny * z];
}

double Tensor::getDeltaE(int x, int y, int z, double J0, double J1, double J2)
{
    /***
    This function computes and returns the energy variation of the system when a spin is switched (from up to down or down to up)
    -----
    There is a '2*' because Einitial = - Efinal, thus DeltaE = Efinal - Einitial = 2 * Efinal.
    There is an other '2*' because it is requiered to count the interaction of i with j, then of j with i (one counts each interaction twice).
    Finally, there is a '4*' in the formula
    -----
    For more information on the use of the 'mod' function, see the description of the 'getEnergy(double J0, double J1, double J2)' function.
    ***/

    return 4 * spins[x + nx * y + nx * ny * z] *
        (J0 * (spins[mod(x + 1, nx) + nx * y + nx * ny * z] + spins[mod(x - 1, nx) + nx * y + nx * ny * z]
        + spins[x + nx * mod(y + 1, ny) + nx * ny * z] + spins[x + nx * mod(y - 1, ny) + nx * ny * z])
        + J1 * (spins[x + nx * y + nx * ny * mod(z + 1, nz)] + spins[x + nx * y + nx * ny * mod(z - 1, nz)])
        + J2 * (spins[x + nx * y + nx * ny * mod(z + 2, nz)] + spins[x + nx * y + nx * ny * mod(z - 2, nz)]));
}

int Tensor::getMagnetization()
{
    /***
    This function computes and returns the total magnetization of the system.
    ***/
    
    int magnetization = 0;

    for(int k = 0; k < nz; k++)
    {
        for(int j = 0; j < ny; j++)
        {
            for(int i = 0; i < nx; i++) // We simply add the value of each spin.
            {
                magnetization += spins[i + nx * j + nx * ny * k];
            }
        }
    }

    return magnetization;
}

int Tensor::getMagnetizationPlane(int z)
{
    /***
    This function computes and returns the value of the magnetisation of a (x,y) plane in the system.
    It is not used anymore, but we used it to check if the spins in an (x,y) plane were aligned or in a paramagnetic state.
    ***/
    
    int magnetizationPlane = 0;

    for(int j = 0; j < ny; j++)
    {
        for(int i = 0; i < nx; i++) // We sum the value of all spins in the (x,y) plane at the position z.
        {
            magnetizationPlane += spins[i + nx * j + nx * ny * z];
        }
    }

    return magnetizationPlane;
}

int Tensor::getSumMagnetizationPlanes()
{
    /***
    This function computes and returns the sum of the absolute values of the magnetizations of each plane of the system.
    It is not used anymore but we used it to check if the spins in each (x,y) planes were aligned or in a paramagnetic state.
    ***/

    int magnetization = 0;

    for(int k = 0; k < nz; k++)
    {
        int magnetizationPlane = 0;

        for(int j = 0; j < ny; j++)
        {
            for(int i = 0; i < nx; i++)
            {
                magnetizationPlane += spins[i + nx * j + nx * ny * k];
            }
        }

        magnetization += abs(magnetizationPlane); // We add the absolute value of each 'magnetizationPlane' value, because we want to know if each plane is in a ferromagnetic state (it is the absolute value that is interesting, not the up or down state (the sign)).
    }

    return magnetization;
}

int Tensor::getAFMCriterium()
{
    /***
    This function conputes and returns an AFM state criterium. It is not used anymore.
    The AFM criterium consists in computing a sum where, when looking at a given spin, if the neighbouring spins (z axis) are antiparallel (AFM state), then we count +1, and if they are parallel (FM state) we count -1.
    ***/
    
    int criterium = 0;

    for(int k = 0; k < nz; k++)
    {
        for(int j = 0; j < ny; j++)
        {
            for(int i = 0; i < nx; i++)
            {
                // If the 'z + 1' and 'z - 1' spins are opposed, the criterium gets '+1', otherwise '-1'.
                criterium += - spins[i + nx * j + nx * ny * mod(k - 1, nz)] * spins[i + nx * j + nx * ny * mod(k + 1, nz)];
            }
        }
    }
    return criterium;
}

Tensor::~Tensor() // The basic destructor of our object.
{
    
}
