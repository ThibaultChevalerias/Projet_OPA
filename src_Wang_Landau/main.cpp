#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <fstream>
#include <string>
#include "fonctionsAnnexes.h"
#include "Tensor.h"

using namespace std;

int main()
{
    /* ================================================================== */
    /* ==================== Declaration of variables ==================== */
    /* ================================================================== */

    int nx = 5;
    int ny = 5;
    int nz = 5;
    
    /* Exchange constants */
    double J0 = 1;
    double J1 = 1;
    double J2 = -1;
    
    /* Wang-Landau parameters */
    double lnf = 1;
    double epsilon = 1/64; // à modifier ?
    
    double flatness = 0; // flatness of the histogram
    double flatness_limit = 0.9; // When we reach this limit, we consider the histogram as flat, and change the value of f.
    
    int step = 0;
    int step_max = 1000000; // Maximum number of steps allowed for the flatness to pass above flatness_limit
    
    /* Variables to choose a random spin */
    int xChosen = 0;
    int yChosen = 0;
    int zChosen = 0;
    
    /* Variables randomly generated */
    double random_double = 0;
    int random_range = 100000000; // float precision (1E8)
    
    /* Energy bins */
    double E_max = nx * ny * nz * (4 * abs(J0) + 2 * abs(J1) + 2 * abs(J2)); // E_max = number of spins * maximum value taken by (4 * J0 + 2 * J1 + 2 * J2). This is an upper boundary.
    double E_min = - E_max;
    int number_bins = 100; // We cut the energy interval in 100 bins
    double deltaE = (E_max - E_min) / number_bins;
    
    double currentEnergy = 0;
    double proposedEnergy = 0;

    /* Histogram and entropy for the Wang-Landau algorithm */
    int DOS[number_bins]; // Density Of energy States
    for(int i = 0; i < number_bins; i++)
    {
        DOS[i] = 0; // Initialization
    }
    double entropy[number_bins]; // Entropy for each energy
    for(int i = 0; i < number_bins; i++)
    {
        entropy[i] = 0; // Initialization
    }

    /* Initialization of the tensor of states */
    
    Tensor States(nx, ny, nz);
    
    States.init();
    
    /* Declaration of streams (data output) */

    string const gE_file("g(E).dat");

    ofstream gE_stream(gE_file.c_str());

    if(gE_stream)
    {
        /* =============================================================== */
        /* ==================== Wang-Landau algorithm ==================== */
        /* =============================================================== */
    
        cout << "Wait while the simulation is running..." << endl;
    
        while(lnf > epsilon)
        {
            currentEnergy = States.getEnergy(J0, J1, J2);
    
            while(flatness < flatness_limit && step < step_max)
            {
            
                /* Changement d'état */
            
            
            
            } // end while(flatness < flatness_limit && step < step_max)
        
            lnf /= 2; // We reduce f until f < epsilon

        } // end  while(lnf > epsilon)



    }//end if(gE_stream)
    else
    {
        cout << "error while opening g(E).dat" << endl;
    }

}
