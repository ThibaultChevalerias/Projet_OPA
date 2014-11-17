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
    int visits[number_bins]; // Histogram of visits of each bin
    for(int i = 0; i < number_bins; i++)
    {
        visits[i] = 0; // Initialization
    }
    double entropy[number_bins]; // Entropy for each energy
    for(int i = 0; i < number_bins; i++)
    {
        entropy[i] = 0; // Initialization
    }

    int currentBin = 0; // in order to know in which bin our current energy is
    int proposedBin = 0; // idem for the proposed energy
    
    /* Initialization of the system */
    
    Tensor System(nx, ny, nz);
    
    System.init();
    
    /* Declaration of streams (data output) */

    string const gE_file("g(E).dat");

    ofstream gE_stream(gE_file.c_str());

    /* Initialization of rand */
    
    srand(time(NULL));
    
    if(gE_stream)
    {
        /* =============================================================== */
        /* ==================== Wang-Landau algorithm ==================== */
        /* =============================================================== */
    
        cout << "Wait while the simulation is running..." << endl;
    
        while(lnf > epsilon)
        {
            cout << "lnf = " << lnf << endl;

            currentEnergy = System.getEnergy(J0, J1, J2);
            currentBin = locateBin(deltaE, currentEnergy);
            
            step = 0;
            
            while(step < (step_max / sqrt(lnf)))
            {
                /* We propose a new state */
                xChosen = rand() % nx;
                yChosen = rand() % ny;
                zChosen = rand() % nz;
                    
                proposedEnergy = currentEnergy + System.getDeltaE(xChosen, yChosen, zChosen, J0, J1, J2);
                proposedBin = locateBin(deltaE, proposedEnergy);
                
                /* Acceptance of the new state */
                random_double = (rand() % random_range + 1) / random_range;
            
                if(log(random_double) < (entropy[currentBin] - entropy[proposedBin])
                {
                    // if accepted, update the energy and the system:
                    currentEnergy = proposedEnergy;
                    System.switchValue(xChosen, yChosen, zChosen);
                }
                // if rejected, nothing changes (same energy, same system)
                
                visits[currentEnergy] ++;
                entropy[currentEnergy] += lnf;
                
                if(step % 10000) // A modifier ??????????????????????????????????
                {
                    if(isFlat(flatness_limit, &visits[0], sizeof(visits)/sizeof(visits[0])))
                    {
                        break;
                    }
                }
                
                step++;
                
            } // end while(flatness < flatness_limit && step < step_max)
            
            for(int i = 0; i < number_bins; i++)
            {
                visits[i] = 0; // We reinitialize the histogram of visits
            }
    
            lnf /= 2; // We reduce f until f < epsilon

        } // end  while(lnf > epsilon)

        /* ======================================================== */
        /* ==================== Writing output ==================== */
        /* ======================================================== */

        gE_stream << "#E_bin_min g(E)" << endl;

        for(int i = 0; i < number_bins; i++)
        {
            gE_stream << E_min + i * deltaE << " " << entropy[i] << endl
        }

        gE_stream.close()

    }//end if(gE_stream)
    else
    {
        cout << "error while opening g(E).dat" << endl;
    }

}
