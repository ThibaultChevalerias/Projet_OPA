#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>
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
    double epsilon = 1.0/8; // à modifier ?????
    
    double flatness_limit = 0.5; // When we reach this limit, we consider the histogram as flat, and change the value of f.
    
    int step = 0;
    int step_max = 1000000; // Maximum number of steps allowed for the flatness to pass above flatness_limit
    
    /* Variables to choose a random spin */
    int xChosen = 0;
    int yChosen = 0;
    int zChosen = 0;
    
    /* Variables randomly generated */
    int random_int = 0;
    double random_double = 0;
    int random_range_int = 100000000;
    double random_range = 100000000; // float precision (1e8)
    
    /* Energy bins */
    double E_max = 525; // E_max = number of spins * maximum value taken by (4 * J0 + 2 * J1 + 2 * J2). This is an upper boundary.
    double E_min = - 725;
    int const number_bins = 50; // We cut the energy interval in number_bins bins
    double deltaE = (E_max - E_min) / number_bins;
    
    double currentEnergy = 0;
    double proposedEnergy = 0;
    
    /* Histogram and entropy for the Wang-Landau algorithm */
    vector<int> visits (number_bins,0); // Histogram of visits of each bin, initialized with zeros
    vector<double> entropy(number_bins,0); // Entropy for each energy, initialized with zeros
    
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

            step = 0;
            
            while(step < (step_max / sqrt(lnf)))
            {
                currentEnergy = System.getEnergy(J0, J1, J2);
                currentBin = locateBin(E_min, deltaE, currentEnergy);

                /* We propose a new state */
                xChosen = rand() % nx;
                yChosen = rand() % ny;
                zChosen = rand() % nz;
                
                proposedEnergy = currentEnergy + System.getDeltaE(xChosen, yChosen, zChosen, J0, J1, J2);
                proposedBin = locateBin(E_min, deltaE, proposedEnergy);
                
                /* Acceptance of the new state */
                random_int = (rand() % random_range_int + 1); // /!\ the algebraic operations with rand() have to be performed with integers.
                random_double = random_int / random_range;
                if(log(random_double) < (entropy[currentBin] - entropy[proposedBin]))
                {
                    // if accepted, update the energy and the system:
                    currentEnergy = proposedEnergy;
                    System.switchValue(xChosen, yChosen, zChosen);
                    visits[proposedBin] ++;
                    entropy[proposedBin] += lnf;
                }
                else
                {
                    visits[currentBin] ++;
                    entropy[currentBin] += lnf;
                }
                
                step++;

                if(step % 10000 == 0) // A modifier ??????????????????????????????????
                {
                    cout << step << endl;
                    if(isFlat(flatness_limit, visits))
                    {
                        cout << "The histogram is flat enough" << endl;
                        break;
                    }
                } 
            
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
        
        gE_stream << "#E_bin_mean g(E)" << endl;
        
        for(int i = 0; i < number_bins; i++)
        {
            gE_stream << E_min + (2 * i + 1) * deltaE / 2 << " " << entropy[i] << endl; //computation of the mean energy of each bin
        }
        
        gE_stream.close();
    }//end if(gE_stream)
    else
    {
        cout << "error while opening g(E).dat" << endl;
    }

    system("PAUSE");
    return O;

}
