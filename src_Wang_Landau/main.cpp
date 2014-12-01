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
    
    /* Size of the spin box */
    int nx = 5;
    int ny = 5;
    int nz = 5;
    
    /* Exchange constants */
    double J0 = 1; // plane (x,y) ferromagnetic here
    double J1 = 1; // nearest neighbours following the z axis (ferro here)
    double J2 = -0.5; // next nearest neighbours following the z axis (may be ferro or antiferro)
    
    /* Wang-Landau parameters */
    double lnf = 1;
    double epsilon = 1.0/8; // à modifier ?????
    
    double flatness_limit = 0.5; // When we reach this limit, we consider the histogram as flat, and change the value of f.
    
    int step = 0;
    int step_max = 1000000; // Maximum number of steps allowed for the flatness to get above flatness_limit
    
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
    double E_max = 550; // E_max < number of spins * maximum value taken by (4 * J0 + 2 * J1 + 2 * J2). This is an upper boundary.
    double E_min = - 750; // E_max et E_min à modifier !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    int const number_bins = 50; // We cut the energy interval in number_bins bins
    double deltaE = (E_max - E_min) / number_bins; // Energy range of each bin
    double rescale = abs(E_min) + deltaE;

    /* Translation of the zero of energy in order to have positive energies */
    if(E_min < 0)
    {
        E_max += rescale; // we add abs(E_min) + deltaE
        E_min += rescale; // we add abs(E_min) + deltaE
        /* This rescale is in order to have strictly (non-zero) positive energies, 
        because negative and zero energies triggers the divergence of Cv in Wang_Landau_computation program */
    }
    
    double currentEnergy = 0; // energy of the current state of the system
    double proposedEnergy = 0; // energy of the system after inversion of a selected spin
    
    /* Histogram and entropy for the Wang-Landau algorithm */
    vector<int> visits (number_bins,0); // Histogram of visits of each energy bin, initialized with zeros
    vector<double> entropy(number_bins,0); // Entropy for each energy bin, initialized with zeros
    
    int currentBin = 0; // in order to know in which bin our current energy is
    int proposedBin = 0; // idem for the proposed energy
    
    /* Initialization of the system */
    
    Tensor System(nx, ny, nz);
    
    System.init();
    
    /* Declaration of streams (data output) */
    
    string const lngE_file("lng(E).dat");
    
    ofstream lngE_stream(lngE_file.c_str());
    
    /* Initialization of rand */
    
    srand(time(NULL));
    
    if(lngE_stream)
    {
        /* =============================================================== */
        /* ==================== Wang-Landau algorithm ==================== */
        /* =============================================================== */
        
        cout << "Wait while the simulation is running..." << endl;
        
        currentEnergy = System.getEnergy(J0, J1, J2) + rescale; // initialisation of the (computation of the) energy
        // For the rescale, see declaration of variables, we have translated the zero of energy
        currentBin = locateBin(E_min, deltaE, currentEnergy); // idem for the current bin

        while(lnf > epsilon)
        {
            cout << "lnf = " << lnf << endl;

            step = 0;
            
            while(step < (step_max / sqrt(lnf)))
            {
                /* We propose a new state */
                xChosen = rand() % nx;
                yChosen = rand() % ny;
                zChosen = rand() % nz;
                
                proposedEnergy = currentEnergy + System.getDeltaE(xChosen, yChosen, zChosen, J0, J1, J2); // no need to rescale again by "-E_min_initial", "currentEnergy" has already been rescaled
                proposedBin = locateBin(E_min, deltaE, proposedEnergy);
                
                /* Acceptance of the new state */
                random_int = (rand() % random_range_int + 1); // /!\ the algebraic operations with rand() have to be performed with integers.
                random_double = random_int / random_range;
                if(log(random_double) < (entropy[currentBin] - entropy[proposedBin]))
                {
                    // if accepted, update the energy, current bin,  and the system:
                    currentEnergy = proposedEnergy; // no need to rescale again by "-E_min_initial", "currentEnergy" has already been rescaled
                    currentBin = proposedBin;
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
            
            } // end while(step < step_max)
            
            for(int i = 0; i < number_bins; i++)
            {
                visits[i] = 0; // We reinitialize the histogram of visits
            }
            
            lnf /= 2; // We reduce f until f < epsilon
        
        } // end  while(lnf > epsilon)
        
        /* ======================================================== */
        /* ==================== Writing output ==================== */
        /* ======================================================== */
        
        lngE_stream << "#number of lines of data in this file :" << endl;
        lngE_stream << number_bins << endl;
        
        lngE_stream << "#Number of sites : " << endl;
        lngE_stream << nx * ny * nz << endl;
        
        lngE_stream << "#E_bin_mean ln(g(E))" << endl;
        
        for (int i = 0; i < number_bins; i++)
        {
            lngE_stream << (E_min + (2 * i + 1) * deltaE / 2) / (nx * ny * nz) << " " << entropy[i] << endl; //computation of the mean energy of each bin
            // We divide by nx * ny * nz, in order to have the energy per site
        }
        
        lngE_stream.close();
    }//end if(lngE_stream)
    else
    {
        cout << "error while opening lng(E).dat" << endl;
    }

    system("PAUSE");
    return 0;

}
