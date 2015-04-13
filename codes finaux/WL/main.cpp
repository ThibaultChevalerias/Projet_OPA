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

/* code giving the value of the density of energy states g(E) in function of the energy E of the system. It writes ln(g(E)) in a configuration file that is used by the main_computation code. */
/* Warning : please take notice of the presence of the logarithme : "lng(E)" is ln(g(E)). It explains the exponnential in the main_computation code. */
/* The presence of this logarithm may vary in the litterature. Same with lnf : it is ln(f), and may be only f in the litterature. Keep in mind that a division by 2 is, in that case, equivalent to a square root if the logarithm is removed */
/* Thus, this code might be slightly different from some litterature because of this choice of logarithm */

int main()
{
    /* ================================================================== */
    /* ==================== Declaration of variables ==================== */
    /* ================================================================== */
    
    /* Size of the system (spin box) */
    int nx = 10; // Size following x
    int ny = 10; // Size following y
    int nz = 56; // Size following z
    
    /* Exchange constants : spin interaction. Determins if spins interactions are ferromagnetic (positiv value) or antiferromagnetic (negativ value) and their respectiv intensity */
    double J0 = 1; // plane (x,y) ferromagnetic here
    double J1 = 1; // nearest neighbours following the z axis (ferro here)
    double J2 = - 0.9; // next nearest neighbours following the z axis (may be ferro or antiferro)
    
    /* Selection of the starting configuration file */
    double T_source = 6.5; // This code doesn't take into account the temperature. However, it needs a starting spin configuration given by the Monte Carlo algorithm.
                           // This double will be used to select the name of the configuration file.
                           // We consider that a spin configuration at a "medium" temperature (not too high, not too low) is a good starting point.
    
    /* Wang-Landau parameters */
    double lnf = 1; // This parameter will be devided by 2 each time the "flatness" criterium (see below) is reached. The Algorithm stops when lnf < epsilon
    double epsilon = 1.0/100000000; // the smaller epsilon is, the more precise the algorithm is, but it also requiers more time to run. Typical value: 1/1e8, must not be bigger or results will be wrong.
    
    double flatness_limit = 0.80; // Criterium determining wether the algorithm has converged or not. This criterium is over the histogram "visits"
                                  // It checks if the histogram "visits" is flat (see fonctionsAnnexes.cpp for more details)
                                  // When we reach this limit, we consider the histogram "visits" as flat, and change the value of f (division by 2), and clean the histogram "visits" (we fill it with 0).
                                  // This criterium must be between 0 and 1.
                                  // The more this criterium is close to 1, the more precise the algorithme is, but it will requier more time to run.
    
    int step = 0;
    double step_max = 10000000000; // Maximum number of steps allowed for the flatness to get above flatness_limit. If this number of steps is reached, we do as if the flatness_limit was reached, but we write a log error to keep track of the issue
    int missed_steps = 0; // see algorithm
    double limit_step = 0; // see algorithm; needs to be a double because it can get too big for an int
    int counter_billions = 0; // see algorithm; counts the billions we remove because the step's number is too high
    
    /* Variables to choose a random spin */
    int x = 0; // coordinate x of the spin
    int y = 0; // coordinate y of the spin
    int z = 0; // coordinate z of the spin
    
    /* Variables randomly generated */
    int random_int = 0;
    double random_double = 0;
    int random_range_int = 100000000;
    double random_range = 100000000; // float precision (1e8)
    
    /************************************/
    /* Getting the min and max energies */
    /************************************/

    // We read the minimum and maximum energy in the energy file from the Monte Carlo algorithm. The goal will be to define the "energy bins" for the WL algorithm, and the extrema over wich we don't want to go.

    /* This string allows us to read the file written by the MC algorithm with the good naming convention */
    string const configWL_file("results/configWL_" + to_string(nx) + "x" + to_string(ny) + "x" + to_string(nz));

    ifstream configWL_stream(configWL_file.c_str()); // reading stream
    
    string line; // will be used to jump some lines and headers in the energy file from the Monte Carlo
    getline(configWL_stream, line); // We jump the first line (it's a comment/header)
    
    /* Some variable initialisation. Their use is explain in the following lines. */
    double J2_read = 0;
    double E_max = 0;
    double E_min = 0;
    
    /* Get the minimum and maximum energies of the system for the choosen J2 parameter */
    while(configWL_stream) // we keep reading the config file until a proper J2 is found
    {
        configWL_stream >> J2_read; // We read the J2 of the current line
        if(J2_read == J2) // if this J2 matches the J2 we're currently using, we take the values of E_min, E_max, and we break the while
        {
            configWL_stream >> E_min;
            configWL_stream >> E_max;
            break;
        }
        else // else, we go to the next line
        {
            getline(configWL_stream, line);
        }
    }
    
    if(E_max == E_min) // if no proper J2 has been found in the config file, we print an error on the screen
    {
        cout << "ERROR WHILE READING configWL FILE !!!" << endl;
    }
    
    /****************************/
    /* Defining the energy bins */
    /****************************/

    int const number_bins = 60000; // We cut the energy interval in number_bins bins
                                   // This parameter must be changed if the size of the system is changed, for there must be more bins if the system is bigger (more intermediary states).
                                   // We haven't found any easy way to determin this parameter except by "trying until it works".
                                   // The goal is to have a different bin for each possible energy that the system may reach (so that, each time its energy changes, the bin is different)
                                   // The easiest way to determin that there are enough bins is when we see bins filled with 0 appear (energy ranged never reached because smaller than the quantums of energy of the system)
                                   // In order to check that, just open the lng(E) file and check if you see a lot of 0 in the middle of the file. It means the bin (energy range) has never been reached.
                                   // Warning : the requiered number of bin also change with the value of the parameters J0, J1 and J2. The smaller they are, the bigger the number of bins.
    
    /**********************************************************************************/
    /* Translation/Rescaling of the zero of energy in order to have positive energies */
    /**********************************************************************************/

    double deltaE = (E_max - E_min) / number_bins; // Energy range of each bin
    double rescale = abs(E_min) + deltaE; // Will be used for the translation of the 0. We add deltaE to be sur and avoid aving the first bin begining at 0.

    if(E_min < 0)
    {
        E_max += rescale; // we add abs(E_min) + deltaE
        E_min += rescale; // we add abs(E_min) + deltaE
        /* This rescale is in order to have strictly (non-zero) positive energies, 
        because negative and zero energies triggers the divergence of Cv in Wang_Landau_computation program (because of an exponential) */
    }
    
    /***************************************/
    /* Some more declarations of variables */
    /***************************************/

    /* keeping track of the energy of the system */
    double currentEnergy = 0; // energy of the current state of the system
    double proposedEnergy = 0; // energy of the system after inversion of a selected spin
    
    /* Histogram and entropy for the Wang-Landau algorithm */
    vector<int> visits (number_bins,0); // Histogram of visits of each energy bin, initialized with zeros
    vector<double> entropy(number_bins,0); // Entropy for each energy bin, initialized with zeros
    
    /* keeping track of the bin corresponding to the energy of the system */
    int currentBin = 0; // in order to know in which bin our current energy is
    int proposedBin = 0; // idem for the proposed energy
    
    /********************************/
    /* Initialization of the system */
    /********************************/

    Tensor System(nx, ny, nz);
    
    System.read_config(T_source, J2); // We have chosen T=6.5 for now
    
    /****************************************/
    /* Declaration of streams (data output) */
    /****************************************/

    /* These strings allow to open files with the correct naming convention */
    string const lngE_file("results/lng(E)_J2=" + to_string(J2) + "_" + to_string(nx) + "x" + to_string(ny) + "x" + to_string(nz) + ".dat");
    string const error_file("results/errors_J2=" + to_string(J2) + "_" + to_string(nx) + "x" + to_string(ny) + "x" + to_string(nz) + ".dat");
    string const log_file("results/log_J2=" + to_string(J2) + "_" + to_string(nx) + "x" + to_string(ny) + "x" + to_string(nz) + ".dat");
    // We put all the output sent to the screen to the 'log' file in order to be able to 're-watch' what happened when the algorithm was running.
    ofstream lngE_stream(lngE_file.c_str());
    ofstream error_stream(error_file.c_str());
    ofstream log_stream(log_file.c_str());
    
    /* Initialization of rand */
    
    srand(time(NULL));
    
    if(lngE_stream && error_stream && log_stream) // Check that all streams are opened without errors
    {
        /* =============================================================== */
        /* ==================== Wang-Landau algorithm ==================== */
        /* =============================================================== */
        
        cout << "Wait while the simulation is running..." << endl << endl;
        log_stream << "Wait while the simulation is running..." << endl << endl;
        
        currentEnergy = System.getEnergy(J0, J1, J2) + rescale; // initialisation of the (computation of the) energy
        // For the rescale, see declaration of variables, we have translated the zero of energy
        currentBin = locateBin(E_min, deltaE, currentEnergy); // idem for the current bin
        
        while(lnf >= epsilon) // As said upward, the algorithm stops when lnf < epsilon
        {
            cout << "lnf = " << lnf << endl << endl;
            log_stream << "lnf = " << lnf << endl << endl; // keeping track of the value of lnf and if the algorithm converges for each value of lnf
            
            step = 0; // reinitialisation
            missed_steps = 0; // reinitialisation
            
            counter_billions = 0; // reinitialisation
            limit_step = step_max / sqrt(lnf); // indeed, when lnf gets smaller, we need more steps
            
            while(step <= limit_step) // As said upward, we stop a "lnf sweep" if there has been too many steps without the convergence (flatness limit) to be reached
            {
                /* We propose a new state (position of the spin selected for a possible switch (up to down or down to up)) */
                x = rand() % nx;
                y = rand() % ny;
                z = rand() % nz;
                
                /* Energy of the system if the spin is switched */
                proposedEnergy = currentEnergy + System.getDeltaE(x, y, z, J0, J1, J2); // no need to rescale again by "-E_min_initial", "currentEnergy" has already been rescaled
                
                if(proposedEnergy >= E_min && proposedEnergy <= E_max) // if we stay in the right energy range...
                {
                    proposedBin = locateBin(E_min, deltaE, proposedEnergy); // Energy bin corresponding to the proposed energy (the proposed energy is within the range of the bin)
                    
                    /* Acceptance of the new state */
                    random_int = (rand() % random_range_int + 1); // /!\ the algebraic operations with rand() have to be performed with integers.
                    random_double = random_int / random_range; // conversion of the integer into a double
                    if(log(random_double) < (entropy[currentBin] - entropy[proposedBin])) // acceptance condition
                    {
                        // if accepted, we update the energy, the current bin, and the system:
                        currentEnergy = proposedEnergy; // no need to rescale again by "-E_min_initial", "currentEnergy" has already been rescaled
                        currentBin = proposedBin;
                        System.switchValue(x, y, z); // we switch the spin
                        visits[proposedBin] ++; // add 1 (bin visited 1 more time)
                        entropy[proposedBin] += lnf; // add lnf
                    }
                    else
                    {
                        visits[currentBin] ++; // we stay at the same energy and bin (currentBin not changed)... which is thus "visited" one more time
                        entropy[currentBin] += lnf;
                    }
                }
                else // if we leave the "right" energy range studied by the algorithm (smaller than Emin or bigger than Emax)
                {
                    missed_steps ++; // nothing happens and we count the number of missed steps. At the end, this number should not be too big compared to the total number of steps
                }
                
                step++;

                if(step % 100000 == 0) // We check the flatness criterium every *** steps. Can be modified if flatness criterium very long or quick to reach.
                {
                    if(isFlat(flatness_limit, visits))
                    {
                        cout << 1000 * counter_billions + step / 1000000 << " million" << endl;
                        cout << "missed steps: " << missed_steps << endl << endl;
                        cout << "**************************************************" << endl;
                        cout << "********** The histogram is flat enough **********" << endl;
                        cout << "**************************************************" << endl << endl;
                        log_stream << 1000 * counter_billions + step / 1000000 << " million" << endl;
                        log_stream << "missed steps: " << missed_steps << endl << endl;
                        log_stream << "**************************************************" << endl;
                        log_stream << "********** The histogram is flat enough **********" << endl;
                        log_stream << "**************************************************" << endl << endl;
                        error_stream << lnf << ": no problem, flatness achieved"; // if the flatness is achieved, we write it in the error_file
                        break;
                    }
                }
                
                if(step % 100000000 == 0) // we show only a few steps in order to follow the program easily and not have too many lines at the screen
                {
                    cout << 1000 * counter_billions + step / 1000000 << " million" << endl;
                    cout << "missed steps: " << missed_steps / 1000000.0 << " million" << endl;
                    log_stream << 1000 * counter_billions + step / 1000000 << " million" << endl;
                    log_stream << "missed steps: " << missed_steps / 1000000.0 << " million" << endl;
                }
                
                if(step % 1000000000 == 0) // above approximately 2 billion, the number of steps will be too big for this int...
                {
                    cout << endl << "rescaling: we remove 1 billion" << endl << "reminder: lnf = " << lnf << endl << endl;
                    log_stream << endl << "rescaling: we remove 1 billion" << endl << "reminder: lnf = " << lnf << endl << endl;
                    counter_billions ++;
                    step = 0; // thus we reinitialize the counter and keep track of the number of billions in another variable
                    limit_step -= 1000000000; // and we update the limit by removing the current number of steps
                }
                
            } // end while(step < step_max)
            
            for(int i = 0; i < number_bins; i++)
            {
                visits[i] = 0; // We reinitialize the histogram of visits
            }
            
            error_stream << "; " << lnf << " --> " << lnf/2 << endl; // if the flatness has not been achieved, there will be a line in the error file without "flatness achieved (see above)
            
            lnf /= 2; // We reduce f until f < epsilon
        
        } // end  while(lnf > epsilon)
        
        /* ======================================================== */
        /* ==================== Writing output ==================== */
        /* ======================================================== */
        
        lngE_stream << "#number of lines of data in this file :" << endl;
        lngE_stream << number_bins << endl;
        
        lngE_stream << "#Number of sites : " << endl;
        lngE_stream << nx * ny * nz << endl;
        
        lngE_stream << "#Rescale added to the energy (to be substracted) : " << endl;
        lngE_stream << rescale << endl;
        
        lngE_stream << "#E_bin_mean ln(g(E))" << endl;
        
        for (int i = 0; i < number_bins; i++)
        {
            lngE_stream << (E_min + (2 * i + 1) * deltaE / 2) - rescale << " " << entropy[i] << endl; // computation of the mean energy of each bin
                                                                                                      // entropy[i] is the value of lng(E) with E the mean energy of the bin i
        }
        
        lngE_stream.close();

    }//end if(lngE_stream)

    else // if there is an error while opening the streams/files
    {
        cout << "error while opening lng(E).dat" << endl;
        log_stream << "error while opening lng(E).dat" << endl;
    }

    return 0;

}
