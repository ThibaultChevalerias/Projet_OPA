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
    /* Code giving the curve of the energy of the system in function of the monte carlo steps to check the convergence of the algorithm  */
    cout << "Wait while the simulation is running..." << endl;

    /*****************************/
    /* Déclaration des variables */
    /*****************************/

    /* Size of the system (spin box) */
    int nx = 10; // size following the x axis (number of spins)
    int ny = 10; // size following the y axis (number of spins)
    int nz = 56; // size following the z axis (number of spins)

    /* Exchange parameter : spin interaction. Determins if spins interactions are ferromagnetic (positiv value) or antiferromagnetic (negativ value) and their respectiv intensity */
    double J0 = 1; // nearest neighbours interaction in the (x,y) plan
    double J1 = 1; // nearest neighbour interaction following the z axis
    double J2 = - 1; // next nearest neighbours interaction following the z axis.
    
    int nombrePas = 5000 * nx * ny * nz; // Number of Monte-Carlo steps
    // it is an order of magnitude, 1000 steps per spin seem good, but more might be needed.
    
    /* Declaration of variables that will contain the position of a selected spin in the spin box. These variables will be used to select a spin for a switch of value (up to down or down to up) */
    int x = 0;
    int y = 0;
    int z = 0;
    
    /* declaration of variables that will be used for the Metropolis algorithm */
    double nombreEntre0Et1 = 0;  //number between O and 1
    double floattemp = 0;
    
    double T = 8; // Temperature for which we want to check the convergence of the algorithm

    /* Variables for output*/
    /* Energy */
    double E = 0; // energy of the system
    double DeltaE = 0; // difference of energie between 2 configurations
    double sum_E = 0; // utile ?
    double sum_E_square = 0; // utile ?

    /* Tensor of the states */
    Tensor States(nx, ny, nz);

    /* Duration estimation */
//    double estimated_time_minutes = (J2init - J2final) * (Tinit - Tinf) * nombrePas / (35000000 *  temperatureStep * J2step);
//    cout << "Estimated time for normal run : " << estimated_time_minutes << " minutes, or " << estimated_time_minutes / 60 << " hours" << endl;

    /* Initialisation of rand */
    srand(time(NULL));
    
    string label;
    
        /* Initialisation of the state of the system */
        States.init();

        /* Declaration of the streams to write in the files */
        bool streams;
        
        /* this string allows to follow the naming convention of the files */
        label = "_J2 = " + to_string(J2) + "_" + to_string(nx) + "x" + to_string(ny) + "x" +to_string(nz);
        // /!\ WANRING: if you modify this label, you also have to change the other label at the end of the code, when the WL_config file is created.
        
        string const convergence_file("convergence" + label + ".dat");
        
        ofstream convergence_stream(convergence_file.c_str());
        
        streams = convergence_stream;

        if(streams) // if no error occured while opening the files...
        {
            /* Headers */
            convergence_stream << "#Step Energy" << endl;

                cout << "J2: " << J2 << "; T: " << T << endl;

                /* We reinitialize the output variables */
                E = 0;
                sum_E = 0;
                sum_E_square = 0;
                DeltaE = 0;

                /*******************************************/
                /********** Monte-Carlo algorithm **********/
                /*******************************************/
                
                for(int step = 0; step < nombrePas; step++)
                {
                    if(step % 10000 == 0)
                    {
                        convergence_stream << step << " " << States.getEnergy(J0, J1, J2) << endl; // writing the energies at each step to check the convergence
                    }
                    
                    /* We choose a random spin */
                    x = rand() % nx;
                    y = rand() % ny;
                    z = rand() % nz;

                    /* Calcul of the generated energy difference DeltaE */
                    DeltaE = States.getDeltaE(x, y, z, J0, J1, J2);

                    /* Decision taken */

                    if(DeltaE <= 0) // If the energy of the system is smaller after the spin switch than before the spin switch
                    {
                        States.switchValue(x, y, z); // We multiply by -1 : we switch the state of the spin (up to down or down to up)
                    }
                    else // If the energy of the system is larger after the spin switch than before the spin switch
                    {
                        floattemp = rand() % 100000000 + 1; // float precision (1e8) of the random number
                        nombreEntre0Et1 = floattemp / 100000000;
                        if(log(nombreEntre0Et1) < - DeltaE/T)
                        {
                            States.switchValue(x, y, z); // We multiply by -1 : we switch the state of the spin
                        } 
                        else // nothing happens
                        {
                            DeltaE = 0;
                        }
                    }
                    
                    // utile ?
                    if(step == nombrePas/2)
                    {
                        E = States.getEnergy(J0, J1, J2);
                        sum_E = E;
                        sum_E_square = E * E;
                    }

                    if(step > nombrePas/2) // We keep the final half of the states to compute the averaged variables
                    {
                        E += DeltaE;
                        sum_E += E;
                        sum_E_square += E * E;
                    }

                }//end for nombrePas

            convergence_stream.close();
            
        }//end if(streams)
        else // if an error occured while opening the file
        {
            cout << "Error while opening files" << endl;
        }

    return 0;
}
