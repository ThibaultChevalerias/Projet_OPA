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
    /* Code giving the curves of <E>, Cv, sigma(E) (standard deviation) in order to use them for determining the parameters of the Wang-Landau algorithm */
    cout << "Wait while the simulation is running..." << endl;

    /* Déclaration des variables */
    int nx = 10;
    int ny = 10;
    int nz = 56;

    double J0 = 1;
    double J1 = 1; // nearest neighbour interaction constant
    double J2 = - 1;
    
    int nombrePas = 5000 * nx * ny * nz; // Number of Monte-Carlo steps
    // it is an order of magnitude, 1000 steps per spin seem good, but more might be needed.
    int x = 0;
    int y = 0;
    int z = 0;
    
    double nombreEntre0Et1 = 0;
    double floattemp = 0;
    double T = 8;

    /* Variables for output*/
    /* Energy */
    double E = 0;
    double DeltaE = 0; // Delta d'énergie entre 2 configurations
    double sum_E = 0;
    double sum_E_square = 0;

    /* Tensor of the states */
    Tensor States(nx, ny, nz);

    /* Duration estimation */
//    double estimated_time_minutes = (J2init - J2final) * (Tinit - Tinf) * nombrePas / (35000000 *  temperatureStep * J2step);
//    cout << "Estimated time for normal run : " << estimated_time_minutes << " minutes, or " << estimated_time_minutes / 60 << " hours" << endl;

    /* Initialisation de rand */
    srand(time(NULL));
    
    string label;
    
        /* Initialisation de l'état du système */
        States.init();

        /* Déclaration des flux pour écrire dans les fichier */
        bool streams;
        
        label = "_J2 = " + to_string(J2) + "_" + to_string(nx) + "x" + to_string(ny) + "x" +to_string(nz);
        // /!\ WANRING: if you modify this label, you also have to change the other label at the end of the code, when the WL_config file is created.
        
        string const convergence_file("convergence" + label + ".dat");
        
        ofstream convergence_stream(convergence_file.c_str());
        
        streams = convergence_stream;

        if(streams)
        {
            /* Headers */
            convergence_stream << "#Step Energy" << endl;
            
            /* Monte-Carlo temperature loop */

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

                    /* Calcul du DeltaE engendré */
                    DeltaE = States.getDeltaE(x, y, z, J0, J1, J2);

                    /* Decision taken */

                    if(DeltaE <= 0) // Si l'énergie finale est plus faible
                    {
                        States.switchValue(x, y, z); // On multiplie par -1 : on change l'état
                    }
                    else // Si l'énergie finale est plus élevée
                    {
                        floattemp = rand() % 100000000 + 1; // float precision 1e8
                        nombreEntre0Et1 = floattemp / 100000000;
                        if(log(nombreEntre0Et1) < - DeltaE/T)
                        {
                            States.switchValue(x, y, z); // On multiplie par -1 : on change l'état
                        } 
                        else // nothing happens
                        {
                            DeltaE = 0;
                        }
                    }
                    
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
        else
        {
            cout << "Error while opening files" << endl;
        }

    return 0;
}
