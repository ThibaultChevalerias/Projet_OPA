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
    
    cout << "Wait while the simulation is running..." << endl;

    /* Déclaration des variables */
    int nx = 4;
    int ny = 4;
    int nz = 8;

    double T = 1; //Temperature
    double J0 = 1;
    double J1 = 1; // constante d'échange entre plus proches voisins
    double J2 = - 0.1; // constante d'échange entre seconds voisins
    
    /* External applied field */
    double H = 0;
    double Hmin = -5;
    double Hmax = 5;
    double Hstep = 0.01;
    
    
    int sweep = 0;
    int sweep_max = 100; // Nombre de sweep de simulation
    int x = 0;
    int y = 0;
    int z = 0;
    
    double nombreEntre0Et1 = 0;
    double floattemp = 0;

    /* Variables for output*/
    double DeltaE = 0;
    /* Magnetization */
    double M = 0;
    int DeltaM = 0; // Delta d'aimantation entre 2 configurations
    double sum_M = 0;

    /* Tensor of the states */
    Tensor States(nx, ny, nz);

    /* Initialisation de rand */
    srand(time(NULL));

    /* Initialisation de l'état du système */
    States.init_all_up();

    /* Déclaration des flux pour écrire dans les fichier */
    string const magnetization_file("magnetization.dat");

    ofstream magnetization_stream(magnetization_file.c_str());

    if (magnetization_stream)
    {
        /* Entêtes */
        magnetization_stream << "#Parameters : J0 = J1 = " << J1 << "; J2 = " << J2 << "; number of sweeps for a given H = " << sweep_max << endl
        << "#external_field_H magnetization" << endl;

        /* Algorithme Monte Carlo */

        /* REMARQUE : FONCTION SEPAREE A FAIRE, CLAIREMENT ! */
        for(H = Hmax; H >= Hmin; H -= Hstep)
        {
            /* We reinitialize the output variables */
            M = 0;
            sum_M = 0;
            DeltaM = 0;
            
            /********************************************/
            /********** Algorithme en lui-même **********/
            /********************************************/
            
            for(sweep = 0; sweep < sweep_max; sweep++)
            {      
                /* One sweep */
                for(x = 0; x < nx; x++)
                {
                    for(y = 0; y < ny; y++)
                    {
                        for(z = 0; z < nz ; z++)
                        {
                            /* Calcul du DeltaE engendré */
                            DeltaE = States.getDeltaE_zeeman(x, y, z, J0, J1, J2, H);
                            
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
                                } // else nothing happens
                            }
                        } //end for z
                    } // end for y
                } // end for x
                /* End of the sweep */
                
                if(sweep >= sweep_max/2) // We keep the final half of the states to compute the averaged variables
                {
                    M = States.getMagnetization();
                    sum_M += M;
                }

            }//end for sweep_max

            /* Writing results */

            magnetization_stream << H << " " << 2 * sum_M / sweep_max << endl;

        }//end for H (max --> min)
        
        for(H = Hmin; H <= Hmax; H += Hstep)
        {
            /* We reinitialize the output variables */
            M = 0;
            sum_M = 0;
            DeltaM = 0;
            
            /********************************************/
            /********** Algorithme en lui-même **********/
            /********************************************/
            
            for(sweep = 0; sweep < sweep_max; sweep++)
            {      
                /* One sweep */
                for(x = 0; x < nx; x++)
                {
                    for(y = 0; y < ny; y++)
                    {
                        for(z = 0; z < nz ; z++)
                        {
                            /* Calcul du DeltaE engendré */
                            DeltaE = States.getDeltaE_zeeman(x, y, z, J0, J1, J2, H);
                            
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
                                } // else nothing happens
                            }
                        } //end for z
                    } // end for y
                } // end for x
                /* End of the sweep */
                
                if(sweep >= sweep_max/2) // We keep the final half of the states to compute the averaged variables
                {
                    M = States.getMagnetization();
                    sum_M += M;
                }

            }//end for sweep_max

            /* Writing results */

            magnetization_stream << H << " " << 2 * sum_M / sweep_max << endl;

        }//end for H (min --> max)
        
        magnetization_stream.close();

    }//end if(streams)
    else
    {
        cout << "Error while opening files" << endl;
    }

    return 0;
}
