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
    /* Choice between FM and AFM study */
    int FM_AFM_choice = 0;

    while(FM_AFM_choice != 1 && FM_AFM_choice != 2) // While the chosen number is not 1 nor 2
    {
        cout << "Please choose between a FM or AFM study" << endl << "1.AFM" << endl << "2.FM" << endl;
        cin >> FM_AFM_choice;
    }

    cout << "Wait while the simulation is running..." << endl;

    /* Déclaration des variables */
    int nx = 4;
    int ny = 4;
    int nz = 8;

    double T = 1; //Temperature
    double J0 = 1;
    double J1 = 1; // constante d'échange entre plus proches voisins
    double J2init = - 0.1; // constante d'échange entre seconds voisins
    double J2final = - 1;
    double J2step = 0.1;
    
    int sweep = 0;
    int sweep_max = 10000; // Nombre de sweep de simulation
    int x = 0;
    int y = 0;
    int z = 0;
    
    double nombreEntre0Et1 = 0;
    double floattemp = 0;
    double temperatureStep = 0.1;
    double Tinit = 15;
    double Tinf = 1;
    int counterT = 0; // see far below where this counter is used

    /* Variables for output*/
    /* Energy */
    double E = 0;
    double DeltaE = 0; // Delta d'énergie entre 2 configurations
    double sum_E = 0;
    double sum_E_square = 0;
    double Cv = 0;
    /* Magnetization */
    double M = 0;
    int DeltaM = 0; // Delta d'aimantation entre 2 configurations
    double sum_M = 0;
    double sum_M_square = 0;
    double chi = 0;
    /* AFM criteria */
    double sumMagnetizationPlanes = 0;
    int DeltaSumMagnetizationPlanes = 0;
    double sum_sumMagnetizationPlanes = 0;
    double AFMcriterium = 0;
    int DeltaAFMcriterium = 0;
    double sum_AFMcriterium = 0;
    /* Configuration */
    vector<double> configuration(nz);
    vector<double> sum_configuration(nz);

    /* Tensor of the states */
    Tensor States(nx, ny, nz);

    

    /* NOUVELLE SECTION : BOUCLE SUR J2/J1 */
    double estimated_time_minutes = (J2init - J2final) * (Tinit - Tinf) * sweep_max/(1000000 * 100 * temperatureStep * J2step);
    cout << "Estimated time : " << estimated_time_minutes << " minutes, or " << estimated_time_minutes / 60 << " hours" << endl;

    for(double J2 = J2init; J2 >= J2final; J2 -= J2step)
    {
        counterT = 0; // see far below where this counter is used
        string label = "_J2=" + to_string(J2);

        /* Initialisation de rand */
        srand(time(NULL));

        /* Initialisation de l'état du système */
        States.init();

        /* Déclaration des flux pour écrire dans les fichier */
        bool streams;

        string const plot_file("results/results" + label);
        string const Cv_file("results/data/Cv" + label + ".dat");
        string const chi_file("results/data/chi" + label + ".dat");
        string const sumMagnetizationPlanes_file("results/data/sumMagnetizationPlanes" + label + ".dat");
        string const AFMcriterium_file("results/data/AFMcriterium" + label + ".dat");
        string const magnetization_file("results/data/magnetization" + label + ".dat");
        string const configuration_file("results/data/configuration" + label + ".dat");

        ofstream plot_stream(plot_file.c_str());
        ofstream Cv_stream(Cv_file.c_str());
        ofstream chi_stream(chi_file.c_str());
        ofstream sumMagnetizationPlanes_stream(sumMagnetizationPlanes_file.c_str());
        ofstream AFMcriterium_stream(AFMcriterium_file.c_str());
        ofstream magnetization_stream(magnetization_file.c_str());
        ofstream configuration_stream(configuration_file.c_str());

        if(FM_AFM_choice == 1) //AFM
        {
            /* We select the useful streams (the streams have to be declared outside the 'if'*/
            magnetization_stream.close();

            streams = sumMagnetizationPlanes_stream && AFMcriterium_stream && Cv_stream && chi_stream && plot_stream && configuration_stream;
        }
        else //FM (we assume FM_AFM_choice ==2)
        {
            sumMagnetizationPlanes_stream.close();
            AFMcriterium_stream.close();

            streams = magnetization_stream && Cv_stream && chi_stream && plot_stream && configuration_stream;
        }

        if(streams)
        {
            /* Entêtes */
            plot_stream << "#Parameters : J0 = J1 = " << J1 << "; J2 = " << J2 << "; number of sweeps = " << sweep_max << endl;
            Cv_stream << "#Temperature Cv" << endl;
            chi_stream << "#Temperature chi" << endl;
            configuration_stream << "#T z sumMagnetization" << endl;

            if(FM_AFM_choice == 1) //AFM
            {
                sumMagnetizationPlanes_stream << "#Temperature sumMagnetizationPlanes" << endl;
                AFMcriterium_stream << "#Temperature AFMcriterium" << endl;
            }
            else //FM (we assume FM_AFM_choice ==2)
            {
                magnetization_stream << "#Temperature magnetization" << endl;
            }

            /* Algorithme Monte Carlo */

            /* REMARQUE : FONCTION SEPAREE A FAIRE, CLAIREMENT ! */
            for(T = Tinit; T >= Tinf; T -= temperatureStep)
            {
                cout << "J2 : " << J2 << "; T : " << T << endl;
                counterT ++;

                /* We reinitialize the output variables */
                E = 0;
                sum_E = 0;
                sum_E_square = 0;
                DeltaE = 0;

                M = 0;
                sum_M = 0;
                sum_M_square = 0;
                DeltaM = 0;

                Cv = 0;
                chi = 0;

                sumMagnetizationPlanes = 0;
                sum_sumMagnetizationPlanes = 0;
                DeltaSumMagnetizationPlanes = 0;
                AFMcriterium = 0;
                sum_AFMcriterium = 0;
                DeltaAFMcriterium = 0;

                for(int i = 0; i < nz; i++)
                {
                    configuration[i] = 0;
                    sum_configuration[i] = 0;
                }
                
                /********************************************/
                /********** Algorithme en lui-même **********/
                /********************************************/
                
                for(sweep = 0; sweep < sweep_max; sweep++)
                {                    
                    for(x = 0; x < nx; x++)
                    {
                        for(y = 0; y < ny; y++)
                        {
                            for(z = 0; z < nz ; z++)
                            {
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
                                    } // else nothing happens
                                }
                            } //end for z
                        } // end for y
                    } // end for x
                        
                    if(sweep >= sweep_max/2) // We keep the final half of the states to compute the averaged variables
                    {
                        E = States.getEnergy(J0, J1, J2);
                        sum_E += E;
                        sum_E_square += E * E;

                        M = States.getMagnetization();
                        sum_M += M;
                        sum_M_square += M * M;

                        if (FM_AFM_choice == 1) //AFM
                        {
                            sumMagnetizationPlanes = States.getSumMagnetizationPlanes();
                            sum_sumMagnetizationPlanes += sumMagnetizationPlanes;

                            AFMcriterium = States.getAFMCriterium();
                            sum_AFMcriterium += AFMcriterium;
                        }

                        for(int i = 0; i < nz; i++)
                        {
                            configuration[i] = States.getMagnetizationPlane(i);
                            sum_configuration[i] += configuration[i];
                        }
                    }

                }//end for sweep_max

                /* Writing results */

                /* Cv = (1/kT^2) * (<E^2> - <E>^2) */
                Cv = (2/(sweep_max*T*T)) * (sum_E_square - sum_E * sum_E * 2 / sweep_max);

                /* chi = (1/kT) * (<M^2> - <M>^2) */
                chi = (2/(sweep_max*T)) * (sum_M_square - sum_M * sum_M * 2 / sweep_max);

                Cv_stream << T << " " << Cv << endl;
                chi_stream << T << " " << chi << endl;

                if(FM_AFM_choice == 1) //AFM
                {
                    sumMagnetizationPlanes_stream << T << " " << 2 * sum_sumMagnetizationPlanes / sweep_max << endl;
                    AFMcriterium_stream << T << " " << 2 * sum_AFMcriterium / sweep_max << endl;
                }
                else //FM (we assume FM_AFM_choice ==2)
                {
                    magnetization_stream << T << " " << abs(2 * sum_M / sweep_max) << endl;
                }

                /* We only plot for a few values of T so that the graph is not too ugly */
                if(counterT % 2 == 0)
                {
                    for(int i = 0; i < nz; i++)
                    {
                        configuration_stream << T << " " << i << " " << abs(2 * sum_configuration[i] / sweep_max) << endl;
                    }
                    configuration_stream << endl;
                }

            }//end for T

            if(FM_AFM_choice == 1) //AFM
            {
                plot_stream << "set xlabel 'Temperature'" << endl << "set ylabel 'sum of the Magnetization of the Planes'" << endl 
                    << "plot 'data/sumMagnetizationPlanes" + label + ".dat'" << endl 
                    << "set ylabel 'AFM criterium'" << endl 
                    << "plot 'data/AFMcriterium" + label + ".dat'" << endl << "pause -1" << endl;

                sumMagnetizationPlanes_stream.close();
                AFMcriterium_stream.close();
            }
            else //FM (we assume FM_AFM_choice == 2)
            {
                plot_stream << "set ylabel 'Total Magnetization'" << endl 
                    << "plot 'data/magnetization" + label + ".dat'" << endl << "pause -1" << endl;

                magnetization_stream.close();
            }

            plot_stream << "set ylabel 'Cv, chi'" << endl 
                << "plot 'data/Cv" + label + ".dat', 'data/chi" + label + ".dat'" << endl << "pause -1" << endl 
                << "set ylabel 'spin position along z'" << endl << "set zlabel 'magnetization of the plane at height z'" << endl
                << "splot [" << Tinf << ":" << Tinit << "] [0:" << nz - 1 << "] [0:" << nx*ny << "] 'data/configuration" + label + ".dat' with lines"
                << endl << "pause -1" << endl;

            Cv_stream.close();
            chi_stream.close();
            configuration_stream.close();
            plot_stream.close();

//            system("start C:\\Users\\Thibault\\Desktop\\Projet_OPA\\gnuplot\\bin\\gnuplot.exe");

        }//end if(streams)
        else
        {
            cout << "Error while opening files" << endl;
        }
    }

    return 0;
}
