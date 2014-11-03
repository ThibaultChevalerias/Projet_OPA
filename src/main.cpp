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
    int nx = 12;
    int ny = 12;
    int nz = 12;

    double T = 1; //Temperature
    double J0 = 1;
    double J1 = 1; // constante d'échange entre plus proches voisins
    double J2init = - 0.3; // constante d'échange entre seconds voisinsq
    double J2final = - 0.5;
    double J2step = 1;
    int nombrePas = 200000; // Nombre de pas de simulation
    int xChosen = 0; //Voir section Monte Carlo
    int yChosen = 0; //Voir section Monte Carlo
    int zChosen = 0; //Voir section Monte Carlo
    double nombreEntre0Et1 = 0;
    double floattemp = 0;
    double temperatureStep = 0.1;
    double Tinit = 15;
    double Tinf = 1;
    int compteurT = 0;

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
    double estimated_time_minutes = (J2init - J2final) * (Tinit - Tinf) * nombrePas/(1000000 * 100 * temperatureStep * J2step);
    cout << "Estimated time : " << estimated_time_minutes << " minutes, or " << estimated_time_minutes / 60 << " hours" << endl;

    for(double J2 = J2init; J2 > J2final; J2 -= J2step)
    {
        compteurT = 0;
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
            plot_stream << "#Parameters : J0 = J1 = " << J1 << "; J2 = " << J2 << "; number of steps = " << nombrePas << endl;
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
            for(T = Tinit; T > Tinf; T -= temperatureStep)
            {
                cout << "J2 : " << J2 << "; T : " << T << endl;
                compteurT ++;

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

                for(int pas = 0; pas < nombrePas; pas++)
                {
                    /* Algorithme en lui-même */

                    xChosen = rand() % nx;
                    yChosen = rand() % ny;
                    zChosen = rand() % nz;

                    /* Calcul du DeltaE engendré */
                    DeltaE = States.getDeltaE(xChosen, yChosen, zChosen, J0, J1, J2);

                    /* DeltaM */
                    DeltaM = - 2 * States.getValue(xChosen, yChosen, zChosen); // Change of +2 or -2 in magnetization

                    if (FM_AFM_choice == 1) //AFM
                    {
                        /* DeltaSumMagnetizationPlanes */
                        DeltaSumMagnetizationPlanes = abs(States.getMagnetizationPlane(zChosen) - 2 * States.getValue(xChosen, yChosen, zChosen))
                            - abs(States.getMagnetizationPlane(zChosen));

                        /* For DeltaAFMcriterium, only the zChosen + 1 and zChosen - 1 spins are affected because the spins around zChosen are the same*/
                        DeltaAFMcriterium = 2 * States.getValue(xChosen, yChosen, zChosen) * (States.getValue(xChosen, yChosen, mod(zChosen + 2, nz))
                            + States.getValue(xChosen, yChosen, mod(zChosen - 2, nz)));
                    }
                    
                    /* Decision taken */

                    if(DeltaE <= 0) // Si l'énergie finale est plus faible
                    {
                        States.switchValue(xChosen, yChosen, zChosen); // On multiplie par -1 : on change l'état
                    }
                    else // Si l'énergie finale est plus élevée
                    {
                        floattemp = rand() % 1000 + 1; // CHANGER LA PRECISION ???
                        nombreEntre0Et1 = floattemp / 1000;
                        if(log(nombreEntre0Et1) < - DeltaE/T)
                        {
                            States.switchValue(xChosen, yChosen, zChosen); // On multiplie par -1 : on change l'état
                        }
                        else
                        {
                            DeltaE = 0; // So that the computation of E leaves it unchanged.
                            DeltaM = 0; // idem for the calculation of M
                            DeltaSumMagnetizationPlanes = 0;
                            DeltaAFMcriterium = 0;                            
                        }
                    }

                    if(pas == nombrePas/2)
                    {
                        E = States.getEnergy(J0, J1, J2);
                        sum_E = E;
                        sum_E_square = E * E;

                        M = States.getMagnetization();
                        sum_M = M;
                        sum_M_square = M * M;

                        if (FM_AFM_choice == 1) //AFM
                        {
                            sumMagnetizationPlanes = States.getSumMagnetizationPlanes();
                            sum_sumMagnetizationPlanes = sumMagnetizationPlanes;

                            AFMcriterium = States.getAFMCriterium();
                            sum_AFMcriterium = AFMcriterium;
                        }

                        for(int i = 0; i < nz; i++)
                        {
                            configuration[i] = States.getMagnetizationPlane(i);
                            sum_configuration[i] = configuration[i];
                        }
                    }

                    if(pas > nombrePas/2) // We keep the final half of the states to compute the averaged variables
                    {
                        E += DeltaE;
                        sum_E += E;
                        sum_E_square += E * E;

                        M += DeltaM;
                        sum_M += M;
                        sum_M_square += M * M;

                        if (FM_AFM_choice == 1) //AFM
                        {
                            sumMagnetizationPlanes += DeltaSumMagnetizationPlanes;
                            sum_sumMagnetizationPlanes += sumMagnetizationPlanes;

                            AFMcriterium += DeltaAFMcriterium;
                            sum_AFMcriterium += AFMcriterium;
                        }

                        configuration[zChosen] += DeltaM;
                        for(int i = 0; i < nz; i++)
                        {
                            sum_configuration[i] += configuration[i];
                        }
                        
                    }
                }//end for nombrePas

                /* Writing results */

                /* Cv = (1/kT^2) * (<E^2> - <E>^2) */
                Cv = (2/(nombrePas*T*T)) * (sum_E_square - sum_E * sum_E * 2 / nombrePas);

                /* chi = (1/kT) * (<M^2> - <M>^2) */
                chi = (2/(nombrePas*T)) * (sum_M_square - sum_M * sum_M * 2 / nombrePas);

                Cv_stream << T << " " << Cv << endl;
                chi_stream << T << " " << chi << endl;

                if(FM_AFM_choice == 1) //AFM
                {
                    sumMagnetizationPlanes_stream << T << " " << 2 * sum_sumMagnetizationPlanes / nombrePas << endl;
                    AFMcriterium_stream << T << " " << 2 * sum_AFMcriterium / nombrePas << endl;
                }
                else //FM (we assume FM_AFM_choice ==2)
                {
                    magnetization_stream << T << " " << abs(2 * sum_M / nombrePas) << endl;
                }

                /* We only plot for a few values of T so that the graph is not too ugly */
                if(compteurT % 5 == 0)
                {
                    for(int i = 0; i < nz; i++)
                    {
                        configuration_stream << T << " " << i << " " << 2 * sum_configuration[i] / nombrePas << endl;
                    }
                    configuration_stream << endl;
                }

            }//end for T

            if(FM_AFM_choice == 1) //AFM
            {
                plot_stream << "plot 'data/sumMagnetizationPlanes" + label + ".dat', 'data/AFMcriterium" + label + ".dat'" << endl << "pause -1" << endl;

                sumMagnetizationPlanes_stream.close();
                AFMcriterium_stream.close();
            }
            else //FM (we assume FM_AFM_choice ==2)
            {
                plot_stream << "plot 'data/magnetization" + label + ".dat'" << endl << "pause -1" << endl;

                magnetization_stream.close();
            }

            plot_stream << "plot 'data/Cv" + label + ".dat', 'data/chi" + label + ".dat'" << endl << "pause -1" << endl
                << "splot [0:" << nz - 1 << "] [-" << nz*nz << ":" << nz*nz << "] 'data/configuration" + label + ".dat' with lines"
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
