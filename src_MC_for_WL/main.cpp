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
    int nz = 40;

    double J0 = 1;
    double J1 = 1; // nearest neighbour interaction constant
    double J2_min = - 1;
    double J2_max = 0.4;
    double J2_step = 0.1;
    
    int nombrePas = 5000 * nx * ny * nz; // Number of Monte-Carlo steps
    // it is an order of magnitude, 1000 steps per spin seem good, but more might be needed.
    int x = 0;
    int y = 0;
    int z = 0;
    
    double nombreEntre0Et1 = 0;
    double floattemp = 0;
    double temperatureStep = 0.1;
    double Tinit = 11;
    double Tinf = 0.1;

    /* Variables for output*/
    /* Energy */
    double E = 0;
    double DeltaE = 0; // Delta d'énergie entre 2 configurations
    double sum_E = 0;
    double sum_E_square = 0;
    double Cv = 0;
    double sigmaE = 0;

    /* Tensor of the states */
    Tensor States(nx, ny, nz);

    /* Duration estimation */
//    double estimated_time_minutes = (J2init - J2final) * (Tinit - Tinf) * nombrePas / (35000000 *  temperatureStep * J2step);
//    cout << "Estimated time for normal run : " << estimated_time_minutes << " minutes, or " << estimated_time_minutes / 60 << " hours" << endl;

    /* Initialisation de rand */
    srand(time(NULL));
    
    string label;
    
    for(double J2 = J2_min; J2 <= J2_max; J2 += J2_step)
    {
        /* Initialisation de l'état du système */
        States.init();

        /* Déclaration des flux pour écrire dans les fichier */
        bool streams;
        
        label = "_J2 = " + to_string(J2) + "_" + to_string(nx) + "x" + to_string(ny) + "x" +to_string(nz);
        // /!\ WANRING: if you modify this label, you also have to change the other label at the end of the code, when the WL_config file is created.
        
        string const plot_file("results/results" + label);
        string const Cv_file("results/data/Cv" + label + ".dat");
        string const energy_file("results/data/energy" + label + ".dat");
        string const sigmaE_file("results/data/sigmaE" + label + ".dat");

        ofstream plot_stream(plot_file.c_str());
        ofstream Cv_stream(Cv_file.c_str());
        ofstream energy_stream(energy_file.c_str());
        ofstream sigmaE_stream(sigmaE_file.c_str());

        streams = plot_stream && Cv_stream && energy_stream && sigmaE_stream;

        if(streams)
        {
            /* Headers */
            plot_stream << "#Parameters : nx = " << nx << "; ny = " << ny << "; nz = " << nz <<
            "; J0 = " << J0 << "; J1 = " << J1 << "; J2 = " << J2 << "; number of steps = " << nombrePas << endl;
            Cv_stream << "#Temperature Cv" << endl;
            energy_stream << "#Temperature Mean_Energy" << endl;
            sigmaE_stream << "#Temperature sigmaE^2" << endl;

            /* Monte-Carlo temperature loop */

            for(double T = Tinit; T >= Tinf; T -= temperatureStep)
            {
                cout << "J2: " << J2 << "; T: " << T << endl;

                /* We reinitialize the output variables */
                E = 0;
                sum_E = 0;
                sum_E_square = 0;
                DeltaE = 0;

                Cv = 0;
                
                /*******************************************/
                /********** Monte-Carlo algorithm **********/
                /*******************************************/
                
                for(int step = 0; step < nombrePas; step++)
                {
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

                /* Writing results */
                
                /* <E> = sum_E / number_of_spins */
                energy_stream << T << " " << sum_E * 2 / nombrePas << endl;
                
                /* Cv = (1/kT^2) * (<E^2> - <E>^2) */
                Cv = (2/(nombrePas*T*T)) * (sum_E_square - sum_E * sum_E * 2 / nombrePas);
                Cv_stream << T << " " << Cv << endl;
                
                /* sigmaE^2 = (<E^2> - <E>^2) */
                sigmaE = ((sum_E_square - sum_E * sum_E * 2 / nombrePas) * 2 / nombrePas);
                sigmaE_stream << T << " " << sigmaE << endl;
                
                /* Writing configurations */
                States.write_config(T, J2);
                
            }//end for T

            plot_stream << "set xlabel 'T'" << endl
                        << "set ylabel 'Cv'" << endl
                        << "plot 'data/Cv" << label << ".dat'" << endl
                        << "pause -1" << endl
                        
                        << "set ylabel '<E>'" << endl
                        << "plot 'data/energy" << label << ".dat'" << endl
                        << "pause -1" << endl
                        
                        << "set ylabel 'sigmaE^2'" << endl
                        << "plot 'data/sigmaE" << label << ".dat'" << endl
                        << "pause -1" << endl;
            
            Cv_stream.close();
            energy_stream.close();
            sigmaE_stream.close();
            plot_stream.close();
            
        }//end if(streams)
        else
        {
            cout << "Error while opening files" << endl;
        }
    } // end for J2



    /*******************************************/
    /*********** Config file for WL ************/
    /*******************************************/

    double Emin = 0;
    double Emax = 0;

    int numberPoints = 0; // the number of points/lines in the energy file over which we want to compute the average min and max energies
    int totalPoints = 0; // the number of points/lines in the energy file
    int jumpPoints = 0; // number of lines we must jump in the middle of the file

    totalPoints = ((Tinit - Tinf) / temperatureStep) + 1;
    numberPoints = totalPoints * 0.05;
    jumpPoints = totalPoints - 2 * numberPoints;

    string const configWL_file("results/configWL_" + to_string(nx) + "x" + to_string(ny) + "x" +to_string(nz));

    ofstream configWL_stream(configWL_file.c_str());

    configWL_stream << "#J2 Emin Emax" << endl;

    if(configWL_stream)
    {
        for(double J2 = J2_min; J2 <= J2_max; J2 += J2_step)
        {
            Emin = 0;
            Emax = 0;

            label = "_J2 = " + to_string(J2) + "_" + to_string(nx) + "x" + to_string(ny) + "x" +to_string(nz);

            string const read_energy_file("results/data/energy" + label + ".dat");

            ifstream read_energy_stream(read_energy_file.c_str());

            string line; // allows to jump lines
            double T_buff; // stores T values temporarily
            double E_buff; // stores E values temporarily

            getline(read_energy_stream, line); // jump the first comment line

            //computation of the mean max Energy
            for(int i = 0; i < numberPoints ; i++)
            {
                read_energy_stream >> T_buff;
                read_energy_stream >> E_buff;

                Emax += E_buff;
            }

            Emax /= numberPoints;

            getline(read_energy_stream, line); // jump to the next line, because after the last >>, we were at the end of the line, and not at the beginning of the following line

            //jump some lines
            for(int i = 0 ; i < jumpPoints ; i++)
            {
                getline(read_energy_stream, line);
            }

            //computation of the mean min Energy
            for(int i = 0; i < numberPoints ; i++)
            {
                read_energy_stream >> T_buff;
                read_energy_stream >> E_buff;

                Emin += E_buff;
            }

            Emin /= numberPoints;

            //writting in configWL file
            configWL_stream << J2 << " " << Emin << " " << Emax << endl;

        } // end for J2

    } // end if(configWL_stream)

    else
    {
        cout << "Error while opening configWL file !" << endl;
    }

    return 0;
}
