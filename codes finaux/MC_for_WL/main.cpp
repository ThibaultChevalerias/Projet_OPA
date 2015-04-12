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
    /* Code giving the curves of <E>, Cv, and sigma(E) (standard deviation) in order to use them for determining the parameters of the Wang-Landau algorithm */
    cout << "Wait while the simulation is running..." << endl;

    /********************************/
    /* Déclaration of the variables */
    /********************************/

    /* Size of the system (spin box) */
    int nx = 10; // size of the system following the x axis
    int ny = 10; // size of the system following the y axis
    int nz = 60; // size of the system following the z axis

    /* Exchange parameter : spin interaction. Determins if spins interactions are ferromagnetic (positiv value) or antiferromagnetic (negativ value) and their respectiv intensity */
    double J0 = 1; // nearest neighbours interaction in the (x,y) plan
    double J1 = 1; // nearest neighbour interaction following the z axis
    double J2_min = - 1; // next nearest neighbours interaction following the z axis. The code will perform a 'for' loop over this parameter, between J2_min and J2_max with a step J2_step
    double J2_max = 0;
    double J2_step = 0.1;
    
    int nombrePas = 5000 * nx * ny * nz; // Number of Monte-Carlo steps
    // it is an order of magnitude, 1000 steps per spin seem good, but more might be needed.

    /* Declaration of variables that will contain the position of a selected spin in the spin box. These variables will be used to select a spin for a switch of value (up to down or down to up) */
    int x = 0;
    int y = 0;
    int z = 0;
    
    /* declaration of variables that will be used for the Metropolis algorithm */
    double nombreEntre0Et1 = 0; //number between O and 1
    double floattemp = 0;

    /* Declaration of variables for a 'for' loop over the temperature. The loop goes from high temperatures to low temperatures */
    // Remark : this variable is without dimension. It is not a real temperature.
    double temperatureStep = 0.1;
    double Tinit = 11;
    double Tinf = 0.1;

    /* Variables for output*/
    /* Energy */
    double E = 0; // Energy
    double DeltaE = 0; // difference of energie between 2 configurations

    /* Thermodynamics variables to compute the specific heat Cv and the standard deviation of the energy sigmaE */
    double sum_E = 0; // sum of the energies of the system for each state
    double sum_E_square = 0; //sum of the squared energies of the system for each state
    double Cv = 0; // scpecific heat
    double sigmaE = 0; // standard deviation of the energy

    /* Tensor of the states (see the object Tensor in Tensor.cpp and Tensor.h) */
    Tensor States(nx, ny, nz);

    /* Duration estimation */
//    double estimated_time_minutes = (J2init - J2final) * (Tinit - Tinf) * nombrePas / (35000000 *  temperatureStep * J2step);
//    cout << "Estimated time for normal run : " << estimated_time_minutes << " minutes, or " << estimated_time_minutes / 60 << " hours" << endl;

    /* Initialisation of the rand function */
    srand(time(NULL));
    
    /* A variable to contain strings */
    string label;


    /***************************/
    /* Beginning of the 'code' */
    /***************************/

    for(double J2 = J2_min; J2 <= J2_max; J2 += J2_step) // loop over J2 to study the phase diagram for different values of -J2/J1 (J1 is fixed)
    {
        if(abs(J2 < 0.00001)) // sometimes, instead of J2=0, we get something e-18 (error on double addition), so we put zero
        {
            J2 = 0;
        }
        
        /* Initialisation of the state of the system */
        States.init();

        /* Declaration of streams to write in files */
        bool streams;
        
        /* This string allows us to write files names with the right convention */
        label = "_J2 = " + to_string(J2) + "_" + to_string(nx) + "x" + to_string(ny) + "x" +to_string(nz);
        // /!\ WANRING: if you modify this label, you also have to change the other label at the end of the code, when the WL_config file is created.
        
        /* Declaration of files */
        string const plot_file("results/results" + label); // This file will contain instructions for gnuplot to plot the graphs quickly
        string const Cv_file("results/data/Cv" + label + ".dat"); // This file will contain a list of values of Cv in function of the temperature, usable directly by Gnuplot
        string const energy_file("results/data/energy" + label + ".dat"); // This file will contain a list of values of the energy of the system in function of the temperature, usable directly by Gnuplot
        string const sigmaE_file("results/data/sigmaE" + label + ".dat"); // This file will contain a list of values of the standard deviation of the energy of the system in function of the temperature, usable directly by Gnuplot

        /* open streams to write in files */
        ofstream plot_stream(plot_file.c_str());
        ofstream Cv_stream(Cv_file.c_str());
        ofstream energy_stream(energy_file.c_str());
        ofstream sigmaE_stream(sigmaE_file.c_str());

        /* A boolean to check that all the streams are opened without an error */
        streams = plot_stream && Cv_stream && energy_stream && sigmaE_stream;

        if(streams)
        {
            /* Headers to write in the files */
            plot_stream << "#Parameters : nx = " << nx << "; ny = " << ny << "; nz = " << nz <<
            "; J0 = " << J0 << "; J1 = " << J1 << "; J2 = " << J2 << "; number of steps = " << nombrePas << endl;
            Cv_stream << "#Temperature Cv" << endl;
            energy_stream << "#Temperature Mean_Energy" << endl;
            sigmaE_stream << "#Temperature sigmaE^2" << endl;

            /********************************/
            /* Monte-Carlo temperature loop */
            /********************************/

            for(double T = Tinit; T >= Tinf; T -= temperatureStep) // loop over the temperature of the system to explore the phase diagram of the system. There is a -= in the loop because we go from high temperatures to low temperatures
            {
                cout << "J2: " << J2 << "; T: " << T << endl; // keeping track of the advancement of the algorythme while running

                /* We reinitialize the output variables */
                E = 0;
                sum_E = 0;
                sum_E_square = 0;
                DeltaE = 0;
                Cv = 0;
                
                /****************************************************/
                /********** Metropolis Algorythm algorithm **********/
                /****************************************************/
                
                for(int step = 0; step < nombrePas; step++) // loop over the Monte-Carlo steps
                {
                    /* We choose a random spin in the spin box */
                    x = rand() % nx; // x position
                    y = rand() % ny; // y position
                    z = rand() % nz; // z position

                    /* Calcul of the difference of energy ofthe system (DeltaE) if the spin is switched */
                    DeltaE = States.getDeltaE(x, y, z, J0, J1, J2);

                    /* Decision making : do we switch the spin or not ? */

                    if(DeltaE <= 0) // If the energy of the system after switching is lower than before switching...
                    {
                        States.switchValue(x, y, z); // We multiply by -1 : spin switched
                    }
                    else // If the energy of the system after switching is higher than before switching...
                    {
                        floattemp = rand() % 100000000 + 1; // float precision 1e8
                        nombreEntre0Et1 = floattemp / 100000000;
                        if(log(nombreEntre0Et1) < - DeltaE/T) // We generate a random number between 0 and 1 with a float precision and check if it is smaller than exp(- DeltaE/T)
                        {
                            States.switchValue(x, y, z); // If it is indeed smaller, we multiply by -1 : spin switched
                        } 
                        else // if not, nothing happens : spin not switched
                        {
                            DeltaE = 0;
                        }
                    }
                    
                    /* computation of Cv and sigmaE */
                    if(step == nombrePas/2) // We begin the computation of Cv and sigmaE on the second half of the Monte Carlo steps (because the system has converged now)
                    {
                        E = States.getEnergy(J0, J1, J2); // Get the energy of the system
                        sum_E = E; // initialisation of the sum
                        sum_E_square = E * E; // initialisation of the sum
                    }

                    if(step > nombrePas/2) // We keep the final half of the states to compute the averaged variables
                    {
                        E += DeltaE; // Update the energy of the system at each Monte-Carlo step
                        sum_E += E; // sum of the energy
                        sum_E_square += E * E; // sum of the squared energy
                    }

                }//end of the loop 'for' over nombrePas


                /*******************/
                /* Writing results */
                /*******************/
                
                /* <E> = sum_E / number_of_spins */
                energy_stream << T << " " << sum_E * 2 / nombrePas << endl;
                
                /* Cv = (1/kT^2) * (<E^2> - <E>^2) */
                Cv = (2/(nombrePas*T*T)) * (sum_E_square - sum_E * sum_E * 2 / nombrePas);
                Cv_stream << T << " " << Cv << endl;
                
                /* sigmaE^2 = (<E^2> - <E>^2) */
                sigmaE = ((sum_E_square - sum_E * sum_E * 2 / nombrePas) * 2 / nombrePas);
                sigmaE_stream << T << " " << sigmaE << endl;
                
                /* Writing configurations */
                // This file will save the curent state of the system for each T and J2 in order to allow to re-use the system for further analysis, or for the WL algorythm
                // It saves the configuration of each spin (up or down)
                States.write_config(T, J2);
                
            }//end of the loop 'for' over T

            /* Writes instructions in the 'plot' file for gnuplot */
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
            
            /* Closes the streams */
            Cv_stream.close();
            energy_stream.close();
            sigmaE_stream.close();
            plot_stream.close();
            
        }//end if(streams)
        else // if there is an error while openning the streams
        {
            cout << "Error while opening files" << endl;
        }
    } // end for J2



    /*******************************************/
    /*********** Config file for WL ************/
    /*******************************************/

    // This will write a configuration file for the Wang Landau Algorithme (in order to run it)

    // We will need the minimum and maximum energy of the system for each value of J2
    double Emin = 0; // minimum energy
    double Emax = 0; // maximum energy

    int numberPoints = 0; // the number of points/lines in the energy file over which we want to compute the average min and max energies (see below for more explanation)
    int totalPoints = 0; // the number of points/lines in the energy file
    int jumpPoints = 0; // number of lines we must jump in the middle of the file

    totalPoints = ((Tinit - Tinf) / temperatureStep) + 1;
    numberPoints = totalPoints * 0.05; // We want 5% of the total number of points/lines in the energy file
    jumpPoints = totalPoints - 2 * numberPoints; // We want to jump the "body" of the energy file, and keep the 5% of points at both extremities (highest temperatures and lowest temperatures, for which the energy of the system is highest or lowest, and almost constant but with fluctuations)

    /* This string allows us to write files names with the right convention */
    string const configWL_file("results/configWL_" + to_string(nx) + "x" + to_string(ny) + "x" +to_string(nz));

    /* Declaration of stream */
    ofstream configWL_stream(configWL_file.c_str());

    /* Header of the file */
    configWL_stream << "#J2 Emin Emax" << endl;

    if(configWL_stream) // check that the stream is opened without errors
    {
        for(double J2 = J2_min; J2 <= J2_max; J2 += J2_step) // We want Emin and Emax for each J2.
        {
            if(abs(J2) < 0.00001) // Sometimes, due to the errors of operations between double type numbers, J2 has the value of around 10^-16 instead of 0.
            {
                J2 = 0; // Therefore, we put J2 to 0 which should be its value.
            }
            
            // Reinitialisation of variables
            Emin = 0;
            Emax = 0;

            /* This string allows us to open the right file following the convention chosen at the writing */
            label = "_J2 = " + to_string(J2) + "_" + to_string(nx) + "x" + to_string(ny) + "x" +to_string(nz);

            string const read_energy_file("results/data/energy" + label + ".dat"); // We read the energies of the system from files written during the main Metropolis Algorythm. These energies depend on J2 and T. We will compute mean values over T.

            /* Reading stream declaration */
            ifstream read_energy_stream(read_energy_file.c_str());

            /* Declaration of temporary variables */
            string line; // allows to jump lines
            double T_buff; // stores T values temporarily
            double E_buff; // stores E values temporarily

            /***************************/
            /* Begining of the reading */
            /***************************/

            getline(read_energy_stream, line); // jump the first comment line

            // Computation of the mean max Energy (we assume that for high temperatures, the system is paramagnetic, thus, its energy is "maximum", and stays "constant" for these high temperatures, with some fluctuations, because we do the computation over the 5% points for which the temperature is the highest.
            // We compute the mean energy of this maximum over these high temperatures to "get ride of" these fluctuations.)
            for(int i = 0; i < numberPoints ; i++)
            {
                read_energy_stream >> T_buff; // Temperature of the system
                read_energy_stream >> E_buff; // Energy of the system

                Emax += E_buff; // Sum of the energyes of the system
            }

            Emax /= numberPoints; // mean value of the maximum over the 5% points corresponding to the highest temperatures

            getline(read_energy_stream, line); // jump to the next line, because after the last >>, we were at the end of the line, and not at the beginning of the following line

            //jump some lines : the "body" of the energy file that doesn't interest us here
            for(int i = 0 ; i < jumpPoints ; i++)
            {
                getline(read_energy_stream, line);
            }

            // Computation of the mean min Energy (the system is magnetic, and its energy is minimum here, and we assume this energy to be almost constant over the 5% points on which we do the computation)
            for(int i = 0; i < numberPoints ; i++)
            {
                read_energy_stream >> T_buff; // Temperature of the system
                read_energy_stream >> E_buff; // Energy of the system

                Emin += E_buff; // Sum of the energyes of the system
            }

            Emin /= numberPoints; // mean value of the minimum over the 5% points corresponding to the lowest temperatures

            // writting in configWL file
            configWL_stream << J2 << " " << Emin << " " << Emax << endl;

        } // end of the loop 'for' over J2

    } // end of the if(configWL_stream)

    else // if an error occured during the opening of the stream
    {
        cout << "Error while opening configWL file !" << endl;
    }

    return 0;
}
