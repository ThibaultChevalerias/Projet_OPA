#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>

using namespace std;

/* This code will compute and write in a file the values of the specific heat Cv in function of the temperature T */
/* It uses the lngE file writen by the WL algorithm in order to compute Cv */

int main()
{
    /* ======================================================== */
    /* ============== Reading ln(g(E)) from file ============== */
    /* ======================================================== */

    vector< pair<double, double> > lngE; // Will contain the energy E (mean energy of a bin) of the system and the value of lng(E) for this energy E
    
    /* Size of the spin box */
    // /!\ /!\ /!\ Beware : it must be the same size as the size of the system used in the WL algorithm.
    // The size of the system is in the name of the different files to help prevent errors. But be carefull. Always check it.
    int nx = 10; // size of the system following the x axis
    int ny = 10; // size of the system following the y axis
    int nz = 56; // size of the system following the z axis
    
    double J2 = - 1; // exchange parameter between next nearest neighbours
                     // /!\ You must change this parameter yourself to explore the phase diagram. There is no loop over J2 in this code.
                     // There is also a protection in the file naming convention to prevent selecting the wrong file
    
    /* reading file following the naming convention */
    string const read_file("results/lng(E)_J2=" + to_string(J2) + "_" + to_string(nx) + "x" + to_string(ny) + "x" + to_string(nz) + ".dat"); // path may be modified to read the appropriate file. This file is the output file of the Wang-landau algorithm code. It is here used as an input to recover the values of lng(E)
                                                                                                                                             // Normaly, the right file to read is selected automaticaly thanks to the naming convention
    ifstream read_stream(read_file.c_str()); // reading stream
    
    /* some more variables declarations */
    int nxnynz = 0; // cannot be declared in the next "if" because it is used later. It will be filled with the value of nx * ny * nz
    double rescale = 0; // This value had been added to the energies in WL, so it is to be substracted.
    double lngE_max = 0; // cannot be declared in the next "if" because it is used later
                         // It will be the maximum value of ln(g(E))
    
    /* begining reading data from file */
    if(read_stream)
    {
        string line; // to jump the comment/header lines

        getline(read_stream, line); // to jump the first line (comment)
        int number_lines = 0;
        read_stream >> number_lines; // read the first important line of the file : the number of data lines in the file (excluding the headers)
        
        getline(read_stream, line); // to jump the comment lines
        getline(read_stream, line); // to jump the comment lines
        read_stream >> nxnynz; // nxnynz = nx * ny * nz (number of sites)
        
        getline(read_stream, line); // to jump the comment lines
        getline(read_stream, line); // to jump the comment lines
        read_stream >> rescale; // rescale of the energy
        
        getline(read_stream, line); // to jump the comment lines
        getline(read_stream, line); // to jump the comment lines
        
        double energy_temp = 0; // temporary variable to transfer data from the file to a vector for easier use
        double lngE_temp = 0; // temporary variable to transfer data from the file to a vector for easier use
        
        for(int i = 0; i < number_lines; i ++) //read from the fourth to the last data line (we don't count headers as data lines and there are no more headers in the body of this file anyway)
        {
            read_stream >> energy_temp;
            read_stream >> lngE_temp;
            lngE.push_back(pair<double, double>(energy_temp, lngE_temp));
            
            /* Computation of the maximum of lngE */
            if(lngE_temp > lngE_max)
            {
                lngE_max = lngE_temp;
            }
        }

        read_stream.close();
        cout << "Data acquired, beginning calculation..." << endl;
    } // end of if(read_stream)
    else // if an error occured while opening the lngE file...
    {
        cout << "Error while reading lng(E) file !!!" << endl;
    }


    /* ========================================================== */
    /* =================== Computation of Cv ==================== */
    /* ========================================================== */

    /* Temperature loop parameters */
    double Tinit = 4; //Initial temperature (implicit kbT with kb=1)
    double Tfinal = 10; // Final temperature (implicit kbT with kb=1)
    double Tstep = 0.001; // Temperature step (implicit kbT with kb=1)

    /* Declaration of some variables */
    long double Cv = 0; // specific heat
    long double Z = 0; // partition function. We are in the Maxwell Boltzman approximation
    long double mean_energy = 0; // for Cv computation
    long double mean_square_energy = 0; // for Cv computation
    long double gEfE = 0; // represents g(E)f(E) = exp(ln(g(E)) - beta * E)

    /* opening the right files thanks to the naming convention */
    string const Cv_file("results/Cv_J2=" + to_string(J2) + "_" + to_string(nx) + "x" + to_string(ny) + "x" + to_string(nz) + ".dat"); // output file. Will allow to plot Cv in function of T, with gnuplot for instance
    ofstream Cv_stream(Cv_file.c_str()); // output stream
    
    string const energy_file("results/energy_J2=" + to_string(J2) + "_" + to_string(nx) + "x" + to_string(ny) + "x" + to_string(nz) + ".dat"); // output file. Will allow to plot Cv in function of T, with gnuplot for instance
    ofstream energy_stream(energy_file.c_str()); // output stream
    
    if(Cv_stream && energy_stream) // checks that no error occured while opening the files
    {
        Cv_stream << "#T Cv" << endl;
        energy_stream << "#T energy" << endl;
        
        for(double T = Tinit; T <= Tfinal; T += Tstep) // loop over the temperature T to explore the phase diagram for a given value of J2
        {
            Cv = 0; // re-initialisation for each temperature
            Z = 0; // re-initialisation for each temperature
            mean_energy = 0; // re-initialisation for each temperature
            mean_square_energy = 0; // re-initialisation for each temperature
            
            for(unsigned int i = 0; i < lngE.size(); i++) // Now we use the values of lngE
            {
                gEfE = exp(lngE[i].second - lngE_max - lngE[i].first / T); // represents g(E)f(E) = exp(ln(g(E)) - ln(g(E))_max - beta * E)
                /* the "- ln(g(E))_max" is to avoid big numbers in the exponential. It does not change anything in the computation of Cv because
                it is cancelled out by dividing by Z later */
                Z += gEfE;
                mean_energy += gEfE * lngE[i].first; // We rescale the energy back here, so that the rescale has not any influence on the calculations (do the calculations)
                mean_square_energy += gEfE * lngE[i].first * lngE[i].first;
            }
            
            mean_energy /= Z; // normalization by the partition function
            mean_square_energy /= Z; // normalization by the partition function
            
            /** Computation of Cv : 
            Cv = (1/k_B*T^2) * (<E^2> - <E>^2)
            **/
            Cv = (mean_square_energy - mean_energy * mean_energy) / (T * T); // The k_B is omitted (arbitrary values)
            
            Cv_stream << T << " " << Cv << endl;
            energy_stream << T << " " << mean_energy << endl;
        } // end of the loop for(int T = Tinit; T < Tfinal; T += Tstep)
        
        Cv_stream.close();
        cout << "Calculation done!" << endl;    
    } // end of if(Cv_stream)
    else // If there is an error while opening the file
    {
        cout << "Error while writing Cv.dat file !!!" << endl;
    }

    return 0;

}
