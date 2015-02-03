#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>

using namespace std;

int main()
{
    /* ======================================================== */
    /* ============== Reading ln(g(E)) from file ============== */
    /* ======================================================== */

    vector< pair<double, double> > lngE; // Will contain the energy (mean energy of a bin) and the value of lng(E) for this energy
    
    double J2 = - 0.4;
    
    string const read_file("results/lng(E)_J2=" + to_string(J2) + ".dat"); // path may be modified to read the appropriate file. This file is the output file of the Wang-landau algorithm code 
    ifstream read_stream(read_file.c_str()); // reading stream
    
    int nxnynz = 0; // cannot be declared in the next "if" because it is used later
    // nxnynz = nx * ny * nz (number of sites)
    double rescale = 0; // This value had been added to the energies in WL, so it is to be substracted.
    double lngE_max = 0; // cannot be declared in the next "if" because it is used later
    // Maximum value of ln(g(E))
    
    if(read_stream)
    {
        string line; // to jump the comment lines

        getline(read_stream, line); // to jump the first line (comment)
        int number_lines = 0;
        read_stream >> number_lines; // read the first important line of the file : the number of data lines
        
        getline(read_stream, line); // to jump the comment lines
        getline(read_stream, line); // to jump the comment lines
        read_stream >> nxnynz; // nxnynz = nx * ny * nz (number of sites)
        
        getline(read_stream, line); // to jump the comment lines
        getline(read_stream, line); // to jump the comment lines
        read_stream >> rescale; // rescale of the energy
        
        getline(read_stream, line); // to jump the comment lines
        getline(read_stream, line); // to jump the comment lines
        
        double energy_temp = 0;
        double lngE_temp = 0;
        
        for(int i = 0; i < number_lines; i ++) //read from the fourth to the last line
        {
            read_stream >> energy_temp;
            read_stream >> lngE_temp;
            lngE.push_back(pair<double, double>(energy_temp, lngE_temp));
            
            /* Calculation of the maximum : */
            if(lngE_temp > lngE_max)
            {
                lngE_max = lngE_temp;
            }
        }

        read_stream.close();
        cout << "Data acquired, beginning calculation..." << endl;
    } // end of if(read_stream)
    else
    {
        cout << "Error while reading lng(E) file !!!" << endl;
    }


    /* ========================================================== */
    /* =================== Computation of Cv ==================== */
    /* ========================================================== */

    /** TEMPERATURES À MIEUX EXPLICITER **/
    double Tinit = 0.001; //Initial temperature (implicit kbT with kb=1)
    double Tfinal = 10; // Final temperature (implicit kbT with kb=1)
    double Tstep = 0.001; // Temperature step (implicit kbT with kb=1)

    long double Cv = 0; // heat capacity
    long double Z = 0; // partition function
    long double mean_energy = 0;
    long double mean_square_energy = 0;
    long double gEfE = 0; // represents g(E)f(E) = exp(ln(g(E)) - beta * E)

    string const Cv_file("results/Cv_J2=" + to_string(J2) + ".dat"); // output file. Will allow to plot Cv in function of T, with gnuplot for instance
    ofstream Cv_stream(Cv_file.c_str()); // output stream
    
    string const energy_file("results/energy_J2=" + to_string(J2) + ".dat"); // output file. Will allow to plot Cv in function of T, with gnuplot for instance
    ofstream energy_stream(energy_file.c_str()); // output stream
    
    if(Cv_stream && energy_stream)
    {
        Cv_stream << "#T Cv" << endl;
        energy_stream << "#T energy" << endl;
        
        for(double T = Tinit; T <= Tfinal; T += Tstep)
        {
            Cv = 0;
            Z = 0;
            mean_energy = 0;
            mean_square_energy = 0;
            
            for(unsigned int i = 0; i < lngE.size(); i++)
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
        } // end of for(int T = Tinit; T < Tfinal; T += Tstep)
        
        Cv_stream.close();
        cout << "Calculation done!" << endl;    
    } // end of if(Cv_stream)
    else
    {
        cout << "Error while writing Cv.dat file !!!" << endl;
    }

    return 0;

}
