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

    string const read_file("lng(E).dat"); // path may be modified to read the appropriate file. This file is the output file of the Wang-landau algorithm code 
    ifstream read_stream(read_file.c_str()); // reading stream
    
    int nxnynz = 0; // cannot be declared in the next "if" because it is used later
    // nxnynz = nx * ny * nz (number of sites)
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

    } // end of if(read_stream)
    else
    {
        cout << "Error while reading lng(E) file !!!" << endl;
    }


    /* ========================================================== */
    /* =================== Computation of Cv ==================== */
    /* ========================================================== */

    /** TEMPERATURES À MIEUX EXPLICITER **/
    double Tinit = 0.01; //Initial temperature (implicit kbT with kb=1)
    double Tfinal = 1; // Final temperature (implicit kbT with kb=1)
    double Tstep = 0.01; // Temperature step (implicit kbT with kb=1)

    vector< pair<double, double> > T_Cv; // Will contain the value of the temperature T and the corresponding value of the specific heat Cv
    double Cv = 0; // heat capacity
    double Z = 0; // partition function
    double mean_energy = 0;
    double mean_square_energy = 0;
    double gEfE = 0; // represents g(E)f(E) = exp(ln(g(E)) - beta * E)

    string const Cv_file("Cv.dat"); // output file. Will allow to plot Cv in function of T, with gnuplot for instance
    ofstream Cv_stream(Cv_file.c_str()); // output stream

    if(Cv_stream)
    {
        Cv_stream << "#T Cv" << endl;

        for(double T = Tinit; T <= Tfinal; T += Tstep)
        {
            Cv = 0;
            Z = 0;
            mean_energy = 0;
            mean_square_energy = 0;

            for(int i = 0; i < lngE.size(); i++)
            {
                gEfE = exp(lngE[i].second - lngE_max - lngE[i].first / T); // represents g(E)f(E) = exp(ln(g(E)) - ln(g(E))_max - beta * E)
                /* the "- ln(g(E))_max" is to avoid big numbers in the exponential. It does not change anything in the computation of Cv because
                it is cancelled out by dividing by Z later */
                Z += gEfE;
                mean_energy += gEfE * lngE[i].first;
                mean_square_energy += gEfE * lngE[i].first * lngE[i].first;
            }
            
            mean_energy /= Z; // normalization by the partition function
            mean_square_energy /= Z; // normalization by the partition function
            
            /** Computation of Cv : 
            Cv = (1/k_B*T^2) * (<E^2> - <E>^2)
            **/
            Cv = (mean_square_energy - mean_energy * mean_energy) / (T * T); // The k_B is omitted (arbitrary values)
            
            Cv *= nxnynz; // We had the energy per site, so the Cv is per site^2 (because there are quadratic terms of the energy in it).
            // So we multiply it by the number of sites.
            
            T_Cv.push_back(pair<double, double>(T, Cv));
            Cv_stream << T << " " << Cv << endl;

        } // end of for(int T = Tinit; T < Tfinal; T += Tstep)
        
        /* ============================================================================================== */
        /* ==== Computation of the maximum of Cv and determination of the transition temperature Tc ===== */
        /* ============================================================================================== */
        double max = T_Cv[0].second;
        int index_max = 0;
        for (int i = 1; i < T_Cv.size(); i++)
        {
            if (T_Cv[i].second > max)
            {
                max = T_Cv[i].second;
                index_max = i;
            }
        }

        Cv_stream << "# The maximum of Cv occurs at the temperature Tc = " << T_Cv[index_max].first << endl;

        Cv_stream.close();
    
    } // end of if(Cv_stream)
    else
    {
        cout << "Error while writing Cv.dat file !!!" << endl;
    }

    system("PAUSE");
    return 0;

}