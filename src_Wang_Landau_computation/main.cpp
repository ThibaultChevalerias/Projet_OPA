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
    /* ================ Reading g(E) from file ================ */
    /* ======================================================== */

    vector< pair<double, double> > gE; // Will contain the energy (mean energy of a bin) and the value of g(E) for this energy

    string const read_file("g(E).dat"); // path may be modified to read the appropriate file. This file is the output file of the Wang-landau algorithm code 
    ifstream read_stream(read_file.c_str()); // reading stream
    if(read_stream)
    {
        string line; // to jump the first and second lines
        getline(read_stream, line); // to jump the first line
        getline(read_stream, line); // to jump the second line

        int number_lines = 0;
        read_stream >> number_lines; // read the first important line of the file : the number of data lines

        double energy_temp = 0;
        double gE_temp = 0;

        for(int i = 0; i < number_lines; i ++) //read from the fourth to the last line
        {
            read_stream >> energy_temp;
            read_stream >> gE_temp;
            gE.push_back(pair<double, double>(energy_temp, gE_temp));
        }

        read_stream.close();

    } // end of if(read_stream)
    else
    {
        cout << "Error while reading g(E) file !!!" << endl;
    }


    /* ========================================================== */
    /* =================== Computation of Cv ==================== */
    /* ========================================================== */

    /** TEMPERATURES À MIEUX EXPLICITER **/
    double Tinit = 1; //Initial temperature (implicit kbT with kb=1)
    double Tfinal = 100; // Final temperature (implicit kbT with kb=1)
    double Tstep = 1; // Temperature step (implicit kbT with kb=1)

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

            for(int i = 0; i < gE.size(); i++)
            {
                gEfE = exp(gE[i].second - gE[i].first / T); // represents g(E)f(E) = exp(ln(g(E)) - beta * E)
                Z += gEfE;
                mean_energy += gEfE * gE[i].first;
                mean_square_energy += gEfE * gE[i].first * gE[i].first;
            }
            
            mean_energy /= Z; // normalization by the partition function
            mean_square_energy /= Z; // normalization by the partition function
            
            /** Computation of Cv : 
            Cv = (1/k_B*T^2) * (<E^2> - <E>^2)
            **/
            Cv = (mean_square_energy - mean_energy * mean_energy) / (T * T); // The k_B is omitted (arbitrary values)
            
            T_Cv.push_back(pair<double, double>(T, Cv));
            Cv_stream << T << " " << Cv << endl;

        } // end of for(int T = Tinit; T < Tfinal; T += Tstep)

        /* Computation of the maximum of Cv and determination of the transition temperature Tc */
        int max = T_Cv[0].second;
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