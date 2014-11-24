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

    vector< pair<double, double> > gE;

    string const read_file("g(E).dat"); // may be modified 
    ifstream read_stream(read_file.c_str());
    if(read_stream)
    {
        string line; // to jump the first and second lines
        getline(read_stream, line); // to jump the first line
        getline(read_stream, line); // to jump the second line

        int number_lines = 0;
        read_stream >> number_lines; // read the first important line : the number of data lines

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

    double Tinit = 0.001; //Initial temperature (implicit kbT with kb=1)
    double Tfinal = 1; // Final temperature (implicit kbT with kb=1)
    double Tstep = 0.001; // Temperature step (implicit kbT with kb=1)

    double Cv = 0;

    string const Cv_file("Cv.dat");
    ofstream Cv_stream(Cv_file.c_str());

    if(Cv_stream)
    {
        Cv_stream << "#T Cv" << endl;

        for(double T = Tinit; T <= Tfinal; T += Tstep)
        {
            Cv = 0;

            for(int i = 0; i < gE.size(); i++)
            {
                Cv += gE[i].first * gE[i].first * gE[i].second * exp(-gE[i].first / T) / (T * T);
            }

            Cv_stream << T << " " << Cv << endl;

        } // end of for(int T = Tinit; T < Tfinal; T += Tstep)

        Cv_stream.close();

    } // end of if(Cv_stream)
    else
    {
        cout << "Error while writing Cv.dat file !!!" << endl;
    }

    system("PAUSE");
    return 0;

}