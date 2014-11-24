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

    }
    else
    {
        cout << "Error while reading g(E) file !!!" << endl;
    }


    /* ========================================================== */
    /* =================== Computation of Cv ==================== */
    /* ========================================================== */

    for(int i = 0; i < 50; i ++)
    {
        cout << gE[i].first << " " << gE[i].second << endl;
    }

    system("PAUSE");
    return 0;

}