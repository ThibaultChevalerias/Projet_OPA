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
    /* Declaration of variables */
    int nx = 5;
    int ny = 5;
    int nz = 5;
    
    /* Exchange constants */
    double J0 = 1;
    double J1 = 1;
    double J2 = -1;
    
    double lnf = 1;
    
    double flatness = 0; // flatness of the histogram
    double flatness_limit = 0.9; // When we reach this limit, we consider the histogram as flat, and change the value of f.
    
    int step = 0;
    int step_max = 1000000; // Maximum number of steps allowed for the flatness to pass above flatness_limit
    
    /* Variables to choose a random spin */
    int xChosen = 0;
    int yChosen = 0;
    int zChosen = 0;
    
    /* Variables randomly generated */
    double random_double = 0;
    int random_range = 100000000; // float precision (1E8)
    
    /* Energy bins */
    /* A COMPLETER E_min = 
    E_max = */
    number_bins = 100; // We cut the energy interval in 100 bins
    
    cout << "Wait while the simulation is running..." << endl;
    
    while(lnf > epsilon)
    {
        while(flatness < flatness_limit && step < step_max)
        {
            
        }
        
        lnf /= 2; // We reduce f until f < epsilon
    }
}
