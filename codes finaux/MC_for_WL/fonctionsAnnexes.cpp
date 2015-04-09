/***
This file contains useful functions for the main code.
***/
#include <iostream>
#include <cmath> // We use the math function 'floor'.
#include <vector>
#include "fonctionsAnnexes.h"

using namespace std;


int mod(int a, int b) // We assume 'b' to be positive, and 'a' to vary between -2 and b+2
{
    /***
    This function returns the number in [0;b], equal to 'a' modulo 'b'.
    ***/
    
    if(a>=0)
    {
        return a%b;
    }
    else
    {
        return b + a%b;
    }
}

int locateBin(double E_min, double deltaE, double energy)
{
    /***
    This function returns the range (or "bin") in energy in which the variable 'energy' is. 'E_min' is the minimum energy and 'deltaE' is the width of an energy range (or "bin").
    ***/
    return floor((energy - E_min) / deltaE); // Euclidian division of energy by deltaE to find in which bin "energy" is. We add (-E_min) in order to have a positive number and divide it by deltaE to get a number between 0 and number_bins - 1
}

int getMax(vector<int>& table) // It only works for integer tables
{
    /***
    This function returns the maximum of a table of integers (vector<int>).
    ***/
    double max = table[0];

    for(unsigned int i = 1; i < table.size(); i++)
    {
        if(max < table[i])
        {
            max = table[i];
        }
    }
    return max;

}

int getMin(vector<int>& table) // It only works for integer tables
{
    /***
    This function returns the *non-zero* minimum of a table of integers (vector<int>).
    ***/
    
    double min = 0;
    int index = 0;
    while (table[index] == 0) // We skip the zeros in the table (we assume there is at least one non-zero number in the table).
    {
        index ++;
    }

    min = table[index];

    for(unsigned int i = index + 1; i < table.size(); i++)
    {
        if(min > table[i] && table[i] != 0)
        {
            min = table[i];
        }
    }

    return min;
}

double getMean(vector<int>& table) // /!\ Eliminates the 0 values in the table from the computation
{
    /***
    This function returns the mean of a table of integers (vector<int>).
    ***/
    
    double mean = 0;
    int count = 0;

    for(unsigned int i = 0; i < table.size(); i++) // Calculation of the mean. We remove the zeros !
    {
        if(table[i] != 0)
        {
            mean += table[i];
            count++;
        }
    }
    if(count != 0) // If there has been at least one non-zero value, everything is fine.
    {
        mean /= count;
        return mean;
    }
    else // If there has been only zeros, count is equal to zero so there is probably a problem with the energy range.
    {
        cout << "Bad energy window, energy out of range." << endl;
        return 0;
    }
}


bool isFlat(double flatness_limit, vector<int>& visits) // It only works for integer tables
{
    /***
    This function returns if a table of integers (vector<int>) is "flat" ('true or 'false'), according to the 'flatness_limit' criterion.
    The limit is given as a number between 0 and 1. For instance, if 'flatness_limit' is equal to 0.8, the criterion for flatness is 80%. Every value has to be within a 20% range around the mean.
    ***/
    
    double mean = 0;
    double min = 0;
    double max = 0;
    double deviation = 0;

    mean = getMean(visits);
    min = getMin(visits);
    max = getMax(visits);

    if ((max - mean) > (mean - min)) // We choose the biggest deviation to be the 'deviation' variable.
    {
        deviation = max - mean;
    }
    else
    {
        deviation = mean - min;
    }
    
    if(1 - (deviation/mean) > flatness_limit) // If we respect the flatness criterion, we return 'true', else it's 'false'.
    // The criterion is that the deviation must be less than mean*(1-flatness_limit).
    // For example, if flatness_limit==0.8, 'deviation' has to be less than 0.2*mean, or 20% of the mean.
    {
        return true;
    }
    else
    {
        return false;
    }

}
