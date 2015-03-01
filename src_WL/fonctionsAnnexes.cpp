#include <iostream>
#include <cmath>
#include "fonctionsAnnexes.h"

using namespace std;

int mod(int a, int b) // We assume b to be positive and a to vary between -2 and b+2
{
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
    return floor((energy - E_min) / deltaE); // euclidian division of energy by deltaE to find in which bin "energy" is. We add (-E_min) in order to have a positive number and divide it by deltaE to get a number between 0 and number_bins - 1
}

int getMax(vector<int>& table) // It only works for integer tables
{
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
    double min = 0;
    int index = 0;
    while (table[index] == 0)
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

double getMean(vector<int>& table) // /!\ eliminates the 0 values in the table from the computation
{
    double mean = 0;
    int count = 0;

    for(unsigned int i = 0; i < table.size(); i++)
    {
        if(table[i] != 0)
        {
            mean += table[i];
            count++;
        }
    }
    if(count != 0)
    {
        mean /= count;
        return mean;
    }
    else
    {
        cout << "bad energy window, energy out of range" << endl;
        return 0;
    }
}


bool isFlat(double flatness_limit, vector<int>& visits) // It only works for integer tables
{
    double mean = 0;
    double min = 0;
    double max = 0;
    double deviation = 0;

    mean = getMean(visits);
    min = getMin(visits);
    max = getMax(visits);

    if ((max - mean) > (mean - min))
    {
        deviation = max - mean;
    }
    else
    {
        deviation = mean - min;
    }
    
    if(1 - (deviation/mean) > flatness_limit)
    {
        return true;
    }
    else
    {
        return false;
    }

}
