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

int locateBin(double deltaE, double energy)
{
    return floor(energy / deltaE); // euclidian division of energy by deltaE to find in which bin "energy" is
}

int getMax(int *table, int size) // It only works for integer tables
{
    double max = *table[0];

    for(int i = 1; i < size; i++)
    {
        if(max < *table[i])
        {
            max = *table[i];
        }
    }
    return max;

}

int getMin(int *table, int size) // It only works for integer tables
{
    double min = *table[0];

    for(int i = 1; i < size; i++)
    {
        if(min > *table[i])
        {
            min = *table[i];
        }
    }
    return min;

}

double getMean(int *table, int size) // /!\ eliminates the 0 values in the table from the computation
{
    double mean = 0;
    int count = 0;

    for(int i = 0; i < size; i++)
    {
        if(*table[i] != 0)
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


bool isFlat(double flatness_limit, int *visits, int size) // It only works for integer tables
{
    double mean = 0;
    double range = 0;

    mean = getMean(visits,size);

    range = getMax(visits,size) - getMin(visits,size);

    if(1-(range/mean) > flatness_limit)
    {
        return true;
    }
    else
    {
        return false;
    }

}