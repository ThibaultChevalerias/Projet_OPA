#include <iostream>
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
