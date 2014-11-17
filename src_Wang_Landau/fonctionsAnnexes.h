#ifndef FONCTIONSANNEXES_H_INCLUDED
#define FONCTIONSANNEXES_H_INCLUDED

#include <vector>

using namespace std;

int mod(int a, int b);
int locateBin(double deltaE, double energy);
int getMax(vector<int>& table);
int getMin(vector<int>& table);
double getMean(vector<int>& table);
bool isFlat(double flatness_limit, vector<int>& visits);

#endif // FONCTIONSANNEXES_H_INCLUDED
