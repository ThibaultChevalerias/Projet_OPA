#ifndef FONCTIONSANNEXES_H_INCLUDED
#define FONCTIONSANNEXES_H_INCLUDED

int mod(int a, int b);
int locateBin(double deltaE, double energy);
int getMax(int *table, int size);
int getMin(int *table, int size)
double getMean(int *table, int size);
bool isFlat(double flatness_limit, int *visits, int size);

#endif // FONCTIONSANNEXES_H_INCLUDED
