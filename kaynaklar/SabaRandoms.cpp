#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <complex>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>

static int wi, wo;
static double a, b, lowerlimit, upperlimit;

double Randf() { //Generating double random number in between (0,1)
	return((double)(rand()) / (double)(RAND_MAX + 1.0));
}

int Randi(int wi){ //Generating an integer
	int wo;
	wo = static_cast<int>(static_cast<double>(rand()) / static_cast<double>(RAND_MAX + 1.0)*wi);

return wo;
}

double SpecRandf(double lowerlimit, double upperlimit){ //Generating double random number in between (lowerlimit,upperlimit)
    double b = Randf();
    double a = (upperlimit - lowerlimit)*b + lowerlimit;

return a;
}
