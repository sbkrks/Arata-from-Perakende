#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <complex>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>

#include "../head/SabaRandoms.h"
using namespace std;

static int j, k, l, m, n, p, q;
static double lowerlimit, upperlimit;
static double a, b, LC, A, SV, NV, sum, u_norm, v_prime_norm;
static double* u, v, w, v_prime, u_new, v_new, rtrn;

double LeviCivita(int k, int l, int m){ //Returning +1.0 or 0.0 or -1.0 depends on the permutation of [klm]
    double LC = (k - l) * (l - m) * (m - k)/2;

return LC;
}

double VecSqr(double* u, int n){ //Calculating the square of given vector of n dimensional.
    double A = 0;

    for(j = 0; j <= n-1; ++j){
        A += u[j] * u[j];
    }

return A;
}

double VecNorm(double* u, int n){ //Calculating the magnitude of given vector of n dimensional.

    double NV = sqrt(VecSqr(u, n));

return NV;
}

double VecDot(double* u, double* v, int n){ //Calculating the dot product of two given vectors of n dimensional.
    double sum = 0;

    for(j = 0; j <= n-1; ++j){
        sum += u[j] * v[j];
    }

return sum;
}

double* VecCross3D(double* u, double* v){ //Calculating the vector product of two 3D vectors and returning the outcome 3D vector.
    double* w = new double[3];

    for(j = 0; j <= 2; ++j){
        w[j] = 0;
        for(k = 0; k <= 2; ++k){
            for(l = 0; l <= 2; ++l){
                w[j] += LeviCivita(j,k,l) * u[k] * v[l];
            }
        }
    }

return w;
}

double* GramSchmidt(int n) { //Creating two orthogonal vectors of n dimensional, returning them as a 1D array, first n element is the first vector, second n elements is the second vector.
	double* u = new double[n];
	double* v = new double[n];
	double* v_prime = new double[n];
	double* u_new = new double[n];
	double* v_new = new double[n];
    double* rtrn = new double[2*n];

    for(j = 0; j <= n-1; ++j){
        u[j] = Randf();
        v[j] = Randf();
    }

    double u_norm = VecNorm(u, n);

    for(j = 0; j <= n - 1; ++j){
        u_new[j]=u[j]/u_norm;
    }

	delete[] u;

    double A = VecDot(v, u_new, n);

    for(j = 0; j <= n - 1 ; ++j){
        v_prime[j] = v[j] - A * u_new[j];
    }

	delete[] v;

    double v_prime_norm = VecNorm(v_prime, n);

    for(j = 0; j <= n - 1; ++j){
        v_new[j]=v_prime[j]/v_prime_norm;
    }

    for(j = 0; j <= n - 1; ++j){
        rtrn[j] = u_new[j];
        rtrn[j+n] = v_new[j];
    }

	delete[] u_new;
	delete[] v_new;
	delete[] v_prime;

return rtrn;
}
