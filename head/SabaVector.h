#pragma once
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <complex>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>

double SpecRandf(double lowerlimit, double upperlimit);

double LeviCivita(int k, int l, int m);

double VecSqr(double* u, int n);

double VecNorm(double* u, int n);

double VecDot(double* u, double* v, int n);

double* VecCross3D(double* u, double* v);

double* GramSchmidt(int n);
