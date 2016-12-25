#pragma once
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <complex>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>

#include "SabaVector.h"
#include "SabaMatrix.h"
#include "SabaRandoms.h"

using namespace std;

complex<double> MatrixDeterminant(complex<double>** A);

complex<double> MatrixTrace(complex<double>** A);

complex<double>** MatrixProduct(double** A, complex<double>** B);
complex<double>** MatrixProduct(complex<double>** A, double** B);
complex<double>** MatrixProduct(complex<double>** A, complex<double>** B);
complex<double>** MatrixProduct(complex<double>** A, complex<double>** B, complex<double>** C);

complex<double>** MatrixSum(complex<double>** A1, complex<double>** A2, complex<double>** A3, complex<double>** A4, complex<double>** A5, complex<double>** A6);

complex<double>** MatrixDifference(complex<double>** A1, complex<double>** A2);

complex<double>** MatrixConjugate(complex<double>** A);

complex<double>** SU3_Rand_Gnrtr();