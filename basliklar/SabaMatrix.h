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

double** MatrixTranspose(double** A);

double MatrixTrace(double** A);

double MatrixDeterminant(double** A);

double** MatrixInverse(double** A);

double** MatrixProduct(double** A, double** B);

double** MatrixConjugate(double** A);