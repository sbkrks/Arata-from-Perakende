#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <complex>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>

#include "../head/SabaVector.h"

static int j, k, l, m, n, p, q;
static double Tr, detA, sum;
static double* ab, a, b, c;
static double** A, B, C;

double** MatrixTranspose(double** A) { //Calculating the transpose of 3x3 matrix.
	double** B;
	B = new double*[3];
	for (j = 0; j <= 2; j++) {
		B[j] = new double[3];
	}
	double x;

	for (j = 0; j <= 2; ++j) {
		for (k = 0; k <= 2; ++k) {
			x = A[j][k];
			B[k][j] = x;
		}
	}
	return B;
}

double MatrixTrace(double** A) {
	double Tr = 0;

	for (j = 0; j <= 2; ++j) {
		Tr += A[j][j];
	}

	return Tr;
}

double MatrixDeterminant(double** A) { //Calculating the determinant of 3x3 matrix.
	double detA = 0;

	for (j = 0; j <= 2; j++) {
		for (k = 0; k <= 2; k++) {
			for (l = 0; l <= 2; l++) {
				detA += LeviCivita(j, k, l) * A[0][j] * A[1][k] * A[2][l];
			}
		}
	}
	return detA;
}

double** MatrixInverse(double** A) { //Calculating the inverse of the matrix if exists.
	double detA = MatrixDeterminant(A);

	double** B = new double*[3];
	for (j = 0; j <= 2; ++j) {
		B[j] = new double[3];
	}
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//                                                                                                              //
	//          B is the inverse matrix of A.                                                                       //
	//                                                                                                              //
	//                          1       ___  ___  ___  ___                                                          //
	//           B    =   _____________ \    \    \    \   Epsilon    Epsilon    A   A                              //
	//            ij        2 det(A)    /__  /__  /__  /__        jmn        ipq  mp  nq                            //
	//                                   m    n    p    q                                                           //
	//                                                                                                              //
	//ref: http://www.continuummechanics.org/cm/matrices.html                                                       //
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	double sum;
	for (j = 0; j <= 2; ++j) {
		for (k = 0; k <= 2; ++k) {
			sum = 0;
			for (m = 0; m <= 2; ++m) {
				for (n = 0; n <= 2; ++n) {
					for (p = 0; p <= 2; ++p) {
						for (q = 0; q <= 2; ++q) {
							sum += LeviCivita(k, m, n) * LeviCivita(j, p, q) * A[m][p] * A[n][q] / (2 * detA);
						}
					}
				}
			}
			B[j][k] = sum;
		}
	}

	return B;
}

double** MatrixProduct(double** A, double** B) {
	double** C = new double*[3];
	for (j = 0; j <= 2; ++j) {
		C[j] = new double[3];
	}

	double sum;

	for (j = 0; j <= 2; ++j) {
		for (k = 0; k <= 2; ++k) {
			sum = 0;
			for (l = 0; l <= 2; ++l) {
				sum += A[j][l] * B[l][k];
			}
			C[j][k] = sum;
		}

	}

	return C;
}

double** MatrixConjugate(double** A) {
	double** B = MatrixTranspose(A);

	return B;
}