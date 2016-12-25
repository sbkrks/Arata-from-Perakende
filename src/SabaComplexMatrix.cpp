#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <complex>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>

#include "../head/SabaVector.h"
#include "../head/SabaMatrix.h"
#include "../head/SabaRandoms.h"
#include "../head/Parameters.h"

using namespace std;

static int j, k, l;
static double alpha1, alpha2, alpha3;
static double* ab, a2, b, c;
static double** P, PI;
static complex<double> detA, sum, z1, z2, z3;
static complex<double>** A, B, C, D, DPI, PDPI, BC, ABC, A1, A2, A3, A4, A5, A6;

complex<double> MatrixDeterminant(complex<double>** A) { //Calculating the determinant of 3x3 matrix.
	complex<double> detA = 0;

	for (j = 0; j <= 2; j++) {
		for (k = 0; k <= 2; k++) {
			for (l = 0; l <= 2; l++) {
				detA += LeviCivita(j, k, l) * A[0][j] * A[1][k] * A[2][l];
			}
		}
	}
	return detA;
}

complex<double> MatrixTrace(complex<double>** A) {
	complex<double> Tr = 0;

	for (j = 0; j <= 2; ++j) {
		Tr += A[j][j];
	}

	return Tr;
}

complex<double>** MatrixProduct(double** A, complex<double>** B) {
	complex<double>** C = new complex<double>*[3];
	for (j = 0; j <= 2; ++j) {
		C[j] = new complex<double>[3];
	}

	complex<double> sum;

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
complex<double>** MatrixProduct(complex<double>** A, double** B) {
	complex<double>** C = new complex<double>*[3];
	for (j = 0; j <= 2; ++j) {
		C[j] = new complex<double>[3];
	}

	complex<double> sum;

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
complex<double>** MatrixProduct(complex<double>** A, complex<double>** B) {
	complex<double>** C = new complex<double>*[3];
	for (j = 0; j <= 2; ++j) {
		C[j] = new complex<double>[3];
	}

	complex<double> sum;

	for (j = 0; j <= 2; ++j) {
		for (k = 0; k <= 2; ++k) {
			sum = 0;
			for (l = 0; l <= 2; ++l) {
				sum = sum + A[j][l] * B[l][k];
			}
			C[j][k] = sum;
		}

	}

	return C;
}
complex<double>** MatrixProduct(complex<double>** A, complex<double>** B, complex<double>** C) {
	complex<double>** BC = new complex<double>*[3];
	for (j = 0; j <= 2; ++j) {
		BC[j] = new complex<double>[3];
	}
	
	complex<double>** ABC = new complex<double>*[3];
	for (j = 0; j <= 2; ++j) {
		ABC[j] = new complex<double>[3];
	}

	complex<double> sum;

	for (j = 0; j <= 2; ++j) {
		for (k = 0; k <= 2; ++k) {
			sum = 0;
			for (l = 0; l <= 2; ++l) {
				sum = sum + B[j][l] * C[l][k];
			}
			BC[j][k] = sum;
		}

	}

	for (j = 0; j <= 2; ++j) {
		for (k = 0; k <= 2; ++k) {
			sum = 0;
			for (l = 0; l <= 2; ++l) {
				sum = sum + A[j][l] * BC[l][k];
			}
			ABC[j][k] = sum;
		}

	}

	for (j = 0; j <= 2; ++j) {
		delete[] BC[j];
	}
	delete[] BC;

	return ABC;
}

complex<double>** MatrixSum(complex<double>** A1, complex<double>** A2) {
	complex<double>** C = new complex<double>*[3];
	for (j = 0; j <= 2; ++j) {
		C[j] = new complex<double>[3];
	}

	for (j = 0; j <= 2; ++j) {
		for (k = 0; k <= 2; ++k) {
			C[j][k] = A1[j][k] + A2[j][k];
		}
	}

	return C;
}
complex<double>** MatrixSum(complex<double>** A1, complex<double>** A2, complex<double>** A3, complex<double>** A4, complex<double>** A5, complex<double>** A6){
	complex<double>** C = new complex<double>*[3];
	for (j = 0; j <= 2; ++j) {
		C[j] = new complex<double>[3];
	}

	for (j = 0; j <= 2; ++j) {
		for (k = 0; k <= 2; ++k) {
			C[j][k] = A1[j][k] + A2[j][k] + A3[j][k] + A4[j][k] + A5[j][k] + A6[j][k];
		}
	}

	return C;
}

complex<double>** MatrixDifference(complex<double>** A1, complex<double>** A2) {
	complex<double>** C = new complex<double>*[3];
	for (j = 0; j <= 2; ++j) {
		C[j] = new complex<double>[3];
	}

	for (j = 0; j <= 2; ++j) {
		for (k = 0; k <= 2; ++k) {
			C[j][k] = A1[j][k] - A2[j][k];
		}
	}

	return C;
}

complex<double>** MatrixConjugate(complex<double>** A) {
	complex<double>** B = new complex<double>*[3];
	for (j = 0; j <= 2; ++j) {
		B[j] = new complex<double>[3];
	}
	for (j = 0; j <= 2; ++j) {
		for (k = 0; k <= 2; ++k) {
			B[k][j] = conj(A[j][k]);
		}
	}

	return B;
}

complex<double>** SU3_Rand_Gnrtr() { //Randomly generates an SU(3) matrix.
	double* ab = GramSchmidt(3);

	double* a2 = new double[3];
	double* b = new double[3];
	for (j = 0; j <= 2; ++j) {
		a2[j] = ab[j];
		b[j] = ab[j + 3];
	}

	delete[] ab;

	double* c = VecCross3D(a2, b);

	double** P = new double*[3];
	for (j = 0; j <= 2; ++j) {
		P[j] = new double[3];
	}
	for (j = 0; j <= 2; ++j) {
		P[0][j] = a2[j];
		P[1][j] = b[j];
		P[2][j] = c[j];
	}

	delete[] a2;
	delete[] b;
	delete[] c;

	double** PI = MatrixInverse(P);

	complex<double>** D = new complex<double>*[3];
	for (j = 0; j <= 2; ++j) {
		D[j] = new complex<double>[3];
		for (k = 0; k <= 2; ++k) {
			D[j][k] = 0;
		}
	}
	double alpha1 = SpecRandf(0., Pi);
	double alpha2 = SpecRandf(0., Pi);
	double alpha3 = 2 * Pi - alpha1 - alpha2;
	//	    z = cos(alpha)+1i*sin(alpha);

	complex<double> z1 = exp(1i*alpha1);
	D[0][0] = z1 / abs(z1);
	complex<double> z2 = exp(1i*alpha2);
	D[1][1] = z2 / abs(z2);
	complex<double> z3 = exp(1i*alpha3);
	D[2][2] = z3 / abs(z3);

	complex<double>** DPI = MatrixProduct(D, PI);

	complex<double>** PDPI = MatrixProduct(P, DPI);

	for (j = 0; j <= 2; ++j)
	{
		delete[] P[j];
		delete[] PI[j];
		delete[] D[j];
		delete[] DPI[j];
	}
	delete[]P;
	delete[]PI;
	delete[]D;
	delete[]DPI;
	
	return PDPI;
}