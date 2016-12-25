#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <complex>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>
#include "../head/Parameters.h"

static int j;
static double*OlcumBar;
static double total;
static double OlcumBarNokta;
static double variance;
static double StaDev;

double JackKnife(double* Olcum) {
	int j;
	double*OlcumBar = new double[N_conf];
	double total = 0;
	double OlcumBarNokta = 0;
	double variance = 0;
	double StaDev;

	for (j = 0; j <= N_conf - 1; ++j) {
		total += Olcum[j];
	}

	for (j = 0; j <= N_conf - 1; ++j) {
		OlcumBar[j] = (total - Olcum[j]) / (N_conf - 1);
	}

	for (j = 0; j <= N_conf - 1; ++j) {
		OlcumBarNokta += OlcumBar[j];
	}
	OlcumBarNokta = OlcumBarNokta / N_conf;

	for (j = 0; j <= N_conf - 1; ++j) {
		variance += (OlcumBar[j] - OlcumBarNokta) * (OlcumBar[j] - OlcumBarNokta);
	}
	variance = (N_conf - 1) * variance / N_conf;

	StaDev = sqrt(variance);

	return StaDev;
}