#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <complex>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>

#include "SabaComplexMatrix.h"
#include "SabaMatrix.h"
#include "SabaRandoms.h"
#include "SabaVector.h"
#include "SabaSiteUpdater.h"
#include "JackKnife.h"
#include "Parameters.h"

//////////////////////////////////////////////////////////
// arxiv.1603.00716v1									//
// Sanjin Benic & Arata Yamamoto 03/03/2016				//
// Quantum Monte Carlo Simulations with A Black Hole	//
//////////////////////////////////////////////////////////


using namespace std;

int n1, n2, n3, n4;		// r, y, z and t indices of a site respectively
int mcs;				// # of Monte Carlo Step index
int row, col;			// several indices for initializing and updating the lattice
int j, k, l, p, q;


int main() {

		ofstream InitField("RandomFile.txt");
        double x = Randf();
		InitField << x << endl;


	return 0;
}

//clock_t start = clock();
// bla bla
//clock_t end = clock();
//float seconds = (float)(end - start) / CLOCKS_PER_SEC;
//cout << "time during bla bla = " << seconds << endl;
