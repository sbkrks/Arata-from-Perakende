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
#include "Parameters.h"

using namespace std;

double PhiUpdater(double f, double PhiIndex, double PhiIndex1, double PhiIndex2, double PhiIndex3, double PhiIndex4, double PhiIndex5, double PhiIndex6, double PhiIndex7, double PhiIndex8);

double* DebugPhiUpdater(double f, double PhiIndex, double PhiIndex1, double PhiIndex2, double PhiIndex3, double PhiIndex4, double PhiIndex5, double PhiIndex6, double PhiIndex7, double PhiIndex8);

int* IndexConfFinder(int n);

int IndexFinder(int* IC);

int* BC(int* siteindex);
