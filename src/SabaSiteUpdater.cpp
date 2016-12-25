#include <complex>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>

#include "../head/SabaComplexMatrix.h"
#include "../head/SabaMatrix.h"
#include "../head/SabaRandoms.h"
#include "../head/SabaVector.h"
#include "../head/Parameters.h"


using namespace std;

static int j, k, n, ni, n1, n2, n3, n4;
static int* IC, siteindex, return_array;
static double f;
static double PhiIndex, PhiIndex1, PhiIndex2, PhiIndex3, PhiIndex4, PhiIndex5, PhiIndex6, PhiIndex7, PhiIndex8;
static double RanPhi, S_init, S_final, DS;
static double x, eta;

double PhiUpdater(double f, double PhiIndex, double PhiIndex1, double PhiIndex2, double PhiIndex3, double PhiIndex4, double PhiIndex5, double PhiIndex6, double PhiIndex7, double PhiIndex8) {
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//	Eylemin degisen parcasi asagidaki gibidir:																													//
	//																																								//
	//																																								//
	//	  a^2	 							 a^2 f								   a^2								    a^2											//
	//	 _____ { PhiIndex - PhiIndex4 }^2 + _______ { PhiIndex - PhiIndex1 }^2 + _______ { PhiIndex - PhiIndex2 }^2 + _______ { PhiIndex - PhiIndex3 }^2			//
	//	  2 f								   2									2									 2											//
	//																																								//
	//	  1					   1																																	//
	// + ___ m^2 PhiIndex^2 + ___ lambda PhiIndex^4																													//
	//	  2					   4																																	//
	//																																								//
	//	  a^2	 							 a^2 f								   a^2								    a^2											//
	// + _____ { PhiIndex8 - PhiIndex }^2 + _______ { PhiIndex5 - PhiIndex }^2 + _______ { PhiIndex6 - PhiIndex }^2 + _______ { PhiIndex7 - PhiIndex }^2			//
	//	  2 f								   2									2									 2											//
	//																																								//
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	double S_init = (a * a * (PhiIndex - PhiIndex4) * (PhiIndex - PhiIndex4)) / (2 * f) + (a * a * f * (PhiIndex - PhiIndex1) * (PhiIndex - PhiIndex1)) / (2)
		+ (a * a * (PhiIndex - PhiIndex2) * (PhiIndex - PhiIndex2)) / (2) + (a * a * (PhiIndex - PhiIndex3) * (PhiIndex - PhiIndex3)) / (2)
		+ (m2 * PhiIndex * PhiIndex) / (2) + (lambda * PhiIndex * PhiIndex * PhiIndex * PhiIndex) / (4)
		+ (a * a * (PhiIndex8 - PhiIndex) * (PhiIndex8 - PhiIndex)) / (2 * f) + (a * a * f * (PhiIndex5 - PhiIndex) * (PhiIndex5 - PhiIndex)) / (2)
		+ (a * a * (PhiIndex6 - PhiIndex) * (PhiIndex6 - PhiIndex)) / (2) + (a * a * (PhiIndex7 - PhiIndex) * (PhiIndex7 - PhiIndex)) / (2);

	double RanPhi = Randf();

	double S_final = (a * a * (RanPhi - PhiIndex4) * (RanPhi - PhiIndex4)) / (2 * f) + (a * a * f * (RanPhi - PhiIndex1) * (RanPhi - PhiIndex1)) / (2)
		+ (a * a * (RanPhi - PhiIndex2) * (RanPhi - PhiIndex2)) / (2) + (a * a * (RanPhi - PhiIndex3) * (RanPhi - PhiIndex3)) / (2)
		+ (m2 * RanPhi * RanPhi) / (2) + (lambda * RanPhi * RanPhi * RanPhi * RanPhi) / (4)
		+ (a * a * (PhiIndex8 - RanPhi) * (PhiIndex8 - RanPhi)) / (2 * f) + (a * a * f * (PhiIndex5 - RanPhi) * (PhiIndex5 - RanPhi)) / (2)
		+ (a * a * (PhiIndex6 - RanPhi) * (PhiIndex6 - RanPhi)) / (2) + (a * a * (PhiIndex7 - RanPhi) * (PhiIndex7 - RanPhi)) / (2);

	double DS = S_final - S_init;
	double x, eta;

	if (DS < 0) {
		x = RanPhi;
	}
	else {
		eta = Randf();
		if (exp(- DS) >= eta) {
			x = RanPhi;
		}
		else {
			x = PhiIndex;
		}
	}
	return x;
}

double* DebugPhiUpdater(double f, double PhiIndex, double PhiIndex1, double PhiIndex2, double PhiIndex3, double PhiIndex4, double PhiIndex5, double PhiIndex6, double PhiIndex7, double PhiIndex8) {
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//	Eylemin degisen parcasi asagidaki gibidir:																													//
	//																																								//
	//																																								//
	//	  a^2	 							 a^2 f								   a^2								    a^2											//
	//	 _____ { PhiIndex - PhiIndex4 }^2 + _______ { PhiIndex - PhiIndex1 }^2 + _______ { PhiIndex - PhiIndex2 }^2 + _______ { PhiIndex - PhiIndex3 }^2			//
	//	  2 f								   2									2									 2											//
	//																																								//
	//	  1					   1																																	//
	// + ___ m^2 PhiIndex^2 + ___ lambda PhiIndex^4																													//
	//	  2					   4																																	//
	//																																								//
	//	  a^2	 							 a^2 f								   a^2								    a^2											//
	// + _____ { PhiIndex8 - PhiIndex }^2 + _______ { PhiIndex5 - PhiIndex }^2 + _______ { PhiIndex6 - PhiIndex }^2 + _______ { PhiIndex7 - PhiIndex }^2			//
	//	  2 f								   2									2									 2											//
	//																																								//
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	double S_init = (a * a * (PhiIndex - PhiIndex4) * (PhiIndex - PhiIndex4)) / (2 * f) + (a * a * f * (PhiIndex - PhiIndex1) * (PhiIndex - PhiIndex1)) / (2)
		+ (a * a * (PhiIndex - PhiIndex2) * (PhiIndex - PhiIndex2)) / (2) + (a * a * (PhiIndex - PhiIndex3) * (PhiIndex - PhiIndex3)) / (2)
		+ (m2 * PhiIndex * PhiIndex) / (2) + (lambda * PhiIndex * PhiIndex * PhiIndex * PhiIndex) / (4)
		+ (a * a * (PhiIndex8 - PhiIndex) * (PhiIndex8 - PhiIndex)) / (2 * f) + (a * a * f * (PhiIndex5 - PhiIndex) * (PhiIndex5 - PhiIndex)) / (2)
		+ (a * a * (PhiIndex6 - PhiIndex) * (PhiIndex6 - PhiIndex)) / (2) + (a * a * (PhiIndex7 - PhiIndex) * (PhiIndex7 - PhiIndex)) / (2);

	double RanPhi = Randf();

	double S_final = (a * a * (RanPhi - PhiIndex4) * (RanPhi - PhiIndex4)) / (2 * f) + (a * a * f * (RanPhi - PhiIndex1) * (RanPhi - PhiIndex1)) / (2)
		+ (a * a * (RanPhi - PhiIndex2) * (RanPhi - PhiIndex2)) / (2) + (a * a * (RanPhi - PhiIndex3) * (RanPhi - PhiIndex3)) / (2)
		+ (m2 * RanPhi * RanPhi) / (2) + (lambda * RanPhi * RanPhi * RanPhi * RanPhi) / (4)
		+ (a * a * (PhiIndex8 - RanPhi) * (PhiIndex8 - RanPhi)) / (2 * f) + (a * a * f * (PhiIndex5 - RanPhi) * (PhiIndex5 - RanPhi)) / (2)
		+ (a * a * (PhiIndex6 - RanPhi) * (PhiIndex6 - RanPhi)) / (2) + (a * a * (PhiIndex7 - RanPhi) * (PhiIndex7 - RanPhi)) / (2);

	double DS = S_final - S_init;
	double x, eta;

	if (DS < 0) {
		x = RanPhi;
	}
	else {
		eta = Randf();
		if (exp(- DS) >= eta) {
			x = RanPhi;
		}
		else {
			x = PhiIndex;
		}
	}
	double* return_array = new double[8];
	return_array[0] = x;
	return_array[1] = PhiIndex;
	return_array[2] = S_init;
	return_array[3] = RanPhi;
	return_array[4] = S_final;
	return_array[5] = DS;
	return_array[6] = eta;
	return_array[7] = exp(- DS);

	return return_array;
}

int* IndexConfFinder(int n) {
	//cout << "This is subroutine IndexConfFinder" << endl;
	int ni = 0;
	int* IC = new int[stdim];
	for (n1 = 0; n1 <= N_r - 1; ++n1) {
		for (n2 = 0; n2 <= N_y - 1; ++n2) {
			for (n3 = 0; n3 <= N_z - 1; ++n3) {
				for (n4 = 0; n4 <= N_t - 1; ++n4) {
					if (ni == n) {
						break;
					}
					++ni;
					//cout << ni << endl;
				}
				if (ni == n) {

					break;
				}
			}
			if (ni == n) {
				break;
			}
		}
		if (ni == n) {
			break;
		}
	}
	if (ni == n) {
		IC[0] = n1;
		IC[1] = n2;
		IC[2] = n3;
		IC[3] = n4;
		if (IC[3] == N_t) {
			IC[3] = 0;
			IC[2] += 1;
		}
		if (IC[2] == N_z) {
			IC[2] = 0;
			IC[1] += 1;
		}
		if (IC[1] == N_y) {
			IC[1] = 0;
			IC[0] += 1;
		}
		//cout << "subr " << "ni = " << ni << " n = " << n << " n1 = " << IC[0] << n1 << " n2 = " << IC[1] << n2 << " n3 = " << IC[2] << n3 << " n4 = " << IC[3] << n4 << endl;
	}

	//cout << "We are exiting the IndexConfFinder now." << endl;
	return IC;
}

//int IndexFinder(int* IC) {
//	int n = IC[0] * N_y * N_z * N_t + IC[1] * N_z * N_t + IC[2] * N_t + IC[3];
//	return n;
//}

int IndexFinder(int* IC) {
	//cout << "This is subroutine IndexFinder" << endl;
	int ni = 0;
	int n = 0;
	for (n1 = 0; n1 <= N_r - 1; ++n1) {
		for (n2 = 0; n2 <= N_y - 1; ++n2) {
			for (n3 = 0; n3 <= N_z - 1; ++n3) {
				for (n4 = 0; n4 <= N_t - 1; ++n4) {
					if (n1 == IC[0] && n2 == IC[1] && n3 == IC[2] && n4 == IC[3]) {
						n = ni;
						break;
					}
					++ni;
				}
				if (n1 == IC[0] && n2 == IC[1] && n3 == IC[2] && n4 == IC[3]) {
					break;
				}
			}
			if (n1 == IC[0] && n2 == IC[1] && n3 == IC[2] && n4 == IC[3]) {
				break;
			}
		}
		if (n1 == IC[0] && n2 == IC[1] && n3 == IC[2] && n4 == IC[3]) {
			break;
		}
	}
	//cout << "We are exiting the IndexFinder now." << endl;
	return n;
}

int* BC(int* siteindex) {
	//siteindex[0] = (siteindex[0] + stdim) % N_r;
	//siteindex[1] = (siteindex[1] + stdim) % N_y;
	//siteindex[2] = (siteindex[2] + stdim) % N_z;
	//siteindex[3] = (siteindex[3] + stdim) % N_t;
	if (siteindex[0] == N_r) { siteindex[0] = 0; }
	if (siteindex[0] == -1) { siteindex[0] = N_r - 1; }
	if (siteindex[1] == N_y) { siteindex[1] = 0; }
	if (siteindex[1] == -1) { siteindex[1] = N_y - 1; }
	if (siteindex[2] == N_z) { siteindex[2] = 0; }
	if (siteindex[2] == -1) { siteindex[2] = N_z - 1; }
	if (siteindex[3] == N_t) { siteindex[3] = 0; }
	if (siteindex[3] == -1) { siteindex[3] = N_t - 1; }
	return siteindex;
}
