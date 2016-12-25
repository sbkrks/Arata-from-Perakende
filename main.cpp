#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <complex>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>

//#include "SabaComplexMatrix.h"
#include "head/SabaMatrix.h"
#include "head/SabaRandoms.h"
#include "head/SabaVector.h"
#include "head/SabaSiteUpdater.h"
#include "head/JackKnife.h"
#include "head/Parameters.h"

//////////////////////////////////////////////////////////
// arxiv.1603.00716v1									//
// Sanjin Benic & Arata Yamamoto 03/03/2016				//
// Quantum Monte Carlo Simulations With A Black Hole	//
//////////////////////////////////////////////////////////

using namespace std;

int n1, n2, n3, n4;		// r, y, z and t indices of a site respectively
int mcs;				// # of Monte Carlo Step index
int row, col;			// several indices for initializing and updating the lattice
int j, k, l, p, q;

//  int to string conversion
/** int to string conversion for file namings */
string itos(int i)
{
	stringstream s;
	s << i;
	return s.str();
}

//  double to string conversion
/** double to string conversion for file namings */
string dtos(double d)
{
	stringstream s;
	s << d;
	return s.str();
}

int main() {

	double f; // The coefficient of dt^2 and inverse coefficient of dr^2 in the invariant line element
			  ////////////////////////////////////////////////////////////////////
			  //							 1									//
			  //	ds ^ 2 = f(r) dt ^ 2 + _____ dr ^ 2 + dy ^ 2 + dz ^ 2		//
			  //							f(r)								//
			  //				R												//
			  //	f(r) = 1 - ___												//
			  //				r												//
			  ////////////////////////////////////////////////////////////////////

			  ////////////////////////////////////////////////////////////
			  //			Initializing the scalar field				//
			  ////////////////////////////////////////////////////////////////////////////////////////////
	double**** Phi = new double***[N_r];				// The real scalar field						//
														//////////////////////////////////////////////////////////
	double***** PhiMeas = new double****[N_r];			// Memory where we store measurements of scalar field.	//
														// Every last index due another measurement.			//
														//////////////////////////////////////////////////////////
	double**** MeanPhi = new double***[N_r];			// Mean value of the scalar function at sites.	//
														// < Phi >										//
														//////////////////////////////////////////////////////////////////
	double**** MeanPhi2 = new double***[N_r];			// Mean value of the square of the scalar function at sites.	//
														// < Phi * Phi >												//
														//////////////////////////////////////////////////////////////////
	double**** G_0 = new double***[N_r];				// Two-point function with |x - x'| = 0			//
														//	G(x, x') = < Phi(x) * Phi(x') >				//
														//	G(0) = < Phi^2 > = < Phi[i] * Phi[i] >		//
														//////////////////////////////////////////////////////////
	double**** G_inf = new double***[N_r];				// Two-point function with |x - x'| = infinity			//
															//	G(x, x') = < Phi(x) * Phi(x') >						//
															//	G(inf) = < Phi >^2 = <Phi[i]> * <Phi[i']>			//
															//////////////////////////////////////////////////////////
	double**** C1 = new double***[N_r];						// Condansate ratio				//
															//		 G(inf)					//
															//	C = ________				//
															//		  G(0)					//
															//////////////////////////////////

	double**** StandartDev = new double***[N_r];

	for (n1 = 0; n1 <= N_r - 1; ++n1) {
		Phi[n1] = new double**[N_y];
		PhiMeas[n1] = new double***[N_y];
		G_0[n1] = new double**[N_y];
		G_inf[n1] = new double**[N_y];
		C1[n1] = new double**[N_y];
		MeanPhi[n1] = new double**[N_y];
		MeanPhi2[n1] = new double**[N_y];
		StandartDev[n1] = new double**[N_y];
		for (n2 = 0; n2 <= N_y - 1; ++n2) {
			Phi[n1][n2] = new double*[N_z];
			PhiMeas[n1][n2] = new double**[N_z];
			G_0[n1][n2] = new double*[N_z];
			G_inf[n1][n2] = new double*[N_z];
			C1[n1][n2] = new double*[N_z];
			MeanPhi[n1][n2] = new double*[N_z];
			MeanPhi2[n1][n2] = new double*[N_z];
			StandartDev[n1][n2] = new double*[N_z];
			for (n3 = 0; n3 <= N_z - 1; ++n3) {
				Phi[n1][n2][n3] = new double[N_t];
				PhiMeas[n1][n2][n3] = new double*[N_t];
				G_0[n1][n2][n3] = new double[N_t];
				G_inf[n1][n2][n3] = new double[N_t];
				C1[n1][n2][n3] = new double[N_t];
				MeanPhi[n1][n2][n3] = new double[N_t];
				MeanPhi2[n1][n2][n3] = new double[N_t];
				StandartDev[n1][n2][n3] = new double[N_t];
				for (n4 = 0; n4 <= N_t - 1; ++n4) {
					PhiMeas[n1][n2][n3][n4] = new double[N_conf];
				}
			}
		}
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////
	// Reading initial values of the scalar field from a text file, just in case.						//
	//																									//
	// You need to create the "OInputField.txt" file for Scalar Field.									//
	// You need to create the "OInputGZero.txt" file for G_0.											//
	// You need to create the "OInputGInf.txt" file for G_inf.											//
	//																									//
	// DO NOT FORGET to close the randomly generation of scalar field below.							//
	// DO NOT FORGET to cancel the process of setting elements of G_0 and G_inf to zero below.			//
	//////////////////////////////////////////////////////////////////////////////////////////////////////
	//cout << "Update edilmis son lattice konfigurasyonu dosyadan okunuyor." << endl;
	//cout << "Hesaplanmis son G_0 elemanlari dosyadan okunuyor." << endl;
	//cout << "Hesaplanmis son G_inf elemanlari dosyadan okunuyor." << endl;
	//fstream InputField("OInputField.txt", ios_base::in);
	//fstream InputGzero("OInputGZero.txt", ios_base::in);
	//fstream InputGinf("OInputGInf.txt", ios_base::in);
	//for (j = 0; j <= N - 1; ++j) {
	//	InputField >> Phi[j];
	//	InputGzero >> G_0[j];
	//	InputGinf >> G_inf[j];
	//}
	//cout << "Okuma islemi tamamlandi." << endl;
	//InputField.close();
	//InputGzero.close();
	//InputGinf.close();

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Updating the lattice until equilibrium																														//
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//	Eylemin degisen parcasi asagidaki gibidir:																													//
	//																																								//
	//																																								//
	//	  a^2	 							  a^2 f								   a^2								    a^2											//
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

	int Index1, Index2, Index3, Index4, Index1m, Index2m, Index3m, Index4m, Index1p, Index2p, Index3p, Index4p;		// indices correspond to IndexConf1, IndexConf2, IndexConf3, IndexConf4, IndexConf5, IndexConf6, IndexConf7, IndexConf8 respectively
//	int aloop;			// index for loop of a
//	double alist[15];
//	double a, m2, R, epsilon;
//	alist[0] = 0.2;
//	alist[1] = 0.5;
//	alist[2] = 0.7;
//	alist[3] = 1.0;
//	alist[4] = 1.5;
//	alist[5] = 2.0;
//	alist[6] = 2.5;
//	alist[7] = 3.0;
//	alist[8] = 3.5;
//	alist[9] = 4.0;
//	alist[10] = 5.0;
//	alist[11] = 6.0;
//	alist[12] = 7.0;
//	alist[13] = 8.0;
//	alist[14] = 9.0;

//		cout << "loop #" << aloop + 1 << endl;
//		a = alist[aloop];				// lattice spacing
        cout << "a = " << a << endl;
		double m2 = -0.2 / (a * a);			// the square of the tachyonic mass
		double R = (N_t * a) / (4 * Pi);		// Schwartzschild radius of the blackhole
		double epsilon = 0.1 * a;				// a small number will be added to r so as to avoid the singularity at r = R

		//char date[9];
        //_strdate_s(date);

//        //string OutputFolder("./Output for a = " + dtos(a) + "/");
//        string InitFieldFolderName("OInitialScalarField_a" + dtos(a) + "_Nconf" + itos(N_conf) + ".txt");
////		cout << InitFieldFolderName << endl;
//		ofstream InitField(InitFieldFolderName);
//		InitField << "date" << "\nCaution!\n This is the one we record the measurements seperately.\nWe measure the mean value of Phi and Phi^2.\nAnd we record Phi during measurements.\n" << endl;
//		InitField << "# of Measurements = " << N_conf << "\t# of inbetween sweeps = " << N_equi << "\n" << endl;
//		InitField << "\na = " << a << "\nm^2 = " << m2 << "\nepsilon = " << epsilon << "\n" << endl;

		string EquiFolderName("OEquilibria_a" + dtos(a) + "_Nconf" + itos(N_conf) + ".txt");
		ofstream Equi(EquiFolderName);
		Equi << "date" << "\nCaution!\n This is the one we record the measurements seperately.\nWe measure the mean value of Phi and Phi^2.\nAnd we record Phi during measurements.\n" << endl;
		Equi << "# of Measurements = " << N_conf << "\t# of inbetween sweeps = " << N_equi << "\n" << endl;
		Equi << "\na = " << a << "\nm^2 = " << m2 << "\nepsilon = " << epsilon << "\n" << endl;

		string CondansateFolderName("OCondensate_a" + dtos(a) + "_Nconf" + itos(N_conf) + ".txt");
		ofstream Cndnst1(CondansateFolderName);
		Cndnst1 << "date" << "\nCaution!\n This is the one we record the measurements seperately.\nWe measure the mean value of Phi and Phi^2.\nAnd we record Phi during measurements.\n" << endl;
		Cndnst1 << "# of Measurements = " << N_conf << "\t# of inbetween sweeps = " << N_equi << "\n" << endl;
		Cndnst1 << "\na = " << a << "\nm^2 = " << m2 << "\nepsilon = " << epsilon << "\n" << endl;
		Cndnst1 << "C = C_inf / C_0" << "\nC_inf[i] = <Phi[i]> * <Phi[i']>" << endl;

//		string FieldMeasFolderName("OFieldMeasurements_a" + dtos(a) + "_Nconf" + itos(N_conf) + ".txt");
//		ofstream FieldMeas(FieldMeasFolderName);
//		FieldMeas << "date" << "\nCaution!\n This is the one we record the measurements seperately.\nWe measure the mean value of Phi and Phi^2.\nAnd we record Phi during measurements.\n" << endl;
//		FieldMeas << "# of Measurements = " << N_conf << "\t# of inbetween sweeps = " << N_equi << "\n" << endl;
//		FieldMeas << "\na = " << a << "\nm^2 = " << m2 << "\nepsilon = " << epsilon << "\n" << endl;
//		FieldMeas << "C = C_inf / C_0" << "\nC_inf[i] = <Phi[i]> * <Phi[i']>" << endl;

		//string EquilibriaFolderName("OEquilibria_a" + dtos(a) + ".txt");
		//ofstream Equi(EquilibriaFolderName);
		//Equi << date << "\nEquilibria data for a = " << a << endl;
		//Equi << "\nCaution! This is the one we don't record the measurements seperately.\nWe measure the mean value of Phi and Phi^2.\nThen we calculate the physical quantities.\n" << endl;
		//Equi << "# of Measurements = " << N_conf << "\t# of inbetween sweeps = " << N_equi << "\n" << endl;
		for (n1 = 0; n1 <= N_r - 1; ++n1) {
			for (n2 = 0; n2 <= N_y - 1; ++n2) {
				for (n3 = 0; n3 <= N_z - 1; ++n3) {
					for (n4 = 0; n4 <= N_t - 1; ++n4) {
						Phi[n1][n2][n3][n4] = Randf();
						//InitField << n1 << "\t" << n2 << "\t" << n3 << "\t" << n4 << "\t" << Phi[n1][n2][n3][n4] << endl;
						C1[n1][n2][n3][n4] = 0;
						G_0[n1][n2][n3][n4] = 0;
						G_inf[n1][n2][n3][n4] = 0;
						MeanPhi[n1][n2][n3][n4] = 0;
						MeanPhi2[n1][n2][n3][n4] = 0;
						StandartDev[n1][n2][n3][n4] = 0;
						for (k = 0; k <= N_conf; ++k) {
							PhiMeas[n1][n2][n3][n4][k] = 0;
						}
					}
				}
			}
		}

		///////////////////////////////////////////////////////
		//		Beginning of The Warm Up Monte Carlo		 //
		///////////////////////////////////////////////////////
		cout << "Warm-Up Monte Carlo basladi." << endl;
		for (mcs = 1; mcs <= 10 * N_equi; ++mcs) {
			for (j = 1; j <= N; ++j) {
				Index1 = Randi(N_r); Index1m = Index1 - 1; Index1p = Index1 + 1;
				Index2 = Randi(N_y); Index2m = Index2 - 1; Index2p = Index2 + 1;
				Index3 = Randi(N_z); Index3m = Index3 - 1; Index3p = Index3 + 1;
				Index4 = Randi(N_t); Index4m = Index4 - 1; Index4p = Index4 + 1;

				if (Index1m == -1) { Index1m = N_r - 1; }
				if (Index2m == -1) { Index2m = N_y - 1; }
				if (Index3m == -1) { Index3m = N_z - 1; }
				if (Index4m == -1) { Index4m = N_t - 1; }

				if (Index1p == N_r) { Index1p = 0; }
				if (Index2p == N_y) { Index2p = 0; }
				if (Index3p == N_z) { Index3p = 0; }
				if (Index4p == N_t) { Index4p = 0; }

				double r = (Index1 * a + R + epsilon);
				f = 1 - R / r;

				Phi[Index1][Index2][Index3][Index4] = PhiUpdater(m2, f, Phi[Index1][Index2][Index3][Index4], Phi[Index1m][Index2][Index3][Index4], Phi[Index1][Index2m][Index3][Index4], Phi[Index1][Index2][Index3m][Index4], Phi[Index1][Index2][Index3][Index4m], Phi[Index1p][Index2][Index3][Index4], Phi[Index1][Index2p][Index3][Index4], Phi[Index1][Index2][Index3p][Index4], Phi[Index1][Index2][Index3][Index4p]);

			}
			//////////////////////////////////////////////
			//		Equilibria Data Part				//
			//////////////////////////////////////////////
			double Stotal = 0;
			for (n1 = 0; n1 <= N_r - 1; ++n1) {
				for (n2 = 0; n2 <= N_y - 1; ++n2) {
					for (n3 = 0; n3 <= N_z - 1; ++n3) {
						for (n4 = 0; n4 <= N_t - 1; ++n4) {
							Index1 = n1; Index1m = Index1 - 1; Index1p = Index1 + 1;
							Index2 = n2; Index2m = Index2 - 1; Index2p = Index2 + 1;
							Index3 = n3; Index3m = Index3 - 1; Index3p = Index3 + 1;
							Index4 = n4; Index4m = Index4 - 1; Index4p = Index4 + 1;

							if (Index1m == -1) { Index1m = N_r - 1; }
							if (Index2m == -1) { Index2m = N_y - 1; }
							if (Index3m == -1) { Index3m = N_z - 1; }
							if (Index4m == -1) { Index4m = N_t - 1; }

							if (Index1p == N_r) { Index1p = 0; }
							if (Index2p == N_y) { Index2p = 0; }
							if (Index3p == N_z) { Index3p = 0; }
							if (Index4p == N_t) { Index4p = 0; }

							double S = (a * a * (Phi[Index1][Index2][Index3][Index4] - Phi[Index1][Index2][Index3][Index4m]) * (Phi[Index1][Index2][Index3][Index4] - Phi[Index1][Index2][Index3][Index4m])) / (2 * f)
								+ (a * a * f * (Phi[Index1][Index2][Index3][Index4] - Phi[Index1m][Index2][Index3][Index4]) * (Phi[Index1][Index2][Index3][Index4] - Phi[Index1m][Index2][Index3][Index4])) / (2)
								+ (a * a * (Phi[Index1][Index2][Index3][Index4] - Phi[Index1][Index2m][Index3][Index4]) *(Phi[Index1][Index2][Index3][Index4] - Phi[Index1][Index2m][Index3][Index4])) / (2)
								+ (a * a * (Phi[Index1][Index2][Index3][Index4] - Phi[Index1][Index2][Index3m][Index4]) *(Phi[Index1][Index2][Index3][Index4] - Phi[Index1][Index2][Index3m][Index4])) / (2)
								+ (m2 * Phi[Index1][Index2][Index3][Index4] * Phi[Index1][Index2][Index3][Index4]) / (2) + (lambda * Phi[Index1][Index2][Index3][Index4] * Phi[Index1][Index2][Index3][Index4] * Phi[Index1][Index2][Index3][Index4] * Phi[Index1][Index2][Index3][Index4]) / (4);
							Stotal += S;
						}
					}
				}
			}
			Equi << Stotal << endl;
		}
		cout << "Warm-Up Monte Carlo bitti." << endl;

		cout << "Olcum basladi." << endl;
		for (k = 0; k <= N_conf - 1; ++k) {
			if ((k + 1) % 100 == 0) { cout << "olcum no = " << k + 1 << endl; }
			for (mcs = 1; mcs <= N_equi; ++mcs) {
				for (j = 1; j <= N; ++j) {
					Index1 = Randi(N_r); Index1m = Index1 - 1; Index1p = Index1 + 1;
					Index2 = Randi(N_y); Index2m = Index2 - 1; Index2p = Index2 + 1;
					Index3 = Randi(N_z); Index3m = Index3 - 1; Index3p = Index3 + 1;
					Index4 = Randi(N_t); Index4m = Index4 - 1; Index4p = Index4 + 1;

					if (Index1m == -1) { Index1m = N_r - 1; }
					if (Index2m == -1) { Index2m = N_y - 1; }
					if (Index3m == -1) { Index3m = N_z - 1; }
					if (Index4m == -1) { Index4m = N_t - 1; }

					if (Index1p == N_r) { Index1p = 0; }
					if (Index2p == N_y) { Index2p = 0; }
					if (Index3p == N_z) { Index3p = 0; }
					if (Index4p == N_t) { Index4p = 0; }

					f = 1 - R / (Index1 * a + R + epsilon);

					Phi[Index1][Index2][Index3][Index4] = PhiUpdater(m2, f, Phi[Index1][Index2][Index3][Index4], Phi[Index1m][Index2][Index3][Index4], Phi[Index1][Index2m][Index3][Index4], Phi[Index1][Index2][Index3m][Index4], Phi[Index1][Index2][Index3][Index4m], Phi[Index1p][Index2][Index3][Index4], Phi[Index1][Index2p][Index3][Index4], Phi[Index1][Index2][Index3p][Index4], Phi[Index1][Index2][Index3][Index4p]);

				}
			}

			for (n1 = 0; n1 <= N_r - 1; ++n1) {
				for (n2 = 0; n2 <= N_y - 1; ++n2) {
					for (n3 = 0; n3 <= N_z - 1; ++n3) {
						for (n4 = 0; n4 <= N_t - 1; ++n4) {
							PhiMeas[n1][n2][n3][n4][k] = Phi[n1][n2][n3][n4];
						}
					}
				}
			}

		}
		cout << "Olcum bitti." << endl;

		cout << "Fiziksel nicelikler hesaplaniyor..." << endl;
		double A;
		double A2;
		for (n1 = 0; n1 <= N_r - 1; ++n1) {
			for (n2 = 0; n2 <= N_y - 1; ++n2) {
				for (n3 = 0; n3 <= N_z - 1; ++n3) {
					for (n4 = 0; n4 <= N_t - 1; ++n4) {
						A = 0;
						A2 = 0;
						for (k = 0; k <= N_conf - 1; ++k) {
							A += PhiMeas[n1][n2][n3][n4][k];
							A2 += PhiMeas[n1][n2][n3][n4][k] * PhiMeas[n1][n2][n3][n4][k];
						}
						MeanPhi[n1][n2][n3][n4] = A / N_conf;
						MeanPhi2[n1][n2][n3][n4] = A2 / N_conf;
						G_0[n1][n2][n3][n4] = MeanPhi2[n1][n2][n3][n4];
					}
				}
			}
		}
		cout << "Fiziksel nicelikler hesaplandi." << endl;

		for (n1 = 0; n1 <= N_r - 1; ++n1) {
			for (n2 = 0; n2 <= N_y - 1; ++n2) {
				for (n3 = 0; n3 <= N_z - 1; ++n3) {
					for (n4 = 0; n4 <= N_t - 1; ++n4) {
						StandartDev[n1][n2][n3][n4] = JackKnife(PhiMeas[n1][n2][n3][n4]);
					}
				}
			}
		}

		for (n1 = 0; n1 <= N_r - 1; ++n1) {
			for (n2 = 0; n2 <= N_y - 1; ++n2) {
				for (n3 = 0; n3 <= N_z - 1; ++n3) {
					for (n4 = 0; n4 <= N_t - 1; ++n4) {
						int n3p = (n3 + (N_z / 2)) % N_z;
						G_inf[n1][n2][n3][n4] = MeanPhi[n1][n2][n3][n4] * MeanPhi[n1][n2][n3p][n4];
						C1[n1][n2][n3][n4] = G_inf[n1][n2][n3][n4] / G_0[n1][n2][n3][n4];
					}
				}
			}
		}

		cout << "Condansate1 yazdiriliyor..." << endl;
		Cndnst1 << "n1 \t n2 \t n3 \t n4 \t Condensate \t Standart Deviation" << endl;
		for (n2 = 0; n2 <= N_y - 1; ++n2) {
			for (n3 = 0; n3 <= N_z - 1; ++n3) {
				for (n4 = 0; n4 <= N_t - 1; ++n4) {
					for (n1 = 0; n1 <= N_r - 1; ++n1) {
						Cndnst1 << n1 << "\t" << n2 << "\t" << n3 << "\t" << n4 << "\t" << C1[n1][n2][n3][n4] << "\t" << StandartDev[n1][n2][n3][n4] << endl;
					}
				}
			}
		}
		cout << "Condansate1 yazimi tamamlandi." << endl;


		//	Time << "Measurements are being written. Time of finishing = " << time(NULL) << endl;
		//InitField.close();
		Cndnst1.close();
		Equi.close();
		//FieldMeas.close();

//
//	/*for (n1 = 0; n1 <= N_r - 1; ++n1) {
//		delete Phi[n1];
//		delete PhiMeas[n1];
//		delete G_0[n1];
//		delete G_inf[n1];
//		delete C1[n1];
//		delete MeanPhi[n1];
//		delete MeanPhi2[n1];
//		delete StandartDev[n1];
//	*/	//for (n2 = 0; n2 <= N_y - 1; ++n2) {
//		//	Phi[n1][n2] = new double*[N_z];
//		//	PhiMeas[n1][n2] = new double**[N_z];
//		//	G_0[n1][n2] = new double*[N_z];
//		//	G_inf[n1][n2] = new double*[N_z];
//		//	C1[n1][n2] = new double*[N_z];
//		//	MeanPhi[n1][n2] = new double*[N_z];
//		//	MeanPhi2[n1][n2] = new double*[N_z];
//		//	StandartDev[n1][n2] = new double*[N_z];
//		//	for (n3 = 0; n3 <= N_z - 1; ++n3) {
//		//		Phi[n1][n2][n3] = new double[N_t];
//		//		PhiMeas[n1][n2][n3] = new double*[N_t];
//		//		G_0[n1][n2][n3] = new double[N_t];
//		//		G_inf[n1][n2][n3] = new double[N_t];
//		//		C1[n1][n2][n3] = new double[N_t];
//		//		MeanPhi[n1][n2][n3] = new double[N_t];
//		//		MeanPhi2[n1][n2][n3] = new double[N_t];
//		//		StandartDev[n1][n2][n3] = new double[N_t];
//		//		for (n4 = 0; n4 <= N_t - 1; ++n4) {
//		//			PhiMeas[n1][n2][n3][n4] = new double[N_conf];
//		//		}
//		//	}
//		//}
	//}

	return 0;
}

//clock_t start = clock();
// bla bla
//clock_t end = clock();
//float seconds = (float)(end - start) / CLOCKS_PER_SEC;
//cout << "time during bla bla = " << seconds << endl;
