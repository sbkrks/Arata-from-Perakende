
#define Pi double( acos(-1.))

#define stdim 4			// dimension of spacetime
#define N_r 10			// # of sites in the direction of r
#define N_y 10			// # of sites in the direction of y
#define N_z 10			// # of sites in the direction of z
#define N_t 60			// # of sites in the ditection of t

//#define a 0.31622776601			// lattice spacing
//#define a 2.176					// lattice spacing by Huseyin Bahtiyar
#define a 9
#define m2 -0.2 / (a * a)			// the square of the tachyonic mass
#define R (N_t * a) / (4 * Pi)		// Schwartzschild radius of the blackhole
#define epsilon 0.1 * a			// a small number added to r so as to avoid the singularity at r = R
#define lambda 0.2					// coupling constant of phi^4 term in action

#define N (N_r * N_y * N_y * N_t)	// total # of sites
//#define N 10
#define N_equi 10					// # of Monte Carlo sweeps to obtain sufficiently independent configuration inbetween every two measurements.
#define N_conf 10000				// # of measurements
#define WU 2						// multiplier to obtain the # of Monte Carlo sweeps for warm up process.
									//We are gonna perform WU * N_equi sweeps for warm up process.

