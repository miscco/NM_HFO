/*
*	Copyright (c) 2014 Michael Schellenberger Costa
*
*	Permission is hereby granted, free of charge, to any person obtaining a copy
*	of this software and associated documentation files (the "Software"), to deal
*	in the Software without restriction, including without limitation the rights
*	to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
*	copies of the Software, and to permit persons to whom the Software is
*	furnished to do so, subject to the following conditions:
*
*	The aboV_p copyright notice and this permission notice shall be included in
*	all copies or substantial portions of the Software.
*
*	THE SOFTWARE IS PROV_sDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
*	IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
*	FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EV_pNT SHALL THE
*	AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
*	LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
*	OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
*	THE SOFTWARE.
*/

/****************************************************************************************************/
/*									Functions of the cortical module								*/
/****************************************************************************************************/
#include "C_Column.h"

/****************************************************************************************************/
/*										 Initialization of RNG 										*/
/****************************************************************************************************/
void C_Column::set_RNG(void) {
	/* Number of independent streams */
	int N = 6;

	/* Create RNG for each stream */
	for (int i=0; i<N; ++i){
		/* Add the RNG */
		MTRands.push_back({ENG(rand()), DIST (mphi, sqrt(dphi))});

		/* Get the random number for the first iteration */
		Rand_vars.push_back(MTRands[i]());
	}
}
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/


/****************************************************************************************************/
/*										 Firing Rate functions 										*/
/****************************************************************************************************/
/* Pyramidal firing rate */
double C_Column::get_Qp	(int N) const{
	double q = Qp_max / (1 + exp(-(Phi_pp[N] - theta_p) / sigma_p));
	return q;
}

/* Slow inhibitory firing rate */
double C_Column::get_Qs	(int N) const{
	double q = Qs_max / (1 + exp(-(Phi_sA[N] - theta_s) / sigma_s));
	return q;
}

/* Fast inhibitory firing rate */
double C_Column::get_Qf	(int N) const{
	double q = Qf_max / (1 + exp(-(Phi_fA[N] - theta_f) / sigma_f));
	return q;
}
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/


/****************************************************************************************************/
/*										 RK noise scaling 											*/
/****************************************************************************************************/
double C_Column::noise_xRK(int N, int M) const{
	extern const double h;
	extern const vector<double> B1, B2;
	double n = 1  / h * (B1[N] * Rand_vars[2*M] + B2[N] * Rand_vars[2*M+1]);
	return n;
}
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/


/****************************************************************************************************/
/*										Calculate the Nth SRK term									*/
/****************************************************************************************************/
void C_Column::set_RK (int N) {
	extern const double dt;
	int M = N+4;
	Phi_pp	[M] = dt*(x_pp[N]);
	Phi_ps	[M] = dt*(x_ps[N]);
	Phi_pf	[M] = dt*(x_pf[N]);
	Phi_sA	[M] = dt*(x_sA[N]);
	Phi_sB	[M] = dt*(x_sB[N]);
	Phi_fA	[M] = dt*(x_fA[N]);
	x_pp  	[M] = dt*(pow(gamma_p,  2) * (N_pp * get_Qp(N) + noise_xRK(N, 0)- Phi_pp[N]) - 2 * gamma_p  * x_pp[N]);
	x_ps  	[M] = dt*(pow(gamma_p,  2) * (N_ps * get_Qp(N) + noise_xRK(N, 1)- Phi_ps[N]) - 2 * gamma_p  * x_ps[N]);
	x_pf  	[M] = dt*(pow(gamma_p,  2) * (N_pf * get_Qp(N) + noise_xRK(N, 2)- Phi_pf[N]) - 2 * gamma_p  * x_pf[N]);
	x_sA  	[M] = dt*(pow(gamma_sA, 2) * (		 get_Qs(N)					- Phi_sA[N]) - 2 * gamma_sA * x_sA[N]);
	x_sB  	[M] = dt*(pow(gamma_sB, 2) * (		 get_Qs(N)					- Phi_sB[N]) - 2 * gamma_sB * x_sB[N]);
	x_fA  	[M] = dt*(pow(gamma_fA, 2) * (		 get_Qs(N)					- Phi_fA[N]) - 2 * gamma_fA * x_fA[N]);

	/* Calculate the intermediate RK terms */
	switch (N){
	case 0:
	case 1:
		set_intermediate_RK(N, 0.5);
		break;
	case 2:
		set_intermediate_RK(N, 1.);
		break;
	case 3:
		break;
	}
}
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/


/****************************************************************************************************/
/*									Update intermediate RK values   								*/
/****************************************************************************************************/
inline void C_Column::set_intermediate_RK (int N, double mult_RK) {
	int M = N+4;
	Phi_pp 	[N+1] = Phi_pp  [0] + mult_RK * Phi_pp  [M];
	Phi_ps 	[N+1] = Phi_ps  [0] + mult_RK * Phi_ps  [M];
	Phi_pf 	[N+1] = Phi_pf  [0] + mult_RK * Phi_pf  [M];
	Phi_sA 	[N+1] = Phi_sA  [0] + mult_RK * Phi_sA  [M];
	Phi_sB 	[N+1] = Phi_sB  [0] + mult_RK * Phi_sB  [M];
	Phi_fA 	[N+1] = Phi_fA  [0] + mult_RK * Phi_fA  [M];
	x_pp 	[N+1] = x_pp    [0] + mult_RK * x_ps    [M];
	x_ps 	[N+1] = x_ps    [0] + mult_RK * x_ps    [M];
	x_pf 	[N+1] = x_pf    [0] + mult_RK * x_pf    [M];
	x_sA 	[N+1] = x_sA    [0] + mult_RK * x_sA    [M];
	x_sB 	[N+1] = x_sB    [0] + mult_RK * x_sB    [M];
	x_fA 	[N+1] = x_fA    [0] + mult_RK * x_fA    [M];
}

/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/


/****************************************************************************************************/
/*									Function that adds all SRK terms								*/
/****************************************************************************************************/
void C_Column::add_RK(void) {
	extern const double h;
	Phi_pp	[0] += (Phi_pp	[4] + Phi_pp[5] * 2 + Phi_pp[6] * 2 + Phi_pp[7])/6;
	Phi_ps	[0] += (Phi_ps	[4] + Phi_ps[5] * 2 + Phi_ps[6] * 2 + Phi_ps[7])/6;
	Phi_pf	[0] += (Phi_pf	[4] + Phi_pf[5] * 2 + Phi_pf[6] * 2 + Phi_pf[7])/6;
	Phi_sA	[0] += (Phi_sA	[4] + Phi_sA[5] * 2 + Phi_sA[6] * 2 + Phi_sA[7])/6;
	Phi_sB	[0] += (Phi_sB	[4] + Phi_sB[5] * 2 + Phi_sB[6] * 2 + Phi_sB[7])/6;
	Phi_fA	[0] += (Phi_fA	[4] + Phi_fA[5] * 2 + Phi_fA[6] * 2 + Phi_fA[7])/6;
	x_pp	[0] += (x_pp	[4] + x_pp  [5] * 2 + x_pp  [6] * 2 + x_pp  [7])/6 + pow(gamma_p, 2) * h * Rand_vars[0];
	x_ps	[0] += (x_ps	[4] + x_ps  [5] * 2 + x_ps  [6] * 2 + x_ps  [7])/6 + pow(gamma_p, 2) * h * Rand_vars[2];
	x_pf	[0] += (x_pf	[4] + x_pf  [5] * 2 + x_pf  [6] * 2 + x_pf  [7])/6 + pow(gamma_p, 2) * h * Rand_vars[4];
	x_sA	[0] += (x_sA	[4] + x_sA  [5] * 2 + x_sA  [6] * 2 + x_sA  [7])/6;
	x_sB	[0] += (x_sB	[4] + x_sB  [5] * 2 + x_sB  [6] * 2 + x_sB  [7])/6;
	x_fA	[0] += (x_fA	[4] + x_fA  [5] * 2 + x_fA  [6] * 2 + x_fA  [7])/6;
	/* Generate noise for the next iteration */
	for (unsigned i=0; i<Rand_vars.size(); ++i) {
		Rand_vars[i] = MTRands[i]() + input;
	}
}
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/
