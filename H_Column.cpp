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
*	THE SOFTWARE IS PROV_fDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
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
#include "H_Column.h"

/****************************************************************************************************/
/*										 Initialization of RNG 										*/
/****************************************************************************************************/
void H_Column::set_RNG(void) {
	/* Number of independent streams */
	int N = 4;

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
double H_Column::get_Qp	(int N) const{
	double q = Qp_max / (1 + exp(-C1 * (V_p[N] - theta_p) / sigma_p));
	return q;
}

/* Inhibitory firing rate */
double H_Column::get_Qf	(int N) const{
	double q = Qf_max / (1 + exp(-C1 * (V_f[N] - theta_f) / sigma_f));
	return q;
}
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/


/****************************************************************************************************/
/*										Synaptic currents											*/
/****************************************************************************************************/
/* Excitatory input to pyramidal population */
double H_Column::I_pp	(int N) const{
	double I = Phi_pp[N] * (V_p[N] - E_AMPA);
	return I;
}

/* Excitatory input to inhibitory population */
double H_Column::I_pf	(int N) const{
	double I = Phi_pf[N] * (V_f[N] - E_AMPA);
	return I;
}

/* Inhibitory input to pyramidal population */
double H_Column::I_fp	(int N) const{
	double I = Phi_fA[N] * N_fp * (V_p[N] - E_GABA);
	return I;
}

/* Inhibitory input to inhibitory population */
double H_Column::I_ff	(int N) const{
	double I = Phi_fA[N] * N_ff * (V_f[N] - E_GABA);
	return I;
}
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/


/****************************************************************************************************/
/*										 Current functions 											*/
/****************************************************************************************************/
/* Leak current of pyramidal population */
double H_Column::I_L_p	(int N) const{
	double I = g_L * (V_p[N] - E_L);
	return I;
}

/* Leak current of inhibitory population */
double H_Column::I_L_f	(int N) const{
	double I = g_L * (V_f[N] - E_L);
	return I;
}
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/


/****************************************************************************************************/
/*										 RK noise scaling 											*/
/****************************************************************************************************/
double H_Column::noise_xRK(int N, int M) const{
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
void H_Column::get_RK (int N) {
	extern const double dt;
	int M = N+4;
	V_p	  	[M] = dt*(-(I_L_p(N) + I_pp(N) + I_fp(N) )/tau_p);
	V_f	  	[M] = dt*(-(I_L_f(N) + I_pf(N) + I_ff(N) )/tau_f);
	Phi_pp	[M] = dt*(x_pp[N]);
	Phi_pf	[M] = dt*(x_pf[N]);
	Phi_fA	[M] = dt*(x_fA[N]);
	x_pp  	[M] = dt*(pow(gamma_p,  2) * (N_pp * get_Qp(N) + noise_xRK(N, 0)- Phi_pp[N]) - 2 * gamma_p  * x_pp[N]);
	x_pf  	[M] = dt*(pow(gamma_p,  2) * (N_pf * get_Qp(N) + noise_xRK(N, 1)- Phi_pf[N]) - 2 * gamma_p  * x_pf[N]);
	x_fA  	[M] = dt*(pow(gamma_fA, 2) * (		 get_Qf(N)					- Phi_fA[N]) - 2 * gamma_fA * x_fA[N]);

	/* Calculate the intermediate RK terms */
	switch (N){
	case 0:
	case 1:
		get_intermediate_RK(N, 0.5);
		break;
	case 2:
		get_intermediate_RK(N, 1.);
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
void H_Column::get_intermediate_RK (int N, double mult_RK) {
	int M = N+4;
	V_p	  	[N+1] = V_p     [0] + mult_RK * V_p     [M];
	V_f	  	[N+1] = V_f     [0] + mult_RK * V_f     [M];
	Phi_pp 	[N+1] = Phi_pp  [0] + mult_RK * Phi_pp  [M];
	Phi_pf 	[N+1] = Phi_pf  [0] + mult_RK * Phi_pf  [M];
	Phi_fA 	[N+1] = Phi_fA  [0] + mult_RK * Phi_fA  [M];
	x_pp 	[N+1] = x_pp    [0] + mult_RK * x_pp    [M];
	x_pf 	[N+1] = x_pf    [0] + mult_RK * x_pf    [M];
	x_fA 	[N+1] = x_fA    [0] + mult_RK * x_fA    [M];
}

/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/


/****************************************************************************************************/
/*									Function that adds all SRK terms								*/
/****************************************************************************************************/
void H_Column::add_RK(void) {
	extern const double h;
	V_p	  	[0] += (V_p		[4] + V_p	[5] * 2 + V_p	[6] * 2 + V_p	[7])/6;
	V_f	  	[0] += (V_f		[4] + V_f	[5] * 2 + V_f	[6] * 2 + V_f	[7])/6;
	Phi_pp	[0] += (Phi_pp	[4] + Phi_pp[5] * 2 + Phi_pp[6] * 2 + Phi_pp[7])/6;
	Phi_pf	[0] += (Phi_pf	[4] + Phi_pf[5] * 2 + Phi_pf[6] * 2 + Phi_pf[7])/6;
	Phi_fA	[0] += (Phi_fA	[4] + Phi_fA[5] * 2 + Phi_fA[6] * 2 + Phi_fA[7])/6;
	x_pp  	[0] += (x_pp	[4] + x_pp	[5] * 2 + x_pp	[6] * 2 + x_pp	[7])/6 + pow(gamma_p, 2) * h * Rand_vars[0];
	x_pf  	[0] += (x_pf	[4] + x_pf	[5] * 2 + x_pf	[6] * 2 + x_pf	[7])/6 + pow(gamma_p, 2) * h * Rand_vars[2];
	x_fA  	[0] += (x_fA	[4] + x_fA	[5] * 2 + x_fA	[6] * 2 + x_fA	[7])/6;

	/* Generate noise for the next iteration */
	for (unsigned i=0; i<Rand_vars.size(); ++i) {
		Rand_vars[i] = MTRands[i]() + input;
	}
}
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/
