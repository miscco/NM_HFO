/*
 *	Copyright (c) 2015 University of Lübeck
 *
 *	Permission is hereby granted, free of charge, to any person obtaining a copy
 *	of this software and associated documentation files (the "Software"), to deal
 *	in the Software without restriction, including without limitation the rights
 *	to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 *	copies of the Software, and to permit persons to whom the Software is
 *	furnished to do so, subject to the following conditions:
 *
 *	The above copyright notice and this permission notice shall be included in
 *	all copies or substantial portions of the Software.
 *
 *	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 *	IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 *	FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 *	AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 *	LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 *	OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 *	THE SOFTWARE.
 *
 *	AUTHORS:	Michael Schellenberger Costa: mschellenbergercosta@gmail.com
 *
 *	Based on:	Computational modeling of high-frequency oscillations at the onset of neocortical
 *				partial seizures: From 'altered structure' to 'dysfunction'
 *				B Molaee-Ardekani, P Benquet, F Bartolomei, F Wendling.
 *				NeuroImage 52(3):1109-1122 (2010)
 */

/****************************************************************************************************/
/*									Functions of the CA3 module										*/
/****************************************************************************************************/
#include "CA3_Column.h"

/****************************************************************************************************/
/*										 Initialization of RNG 										*/
/****************************************************************************************************/
void CA3_Column::set_RNG(void) {
	extern const double dt;
	/* Number of independent random variables */
	int N = 2;

	/* Create RNG for each stream */
	for (int i=0; i<N; ++i){
		/* Add the RNG for I_{l}*/
		MTRands.push_back(random_stream_normal(0.0, dphi*dt));

		/* Add the RNG for I_{l,0} */
		MTRands.push_back(random_stream_normal(0.0, dt));

		/* Get the random number for the first iteration */
		Rand_vars.push_back(MTRands[2*i]());
		Rand_vars.push_back(MTRands[2*i+1]());
	}
}
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/


/****************************************************************************************************/
/*										 RK noise scaling 											*/
/****************************************************************************************************/
double CA3_Column::noise_xRK(int N, int M) const{
	return gamma_p * gamma_p * (Rand_vars[2*M] + Rand_vars[2*M+1]/std::sqrt(3))*B[N];
}

double CA3_Column::noise_aRK(int M) const{
	return gamma_p * gamma_p * (Rand_vars[2*M] - Rand_vars[2*M+1]*std::sqrt(3))/4;
}
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/


/****************************************************************************************************/
/*										 Firing Rate functions 										*/
/****************************************************************************************************/
/* Pyramidal firing rate */
double CA3_Column::get_Qp	(int N) const{
	double q = Qp_max / (1 + exp(-C1 * (N_pp * y_pp[N] - N_fp * y_fA[N] - theta_p) / sigma_p));
	return q;
}

/* Inhibitory firing rate */
double CA3_Column::get_Qf	(int N) const{
	double q = Qf_max / (1 + exp(-C1 * (N_pf * y_pf[N] - N_ff * y_fA[N]- theta_f) / sigma_f));
	return q;
}
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/


/****************************************************************************************************/
/*										Synaptic currents											*/
/****************************************************************************************************/
/* Excitatory input to pyramidal population */
double CA3_Column::I_pp	(int N) const{
	double I = y_pp[N] * (V_p[N] - E_AMPA);
	return I;
}

/* Excitatory input to inhibitory population */
double CA3_Column::I_pf	(int N) const{
	double I = y_pf[N] * (V_f[N] - E_AMPA);
	return I;
}

/* Inhibitory input to pyramidal population */
double CA3_Column::I_fp	(int N) const{
	double I = y_fA[N] * N_fp * (V_p[N] - E_GABA);
	return I;
}

/* Inhibitory input to inhibitory population */
double CA3_Column::I_ff	(int N) const{
	double I = y_fA[N] * N_ff * (V_f[N] - E_GABA);
	return I;
}
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/


/****************************************************************************************************/
/*										 Current functions 											*/
/****************************************************************************************************/
/* Leak current of pyramidal population */
double CA3_Column::I_L_p	(int N) const{
	double I = g_L * (V_p[N] - E_L);
	return I;
}

/* Leak current of inhibitory population */
double CA3_Column::I_L_f	(int N) const{
	double I = g_L * (V_f[N] - E_L);
	return I;
}
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/


/****************************************************************************************************/
/*										Calculate the Nth SRK term									*/
/****************************************************************************************************/
void CA3_Column::get_RK (int N) {
	extern const double dt;
	V_p	[N+1] = V_p [0] + A[N]*dt*(-(I_L_p(N) + I_pp(N) + I_fp(N) )/tau_p);
	V_f	[N+1] = V_f [0] + A[N]*dt*(-(I_L_f(N) + I_pf(N) + I_ff(N) )/tau_f);
	y_pp[N+1] = y_pp[0] + A[N]*dt*(x_pp[N]);
	y_pf[N+1] = y_pf[0] + A[N]*dt*(x_pf[N]);
	y_fA[N+1] = y_fA[0] + A[N]*dt*(x_fA[N]);
	x_pp[N+1] = x_pp[0] + A[N]*dt*gamma_p *(G_p * (get_Qp(N) - y_pp[N]) - 2 * x_pp[N]) + noise_xRK(N, 0);
	x_pf[N+1] = x_pf[0] + A[N]*dt*gamma_p *(G_p * (get_Qp(N) - y_pf[N]) - 2 * x_pf[N]) + noise_xRK(N, 1);
	x_fA[N+1] = x_fA[0] + A[N]*dt*gamma_fA*(G_fA* (get_Qf(N) - y_fA[N]) - 2 * x_fA[N]);
}
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/


/****************************************************************************************************/
/*									Function that adds all SRK terms								*/
/****************************************************************************************************/
void CA3_Column::add_RK(void) {
	V_p	[0] = (-3*V_p [0] + 2*V_p [1] + 4*V_p [2] + 2*V_p [3] + V_p	[4])/6;
	V_f	[0] = (-3*V_f [0] + 2*V_f [1] + 4*V_f [2] + 2*V_f [3] + V_f	[4])/6;
	y_pp[0] = (-3*y_pp[0] + 2*y_pp[1] + 4*y_pp[2] + 2*y_pp[3] + y_pp[4])/6;
	y_pf[0] = (-3*y_pf[0] + 2*y_pf[1] + 4*y_pf[2] + 2*y_pf[3] + y_pf[4])/6;
	y_fA[0] = (-3*y_fA[0] + 2*y_fA[1] + 4*y_fA[2] + 2*y_fA[3] + y_fA[4])/6;
	x_pp[0] = (-3*x_pp[0] + 2*x_pp[1] + 4*x_pp[2] + 2*x_pp[3] + x_pp[4])/6 + noise_aRK(0);
	x_pf[0] = (-3*x_pf[0] + 2*x_pf[1] + 4*x_pf[2] + 2*x_pf[3] + x_pf[4])/6 + noise_aRK(1);
	x_fA[0] = (-3*x_fA[0] + 2*x_fA[1] + 4*x_fA[2] + 2*x_fA[3] + x_fA[4])/6;

	/* Generate noise for the next iteration */
	for (unsigned i=0; i<Rand_vars.size(); ++i) {
		Rand_vars[i] = MTRands[i]() + input;
	}
}
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/
