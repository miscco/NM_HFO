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
/*									Functions of the cortical module								*/
/****************************************************************************************************/
#include "Cortical_Column.h"

/****************************************************************************************************/
/*										 Initialization of RNG 										*/
/****************************************************************************************************/
void Cortical_Column::set_RNG(void) {
	extern const double dt;
	/* Number of independent random variables */
	int N = 3;

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
double Cortical_Column::noise_xRK(int N, int M) const{
	return gamma_p * gamma_p * (Rand_vars[2*M] + Rand_vars[2*M+1]/std::sqrt(3))*B[N];
}

double Cortical_Column::noise_aRK(int M) const{
	return gamma_p * gamma_p * (Rand_vars[2*M] - Rand_vars[2*M+1]*std::sqrt(3))/4;
}
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/


/****************************************************************************************************/
/*										 Firing Rate functions 										*/
/****************************************************************************************************/
/* Pyramidal firing rate */
double Cortical_Column::get_Qp	(int N) const{
	double q = Qp_max / (1 + exp(-(N_pp * y_pp[N] - N_sp * y_sA[N] - N_fp * y_fA[N] - theta_p) / sigma_p));
	return q;
}

/* Slow inhibitory firing rate */
double Cortical_Column::get_Qs	(int N) const{
	double q = Qs_max / (1 + exp(-(N_ps * y_ps[N] - N_ss * y_sA[N] - theta_s) / sigma_s));
	return q;
}

/* Fast inhibitory firing rate */
double Cortical_Column::get_Qf	(int N) const{
	double q = Qf_max / (1 + exp(-(N_pf * y_pf[N] - N_sf * y_sA[N] - N_ff * y_fA[N] - theta_f) / sigma_f));
	return q;
}
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/


/****************************************************************************************************/
/*										Calculate the Nth SRK term									*/
/****************************************************************************************************/
void Cortical_Column::set_RK (int N) {
	extern const double dt;
	y_pp	[N+1] = y_pp[0] + A[N]*dt*(x_pp[N]);
	y_ps	[N+1] = y_ps[0] + A[N]*dt*(x_ps[N]);
	y_pf	[N+1] = y_pf[0] + A[N]*dt*(x_pf[N]);
	y_sA	[N+1] = y_sA[0] + A[N]*dt*(x_sA[N]);
	y_sB	[N+1] = y_sB[0] + A[N]*dt*(x_sB[N]);
	y_fA	[N+1] = y_fA[0] + A[N]*dt*(x_fA[N]);
	x_pp  	[N+1] = x_pp[0] + A[N]*dt*gamma_p*(G_p * (get_Qp(N) - y_pp[N]) - 2 * x_pp[N]) + noise_xRK(N, 0);
	x_ps  	[N+1] = x_ps[0] + A[N]*dt*gamma_p*(G_p * (get_Qp(N) - y_ps[N]) - 2 * x_ps[N]) + noise_xRK(N, 1);
	x_pf  	[N+1] = x_pf[0] + A[N]*dt*gamma_p*(G_p * (get_Qp(N) - y_pf[N]) - 2 * x_pf[N]) + noise_xRK(N, 2);
	x_sA  	[N+1] = x_sA[0] + A[N]*dt*gamma_p*(G_sA* (get_Qs(N) - y_sA[N]) - 2 * x_sA[N]);
	x_sB  	[N+1] = x_sB[0] + A[N]*dt*gamma_p*(G_sB* (get_Qs(N) - y_sB[N]) - 2 * x_sB[N]);
	x_fA  	[N+1] = x_fA[0] + A[N]*dt*gamma_p*(G_fA* (get_Qs(N) - y_fA[N]) - 2 * x_fA[N]);
}
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/


/****************************************************************************************************/
/*									Function that adds all SRK terms								*/
/****************************************************************************************************/
void Cortical_Column::add_RK(void) {
	y_pp[0] = (-3*y_pp[0] + 2*y_pp[1] + 4*y_pp[2] + 2*y_pp[3] + y_pp[4])/6;
	y_ps[0] = (-3*y_ps[0] + 2*y_ps[1] + 4*y_ps[2] + 2*y_ps[3] + y_ps[4])/6;
	y_pf[0] = (-3*y_pf[0] + 2*y_pf[1] + 4*y_pf[2] + 2*y_pf[3] + y_pf[4])/6;
	y_sA[0] = (-3*y_sA[0] + 2*y_sA[1] + 4*y_sA[2] + 2*y_sA[3] + y_sA[4])/6;
	y_sB[0] = (-3*y_sB[0] + 2*y_sB[1] + 4*y_sB[2] + 2*y_sB[3] + y_sB[4])/6;
	y_fA[0] = (-3*y_fA[0] + 2*y_fA[1] + 4*y_fA[2] + 2*y_fA[3] + y_fA[4])/6;
	x_pp[0] = (-3*x_pp[0] + 2*x_pp[1] + 4*x_pp[2] + 2*x_pp[3] + x_pp[4])/6 + noise_aRK(0);
	x_ps[0] = (-3*x_ps[0] + 2*x_ps[1] + 4*x_ps[2] + 2*x_ps[3] + x_ps[4])/6 + noise_aRK(1);
	x_pf[0] = (-3*x_pf[0] + 2*x_pf[1] + 4*x_pf[2] + 2*x_pf[3] + x_pf[4])/6 + noise_aRK(2);
	x_sA[0] = (-3*x_sA[0] + 2*x_sA[1] + 4*x_sA[2] + 2*x_sA[3] + x_sA[4])/6;
	x_sB[0] = (-3*x_sB[0] + 2*x_sB[1] + 4*x_sB[2] + 2*x_sB[3] + x_sB[4])/6;
	x_fA[0] = (-3*x_fA[0] + 2*x_fA[1] + 4*x_fA[2] + 2*x_fA[3] + x_fA[4])/6;
	/* Generate noise for the next iteration */
	for (unsigned i=0; i<Rand_vars.size(); ++i) {
		Rand_vars[i] = MTRands[i]() + input;
	}
}
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/
