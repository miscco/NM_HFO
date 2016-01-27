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

/************************************************************************************************/
/*									Header file of a CA3 module									*/
/************************************************************************************************/
#pragma once
#include <cmath>
#include <vector>
#include "Random_Stream.h"
using std::vector;

/****************************************************************************************************/
/*									Macro for vector initialization									*/
/****************************************************************************************************/
#ifndef _INIT
#define _INIT(x)	{x, 0.0, 0.0, 0.0, 0.0}
#endif
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/


/****************************************************************************************************/
/*									Implementation of the CA3 module 								*/
/****************************************************************************************************/
class CA3_Column {
public:
	/* Constructors */
	CA3_Column(void)
	{set_RNG();}

	/* Initialize the RNGs */
	void 	set_RNG		(void);

	/* Set strength of input */
	void	set_input	(double I) {input = I;}

	/* Firing rates */
	double 	get_Qp		(int) const;
	double 	get_Qf		(int) const;

	/* Currents */
	double 	I_pp		(int) const;
	double 	I_pf		(int) const;
	double 	I_fp		(int) const;
	double 	I_ff		(int) const;
	double 	I_L_p		(int) const;
	double 	I_L_f		(int) const;

	/* Noise function */
	double 	noise_xRK 	(int, int) const;
	double 	noise_aRK 	(int) const;

	/* ODE functions */
	void 	get_RK		(int);
	void 	add_RK		(void);

	/* Data storage  access */
	void	get_data (int N, double* V, double * Y) {V[N] = V_p[0]; Y[N] = N_pp*y_pp[0] - N_fp*y_fA[0];}

private:
	/* Random number generators */
	vector<random_stream_normal> MTRands;

	/* Container for noise */
	vector<double>	Rand_vars;

	/* Declaration and Initialization of parameters */
	/* Membrane time in ms */
	const double 	tau_p 		= 1.;
	const double 	tau_f 		= 1.;

	/* Maximum firing rate in ms^-1 */
	const double 	Qp_max		= 30.E-3;
	const double 	Qf_max		= 60.E-3;

	/* Sigmoid threshold in mV */
	const double 	theta_p		= -58.5;
	const double 	theta_f		= -58.5;

	/* Sigmoid gain in mV */
	const double 	sigma_p		= 4;
	const double 	sigma_f		= 6;

	/* Scaling parameter for sigmoidal mapping (dimensionless) */
	const double 	C1          = (3.14159265/sqrt(3));

	/* PSP rise time in ms^-1 */
	const double 	gamma_p		= 180E-3;
	const double 	gamma_fA	= 220E-3;

	/* PSP amplitude in mV */
	const double 	G_p         = 18;
	const double 	G_fA        = 30;

	/* Conductivities */
	/* Leak */
	const double 	g_L    		= 1.;

	/* Reversal potentials in mV */
	/* synaptic */
	const double 	E_AMPA  	= 0;
	const double 	E_GABA  	= -70;

	/* Leak */
	const double 	E_L 		= -60;

	/* Noise parameters in ms^-1 */
	const double	dphi		= 5E-3;
	double			input		= 0.0;

	/* Connectivities (dimensionless) */
	const double 	N_pp		= 280;
	const double 	N_pf		= 600;
	const double 	N_fp		= 280;
	const double 	N_ff		= 400;

	/* Parameters for SRK4 iteration */
	const vector<double> A = {0.5,  0.5,  1.0, 1.0};
	const vector<double> B = {0.75, 0.75, 0.0, 0.0};

	/* Population variables */
	vector<double> 	V_p		= _INIT(E_L),		/* pyramidal membrane voltage       */
					V_f     = _INIT(E_L),		/* fast inhibitory membrane voltage */
					y_pp	= _INIT(0.0),		/* PostSP p to p                    */
					y_pf	= _INIT(0.0),		/* PostSP p to f                    */
					y_fA	= _INIT(0.0),		/* PostSP f to p                    */
					x_pp	= _INIT(0.0),		/* derivative of y_pp				*/
					x_pf	= _INIT(0.0),		/* derivative of y_pf				*/
					x_fA	= _INIT(0.0);		/* derivative of y_ff				*/
};
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/

