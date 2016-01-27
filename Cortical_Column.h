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
/*									Header file of a cortical module							*/
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
/*									Implementation of the CA1 module 								*/
/****************************************************************************************************/
class Cortical_Column {
public:
	/* Constructors */
	Cortical_Column(void)
	{set_RNG();}

	/* Initialize the RNGs */
	void 	set_RNG		(void);

	/* Set strength of input */
	void	set_input	(double I) {input = I;}

	/* Firing rates */
	double 	get_Qp		(int) const;
	double 	get_Qs		(int) const;
	double 	get_Qf		(int) const;

	/* Currents */
	double 	I_pp		(int) const;
	double 	I_ps		(int) const;
	double 	I_pf		(int) const;
	double 	I_sp		(int) const;
	double 	I_ss		(int) const;
	double 	I_sf		(int) const;
	double 	I_fp		(int) const;
	double 	I_ff		(int) const;
	double 	I_L_p		(int) const;
	double 	I_L_s		(int) const;
	double 	I_L_f		(int) const;

	/* Noise function */
	double 	noise_xRK 	(int, int) const;
	double 	noise_aRK 	(int) const;

	/* ODE functions */
	void 	set_RK		(int);
	void 	add_RK		(void);

	/* Data storage  access */
	void	get_data (int N, double* V) {V[N] = N_pp * y_pp[0] - N_fp * y_fA[0] - N_sp * (y_sA[0] + y_sB[0]);}
	/* Stimulation protocoll acces */
	friend class Stim;

private:
	/* Random number generators */
	vector<random_stream_normal> MTRands;

	/* Container for noise */
	vector<double>	Rand_vars;

	/* Declaration and Initialization of parameters */
	/* Membrane time in ms */
	const double 	tau_p 		= 3;
	const double 	tau_s 		= 3;
	const double 	tau_f 		= 3;

	/* Maximum firing rate in ms^-1 */
	const double 	Qp_max		= 5.E-3;
	const double 	Qs_max		= 5.E-3;
	const double 	Qf_max		= 5.E-3;

	/* Sigmoid threshold in mV */
	const double 	theta_p		= 1;
	const double 	theta_s		= 6;
	const double 	theta_f		= 6;

	/* Sigmoid gain in mV */
	const double 	sigma_p		= 0.56;
	const double 	sigma_s		= 0.56;
	const double 	sigma_f		= 0.56;

	/* PSP rise time in ms^-1 */
	const double 	gamma_p		= 180E-3;
	const double 	gamma_sA	= 33E-3;
	const double 	gamma_fA	= 220E-3;
	const double 	gamma_sB	= 3.3E-3;

	/* PSP amplitudes in mV */
	const double 	G_p         = 5;
	const double 	G_sA        = 50;
	const double 	G_fA        = 20;
	const double 	G_sB        = 3;

	/* Noise parameters in ms^-1 */
	const double	dphi		= 5E-3;
	double			input		= 0.0;

	/* Connectivities (dimensionless) */
	const double 	N_pp		= 200;
	const double 	N_ps		= 200;
	const double 	N_pf		= 200;
	const double 	N_sp		= 240;
	const double 	N_ss		= 400;
	const double 	N_sf		= 400;
	const double 	N_fp		= 100;
	const double 	N_ff		= 100;

	/* Parameters for SRK4 iteration */
	const vector<double> A = {0.5,  0.5,  1.0, 1.0};
	const vector<double> B = {0.75, 0.75, 0.0, 0.0};

	/* Population variables															*/
	/* Excitatory PSPs have to be treated individually as they are noisy.           */
	/* In contrast inhibitory PSP are noise free and therewith only computed once   */
	vector<double> 	y_pp	= _INIT(0.0),		/* Pyramidal to p  AMPA  PSP        */
					y_ps	= _INIT(0.0),		/* Pyramidal to s  AMPA  PSP        */
					y_pf	= _INIT(0.0),		/* Pyramidal to f  AMPA  PSP        */
					y_sA	= _INIT(0.0),		/* Slow inhibitory GABAA PSP        */
					y_sB	= _INIT(0.0),		/* Slow inhibitory GABAB PSP        */
					y_fA	= _INIT(0.0),		/* Fast inhibitory GABAA PSP        */
					x_pp    = _INIT(0.0),		/* derivative of y_pp				*/
					x_ps    = _INIT(0.0),		/* derivative of y_ps				*/
					x_pf    = _INIT(0.0),		/* derivative of y_pf				*/
					x_sA	= _INIT(0.0),		/* derivative of y_sA				*/
					x_sB	= _INIT(0.0),		/* derivative of y_sB				*/
					x_fA	= _INIT(0.0);		/* derivative of y_fB				*/
};
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/

