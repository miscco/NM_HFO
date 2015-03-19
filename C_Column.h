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
*/

/************************************************************************************************/
/*									Header file of a cortical module							*/
/************************************************************************************************/
#pragma once
#include <cmath>
#include <vector>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>
using std::vector;

/****************************************************************************************************/
/*										Typedefs for RNG											*/
/****************************************************************************************************/
typedef boost::mt11213b                    	ENG;    /* Mersenne Twister		*/
typedef boost::normal_distribution<double>	DIST;   /* Normal Distribution	*/
typedef boost::variate_generator<ENG,DIST> 	GEN;    /* Variate generator	*/
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/


/****************************************************************************************************/
/*									Macro for vector initialization									*/
/****************************************************************************************************/
#ifndef _INIT
#define _INIT(x)	{x, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
#endif
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/


/****************************************************************************************************/
/*									Implementation of the CA1 module 								*/
/****************************************************************************************************/
class C_Column {
public:
	/* Constructors */
	C_Column(void)
	{set_RNG();}

	/*
	C_Column(double* Par)
	 //:sigma_p 	(Par[0]),		dphi		(Par[1])
	{set_RNG();}
	*/

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

	/* ODE functions */
	void 	set_RK		(int);
	void 	add_RK		(void);
	inline void set_intermediate_RK	(int, double);

	/* Data storage  access */
	void	get_data (int N, double* V) {V[N] = N_pp * Phi_pp[0] - N_fp * Phi_fA[0] - N_sp * (Phi_sA[0] + Phi_sB[0]);}
	/* Stimulation protocoll acces */
	friend class Stim;

private:
	/* Population variables															*/
	/* Excitatory PSPs have to be treated individually as they are noisy.           */
	/* In contrast inhibitory PSP are noise free and therewith only computed once   */

	vector<double> 	Phi_pp	= _INIT(0.0),		/* Pyramidal to p  AMPA  PSP        */
					Phi_ps	= _INIT(0.0),		/* Pyramidal to s  AMPA  PSP        */
					Phi_pf	= _INIT(0.0),		/* Pyramidal to f  AMPA  PSP        */
					Phi_sA	= _INIT(0.0),		/* Slow inhibitory GABAA PSP        */
					Phi_sB	= _INIT(0.0),		/* Slow inhibitory GABAB PSP        */
					Phi_fA	= _INIT(0.0),		/* Fast inhibitory GABAA PSP        */
					x_pp    = _INIT(0.0),		/* derivative of Phi_pp				*/
					x_ps    = _INIT(0.0),		/* derivative of Phi_ps				*/
					x_pf    = _INIT(0.0),		/* derivative of Phi_pf				*/
					x_sA	= _INIT(0.0),		/* derivative of Phi_sA				*/
					x_sB	= _INIT(0.0),		/* derivative of Phi_sB				*/
					x_fA	= _INIT(0.0);		/* derivative of Phi_fB				*/

	/* Random number generators */
	vector<GEN>		MTRands;

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

	/* Noise parameters in ms^-1 */
	const double 	mphi		= 0E-3;
	const double	dphi		= 10E-3;
	double			input		= 0.0;

	/* Connectivities (dimensionless) */
	const double 	N_pp		= 20;
	const double 	N_ps		= 20;
	const double 	N_pf		= 20;
	const double 	N_sp		= 10;
	const double 	N_ss		= 10;
	const double 	N_sf		= 10;
	const double 	N_fp		= 10;
	const double 	N_ff		= 10;
};
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/

