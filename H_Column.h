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
/*									Header file of a CA3 module									*/
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
/*									Implementation of the CA3 module 								*/
/****************************************************************************************************/
class H_Column {
public:
	/* Constructors */
	H_Column(void)
	{set_RNG();}

	/*
	H_Column(double* Par)
	 //:sigma_p 	(Par[0]),		dphi		(Par[1])
	{set_RNG();}
	*/

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

	/* ODE functions */
	void 	get_RK              (int);
	void 	get_intermediate_RK	(int, double);
	void 	add_RK              (void);

	/* Data storage  access */
	void	get_data (int N, double* V) {V[N] = N_pp * V_p[0] - N_fp * V_f[0];}

	/* Stimulation protocoll acces */
	friend class Stim;

private:
	/* Population variables */
	vector<double> 	V_p		= _INIT(E_L),		/* pyramidal membrane voltage       */
					V_f     = _INIT(E_L),		/* fast inhibitory membrane voltage */
					Phi_pp	= _INIT(0.0),		/* PostSP p to p                    */
					Phi_pf	= _INIT(0.0),		/* PostSP p to f                    */
					Phi_fA	= _INIT(0.0),		/* PostSP f to p                    */
					x_pp	= _INIT(0.0),		/* derivative of Phi_pp				*/
					x_pf	= _INIT(0.0),		/* derivative of Phi_pf				*/
					x_fA	= _INIT(0.0);		/* derivative of Phi_ff				*/

	/* Random number generators */
	vector<GEN>		MTRands;

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

	/* Conductivities in mS/cm^-2 */
	/* Leak */
	const double 	g_L    		= 1.;

	/* Reversal potentials in mV */
	/* synaptic */
	const double 	E_AMPA  	= 0;
	const double 	E_GABA  	= -70;

	/* Leak */
	const double 	E_L 		= -60;

	/* Noise parameters in ms^-1 */
	const double 	mphi		= 0E-3;
	const double	dphi		= 1E-3;
	double			input		= 0.0;

	/* Connectivities (dimensionless) */
	const double 	N_pp		= 40;
	const double 	N_pf		= 60;
	const double 	N_fp		= 40;
	const double 	N_ff		= 55;
};
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/

