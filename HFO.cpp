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

/****************************************************************************************************/
/* 		Implementation of the simulation as MATLAB routine (mex compiler)							*/
/* 		mex command is given by:																	*/
/* 		mex CXXFLAGS="\$CXXFLAGS -std=c++11 -O3" HFO.cpp C_Column.cpp H_Column.cpp		            */
/*		The Simulation requires the following boost libraries:	Random								*/
/****************************************************************************************************/
#include "mex.h"
#include "matrix.h"
#include "C_Column.h"
#include "H_Column.h"
#include "saves.h"
#include "ODE.h"

/****************************************************************************************************/
/*										Fixed simulation settings									*/
/****************************************************************************************************/
extern const int onset	= 10;								/* time until data is stored in  s		*/
extern const int res 	= 1E4;								/* number of iteration steps per s		*/
extern const int red 	= 1;								/* number of iterations that is saved	*/
extern const double dt 	= 1E3/res;							/* duration of a timestep in ms			*/
extern const double h	= sqrt(dt);							/* squareroot of dt for SRK iteration	*/
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/


/****************************************************************************************************/
/*										Simulation routine	 										*/
/*										lhs defines outputs											*/
/*										rhs defines inputs											*/
/****************************************************************************************************/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	/* Initialize the seeder */
	srand(time(NULL));

	/* Fetch inputs */
	const int T				= (int) (mxGetScalar(prhs[0]));	/* Duration of simulation in s			*/
	const int Time 			= (T+onset)*res;				/* Total number of iteration steps		*/
	double* Param_C         = mxGetPr (prhs[1]);			/* Parameters of C module				*/
	double* Param_H         = mxGetPr (prhs[2]);			/* Parameters of H module				*/

	/* Initialize the populations */
	C_Column Cortex;
	H_Column HFO;

	/* Create data containers */
	mxArray* V_C		= SetMexArray(1, T*res/red);
	mxArray* V_H		= SetMexArray(1, T*res/red);

	/* Pointer to the actual data block */
	double* Pr_V_C	= mxGetPr(V_C);
	double* Pr_V_H	= mxGetPr(V_H);

	/* Simulation */
	int count = 0;
	for (int t=0; t<Time; ++t) {
		ODE (Cortex, HFO);
		if(t>=onset*res && t%red==0){
			get_data(count, Cortex, HFO, Pr_V_C, Pr_V_H);
			++count;
		}
	}

	/* Output of the simulation */
	plhs[0] = V_C;
	plhs[1] = V_H;
return;
}
/****************************************************************************************************/
/*												end													*/
/****************************************************************************************************/
