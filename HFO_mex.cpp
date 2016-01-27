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
 *				partial seizures: From 'altered structure' to 'dysfunction'.
 *				B Molaee-Ardekani, P Benquet, F Bartolomei, F Wendling.
 *				NeuroImage 52(3):1109-1122 (2010)
 */

/****************************************************************************************************/
/* 		Implementation of the simulation as MATLAB routine (mex compiler)							*/
/* 		mex command is given by:																	*/
/* 		mex CXXFLAGS="\$CXXFLAGS -std=c++11 -O3" HFO_mex.cpp CA3_Column.cpp Cortical_Column.cpp     */
/****************************************************************************************************/
#include "mex.h"
#include "matrix.h"
#include "Data_Storage.h"
#include "ODE.h"
mxArray* SetMexArray(int N, int M);

/****************************************************************************************************/
/*										Fixed simulation settings									*/
/****************************************************************************************************/
extern const int onset	= 10;								/* time until data is stored in  s		*/
extern const int res 	= 1E4;								/* number of iteration steps per s		*/
extern const double dt 	= 1E3/res;							/* duration of a timestep in ms			*/
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
	mxArray* V_C		= SetMexArray(1, T*res);
	mxArray* V_H		= SetMexArray(1, T*res);
	mxArray* Y_H		= SetMexArray(1, T*res);

	/* Pointer to the actual data block */
	double* Pr_V_C	= mxGetPr(V_C);
	double* Pr_V_H	= mxGetPr(V_H);
	double* Pr_Y_H	= mxGetPr(Y_H);

	/* Simulation */
	int count = 0;
	for (int t=0; t<Time; ++t) {
		ODE (Cortex, HFO);
		if(t>=onset*res){
			get_data(count, Cortex, HFO, Pr_V_C, Pr_V_H, Pr_Y_H);
			++count;
		}
	}

	/* Output of the simulation */
	plhs[0] = V_C;
	plhs[1] = V_H;
	plhs[2] = Y_H;
return;
}
/****************************************************************************************************/
/*												end													*/
/****************************************************************************************************/


/****************************************************************************************************/
/*									Create MATLAB data container									*/
/****************************************************************************************************/
mxArray* SetMexArray(int N, int M) {
	mxArray* Array	= mxCreateDoubleMatrix(0, 0, mxREAL);
	mxSetM(Array, N);
	mxSetN(Array, M);
	mxSetData(Array, mxMalloc(sizeof(double)*M*N));
	return Array;
}
/****************************************************************************************************/
/*										 		end													*/
/****************************************************************************************************/
