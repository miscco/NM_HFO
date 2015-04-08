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
/*									Functions for data storage										*/
/****************************************************************************************************/
#pragma once
#include "C_Column.h"
#include "H_Column.h"

/****************************************************************************************************/
/*											Save data												*/
/****************************************************************************************************/
inline void get_data(int counter, C_Column& C,  H_Column& H, double* V_C, double* V_H) {
	C.get_data(counter, V_C);
	H.get_data(counter, V_H);
}
/****************************************************************************************************/
/*										 		end													*/
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
