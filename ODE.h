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
/*									Implementation of the ODE solver								*/
/****************************************************************************************************/
#pragma once
#include "CA3_Column.h"
#include "Cortical_Column.h"

/****************************************************************************************************/
/*										Evaluation of SRK4											*/
/****************************************************************************************************/
void ODE(Cortical_Column& Cortex, CA3_Column& CA3) {
	/* First calculating every ith RK moment. Has to be in order, 1th moment first */
	for (int i=0; i<4; ++i) {
		Cortex.set_RK(i);
		CA3.get_RK(i);
	}

	/* Add all moments */
	Cortex.add_RK();
	CA3.add_RK();
}
/****************************************************************************************************/
/*										 		end													*/
/****************************************************************************************************/
