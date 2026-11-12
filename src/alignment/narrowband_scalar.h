/*
Tandem Repeats Finder 
Copyright (C) 1999-2020 Gary Benson

This file is part of the Tandem Repeats Finder (TRF) program.

TRF is free software: you can redistribute it and/or modify 
it under the terms of the GNU Affero General Public License as 
published by the Free Software Foundation, either version 3 of 
the License, or (at your option) any later version.

TRF is distributed in the hope that it will be useful, but 
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See 
the GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public 
License along with TRF.  If not, see <https://www.gnu.org/licenses/>.
*/

/* Scalar narrowband kernel, and the row bookkeeping it shares with the vector kernel. */

#ifndef TRF_ALIGNMENT_NARROWBAND_SCALAR_H
#define TRF_ALIGNMENT_NARROWBAND_SCALAR_H

/* Fast C remainder for the common [0, 2*size) row-coordinate range.
 * Negative sentinel expressions and larger values retain C's signed remainder.
 * size is positive; value-size cannot overflow after the nonnegative check. */
static inline int trf_narrowband_remainder(int value, int size)
{
    if(value>=0 && value-size<size)
        return value<size?value:value-size;
    return value%size;
}

/*******************************************************************/

/***************************** narrowbandwrap()  macros **************************/

/*******************************************************************/

/* *pcurr>=maxscore added 2.17.05 gary benson -- to extend alignment as far as possible */
#define test_trace_and_backwards_maxscore \
	if((realr<=start-size)&&(*pcurr==0))\
pleft=(*pcurr=-1000);\
else end_of_trace=FALSE;\
if(*pcurr>=maxscore)\
{\
	maxscore=*pcurr;\
	minrealrow=realr;\
	mincol=c;\
	/***6/9/05 G. Benson***/ mincolbandcenter=Bandcenter[r];\
	/*** 6/14/05 G. Benson ***/ mincolposition=i;\
}


#define test_maxrowscore_with_match \
	if(*pcurr>maxrowscore)\
{\
	maxrowscore=*pcurr;\
	if((*pcurr==diag_score)&&(match_yes_no==Alpha))\
	matchatmax_col=c;\
	else matchatmax_col=-2;\
}

#define test_maxrowscore_without_match \
	if(*pcurr>maxrowscore)\
matchatmax_col=-2;

/* *pcurr>=maxscore added 2.17.05 gary benson -- to extend alignment as far as possible */     
#define test_trace_and_forward_maxscore \
	if((realr>=start)&&(*pcurr==0))\
pleft=(*pcurr=-1000);\
else end_of_trace=FALSE;\
if(*pcurr>=maxscore)\
{\
	maxscore=*pcurr;\
	maxrealrow=realr;\
	maxrow=r;\
	maxcol=c;\
}

void trf_narrowband_scalar(int start, int size, int bandradius, int bandradiusforward,
                           int option, int tuplesize);

#endif
