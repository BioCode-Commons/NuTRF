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

/* Compile-time limits, option codes, and ordered maximum macros. */

#ifndef TRF_CORE_CONSTANTS_H
#define TRF_CORE_CONSTANTS_H

#define WEIGHTCONSENSUS 0
#define REPEATCONSENSUS 0
#define MAXDISTANCECONSTANT 2000 /* largest period size */
#define MINDISTANCE 10
#define MINBANDRADIUS 6
#define RECENTERCRITERION 3
#define MAXBANDWIDTH 150
#define MAXTUPLESIZES 10

#define MAXPATTERNSIZECONSTANT MAXDISTANCECONSTANT
#define DASH '-'
#define BLANK ' '

#define WITHCONSENSUS 1
#define WITHOUTCONSENSUS 0

#define MULTIPLES 3 /* We keep MULTIPLES*d distance entries */
#define FILLMULTIPLE 3
#define LOWERENDMULTIPLE                     \
    3 /* how close the lower end of distance \
         indices has to be to the end of one copy */
#define ACCEPTED 1
#define NOTACCEPTED 0

#define SMALLDISTANCE 20 /* these are sizes for which we do a */
/* full array wraparound dynamic */
/* programming and for which we set d_range */
/* by hand */

#define TAGSEP 50 /* index separation for tags table for linking */
/* active distances for distance range addition of */
/* matches */

#define Lookratio .4

#define NUMBER_OF_PERIODS 5 /* determines 5 best periods for a repeat */
#define NUMBER_OF_PERIODS_TO_TEST \
    3 /* only test 3 best periods for multiples test */
#define NUMBER_OF_PERIODS_INTO_SORTMULTIPLES \
    5 /* for compatibility with bestperiodslist */

/* alignment options */
#define GLOBAL 0
#define LOCAL 1

/* control routine status codes */
#define CTRL_SUCCESS 0
#define CTRL_BADFNAME -1
#define CTRL_BADFORMAT -2
#define CTRL_NOTHINGPROCESSED -3

#define MAXSEQNAMELEN 200

#define max4( a, b, c, d )                                                 \
    ( ( a >= b )                                                           \
        ? ( ( a >= c ) ? ( ( a >= d ) ? a : d ) : ( ( c >= d ) ? c : d ) ) \
        : ( ( b >= c ) ? ( ( b >= d ) ? b : d ) : ( ( c >= d ) ? c : d ) ) )
/* returns max of 4 in order a,b,c,d */

#define max3( a, b, c ) \
    ( ( a >= b ) ? ( ( a >= c ) ? a : c ) : ( ( b >= c ) ? b : c ) )
/* returns max of 3 in order a,b,c */

#define max2( a, b ) ( ( a >= b ) ? a : b )
/* returns max of 2 in order a,b */

#define max( a, b ) ( ( ( a ) >= ( b ) ) ? ( a ) : ( b ) )
#define min( a, b ) ( ( ( a ) <= ( b ) ) ? ( a ) : ( b ) )

#endif
