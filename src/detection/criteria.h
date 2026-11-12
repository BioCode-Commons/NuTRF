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

/* Statistical criteria that decide whether a candidate distance is aligned and reported. */

#ifndef TRF_DETECTION_CRITERIA_H
#define TRF_DETECTION_CRITERIA_H

#include <math.h>
#include "core/globals.h"

/* half-width of the range of distances that a repeat of period d spreads over */
static inline int d_range(int d)
{
	return((int)floor(2.3*sqrt(Pindel*d)));
}

int GetTopPeriods(unsigned char *pattern, int length, int *toparray);
int multiples_criteria_4(int found_d);
int new_meet_criteria_3(int d, int location, int tuplesize);

#endif
