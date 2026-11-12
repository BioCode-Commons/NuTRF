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

/* Removal of repeats by period size and of repeats made redundant by others. */

#ifndef TRF_REDUNDANCY_REDUNDANCY_H
#define TRF_REDUNDANCY_REDUNDANCY_H

#include "core/types.h"

IL*     RemoveBySize(IL * headptr, int maxsize);
IL*     RemoveRedundancy(IL * headptr);
int     IntervalOverlap(IL* iptr, IL* jptr);
int     IsRedundant(IL* iptr, IL* jptr);

#endif
