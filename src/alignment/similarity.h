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

/* Match and mismatch scores for every pair of byte values. */

#ifndef TRF_ALIGNMENT_SIMILARITY_H
#define TRF_ALIGNMENT_SIMILARITY_H

#include "core/globals.h"

void init_sm(int match, int mismatch);

/* returns the match or mismatch score; N does not match itself */
#define match( a, b ) ( SM[256 * ( ( a ) ) + ( b )] )

#endif
