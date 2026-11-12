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

/* Resolution of the repeat list of one sequence into its final output. */

#ifndef TRF_REDUNDANCY_CLEAN_H
#define TRF_REDUNDANCY_CLEAN_H

void TRFClean( char * datafile, char* alignmentfile, char* tablefile,
		int maxsize, int data, int masked, char* Sequence,
		char* maskfile);

#endif
