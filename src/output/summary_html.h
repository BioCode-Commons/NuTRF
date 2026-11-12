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

/* Summary page that lists the sequences of a multi-sequence input. */

#ifndef TRF_OUTPUT_SUMMARY_HTML_H
#define TRF_OUTPUT_SUMMARY_HTML_H

#include <stdio.h>

void    WriteSummaryHeading(FILE* fp);
void    WriteSummaryRow(FILE* fp, int ordinal, char* tablefile, char* name, int repeats);
void    WriteSummaryEnd(FILE* fp, int foundsome);

#endif
