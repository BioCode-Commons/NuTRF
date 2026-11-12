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

/* Best periods of recent alignments, used to filter large candidate distances. */

#ifndef TRF_PERIODS_BEST_PERIOD_LIST_H
#define TRF_PERIODS_BEST_PERIOD_LIST_H

void init_bestperiodlist(void);
void free_bestperiodlist(void);
void add_to_bestperiodlist(int d);
void adjust_bestperiod_entry(int d);
int search_for_range_in_bestperiodlist(int start, int distance);

#endif
