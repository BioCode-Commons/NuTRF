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

/* Windows of k-tuple matches for each distance, and the links between active distances. */

#ifndef TRF_DETECTION_DISTANCE_LIST_H
#define TRF_DETECTION_DISTANCE_LIST_H

#include "core/types.h"

struct distancelist *new_distancelist(void);
void clear_distancelist(struct distancelist *objptr);
void init_links(void);
void add_tuple_match_to_Distance_entry(int location, int size, int d,
		struct distancelist *objptr);
void link_Distance_window(int d);
void untag_Distance_window(int d, int linkdown);
int no_matches_so_unlink_Distance(int d, int location,
		struct distancelist *objptr);

#endif
