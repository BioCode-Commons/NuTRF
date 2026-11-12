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

/* Extent of earlier alignments at each distance. */

#ifndef TRF_DETECTION_DISTANCE_SEEN_H
#define TRF_DETECTION_DISTANCE_SEEN_H

void init_distanceseenarray(void);
void free_distanceseenarray(void);
void add_to_distanceseenarray(int location,int distance,int end,int score);
int search_for_distance_match_in_distanceseenarray(int distance,int start);

void init_distanceseenlist(void);
void free_distanceseenlist(void);
void add_to_distanceseenlist(int location,int distance,int end,int  score, int acceptstatus);
int search_for_distance_match_in_distanceseenlist(int distance,int start);

#endif
