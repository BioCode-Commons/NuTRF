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

/* Program state shared by all modules. Every variable is defined in globals.c. */

#ifndef TRF_CORE_GLOBALS_H
#define TRF_CORE_GLOBALS_H

#include <stdio.h>
#include <stddef.h> /* has size_t definition */
#include <stdlib.h> /* has calloc definition */
#include <ctype.h>  /* includes toupper(c) */
#include <string.h> /* includes strncat() */
#include <math.h>   /* for ceil function */
#include <stdarg.h> /* for trf_message() function */
#include "core/platform.h"
#include "core/version.h"
#include "core/constants.h"
#include "core/types.h"

extern int counterInSeq;  /* number of repeats recorded for the current sequence */
extern unsigned int maxwraplength;  /* rows available to one wraparound alignment */
extern int Min_Distance_Entries;  /* minimum number of places to store a tuple match; larger than the distance for small distances */
extern int Min_Distance_Window;  /* minimum size of distance window; larger than the distance for small distances */
extern int PM;  /* match probability, percent */
extern int PI;  /* indel probability, percent */
extern double Pindel;  /* expected probability of a single character indel in the worst case tandem repeat */
extern int MAXDISTANCE;
extern int MAXPATTERNSIZE;
extern char debugbuffer[500];
extern unsigned char EC[2 * ( MAXPATTERNSIZECONSTANT + 1 )];  /* pattern used by the current alignment */
extern int *Index;  /* maps A, C, G, T to 0..3 */
extern int *ACGTcount;
extern unsigned char *Sequence;  /* current sequence, 1-based */
extern int Length;
extern int Delta;  /* indel penalty */
extern int Alpha;  /* match bonus */
extern int Beta;  /* mismatch penalty */
extern int AFDelta;  /* affine gap initiation penalty */
extern int AFGamma;  /* affine gap extension penalty */
extern int pwidth;
extern int Reportmin;
extern int Heading;
extern int Classlength;
extern int Test;
extern double Rows;
extern double Totalcharacters;
extern int Wrapend;
extern int Maxrealrow;
extern int Maxrow;
extern int Maxcol;
extern int Maxscore;
extern int ConsClasslength;
extern int *Tag;  /* list of tags for linking active distances */
extern int Toptag;  /* last tag in list */
extern struct pairalign AlignPair;
extern struct cons_data Consensus;
extern struct bestperiodlistelement Bestperiodlist[1];
extern struct distancelist *Distance;
extern struct distanceentry *_DistanceEntries;  /* storage shared by all Distance[].entry arrays */
extern struct distanceseenarrayelement *Distanceseenarray;  /* extent of the last alignment at each distance; blocks a repeated alignment in the same region */
extern struct distancelistelement Distanceseenlist[1];
extern int *SM;  /* similarity matrix, 256 x 256; written only by init_sm */
extern double Copynumber;
extern double WDPcount;
extern double OUTPUTcount;
extern int *Criteria_count;
extern int *Consensus_count;
extern int *Outputsize_count;
extern double *Cell_count;
extern double Try_waiting_time_count;
extern double Fail_waiting_time_count;
extern double Cell_total;
extern double Wasted_total;
extern TRFPARAMSET paramset;  /* this global controls the algorithm */
extern int *Bandcenter;  /* band center column of each narrowband row */
extern int **S;  /* alignment score matrix */
extern int Up[MAXBANDWIDTH + 1];
extern int Diag[MAXBANDWIDTH + 1];
extern int Maxrealcol;
extern int four_to_the[];
extern int NTS;  /* number of different tuple sizes to use; preset for all distances */
extern int Tuplesize[MAXTUPLESIZES + 1];
extern int Tuplemaxdistance[MAXTUPLESIZES + 1];  /* upper distance for each tuple size */
extern int Tuplecode[MAXTUPLESIZES + 1];  /* tuple codes encountered at a sequence location */
extern int *Tuplehash[MAXTUPLESIZES + 1];  /* points to last location of code in history list */
extern int Historysize[MAXTUPLESIZES + 1];  /* size of history lists */
extern int Nextfreehistoryindex[MAXTUPLESIZES + 1];  /* next free location in history index */
extern struct historyentry *History[MAXTUPLESIZES + 1];
extern int Criteria_print;
extern int Meet_criteria_print;
extern int *Sortmultiples;
extern int F[MAXBANDWIDTH + 1];
extern int Fdistance[MAXBANDWIDTH + 1];
extern int ldong;
extern int *Statistics_Distance;
extern FILE *Fptxt;  /* alignment file of the current sequence */
extern FILE *Fpdat;
extern int Minsize;
extern int Minscore;
extern int MaxPeriod;
extern int Period;
extern int print_flanking;
extern IL *GlobalIndexList;  /* repeats found in the current sequence */
extern IL *GlobalIndexListTail;
extern char hsequence[256];  /* heading lines of the HTML tables */
extern char hparameters[256];
extern char hlength[256];

#endif
