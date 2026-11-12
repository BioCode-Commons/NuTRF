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

/* Shared structure types. */

#ifndef TRF_CORE_TYPES_H
#define TRF_CORE_TYPES_H

#include <stdio.h>
#include "core/platform.h"
#include "core/constants.h"

/* alignment of the current primary and secondary sequences */
struct pairalign {
    int   length;
    int   score;
    char *textprime, *textsecnd;
    int * indexprime, *indexsecnd;
};

struct cons_data {
    char pattern[2 * ( MAXPATTERNSIZECONSTANT + 1 )];
    int  A[2 * ( MAXPATTERNSIZECONSTANT + 1 )],
      C[2 * ( MAXPATTERNSIZECONSTANT + 1 )],
      G[2 * ( MAXPATTERNSIZECONSTANT + 1 )],
      T[2 * ( MAXPATTERNSIZECONSTANT + 1 )],
      dash[2 * ( MAXPATTERNSIZECONSTANT + 1 )],
      insert[2 * ( MAXPATTERNSIZECONSTANT + 1 )],
      letters[2 * ( MAXPATTERNSIZECONSTANT + 1 )],
      total[2 * ( MAXPATTERNSIZECONSTANT + 1 )];
};

/*******************************************/
/*** new program started 11-29-95 **********/

struct bestperiodlistelement {
    int                           indexhigh;
    int                           indexlow;
    int                           best1;
    int                           best2;
    int                           best3;
    int                           best4;
    int                           best5;
    struct bestperiodlistelement *next;
};

struct distanceentry {
    int location;
    int size;
};

struct distancelist {
    int k_run_sums_criteria, waiting_time_criteria, lo_d_range, hi_d_range;
    int numentries, nummatches;
    int lowindex, highindex;
    int linked;
    int linkdown, linkup;
    struct distanceentry *entry;
};

/* created 5/23/05 G. Benson */
/* this array replaces the distanceseenlist.  It stores the extent of
   alignments, both pre- or post-consensus, and is used to block an
   alignment in the same region with the same distance */

struct distanceseenarrayelement {
    int index;
    int end;
    int score;
};

struct distancelistelement {
    int index;
    int distance;
    int changed_from_distance; /* use for test in
                                  search_for_distance_match_in_distanceseenlist
                                  3/10/05 */
    int end;
    int score;
    int best_possible_score; /* number of copies X length X match weight */
    int accepted;
    struct distancelistelement *next;
};

typedef struct {
    unsigned int match;
    unsigned int mismatch;
    unsigned int indel;
    unsigned int minscore;
    unsigned int maxperiod;
    unsigned int PM;
    unsigned int PI;
    int          datafile;
    int          maskedfile;
    int          flankingsequence;
    unsigned int flankinglength;
    int          HTMLoff;
    int          redundoff;
    int          ngs;
    int          use_stdin;
    unsigned int maxwraplength;

    char  inputfilename[_MAX_PATH]; /* constant defined in stdlib */
    char  outputprefix[_MAX_PATH];
    char  outputdirectory[_MAX_PATH];
    char  outputfilename[_MAX_PATH];
    int   multisequencefile; /* flags if file has more than one sequence */
    int   sequenceordinal;   /* holds seq. index starting on 1 */
    int   outputcount;       /* repeats found */
    int   guihandle;         /* this variable is only used in the GUI version */
    int   running;
    char *endstatus;
    int   percent;
} TRFPARAMSET;

struct historyentry {
    int location, previous, code;
};

struct distribution_parameters {
    double exp;
    double var;
};

struct MDDtype {
    int   distance;
    char *direction;
}; /* MDD[MAXWRAPLENGTH+1][MAXBANDWIDTH+1];*/

/* the following structure is used to pass a sequence to the algorithm */
typedef struct {
    /* Changed to unsigned Feb 16, 2016 Yozen */
    unsigned int length;
    int          composition[26];
    int          nucleotides;
    char         name[MAXSEQNAMELEN];
    char *       sequence;

} FASTASEQUENCE;

/* one reported repeat; the list is kept in order of detection */
struct index_list
{
	int     count;      /* indicates order in original file */
	char    ref[45];    /* records label for linking */
	int     first;      /* first index */
	int     last;       /* last index */
	int     period;     /* period size */
	float   copies;     /* number of copies */
	int     size;  /* consensus size */
	int     matches;
	int     indels;
	int     score;
	int     acount;
	int     ccount;
	int     gcount;
	int     tcount;
	float   entropy;
	char*   pattern;

	struct index_list* next;
};
typedef struct index_list IL;

#endif
