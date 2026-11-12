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

/* Definitions of the shared program state declared in globals.h. */

#include "core/globals.h"

int counterInSeq = 0;
unsigned int maxwraplength = 0;
int Min_Distance_Entries = 20;
int Min_Distance_Window = 20;
int PM;
int PI;
double Pindel;
int MAXDISTANCE = 500;
int MAXPATTERNSIZE = 500;
char debugbuffer[500];
unsigned char EC[2 * ( MAXPATTERNSIZECONSTANT + 1 )];
int *Index;
int *ACGTcount;
unsigned char *Sequence;
int Length;
int Delta;
int Alpha;
int Beta;
int AFDelta;
int AFGamma;
int pwidth = 75;
int Reportmin;
int Heading;
int Classlength;
int Test;
double Rows;
double Totalcharacters;
int Wrapend;
int Maxrealrow;
int Maxrow;
int Maxcol;
int Maxscore;
int ConsClasslength;
int *Tag;
int Toptag;
struct pairalign AlignPair;
struct cons_data Consensus;
struct bestperiodlistelement Bestperiodlist[1];
struct distancelist *Distance;
struct distanceentry *_DistanceEntries;
struct distanceseenarrayelement *Distanceseenarray;
struct distancelistelement Distanceseenlist[1];
int *SM = NULL;
double Copynumber;
double WDPcount;
double OUTPUTcount;
int *Criteria_count;
int *Consensus_count;
int *Outputsize_count;
double *Cell_count;
double Try_waiting_time_count;
double Fail_waiting_time_count;
double Cell_total;
double Wasted_total;
TRFPARAMSET paramset;
int *Bandcenter = NULL;
int **S;
int Up[MAXBANDWIDTH + 1];
int Diag[MAXBANDWIDTH + 1];
int Maxrealcol;
int four_to_the[] = {
  1, 4, 16, 64, 256, 1024, 4096, 16384, 65536, 262144, 1048576 };
int NTS;
int Tuplesize[MAXTUPLESIZES + 1];
int Tuplemaxdistance[MAXTUPLESIZES + 1];
int Tuplecode[MAXTUPLESIZES + 1];
int *Tuplehash[MAXTUPLESIZES + 1];
int Historysize[MAXTUPLESIZES + 1];
int Nextfreehistoryindex[MAXTUPLESIZES + 1];
struct historyentry *History[MAXTUPLESIZES + 1];
int Criteria_print = 0;
int Meet_criteria_print = 0;
int *Sortmultiples;
int F[MAXBANDWIDTH + 1];
int Fdistance[MAXBANDWIDTH + 1];
int ldong;
int *Statistics_Distance;
FILE *Fptxt;
FILE *Fpdat;
int Minsize = 1;
int Minscore;
int MaxPeriod;
int Period;
int print_flanking = 0;
IL *GlobalIndexList = NULL;
IL *GlobalIndexListTail = NULL;
char hsequence[256];
char hparameters[256];
char hlength[256];
