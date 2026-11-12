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

/* Analysis of one sequence: allocation, scan, and resolution of the output. */

#include <errno.h>
#include "core/globals.h"
#include "core/report.h"
#include "core/arrays.h"
#include "detection/coin_toss.h"
#include "detection/criteria.h"
#include "detection/distance_list.h"
#include "detection/distance_seen.h"
#include "detection/tuple_scan.h"
#include "output/html.h"
#include "periods/best_period_list.h"
#include "redundancy/clean.h"
#include "run/sequence.h"

/*************************************************
 *   This routine acts on single-sequence files and
 *   is used by the control routine.
 **************************************************/
void TRF(FASTASEQUENCE* pseq)
{
	unsigned int i;  /* used at the end to free memory */
	int *stemp;
	char htmlstring[_MAX_PATH], txtstring[_MAX_PATH],
	     paramstring[_MAX_PATH],datstring[_MAX_PATH],
	     maskstring[_MAX_PATH],messagebuffer[100];

	/* added Nov 1, 2012 by Y. Gelfand */
	init_bestperiodlist();


	/*  Set global print_flanking that controls the generation of flanking */
	print_flanking = paramset.flankingsequence;

	/* allocate memory for file names */
	if (paramset.ngs != 1) PrintProgress("Allocating Memory...");


	/* change made for NGS data analysis */
	/* G. Benson */
	/* 1.26.10 */
	/* make MAXWRAPLENGTH = 1000 for smaller for small sequences */
	maxwraplength = min(paramset.maxwraplength,pseq->length);

	/* allocate memory */
	S = (int **) malloc((maxwraplength+1) *  sizeof(int *));
	if(S==NULL)
	{
		PrintError("Unable to allocate memory for S array");
		exit(-1);
	}
	// stemp = (int *) malloc(((MAXWRAPLENGTH+1)*(MAXBANDWIDTH+1)) * sizeof(int));
	// memset(stemp,0,((MAXWRAPLENGTH+1)*(MAXBANDWIDTH+1)) * sizeof(int));
	/* Yozen Jan 26, 2016: We control the compilation and we're going to be using C99
	or greater standard C; don't need to cast, and we can use the pointer to determine the size.
	Also, use calloc instead of malloc+memset. */
	stemp = calloc(((size_t)(maxwraplength+1)*(MAXBANDWIDTH+1)), sizeof(*stemp));
	if(stemp==NULL)
	{
		char errmsg[255];
		#if (__x86_64__ + __x86_64 + __amd64 + __amd64__ + _M_AMD64 + _M_X64) > 1
		snprintf(errmsg, 255, "Unable to allocate %lu bytes for stemp array. Please set a lower value for the longest TR length. (%s:%d)\n", ((maxwraplength+1)*(MAXBANDWIDTH+1)) * sizeof(*stemp), __FILE__, __LINE__);
		#else
		snprintf(errmsg, 255, "Unable to allocate %u bytes for stemp array. Please set a lower value for the longest TR length. (%s:%d)\n", ((maxwraplength+1)*(MAXBANDWIDTH+1)) * sizeof(*stemp), __FILE__, __LINE__);
		#endif
		PrintError(errmsg);
		exit(-1);
	}
	for(i=0;i<=maxwraplength;i++)
	{
		S[i]= stemp;
		stemp+=MAXBANDWIDTH+1;
	}
	S[0][0]=1;

	/* AlignPair holds the characters and alignments of the current */
	/* primary and secondary sequences  */
	AlignPair.textprime=newAlignPairtext(2*maxwraplength);
	if (paramset.endstatus) {
		return;
	}
	AlignPair.textsecnd=newAlignPairtext(2*maxwraplength);
	if (paramset.endstatus) {
		return;
	}
	AlignPair.indexprime=newAlignPairindex(2*maxwraplength);
	if (paramset.endstatus) {
		return;
	}
	AlignPair.indexsecnd=newAlignPairindex(2*maxwraplength);
	if (paramset.endstatus) {
		return;
	}

	/* set algorithm's parameters */
	Alpha = paramset.match;
	Beta  = -paramset.mismatch;
	Delta = -paramset.indel;
	PM    = paramset.PM;
	PI    = paramset.PI;
	Minscore = paramset.minscore;
	MaxPeriod= paramset.maxperiod;
	MAXDISTANCE=MAXPATTERNSIZE=paramset.maxperiod;
	if(MAXDISTANCE<500)
	{
		MAXDISTANCE = MAXPATTERNSIZE=500;
	}

	MAXDISTANCE = MAXPATTERNSIZE = min( MAXDISTANCE, (int) (pseq->length * .6));
	MAXDISTANCE = MAXPATTERNSIZE = max( MAXDISTANCE, 200);


	/* generate the parameter string to be used in file names */
	sprintf(paramstring,"%d.%d.%d.%d.%d.%d.%d",
			paramset.match,paramset.mismatch,paramset.indel,
			paramset.PM,paramset.PI,paramset.minscore,paramset.maxperiod);


	Reportmin= 0;
	ldong=0;
	Rows=0;
	Totalcharacters=0;
	Test=1;

	/* print the names of the files */
	sprintf(htmlstring,"%s.%s.html",paramset.outputprefix,paramstring);
	sprintf(txtstring ,"%s.%s.txt.html",paramset.outputprefix,paramstring);
	sprintf(datstring ,"%s.%s.dat",paramset.outputprefix,paramstring);
	sprintf(maskstring,"%s.%s.mask",paramset.outputprefix,paramstring);

	/* start txt file */
	if (!paramset.HTMLoff) {
		Fptxt=fopen(txtstring,"w");
		if(Fptxt==NULL)
		{
			PrintError("Unable to open alignment file for writing in TRF() routine!");
			exit(-1);
		}



		fprintf(Fptxt,"<HTML>");
		fprintf(Fptxt,"<HEAD>");
		fprintf(Fptxt,"<TITLE>%s</TITLE>",txtstring);
		fprintf(Fptxt,"</HEAD>");
		fprintf(Fptxt,"<BODY bgcolor=\"#FBF8BC\">");
		fprintf(Fptxt,"<PRE>");

		fprintf(Fptxt,"\nTandem Repeats Finder Program written by:");
		fprintf(Fptxt,"\n\n                 Gary Benson");
		fprintf(Fptxt,"\n      Program in Bioinformatics");
		fprintf(Fptxt,"\n          Boston University");
		fprintf(Fptxt,"\n\nVersion %s", versionstring);

		fprintf(Fptxt,"\n\nSequence: %s\n\nParameters: %d %d %d %d %d %d %d\n",
				pseq->name,paramset.match,paramset.mismatch,paramset.indel,
				paramset.PM,paramset.PI,paramset.minscore,paramset.maxperiod);
	}



	if (paramset.ngs != 1) PrintProgress("Initializing data structures...");



	Distance=new_distancelist();
	clear_distancelist(Distance);
	Tag=newTags(MAXDISTANCE/TAGSEP+1);
	Toptag=(int)ceil(MAXDISTANCE/TAGSEP);
	init_links();

	init_index();
	/* init_distanceseenlist(); */
	/* modified 5/23/05 G. Benson */
	init_distanceseenarray();


	if (paramset.ngs != 1) PrintProgress("Computing TR Model Statistics...");

	init_and_fill_coin_toss_stats2000_with_4tuplesizes();

	/* over allocate statistics_distance array to prevent spill in alignments
	   with execive insertion counts Jan 07, 2003 */
	Statistics_Distance=(int *)calloc(4*MAXDISTANCE,sizeof(int));
	if(Statistics_Distance==NULL)
	{
		PrintError("Unable to allocate memory for Statistics_Distance array");
		exit(-3);
	}



	/* set the sequence pointer. more global vars! */
	Sequence = pseq->sequence-1; /* start one character before */
	Length=pseq->length;

	if (!paramset.HTMLoff) {
		fprintf(Fptxt,"\n\nLength: %d",Length);
		fprintf(Fptxt,"\nACGTcount: A:%3.2f, C:%3.2f, G:%3.2f, T:%3.2f\n\n",
				(double)pseq->composition['A'-'A']/Length,
				(double)pseq->composition['C'-'A']/Length,
				(double)pseq->composition['G'-'A']/Length,
				(double)pseq->composition['T'-'A']/Length);

		if((pseq->length-pseq->nucleotides)>0)
		{
			fprintf(Fptxt,"Warning! %d characters in sequence are not A, C, G, or T\n\n",
					(pseq->length-pseq->nucleotides));
		}
	}

	Totalcharacters+=Length;
	WDPcount=0;

	/* G. Benson 1/28/2004 */
	/* following four memory allocations increased to avoid memory error when
	   consensus length exceeds MAXDISTANCE after returning from get_consensus(d) */

	Criteria_count=(int *)calloc(2*(MAXDISTANCE+1),sizeof(int));
	if(Criteria_count==NULL)
	{
		PrintError("Unable to allocate Criteria_count");
		exit(-4);
	}
	Consensus_count=(int *)calloc(2*(MAXDISTANCE+1),sizeof(int));
	if(Consensus_count==NULL)
	{
		PrintError("Unable to allocate memory for Consensus_count");
		exit(-5);
	}
	Cell_count=(double *)calloc(2*(MAXDISTANCE+1),sizeof(double));
	if(Cell_count==NULL)
	{
		PrintError("Unable to allocate memory for Cell_count");
		exit(-6);
	}
	Outputsize_count=(int *)calloc(2*(MAXDISTANCE+1),sizeof(int));
	if(Outputsize_count==NULL)
	{
		PrintError("Unable to allocate memory for Outputsize_count");
		exit(-7);
	}


	if(paramset.multisequencefile)
	{
		sprintf(messagebuffer,"Scanning Sequence %d...",
				paramset.sequenceordinal);
		if (paramset.ngs != 1) PrintProgress(messagebuffer);
	}
	else if (paramset.ngs != 1)
	{
		PrintProgress("Scanning...");
	}


	clear_distancelist(Distance);

	newtupbo(); /* this is the main function of the algorithm */


	Cell_total=0;
	for(i=MAXDISTANCE;i>=1;i--)
	{
		Cell_total+=Cell_count[i];
		if(i<=SMALLDISTANCE)
			Wasted_total+=(Criteria_count[i]+Consensus_count[i])*i*i*2;
		else
			Wasted_total+=(Criteria_count[i]+Consensus_count[i])*i
				*(2*d_range(i)+1)*2;
	}

	if (!paramset.HTMLoff) {
		fprintf(Fptxt,"Done.\n");
		fprintf(Fptxt,"</PRE>");
		fprintf(Fptxt,"</BODY>");
		fprintf(Fptxt,"</HTML>");
		fclose(Fptxt);
	}


	/****************************************************************
	 * The following memory deallocations where not originally
	 * implemented. There may still exist some marginal memory leaks
	 * but due messy code structure and extensive use of macros and
	 * global variables this could be impossible to fix.
	 *****************************************************************/

	if (paramset.ngs != 1) PrintProgress("Freeing Memory...");

	free(S[0]);
	free(S);
	free(Statistics_Distance);
	free(Criteria_count);
	free(Consensus_count);
	free(Cell_count);
	free(Outputsize_count);
	free(AlignPair.textprime);
	free(AlignPair.textsecnd);
	free(AlignPair.indexprime);
	free(AlignPair.indexsecnd);
	free(Tag);
	free(Index);
	/* free_distanceseenlist(); */
	/* modified 5/23/05 G. Benson */
	free_distanceseenarray();

	/* free distance list and all its entries */
	//for(i=1;i<=MAXDISTANCE;i++) free(Distance[i].entry);
	free(_DistanceEntries);
	free(Distance);

	for(i=1; i<=NTS; i++)
	{
		free(Tuplehash[i]);
		free(History[i]);
	}

	free(Sortmultiples);

	if (paramset.ngs != 1) PrintProgress("Resolving output...");

	/* sort, remove redundancy, and write the HTML output */
	TRFClean( datstring, txtstring, htmlstring,MaxPeriod,
			paramset.datafile,paramset.maskedfile,Sequence,maskstring);

	/* Set the name of the outputfilename global to name given to
	   file by the HTML output routines */
	MakeFileName(paramset.outputfilename,htmlstring,1);

	/* added Nov 1, 2012 by Y. Gelfand */
	free_bestperiodlist();

	if (paramset.ngs != 1) PrintProgress("Done.");

	return;
}
