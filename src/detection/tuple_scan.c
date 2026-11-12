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

/* k-tuple scan of one sequence; drives alignment and reporting of candidates. */

#include "core/globals.h"
#include "core/report.h"
#include "alignment/similarity.h"
#include "alignment/newwrap.h"
#include "alignment/narrowband.h"
#include "alignment/pair_alignment.h"
#include "consensus/consensus.h"
#include "consensus/statistics.h"
#include "detection/criteria.h"
#include "detection/distance_list.h"
#include "detection/distance_seen.h"
#include "detection/tuple_scan.h"
#include "periods/best_period_list.h"

/* Exact, division-free form of a hot remainder. It falls back to the
 * original expression outside the range in which the shortcut is exact. */
static inline int trf_mod_power_of_two(int value, int modulus)
{
	/* value%modulus; four_to_the[] holds powers of two */
	if(value>=0 && modulus>0 && (modulus&(modulus-1))==0) return value&(modulus-1);
	return value%modulus;
}

/* strchr("acgtACGT",c)!=NULL for every byte value, including the terminator */
static unsigned char trf_acgt_table[256];
static int trf_acgt_table_ready = 0;
static void trf_init_acgt_table(void)
{
	int c;
	for(c=0;c<256;c++) trf_acgt_table[c]=(strchr("acgtACGT",c)!=NULL);
	trf_acgt_table_ready=1;
}

/*******************************************************************/

/****************************  newtupbo()  *****************************/

/*******************************************************************/

/*??*/

/* last update May 19, 1997 */
/* started may 13, 1997 */
/* uses different tuple sizes for different distances */

void newtupbo(void)
{
	int mintuplesize,maxtuplesize,build_entire_code,g,badcharindex;
	int code,y,i,h,d,yy,j;
	int found,progbarpos,percentincrease,onepercent;
	double oldcopynumber;
	int pass_multiples_test;

	/* Moved here by Yozen on Feb 16, 2016.
	Always initialize here. Freed at the end of this function.
	maxwraplength has already been initialized in the calling
	function, TRF */
	Bandcenter = calloc(maxwraplength+1, sizeof(*Bandcenter));

	/* Jan 27, 2006, Gelfand, changed to use Similarity Matrix to avoid N matching itself */
	/* This function may be called multiple times (for different match/mismatch scores) */
	init_sm(Alpha,Beta);
	if(!trf_acgt_table_ready) trf_init_acgt_table();


	/* set progress indicator to zero  */
	paramset.percent=0;
	if (paramset.ngs != 1) SetProgressBar();


	OUTPUTcount = 0; /* needed to make browser label unique */

	mintuplesize=Tuplesize[1];
	maxtuplesize=Tuplesize[NTS];
	for(g=1;g<=NTS;g++)
	{
		Tuplehash[g]=(int *)farcalloc(four_to_the[Tuplesize[g]],sizeof(int));
		Historysize[g]=2*(Tuplemaxdistance[g]+1)+2; /* The idea here is that no */
		/* previous history pointer points back */
		/* more than Tuplemaxdistance.  Then, when */
		/* History entry is reused, following */
		/* links from the current will exceed the */
		/* maxdistance before reaching the reused */
		/* entry. */
		History[g]=(struct historyentry *)
			farcalloc(Historysize[g], sizeof(struct historyentry));
		Nextfreehistoryindex[g]=1;  /* set all to 1 because 0 indicates */
		/* Tuplehash points to nothing */

	}

	Sortmultiples=(int *)calloc(MAXDISTANCE+1,sizeof(int));

	build_entire_code=1;

	onepercent = Length/100;
	percentincrease = 0;
	progbarpos = 0;
	for(i=0;i<=Length;i++)
	{
		/* if percent changed then set indicator */
		percentincrease++;
		if(percentincrease==onepercent)
		{
			percentincrease = 0;
			progbarpos++;
			paramset.percent = progbarpos;
			if (paramset.ngs != 1) SetProgressBar();
		}
		if((i==0)  /* before start of sequence or */
				||(!trf_acgt_table[Sequence[i]]))  /* not one of A,C,G,T */
		{
			badcharindex=i;
			build_entire_code=1;
			/* find first good string of mintupsize characters */
			g=0;
			while((g<mintuplesize)&&(i<Length))
			{
				i++;
				if(!trf_acgt_table[Sequence[i]])
				{
					badcharindex=i;
					g=0;
				}
				else
					g++;
			}
			if(g<mintuplesize) break; /* i=Length and minimum tuple not found */
		}
		if(build_entire_code)
		{
			code=0;
			for(g=badcharindex+1;g<=i;g++)
			{
				code=code*4+Index[Sequence[g]];
			}
			if(i-badcharindex>=maxtuplesize){build_entire_code=0;}
		}
		else
		{
			code=trf_mod_power_of_two(code,four_to_the[Tuplesize[NTS]-1])
				*4+Index[Sequence[i]];
		}
		Tuplecode[NTS]=code;
		for(h=NTS-1;h>=1;h--)
		{
			Tuplecode[h]=trf_mod_power_of_two(code,four_to_the[Tuplesize[h]]);
		}

		/* process index i using all the tuplesizes */

		g=1;
		while((g<=NTS)&&(i-badcharindex>=Tuplesize[g]))
		{

			/* change 5/25/99 ends here */

			y=Tuplehash[g][Tuplecode[g]]; /* index in history list of last
							 occurrence of code */
			h=Nextfreehistoryindex[g]; /* next free index in history list */
			j=h+1;            /* advance next free index */
			if(j==Historysize[g]) j=1; /* we use a circular history list */

			if((History[g][j].location!=0) /* if the next entry has 
							  already been used */
					&&(j==Tuplehash[g][History[g][j].code])) /* check Tuplehash.  
										    If it still */
			{             /* points here, */
				Tuplehash[g][History[g][j].code]=0; /* zero it out.   */

			}
			Nextfreehistoryindex[g]=j;

			Tuplehash[g][Tuplecode[g]]=h; /* store index of current tuple */
			History[g][h].location=i; /* store info about current tuple */
			History[g][h].previous=y;
			History[g][h].code=Tuplecode[g];

			yy=h;         /* yy holds entry which points to y */
			while(y!=0)
			{
				d=i-History[g][y].location; /* d=distance between matching 
							       tuples */
				if(d>Tuplemaxdistance[g]) /* if d exceeds Tuplemaxdistance, 
							     then */
				{           /* make the previous location 0.  We */
					History[g][yy].previous=0; /* are no longer interested 
								      in the y */
					y=0;          /* entry.  It will be zeroed out when */
				} /* reused */
				else
				{
					yy=y;
					y=History[g][y].previous; /* get next matching tuple */
					/* process */
					/* is this a distance that is too small for the tuplesize? */
					/* recall that the History lists do not exclude distance */
					/* that are too short, only those that are too long  */
					if(d>Tuplemaxdistance[g-1])
					{
						/* add tuplematch to Distance */
						add_tuple_match_to_Distance_entry(i,Tuplesize[g],d,Distance);


						/* check if this distance has already been processed */

						/* found=search_for_distance_match_in_distanceseenlist(d,i); */
						/* modified 5/23/05 G. Benson */
						found=search_for_distance_match_in_distanceseenarray(d,i);


						if(!found)
						{


							/* is distance d linked into other nonzero distances? */
							/* if not, link it in */
							if(!Distance[d].linked) link_Distance_window(d);

							/* test criteria for candidate */


							if((new_meet_criteria_3(d,i,Tuplesize[g]))
									/* change 2: changed 250 into 100 */
									&&((d<=250)||(search_for_range_in_bestperiodlist(i,d)))) 
								/* use bestperiod list only for distances greater than 500 */           
							{
								/* align sequence against candidate */
								/* and get alignment */
								WDPcount++;
								Criteria_count[d]++;
								Rows=0;
								if(d<=SMALLDISTANCE)
								{
									newwrap(i,d,WITHOUTCONSENSUS);
									Cell_count[d]+=(Rows*d);
									get_pair_alignment_with_copynumber(d);
								}
								else        /* d is a large distance */
								{
									if (1%100==0) fprintf(stderr,"\ni=%d  d=%d",i,d);
									narrowbandwrap(i,d,max(MINBANDRADIUS,d_range(d)), min(2*max(MINBANDRADIUS,d_range(d)), (d/3) ), 
											WITHOUTCONSENSUS,RECENTERCRITERION);
									Cell_count[d]+=(Rows*(2*max(MINBANDRADIUS,d_range(d))+1));
									get_narrowband_pair_alignment_with_copynumber( d, min(2*max(MINBANDRADIUS,d_range(d)), (d/3) ), LOCAL );
								} 
								if(Meet_criteria_print)trf_message("\nFrom:%d, To:%d,  Copynumber:%f",
										AlignPair.indexprime[AlignPair.length],
										AlignPair.indexprime[1],Copynumber);
								/* add_to_distanceseenlist(i,d,Maxrealrow,Maxscore,WITHOUTCONSENSUS); */
								/* modified 5/23/05 G. Benson */
								add_to_distanceseenarray(i,d,Maxrealrow,Maxscore);




								/* changed 2/17/05 gary benson
								   change to make number of copies required less restrictive for pattern sizes >= 50
								   this ramps from 1.9 at pattern size = 50 down to 1.8 at pattern size = 100 or above */
								if((d<=50 && Copynumber<1.9)
										|| (d>50 && d<=100 && Copynumber < 1.9 - 0.002*(d-50))
										|| (d>100 && Copynumber < 1.8))
									/* if(Copynumber<1.9) */
									/* max(1.9,(double)Min_Distance_Window/(double)d)) */
									/*  ceil((double)Min_Distance_Window/(double)d)))*/
									/* rationale for max here; Min_Distance_Window is the */
									/* smallest actual tandem repeat we want to see, so */
									/* divide by d to get minimum number of copies */ 
								{

								}
								else
								{

									pass_multiples_test=multiples_criteria_4(d);
									add_to_bestperiodlist(d);
									if(pass_multiples_test)
									{

										/* get consensus */
										Classlength=d;
										get_consensus(d);
										if(ConsClasslength!=Classlength)
											Classlength=ConsClasslength;
										oldcopynumber=Copynumber;
										/* repeat alignment using consensus */
										Consensus_count[Classlength]++;
										Rows=0;
										if(Classlength<=SMALLDISTANCE)
										{
											newwrap(i,Classlength,WITHCONSENSUS);
											Cell_count[Classlength]+=(Rows*d);
											get_pair_alignment_with_copynumber(Classlength);
										}
										else    /* d is a large distance */
										{
											narrowbandwrap(i,Classlength,
													max(MINBANDRADIUS,d_range(Classlength)), min(2*max(MINBANDRADIUS,d_range(Classlength)), (Classlength/3) ), 
													WITHCONSENSUS,RECENTERCRITERION);
											Cell_count[Classlength]+=(Rows*(2*max(MINBANDRADIUS,
															d_range(Classlength))+1));
											get_narrowband_pair_alignment_with_copynumber( Classlength, min(2*max(MINBANDRADIUS,d_range(Classlength)), (Classlength/3) ), LOCAL );
										}
										/*add_to_distanceseenlist(i,Classlength,Maxrealrow,
										  Maxscore,WITHCONSENSUS); */
										/* modified 5/23/05 G. Benson */
										add_to_distanceseenarray(i,d,Maxrealrow,Maxscore);

										/* modified 6/10/05 G. Benson */
										adjust_bestperiod_entry(d);


										/* changed 2/17/05 gary
										   change to make number of copies required less restrictive for pattern sizes >= 50
										   this ramps from 1.9 at pattern size = 50 down to 1.8 at pattern size = 100 or above */
										if((Classlength<=50 && Copynumber<1.9)
												|| (Classlength>50 && Classlength<=100 && Copynumber < 1.9 - 0.002*(d-50))
												|| (Classlength>100 && Copynumber < 1.8))

											/* if(Copynumber<1.9) */
											/*  max(1.9,(double)Min_Distance_Window/(double)Classlength)) */
										{
										}
										else
										{
											if((Classlength>=Minsize)&&(AlignPair.score>=Minscore))
											{

												/* output repeat */
												OUTPUTcount++;          Outputsize_count[Classlength]++;
												/* sequence */
												trf_message("\nFound at i:%d original size:%d final size:%d",
														i,d,Classlength);
												get_statistics(Classlength);
											}


										}

									}



								}
							}
						}
					}
				}
			}
			g++;
		}
	}

	/* close progress indicator */
	paramset.percent = -1;
	if (paramset.ngs != 1) SetProgressBar();

	/* Added by Yozen on Feb 16, 2016. */
	free(Bandcenter);
}
