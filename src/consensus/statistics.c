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

/* Statistics of one repeat; records it in the repeat list and the alignment file. */

#include "core/globals.h"
#include "core/report.h"
#include "core/repeat_list.h"
#include "alignment/pair_alignment.h"
#include "consensus/statistics.h"
#include "detection/criteria.h"
#include "output/alignment_text.h"

/*******************************************************************/

/***************************** statistics() **************************/

/*******************************************************************/

void get_statistics(int consensussize)
{
	int g,x,lp,rp,match,mismatch,indel,mindistance,maxdistance;
	int d;
	int best_match_distance,best_match_count;
	int i,count;
	double diversity[4];
	double  temp, entropy;
	int e,f;
	int size;
	int startECpos; /* to print consensus pattern */
	int ACGTcount[26];

	size=consensussize;


	for(g=1;g<=MAXDISTANCE+d_range(MAXDISTANCE);g++)
		Statistics_Distance[g]=0;

	match=0;
	mismatch=0;
	indel=0;



	lp=1;
	while(AlignPair.textsecnd[lp]=='-')
	{
		lp++;
		if(lp>AlignPair.length)
		{
			trf_message("\nError in statistics.");
			trf_message("\nInitial left pointer exceeds AlignPair.length while");
			trf_message("\nlooking for first non -");
			exit(-10); 
		}
	}
	rp=lp+1;
	if(rp>AlignPair.length)
	{
		trf_message("\nError in statistics.");
		trf_message("\nInitial right pointer exceeds AlignPair.length");
		exit(-10);
	} 
	while(AlignPair.indexsecnd[rp]!=AlignPair.indexsecnd[lp])
	{ 
		rp++;
		if(rp>AlignPair.length)
		{
			trf_message("\nError in statistics.");
			trf_message("\nInitial right pointer exceeds AlignPair.length while");
			trf_message("\nlooking for AlignPair.indexsecnd[lp]");
			exit(-10); 
		}
	}
	while(AlignPair.textsecnd[rp]=='-')
	{
		rp++;
		if(rp>AlignPair.length)
		{
			trf_message("\nError in statistics.");
			trf_message("\nInitial right pointer exceeds AlignPair.length while");
			trf_message("\nlooking for first non -");
			exit(-10); 
		}
	}
	if(AlignPair.indexsecnd[lp]!=AlignPair.indexsecnd[rp])
	{
		trf_message("\nError in statistics.");
		trf_message("\nInitial left pointer index not the same as");
		trf_message("\ninitial right pointer index");
		exit(-10);
	}

	/* now do the tests */
	mindistance=size;
	maxdistance=size;

	while((rp<=AlignPair.length)&&(lp<rp))
	{

		if((AlignPair.textsecnd[lp]!='-')&&(AlignPair.textsecnd[rp]!='-'))
		{
			if((AlignPair.textprime[lp]!='-')&&(AlignPair.textprime[rp]!='-'))
			{
				if(AlignPair.indexsecnd[lp]!=AlignPair.indexsecnd[rp])
				{
					trf_message("\nError in statistics.");
					trf_message("\nLeft pointer index not the same as right pointer index");
					trf_message("\non match or mismatch when nothing is dash.");
					trf_message("\nAlignPair.indexsecnd[%d]: %d, AlignPair.indexprime[%d]: %d"
							,lp,AlignPair.indexsecnd[lp],lp,AlignPair.indexprime[lp]);
					trf_message("\nAlignPair.indexsecnd[%d]: %d, AlignPair.indexprime[%d]: %d"
							,rp,AlignPair.indexsecnd[rp],rp,AlignPair.indexprime[rp]);

					exit(-10);
				}
				if(AlignPair.textprime[lp]==AlignPair.textprime[rp])
				{
					match++;
					d=AlignPair.indexprime[rp]-AlignPair.indexprime[lp];


					if(d<0)d=-d;

					/* protect for memory override January 08, 2003 */
					if(d<(4*MAXDISTANCE))
					{
						Statistics_Distance[d]++;
						if (d<mindistance)mindistance=d;
						if (d>maxdistance)maxdistance=d;



					}	  
				}
				else 
				{
					mismatch++;

				}
				lp++;
				rp++;


			}
			else if((AlignPair.textprime[lp]=='-')&&(AlignPair.textprime[rp]=='-'))
			{
				/* do nothing */
				if(AlignPair.indexsecnd[lp]!=AlignPair.indexsecnd[rp])
				{
					trf_message("\nError in statistics.");
					trf_message("\nLeft pointer index not the same as right pointer index");
					trf_message("\non do nothing.\nlp: %d,  rp: %d",lp,rp);
					exit(-10);
				}

				lp++;
				rp++;
			}
			else if((AlignPair.textprime[lp]=='-')||(AlignPair.textprime[rp]=='-'))
			{
				if(AlignPair.indexsecnd[lp]!=AlignPair.indexsecnd[rp])
				{
					trf_message("\nError in statistics.");
					trf_message("\nLeft pointer index not the same as right pointer index");
					trf_message("\non indel caused by dash in textprime.");
					trf_message("\nAlignPair.indexsecnd[%d]: %d, AlignPair.indexprime[%d]: %d"
							,lp,AlignPair.indexsecnd[lp],lp,AlignPair.indexprime[lp]);
					trf_message("\nAlignPair.indexsecnd[%d]: %d, AlignPair.indexprime[%d]: %d"
							,rp,AlignPair.indexsecnd[rp],rp,AlignPair.indexprime[rp]);

					exit(-10);
				}
				indel++;

				lp++;
				rp++;
			}
		}
		else if ((AlignPair.textsecnd[lp]=='-')&&(AlignPair.textsecnd[rp]=='-'))
		{
			if(AlignPair.indexsecnd[lp]!=AlignPair.indexsecnd[rp])
			{
				trf_message("\nError in statistics.");
				trf_message("\nLeft pointer index not the same as right pointer index");
				trf_message("\non match or mismatch when both secnd are dash.");
				trf_message("\nlp: %d,  rp: %d",lp,rp);
				exit(-10);
			}
			if (AlignPair.textprime[lp]==AlignPair.textprime[rp])
			{
				match++;
				d=AlignPair.indexprime[rp]-AlignPair.indexprime[lp];
				if(d<0) d=-d;

				/* protect for memory override January 08, 2003 */
				if(d<(4*MAXDISTANCE))
				{
					Statistics_Distance[d]++;
					if (d<mindistance)mindistance=d;
					if (d>maxdistance)maxdistance=d;
				}

			}
			else 
			{
				mismatch++;

			}
			lp++;
			rp++;
		}
		else if(AlignPair.textsecnd[lp]=='-')
		{

			indel++;
			lp++;
		}
		else if(AlignPair.textsecnd[rp]=='-') 
		{

			indel++;
			rp++;
		}
	}
	if(lp>=rp)
	{
		trf_message("\nError in statistics.");
		trf_message("\nLeft pointer >= right pointer");
		exit(-10); 
	}

	x=match+mismatch+indel;

	best_match_distance=0;
	best_match_count=0;
	for(g=mindistance;g<=maxdistance;g++)
		if(Statistics_Distance[g]!=0)
		{
			if(Statistics_Distance[g]>best_match_count)
			{
				best_match_count=Statistics_Distance[g];
				best_match_distance=g;
			}
		}


	reverse();
	Period=best_match_distance;

	if (!paramset.HTMLoff)
		print_alignment_headings(Classlength);

	/* save the starting position of consensus in EC */
	startECpos = (int) AlignPair.indexsecnd[1];

	if (!paramset.HTMLoff)
		alt3_print_alignment(Classlength);

	reverse();

	if (!paramset.HTMLoff) {	
		fprintf(Fptxt,"\nStatistics");
		fprintf(Fptxt,"\nMatches: %d,  Mismatches: %d, Indels: %d",match, mismatch, indel);
		fprintf(Fptxt,"\n        %0.2f            %0.2f        %0.2f",
				(float)match/x, (float)mismatch/x, (float)indel/x);
		fprintf(Fptxt,"\n");
		fprintf(Fptxt,"\nMatches are distributed among these distances:");
		for(g=mindistance;g<=maxdistance;g++)
			if(Statistics_Distance[g]!=0)
			{
				fprintf(Fptxt,"\n %3d  %3d  %0.2f",
						g,Statistics_Distance[g],(float)Statistics_Distance[g]/match);
			}
	}

	ACGTcount['A'-'A']=0;
	ACGTcount['C'-'A']=0;
	ACGTcount['G'-'A']=0;
	ACGTcount['T'-'A']=0;
	i=1;
	count=0;
	while (i<=AlignPair.length)
	{
		if(AlignPair.textprime[i]!='-')
		{
			ACGTcount[AlignPair.textprime[i]-'A']++;
			count++;
		}
		i++;
	}

	if(count!=0)
		trf_message("\n\nACGTcount: A:%3.2f, C:%3.2f, G:%3.2f, T:%3.2f",
				(double)ACGTcount['A'-'A']/count,(double)ACGTcount['C'-'A']/count,
				(double)ACGTcount['G'-'A']/count,(double)ACGTcount['T'-'A']/count); 
	else trf_message("\nError in statistics: ACGTcount=0");

	diversity[0]=(double)ACGTcount['A'-'A']/count;
	diversity[1]=(double)ACGTcount['C'-'A']/count;
	diversity[2]=(double)ACGTcount['G'-'A']/count;
	diversity[3]=(double)ACGTcount['T'-'A']/count;

	for(e=2;e>=0;e--)
		for(f=0;f<=e;f++)
		{
			if(diversity[f]<diversity[f+1])
			{
				temp=diversity[f+1];
				diversity[f+1]=diversity[f];
				diversity[f]=temp;
			}
		}

	entropy=(((diversity[0]==0)?0:(diversity[0]*(log(diversity[0])/log(2))))+
			((diversity[1]==0)?0:(diversity[1]*(log(diversity[1])/log(2))))+
			((diversity[2]==0)?0:(diversity[2]*(log(diversity[2])/log(2))))+
			((diversity[3]==0)?0:(diversity[3]*(log(diversity[3])/log(2)))));
	if(entropy<0)entropy=-entropy;

	if (!paramset.HTMLoff) {
		fprintf(Fptxt,"\n");
		fprintf(Fptxt,"\n");
	}

	/*  changed by Gary Benson, 6/1/99, to remove summary
	    information and add flanking sequence */

	/*  fprintf(Fptxt,"\n\n %d--%d & %d & %3.1f & %d & %d & %d & %d & %d & %d & %d & %d & %3.2f \\\\\n\n",
	    AlignPair.indexprime[AlignPair.length],
	    AlignPair.indexprime[1],best_match_distance,Copynumber,Classlength,
	    (int)(100*(float)match/x),(int)(100*(float)indel/x),AlignPair.score,
	    (int)(100*(double)ACGTcount['A'-'A']/count),
	    (int)(100*(double)ACGTcount['C'-'A']/count), 
	    (int)(100*(double)ACGTcount['G'-'A']/count),
	    (int)(100*(double)ACGTcount['T'-'A']/count),
	    entropy); 
	 */

	/* prints line showing the consensus pattern */

	if (!paramset.HTMLoff) {

		printECtoAlignments(Fptxt, startECpos, consensussize);


		if(print_flanking)
		{
			reverse();
			print_flanking_sequence(paramset.flankinglength);
			reverse();
		}       

	}
	/* end of change, 6/1/99 */

	/* Added by Yevgeniy Gelfand on Jan 27, 2010  */
	/* To have smaller sequences not send results */
	/* to disc to improve performance             */
	{

		IL      *newptr;

		/* create new index list element */
		newptr = (IL*) malloc(sizeof(IL));
		if(newptr==NULL)
		{
			FreeList(GlobalIndexList);
			GlobalIndexList=NULL;
			GlobalIndexListTail=NULL;
			return;
		}
		counterInSeq++;
		newptr->count = counterInSeq;



		/* assign data to fields */
		sprintf(newptr->ref,"%d--%d,%d,%3.1f,%d,%d",
				AlignPair.indexprime[AlignPair.length],AlignPair.indexprime[1],best_match_distance,Copynumber,Classlength,(int ) OUTPUTcount);

		newptr->first = AlignPair.indexprime[AlignPair.length];
		newptr->last = AlignPair.indexprime[1];
		newptr->period = best_match_distance;
		newptr->copies = Copynumber;
		newptr->size = Classlength;
		newptr->matches = (int)(100*(float)match/x);
		newptr->indels = (int)(100*(float)indel/x);
		newptr->score = AlignPair.score;
		newptr->acount = (int)(100*(double)ACGTcount['A'-'A']/count);
		newptr->ccount = (int)(100*(double)ACGTcount['C'-'A']/count); 
		newptr->gcount = (int)(100*(double)ACGTcount['G'-'A']/count);
		newptr->tcount = (int)(100*(double)ACGTcount['T'-'A']/count);
		newptr->entropy = entropy;


		/* allocate memory to place the pattern and copy data into it*/
		newptr->pattern = (char*) malloc((consensussize+1)*sizeof(char));
		if(newptr->pattern==NULL)
		{
			free(newptr);
			FreeList(GlobalIndexList);
			GlobalIndexList=NULL;
			GlobalIndexListTail=NULL;
			return;
		}
		printECtoBuffer(newptr->pattern, startECpos, consensussize);

		if (GlobalIndexList == NULL) /* first element */
		{
			GlobalIndexList = GlobalIndexListTail = newptr;
			GlobalIndexListTail->next = NULL;
		}
		else /* add new element to end of list */
		{
			GlobalIndexListTail->next = newptr;
			GlobalIndexListTail = newptr;
			GlobalIndexListTail->next = NULL;
		}

	}



}
