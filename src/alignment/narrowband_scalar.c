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

/* Scalar narrowband kernel. */

#include "core/globals.h"
#include "core/report.h"
#include "alignment/narrowband_scalar.h"

/*******************************************************************/

/***************************** narrowbandwrap()  **************************/

/*******************************************************************/
void trf_narrowband_scalar(int start, int size,int bandradius, int bandradiusforward, 
		int option,int tuplesize)

/* start is end of pattern in text */
/* tuplesize is the size of tuple used for this pattern size */

/* started Feb. 7 1997 */
{
	int g;
	register int *pup, *pdiag, *pcurr,pleft;
	int c,realr,end_of_trace,
	    maxscore,maxrow,maxcol,mincol;
	int maxrealrow, minrealrow;
	char currchar;
	const int *similarity_row;
	/* Diag/Up originally copy the preceding S row and S+Delta. Each
	 * diagonal index is consumed once, so the read-only row plus local
	 * diag_score preserves the immediate match observer without stores. */
	const int *previous, *previous_diag, *previous_up;
	int diag_score;
	int w,matches_in_diagonal,matchatmax_col,i,k,maxrowscore,
	    lastmatchatmax_col,match_yes_no;
	/*** 6/9/05 G. Benson ***/ 
	int mincolbandcenter,zeroat,mincolposition;
	/* Feb 16, 2016 Yozen */
	unsigned int r;

	/* G. Benson */
	/* 1/26/10 */
	/* change MAXWRAPLENGTH to MAXWRAPLENGTHCONST so MAXWRAPLENGTH can be used as an int */
	/* int Bandcenter[MAXWRAPLENGTH+1]; */
	// if (Bandcenter == NULL) {
	// 	Bandcenter = calloc(maxwraplength+1, sizeof(*Bandcenter));
	// }
	/* Bandcenter init removed by Yozen on Feb 16, 2016. */

	w=bandradius;
	if(MAXBANDWIDTH<2*w+1)
	{
		trf_message("\nIn narrowbandwrap, MAXBANDWIDTH: %d exceeded by 2*w+1: %d\n",
				MAXBANDWIDTH,2*w+1);
		exit(-1);
	}

	/* fill EC */
	if(option==WITHCONSENSUS)
		for(g=0;g<size;g++) EC[g]=Consensus.pattern[g];
	else
		for(g=0;g<size;g++) EC[g]=Sequence[start-size+g+1];

	/* backward wdp */
	maxscore=0;
	realr=start+1;
	r=maxwraplength;
	Bandcenter[r]=0;
	matches_in_diagonal=0;
	matchatmax_col=-2;

	pcurr=&S[r][0];
	pdiag=&Diag[0];
	pup=&Up[0];




	/* 3/14/05 gary benson -- reverse direction */
	/* change zeroth row values to put in -1000 in unreachable cells
	   and gap penalty in cells beyond start location */


	for (i=0;i<=w;i++)
	{
		*pup=(*pdiag=(*pcurr=0+Delta*(w-i)))+ Delta;
		pup++;
		pdiag++;
		pcurr++;
	}

	for (i=w+1;i<=2*w;i++)
	{
		*pup=(*pdiag=(*pcurr=-1000))+ Delta;
		pup++;
		pdiag++;
		pcurr++;
	}





	/* backwards */

	end_of_trace=FALSE;  
	while ((!end_of_trace) && (realr>1) && (r>0))
	{
		r--;
		realr--;
		Rows++;
		maxrowscore=-1;
		lastmatchatmax_col=matchatmax_col;
		end_of_trace=TRUE;
		currchar=Sequence[realr];
		similarity_row=SM+256*currchar;
		previous=S[r+1];
		pcurr=&S[r][2*w];
		pleft=-1000;

		if(matches_in_diagonal>=tuplesize)
		{
			/* recenter band */
			Bandcenter[r]=trf_narrowband_remainder(matchatmax_col-1+size,size);
		}
		else
		{

			/* don't recenter */
			Bandcenter[r]=trf_narrowband_remainder(Bandcenter[r+1]-1+size,size);
		}

		/* change of bandcenter determines which inputs go into which
		   cells */

		k=trf_narrowband_remainder(Bandcenter[r]-Bandcenter[r+1]+size,size);
		if(size-k<=k) k=-(size-k);
		c=trf_narrowband_remainder(Bandcenter[r]+w,size);
		if(k<=-1)           /* band shifts left */
		{
			k=-k;
			previous_diag=&previous[2*w-k+1];
			previous_up=&previous[2*w-k];
			for(i=2*w;i>=k;i--)
			{
				diag_score=*previous_diag+(match_yes_no=similarity_row[EC[c]]);
				pleft=(*pcurr=max4(0,diag_score,*previous_up+Delta,pleft))+Delta;
				test_trace_and_backwards_maxscore;
				test_maxrowscore_with_match;
				pcurr--;
				previous_diag--;
				previous_up--;
				c=(c==0)?size-1:c-1;
			}
			i=k-1;
			diag_score=*previous_diag+(match_yes_no=similarity_row[EC[c]]);
			pleft=(*pcurr=max3(0,diag_score,pleft))+Delta;
			test_trace_and_backwards_maxscore;
			test_maxrowscore_with_match;
			pcurr--;
			c=(c==0)?size-1:c-1;
			for(i=k-2;i>=0;i--)
			{
				pleft=(*pcurr=max2(0,pleft))+Delta;
				test_trace_and_backwards_maxscore;
				test_maxrowscore_without_match;
				pcurr--;
				c=(c==0)?size-1:c-1;
			}
		}
		else            /* band shifts right */
		{
			previous_diag=&previous[2*w];
			previous_up=&previous[2*w];
			for(i=2*w;i>=2*w-k+1;i--)
			{
				pleft=(*pcurr=max2(0,pleft))+Delta;
				test_trace_and_backwards_maxscore;
				test_maxrowscore_without_match;
				pcurr--;
				c=(c==0)?size-1:c-1;
			}
			i=2*w-k;
			pleft=(*pcurr=max3(0,*previous_up+Delta,pleft))+Delta;
			test_trace_and_backwards_maxscore;
			test_maxrowscore_without_match;
			pcurr--;
			previous_up--;
			c=(c==0)?size-1:c-1;
			for(i=2*w-k-1;i>=0;i--)
			{
				diag_score=*previous_diag+(match_yes_no=similarity_row[EC[c]]);
				pleft=(*pcurr=max4(0,diag_score,*previous_up+Delta,pleft))+Delta;
				test_trace_and_backwards_maxscore;
				test_maxrowscore_with_match;
				pcurr--;
				previous_diag--;
				previous_up--;
				c=(c==0)?size-1:c-1;
			}
		}
		if(trf_narrowband_remainder(matchatmax_col-lastmatchatmax_col+size,size)==size-1)
			matches_in_diagonal++;
		else matches_in_diagonal=0;
	}
	/* Restore backward scratch once; forward may initialize a narrower band. */
	pcurr=&S[r][0];
	pdiag=&Diag[0];
	pup=&Up[0];
	for (i=0;i<=2*w;i++)
	{
		*pup=(*pdiag=*pcurr)+Delta;
		pcurr++;
		pdiag++;
		pup++;
	}


	/**********************************************************************/
	/**********************************************************************/
	/**********************************************************************/
	/**********************************************************************/

	/* forward */

	if (ldong==1) /* null statement */;
	else
	{

		r=0;
		realr=minrealrow-1;
		/* mincol matches minrealrow in previous alignment above */

		/*** 6/14/05 G. Benson ***/
		/* modification version g */
		Bandcenter[0]=trf_narrowband_remainder(mincolbandcenter-1+size,size);
		matches_in_diagonal=0;
		matchatmax_col=-2;
		maxscore=0;
		pup=&Up[0];
		pdiag=&Diag[0];
		pcurr=&S[0][0];




		/* 3/14/05 gary benson -- forward direction */
		/* change zeroth row values to put in -1000 in unreachable cells
		   and gap penalty in cells beyond start location */

		/* initialize top row */

		/*** 6/14/05 G. Benson ***/
		/* set zero at mincol, not at bandcenter, which is now mincolbandcenter, 
		   so that forward alignment starts out the same as backwards alignment ended. */ 
		zeroat=mincolposition-w;



		/* July 9, 20009, G. Benson */
		/* widening the narrowband in the forward direction to catch alignments that don't return to correct start column */
		w = bandradiusforward;

		if(MAXBANDWIDTH<2*w+1)
		{
			trf_message("\nIn narrowbandwrap, MAXBANDWIDTH: %d exceeded by bandradiusforward 2*w+1: %d\n",
					MAXBANDWIDTH,2*w+1);
			exit(-1);
		}

		//	if((zeroat>w)||(zeroat<-w)) fprintf(stderr,"\nfound zeroat out of bounds");
		/*** 6/14/05 G. Benson restore what was in all earlier versions ***/
		for(i=0;i<w+zeroat;i++)
		{
			*pup=(*pdiag=(*pcurr=-1000))+Delta;
			pup++;
			pdiag++;
			pcurr++;
		}

		/*** 6/9/05 G. Benson ***/
		for(i=w+zeroat;i<=2*w;i++)
		{
			/*** 6/9/05 G. Benson ***/
			*pup=(*pdiag=(*pcurr=0+Delta*(i-(w+zeroat))))+Delta;
			pup++;
			pdiag++;
			pcurr++;
		}




		/* compute until end of trace */      
		end_of_trace=FALSE;  
		while ((!end_of_trace) && (realr<Length) && (r<maxwraplength))
		{
			r++;
			realr++;
			Rows++;
			end_of_trace=TRUE;
			maxrowscore=-1;
			lastmatchatmax_col=matchatmax_col;
			currchar=Sequence[realr];
			similarity_row=SM+256*currchar;
			previous=S[r-1];
			pcurr=&S[r][0];
			pleft=-1000;      /* don't use pleft for first entry */   

			if(matches_in_diagonal>=tuplesize)
			{
				/* recenter band */
				Bandcenter[r]=trf_narrowband_remainder(matchatmax_col+1,size);
			}
			else
			{

				/* don't recenter */
				Bandcenter[r]=trf_narrowband_remainder(Bandcenter[r-1]+1,size);
			}

			/* change of bandcenter determines which inputs go into which
			   cells */

			k=trf_narrowband_remainder(Bandcenter[r]-Bandcenter[r-1]+size,size);
			if(size-k<=k) k=-(size-k);
			if(k>=1)          /* band shifts right */
			{
				previous_diag=&previous[k-1];
				previous_up=&previous[k];
				c=trf_narrowband_remainder(Bandcenter[r]-w+size,size);
				for(i=0;i<=2*w-k;i++)
				{
					diag_score=*previous_diag+(match_yes_no=similarity_row[EC[c]]);
					pleft=(*pcurr=max4(0,diag_score,*previous_up+Delta,pleft))+Delta;
					test_trace_and_forward_maxscore;
					test_maxrowscore_with_match;
					pcurr++;
					previous_diag++;
					previous_up++;
					c=(c+1==size)?0:c+1;
				}
				i=2*w-k+1;
				diag_score=*previous_diag+(match_yes_no=similarity_row[EC[c]]);
				pleft=(*pcurr=max3(0,diag_score,pleft))+Delta;
				test_trace_and_forward_maxscore;
				test_maxrowscore_with_match;
				pcurr++;
				c=(c+1==size)?0:c+1; 
				for(i=2*w-k+2;i<=2*w;i++)
				{
					pleft=(*pcurr=max2(0,pleft))+Delta;
					test_trace_and_forward_maxscore;
					test_maxrowscore_without_match;
					pcurr++;
					c=(c+1==size)?0:c+1;
				}
			}
			else          /* band shifts left */
			{
				k=-k;
				c=trf_narrowband_remainder(Bandcenter[r]-w+size,size);
				previous_up=&previous[0];
				previous_diag=&previous[0];
				for(i=0;i<=k-1;i++)
				{
					pleft=(*pcurr=max2(0,pleft))+Delta;
					test_trace_and_forward_maxscore;
					test_maxrowscore_without_match;
					pcurr++;
					c=(c+1==size)?0:c+1;
				}
				i=k;
				pleft=(*pcurr=max3(0,*previous_up+Delta,pleft))+Delta;
				test_trace_and_forward_maxscore;
				test_maxrowscore_without_match;
				pcurr++;
				previous_up++;
				c=(c+1==size)?0:c+1;
				for(i=k+1;i<=2*w;i++)
				{
					diag_score=*previous_diag+(match_yes_no=similarity_row[EC[c]]);
					pleft=(*pcurr=max4(0,diag_score,*previous_up+Delta,pleft))+Delta;
					test_trace_and_forward_maxscore;
					test_maxrowscore_with_match;
					pcurr++;
					previous_diag++;
					previous_up++;
					c=(c+1==size)?0:c+1;
				}
			}

			if(trf_narrowband_remainder(matchatmax_col-lastmatchatmax_col+size,size)==1)
				matches_in_diagonal++;
			else matches_in_diagonal=0;
		}
		/* Restore the caller-visible scratch from the last forward row. */
		pcurr=&S[r][0];
		pdiag=&Diag[0];
		pup=&Up[0];
		for (i=0;i<=2*w;i++)
		{
			*pup=(*pdiag=*pcurr)+Delta;
			pcurr++;
			pdiag++;
			pup++;
		}
		/* store last position scanned with this pattern */
		Wrapend=realr;

		if (Criteria_print) trf_message("   DWPlength:%d",Wrapend-(minrealrow-1)+1);
		/* test for report */
		if (maxscore>=Reportmin)
		{
			Maxrealrow=maxrealrow;
			Maxrow=maxrow;
			Maxcol=maxcol;
			Maxscore=maxscore;

		}

	}

}
