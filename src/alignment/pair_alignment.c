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

/* Traceback of a wraparound alignment into AlignPair, with the copy number. */

#include "core/globals.h"
#include "core/report.h"
#include "alignment/similarity.h"
#include "alignment/pair_alignment.h"

#define fill_align_pair( c1, c2, l, i, j ) \
    AlignPair.textprime[l]  = c1;          \
    AlignPair.textsecnd[l]  = c2;          \
    AlignPair.indexprime[l] = i;           \
    AlignPair.indexsecnd[l] = j

/*******************************************************************/

/*********************  get_narrowband_pair_alignment_with_copynumber()  **********************/

/*******************************************************************/


void get_narrowband_pair_alignment_with_copynumber(int size,
		int bandradius,int option)


/* for a repeat of EC, do a traceback alignment, */
/* ending at row Maxrow and column Maxcol */




#define test_match_mismatch \
	if (S[r][i]==S[r-1][upi-1]+match(x[realr], y[c]))\
{\
	length++;\
	if(c==fullcopy)Copynumber++;\
	if(option==LOCAL){fill_align_pair(x[realr],y[c],length,realr,c);c=(c-1+size)%size;}\
	else {fill_align_pair(x[realr],y[c],length,realr,c+Maxrealcol-Maxcol);c=c-1;}\
	realr--;\
	r--;\
	i=upi-1;\
}\
else

#define test_up \
	if (S[r][i]==S[r-1][upi]+Delta)\
{\
	length++;\
	if(option==LOCAL) {fill_align_pair(x[realr],'-',length,realr,(c+1)%size);}\
	else {fill_align_pair(x[realr],'-',length,realr,c+1+Maxrealcol-Maxcol);}\
	realr--;\
	r--;\
	i=upi;\
}\
else

#define test_left \
	if(i==0)\
{\
	trf_message("\nget_pair_alignment_with_copynumber: error in trace back");\
	trf_message("\nattempted to compute left branch when i==0");\
	trf_message("\nS[%d][%d]=%d",r,i,S[r][i]);\
	break;\
}\
else if (S[r][i]==S[r][i-1]+Delta)\
{\
	length++;\
	if(c==fullcopy)Copynumber++;\
	if(option==LOCAL){fill_align_pair('-',y[c],length,realr+1,c);c=(c-1+size)%size;}\
	else {fill_align_pair('-',y[c],length,realr+1,c+Maxrealcol-Maxcol);c=c-1;}\
	i=i-1;\
}\
else

#define report_error_match_up_left \
{\
	trf_message("\nget_pair_alignment_with_copynumber: error in trace back");\
	trf_message("\nS: row=%d  column=%d  upi=%d  Sequence: realrow=%d  EC: realcol=%d",\
			r,i,upi,realr,c);\
	trf_message("\nS=%d  Sleft=%d  Sup=%d  Sdiag=%d  match=%d",\
			S[r][i],S[r][i-1],S[r-1][upi],\
			S[r-1][upi-1],match(x[realr], y[c]));\
	break;\
}
#define report_error_match_left \
{\
	trf_message("\nget_pair_alignment_with_copynumber: error in trace back");\
	trf_message("\nS: row=%d  column=%d  Sequence: realrow=%d  EC: realcol=%d",\
			r,i,realr,c);\
	trf_message("\nS=%d  Sleft=%d  Sdiag=%d  match=%d",\
			S[r][i],S[r][i-1],S[r-1][upi-1],match(x[realr], y[c]));\
	break;\
}
#define report_error_up_left \
{\
	trf_message("\nget_pair_alignment_with_copynumber: error in trace back");\
	trf_message("\nS: row=%d  column=%d  Sequence: realrow=%d  EC: realcol=%d",\
			r,i,realr,c);\
	trf_message("\nS=%d  Sleft=%d  Sup=%d",\
			S[r][i],S[r][i-1],S[r-1][upi]);\
	break;\
}
#define report_error_left \
{\
	trf_message("\nget_pair_alignment_with_copynumber: error in trace back");\
	trf_message("\nS: row=%d  column=%d  Sequence: realrow=%d  EC: realcol=%d",\
			r,i,realr,c);\
	trf_message("\nS=%d  Sleft=%d",\
			S[r][i],S[r][i-1]);\
	break;\
}


#define report_error_up \
{\
	trf_message("\nget_pair_alignment_with_copynumber: error in trace back");\
	trf_message("\nfailed to go up when c=-1");\
	trf_message("\nS: row=%d  column=%d  upi=%d  Sequence: realrow=%d  EC: realcol=%d",\
			r,i,upi,realr,c);\
	trf_message("\nS=%d  up=%d",\
			S[r][i],S[r-1][upi]);\
	break;\
}




{

	int i,length,fullcopy;
	unsigned char *x,*y;
	int realr,r,c,w;
	int upi,k;
	int legitimateZero;



	if(ldong) /* null statement */;

	else
	{
		w=bandradius;
		x=Sequence;
		y=EC;


		realr=Maxrealrow;
		r=Maxrow;
		c=Maxcol;
		k=(Maxcol-Bandcenter[r]+size)%size;
		if(size-k<=k) k=-(size-k);
		i=w+k;          
		fullcopy=(Maxcol+1)%size;

		AlignPair.score=S[r][i]; 
		length=0;
		Copynumber=0;
		for (;;)
		{      
			/* stop at zeros or -1000 for local */
			/* stop at r=0 for global */
			
			//change 8.20.12 Gary Benson
			//check for legitimate zero in alignment score (not from local alignment zeroing)
			//and allow to go back past that zero in alignment array
			//used to capture some alignments with mutations near the beginning
			//check to see if LOCAL zero score is a legitimate zero not a local zero
			legitimateZero=0;
			if((option==LOCAL)&&(S[r][i]==0)&&r!=0)
			{
				//get upi
				k=(c-Bandcenter[r-1]+size)%size;
				if(size-k<=k) k=-(size-k);
				upi=w+k;
				//get k
				k=(Bandcenter[r]-Bandcenter[r-1]+size)%size;
				if(size-k<=k) k=-(size-k);
				if(k>=1)      /* band shifts right */
				{
					if(i<=2*w-k)
					{
						if((S[r][i]==S[r-1][upi-1]+match(x[realr], y[c]))
						   ||(S[r][i]==S[r-1][upi]+Delta)
						   ||(S[r][i]==S[r][i-1]+Delta)) legitimateZero=1;
					}
					else if(i==2*w-k+1)
					{
						if((S[r][i]==S[r-1][upi-1]+match(x[realr], y[c]))
						   ||(S[r][i]==S[r][i-1]+Delta)) legitimateZero=1;
					}
					else 
					{
						if(S[r][i]==S[r][i-1]+Delta) legitimateZero=1;						
					}
				}
				else          /* (k<=0) band shifts left */
				{
					k=-k;
					if(i<=k-1)
					{
						if(S[r][i]==S[r][i-1]+Delta) legitimateZero=1;
					}
					else if(i==k)
					{
						if((S[r][i]==S[r-1][upi]+Delta)
						   ||(S[r][i]==S[r][i-1]+Delta)) legitimateZero=1;	
					}
					else
					{
						if((S[r][i]==S[r-1][upi-1]+match(x[realr], y[c]))
						   ||(S[r][i]==S[r-1][upi]+Delta)
						   ||(S[r][i]==S[r][i-1]+Delta)) legitimateZero=1;
					}
				}
				//debug
				//if (legitimateZero) {
					//report where
				//	printf("\nlegitimate zero at r:%d,i:%d",r,i);
				//}
			}
			
			
			/* stop at zeros or -1000 for local */
			/* stop at r=0 for global */
			if (r==0
				||((option==LOCAL)&&(S[r][i]<=0)&&!legitimateZero)
				||((option==GLOBAL)&&(r==0)&&(c==-1)))
				//debug (uncomment to restore)	
			/*
			 if (((option==LOCAL)&&(S[r][i]<=0))||((option==GLOBAL)&&(r==0)&&(c==-1)))
			 */
			{
				legitimateZero=0;
				AlignPair.length=length;
				if(Maxcol>=c){Copynumber+=((double)(Maxcol-c))/size;}
				else {Copynumber+=((double)(Maxcol+size-c))/size;}return;
			}
			else
			{
				k=(c-Bandcenter[r-1]+size)%size;
				if(size-k<=k) k=-(size-k);
				upi=w+k;
				/* in case global goes to c==-1 */
				if((option==GLOBAL)&&(c==-1))
				{
					test_up
						report_error_up
				}
				else
				{
					k=(Bandcenter[r]-Bandcenter[r-1]+size)%size;
					if(size-k<=k) k=-(size-k);
					if(k>=1)      /* band shifts right */
					{
						if(i<=2*w-k)
						{
							test_match_mismatch
								test_up
								test_left
								report_error_match_up_left;
						}
						else if(i==2*w-k+1)
						{
							test_match_mismatch
								test_left
								report_error_match_left;
						}
						else 
						{
							test_left
								report_error_left;
						}
					}
					else          /* (k<=0) band shifts left */
					{
						k=-k;
						if(i<=k-1)
						{
							test_left
								report_error_left;
						}
						else if(i==k)
						{
							test_up
								test_left
								report_error_up_left;
						}
						else
						{
							test_match_mismatch
								test_up
								test_left
								report_error_match_up_left;
						}
					}
				}
			}
		}
	}
}

/************************************************************/ 



/*******************************************************************/

/*********************  get_pair_alignment_with_copynumber()  **********************/

/*******************************************************************/


void get_pair_alignment_with_copynumber(int size)


	/* for a repeat of EC, do a traceback alignment, */
	/* ending at row Maxrow and column Maxcol */






{

	int i,j,si,length,adjlength,fullcopy;
	unsigned char *x,*y;



	x=Sequence;
	y=EC;


	i=Maxrealrow;
	si=Maxrow;
	j=Maxcol;
	adjlength=size-1;
	fullcopy=(j+1)%size;

	AlignPair.score=S[si][j]; 

	length=0;
	Copynumber=0;
	for (;;)
	{      
		/* stop at zeros or -1000 */
		if (S[si][j]<=0)
		{
			AlignPair.length=length;
			if(Maxcol>=j){Copynumber+=((double)(Maxcol-j))/size;}
			else {Copynumber+=((double)(Maxcol+size-j))/size;}return;
		}

		/* check match/mismatch branch */

		else if (S[si][j]==S[si-1][(j+adjlength)%size]+match(x[i], y[j]))
		{
			length++;
			if(j==fullcopy)Copynumber++;
			fill_align_pair(x[i],y[j],length,i,j);
			i--;
			si--;
			j=(j+adjlength)%size;
		}

		/* check deletion branch */
		else if (S[si][j]==S[si-1][j]+Delta)
		{
			length++;
			fill_align_pair(x[i],'-',length,i,(j+1)%size);
			i--;
			si--;
		}

		/* check other deletion branch */
		else if (S[si][j]==S[si][(j+adjlength)%size]+Delta)
		{
			length++;
			if(j==fullcopy)Copynumber++;
			fill_align_pair('-',y[j],length,i+1,j);
			j=(j+adjlength)%size;
		}
		else
		{
			trf_message("\nget_pair_alignment_with_copynumber: error in trace back");
			trf_message("\nrow=%d  column=%d",i,j);
			trf_message("\nS=%d  Sleft=%d  Sup=%d  Sdiag=%d  match=%d",S[si][j],
					S[si][(j+adjlength)%size],S[si-1][j],
					S[si-1][(j+adjlength)%size],match(x[i], y[j]));

			break;
		}
	}


}










/*******************************************************************/

/**************************  reverse()  ****************************/

/*******************************************************************/


void reverse()

	/* reverses the alignment in AlignPair */




{
	int j, tempi,ml;
	char temp;

	ml=AlignPair.length;

	for (j=1;j<=ml/2;j++)
	{

		temp=AlignPair.textprime[j];
		AlignPair.textprime[j]=AlignPair.textprime[ml-j+1];
		AlignPair.textprime[ml-j+1]=temp;

		temp=AlignPair.textsecnd[j];
		AlignPair.textsecnd[j]=AlignPair.textsecnd[ml-j+1];
		AlignPair.textsecnd[ml-j+1]=temp;

		tempi=AlignPair.indexprime[j];
		AlignPair.indexprime[j]=AlignPair.indexprime[ml-j+1];
		AlignPair.indexprime[ml-j+1]=tempi;

		tempi=AlignPair.indexsecnd[j];
		AlignPair.indexsecnd[j]=AlignPair.indexsecnd[ml-j+1];
		AlignPair.indexsecnd[ml-j+1]=tempi;
	}
}









void shift_pattern_indices (int patternsize)
{
	int downshift, upshift, l;

	downshift = AlignPair.indexsecnd[1];
	upshift = patternsize-downshift;
	for(l=1;l<=AlignPair.length;l++)
	{
		if(AlignPair.indexsecnd[l]>=downshift)
			AlignPair.indexsecnd[l]-=downshift;
		else
			AlignPair.indexsecnd[l]+=upshift;
	}
} 
