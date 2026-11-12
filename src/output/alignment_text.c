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

/* Text of one alignment in the alignment file. */

#include "core/globals.h"
#include "alignment/pair_alignment.h"
#include "output/alignment_text.h"

/*******************************************************************/

/***********************  alt_print_alignment()  ***********************/

/*******************************************************************/

void alt3_print_alignment(int patternwidth) /* change added to make indices run
					       from 1 to patternsize in the 
					       output */

	/* prints out the alignment in AlignPair */

{

	extern int pwidth;
	int i, j, g, h, first, m;

	shift_pattern_indices(patternwidth);

	if (pwidth>0)
	{
		if(AlignPair.indexprime[1]!=1)
		{
			m=AlignPair.indexprime[1]-10;
			if (m<1) m=1;
			j=m;
			fprintf(Fptxt,"  %9d ",j);
			for(i=1;i<=10;i++)
			{
				fputc(Sequence[j],Fptxt);
				j++;
				if(j==AlignPair.indexprime[1]) break;
			}
			fprintf(Fptxt,"\n\n");
		}
		g=1;
		for (;;)
		{
			j=g;
			fprintf(Fptxt,"            ");
			first=TRUE;
			h=0;
			for (i=0;i<pwidth-10; ++i)
			{
				if (j==AlignPair.length+1)
					break;
				if ((!first)&&(AlignPair.indexsecnd[j]==0)
						&&(AlignPair.indexsecnd[j-1]!=0))
				{
					if ((patternwidth>6)||(pwidth-10-h<(2*patternwidth)))
						break;
					fprintf(Fptxt," ");
					h++;
				}
				first=FALSE;
				if ((AlignPair.textprime[j]!=AlignPair.textsecnd[j])&&
						(AlignPair.textprime[j]!='-') &&
						(AlignPair.textsecnd[j]!='-'))
				{fputc('*',Fptxt);}
				else fputc(' ',Fptxt);
				j++;
				h++;
			}  
			fputc('\n',Fptxt);
			j=g;
			fprintf(Fptxt,"  %9d ",AlignPair.indexprime[j]);
			first=TRUE;
			h=0;
			for (i=0;i<pwidth-10; ++i)
			{
				if (j==AlignPair.length+1) break;
				if ((!first)&&(AlignPair.indexsecnd[j]==0)
						&&(AlignPair.indexsecnd[j-1]!=0))
				{
					if ((patternwidth>6)||(pwidth-10-h<(2*patternwidth)))
						break;
					fprintf(Fptxt," ");
					h++;
				}
				first=FALSE;
				fputc(AlignPair.textprime[j++],Fptxt);
				h++;
			}  
			fputc('\n',Fptxt);
			j=g;
			fprintf(Fptxt,"  %9d ",AlignPair.indexsecnd[j]+1); /* +1 added to make 
									      indices run from 
									      1 to 
									      patternwidth in 
									      output */
			first=TRUE;
			h=0;
			for (i=0;i<pwidth-10; ++i)
			{
				if (j==AlignPair.length+1)
					break;
				if ((!first)&&(AlignPair.indexsecnd[j]==0)
						&&(AlignPair.indexsecnd[j-1]!=0))
				{
					if ((patternwidth>6)||(pwidth-10-h<(2*patternwidth)))
						break;
					fprintf(Fptxt," ");
					h++;
				}
				first=FALSE;
				fputc(AlignPair.textsecnd[j++],Fptxt);
				h++;
			}
			fprintf(Fptxt,"\n\n");
			g=j;
			if (j==AlignPair.length+1) break;
		}
		if(AlignPair.indexprime[AlignPair.length]!=Length)
		{
			m=AlignPair.indexprime[AlignPair.length]+10;
			if (m>Length) m=Length;
			j=AlignPair.indexprime[AlignPair.length]+1;
			fprintf(Fptxt,"  %9d ",j);
			for(i=1;i<=10;i++) 
			{
				fputc(Sequence[j],Fptxt); 
				j++;
				if(j>m) break;
			}
			fprintf(Fptxt,"\n\n");
		}
	}
	else
		fprintf(Fptxt,"Error, pwidth<=0, can't report alignments.\n");
}



void print_alignment_headings(int consensuslength)
{
	/* headings */
	if (Heading==0)
	{
		Heading=1;
	}
	fprintf(Fptxt,"\n\n<A NAME=\"%d--%d,%d,%3.1f,%d,%d\">",
			AlignPair.indexprime[1],
			AlignPair.indexprime[AlignPair.length],Period,Copynumber,consensuslength,(int) OUTPUTcount);
	fprintf(Fptxt,"</A>");


#if defined(WINDOWSGUI)
	fprintf(Fptxt,"<P>See <FONT COLOR=\"#0000FF\">Alignment Explanation</FONT> in Tandem Repeats Finder Help</P><BR>\n");
#elif defined(WINDOWSCONSOLE)
	fprintf(Fptxt,"<A HREF=\"http://tandem.bu.edu/trf/trf.definitions.html#alignment\" target =\"explanation\">Alignment explanation</A><BR><BR>\n");
#elif defined(UNIXGUI)
	fprintf(Fptxt,"<P>See <FONT COLOR=\"#0000FF\">Alignment Explanation</FONT> in Tandem Repeats Finder Help</P><BR>\n");
#elif defined(UNIXCONSOLE)
	fprintf(Fptxt,"<A HREF=\"http://tandem.bu.edu/trf/trf.definitions.html#alignment\" target =\"explanation\">Alignment explanation</A><BR><BR>\n");
#endif




	fprintf(Fptxt,"    Indices: %d--%d",
			AlignPair.indexprime[1],
			AlignPair.indexprime[AlignPair.length]);
	fprintf(Fptxt,"  Score: %d",Maxscore);
	fprintf(Fptxt,"\n    Period size: %d  Copynumber: %3.1f  Consensus size: %d\n\n",
			Period,Copynumber,consensuslength);
}  

/******************** flanking sequence *********************/

void print_flanking_sequence(int flank_length)
{

	int m,n,k,i,j;

	m=AlignPair.indexprime[1]-flank_length;
	if (m<1) m=1;
	n=AlignPair.indexprime[AlignPair.length]+flank_length;
	if (n>Length) n=Length;

	if (m==AlignPair.indexprime[1])
	{
		fprintf(Fptxt,"\nLeft flanking sequence: None");
	}
	else
	{
		fprintf(Fptxt,"\nLeft flanking sequence: Indices %d -- %d\n",
				m,AlignPair.indexprime[1]-1);
		k=AlignPair.indexprime[1];
		j=m;
		for(;;)
		{
			for(i=1;i<=pwidth-10;i++)
			{
				fputc(Sequence[j],Fptxt);
				j++;
				if(j>=k) break;
			}
			fprintf(Fptxt,"\n");
			if(j>=k) break;
		} 
	}
	if (n==AlignPair.indexprime[AlignPair.length])
	{
		fprintf(Fptxt,"\n\nRight flanking sequence: None");
	}
	else
	{
		fprintf(Fptxt,"\n\nRight flanking sequence: Indices %d -- %d\n",
				AlignPair.indexprime[AlignPair.length]+1,n);
		j=AlignPair.indexprime[AlignPair.length]+1;
		for(;;)
		{
			for(i=1;i<=pwidth-10;i++)
			{
				fputc(Sequence[j],Fptxt);
				j++;
				if(j>n) break;
			}
			fprintf(Fptxt,"\n");
			if(j>n) break;
		} 
	}
	fprintf(Fptxt,"\n\n");

	return;

}


void printECtoAlignments(FILE* fp, int start, int width)
{
	int counter;
	int chars=0;

	fprintf(Fptxt,"\nConsensus pattern (%d bp):   ",width);
	for(counter=start;counter<width;counter++,chars++)
	{
		if(chars%65==0) fprintf(fp,"\n");
		fprintf(fp,"%c",EC[counter]);
	}
	for(counter=0;counter<start;counter++,chars++)
	{
		if(chars%65==0) fprintf(fp,"\n");
		fprintf(fp,"%c",EC[counter]);
	}

	fprintf(fp,"\n");

	return;
}

void printECtoBuffer(char *trg, int start, int width)
{
	int counter;

	for(counter=start;counter<width;counter++,trg++)
	{
		*trg=EC[counter];
	}
	for(counter=0;counter<start;counter++,trg++)
	{
		*trg=EC[counter];
	}
	*trg='\0';
	return;
}

void printECtoData(FILE* fp, int start, int width)
{
	int counter;

	for(counter=start;counter<width;counter++)
	{
		fprintf(fp,"%c",EC[counter]);
	}
	for(counter=0;counter<start;counter++)
	{
		fprintf(fp,"%c",EC[counter]);
	}
	fprintf(fp,"\n");
	return;
}
