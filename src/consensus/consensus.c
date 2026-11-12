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

/* Consensus pattern of the aligned copies in AlignPair. */

#include "core/globals.h"
#include "consensus/consensus.h"

#if WEIGHTCONSENSUS
/*******************************************************************/

/****************************  get_consensus()  ********************/

/*******************************************************************/

void get_consensus(int patternsize)

{
	int c,lastindex,j,i,max,letters,newinsert,insertA,insertC,insertG,insertT;
	char maxchar;

	/* initialize counts */

	for (c=0;c<=2*(MAXPATTERNSIZE);c++)
	{
		Consensus.A[c]=0;
		Consensus.C[c]=0;
		Consensus.G[c]=0;
		Consensus.T[c]=0;
		Consensus.dash[c]=0;
		Consensus.insert[c]=0;  /* times insert occurs */
		Consensus.letters[c]=0; /* number of letters in all inserts */
		Consensus.total[c]=0;   /* occurrences of a position */
		Consensus.pattern[c]=DASH;
	}

	/* start consensus */

	lastindex=-1;

	i=1;
	while (i<=AlignPair.length)
	{
		if (AlignPair.indexsecnd[i]!=lastindex)
		{
			switch(AlignPair.textprime[i]){
				case 'A':
					Consensus.A[2*AlignPair.indexsecnd[i]+1]++;
					break;
				case 'C':
					Consensus.C[2*AlignPair.indexsecnd[i]+1]++;
					break;
				case 'G':
					Consensus.G[2*AlignPair.indexsecnd[i]+1]++;
					break;
				case 'T':
					Consensus.T[2*AlignPair.indexsecnd[i]+1]++;
					break;
				case DASH:
					Consensus.dash[2*AlignPair.indexsecnd[i]+1]++;
					break;
			}
			Consensus.total[2*AlignPair.indexsecnd[i]+1]++;
			if (lastindex!=-1)
			{
				if (AlignPair.indexsecnd[i]==patternsize-1)
				{
					Consensus.total[0]++;
				}
				else
				{
					Consensus.total[2*AlignPair.indexsecnd[i]+2]++;
				}
			}
			lastindex=AlignPair.indexsecnd[i];
			i++;

		}
		else            /* AlignPair.indexsecnd[i]==lastindex */
		{
			Consensus.insert[2*AlignPair.indexsecnd[i]]++;
			insertA=0;
			insertC=0;
			insertG=0;
			insertT=0;
			while((AlignPair.indexsecnd[i]==lastindex)
					&&
					(i<=AlignPair.length))
			{

				switch(AlignPair.textprime[i]){
					case 'A':
						if (insertA==0)
						{
							Consensus.A[2*AlignPair.indexsecnd[i]]++;
							insertA=1;
							break;
						}
					case 'C':
						if (insertC==0)
						{
							Consensus.C[2*AlignPair.indexsecnd[i]]++;
							insertC=1;
							break;
						}
					case 'G':
						if (insertG==0)
						{
							Consensus.G[2*AlignPair.indexsecnd[i]]++;
							insertG=1;
							break;
						}
					case 'T':
						if (insertT==0)
						{
							Consensus.T[2*AlignPair.indexsecnd[i]]++;
							insertT=1;
							break;
						}
				}
				Consensus.letters[2*AlignPair.indexsecnd[i]]++;
				i++;
			}
		}

	} 

	/* get consensus for letters that exist */
	for(i=1;i<=2*(patternlength);i+=2)
	{
		letters=Consensus.total[i]-Consensus.dash[i];
		max=Consensus.A[i];maxchar='A';
		if(max<Consensus.C[i]){max=Consensus.C[i];maxchar='C';}
		if(max<Consensus.G[i]){max=Consensus.G[i];maxchar='G';}
		if(max<Consensus.T[i]){max=Consensus.T[i];maxchar='T';}

		/* weighted consensus */
		if (((max*Alpha)+((letters-max)*Beta)+(Consensus.dash[i]*Delta))
				>=(letters*Delta))
			Consensus.pattern[i]=maxchar;
		else Consensus.pattern[i]=DASH;

	}

	/* get consensus for inserted letters */
	for(i=0;i<=2*(patternlength);i+=2)
	{
		if (Consensus.total[i]!=0)
		{
			max=Consensus.A[i];maxchar='A';
			if(max<Consensus.C[i]){max=Consensus.C[i];maxchar='C';}
			if(max<Consensus.G[i]){max=Consensus.G[i];maxchar='G';}
			if(max<Consensus.T[i]){max=Consensus.T[i];maxchar='T';}
			if (((max*Alpha)+
						((Consensus.insert[i]-max)*Beta)+
						((Consensus.letters[i]-Consensus.insert[i])*Delta)+
						((Consensus.total[i]-Consensus.insert[i])*Delta))
					>
					(Consensus.letters[i]*Delta))
				Consensus.pattern[i]=maxchar;
			else Consensus.pattern[i]=DASH;

		}
	}

	/* compress consensus */

	j=0;
	for(i=0;i<=2*patternlength;i++)
		if (Consensus.pattern[i]!=DASH)
		{
			Consensus.pattern[j]=Consensus.pattern[i];
			j++;
		}
	ConsClasslength=j;

}





#else


/*******************************************************************/

/****************************  get_consensus()  ********************/

/*******************************************************************/

void get_consensus(int patternlength)

{
	int c,lastindex,j,i,max;
	char maxchar;

	/* initialize counts */

	/* Y. Hernandez Oct 15, 2018 */
	/* Changed from MAXPATTERNSIZE to MAXPATTERNSIZECONSTANT since that
	 * is what is used to init the arrays in Consensus. Using
	 * MAXPATTERNSIZE will result in an out-of-bounds error if > 2000.
	 */
	for (c=0;c<=2*(MAXPATTERNSIZECONSTANT);c++)
	{
		Consensus.A[c]=0;
		Consensus.C[c]=0;
		Consensus.G[c]=0;
		Consensus.T[c]=0;
		Consensus.dash[c]=0;
		Consensus.insert[c]=0;
		Consensus.total[c]=0;
		Consensus.pattern[c]=DASH;
	}

	/* start consensus */

	lastindex=-1;

	i=1;
	while (i<=AlignPair.length)
	{
		if (AlignPair.indexsecnd[i]!=lastindex)
		{
			switch(AlignPair.textprime[i]){
				case 'A':
					Consensus.A[2*AlignPair.indexsecnd[i]+1]++;
					break;
				case 'C':
					Consensus.C[2*AlignPair.indexsecnd[i]+1]++;
					break;
				case 'G':
					Consensus.G[2*AlignPair.indexsecnd[i]+1]++;
					break;
				case 'T':
					Consensus.T[2*AlignPair.indexsecnd[i]+1]++;
					break;
				case DASH:
					Consensus.dash[2*AlignPair.indexsecnd[i]+1]++;
					break;
			}
			if (lastindex!=-1)
			{
				if (AlignPair.indexsecnd[i]==patternlength-1)
				{
					Consensus.total[0]++;
				}
				else
				{
					Consensus.total[2*AlignPair.indexsecnd[i]+2]++;
				}
			}
			lastindex=AlignPair.indexsecnd[i];
			i++;

		}
		else            /* AlignPair.indexsecnd[i]==lastindex */
		{
			Consensus.insert[2*AlignPair.indexsecnd[i]]++;

			while((AlignPair.indexsecnd[i]==lastindex)
					&&
					(i<=AlignPair.length))
			{

				switch(AlignPair.textprime[i]){
					case 'A':
						Consensus.A[2*AlignPair.indexsecnd[i]]++;
						break;
					case 'C':
						Consensus.C[2*AlignPair.indexsecnd[i]]++;
						break;
					case 'G':
						Consensus.G[2*AlignPair.indexsecnd[i]]++;
						break;
					case 'T':
						Consensus.T[2*AlignPair.indexsecnd[i]]++;
						break;
					case DASH:
						Consensus.dash[2*AlignPair.indexsecnd[i]]++;
						break;
				}
				i++;
			}
		}

	} 

	for(i=1;i<=2*(patternlength);i+=2)
	{
		max=Consensus.dash[i];
		maxchar=DASH;
		if(max<Consensus.A[i]){max=Consensus.A[i];maxchar='A';}
		if(max<Consensus.C[i]){max=Consensus.C[i];maxchar='C';}
		if(max<Consensus.G[i]){max=Consensus.G[i];maxchar='G';}
		if(max<Consensus.T[i]){max=Consensus.T[i];maxchar='T';}
		Consensus.pattern[i]=maxchar;
	}
	for(i=0;i<=2*(patternlength);i+=2)
	{
		if ((Consensus.total[i]!=0)
				&&
				(((float)Consensus.insert[i]/Consensus.total[i])>=0.5))
		{
			max=Consensus.A[i];maxchar='A';
			if(max<Consensus.C[i]){max=Consensus.C[i];maxchar='C';}
			if(max<Consensus.G[i]){max=Consensus.G[i];maxchar='G';}
			if(max<Consensus.T[i]){max=Consensus.T[i];maxchar='T';}
			/*      if (max<Consensus.total[i]/2){maxchar='X';} */
		}


		else maxchar=DASH;
		Consensus.pattern[i]=maxchar;

	}
	j=0;
	for(i=0;i<=2*patternlength;i++)
		if (Consensus.pattern[i]!=DASH)
		{
			Consensus.pattern[j]=Consensus.pattern[i];
			j++;
		}
	ConsClasslength=j;

}


#endif
