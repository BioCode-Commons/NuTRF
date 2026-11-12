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

/* Repeat table as plain text: the .dat file and the -ngs stream. */

#include "core/globals.h"
#include "core/report.h"
#include "core/repeat_list.h"
#include "output/data_file.h"

/* program heading of a .dat file; the -ngs stream has none */
void    WriteDataHeading(FILE* fp)
{
	if (paramset.ngs != 1) {
		fprintf(fp,"Tandem Repeats Finder Program written by:\n\n");
		fprintf(fp,"Gary Benson\n");
		fprintf(fp,"Program in Bioinformatics\n");
		fprintf(fp,"Boston University\n");
		fprintf(fp,"Version %s\n", versionstring);
	}
}

/* heading of the current sequence, then one line for each repeat in
   GlobalIndexList; -ngs lines also carry 50 bases of flanking sequence */
void    WriteDataRecords(FILE* fp, char* name)
{
	IL* lpointer;
	int charcount;

	if (paramset.ngs) {

		/* only print if we have at least 1 record */
		if (NULL!=GlobalIndexList) {
			fprintf(fp,"@%s\n",name);
		}

	} else {
		fprintf(fp,"\n\nSequence: %s\n\n\n\nParameters: %d %d %d %d %d %d %d\n\n\n",
				name,paramset.match,paramset.mismatch,paramset.indel,paramset.PM,paramset.PI,paramset.minscore,paramset.maxperiod);
	}

	for(lpointer=GlobalIndexList;lpointer!=NULL;lpointer=lpointer->next)
	{
		fprintf(fp,"%d %d %d %.1f %d %d %d %d %d %d %d %d %.2f %s ",
				lpointer->first, lpointer->last, lpointer->period,
				lpointer->copies, lpointer->size, lpointer->matches,
				lpointer->indels, lpointer->score, lpointer->acount,
				lpointer->ccount, lpointer->gcount, lpointer->tcount,
				lpointer->entropy, lpointer->pattern );
		for(charcount=lpointer->first; charcount<=lpointer->last;charcount++)
			fprintf(fp,"%c", Sequence[charcount]);

		/* print short flanks to .dat file */
		if (paramset.ngs) {
			int flankstart,flankend;

	  					flankstart = lpointer->first - 50; 
			flankstart=max(1,flankstart);
			flankend = lpointer->last + 50; 
			flankend=min(Length,flankend);

			fprintf(fp," ");				
			if (lpointer->first == 1) {
				fprintf(fp,".");
			} else {
				for(charcount=flankstart; charcount<lpointer->first;charcount++)
					fprintf(fp,"%c", Sequence[charcount]);
			}

			fprintf(fp," ");	
			if (lpointer->last == Length) {
				fprintf(fp,".");
			} else {
				for(charcount=lpointer->last+1; charcount<=flankend;charcount++)
					fprintf(fp,"%c", Sequence[charcount]);
			}

		}
		
		fprintf(fp,"\n");
	}
}


/***********GetList()*******************************************************/

IL*     GetList(char * datafile)
{
	FILE    *fp;
	IL      *headptr, *newptr, *lastptr;
	int     counter, i;
	char patbuffer[MAXDISTANCECONSTANT+1000];

	headptr=newptr=lastptr=NULL;

	/* open file */
	fp = fopen(datafile, "r");
	if ( NULL==fp ) return NULL;

	/* get hsequence line on ninth line of data file*/
	for (counter =0; counter<9; counter++) fgets(hsequence, 255, fp);

	/* get hparameters line on twelvth line of data file*/
	for (counter =0; counter<3; counter++) fgets(hparameters, 255, fp);

	/* get hlength from another global variable (bad practice)*/
	sprintf(hlength, "Length:  %d", Length);

	/* loop to fill out list from buffer */
	counter = 1; /* keeps track of order they are found */
	while(1)
	{
		/* create new index list element */
		newptr = (IL*) malloc(sizeof(IL));
		if(newptr==NULL)
		{
			FreeList(headptr);
			return NULL;
		}
		newptr->count = counter++;

		/* get data from file */
		i = fscanf(fp, "%s %d %d %d %f %d %d %d %d %d %d %d %d %f %s",
				newptr->ref, &newptr->first, &newptr->last,
				&newptr->period, &newptr->copies, &newptr->size,
				&newptr->matches, &newptr->indels, &newptr->score,
				&newptr->acount, &newptr->ccount, &newptr->gcount,
				&newptr->tcount, &newptr->entropy,patbuffer);

		if(i==EOF)
		{
			free(newptr);
			break;
		}

		/* allocate memory to place the pattern and copy data into it*/
		newptr->pattern = (char*) malloc((strlen(patbuffer)+1)*sizeof(char));
		if(newptr->pattern==NULL)
		{
			free(newptr);
			FreeList(headptr);
			return NULL;
		}
		strcpy(newptr->pattern,patbuffer);

		if (headptr == NULL) /* first element */
		{
			headptr = lastptr = newptr;
			lastptr->next = NULL;
		}
		else /* add new element to end of list */
		{
			lastptr->next = newptr;
			lastptr = newptr;
			lastptr->next = NULL;
		}
	}

	fclose(fp);
	return  headptr;
}

void    MakeDataFile(IL * headptr,char * datafile,int data)
{
	FILE* fp;
	IL* lpointer;
	int charcount;

	/* if data = 1 then produce new datafile overwriting old one */
	/* otherwise delete old file */
	if(data)
	{
		fp = fopen(datafile,"w");
		if(fp==NULL)
		{
			PrintError("Unable to open output file for writing in MakeDataFile routine!");
			exit(-1);
		}

		if (paramset.ngs != 1) {
			fprintf(fp,"Tandem Repeats Finder Program writen by:\n\nGary Benson\nProgram in Bioinformatics\nBoston University\nVersion %s\n\n\n%s\n\n\n%s\n\n",versionstring, hsequence, hparameters);
		}
		for(lpointer=headptr;lpointer!=NULL;lpointer=lpointer->next)
		{
			fprintf(fp,"%d %d %d %.1f %d %d %d %d %d %d %d %d %.2f %s ",
					lpointer->first, lpointer->last, lpointer->period,
					lpointer->copies, lpointer->size, lpointer->matches,
					lpointer->indels, lpointer->score, lpointer->acount,
					lpointer->ccount, lpointer->gcount, lpointer->tcount,
					lpointer->entropy, lpointer->pattern );
			for(charcount=lpointer->first; charcount<=lpointer->last;charcount++)
				fprintf(fp,"%c", Sequence[charcount]);
			fprintf(fp,"\n");
		}
		fclose(fp);
	}
	else remove(datafile);
	return;
}
