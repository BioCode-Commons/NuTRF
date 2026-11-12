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

/* Best periods of recent alignments, used to filter large candidate distances. */

#include "core/globals.h"
#include "core/report.h"
#include "periods/best_period_list.h"

/**********************************************************/
/**********************************************************/

void init_bestperiodlist(void)
{
	Bestperiodlist->next=NULL;
}

/**********************************************************/
/**********************************************************/

void free_bestperiodlist(void)
{
	struct bestperiodlistelement *entry, *entrylast;
	entry=Bestperiodlist->next;
	Bestperiodlist->next=NULL;
	while (entry!=NULL)
	{
		entrylast=entry;
		entry=entry->next;
		free(entrylast);
	}
}

/**********************************************************/
/**********************************************************/

/*** modified 6/2/05 G. Benson ***/
void add_to_bestperiodlist(int d)

{
	struct bestperiodlistelement *ptr;

	/*** added 6/2/05 G. Benson***/
	/* when distance is 1, no Sortmultiples are defined due to change in multiples_criteria_4*/
	if(d!=1)
	{

		ptr=(struct bestperiodlistelement *)
			calloc(1,sizeof(struct bestperiodlistelement));
		if(ptr==NULL){trf_message("\nAdd_to_bestperiodlist: Out of memory!");exit(-1);}

		ptr->indexlow=AlignPair.indexprime[AlignPair.length]; 
		ptr->indexhigh=AlignPair.indexprime[1];
		/*** G. Benson 6/2/05 fixed bug using Sortmultiples 1-5 instead of 0-4  ***/
		ptr->best1=Sortmultiples[0];
		ptr->best2=Sortmultiples[1];
		ptr->best3=Sortmultiples[2];
		ptr->best4=Sortmultiples[3];
		ptr->best5=Sortmultiples[4];
		ptr->next=Bestperiodlist->next;
		Bestperiodlist->next=ptr;
	}
}


/*** modified 6/2/05 G. Benson ***/
void adjust_bestperiod_entry(int d)
	/* shortens length of best period entry if the consensus alignment 
	   turns out to be shorter than the first alignment; similar 
	   to distanceseen */
{
	struct bestperiodlistelement *ptr;

	/*** added 6/2/05 G. Benson***/
	/* when distance is 1, no Sortmultiples are defined due to change in multiples_criteria_4*/
	if(d!=1)
	{

		ptr=Bestperiodlist->next;
		if(ptr->indexhigh>AlignPair.indexprime[1])
			ptr->indexhigh=AlignPair.indexprime[1];

	}
}
/**********************************************************/
/**********************************************************/

/* difference plus/minus 5 percent of smaller */
/* +/- 1 percent*/
#define XYZ(a,b) ((a>b)?(a-b)<=(0.01*b):(b-a)<=(0.01*a))

int search_for_range_in_bestperiodlist(int start, int distance)

	/* returns true if
	   1) an entry brackets the range start-distance to start and
	   the distance matches one of the first three best periods
	   2) no entry brackets the range
	   returns false otherwise */

{
	struct bestperiodlistelement *entry, *entrylast, *temp;
	int range_covered;

	entry=Bestperiodlist->next;
	entrylast=Bestperiodlist;
	range_covered=FALSE;
	while (entry!=NULL)
	{
		if(entry->indexhigh<start-2*MAXDISTANCE)
		{
			/* remove current entry, too far back */

			entrylast->next=entry->next;
			temp=entry;
			entry=entrylast;
			free(temp);
		}
		else 
		{
			/* this specifies how much must be bracketed */
			if((entry->indexlow<=start-2*distance+1+Distance[distance].waiting_time_criteria)
					&&(entry->indexhigh>=start)/*&&(entry->indexhigh-entry->indexlow<=3*distance)*/)
			{
				range_covered=TRUE;

				if((entry->best1==distance)
						||(entry->best2==distance)
						||(entry->best3==distance)
						||(entry->best4==distance)
						||(entry->best5==distance)) return(TRUE);


				/* change 1 */
				/*

				   if(XYZ(entry->best1,distance)
				   ||XYZ(entry->best2,distance)
				   ||XYZ(entry->best3,distance)
				   ||XYZ(entry->best4,distance)
				   ||XYZ(entry->best5,distance)) return(TRUE);
				 */


			}
		}
		/* go to next entry */
		entrylast=entry;  
		entry=entry->next;
	}
	if(!range_covered) return(TRUE);
	else return(FALSE);


}
