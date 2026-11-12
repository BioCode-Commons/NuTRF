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

/* Extent of earlier alignments at each distance. */

#include "core/globals.h"
#include "core/report.h"
#include "detection/distance_seen.h"

/*******************************************************************/
/*******************************************************************/
void init_distanceseenarray(void)
	/* created 5/23/05 G. Benson */
{
	Distanceseenarray=(struct distanceseenarrayelement *)
		calloc(MAXDISTANCECONSTANT+1,sizeof(struct distanceseenarrayelement));
	if(Distanceseenarray==NULL){trf_message("\nInit Distanceseenarray: Out of memory!");
		exit(-1);}
}
/*******************************************************************/
/*******************************************************************/
void free_distanceseenarray(void)
	/* created 5/23/05 G. Benson */
{
	free(Distanceseenarray);
}
/*******************************************************************/
/*******************************************************************/
void add_to_distanceseenarray(int location,int distance,int end,int score)
	/* created 5/23/05 G. Benson */
	/* adds the extent of an alignment (end) at the given distance 
	   location and score are not currently used */
{
	struct distanceseenarrayelement *ptr;

	ptr = &Distanceseenarray[distance];

	ptr->index=location;
	ptr->end=end; 
	ptr->score=score;

}
/*******************************************************************/
/*******************************************************************/
int search_for_distance_match_in_distanceseenarray(int distance,int start)
	/* created 5/23/05 G. Benson */
	/* searches for an alignment with patternsize of distance in the region 
	   including start.  True means found and alignment should be blocked.*/
{
	struct distanceseenarrayelement ptr;

	ptr = Distanceseenarray[distance];

	if (ptr.end>=start) return(TRUE);
	else return(FALSE);
}
/*******************************************************************/
/*******************************************************************/


/*******************************************************************/
/*******************************************************************/
/*******************************************************************/
/*******************************************************************/
void init_distanceseenlist()
{
	Distanceseenlist->next=NULL;
}

/*******************************************************************/
/*******************************************************************/

void free_distanceseenlist()
{
	struct distancelistelement *entry, *entrylast;
	entry=Distanceseenlist->next;
	Distanceseenlist->next=NULL;
	while (entry!=NULL)
	{
		entrylast=entry;
		entry=entry->next;
		free(entrylast);
	}
}

/*******************************************************************/
/*******************************************************************/

void add_to_distanceseenlist(int location,int distance,int end,int  score, int acceptstatus)



{
	int changed_from_distance; /* 3/14/05 */
	struct distancelistelement *ptr, *entry;
	ptr=(struct distancelistelement *)
		calloc(1,sizeof(struct distancelistelement));
	if(ptr==NULL){trf_message("\nAdd_to_distanceseenlist: Out of memory!");exit(-1);}

	/* 3/15/02 Gary Benson
	   fix to remove from Distanceseenlist the first report of a repeat when
	   it blocks further along in the sequence than reported by the consensus 
	   sequence */
	entry=Distanceseenlist->next;
	/*
	   if((acceptstatus==WITHCONSENSUS)&&(distance==entry->distance)
	   &&(location==entry->index)&&(end<=entry->end))
	 */
	if((acceptstatus==WITHCONSENSUS) /* removed matching d because consensussize may be different */
			&&(location==entry->index)&&(end<=entry->end)&&
			((distance<=250)||(distance==entry->distance))) /* put back d match unless d under 250 */
	{
		/* don't remove Distanceseenlist element,
		   but make end=0 so it will be removed by 
		   search_for_distance_in_distanceseenlist
		   procedure */
		entry->end=0;

	}

	/* 3/14/05 */
	if(acceptstatus==WITHCONSENSUS) changed_from_distance=entry->distance;

	ptr->index=location;
	ptr->distance=distance;
	/* 3/14/05 gary benson -- new field to allow blocking of a distance that always shifts
	   to a nearby distance when consensus is computed */
	if(acceptstatus==WITHCONSENSUS) ptr->changed_from_distance=changed_from_distance;
	else ptr->changed_from_distance=distance;
	ptr->end=end;
	ptr->score=score;
	ptr->accepted=acceptstatus;
	ptr->next=Distanceseenlist->next;
	Distanceseenlist->next=ptr;
}

/*******************************************************************/
/*******************************************************************/

int search_for_distance_match_in_distanceseenlist(int distance,int start)

	/* tests for exact same distance or a difference in distance 
	   of less equal 2; finding a match means ruling 
	   out redoing the alignment with that distance */

	/* later, should include a test based on multiples where the actual
	   score and the theoretical best score (all matches) are compared; if
	   the actual score is close to the theoretical score, then don't redo
	   the alignment */     

{
	struct distancelistelement *entry, *entrylast, *temp;
	int absdiff;
	entry=Distanceseenlist->next;
	entrylast=Distanceseenlist;
	while (entry!=NULL)
	{

		/*	if (entry->end<min(start-distance,start-Min_Distance_Window))*/ /***???***/
		/* This change added 11/1/01 to allow repeats that are smaller than Min_Distance_Window to
		   be detected *following* another repeat.  The problem was that the first repeat could block
		   any repeat that followed it, but ended within Min_Distance_Window of the end of the first
		   repeat.  By using Minscore/Alpha, the first repeat can't block anything past its own end. */

		/*    if (entry->end<min(start-distance,start-Minscore/Alpha)+1)  */

		if (entry->end<min(start-distance,max(start-Minscore/Alpha,start-Min_Distance_Window))+1) 

			/* +1 added to correct for the < instead of <= */

			/* changed start to start-distance */
		{
			entrylast->next=entry->next;
			temp=entry;
			entry=entrylast;
			free(temp);
		}
		else 
		{
			if(distance>=entry->distance) absdiff=distance-entry->distance;
			else absdiff=entry->distance-distance;
			if(absdiff==0)
			{
				return(TRUE);
			}
			/* new idea 2/11/05 use an extra field called changed_from_distance to store the old distance
			   when you have a consensus.  In search_for_distance_on_distanceseen_list, if distance doesn't
			   match, but is close, say within 5, then also look at changed_from_distance for a match.
			 */  
			else if((absdiff<=5) && (distance==entry->changed_from_distance))
			{
				return(TRUE);
			}



		}
		/* go to next entry */
		entrylast=entry;  
		entry=entry->next;
	}
	return(FALSE);

}
