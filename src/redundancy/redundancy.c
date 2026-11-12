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

/* Removal of repeats by period size and of repeats made redundant by others. */

#include <stdlib.h>
#include "redundancy/redundancy.h"

/************ RemoveBySize() **********************************************/

IL*     RemoveBySize(IL * headptr, int maxsize)
{
	IL* currptr;
	IL* prevptr;

	/* loop thru list removing all elements with period > maxsize */
	for(currptr=headptr; currptr!=NULL;)
	{
		if(currptr->period>maxsize) /* remove */
		{
			if(currptr==headptr)
			{
				headptr=headptr->next;
				free(currptr->pattern);
				free(currptr);
				currptr=headptr; 
			}
			else
			{
				prevptr->next = currptr->next;
				free(currptr->pattern);
				free(currptr);
				currptr = prevptr->next;
			}

		}
		else
		{
			prevptr= currptr;
			currptr= currptr->next;
		}
	}

	return headptr;

}

IL*     RemoveRedundancy(IL * headptr)
{
	IL* iptr;       /* first pointer of a pair being examined */
	IL* jptr;       /* second pointer of a pair being examined */
	IL* previptr;   /* points to element before iptr */
	IL* prevjptr;   /* points to element before jptr */
	int overlap;    /* overlap of two intervals */
	int iinterval;  /* first pointer index interval */
	int jinterval;  /* second pointer index interval */

	if(headptr==NULL) return headptr; /* return if no elements */
	if(headptr->next==NULL) return headptr; /* return if one element only */

	iptr = headptr;       /* initialize to start at head of list */

	while(iptr!=NULL) /* loop from start to end of list*/
	{
		//iinterval = iptr->last - iptr->first + 1;

		jptr=iptr->next;
		prevjptr=iptr;
		while(jptr!=NULL) /* loop until not enough overlap */
		{
			iinterval = iptr->last - iptr->first + 1;
			jinterval = jptr->last - jptr->first + 1;
			overlap = IntervalOverlap(iptr, jptr);
			if(overlap==0) break;

			/* if neither overlap satisfies minimum of 90% break loop*/
			/*if( (overlap/(double)iinterval<0.9) && */
			/*  (overlap/(double)jinterval<0.9) ) break; */

			/* overlap in relation to iinterval */
			if( !(overlap/(double)iinterval<0.9) && IsRedundant(iptr,jptr))
			{
				/* remove iptr and break from inner loop*/
				if(iptr==headptr)
				{
					headptr=headptr->next;
					free(iptr->pattern);
					free(iptr);
					iptr=headptr;
					jptr=iptr->next;
					prevjptr=iptr;
					continue;
				}
				else
				{
					previptr->next = iptr->next;
					free(iptr->pattern);
					free(iptr);
					iptr = previptr->next;
					jptr=iptr->next;
					prevjptr=iptr;
					continue;
				}
			}

			/* overlap in relation to jinterval */
			if( !(overlap/(double)jinterval<0.9) && IsRedundant(jptr,iptr))
			{
				/* remove jptr and continue next iteration of inner loop*/
				prevjptr->next = jptr->next;
				free(jptr->pattern);
				free(jptr);
				jptr = prevjptr->next;
				continue;
			}

			/* update and continue to next iteration of inner loop */
			prevjptr = jptr;
			jptr = jptr->next;

		}
		previptr = iptr;
		iptr = iptr->next;
	}
	return headptr;
}

int     IntervalOverlap(IL* iptr, IL* jptr)
{
	int beg,end, overlap;
	beg = ( iptr->first > jptr->first ) ? iptr->first : jptr->first;
	end = ( iptr->last  < jptr->last  ) ? iptr->last  : jptr->last;

	overlap = end - beg + 1;    
	return (( overlap > 0 ) ? overlap:0);
}

/* this fuction retuns 1 if iptr is redundant in reference to jptr */
int     IsRedundant(IL* iptr,IL* jptr)
{

	if( (iptr->period>jptr->period)&&
			(iptr->period%jptr->period==0)&&
			(iptr->score<= 1.1*jptr->score)  ) return 1;
	if( (iptr->period==jptr->period)&&
			(iptr->score<= jptr->score)  ) return 1;

	return 0;
}    
