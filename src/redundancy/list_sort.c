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

/* Stable sorts of the repeat list. */

#include <stddef.h>
#include "redundancy/list_sort.h"

/* Stable bottom-up mergesort, relinking only: no allocation or payload writes.
 * by_count selects the original single-field comparator. Equal keys always
 * take the left node, retaining the original adjacent-swap sorts' ordering.
 */
static IL* trf_sort_list(IL *headptr, int by_count)
{
	IL *left, *right, *node, *tail, *remaining;
	size_t width, leftsize, rightsize, merges, i;

	if(headptr==NULL || headptr->next==NULL) return headptr;

	/* Retain the linear best case for already ordered and all-equal lists. */
	for(left=headptr; left->next!=NULL; left=left->next)
	{
		if((by_count ? left->count : left->first) >
			(by_count ? left->next->count : left->next->first)) break;
	}
	if(left->next==NULL) return headptr;

	for(width=1;;width*=2)
	{
		remaining=headptr;
		headptr=tail=NULL;
		merges=0;
		while(remaining!=NULL)
		{
			merges++;
			left=right=remaining;
			leftsize=0;
			for(i=0; i<width && right!=NULL; i++)
			{
				leftsize++;
				right=right->next;
			}
			rightsize=width;
			while(leftsize>0 || (rightsize>0 && right!=NULL))
			{
				if(leftsize==0)
				{
					node=right;
					right=right->next;
					rightsize--;
				}
				else if(rightsize==0 || right==NULL ||
					(by_count ? left->count<=right->count :
						left->first<=right->first))
				{
					node=left;
					left=left->next;
					leftsize--;
				}
				else
				{
					node=right;
					right=right->next;
					rightsize--;
				}
				if(tail!=NULL) tail->next=node;
				else headptr=node;
				tail=node;
			}
			remaining=right;
		}
		tail->next=NULL;
		/* A single merged pair covers the entire list. This stops before
		 * width can overflow for any realizable list of IL nodes.
		 */
		if(merges<=1) return headptr;
	}
}

IL*     SortByIndex(IL * headptr)
{
	return trf_sort_list(headptr, 0);
}

IL*     SortByCount(IL * headptr)
{
	return trf_sort_list(headptr, 1);
}
