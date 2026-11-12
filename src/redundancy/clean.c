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

/****************************************************************
 *   Resolution of the repeat list of one sequence
 *****************************************************************
 *   The repeats found in a sequence are kept in a linked list.
 *   The general procedure is the following :
 *   i.      Remove pattern sizes greater than the maximum period.
 *   ii.     Get the list sorted by first index, keeping track of
 *           the original order.
 *   iii.    Remove redundant elements based on given criteria.
 *   iv.     Sort remaining elements back to original order.
 *   v.      In the alignment file remove sections that correspond
 *           to records deleted in i and iii.
 *   vi.     Break the alignment file into as many files as needed
 *           for the maximum number of alignments per file.
 *   vii.    Save all remaining elements in HTML format, creating
 *           as many files as necessary.
 *   The data file and the masked file are written afterwards by
 *   the control routine from the resolved list.
 *****************************************************************/

#include "core/globals.h"
#include "output/html.h"
#include "redundancy/clean.h"
#include "redundancy/list_sort.h"
#include "redundancy/redundancy.h"

void TRFClean( char * datafile, char* alignmentfile, char* tablefile,
		int maxsize, int data, int masked, char* Sequence,
		char* maskfile)
{
	IL *headptr=NULL,*currptr;
	int i;

	/* Added by Yevgeniy Gelfand on Jan 27, 2010  */
	/* To have smaller sequences not send results */
	/* to disc to improve performance             */
	headptr = GlobalIndexList;

	headptr = RemoveBySize(headptr, maxsize);

	headptr = SortByIndex(headptr);

	if (!paramset.redundoff) {
		headptr = RemoveRedundancy(headptr);
	}

	headptr = SortByCount(headptr);

	if(!paramset.HTMLoff)
	{

		CleanAlignments(headptr, alignmentfile);

		BreakAlignments(headptr, alignmentfile);

		OutputHTML(headptr, tablefile, alignmentfile);
	}

	/* Added by Yevgeniy Gelfand on Jan 27, 2010  */
	/* To have smaller sequences not send results */
	/* to disc to improve performance             */
	//MakeDataFile(headptr,datafile,data);

	//MakeMaskedFile(headptr, masked, Sequence, maskfile);

	/* update the global result */
	for(i=0,currptr=headptr;currptr!=NULL;i++,currptr=currptr->next);
	paramset.outputcount= i;

	/* Added by Yevgeniy Gelfand on Jan 27, 2010  */
	/* To have smaller sequences not send results */
	/* to disc to improve performance             */
	GlobalIndexList = headptr;

	return;
}
