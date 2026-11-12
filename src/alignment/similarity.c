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

/* Match and mismatch scores for every pair of byte values. */

#include "core/globals.h"
#include "core/report.h"
#include "alignment/similarity.h"

/* Jan 27, 2006, Gelfand, changed to use Similarity Matrix to avoid N matching itself */
/* This function may be called multiple times (for different match/mismatch scores) */
/*******************************************************************/
/**************************    init_sm()   *************************/
/*******************************************************************/
void init_sm(int match, int mismatch)

{
	static int sm_allocated = 0;
	static int previous_match, previous_mismatch;
	int i,j,*currint;

	/* The matrix is read-only between calls. WGS records normally reuse the
	 * same scoring parameters; rebuild only when either parameter changes. */
	if (sm_allocated && match == previous_match && mismatch == previous_mismatch)
		return;

	/* SM has 256*w56 entries to map MATCH-MISMATCH matrix */
	if (!sm_allocated) 
	{
		SM=(int *)calloc(256*256,sizeof(int));
		if(SM==NULL){trf_message("\nInit_sm: Out of memory!");exit(-1);}
		sm_allocated = 1;
	}

	/* generate 256x256 into matrix */
	for(i=0,currint=SM;i<=255;i++)
	{
		for(j=0;j<=255;j++,currint++)
		{
			*currint=mismatch;
		}
	}

	SM['A'*256+'A']=match;
	SM['C'*256+'C']=match;
	SM['G'*256+'G']=match;
	SM['T'*256+'T']=match;
	previous_match=match;
	previous_mismatch=mismatch;
}
