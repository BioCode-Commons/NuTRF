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

/* Allocation of the shared work arrays. */

#include "core/globals.h"
#include "core/report.h"
#include "core/arrays.h"

#define new1Darrayfunc(type,functionname,length)\
	type *functionname(int length)\
{\
	type *objptr=calloc(length,sizeof(*objptr));\
	if(objptr==NULL)\
	{\
		paramset.endstatus = "functionname: Out of memory!";\
		return NULL;\
	}\
	return(objptr);\
}
new1Darrayfunc(char,newAlignPairtext,length)

new1Darrayfunc(char,newLine,length)

new1Darrayfunc(int,newAlignPairindex,length)

new1Darrayfunc(int,newTags,length)

/*******************************************************************/

/**************************  init_index()  *************************/

/*******************************************************************/

void init_index()

{
	/* index has 256 entries so that finding the entries for A, C, G and T */
	/* require no calculation */

	Index=(int *)calloc(256,sizeof(int));
	if(Index==NULL){trf_message("\nInit_index: Out of memory!");exit(-1);}
	Index['A']=0;
	Index['C']=1;
	Index['G']=2;
	Index['T']=3;
}
