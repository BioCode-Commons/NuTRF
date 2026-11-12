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

/* Masked sequence file. */

#include "core/globals.h"
#include "core/report.h"
#include "output/mask.h"

void    MakeMaskedFile(IL* headptr,int masked,char*  Sequence,char* maskfile)
{
	int count,printcr;
	int masker;
	FILE *fp;
	IL* lpointer;

	if(masked)
	{
		fp = fopen(maskfile,"w");
		if(fp==NULL)
		{
			PrintError("Unable to open output file for writing in MakeMaskedFile routine!");
			exit(-1);
		}

		/* Ouput sequence description from global variable to file*/
		fprintf(fp,">%s",&hsequence[10]);

		for(lpointer=headptr;lpointer!=NULL;lpointer=lpointer->next)
		{
			for(masker=lpointer->first; masker<=lpointer->last; masker++)
				Sequence[masker]='N';
		}
		printcr=0;
		for(count=1;Sequence[count]!='\0';count++)
		{
			fputc(Sequence[count], fp);
			printcr++;
			if(printcr>=60)
			{
				printcr=0;
				fputc('\n', fp);
			}
		}
		fputc('\n', fp);
		fputc('\n', fp);
		fclose(fp);
	}

	return;
}
