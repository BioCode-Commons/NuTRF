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

/* Control of one run: loops over the FASTA records of the input, analyzes
   each sequence, and assembles the data, masked, and summary files. */

#include <errno.h>
#include "core/globals.h"
#include "core/report.h"
#include "core/repeat_list.h"
#include "output/data_file.h"
#include "output/mask.h"
#include "output/summary_html.h"
#include "run/control.h"
#include "run/fasta.h"
#include "run/sequence.h"

/***********************************************
 *   This routine can act on a multiple-sequence
 *   file and calls TRF() routine as many times
 *   as it needs to.
 ************************************************/

void TRFControlRoutine(void)
{
	FILE *srcfp,*outmfp,*destmfp,*destdfp = NULL;
	char  source[_MAX_PATH],input[_MAX_PATH],outm[_MAX_PATH],
	      prefix[_MAX_PATH],destm[_MAX_PATH],destd[_MAX_PATH],
	      desth[_MAX_PATH],paramstring[_MAX_PATH],outh[_MAX_PATH];
	int a,i,loadstatus,foundsome=0;
	char  line[1000];
	FILE  *desthfp;
	FASTASEQUENCE seq;


	/* save names locally so they can be replaced later */
	strcpy(source,paramset.inputfilename);
	strcpy(prefix,paramset.outputprefix);


	/* open input file for reading */
	if (paramset.use_stdin) {
		srcfp = stdin;
	} else {
		srcfp = fopen(source,"rb");
		if(srcfp==NULL)
		{
			sprintf(line,"Error opening %s: %s\n",source, strerror(errno));
			PrintError(line);
			paramset.endstatus = "Bad filename";
			paramset.running = 0;
			// return CTRL_BADFNAME;
			return;
		}
	}

	/* get the first sequence */
	if (paramset.ngs != 1) PrintProgress("Loading sequence...");


	loadstatus = LoadSequenceFromFileEugene(&seq,srcfp);


	if(loadstatus<0)
	{
		PrintError("Could not load sequence. Empty file or bad format.");
		paramset.endstatus = "Bad format."; /* ok for now */
		paramset.running = 0;
		fclose(srcfp);
		// return CTRL_BADFORMAT;
		return;
	}


	/* generate the parameter string to be used in file names */
	sprintf(paramstring,"%d.%d.%d.%d.%d.%d.%d",
			paramset.match,paramset.mismatch,paramset.indel,
			paramset.PM,paramset.PI,paramset.minscore,paramset.maxperiod);


	sprintf(hparameters,"Parameters: %d %d %d %d %d %d %d\n",
			paramset.match,paramset.mismatch,paramset.indel,
			paramset.PM,paramset.PI,paramset.minscore,paramset.maxperiod);



	/* based on number of sequences in file use different approach */
	if(loadstatus==0) /* only one sequence in file */
	{

		sprintf(hsequence,"Sequence: %s\n",seq.name);
		sprintf(hlength,"Length:  %d",seq.length);

		paramset.multisequencefile = 0;
		paramset.sequenceordinal = 1;
		/* call trf and return */
		counterInSeq=0;
		TRF(&seq);

		if (paramset.endstatus) {
			return;
		}

		if(paramset.datafile)
		{


			/* Added by Yevgeniy Gelfand on Jan 27, 2010  */
			/* To have smaller sequences not send results */
			/* to disc to improve performance             */



	                if (paramset.ngs) {
        	                destdfp = stdout;
	                } else {
        	                sprintf(destd,"%s.%s.dat",prefix,paramstring);
	                        destdfp = fopen(destd,"w");
        	                if(destdfp==NULL)
                	        {
                        	        PrintError("Unable to open data file for writing in TRFControlRoutine routine!");
	                                exit(-1);
        	                }
	                }


			WriteDataHeading(destdfp);
			WriteDataRecords(destdfp,seq.name);

		}


		/* masked file moved here so Sequence is not "ruined" by Ns for .dat output */
		{		
			char maskstring[_MAX_PATH];
			sprintf(maskstring,"%s.%s.mask",paramset.outputprefix,paramstring);
			MakeMaskedFile(GlobalIndexList, paramset.maskedfile, Sequence, maskstring);
		}


		FreeList(GlobalIndexList);
		GlobalIndexList = NULL;
		GlobalIndexListTail = NULL;


		free(seq.sequence);
		fclose(srcfp);

		if(destdfp) { fclose(destdfp); destdfp = NULL; }

		paramset.endstatus = NULL;
		paramset.running = 0;
		// return CTRL_SUCCESS;
		return;
	}
	paramset.multisequencefile = 1;
	paramset.sequenceordinal = 1;

	/*********************************************************
	 *   if there are more files need to produce sumary-style
	 *   output.
	 **********************************************************/

	/* generate the parameter string to be used in file names */
	sprintf(paramstring,"%d.%d.%d.%d.%d.%d.%d",
			paramset.match,paramset.mismatch,paramset.indel,
			paramset.PM,paramset.PI,paramset.minscore,paramset.maxperiod);

	/* open sumary table file */
	sprintf(desth,"%s.%s.summary.html",prefix,paramstring);
	if (!paramset.HTMLoff) {
		desthfp = fopen(desth,"w");
		if(desthfp==NULL)
		{
			PrintError("Unable to open summary file for writing in TRFControlRoutine routine!");
			exit(-1);
		}
	}

	/* open masked file if requested */
	if(paramset.maskedfile)
	{
		sprintf(destm,"%s.%s.mask",prefix,paramstring);
		destmfp = fopen(destm,"w");
		if(destmfp==NULL)
		{
			PrintError("Unable to open masked file for writing in TRFControlRoutine routine!");
			exit(-1);
		}
	}
	/* open datafile if requested */
	if(paramset.datafile)
	{
		if (paramset.ngs) {
			destdfp = stdout;
		} else {
			sprintf(destd,"%s.%s.dat",prefix,paramstring);
			destdfp = fopen(destd,"w");
			if(destdfp==NULL)
			{
				PrintError("Unable to open data file for writing in TRFControlRoutine routine!");
				exit(-1);
			}
		}
	}

	if (!paramset.HTMLoff) {
		WriteSummaryHeading(desthfp);
	}

	/******************************************
	 *   process every sequence in file
	 *******************************************/
	i=1;
	for(;;)
	{

		sprintf(hsequence,"Sequence: %s\n",seq.name);
		sprintf(hlength,"Length:  %d",seq.length);

		/* set the prefix to be used for naming of output */
		sprintf(input,"%s.s%d",prefix,i);
		strcpy(paramset.inputfilename,input);
		strcpy(paramset.outputprefix,input);

		/* call the tandem repeats finder routine */
		counterInSeq=0;
		TRF(&seq);

		if(paramset.datafile)
		{

			/* only for the first one write the header */
			if (i==1) {
				WriteDataHeading(destdfp);
			}
			WriteDataRecords(destdfp,seq.name);

		}

		if (!paramset.HTMLoff) {
			/* print table rows based on repeat count */
			sprintf(outh,"%s.%s.1.html",input,paramstring);
			if(paramset.outputcount>0)
			{
				/* print a table raw to the summary table */
				WriteSummaryRow(desthfp,i,outh,seq.name,paramset.outputcount);
				foundsome=1;
			}
			else
			{
				/* remove html files if no output in it */
				remove(outh);
				sprintf(line,"%s.%s.1.txt.html",input,paramstring);
				remove(line);
			}
		}


		/* masked file moved here so Sequence is not "ruined" by Ns for .dat output */
		{		
			char maskstring[_MAX_PATH];
			sprintf(maskstring,"%s.s%d.%s.mask",prefix,i,paramstring);
			MakeMaskedFile(GlobalIndexList, paramset.maskedfile, Sequence, maskstring);
		}


		/* append new output to destination files */
		if(paramset.maskedfile)
		{
			/* recreate the name of the masked sequence file */
			sprintf(outm,"%s.s%d.%s.mask",prefix,i,paramstring);
			outmfp = fopen(outm,"r");
			if(outmfp==NULL)
			{
				PrintError("Unable to open masked file for reading in TRFControlRoutine routine!");
				exit(-1);
			}
			/* copy until end of file */
			while(1)
			{
				a = getc(outmfp);
				if(a==EOF) break;
				putc(a,destmfp);
			}
			fclose(outmfp);

			/* remove intermediary file */
			remove(outm);
		}



		FreeList(GlobalIndexList);
		GlobalIndexList = NULL;
		GlobalIndexListTail = NULL;


		/* free the data associated with the sequence */
		free(seq.sequence);

		/* if more sequences load and repeat */
		if(loadstatus>0)
		{
			if (paramset.ngs != 1) PrintProgress("Loading sequence file...");

			loadstatus = LoadSequenceFromFileEugene(&seq,srcfp);

			paramset.sequenceordinal++;
			i++;
		}
		else
		{
			break;
		}
	}

	if (!paramset.HTMLoff) {
		WriteSummaryEnd(desthfp,foundsome);
	}

	/* close files */
	fclose(srcfp);
	if(paramset.maskedfile) fclose(destmfp);
	if(paramset.datafile) fclose(destdfp);

	if (!paramset.HTMLoff)
		fclose(desthfp);

	/* set output file name to the summary table */
	strcpy(paramset.outputfilename,desth);

	paramset.endstatus = NULL;
	paramset.running = 0;
	// return CTRL_SUCCESS;
	return;
}
