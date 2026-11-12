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

/* HTML repeat tables and alignment pages of one sequence. */

#include "core/globals.h"
#include "core/report.h"
#include "output/html.h"

/* max # of items in tables for extended output format*/
#define EO_MAX_TBL 120

void   CleanAlignments(IL * headptr, char * alignmentfile)
{
	char string1[260];
	char string2[260];
	char string3[260];
	char tempfile[264];
	FILE *al_fp, *tmp_fp;
	IL*  currptr;
	int moving; /* in loop moving=1 to copy strings to tempfile */

	/* execute even if no repeats in list */

	/* make name of temporary file */
	strcpy(tempfile,alignmentfile);
	strcat(tempfile,".tmp");

	/* open files */
	al_fp = fopen(alignmentfile, "r");
	if(al_fp==NULL)
	{
		PrintError("Unable to open alignment file for reading in CleanAlignments routine!");
		exit(-1);
	}
	tmp_fp= fopen(tempfile, "w");
	if(tmp_fp==NULL)
	{
		PrintError("Unable to open temp file for writing in CleanAlignments routine!");
		exit(-1);
	}

	moving = 1; /* starts by moving lines to temporary file */
	currptr=headptr;
	while( fgets(string1, 260, al_fp) != NULL )
	{
		if( string1[0]=='F' /* Gelfand Dec 15, 2013 */&& string1[1]=='o') /* Ocurrs when "Found at" is encountered */
		{
			fgets(string2, 260, al_fp);
			fgets(string3, 260, al_fp);
			if(currptr!=NULL)
			{
				if( strstr(string3, currptr->ref)!=NULL )
				{
					fputs(string1, tmp_fp);
					fputs(string2, tmp_fp);
					fputs(string3, tmp_fp);
					currptr = currptr->next;
					moving = 1;
					continue;
				}
			}
			moving = 0;
		}

		if(string1[0]=='D' /* Gelfand Dec 15, 2013 */&& string1[1]=='o') /* Occurs when "Done." is encountered */
		{
			fgets(string2, 260, al_fp);
			fputs(string1, tmp_fp);
			fputs(string2, tmp_fp);
			fputs("\n", tmp_fp);
			break;
		}

		if(moving) fputs(string1, tmp_fp);
	}

	fclose(al_fp);
	fclose(tmp_fp);

	remove(alignmentfile);
	rename(tempfile,alignmentfile);

	return;

}

void   BreakAlignments(IL * headptr, char * alignmentfile)
{
	FILE * al_fp;
	FILE * out_fp;
	char outfile[260];
	int nfiles;
	IL* currptr;
	int alignments;
	char headlns[30][200]; /* to hold the heading section of the file */
	char buffer[200]; /* holds one line of alignment data */
	int headcnt;
	char nextchar;
	int i,j;


	/* Find out how many alignments there are and how many files will
	   be needed */
	for( alignments=0,currptr=headptr; currptr!=NULL; currptr=currptr->next, alignments++);
	nfiles = alignments / EO_MAX_TBL ;
	if( (alignments % EO_MAX_TBL )> 0 ) nfiles++ ;
	if(nfiles==0) nfiles = 1; /* make sure at least one file is generated */

	/* if only one file will be needed just rename file */
	if(nfiles==1)
	{
		MakeFileName( outfile, alignmentfile,1);
		remove(outfile);
		rename(alignmentfile,outfile);
		return;    
	}

	/* get heading lines */
	al_fp = fopen(alignmentfile,"r");
	if(al_fp==NULL)
	{
		PrintError("Unable to open alignment file for reading in BreakAlignments routine!");
		exit(-1);
	}

	for (i=0,headcnt=0; i<30; i++)
	{
		/* find if next line is the "Found at" line */
		nextchar = getc(al_fp); 
		ungetc(nextchar,al_fp);
		if(nextchar=='F') break;
		fgets( headlns[i], 200, al_fp);
		headcnt++;
	}

	/* loop creating files */
	for(i=1;i<=nfiles;i++)
	{
		/* create name of ith file */
		MakeFileName( outfile, alignmentfile,i);
		/* open the file for writing */
		out_fp = fopen(outfile,"w");
		if(out_fp==NULL)
		{
			PrintError("Unable to open output file for writing in BreakAlignments routine!");
			exit(-1);
		}
		/* output heading */
		for(j=0;j<headcnt;j++) fputs(headlns[j],out_fp);
		/* output n of N identifier */
		//fprintf(out_fp,"File %d of %d\n\n",i+1,nfiles);
		fprintf(out_fp,"File %d of %d\n\n",i,nfiles);

		/*move alignments to the current output file*/
		for(j=0; j<EO_MAX_TBL; j++)
		{
			/* copy the first line of the alignment */
			fgets(buffer, 200, al_fp);
			fputs(buffer, out_fp);
			/* copy successive lines checking for the "Found at" and
			   the "Done." lines */
			while(1)
			{
				nextchar = getc(al_fp);
				ungetc(nextchar,al_fp);
				if(nextchar=='F'||nextchar=='D') break;

				// added check for NULL, Gelfand Dec 15, 2013 because of suspicion this may under 
				// certain conditions enter an infinite loop. 
				if (NULL==fgets(buffer, 200, al_fp)) { nextchar='D'; break; }

				fputs(buffer, out_fp);
			}
			if(nextchar=='F') continue;
			if(nextchar=='D') break;
		}

		/* Output closing Lines */
		fprintf(out_fp,"\nDone.\n</PRE></BODY></HTML>\n");

		/*close file*/
		fclose(out_fp);
	}

	/*close original file*/
	fclose(al_fp);

	/* delete original file */
	remove(alignmentfile);

	return;
}

void    MakeFileName(char* newname, char* oldname, int tag)
{
	char newext[20];
	char oldext[10];
	int numindex;
	int count;

	/*copy oldname to newname buffer*/
	strcpy(newname,oldname);

	/*find position of last number in name*/
	for(count=0, numindex=0; newname[count]!='\0';count++)
	{
		if(isdigit((int) newname[count])) numindex = count;
	}


	/*get old extension */
	strcpy(oldext,&newname[numindex+1]);

	/* create new extension based on tag*/
	sprintf(newext,".%d%s",tag,oldext);

	/* copy newext in place over old newname */
	strcpy(&newname[numindex+1], newext);

	return;
}



void    OutputHTML(IL * headptr, char * tablefile, char *alignmentfile)
{
	FILE* fp;
	IL* currptr;
	int i,j;
	int alignments;
	int nfiles;
	char outfile[260];
	char linkfile[260];
	char namebuffer[260];


	/* find out how many elements there are and how
	   many files will be needed */
	for(alignments=0,currptr=headptr;
			currptr!=NULL;
			currptr=currptr->next,alignments++);
	nfiles = alignments/EO_MAX_TBL;
	if((alignments%EO_MAX_TBL)>0)nfiles++;
	if(nfiles==0) nfiles=1; /* make sure at least one file is generated */

	/* loop creating files */
	currptr=headptr;
	for(i=1;i<=nfiles;i++)
	{
		/* create name of ith file */
		MakeFileName( outfile, tablefile, i);

		/* create name of file where the alignments are found */
		MakeFileName( linkfile, alignmentfile, i);

		/* open the file for writing */
		fp = fopen(outfile,"w");
		if(fp==NULL)
		{
			PrintError("Unable to open output file for writing in OutputHTML routine!");
			exit(-1);
		}

		/* output heading */
		OutputHeading(fp, outfile,alignmentfile);

		/* print links to other tables */
		fprintf(fp,"\n<P><PRE>Tables:   ");
		for(j=1;j<=nfiles;j++)
		{
			/* output a link for all but the current page */ 
			if(j!=i)
			{
				MakeFileName(namebuffer,tablefile,j);
				fprintf(fp,"<A HREF=\"%s\" target=\"_self\">%d</A>   ", namebuffer, j);
			}
			else
			{
				fprintf(fp,"%d   ", j);
			}
			if(j%16==0&&j<nfiles) fprintf(fp,"\n          ");
		}

		/* printf "Table n of N" line */
		fprintf(fp,"\n\nThis is table  %d  of  %d  ( %d repeats found )\n</PRE>",i,nfiles,alignments);

		/*print help lines */
		fprintf(fp,"<PRE>\nClick on indices to view alignment\n");

#if defined(WINDOWSGUI)
		fprintf(fp,"</PRE><P>See <FONT COLOR=\"#0000FF\">Table Explanation</FONT> in Tandem Repeats Finder Help</P><BR>\n");
#elif defined(WINDOWSCONSOLE)
		fprintf(fp,"</PRE><A HREF=\"http://tandem.bu.edu/trf/trf.definitions.html#table\" target = \"explanation\">Table Explanation</A><BR><BR>\n");
#elif defined(UNIXGUI)
		fprintf(fp,"</PRE><P>See <FONT COLOR=\"#0000FF\">Table Explanation</FONT> in Tandem Repeats Finder Help</P><BR>\n");
#elif defined(UNIXCONSOLE)
		fprintf(fp,"</PRE><A HREF=\"http://tandem.bu.edu/trf/trf.definitions.html#table\" target = \"explanation\">Table Explanation</A><BR><BR>\n");
#endif


		/* print beginning of table */
		fprintf(fp,"<TABLE BORDER=1 CELLSPACING=0 CELLPADDING=0>\n");

		/* output rows of data */
		for(j=0;j<EO_MAX_TBL && currptr!=NULL;currptr=currptr->next,j++)
		{

			if(j%22==0) fprintf(fp,"<TR><TD WIDTH=140><CENTER>Indices</CENTER></TD><TD WIDTH=80><CENTER>Period<BR>Size </CENTER></TD><TD WIDTH=70><CENTER>Copy<BR>Number</CENTER></TD><TD WIDTH=70><CENTER>Consensus<BR>Size</CENTER></TD><TD WIDTH=70><CENTER>Percent<BR>Matches</CENTER></TD><TD WIDTH=70><CENTER>Percent<BR>Indels</CENTER></TD><TD WIDTH=60><CENTER>Score</CENTER></TD><TD WIDTH=40><CENTER>A</CENTER></TD><TD WIDTH=40><CENTER>C</CENTER></TD><TD WIDTH=40><CENTER>G</CENTER></TD><TD WIDTH=40><CENTER>T</CENTER></TD><TD WIDTH=70><CENTER>Entropy<BR>(0-2)</CENTER></TD></TR>\n");
			fprintf(fp,"<TR><TD><CENTER><A HREF=\"%s#%s\">%d--%d</A></CENTER></TD><TD><CENTER>%d</CENTER></TD><TD><CENTER>%1.1f</CENTER></TD><TD><CENTER>%d</CENTER></TD><TD><CENTER>%d</CENTER></TD><TD><CENTER>%d</CENTER></TD><TD><CENTER>%d</CENTER></TD><TD><CENTER>%d</CENTER></TD><TD><CENTER>%d</CENTER></TD><TD><CENTER>%d</CENTER></TD><TD><CENTER>%d</CENTER></TD><TD><CENTER>%1.2f</CENTER></TD></TR>\n",
					linkfile,currptr->ref,currptr->first,currptr->last,currptr->period,currptr->copies,currptr->size,currptr->matches,currptr->indels,currptr->score,currptr->acount,currptr->ccount,currptr->gcount,currptr->tcount,currptr->entropy);
		}

		/* if empty list print at least one heading */
		if(headptr==NULL) fprintf(fp,"<TR><TD WIDTH=140><CENTER>Indices</CENTER></TD><TD WIDTH=80><CENTER>Period<BR>Size </CENTER></TD><TD WIDTH=70><CENTER>Copy<BR>Number</CENTER></TD><TD WIDTH=70><CENTER>Consensus<BR>Size</CENTER></TD><TD WIDTH=70><CENTER>Percent<BR>Matches</CENTER></TD><TD WIDTH=70><CENTER>Percent<BR>Indels</CENTER></TD><TD WIDTH=60><CENTER>Score</CENTER></TD><TD WIDTH=40><CENTER>A</CENTER></TD><TD WIDTH=40><CENTER>C</CENTER></TD><TD WIDTH=40><CENTER>G</CENTER></TD><TD WIDTH=40><CENTER>T</CENTER></TD><TD WIDTH=70><CENTER>Entropy<BR>(0-2)</CENTER></TD></TR>\n");

		/* close table */
		fprintf(fp,"\n</TABLE>\n");

		/* if no repeats print message */
		if(headptr==NULL) fprintf(fp, "\nNo Repeats Found!<BR>");/*for empty list*/

		/* print links to other tables (again)*/
		fprintf(fp,"\n<P><PRE>Tables:   ");
		for(j=1;j<=nfiles;j++)
		{
			/* output a link for all but the current page */ 
			if(j!=i)
			{
				MakeFileName(namebuffer,tablefile,j);
				fprintf(fp,"<A HREF=\"%s\" target=\"_self\">%d</A>   ", namebuffer, j);
			}
			else
			{
				fprintf(fp,"%d   ", j);
			}
			if(j%16==0&&j<nfiles) fprintf(fp,"\n          ");
		}
		fprintf(fp,"\n</PRE>");

		if(i==nfiles) fprintf(fp,"<P>The End!\n");

		fprintf(fp,"\n</BODY></HTML>\n");
		fclose(fp);
	}

	return;
}

void    OutputHeading(FILE* fp, char * tablefile, char *alignmentfile)
{
	/* output fixed (old) heading */
	fprintf(fp,"<HTML><HEAD><TITLE>%s</TITLE><BASE TARGET=\"%s\"></HEAD><BODY bgcolor=\"#FBF8BC\"><BR><PRE>Tandem Repeats Finder Program written by:</PRE><PRE><CENTER>Gary Benson<BR>Program in Bioinformatics<BR>Boston University<BR>Version %s<BR></CENTER>",tablefile, alignmentfile, versionstring);
	fprintf(fp,"\nPlease cite:\nG. Benson,\n\"Tandem repeats finder: a program to analyze DNA sequences\"\nNucleic Acid Research(1999)\nVol. 27, No. 2, pp. 573-580.\n");

	fprintf(fp,"\n%s",hsequence);
	fprintf(fp,"%s",hparameters);
	fprintf(fp,"%s</PRE>\n",hlength);

	return;
}
