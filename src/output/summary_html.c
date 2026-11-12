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

/* Summary page that lists the sequences of a multi-sequence input. */

#include "core/globals.h"
#include "output/summary_html.h"

void    WriteSummaryHeading(FILE* fp)
{
	/* start output of sumary file */
	fprintf(fp,"<HTML>");
	fprintf(fp,"<HEAD>");
	fprintf(fp,"<TITLE>Output Summary</TITLE>");
	/* fprintf(fp,"<BASE TARGET=\"Table\">"); */
	fprintf(fp,"</HEAD>");
	fprintf(fp,"<BODY bgcolor=\"#FBF8BC\">");
	fprintf(fp,"<PRE>");

	fprintf(fp,"\nTandem Repeats Finder Program written by:<CENTER>");
	fprintf(fp,"\nGary Benson");
	fprintf(fp,"\nProgram in Bioinformatics");
	fprintf(fp,"\nBoston University");
	fprintf(fp,"\nVersion %s</CENTER>\n", versionstring);
	fprintf(fp,"\nPlease cite:\nG. Benson,\n\"Tandem repeats finder: a program to analyze DNA sequences\"\nNucleic Acid Research(1999)\nVol. 27, No. 2, pp. 573-580.\n");


	fprintf(fp,"\n\n<B>Multiple Sequence Summary</B>\n\n");
	fprintf(fp,"Only sequences containing repeats are shown!\n\n");
	fprintf(fp,"Click on sequence description to view repeat table.\n\n");

	/* print beginning of table */
	fprintf(fp,"<TABLE BORDER=1 CELLSPACING=0 CELLPADDING=0>\n");
	fprintf(fp,"<TR><TD WIDTH=80><CENTER>Sequence\nIndex</CENTER></TD>"
			"<TD WIDTH=400><CENTER>Sequence\nDescription</CENTER></TD>"
			"<TD WIDTH=80><CENTER>Number of\nRepeats</CENTER></TD>"
			"</TR>\n");
}

/* one table row for a sequence that contains repeats */
void    WriteSummaryRow(FILE* fp, int ordinal, char* tablefile, char* name, int repeats)
{
	fprintf(fp,"<TR><TD><CENTER>%d</CENTER></TD>"
			"<TD><CENTER><A TARGET=\"%s\" HREF=\"%s\">%s</A>"
			"</CENTER></TD><TD><CENTER>%d</CENTER></TD></TR>",
			ordinal,tablefile,tablefile,name,repeats);
}

void    WriteSummaryEnd(FILE* fp, int foundsome)
{
	/* close table and html body */
	fprintf(fp,"\n</TABLE>\n");
	if(!foundsome)
	{
		fprintf(fp, "\nNo Repeats Found!<BR>");
	}
	fprintf(fp,"\n</BODY></HTML>\n");
}
