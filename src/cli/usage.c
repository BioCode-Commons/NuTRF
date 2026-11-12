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

/* Usage text and version banner. */

#include <stdio.h>
#include "core/platform.h"
#include "core/version.h"
#include "cli/usage.h"

const char *usage = "\n\nPlease use: %s File Match Mismatch Delta PM PI Minscore MaxPeriod [options]\n"
"\nWhere: (all weights, penalties, and scores are positive)"
"\n  File = sequences input file"
"\n  Match  = matching weight"
"\n  Mismatch  = mismatching penalty"
"\n  Delta = indel penalty"
"\n  PM = match probability (whole number)"
"\n  PI = indel probability (whole number)"
"\n  Minscore = minimum alignment score to report"
"\n  MaxPeriod = maximum period size to report. Must be between 1 and 2000, inclusive"
"\n  [options] = one or more of the following:"
"\n        -m        masked sequence file"
"\n        -f        flanking sequence"
"\n        -d        data file"
"\n        -h        suppress html output"
"\n        -r        no redundancy elimination"
"\n        -l <n>    maximum TR length expected (in millions) (eg, -l 3 or -l=3 for 3 million)"
"\n                  Human genome HG38 would need -l 6"
#if (defined(UNIXGUI)+defined(UNIXCONSOLE))>=1
"\n        -ngs      more compact .dat output on multisequence files, returns 0 on success."
"\n                  Output is printed to the screen, not a file. You may pipe input in with"
"\n                  this option using - for file name. Short 50 flanks are appended to .dat"
"\n                  output."
#endif
"\n"
"\nSee more information on the TRF Unix Help web page: https://tandem.bu.edu/trf/trf.unix.help.html"
"\n"
"\nNote the sequence file should be in FASTA format:"
"\n"
"\n>Name of sequence"
"\naggaaacctgccatggcctcctggtgagctgtcctcatccactgctcgctgcctctccag"
"\natactctgacccatggatcccctgggtgcagccaagccacaatggccatggcgccgctgt"
"\nactcccacccgccccaccctcctgatcctgctatggacatggcctttccacatccctgtg"
"\n"
"\n"
"\nTandem Repeats Finder"
"\nCopyright (C) 1999-2020 Gary Benson"
"\n"
"\nThis program is free software: you can redistribute it and/or modify"
"\nit under the terms of the GNU Affero General Public License as"
"\npublished by the Free Software Foundation, either version 3 of"
"\nthe License, or (at your option) any later version."
"\n";

void    PrintBanner(void)
{
	fprintf(stderr,"\nTandem Repeats Finder, Version %s", versionstring);
	fprintf(stderr,"\nCopyright (C) Dr. Gary Benson 1999-2012. All rights reserved.\n");

	return;
}
