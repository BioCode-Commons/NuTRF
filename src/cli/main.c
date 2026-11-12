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


/**************************************************************
 * main.c :   Command line interface to the Tandem Repeats Finder
 *
 ***************************************************************/


#include <stdio.h>
#include "core/globals.h"
#include "cli/arguments.h"
#include "run/control.h"

int main(int ac, char** av)
{
	ParseArguments(ac, av);

	/* call the routine that controls execution */
	TRFControlRoutine();
	
	/* Check the status by looking at the output count and the endstatus
	   members of the paramset struct */
	if (paramset.outputcount == 0 && !paramset.endstatus) {
		/* If no TRs were found, print informative message to STDERR */
		printf("No TRs found. Exiting...\n");
		return 0;
	}
	else if (paramset.endstatus) {
		printf("Error processing input: %s\n", paramset.endstatus);
		return 1;
	}
	else {
		return 0;
	}
}
