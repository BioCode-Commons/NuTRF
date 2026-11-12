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

/* Command-line parsing into the global parameter set. */

#ifndef TRF_CLI_ARGUMENTS_H
#define TRF_CLI_ARGUMENTS_H

/* Fills paramset from the command line. Prints a message and exits on a
   usage request, a version request, or an invalid argument. */
void ParseArguments(int ac, char** av);

char* GetNamePartAddress(char* name);

#endif
