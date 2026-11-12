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

/* Messages to the alignment file, the error stream, and the progress indicator. */

#ifndef TRF_CORE_REPORT_H
#define TRF_CORE_REPORT_H

void trf_message( char *format, ... );
void PrintError(char* errortext);
void PrintProgress(char* progresstext);
void SetProgressBar(void);

#ifdef UNIXGUI
extern void set_progress_bar( double fraction );
#endif

#endif
