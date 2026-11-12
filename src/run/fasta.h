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

/* FASTA record loading. */

#ifndef TRF_RUN_FASTA_H
#define TRF_RUN_FASTA_H

#include <stdio.h>
#include "core/types.h"

int LoadSequenceFromFileEugene(FASTASEQUENCE *pseq,FILE* fp); /* may use stdin and better file reading for handling large files over 2GB */
int LoadSequenceFromFileBenson(FASTASEQUENCE *pseq,FILE* fp); /* old function, uses filepos, 32bit version of this would not process a file over 2GB properly */

#endif
