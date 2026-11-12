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

#include "core/globals.h"
#include "core/report.h"
#include "run/fasta.h"

// how much memory to allocate at first when loading a sequence
// (added by Eugene Scherba, 2010-02-16)
#define MEMBLOCK (10*1024*1024)

/*******************************************************
 * this is an older version of LoadSequenceFromFile
 * it was replaced because the older version bactracked in the
 * input file after determining sequence length which made the
 * program fail when input stream from stdin had to be read.
 *
 * new version written on 2010-02-18 by Eugene Scherba
 ******************************************************/
int LoadSequenceFromFileBenson(FASTASEQUENCE *pseq,FILE* fp)
{
	int letter,i,pos1,length,next;
	char *ptext;

	// read the FASTA '>' symbol
	letter = getc(fp);
	if(letter!='>') return -1; // invalid format

	// read name and description text
	for(i=0,ptext=pseq->name;i<(MAXSEQNAMELEN-1);i++,ptext++)
	{
		letter = getc(fp);
		if(letter==EOF) break;
		*ptext = letter;
		if(*ptext==10||*ptext==13)
		{
			break;
		}
	}
	*ptext='\0';
	// if line was not read completely flush the rest
	if(i==(MAXSEQNAMELEN-1))
	{
		letter = 0;
		while(letter!=13&&letter!=10&&letter!=EOF) letter = getc(fp);
	}

	// get the length of the sequence
	pos1 = ftell(fp);
	length=0;
	for(;;)
	{
		letter = getc(fp);
		if(letter==EOF||letter=='>') break;
		length++;
	}

	// allocate memory including white space
	pseq->sequence = (char*) malloc(sizeof(char)*(length+1));
	if(pseq->sequence==NULL) return -1;

	// read sequence into buffer
	fseek(fp,pos1,SEEK_SET);
	pseq->length=0;
	ptext = pseq->sequence;
	for(i=0;i<26;i++) pseq->composition[i]=0;
	for(;;)
	{
		// get a character from file
		letter = getc(fp);
		// break if end of file
		if(letter==EOF) { next = 0; break; }
		// break if another sequence found
		if(letter=='>') { next = 1; ungetc('>',fp); break; }
		*ptext = letter;
		// if character is in range of alpha characters
		if(*ptext>='A'&&*ptext<='z') // in alpha range
		{
			if(*ptext<='Z') // in upper case range
			{
				pseq->composition[*ptext-'A']++;
				ptext++;
				pseq->length++;
			}
			else if(*ptext>='a') // in lower case range
			{
				// make upper case
				(*ptext)+=('A'-'a');

				pseq->composition[*ptext-'A']++;
				ptext++;
				pseq->length++;
			}
		}
	}
	// terminate sequence text as a string
	*ptext='\0';

	// compute member
	pseq->nucleotides =
		pseq->composition['A'-'A'] + pseq->composition['C'-'A'] +
		pseq->composition['G'-'A'] + pseq->composition['T'-'A'];

	return next;
}

/*******************************************************
 *   loads a sequence from an open file. If the routine
 *   comes to another sequence header while reading a
 *   sequence the sequence returns 1 to indicate that
 *   more sequences are present. Otherwise the routine
 *   returns 0 to indicate EOF or -1 to indicate error.
 *   The member sequence must be NULL before the routine
 *   iscalled. The calling function must free the allo-
 *   cated memory after use.
 ********************************************************/

int LoadSequenceFromFileEugene(FASTASEQUENCE *pseq, FILE* fp)
{
	int i, j, c;
	int next = -1;	// whether a next sequence was encountered
	char *ptemp;
	char to_upper;

	// read the FASTA '>' symbol
	c = getc(fp);
	if ((char)c != '>' || c == EOF) return -1; // invalid format

	// read name and description text
	for(i = 0; i < MAXSEQNAMELEN - 1; i++) {
		c = getc(fp);
		if (c == 10 || c == 13) {
			break;
		} else if (c == EOF) {
			PrintError("FASTA input terminated too early");
			return -1;
		} else {
			pseq->name[i] = (char)c;
		}
	}
	pseq->name[i] = '\0';

	// if line was not read completely flush the rest
	if (i == MAXSEQNAMELEN - 1) {
		c = 0;
		while (c != 13 && c != 10 && c != EOF) c = getc(fp);
	}

	// read sequence into buffer
	pseq->sequence = NULL;
	for (i = 0; i < 26; i++) pseq->composition[i] = 0;
	to_upper = 'A' - 'a';

	i = 0;
	for (j = 0; c != EOF && (char)c != '>'; j += MEMBLOCK) {
		if ((ptemp = realloc(pseq->sequence, sizeof(char)*(i + MEMBLOCK + 1))) == NULL) {
			PrintError("Insufficient memory");
			return -1;
		}
		pseq->sequence = ptemp;
		for(i = j; i < j + MEMBLOCK;) {
			c = getc(fp);			// get a character from file
			if ((char)c >= 'A' && (char)c <= 'Z') { // in upper-case range of alpha characters
				pseq->sequence[i] = (char)c;
				pseq->composition[c - (int)'A']++;
				i++;
			} else if ((char)c >= 'a' && (char)c <= 'z') {	// in lower-case range of alpha characters
				pseq->sequence[i] = (char)c + to_upper;	// make upper-case
				pseq->composition[c - (int)'a']++;
				i++;
			} else if ((char)c == '>') {	// break if another sequence found
				next = 1;
				ungetc('>', fp);
				break;
			} else if (c == EOF) {		// break if end of file
				next = 0;
				break;
			}
		}
	}
	pseq->length = i;		// set sequence length
	if (i > 0) {
		pseq->sequence[i] = '\0';	// terminate sequence text as a string
	}

	// compute member
	pseq->nucleotides =
		pseq->composition['A'-'A'] + pseq->composition['C'-'A'] +
		pseq->composition['G'-'A'] + pseq->composition['T'-'A'];

	return next;
}
