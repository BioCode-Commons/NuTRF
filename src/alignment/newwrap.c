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

/* Wraparound alignment over the whole pattern width, used for small periods. */

#include "core/globals.h"
#include "core/report.h"
#include "alignment/similarity.h"
#include "alignment/newwrap.h"


/* For nonpositive gap scores, the second circular pass can reuse each
 * first-pass cell. Its incoming score dominates the first-pass incoming
 * score; after a zero sentinel reset, the next first-pass cell has no
 * positive contribution from that reset predecessor. Keep the original
 * implementation for positive gap scores. Integer-overflow behavior is
 * subject to the same inherited score bounds as the original recurrence. */

/*******************************************************************/

/***************************** newwrap()  **************************/

/*******************************************************************/
static void newwrap_original(int start, int size, int consensuspresent)

{
	int g;
	register int *pup, *pdiag, *pcurr,pleft;
	int adjlength,adjmone,c,realr,end_of_trace,
	    maxscore,minrow,maxrow,maxcol,modstart,maxrealrow;
	char currchar;
	/* Feb 16, 2016 Yozen */
	unsigned int r;

	/* fill EC */
	if(consensuspresent)
		for(g=0;g<size;g++) EC[g]=Consensus.pattern[g];
	else
		for(g=0;g<size;g++) EC[g]=Sequence[start-size+g+1];


	/* backward wdp */
	maxscore=0;
	realr=start+1;
	r=maxwraplength;

	adjlength=size-1;
	adjmone=adjlength-1;
	pup=Up;
	pdiag=Diag;
	pcurr=&S[r][0];

	for (c=0;c<size;c++)
	{
		*pup=(*pdiag=(*pcurr=maxscore))+ Delta;
		pup++;
		pdiag++;
		pcurr++;
	}
	/* backwards */
	end_of_trace=FALSE;  
	while ((!end_of_trace) && (realr>1) && (r>0))
	{
		r--;
		realr--;
		Rows++;
		currchar=Sequence[realr];

		pcurr=&S[r][adjlength];
		pleft=Delta;        /* first pass. set S[r][0]=0 */
		pdiag=&Diag[adjlength];
		pup=&Up[adjlength];
		for (c=adjlength;c>=0;c--)
		{
			*pdiag+=match(currchar,EC[c]);
			pleft=(*pcurr=max4(0,*pdiag,*pup,pleft))+Delta;
			pcurr--;
			pdiag--;
			pup--;
		}

		end_of_trace=TRUE;      /* setup for encountering a break in the trace */
		/*second pass */
		pcurr=&S[r][adjlength];
		/* pleft set from first pass */
		pdiag=&Diag[adjlength];
		pup=&Up[adjlength];
		for (c=adjlength;c>=0;c--)
		{
			*pup=(pleft=(*pcurr=max4(0,*pdiag,*pup,pleft))+Delta);
			if((realr<=start-max(size,Min_Distance_Window))&&(*pcurr==0)) 
				/* if (*pcurr==0) doesn't work with consensus */
				*pup=pleft=*pcurr=-1000;
			else
				end_of_trace=FALSE;


			/* 3/14/05 gary benson -- made >= to extend alignments as far as possible */ 
			if (*pcurr>=maxscore)  /* test for maximum */
			{
				maxscore=*pcurr;
				minrow=realr;
			}
			pcurr--;
			pdiag--;
			pup--;

		}  
		pcurr=&S[r][adjlength];
		pdiag=&Diag[adjmone];
		for (c=adjmone;c>=0;c--)
		{
			*pdiag=*pcurr;
			pcurr--;
			pdiag--;
		}
		Diag[adjlength]=*pcurr;


	}
	/***********************************/
	r=0;
	realr=start-1;

	pup=Up;
	pdiag=&Diag[0];
	pcurr=&S[r][0];

	/* initialize_scoring_array top row (*pcurr)*/
	/* initialize diagonal branch (*pdiag) and up branch (*pup) */
	for (c=0;c<size;c++)
	{
		*pup=(*pdiag=(*pcurr=maxscore))+ Delta;
		pup++;
		pdiag++;
		pcurr++;
	}

	end_of_trace=FALSE;  
	while ((!end_of_trace) && (realr<Length) && (r<maxwraplength))
	{
		r++;
		realr++;
		Rows++;
		currchar=Sequence[realr];

		pcurr=&S[r][0];
		pleft=Delta;        /* first pass. S[r][adjlength]=0  */
		pdiag=&Diag[0];
		pup=&Up[0];
		for (c=0;c<size;c++)
		{
			*pdiag+=match(currchar,EC[c]);
			pleft=(*pcurr=max4(0,*pdiag,*pup,pleft))+Delta;
			pcurr++;
			pdiag++;
			pup++;
		}
		end_of_trace=TRUE;      /* setup for encountering a break in the trace */
		/*second pass */
		pcurr=&S[r][0];
		/* pleft set from first pass */
		pdiag=&Diag[0];
		pup=&Up[0];
		for (c=0;c<size;c++)
		{
			*pup=(pleft=(*pcurr=max4(0,*pdiag,*pup,pleft))+Delta);
			if((realr>=start)&&(*pcurr==0))
				/* if (*pcurr==0) doesn't work with consensus */
				*pup=pleft=*pcurr=-1000;
			else
				end_of_trace=FALSE;


			/* 3/14/05 gary benson -- made >= to extend alignments as far as possible */ 
			if (*pcurr>=maxscore)  /* test for maximum */
			{
				maxscore=*pcurr;
				maxrow=realr;       /* ??? do we need this? */
			}
			pcurr++;
			pdiag++;
			pup++;

		}  
		pcurr=&S[r][0];
		pdiag=&Diag[1];
		for (c=0;c<size;c++)
		{
			*pdiag=*pcurr;
			pcurr++;
			pdiag++;
		}
		Diag[0]=Diag[size];

	}    

	/* store last position scanned with this pattern */

	Wrapend=realr;

	/* redo for an accurate alignment if score is big enough to report */
	Reportmin=0;
	if (maxscore>=Reportmin)
	{

		r=0;
		modstart=start+1;
		realr=minrow-size-1;    /* go back at least one additional pattern */
		if (realr<0) realr=0;   /* length in case of misalignment in */
		/* backwards scan due to starting with */
		/* maxscore in every column */
		maxscore=0;
		/****/
		pup=Up;
		pdiag=&Diag[0];
		pcurr=&S[r][0];

		/* initialize_scoring_array top row (*pcurr)*/
		/* initialize diagonal branch (*pdiag) and up branch (*pup) */
		for (c=0;c<size;c++)
		{
			*pup=(*pdiag=(*pcurr=maxscore))+ Delta;
			pup++;
			pdiag++;
			pcurr++;
		}

		end_of_trace=FALSE;  
		while ((!end_of_trace) && (realr<Length) && (r<maxwraplength))
		{
			r++;
			realr++;
			Rows++;
			currchar=Sequence[realr];

			pcurr=&S[r][0];
			pleft=Delta;      /* first pass. S[r][adjlength]=0  */
			pdiag=&Diag[0];
			pup=&Up[0];
			for (c=0;c<size;c++)
			{
				*pdiag+=match(currchar,EC[c]);
				pleft=(*pcurr=max4(0,*pdiag,*pup,pleft))+Delta;
				pcurr++;
				pdiag++;
				pup++;
			}
			end_of_trace=TRUE;    /* setup for encountering a break in the trace */
			/*second pass */
			pcurr=&S[r][0];
			/* pleft set from first pass */
			pdiag=&Diag[0];
			pup=&Up[0];
			for (c=0;c<size;c++)
			{
				*pup=(pleft=(*pcurr=max4(0,*pdiag,*pup,pleft))+Delta);
				if ((realr>=modstart)
						&&(*pcurr==0))
				{*pup=pleft=*pcurr=-1000;}
				else
					end_of_trace=FALSE;

				if (*pcurr>maxscore)    /* test for maximum */
				{
					maxscore=*pcurr;
					maxrealrow=realr; 
					maxrow=r;
					maxcol=c;
				}
				pcurr++;
				pdiag++;
				pup++;

			}  
			pcurr=&S[r][0];
			pdiag=&Diag[1];
			for (c=0;c<size;c++)
			{
				*pdiag=*pcurr;
				pcurr++;
				pdiag++;
			}
			Diag[0]=Diag[size];

		}          
	}

	if(Criteria_print)trf_message("   DWPlength:%d",realr-(minrow-size+1)+1);

	/* test for report */
	if (maxscore>=Reportmin)
	{
		Maxrealrow=maxrealrow;
		Maxrow=maxrow;
		Maxcol=maxcol;
		Maxscore=maxscore;

	}

}

/*******************************************************************/

/***************************** newwrap()  **************************/

/*******************************************************************/
#ifdef TRF_TEST_FULLSTATE
static void newwrap_fullstate(int start, int size, int consensuspresent)

{
	int g;
	int diag_score;
	const int *similarity_row;
	if (Delta>0)
	{
		newwrap_original(start,size,consensuspresent);
		return;
	}
	register int *pup, *pdiag, *pcurr,pleft;
	int adjlength,adjmone,c,realr,end_of_trace,
	    maxscore,minrow,maxrow,maxcol,modstart,maxrealrow;
	char currchar;
	/* Feb 16, 2016 Yozen */
	unsigned int r;

	/* fill EC */
	if(consensuspresent)
		for(g=0;g<size;g++) EC[g]=Consensus.pattern[g];
	else
		for(g=0;g<size;g++) EC[g]=Sequence[start-size+g+1];


	/* backward wdp */
	maxscore=0;
	realr=start+1;
	r=maxwraplength;

	adjlength=size-1;
	adjmone=adjlength-1;
	pup=Up;
	pdiag=Diag;
	pcurr=&S[r][0];

	for (c=0;c<size;c++)
	{
		*pup=(*pdiag=(*pcurr=maxscore))+ Delta;
		pup++;
		pdiag++;
		pcurr++;
	}
	/* backwards */
	end_of_trace=FALSE;  
	while ((!end_of_trace) && (realr>1) && (r>0))
	{
		r--;
		realr--;
		Rows++;
		currchar=Sequence[realr];
		similarity_row=SM+256*currchar;

		pcurr=&S[r][adjlength];
		pleft=Delta;        /* first pass. set S[r][0]=0 */
		pdiag=&Diag[adjlength];
		pup=&Up[adjlength];
		for (c=adjlength;c>=0;c--)
		{
			diag_score=*pdiag+similarity_row[EC[c]];
			pleft=(*pcurr=max4(0,diag_score,*pup,pleft))+Delta;
			if(c>0) { pcurr--; pdiag--; pup--; }
		}

		end_of_trace=TRUE;      /* setup for encountering a break in the trace */
		/*second pass */
		pcurr=&S[r][adjlength];
		/* pleft set from first pass */
		pdiag=&Diag[adjlength];
		pup=&Up[adjlength];
		for (c=adjlength;c>=0;c--)
		{
			*pup=(pleft=(*pcurr=max2(*pcurr,pleft))+Delta);
			if((realr<=start-max(size,Min_Distance_Window))&&(*pcurr==0)) 
				/* if (*pcurr==0) doesn't work with consensus */
				*pup=pleft=*pcurr=-1000;
			else
				end_of_trace=FALSE;


			/* 3/14/05 gary benson -- made >= to extend alignments as far as possible */ 
			if (*pcurr>=maxscore)  /* test for maximum */
			{
				maxscore=*pcurr;
				minrow=realr;
			}
			if(c>0) Diag[c-1]=*pcurr;
			if(c>0) { pcurr--; pup--; }

		}
		Diag[adjlength]=S[r][0];


	}
	/***********************************/
	r=0;
	realr=start-1;

	pup=Up;
	pdiag=&Diag[0];
	pcurr=&S[r][0];

	/* initialize_scoring_array top row (*pcurr)*/
	/* initialize diagonal branch (*pdiag) and up branch (*pup) */
	for (c=0;c<size;c++)
	{
		*pup=(*pdiag=(*pcurr=maxscore))+ Delta;
		pup++;
		pdiag++;
		pcurr++;
	}

	end_of_trace=FALSE;  
	while ((!end_of_trace) && (realr<Length) && (r<maxwraplength))
	{
		r++;
		realr++;
		Rows++;
		currchar=Sequence[realr];
		similarity_row=SM+256*currchar;

		pcurr=&S[r][0];
		pleft=Delta;        /* first pass. S[r][adjlength]=0  */
		pdiag=&Diag[0];
		pup=&Up[0];
		for (c=0;c<size;c++)
		{
			diag_score=*pdiag+similarity_row[EC[c]];
			pleft=(*pcurr=max4(0,diag_score,*pup,pleft))+Delta;
			pcurr++;
			pdiag++;
			pup++;
		}
		end_of_trace=TRUE;      /* setup for encountering a break in the trace */
		/*second pass */
		pcurr=&S[r][0];
		/* pleft set from first pass */
		pdiag=&Diag[0];
		pup=&Up[0];
		for (c=0;c<size;c++)
		{
			*pup=(pleft=(*pcurr=max2(*pcurr,pleft))+Delta);
			if((realr>=start)&&(*pcurr==0))
				/* if (*pcurr==0) doesn't work with consensus */
				*pup=pleft=*pcurr=-1000;
			else
				end_of_trace=FALSE;


			/* 3/14/05 gary benson -- made >= to extend alignments as far as possible */ 
			if (*pcurr>=maxscore)  /* test for maximum */
			{
				maxscore=*pcurr;
				maxrow=realr;       /* ??? do we need this? */
			}
			Diag[c+1]=*pcurr;
			pcurr++;
			pup++;

		}
		Diag[0]=Diag[size];

	}    

	/* store last position scanned with this pattern */

	Wrapend=realr;

	/* redo for an accurate alignment if score is big enough to report */
	Reportmin=0;
	if (maxscore>=Reportmin)
	{

		r=0;
		modstart=start+1;
		realr=minrow-size-1;    /* go back at least one additional pattern */
		if (realr<0) realr=0;   /* length in case of misalignment in */
		/* backwards scan due to starting with */
		/* maxscore in every column */
		maxscore=0;
		/****/
		pup=Up;
		pdiag=&Diag[0];
		pcurr=&S[r][0];

		/* initialize_scoring_array top row (*pcurr)*/
		/* initialize diagonal branch (*pdiag) and up branch (*pup) */
		for (c=0;c<size;c++)
		{
			*pup=(*pdiag=(*pcurr=maxscore))+ Delta;
			pup++;
			pdiag++;
			pcurr++;
		}

		end_of_trace=FALSE;  
		while ((!end_of_trace) && (realr<Length) && (r<maxwraplength))
		{
			r++;
			realr++;
			Rows++;
			currchar=Sequence[realr];
			similarity_row=SM+256*currchar;

			pcurr=&S[r][0];
			pleft=Delta;      /* first pass. S[r][adjlength]=0  */
			pdiag=&Diag[0];
			pup=&Up[0];
			for (c=0;c<size;c++)
			{
				diag_score=*pdiag+similarity_row[EC[c]];
				pleft=(*pcurr=max4(0,diag_score,*pup,pleft))+Delta;
				pcurr++;
				pdiag++;
				pup++;
			}
			end_of_trace=TRUE;    /* setup for encountering a break in the trace */
			/*second pass */
			pcurr=&S[r][0];
			/* pleft set from first pass */
			pdiag=&Diag[0];
			pup=&Up[0];
			for (c=0;c<size;c++)
			{
				*pup=(pleft=(*pcurr=max2(*pcurr,pleft))+Delta);
				if ((realr>=modstart)
						&&(*pcurr==0))
				{*pup=pleft=*pcurr=-1000;}
				else
					end_of_trace=FALSE;

				if (*pcurr>maxscore)    /* test for maximum */
				{
					maxscore=*pcurr;
					maxrealrow=realr; 
					maxrow=r;
					maxcol=c;
				}
				Diag[c+1]=*pcurr;
				pcurr++;
				pup++;

			}
			Diag[0]=Diag[size];

		}          
	}

	if(Criteria_print)trf_message("   DWPlength:%d",realr-(minrow-size+1)+1);

	/* test for report */
	if (maxscore>=Reportmin)
	{
		Maxrealrow=maxrealrow;
		Maxrow=maxrow;
		Maxcol=maxcol;
		Maxscore=maxscore;

	}

}
#endif

/* AVX2 rows for the fast path below.
 *
 * For nonpositive Delta one circular row is
 *   v1[c] = max(0, diag[c]+sim[c], prev[c]+Delta, v1[c-1]+Delta)   (first pass)
 *   S[c]  = max(v1[c], v1[last]+(position of c after last)*Delta)  (second pass)
 * followed by the zero -> -1000 trace sentinel. The second line is the scalar
 * second pass in closed form: every first-pass predecessor is already
 * dominated by v1[c], so only the value carried around the circle remains,
 * and a sentinel reset can only lower a carried value that was not positive.
 * The first pass is a prefix maximum in ramp coordinates, as in the
 * narrowband kernels (no packed reductions are used here, so there is no
 * score-range restriction). diag is the preceding row rotated by one column and the
 * upper neighbour is the preceding row itself, so Diag/Up are not maintained;
 * a reset upper neighbour contributes -1000+Delta instead of -1000, both
 * negative. All arithmetic is exact 32-bit integer arithmetic under
 * trf_newwrap_avx2_domain, so stored cells, endpoints and maxima equal the
 * scalar ones. Cells are handled in at most three eight-int blocks; lanes at
 * or beyond size are excluded and written back unchanged.
 */
#include <limits.h>
#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__)) && \
    INT_MAX == 2147483647 && !defined(TRF_NEWWRAP_DISABLE_SIMD)
#define TRF_NEWWRAP_HAVE_AVX2 1
#include <immintrin.h>
#include <stdint.h>

static const int32_t trf_newwrap_lanes[16]={-1,-1,-1,-1,-1,-1,-1,-1,0,0,0,0,0,0,0,0};

static int trf_newwrap_avx2_domain(int size)
{
    int64_t a=Alpha,b=-(int64_t)Beta,d=-(int64_t)Delta;
    int64_t top=a*(int64_t)maxwraplength;
#ifdef TRF_NEWWRAP_FORCE_SCALAR
    return 0;
#endif
    if(size<1 || size>20 || a<0 || b<0 || d<0 || maxwraplength<1) return 0;
    /* scores stay within top; ramp coordinates add at most 32*d */
    if(top>INT_MAX/2 || b>(1<<28) || d>(1<<20)) return 0;
    return __builtin_cpu_supports("avx2")!=0;
}

__attribute__((target("avx2")))
static inline __attribute__((always_inline)) int trf_newwrap_hmax(__m256i x)
{
    __m128i m=_mm_max_epi32(_mm256_castsi256_si128(x),_mm256_extracti128_si256(x,1));
    m=_mm_max_epi32(m,_mm_shuffle_epi32(m,_MM_SHUFFLE(1,0,3,2)));
    m=_mm_max_epi32(m,_mm_shuffle_epi32(m,_MM_SHUFFLE(2,3,0,1)));
    return _mm_cvtsi128_si32(m);
}

/* Returns 0 if the call was not handled (nothing observable was changed
 * beyond what the scalar path rewrites); otherwise stores the endpoints of
 * the accurate forward scan and returns 1. *score receives its maximum. */
__attribute__((target("avx2"),noinline))
static int trf_newwrap_avx2(int start,int size,int *score,int *endrealrow,
                            int *endrow,int *endcol)
{
    const __m256i zero=_mm256_setzero_si256(),delta=_mm256_set1_epi32(Delta);
    const __m256i alpha=_mm256_set1_epi32(Alpha),beta=_mm256_set1_epi32(Beta);
    const __m256i iota=_mm256_setr_epi32(0,1,2,3,4,5,6,7);
    const __m256i sevens=_mm256_set1_epi32(7);
    const __m256i slope=_mm256_set1_epi32(-Delta);
    const int blocks=(size+7)>>3;                   /* 1..3 */
    const int lastblock=blocks-1,lastlane=(size-1)&7;
    const __m256i valid=_mm256_loadu_si256((const __m256i *)(trf_newwrap_lanes+7-lastlane));
    const __m256i lastmask=_mm256_cmpeq_epi32(iota,_mm256_set1_epi32(lastlane));
    const __m256i lastindex=_mm256_set1_epi32(lastlane);
    const int backgate=start-max(size,Min_Distance_Window);
    const int length=Length;
    const unsigned int limit=maxwraplength;
    __m256i letters[3],index[3],ramp[3],around[3],b[3],v1[3],v2[3];
    __m256i ch,gate,carry,diag,sim,u,p,v,best,stored;
    unsigned int r,rows=0;
    int realr,end_of_trace,maxscore,minrow=0,maxrealrow=0,maxrow=0,maxcol=0;
    int k,character,canonical,trace_gate,q,dead,better;
    unsigned int where;
    const int *previous;
    int *current;

    for(k=0;k<blocks;k++)
        letters[k]=_mm256_cvtepu8_epi32(_mm_loadl_epi64((const __m128i *)(EC+8*k)));

    /* ---- backward scan: cells descend, position of c is size-1-c ---- */
    for(k=0;k<blocks;k++) {
        index[k]=_mm256_sub_epi32(_mm256_set1_epi32(size-1-8*k),iota);
        if(k==lastblock) index[k]=_mm256_and_si256(index[k],valid);
        ramp[k]=_mm256_mullo_epi32(index[k],slope);
        /* carried value reaches position t after t+1 steps */
        around[k]=_mm256_sub_epi32(zero,_mm256_add_epi32(ramp[k],slope));
    }
    maxscore=0;
    realr=start+1;
    r=limit;
    current=S[r];
    for(k=0;k<blocks;k++) {
        stored=zero;
        if(k==lastblock)
            stored=_mm256_blendv_epi8(_mm256_loadu_si256((const __m256i *)(current+8*k)),stored,valid);
        _mm256_storeu_si256((__m256i *)(current+8*k),stored);
    }
    end_of_trace=FALSE;
    while((!end_of_trace) && (realr>1) && (r>0))
    {
        r--;
        realr--;
        rows++;
        character=(char)Sequence[realr];
        canonical=(character=='A'||character=='C'||character=='G'||character=='T');
        ch=_mm256_set1_epi32(canonical?character:-1);
        trace_gate=(realr<=backgate);
        gate=_mm256_set1_epi32(trace_gate?-1000:0);
        previous=S[r+1];
        current=S[r];
        for(k=0;k<blocks;k++) b[k]=_mm256_loadu_si256((const __m256i *)(previous+8*k));
        carry=zero;
        for(k=lastblock;k>=0;k--) {
            /* diag[c]=previous[c+1], and previous[0] for the last cell */
            diag=_mm256_alignr_epi8(_mm256_permute2x128_si256(b[k],k<lastblock?b[k+1]:zero,0x21),b[k],4);
            if(k==lastblock)
                diag=_mm256_blendv_epi8(diag,_mm256_broadcastd_epi32(_mm256_castsi256_si128(b[0])),lastmask);
            sim=_mm256_blendv_epi8(beta,alpha,_mm256_cmpeq_epi32(letters[k],ch));
            u=_mm256_max_epi32(zero,_mm256_max_epi32(_mm256_add_epi32(diag,sim),_mm256_add_epi32(b[k],delta)));
            if(k==lastblock) u=_mm256_and_si256(u,valid);
            u=_mm256_add_epi32(u,ramp[k]);
            u=_mm256_max_epi32(u,_mm256_srli_si256(u,4));
            u=_mm256_max_epi32(u,_mm256_srli_si256(u,8));
            u=_mm256_max_epi32(u,_mm256_shuffle_epi32(_mm256_permute2x128_si256(u,u,0x81),_MM_SHUFFLE(0,0,0,0)));
            p=_mm256_max_epi32(u,carry);
            carry=_mm256_broadcastd_epi32(_mm256_castsi256_si128(p));
            v1[k]=_mm256_sub_epi32(p,ramp[k]);
        }
        /* value carried around the circle: first-pass score of cell 0 */
        carry=_mm256_broadcastd_epi32(_mm256_castsi256_si128(v1[0]));
        best=zero;
        for(k=0;k<blocks;k++) {
            v=_mm256_max_epi32(v1[k],_mm256_add_epi32(carry,around[k]));
            stored=_mm256_blendv_epi8(v,gate,_mm256_cmpeq_epi32(v,zero));
            if(k==lastblock) {
                v=_mm256_and_si256(v,valid);
                stored=_mm256_blendv_epi8(_mm256_loadu_si256((const __m256i *)(current+8*k)),stored,valid);
            }
            best=_mm256_max_epi32(best,v);
            _mm256_storeu_si256((__m256i *)(current+8*k),stored);
        }
        q=trf_newwrap_hmax(best);
        dead=trace_gate & (q==0);
        end_of_trace=dead;
        /* >= : extend alignments as far as possible (cells of a dead row are -1000) */
        better=(!dead) & (q>=maxscore);
        maxscore=better?q:maxscore;
        minrow=better?realr:minrow;
    }
    if(rows==0) return 0;       /* minrow would be undefined */

    /* ---- accurate forward scan: cells ascend, position of c is c ---- */
    for(k=0;k<blocks;k++) {
        index[k]=_mm256_add_epi32(_mm256_set1_epi32(8*k),iota);
        ramp[k]=_mm256_mullo_epi32(index[k],slope);
        around[k]=_mm256_sub_epi32(zero,_mm256_add_epi32(ramp[k],slope));
    }
    r=0;
    realr=minrow-size-1;    /* go back at least one additional pattern */
    if(realr<0) realr=0;
    maxscore=0;
    current=S[0];
    for(k=0;k<blocks;k++) {
        stored=zero;
        if(k==lastblock)
            stored=_mm256_blendv_epi8(_mm256_loadu_si256((const __m256i *)(current+8*k)),stored,valid);
        _mm256_storeu_si256((__m256i *)(current+8*k),stored);
    }
    end_of_trace=FALSE;
    while((!end_of_trace) && (realr<length) && (r<limit))
    {
        r++;
        realr++;
        rows++;
        character=(char)Sequence[realr];
        canonical=(character=='A'||character=='C'||character=='G'||character=='T');
        ch=_mm256_set1_epi32(canonical?character:-1);
        trace_gate=(realr>=start+1);
        gate=_mm256_set1_epi32(trace_gate?-1000:0);
        previous=S[r-1];
        current=S[r];
        for(k=0;k<blocks;k++) b[k]=_mm256_loadu_si256((const __m256i *)(previous+8*k));
        carry=zero;
        for(k=0;k<blocks;k++) {
            /* diag[c]=previous[c-1], and previous[size-1] for cell 0 */
            if(k==0)
                diag=_mm256_permutevar8x32_epi32(b[lastblock],lastindex);
            else
                diag=b[k-1];
            diag=_mm256_alignr_epi8(b[k],_mm256_permute2x128_si256(diag,b[k],0x21),12);
            sim=_mm256_blendv_epi8(beta,alpha,_mm256_cmpeq_epi32(letters[k],ch));
            u=_mm256_max_epi32(zero,_mm256_max_epi32(_mm256_add_epi32(diag,sim),_mm256_add_epi32(b[k],delta)));
            u=_mm256_add_epi32(u,ramp[k]);
            u=_mm256_max_epi32(u,_mm256_slli_si256(u,4));
            u=_mm256_max_epi32(u,_mm256_slli_si256(u,8));
            u=_mm256_max_epi32(u,_mm256_shuffle_epi32(_mm256_permute2x128_si256(u,u,0x08),_MM_SHUFFLE(3,3,3,3)));
            p=_mm256_max_epi32(u,carry);
            carry=_mm256_permutevar8x32_epi32(p,sevens);
            v1[k]=_mm256_sub_epi32(p,ramp[k]);
        }
        /* value carried around the circle: first-pass score of the last cell */
        carry=_mm256_permutevar8x32_epi32(v1[lastblock],lastindex);
        best=zero;
        for(k=0;k<blocks;k++) {
            v=_mm256_max_epi32(v1[k],_mm256_add_epi32(carry,around[k]));
            stored=_mm256_blendv_epi8(v,gate,_mm256_cmpeq_epi32(v,zero));
            if(k==lastblock) {
                v=_mm256_and_si256(v,valid);
                stored=_mm256_blendv_epi8(_mm256_loadu_si256((const __m256i *)(current+8*k)),stored,valid);
            }
            v2[k]=v;
            best=_mm256_max_epi32(best,v);
            _mm256_storeu_si256((__m256i *)(current+8*k),stored);
        }
        q=trf_newwrap_hmax(best);
        /* first column holding the row maximum (cells are nonnegative, and
         * lanes at or beyond size are excluded below) */
        best=_mm256_set1_epi32(q);
        where=0;
        for(k=0;k<blocks;k++)
            where|=(unsigned int)_mm256_movemask_ps(_mm256_castsi256_ps(
                       _mm256_cmpeq_epi32(v2[k],best)))<<(8*k);
        where&=(size==32)?~0u:((1u<<size)-1u);
        dead=trace_gate & (q==0);
        end_of_trace=dead;
        better=(q>maxscore);
        maxscore=better?q:maxscore;
        maxrealrow=better?realr:maxrealrow;
        maxrow=better?(int)r:maxrow;
        maxcol=better?__builtin_ctz(where):maxcol;
    }
    /* Leave the scratch arrays as the scalar passes do after their last row
     * (r) of the accurate scan: a reset cell's upper score is the sentinel. */
    current=S[r];
    for(k=0;k<size;k++) Up[k]=(current[k]==-1000)?-1000:current[k]+Delta;
    if(r==0) for(k=0;k<size;k++) Diag[k]=0;
    else {
        for(k=0;k<size;k++) Diag[k+1]=current[k];
        Diag[0]=Diag[size];
    }
    Rows+=rows;
    *score=maxscore;
    *endrealrow=maxrealrow;
    *endrow=maxrow;
    *endcol=maxcol;
    return 1;
}
#else
#define TRF_NEWWRAP_HAVE_AVX2 0
#endif

/* The preliminary forward scan never changes minrow. Its matrix is
 * overwritten by the accurate forward scan, which is always entered because
 * Reportmin is reset to zero. For a positive accurate maximum, endpoints and
 * all matrix cells used by traceback are therefore independent of that scan.
 * Rows, Wrapend and trailing scratch rows have no observable consumer on the
 * newwrap path. Retain the full-state implementation above for differential
 * testing, and the original path for unsupported or degenerate domains. */
void newwrap(int start, int size, int consensuspresent)

{
	int g;
	double entry_rows=Rows;
	int diag_score;
	const int *similarity_row;
	if (Delta>0 || size<1 || size>20 || maxwraplength==0 || Criteria_print)
	{
		newwrap_original(start,size,consensuspresent);
		return;
	}
	register int *pup, *pdiag, *pcurr,pleft;
	int adjlength,adjmone,c,realr,end_of_trace,
	    maxscore,minrow,maxrow,maxcol,modstart,maxrealrow;
	char currchar;
	/* Feb 16, 2016 Yozen */
	unsigned int r;

	/* fill EC */
	if(consensuspresent)
		for(g=0;g<size;g++) EC[g]=Consensus.pattern[g];
	else
		for(g=0;g<size;g++) EC[g]=Sequence[start-size+g+1];

#if TRF_NEWWRAP_HAVE_AVX2
	if(trf_newwrap_avx2_domain(size) &&
	   trf_newwrap_avx2(start,size,&maxscore,&maxrealrow,&maxrow,&maxcol))
	{
		Reportmin=0;
		/* A positive accurate maximum defines all reported endpoints. Preserve
		 * the original path for degenerate scoring where they may be undefined. */
		if(maxscore<=0)
		{
			Rows=entry_rows;
			newwrap_original(start,size,consensuspresent);
			return;
		}
		Maxrealrow=maxrealrow;
		Maxrow=maxrow;
		Maxcol=maxcol;
		Maxscore=maxscore;
		return;
	}
#endif

	/* backward wdp */
	maxscore=0;
	realr=start+1;
	r=maxwraplength;

	adjlength=size-1;
	adjmone=adjlength-1;
	pup=Up;
	pdiag=Diag;
	pcurr=&S[r][0];

	for (c=0;c<size;c++)
	{
		*pup=(*pdiag=(*pcurr=maxscore))+ Delta;
		pup++;
		pdiag++;
		pcurr++;
	}
	/* backwards */
	end_of_trace=FALSE;  
	while ((!end_of_trace) && (realr>1) && (r>0))
	{
		r--;
		realr--;
		Rows++;
		currchar=Sequence[realr];
		similarity_row=SM+256*currchar;

		pcurr=&S[r][adjlength];
		pleft=Delta;        /* first pass. set S[r][0]=0 */
		pdiag=&Diag[adjlength];
		pup=&Up[adjlength];
		for (c=adjlength;c>=0;c--)
		{
			diag_score=*pdiag+similarity_row[EC[c]];
			/* Only the numeric score survives; pairwise max avoids redundant comparisons. */
			diag_score=max2(diag_score,*pup);
			diag_score=max2(diag_score,pleft);
			diag_score=max2(0,diag_score);
			pleft=(*pcurr=diag_score)+Delta;
			if(c>0) { pcurr--; pdiag--; pup--; }
		}

		end_of_trace=TRUE;      /* setup for encountering a break in the trace */
		/*second pass */
		pcurr=&S[r][adjlength];
		/* pleft set from first pass */
		pdiag=&Diag[adjlength];
		pup=&Up[adjlength];
		for (c=adjlength;c>=0;c--)
		{
			*pup=(pleft=(*pcurr=max2(*pcurr,pleft))+Delta);
			if((realr<=start-max(size,Min_Distance_Window))&&(*pcurr==0)) 
				/* if (*pcurr==0) doesn't work with consensus */
				*pup=pleft=*pcurr=-1000;
			else
				end_of_trace=FALSE;


			/* 3/14/05 gary benson -- made >= to extend alignments as far as possible */ 
			if (*pcurr>=maxscore)  /* test for maximum */
			{
				maxscore=*pcurr;
				minrow=realr;
			}
			if(c>0) Diag[c-1]=*pcurr;
			if(c>0) { pcurr--; pup--; }

		}
		Diag[adjlength]=S[r][0];


	}
	/***********************************/
	/* The accurate scan below reinitializes its full live matrix region. */

	/* redo for an accurate alignment if score is big enough to report */
	Reportmin=0;
	if (maxscore>=Reportmin)
	{

		r=0;
		modstart=start+1;
		realr=minrow-size-1;    /* go back at least one additional pattern */
		if (realr<0) realr=0;   /* length in case of misalignment in */
		/* backwards scan due to starting with */
		/* maxscore in every column */
		maxscore=0;
		/****/
		pup=Up;
		pdiag=&Diag[0];
		pcurr=&S[r][0];

		/* initialize_scoring_array top row (*pcurr)*/
		/* initialize diagonal branch (*pdiag) and up branch (*pup) */
		for (c=0;c<size;c++)
		{
			*pup=(*pdiag=(*pcurr=maxscore))+ Delta;
			pup++;
			pdiag++;
			pcurr++;
		}

		end_of_trace=FALSE;  
		while ((!end_of_trace) && (realr<Length) && (r<maxwraplength))
		{
			r++;
			realr++;
			Rows++;
			currchar=Sequence[realr];
			similarity_row=SM+256*currchar;

			pcurr=&S[r][0];
			pleft=Delta;      /* first pass. S[r][adjlength]=0  */
			pdiag=&Diag[0];
			pup=&Up[0];
			for (c=0;c<size;c++)
			{
				diag_score=*pdiag+similarity_row[EC[c]];
				diag_score=max2(diag_score,*pup);
				diag_score=max2(diag_score,pleft);
				diag_score=max2(0,diag_score);
				pleft=(*pcurr=diag_score)+Delta;
				pcurr++;
				pdiag++;
				pup++;
			}
			end_of_trace=TRUE;    /* setup for encountering a break in the trace */
			/*second pass */
			pcurr=&S[r][0];
			/* pleft set from first pass */
			pdiag=&Diag[0];
			pup=&Up[0];
			for (c=0;c<size;c++)
			{
				*pup=(pleft=(*pcurr=max2(*pcurr,pleft))+Delta);
				if ((realr>=modstart)
						&&(*pcurr==0))
				{*pup=pleft=*pcurr=-1000;}
				else
					end_of_trace=FALSE;

				if (*pcurr>maxscore)    /* test for maximum */
				{
					maxscore=*pcurr;
					maxrealrow=realr; 
					maxrow=r;
					maxcol=c;
				}
				Diag[c+1]=*pcurr;
				pcurr++;
				pup++;

			}
			Diag[0]=Diag[size];

		}          
	}

	/* A positive accurate maximum defines all reported endpoints. Preserve
	 * the original path for degenerate scoring where they may be undefined. */
	if(maxscore<=0)
	{
		Rows=entry_rows;
		newwrap_original(start,size,consensuspresent);
		return;
	}

	if(Criteria_print)trf_message("   DWPlength:%d",realr-(minrow-size+1)+1);

	/* test for report */
	if (maxscore>=Reportmin)
	{
		Maxrealrow=maxrealrow;
		Maxrow=maxrow;
		Maxcol=maxcol;
		Maxscore=maxscore;

	}

}
