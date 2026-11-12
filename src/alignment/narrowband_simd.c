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

/* Exact SSE4.1 prefix scan for the ordinary one-column band shifts.
 * The scalar kernel remains the baseline ISA and unusual-domain fallback.
 * SM follows init_sm's production invariant: canonical self-pairs are Alpha,
 * every other byte pair is Beta. No matrix writer exists outside init_sm.
 */
#include <stdint.h>
#include <limits.h>
#include <string.h>
#include "core/globals.h"
#include "core/report.h"
#include "alignment/narrowband_scalar.h"
#include "alignment/narrowband_simd.h"

#if TRF_NARROWBAND_HAVE_SIMD
#include <smmintrin.h>

struct trf_narrowband_row { int score,last,matchcol,dead; };

/* Hide a value's provenance / forbid if-conversion of a rare branch. */
#define TRF_NARROWBAND_OPAQUE(x) __asm__ __volatile__("" : "+r"(x))
#define TRF_NARROWBAND_BARRIER() __asm__ __volatile__("")

static inline int trf_narrowband_column(int c,int offset,int direction,int size)
{
    int value=c+direction*offset;
    if(value<0) value+=size;
    else if(value>=size) value-=size;
    return value;
}

int trf_narrowband_vector_domain(int start,int size,int wb,int wf)
{
    int64_t a=Alpha,b=-(int64_t)Beta,d=-(int64_t)Delta;
    if(size<3 || size>4000 || wb<1 || wf<1 || wb>(MAXBANDWIDTH-1)/2 ||
       wf>(MAXBANDWIDTH-1)/2 || 2*wb>=size || 2*wf>=size || start<size ||
       maxwraplength<1 || a<0 || b<0 || d<0) return 0;
    int64_t low=((int64_t)wb+wf)*d;
    if(low<1000) low=1000;
    if(
       a*maxwraplength>INT_MAX || 8*d>INT_MAX ||
       low+(b>d?b:d)>-(int64_t)INT_MIN) return 0;
#ifdef TRF_NARROWBAND_FORCE_SCALAR
    return 0;
#else
    return __builtin_cpu_supports("sse4.1")!=0;
#endif
}

/* Each block is in traversal order, including reversed backward loads/stores.
 * With nonpositive Delta, deferred zero -> -1000 stores cannot affect the
 * following raw maximum: both outgoing values are dominated by its zero.
 * Track FIRST and LAST row maxima separately, matching the original > / >=.
 */
__attribute__((target("sse4.1"),noinline))
static void trf_narrowband_scan_unpacked(int *current,const int *previous,
    const unsigned char *pattern,int size,int c,int w,int direction,
    int trace_gate,int character,struct trf_narrowband_row *out)
{
    __m128i zero=_mm_setzero_si128(),delta=_mm_set1_epi32(Delta);
    __m128i delta2=_mm_set1_epi32(2*Delta), alpha=_mm_set1_epi32(Alpha);
    __m128i beta=_mm_set1_epi32(Beta), ch=_mm_set1_epi32(character);
    __m128i lane_best=_mm_set1_epi32(-1), first=zero,last=zero;
    __m128i index=_mm_setr_epi32(0,1,2,3),step=_mm_set1_epi32(4);
    int t=0,carry=0,canonical=(character=='A'||character=='C'||character=='G'||character=='T');
    for(;t+4<=2*w;t+=4) {
        int c0=c,c1=trf_narrowband_column(c,1,direction,size);
        int c2=trf_narrowband_column(c,2,direction,size);
        int c3=trf_narrowband_column(c,3,direction,size);
        __m128i letters=_mm_setr_epi32(pattern[c0],pattern[c1],pattern[c2],pattern[c3]);
        __m128i sim=canonical?_mm_blendv_epi8(beta,alpha,_mm_cmpeq_epi32(letters,ch)):beta;
        __m128i diag,up;
        if(direction>0) {
            diag=_mm_loadu_si128((const __m128i *)(previous+t));
            up=_mm_loadu_si128((const __m128i *)(previous+t+1));
        } else {
            diag=_mm_shuffle_epi32(_mm_loadu_si128((const __m128i *)(previous+2*w-t-3)),_MM_SHUFFLE(0,1,2,3));
            up=_mm_shuffle_epi32(_mm_loadu_si128((const __m128i *)(previous+2*w-t-4)),_MM_SHUFFLE(0,1,2,3));
        }
        __m128i v=_mm_max_epi32(zero,_mm_max_epi32(_mm_add_epi32(diag,sim),_mm_add_epi32(up,delta)));
        v=_mm_insert_epi32(v,max2(_mm_extract_epi32(v,0),carry+Delta),0);
        v=_mm_max_epi32(v,_mm_add_epi32(_mm_slli_si128(v,4),delta));
        v=_mm_max_epi32(v,_mm_add_epi32(_mm_slli_si128(v,8),delta2));
        carry=_mm_extract_epi32(v,3);
        __m128i greater=_mm_cmpgt_epi32(v,lane_best);
        __m128i ge=_mm_or_si128(greater,_mm_cmpeq_epi32(v,lane_best));
        first=_mm_blendv_epi8(first,index,greater);
        last=_mm_blendv_epi8(last,index,ge);
        lane_best=_mm_max_epi32(lane_best,v);
        __m128i stored=trace_gate?_mm_blendv_epi8(v,_mm_set1_epi32(-1000),_mm_cmpeq_epi32(v,zero)):v;
        if(direction>0) _mm_storeu_si128((__m128i *)(current+t),stored);
        else _mm_storeu_si128((__m128i *)(current+2*w-t-3),_mm_shuffle_epi32(stored,_MM_SHUFFLE(0,1,2,3)));
        index=_mm_add_epi32(index,step);
        c=trf_narrowband_column(c,4,direction,size);
    }
    int best[4],fi[4],la[4],q=-1,first_at=0,last_at=0;
    _mm_storeu_si128((__m128i *)best,lane_best);
    _mm_storeu_si128((__m128i *)fi,first);
    _mm_storeu_si128((__m128i *)la,last);
    for(int lane=0;lane<4;++lane) {
        if(best[lane]>q) { q=best[lane]; first_at=fi[lane]; last_at=la[lane]; }
        else if(best[lane]==q) {
            if(fi[lane]<first_at) first_at=fi[lane];
            if(la[lane]>last_at) last_at=la[lane];
        }
    }
    for(;t<=2*w;++t) {
        int i=direction>0?t:2*w-t;
        int sim=(canonical && pattern[c]==character)?Alpha:Beta;
        int raw=max3(0,previous[i]+sim,carry+Delta);
        if(t<2*w) raw=max2(raw,previous[i+direction]+Delta);
        current[i]=(trace_gate && raw==0)?-1000:raw;
        carry=raw;
        if(raw>q) { q=raw; first_at=t; last_at=t; }
        else if(raw==q) last_at=t;
        c=trf_narrowband_column(c,1,direction,size);
    }
    out->dead=trace_gate && q==0;
    out->score=out->dead?-1000:q;
    out->last=last_at;
    int first_col=trf_narrowband_column(c,first_at-(2*w+1),direction,size);
    /* c has advanced one full row; recompute the winning diagonal only. */
    int first_i=direction>0?first_at:2*w-first_at;
    int sim=(canonical && pattern[first_col]==character)?Alpha:Beta;
    out->matchcol=(q==previous[first_i]+sim && sim==Alpha)?first_col:-2;
}

__attribute__((target("sse4.1")))
static inline __attribute__((always_inline)) void trf_narrowband_scan(int *current,const int *previous,
    const unsigned char *pattern,int size,int c,int w,int direction,
    int trace_gate,int character,struct trf_narrowband_row *out)
{
    if((int64_t)Alpha*maxwraplength>(INT_MAX-255)/256) {
        trf_narrowband_scan_unpacked(current,previous,pattern,size,c,w,direction,
                                    trace_gate,character,out);
        return;
    }
    __m128i zero=_mm_setzero_si128(),delta=_mm_set1_epi32(Delta);
    __m128i delta2=_mm_set1_epi32(2*Delta), alpha=_mm_set1_epi32(Alpha);
    __m128i beta=_mm_set1_epi32(Beta), ch=_mm_set1_epi32(character);
    __m128i best_first=_mm_set1_epi32(INT_MIN),best_last=best_first;
    __m128i index=_mm_setr_epi32(0,1,2,3),step=_mm_set1_epi32(4);
    int t=0,carry=0,canonical=(character=='A'||character=='C'||character=='G'||character=='T');
    for(;t+4<=2*w;t+=4) {
        __m128i letters;
        if(direction>0 && c+3<size) {
            uint32_t bytes;
            memcpy(&bytes,pattern+c,sizeof bytes);
            letters=_mm_cvtepu8_epi32(_mm_cvtsi32_si128(bytes));
        } else if(direction<0 && c>=3) {
            uint32_t bytes;
            memcpy(&bytes,pattern+c-3,sizeof bytes);
            letters=_mm_shuffle_epi32(_mm_cvtepu8_epi32(_mm_cvtsi32_si128(bytes)),_MM_SHUFFLE(0,1,2,3));
        } else {
            int c1=trf_narrowband_column(c,1,direction,size);
            int c2=trf_narrowband_column(c,2,direction,size);
            int c3=trf_narrowband_column(c,3,direction,size);
            letters=_mm_setr_epi32(pattern[c],pattern[c1],pattern[c2],pattern[c3]);
        }
        __m128i sim=canonical?_mm_blendv_epi8(beta,alpha,_mm_cmpeq_epi32(letters,ch)):beta;
        __m128i diag,up;
        if(direction>0) {
            diag=_mm_loadu_si128((const __m128i *)(previous+t));
            up=_mm_loadu_si128((const __m128i *)(previous+t+1));
        } else {
            diag=_mm_shuffle_epi32(_mm_loadu_si128((const __m128i *)(previous+2*w-t-3)),_MM_SHUFFLE(0,1,2,3));
            up=_mm_shuffle_epi32(_mm_loadu_si128((const __m128i *)(previous+2*w-t-4)),_MM_SHUFFLE(0,1,2,3));
        }
        __m128i v=_mm_max_epi32(zero,_mm_max_epi32(_mm_add_epi32(diag,sim),_mm_add_epi32(up,delta)));
        v=_mm_insert_epi32(v,max2(_mm_extract_epi32(v,0),carry+Delta),0);
        v=_mm_max_epi32(v,_mm_add_epi32(_mm_slli_si128(v,4),delta));
        v=_mm_max_epi32(v,_mm_add_epi32(_mm_slli_si128(v,8),delta2));
        carry=_mm_extract_epi32(v,3);
        /* Under the dispatch bound the low eight bits encode exact ties:
         * FIRST maximizes 255-index and LAST maximizes index. */
        __m128i packed=_mm_slli_epi32(v,8);
        best_first=_mm_max_epi32(best_first,_mm_add_epi32(packed,_mm_sub_epi32(_mm_set1_epi32(255),index)));
        best_last=_mm_max_epi32(best_last,_mm_add_epi32(packed,index));
        __m128i stored=trace_gate?_mm_blendv_epi8(v,_mm_set1_epi32(-1000),_mm_cmpeq_epi32(v,zero)):v;
        if(direction>0) _mm_storeu_si128((__m128i *)(current+t),stored);
        else _mm_storeu_si128((__m128i *)(current+2*w-t-3),_mm_shuffle_epi32(stored,_MM_SHUFFLE(0,1,2,3)));
        index=_mm_add_epi32(index,step);
        c=trf_narrowband_column(c,4,direction,size);
    }
    best_first=_mm_max_epi32(best_first,_mm_shuffle_epi32(best_first,_MM_SHUFFLE(1,0,3,2)));
    best_first=_mm_max_epi32(best_first,_mm_shuffle_epi32(best_first,_MM_SHUFFLE(2,3,0,1)));
    best_last=_mm_max_epi32(best_last,_mm_shuffle_epi32(best_last,_MM_SHUFFLE(1,0,3,2)));
    best_last=_mm_max_epi32(best_last,_mm_shuffle_epi32(best_last,_MM_SHUFFLE(2,3,0,1)));
    int first_packed=_mm_cvtsi128_si32(best_first),last_packed=_mm_cvtsi128_si32(best_last);
    int q=first_packed==INT_MIN?-1:first_packed>>8;
    int first_at=255-(first_packed&255),last_at=last_packed&255;
    for(;t<=2*w;++t) {
        int i=direction>0?t:2*w-t;
        int sim=(canonical && pattern[c]==character)?Alpha:Beta;
        int raw=max3(0,previous[i]+sim,carry+Delta);
        if(t<2*w) raw=max2(raw,previous[i+direction]+Delta);
        current[i]=(trace_gate && raw==0)?-1000:raw;
        carry=raw;
        if(raw>q) { q=raw; first_at=t; last_at=t; }
        else if(raw==q) last_at=t;
        c=trf_narrowband_column(c,1,direction,size);
    }
    out->dead=trace_gate && q==0;
    out->score=out->dead?-1000:q;
    out->last=last_at;
    int first_col=trf_narrowband_column(c,first_at-(2*w+1),direction,size);
    /* c has advanced one full row; recompute the winning diagonal only. */
    int first_i=direction>0?first_at:2*w-first_at;
    int sim=(canonical && pattern[first_col]==character)?Alpha:Beta;
    out->matchcol=(q==previous[first_i]+sim && sim==Alpha)?first_col:-2;
}

/* Row kernels for the same one-column shifts.
 *
 * The left recurrence raw[t]=max(base[t],raw[t-1]+Delta) is evaluated in ramp
 * coordinates: with D=-Delta>=0 and u[t]=base[t]+t*D, raw[t]+t*D is the prefix
 * maximum of u. All quantities are exact 32-bit integers under the dispatch
 * bounds (trf_narrowband_ramp_domain, trf_narrowband_packed_limit), so every stored cell equals the scalar
 * value. The only loop-carried vector dependency is one broadcast and one
 * maximum per block, instead of an extract/insert round trip.
 *
 * Cells live on a fixed grid of four-int blocks starting at cell 0. The upper
 * neighbour is assembled from the adjacent block (palignr), so every load of
 * the preceding row coincides with one of its stores. The band has 2w+1
 * cells, hence the highest block holds one or three band cells. Its other
 * lanes lie inside the row (2w+3<=MAXBANDWIDTH is part of the dispatch), are
 * excluded from the recurrence and the reductions, and are written back
 * unchanged. There is no scalar remainder loop.
 * Row reductions use the same packed FIRST/LAST encoding as above.
 */
/* The packed reductions need every cell of a row to be at most
 * TRF_NARROWBAND_PACKED_MAX. A cell exceeds the largest cell of the preceding
 * row by at most Alpha (Alpha>=0, Delta<=0), and the running maxscore of
 * either scan is at least every cell written so far in that scan (the
 * initial rows are nonpositive, maxscore starts at 0). A row may therefore
 * use these kernels while maxscore<=trf_narrowband_packed_limit(); the row
 * loops test that per row and otherwise use the unencoded scan. This keeps
 * the kernels available for large maximum TR lengths (-l), where
 * Alpha*maxwraplength itself exceeds the packed range. */
#define TRF_NARROWBAND_PACKED_MAX ((INT_MAX-255)/256)
static int trf_narrowband_packed_limit(void)
{
    int64_t limit=(int64_t)TRF_NARROWBAND_PACKED_MAX-(int64_t)Alpha;
    return limit<0?-1:(int)limit;
}

static int trf_narrowband_ramp_domain(int wb,int wf)
{
    int64_t a=Alpha;
    int64_t d=-(int64_t)Delta;
#ifdef TRF_NARROWBAND_DISABLE_RAMP
    (void)wb; (void)wf; (void)a; (void)d;
    return 0;
#else
    /* ramp coordinates add at most (MAXBANDWIDTH+8)*d to a packed-range cell */
    return a>=0 && a<=TRF_NARROWBAND_PACKED_MAX && d>=0 &&
           (int64_t)TRF_NARROWBAND_PACKED_MAX+(int64_t)(MAXBANDWIDTH+8)*d<=INT_MAX &&
           2*wb+3<=MAXBANDWIDTH && 2*wf+3<=MAXBANDWIDTH;
#endif
}

/* lanes below n are selected by loading at trf_narrowband_lanes+4-n */
static const int32_t trf_narrowband_lanes[8]={-1,-1,-1,-1,0,0,0,0};

#define TRF_NARROWBAND_REDUCE(x) do { \
    x=_mm_max_epi32(x,_mm_shuffle_epi32(x,_MM_SHUFFLE(1,0,3,2))); \
    x=_mm_max_epi32(x,_mm_shuffle_epi32(x,_MM_SHUFFLE(2,3,0,1))); } while(0)

__attribute__((target("sse4.1")))
static inline __attribute__((always_inline)) void trf_narrowband_scan_forward(
    int *current,const int *previous,const unsigned char *pattern,int size,
    int c,int w,int trace_gate,int character,struct trf_narrowband_row *out)
{
    const int direction=1;
    __m128i zero=_mm_setzero_si128(),delta=_mm_set1_epi32(Delta);
    __m128i alpha=_mm_set1_epi32(Alpha),beta=_mm_set1_epi32(Beta);
    __m128i ch=_mm_set1_epi32(character);
    __m128i ramp=_mm_setr_epi32(0,-Delta,-2*Delta,-3*Delta);
    __m128i rampstep=_mm_set1_epi32(-4*Delta);
    __m128i gate=_mm_set1_epi32(trace_gate?-1000:0);
    __m128i best_first=_mm_set1_epi32(INT_MIN),best_last=best_first;
    __m128i index=_mm_setr_epi32(0,1,2,3),step=_mm_set1_epi32(4);
    __m128i rindex=_mm_setr_epi32(255,254,253,252);
    __m128i carry=zero,block,next,letters,sim,up,u,p,v,packed,stored,valid;
    const int full=(2*w)&~3;        /* cells in complete blocks */
    const int n=2*w+1-full;         /* band cells in the highest block: 1 or 3 */
    const int c0=c;
    int t,canonical=(character=='A'||character=='C'||character=='G'||character=='T');
    if(!canonical) alpha=beta;
    block=_mm_loadu_si128((const __m128i *)previous);
    for(t=0;t<full;t+=4) {
        if(c+3<size) {
            uint32_t bytes;
            memcpy(&bytes,pattern+c,sizeof bytes);
            letters=_mm_cvtepu8_epi32(_mm_cvtsi32_si128((int)bytes));
        } else {
            int c1=trf_narrowband_column(c,1,direction,size);
            int c2=trf_narrowband_column(c,2,direction,size);
            int c3=trf_narrowband_column(c,3,direction,size);
            letters=_mm_setr_epi32(pattern[c],pattern[c1],pattern[c2],pattern[c3]);
        }
        sim=_mm_blendv_epi8(beta,alpha,_mm_cmpeq_epi32(letters,ch));
        next=_mm_loadu_si128((const __m128i *)(previous+t+4));
        up=_mm_alignr_epi8(next,block,4);
        u=_mm_max_epi32(zero,_mm_max_epi32(_mm_add_epi32(block,sim),_mm_add_epi32(up,delta)));
        u=_mm_add_epi32(u,ramp);
        u=_mm_max_epi32(u,_mm_slli_si128(u,4));
        u=_mm_max_epi32(u,_mm_slli_si128(u,8));
        p=_mm_max_epi32(u,carry);
        carry=_mm_shuffle_epi32(p,_MM_SHUFFLE(3,3,3,3));
        v=_mm_sub_epi32(p,ramp);
        packed=_mm_slli_epi32(v,8);
        best_first=_mm_max_epi32(best_first,_mm_add_epi32(packed,rindex));
        best_last=_mm_max_epi32(best_last,_mm_add_epi32(packed,index));
        stored=_mm_add_epi32(v,_mm_and_si128(gate,_mm_cmpeq_epi32(v,zero)));
        _mm_storeu_si128((__m128i *)(current+t),stored);
        index=_mm_add_epi32(index,step);
        rindex=_mm_sub_epi32(rindex,step);
        ramp=_mm_add_epi32(ramp,rampstep);
        block=next;
        c=trf_narrowband_column(c,4,direction,size);
    }
    /* highest block: n band cells, the last of which has no upper neighbour */
    {
        int c1=trf_narrowband_column(c,1,direction,size);
        int c2=trf_narrowband_column(c,2,direction,size);
        letters=_mm_setr_epi32(pattern[c],pattern[c1],pattern[c2],0);
        valid=_mm_loadu_si128((const __m128i *)(trf_narrowband_lanes+4-n));
        sim=_mm_blendv_epi8(beta,alpha,_mm_cmpeq_epi32(letters,ch));
        up=_mm_add_epi32(_mm_srli_si128(block,4),delta);
        up=_mm_and_si128(up,_mm_loadu_si128((const __m128i *)(trf_narrowband_lanes+5-n)));
        u=_mm_max_epi32(zero,_mm_max_epi32(_mm_add_epi32(block,sim),up));
        u=_mm_add_epi32(u,ramp);
        u=_mm_max_epi32(u,_mm_slli_si128(u,4));
        u=_mm_max_epi32(u,_mm_slli_si128(u,8));
        p=_mm_max_epi32(u,carry);
        v=_mm_sub_epi32(p,ramp);
        packed=_mm_blendv_epi8(_mm_set1_epi32(INT_MIN),_mm_slli_epi32(v,8),valid);
        best_first=_mm_max_epi32(best_first,_mm_add_epi32(packed,rindex));
        best_last=_mm_max_epi32(best_last,_mm_add_epi32(packed,index));
        stored=_mm_add_epi32(v,_mm_and_si128(gate,_mm_cmpeq_epi32(v,zero)));
        stored=_mm_blendv_epi8(_mm_loadu_si128((const __m128i *)(current+t)),stored,valid);
        _mm_storeu_si128((__m128i *)(current+t),stored);
    }
    TRF_NARROWBAND_REDUCE(best_first);
    TRF_NARROWBAND_REDUCE(best_last);
    {
        int first_packed=_mm_cvtsi128_si32(best_first),last_packed=_mm_cvtsi128_si32(best_last);
        int q=first_packed>>8;
        int first_at=255-(first_packed&255),last_at=last_packed&255;
        int first_col=trf_narrowband_column(c0,first_at,direction,size);
        int s1=(canonical & (pattern[first_col]==character))?Alpha:Beta;
        out->dead=trace_gate & (q==0);
        out->score=out->dead?-1000:q;
        out->last=last_at;
        out->matchcol=((q==previous[first_at]+s1) & (s1==Alpha))?first_col:-2;
    }
}

/* Backward rows traverse memory downwards on the same grid. Lanes stay in
 * memory order, so the prefix maximum runs from lane 3 to lane 0 and no lane
 * reversal is needed. The partial highest block is traversed first. */
__attribute__((target("sse4.1")))
static inline __attribute__((always_inline)) void trf_narrowband_scan_backward(
    int *current,const int *previous,const unsigned char *pattern,int size,
    int c,int w,int trace_gate,int character,struct trf_narrowband_row *out)
{
    const int direction=-1;
    __m128i zero=_mm_setzero_si128(),delta=_mm_set1_epi32(Delta);
    __m128i alpha=_mm_set1_epi32(Alpha),beta=_mm_set1_epi32(Beta);
    __m128i ch=_mm_set1_epi32(character);
    __m128i rampstep=_mm_set1_epi32(-4*Delta);
    __m128i gate=_mm_set1_epi32(trace_gate?-1000:0);
    __m128i best_first,best_last,step=_mm_set1_epi32(4);
    __m128i carry,block,next,letters,sim,up,u,p,v,packed,stored,valid,ramp,index,rindex;
    const int full=(2*w)&~3;        /* cells in complete blocks */
    const int n=2*w+1-full;         /* band cells in the highest block: 1 or 3 */
    const int c0=c;
    int low,canonical=(character=='A'||character=='C'||character=='G'||character=='T');
    if(!canonical) alpha=beta;
    /* highest block: lane L is traversal position n-1-L for L<n */
    {
        int c1=trf_narrowband_column(c,1,direction,size);
        int c2=trf_narrowband_column(c,2,direction,size);
        if(n==3) letters=_mm_setr_epi32(pattern[c2],pattern[c1],pattern[c],0);
        else letters=_mm_setr_epi32(pattern[c],0,0,0);
        valid=_mm_loadu_si128((const __m128i *)(trf_narrowband_lanes+4-n));
        index=_mm_and_si128(valid,_mm_setr_epi32(n-1,n-2,n-3,0));
        ramp=_mm_and_si128(valid,_mm_setr_epi32(-(n-1)*Delta,-(n-2)*Delta,-(n-3)*Delta,0));
        rindex=_mm_sub_epi32(_mm_set1_epi32(255),index);
        block=_mm_loadu_si128((const __m128i *)(previous+full));
        if(full>0) next=_mm_loadu_si128((const __m128i *)(previous+full-4));
        else next=zero;
        sim=_mm_blendv_epi8(beta,alpha,_mm_cmpeq_epi32(letters,ch));
        up=_mm_alignr_epi8(block,next,12);
        u=_mm_max_epi32(zero,_mm_max_epi32(_mm_add_epi32(block,sim),_mm_add_epi32(up,delta)));
        u=_mm_add_epi32(_mm_and_si128(u,valid),ramp);
        u=_mm_max_epi32(u,_mm_srli_si128(u,4));
        u=_mm_max_epi32(u,_mm_srli_si128(u,8));
        p=u;
        carry=_mm_shuffle_epi32(p,_MM_SHUFFLE(0,0,0,0));
        v=_mm_sub_epi32(p,ramp);
        packed=_mm_blendv_epi8(_mm_set1_epi32(INT_MIN),_mm_slli_epi32(v,8),valid);
        best_first=_mm_add_epi32(packed,rindex);
        best_last=_mm_add_epi32(packed,index);
        stored=_mm_add_epi32(v,_mm_and_si128(gate,_mm_cmpeq_epi32(v,zero)));
        stored=_mm_blendv_epi8(_mm_loadu_si128((const __m128i *)(current+full)),stored,valid);
        _mm_storeu_si128((__m128i *)(current+full),stored);
        index=_mm_setr_epi32(n+3,n+2,n+1,n);
        rindex=_mm_sub_epi32(_mm_set1_epi32(255),index);
        ramp=_mm_setr_epi32(-(n+3)*Delta,-(n+2)*Delta,-(n+1)*Delta,-n*Delta);
        block=next;
        c=trf_narrowband_column(c,n,direction,size);
    }
    for(low=full-4;low>=0;low-=4) {
        if(c>=3) {
            uint32_t bytes;
            memcpy(&bytes,pattern+c-3,sizeof bytes);
            letters=_mm_cvtepu8_epi32(_mm_cvtsi32_si128((int)bytes));
        } else {
            int c1=trf_narrowband_column(c,1,direction,size);
            int c2=trf_narrowband_column(c,2,direction,size);
            int c3=trf_narrowband_column(c,3,direction,size);
            letters=_mm_setr_epi32(pattern[c3],pattern[c2],pattern[c1],pattern[c]);
        }
        sim=_mm_blendv_epi8(beta,alpha,_mm_cmpeq_epi32(letters,ch));
        if(low>0) next=_mm_loadu_si128((const __m128i *)(previous+low-4));
        else next=zero;         /* cell 0 has no upper neighbour: 0+Delta<=0 */
        up=_mm_alignr_epi8(block,next,12);
        u=_mm_max_epi32(zero,_mm_max_epi32(_mm_add_epi32(block,sim),_mm_add_epi32(up,delta)));
        u=_mm_add_epi32(u,ramp);
        u=_mm_max_epi32(u,_mm_srli_si128(u,4));
        u=_mm_max_epi32(u,_mm_srli_si128(u,8));
        p=_mm_max_epi32(u,carry);
        carry=_mm_shuffle_epi32(p,_MM_SHUFFLE(0,0,0,0));
        v=_mm_sub_epi32(p,ramp);
        packed=_mm_slli_epi32(v,8);
        best_first=_mm_max_epi32(best_first,_mm_add_epi32(packed,rindex));
        best_last=_mm_max_epi32(best_last,_mm_add_epi32(packed,index));
        stored=_mm_add_epi32(v,_mm_and_si128(gate,_mm_cmpeq_epi32(v,zero)));
        _mm_storeu_si128((__m128i *)(current+low),stored);
        index=_mm_add_epi32(index,step);
        rindex=_mm_sub_epi32(rindex,step);
        ramp=_mm_add_epi32(ramp,rampstep);
        block=next;
        c=trf_narrowband_column(c,4,direction,size);
    }
    TRF_NARROWBAND_REDUCE(best_first);
    TRF_NARROWBAND_REDUCE(best_last);
    {
        int first_packed=_mm_cvtsi128_si32(best_first),last_packed=_mm_cvtsi128_si32(best_last);
        int q=first_packed>>8;
        int first_at=255-(first_packed&255),last_at=last_packed&255;
        int first_col=trf_narrowband_column(c0,first_at,direction,size);
        int s1=(canonical & (pattern[first_col]==character))?Alpha:Beta;
        out->dead=trace_gate & (q==0);
        out->score=out->dead?-1000:q;
        out->last=last_at;
        out->matchcol=((q==previous[2*w-first_at]+s1) & (s1==Alpha))?first_col:-2;
    }
}

/* AVX2 row runs.
 *
 * The two functions below execute consecutive ordinary rows (band shifted by
 * exactly one column, canonical band center) of the backward and forward
 * loops of trf_narrowband_vector. They use the same exact integer recurrence,
 * stored cells and FIRST/LAST reductions as the SSE4.1 kernels above, with
 * eight cells per block on a grid of eight-int blocks, and return to the
 * caller's loop, having changed nothing for that row, as soon as the next row
 * is not ordinary. The caller's loop then handles that row as before.
 * Dispatch additionally requires 2w+7<=MAXBANDWIDTH so that whole blocks stay
 * inside a row. The pattern is read from a circularly extended copy, so
 * column wrap-around needs no per-block handling.
 */
#include <immintrin.h>
#define TRF_NARROWBAND_PAD 176      /* >= MAXBANDWIDTH+16 */
static const int32_t trf_narrowband_lanes8[16]={-1,-1,-1,-1,-1,-1,-1,-1,0,0,0,0,0,0,0,0};

/* Loop state shared with trf_narrowband_vector. best_* are the endpoint
 * variables of the respective direction. */
struct trf_narrowband_run {
    unsigned int r;
    int realr,end_of_trace,maxscore,matches_in_diagonal,matchatmax_col,center;
    int best_realr,best_row,best_col,best_center,best_position;
    unsigned int rows;
};

/* Horizontal maxima of two vectors at once: low and high halves of the result. */
__attribute__((target("avx2")))
static inline __attribute__((always_inline)) uint64_t trf_narrowband_hmax8x2(__m256i first,__m256i last)
{
    __m128i f=_mm_max_epi32(_mm256_castsi256_si128(first),_mm256_extracti128_si256(first,1));
    __m128i l=_mm_max_epi32(_mm256_castsi256_si128(last),_mm256_extracti128_si256(last,1));
    __m128i z=_mm_max_epi32(_mm_unpacklo_epi32(f,l),_mm_unpackhi_epi32(f,l));
    z=_mm_max_epi32(z,_mm_shuffle_epi32(z,_MM_SHUFFLE(1,0,3,2)));
    return (uint64_t)_mm_cvtsi128_si64(z);
}

#define TRF_NARROWBAND_CANONICAL(ch) ((ch)=='A'||(ch)=='C'||(ch)=='G'||(ch)=='T')

__attribute__((target("avx2"),noinline))
static void trf_narrowband_run_forward_avx2(struct trf_narrowband_run *loop,
    int start,int size,int w,int tuplesize,const unsigned char *circular,
    int packed_limit)
{
    const __m256i zero=_mm256_setzero_si256(),delta=_mm256_set1_epi32(Delta);
    const __m256i alpha=_mm256_set1_epi32(Alpha),beta=_mm256_set1_epi32(Beta);
    const __m256i iota=_mm256_setr_epi32(0,1,2,3,4,5,6,7);
    const __m256i eight=_mm256_set1_epi32(8),sevens=_mm256_set1_epi32(7);
    const __m256i floor=_mm256_set1_epi32(-(1<<30));
    const __m256i rampstep=_mm256_set1_epi32(-8*Delta);
    const __m256i ramp0=_mm256_mullo_epi32(iota,_mm256_set1_epi32(-Delta));
    const int cells=2*w+1,full=cells&~7,n=cells-full;   /* n is odd: 1..7 */
    const __m256i valid=_mm256_loadu_si256((const __m256i *)(trf_narrowband_lanes8+8-n));
    const __m256i upvalid=_mm256_loadu_si256((const __m256i *)(trf_narrowband_lanes8+9-n));
    const int alpha_s=Alpha,beta_s=Beta,length=Length;
    const unsigned int limit=maxwraplength;
    unsigned int r=loop->r,rows=0;
    int realr=loop->realr,end_of_trace=loop->end_of_trace,maxscore=loop->maxscore;
    int matches_in_diagonal=loop->matches_in_diagonal,matchatmax_col=loop->matchatmax_col;
    int center=loop->center;
    int maxrealrow=loop->best_realr,maxrow=loop->best_row,maxcol=loop->best_col;

    while((!end_of_trace) && (realr<length) && (r<limit))
    {
        int predicted,hidden,recenter,fastcenter,c,t,character,canonical,trace_gate;
        int lastmatchatmax_col,q,first_at,last_at,first_col,s1,dead,score,better,bestcol;
        const unsigned char *letters;
        const int *previous;
        int *current;
        __m256i ch,gate,ramp,index,best_first,best_last,carry;
        __m256i block,next,sim,up,u,p,v,packed,stored;
        uint64_t both;

        /* Is the next row ordinary? Nothing is modified before this is known. */
        if(maxscore>packed_limit) break;    /* row may leave the packed range */
        if(center<0 || center>=size) break;
        predicted=(center+1==size)?0:center+1;
        recenter=(matches_in_diagonal>=tuplesize);
        if(recenter & (matchatmax_col<0)) break;
        fastcenter=(matchatmax_col+1==size)?0:matchatmax_col+1;
        hidden=predicted;
        TRF_NARROWBAND_OPAQUE(hidden);
        if((recenter?fastcenter:hidden)!=hidden) break;

        r++;
        realr++;
        rows++;
        lastmatchatmax_col=matchatmax_col;
        character=(char)Sequence[realr];
        Bandcenter[r]=predicted;
        center=predicted;
        c=predicted-w;
        if(c<0) c+=size;
        letters=circular+c;
        previous=S[r-1];
        current=S[r];
        trace_gate=(realr>=start);
        canonical=TRF_NARROWBAND_CANONICAL(character);
        ch=_mm256_set1_epi32(canonical?character:-1);
        gate=_mm256_set1_epi32(trace_gate?-1000:0);
        ramp=ramp0;
        index=iota;
        best_first=best_last=floor;
        carry=zero;

        block=_mm256_loadu_si256((const __m256i *)previous);
        for(t=0;t<full;t+=8) {
            sim=_mm256_cvtepu8_epi32(_mm_loadl_epi64((const __m128i *)(letters+t)));
            sim=_mm256_blendv_epi8(beta,alpha,_mm256_cmpeq_epi32(sim,ch));
            next=_mm256_loadu_si256((const __m256i *)(previous+t+8));
            up=_mm256_alignr_epi8(_mm256_permute2x128_si256(block,next,0x21),block,4);
            u=_mm256_max_epi32(zero,_mm256_max_epi32(_mm256_add_epi32(block,sim),_mm256_add_epi32(up,delta)));
            u=_mm256_add_epi32(u,ramp);
            u=_mm256_max_epi32(u,_mm256_slli_si256(u,4));
            u=_mm256_max_epi32(u,_mm256_slli_si256(u,8));
            u=_mm256_max_epi32(u,_mm256_shuffle_epi32(_mm256_permute2x128_si256(u,u,0x08),_MM_SHUFFLE(3,3,3,3)));
            p=_mm256_max_epi32(u,carry);
            carry=_mm256_permutevar8x32_epi32(p,sevens);
            v=_mm256_sub_epi32(p,ramp);
            packed=_mm256_slli_epi32(v,8);
            best_first=_mm256_max_epi32(best_first,_mm256_sub_epi32(packed,index));
            best_last=_mm256_max_epi32(best_last,_mm256_add_epi32(packed,index));
            stored=_mm256_add_epi32(v,_mm256_and_si256(gate,_mm256_cmpeq_epi32(v,zero)));
            _mm256_storeu_si256((__m256i *)(current+t),stored);
            index=_mm256_add_epi32(index,eight);
            ramp=_mm256_add_epi32(ramp,rampstep);
            block=next;
        }
        /* highest block: n band cells, the last of which has no upper neighbour */
        sim=_mm256_cvtepu8_epi32(_mm_loadl_epi64((const __m128i *)(letters+t)));
        sim=_mm256_blendv_epi8(beta,alpha,_mm256_cmpeq_epi32(sim,ch));
        up=_mm256_alignr_epi8(_mm256_permute2x128_si256(block,block,0x81),block,4);
        up=_mm256_and_si256(_mm256_add_epi32(up,delta),upvalid);
        u=_mm256_max_epi32(zero,_mm256_max_epi32(_mm256_add_epi32(block,sim),up));
        u=_mm256_add_epi32(u,ramp);
        u=_mm256_max_epi32(u,_mm256_slli_si256(u,4));
        u=_mm256_max_epi32(u,_mm256_slli_si256(u,8));
        u=_mm256_max_epi32(u,_mm256_shuffle_epi32(_mm256_permute2x128_si256(u,u,0x08),_MM_SHUFFLE(3,3,3,3)));
        p=_mm256_max_epi32(u,carry);
        v=_mm256_sub_epi32(p,ramp);
        packed=_mm256_blendv_epi8(floor,_mm256_slli_epi32(v,8),valid);
        best_first=_mm256_max_epi32(best_first,_mm256_sub_epi32(packed,index));
        best_last=_mm256_max_epi32(best_last,_mm256_add_epi32(packed,index));
        stored=_mm256_add_epi32(v,_mm256_and_si256(gate,_mm256_cmpeq_epi32(v,zero)));
        stored=_mm256_blendv_epi8(_mm256_loadu_si256((const __m256i *)(current+t)),stored,valid);
        _mm256_storeu_si256((__m256i *)(current+t),stored);

        both=trf_narrowband_hmax8x2(best_first,best_last);
        q=((int)(uint32_t)both+255)>>8;
        first_at=255-(((int)(uint32_t)both+255)&255);
        last_at=(int)((both>>32)&255);
        first_col=c+first_at;
        if(first_col>=size) first_col-=size;
        s1=(canonical & (letters[first_at]==character))?alpha_s:beta_s;
        dead=trace_gate & (q==0);
        score=dead?-1000:q;

        end_of_trace=dead;
        better=(score>=maxscore);
        bestcol=c+last_at;
        if(bestcol>=size) bestcol-=size;
        maxscore=better?score:maxscore;
        maxrealrow=better?realr:maxrealrow;
        maxrow=better?(int)r:maxrow;
        maxcol=better?bestcol:maxcol;
        if(!dead) matchatmax_col=((q==previous[first_at]+s1) & (s1==alpha_s))?first_col:-2;
        matches_in_diagonal=(trf_narrowband_remainder(matchatmax_col-lastmatchatmax_col+size,size)==1)?
            matches_in_diagonal+1:0;
    }

    loop->r=r; loop->realr=realr; loop->end_of_trace=end_of_trace;
    loop->maxscore=maxscore; loop->matches_in_diagonal=matches_in_diagonal;
    loop->matchatmax_col=matchatmax_col; loop->center=center;
    loop->best_realr=maxrealrow; loop->best_row=maxrow; loop->best_col=maxcol;
    loop->rows=rows;
}

/* Backward rows start at cell 2w and descend. Lanes stay in memory order:
 * lane L of the block at cell B is traversal position 2w-B-L. */
__attribute__((target("avx2"),noinline))
static void trf_narrowband_run_backward_avx2(struct trf_narrowband_run *loop,
    int start,int size,int w,int tuplesize,const unsigned char *circular,
    int packed_limit)
{
    const __m256i zero=_mm256_setzero_si256(),delta=_mm256_set1_epi32(Delta);
    const __m256i alpha=_mm256_set1_epi32(Alpha),beta=_mm256_set1_epi32(Beta);
    const __m256i iota=_mm256_setr_epi32(0,1,2,3,4,5,6,7);
    const __m256i eight=_mm256_set1_epi32(8);
    const __m256i floor=_mm256_set1_epi32(-(1<<30));
    const __m256i rampstep=_mm256_set1_epi32(-8*Delta);
    const __m256i slope=_mm256_set1_epi32(-Delta);
    const int cells=2*w+1,full=cells&~7,n=cells-full;   /* n is odd: 1..7 */
    const __m256i valid=_mm256_loadu_si256((const __m256i *)(trf_narrowband_lanes8+8-n));
    const __m256i index_high=_mm256_and_si256(valid,_mm256_sub_epi32(_mm256_set1_epi32(n-1),iota));
    const __m256i ramp_high=_mm256_mullo_epi32(index_high,slope);
    const __m256i index_next=_mm256_sub_epi32(_mm256_set1_epi32(n+7),iota);
    const __m256i ramp_next=_mm256_mullo_epi32(index_next,slope);
    const int alpha_s=Alpha,beta_s=Beta;
    unsigned int r=loop->r,rows=0;
    int realr=loop->realr,end_of_trace=loop->end_of_trace,maxscore=loop->maxscore;
    int matches_in_diagonal=loop->matches_in_diagonal,matchatmax_col=loop->matchatmax_col;
    int center=loop->center;
    int minrealrow=loop->best_realr,mincol=loop->best_col;
    int mincolbandcenter=loop->best_center,mincolposition=loop->best_position;

    while((!end_of_trace) && (realr>1) && (r>0))
    {
        int predicted,hidden,recenter,fastcenter,c,low,character,canonical,trace_gate;
        int lastmatchatmax_col,q,first_at,last_at,first_col,s1,dead,score,better,bestcol;
        const unsigned char *letters;
        const int *previous;
        int *current;
        __m256i ch,gate,ramp,index,best_first,best_last,carry;
        __m256i block,next,sim,up,u,p,v,packed,stored;
        uint64_t both;

        /* Is the next row ordinary? Nothing is modified before this is known. */
        if(maxscore>packed_limit) break;    /* row may leave the packed range */
        if(center<0 || center>=size) break;
        predicted=(center==0)?size-1:center-1;
        recenter=(matches_in_diagonal>=tuplesize);
        if(recenter & (matchatmax_col<0)) break;
        fastcenter=(matchatmax_col==0)?size-1:matchatmax_col-1;
        hidden=predicted;
        TRF_NARROWBAND_OPAQUE(hidden);
        if((recenter?fastcenter:hidden)!=hidden) break;

        r--;
        realr--;
        rows++;
        lastmatchatmax_col=matchatmax_col;
        character=(char)Sequence[realr];
        Bandcenter[r]=predicted;
        center=predicted;
        c=predicted+w;
        if(c>=size) c-=size;
        letters=circular+c-2*w;         /* letters[i] belongs to band cell i */
        previous=S[r+1];
        current=S[r];
        trace_gate=(realr<=start-size);
        canonical=TRF_NARROWBAND_CANONICAL(character);
        ch=_mm256_set1_epi32(canonical?character:-1);
        gate=_mm256_set1_epi32(trace_gate?-1000:0);

        /* highest block first */
        block=_mm256_loadu_si256((const __m256i *)(previous+full));
        if(full>0) next=_mm256_loadu_si256((const __m256i *)(previous+full-8));
        else next=zero;
        sim=_mm256_cvtepu8_epi32(_mm_loadl_epi64((const __m128i *)(letters+full)));
        sim=_mm256_blendv_epi8(beta,alpha,_mm256_cmpeq_epi32(sim,ch));
        up=_mm256_alignr_epi8(block,_mm256_permute2x128_si256(next,block,0x21),12);
        u=_mm256_max_epi32(zero,_mm256_max_epi32(_mm256_add_epi32(block,sim),_mm256_add_epi32(up,delta)));
        u=_mm256_add_epi32(_mm256_and_si256(u,valid),ramp_high);
        u=_mm256_max_epi32(u,_mm256_srli_si256(u,4));
        u=_mm256_max_epi32(u,_mm256_srli_si256(u,8));
        u=_mm256_max_epi32(u,_mm256_shuffle_epi32(_mm256_permute2x128_si256(u,u,0x81),_MM_SHUFFLE(0,0,0,0)));
        p=u;
        carry=_mm256_broadcastd_epi32(_mm256_castsi256_si128(p));
        v=_mm256_sub_epi32(p,ramp_high);
        packed=_mm256_blendv_epi8(floor,_mm256_slli_epi32(v,8),valid);
        best_first=_mm256_sub_epi32(packed,index_high);
        best_last=_mm256_add_epi32(packed,index_high);
        stored=_mm256_add_epi32(v,_mm256_and_si256(gate,_mm256_cmpeq_epi32(v,zero)));
        stored=_mm256_blendv_epi8(_mm256_loadu_si256((const __m256i *)(current+full)),stored,valid);
        _mm256_storeu_si256((__m256i *)(current+full),stored);
        index=index_next;
        ramp=ramp_next;
        block=next;
        for(low=full-8;low>=0;low-=8) {
            sim=_mm256_cvtepu8_epi32(_mm_loadl_epi64((const __m128i *)(letters+low)));
            sim=_mm256_blendv_epi8(beta,alpha,_mm256_cmpeq_epi32(sim,ch));
            if(low>0) next=_mm256_loadu_si256((const __m256i *)(previous+low-8));
            else next=zero;         /* cell 0 has no upper neighbour: 0+Delta<=0 */
            up=_mm256_alignr_epi8(block,_mm256_permute2x128_si256(next,block,0x21),12);
            u=_mm256_max_epi32(zero,_mm256_max_epi32(_mm256_add_epi32(block,sim),_mm256_add_epi32(up,delta)));
            u=_mm256_add_epi32(u,ramp);
            u=_mm256_max_epi32(u,_mm256_srli_si256(u,4));
            u=_mm256_max_epi32(u,_mm256_srli_si256(u,8));
            u=_mm256_max_epi32(u,_mm256_shuffle_epi32(_mm256_permute2x128_si256(u,u,0x81),_MM_SHUFFLE(0,0,0,0)));
            p=_mm256_max_epi32(u,carry);
            carry=_mm256_broadcastd_epi32(_mm256_castsi256_si128(p));
            v=_mm256_sub_epi32(p,ramp);
            packed=_mm256_slli_epi32(v,8);
            best_first=_mm256_max_epi32(best_first,_mm256_sub_epi32(packed,index));
            best_last=_mm256_max_epi32(best_last,_mm256_add_epi32(packed,index));
            stored=_mm256_add_epi32(v,_mm256_and_si256(gate,_mm256_cmpeq_epi32(v,zero)));
            _mm256_storeu_si256((__m256i *)(current+low),stored);
            index=_mm256_add_epi32(index,eight);
            ramp=_mm256_add_epi32(ramp,rampstep);
            block=next;
        }

        both=trf_narrowband_hmax8x2(best_first,best_last);
        q=((int)(uint32_t)both+255)>>8;
        first_at=255-(((int)(uint32_t)both+255)&255);
        last_at=(int)((both>>32)&255);
        first_col=c-first_at;
        if(first_col<0) first_col+=size;
        s1=(canonical & (letters[2*w-first_at]==character))?alpha_s:beta_s;
        dead=trace_gate & (q==0);
        score=dead?-1000:q;

        end_of_trace=dead;
        better=(score>=maxscore);
        bestcol=c-last_at;
        if(bestcol<0) bestcol+=size;
        maxscore=better?score:maxscore;
        minrealrow=better?realr:minrealrow;
        mincol=better?bestcol:mincol;
        mincolbandcenter=better?predicted:mincolbandcenter;
        mincolposition=better?2*w-last_at:mincolposition;
        if(!dead) matchatmax_col=((q==previous[2*w-first_at]+s1) & (s1==alpha_s))?first_col:-2;
        matches_in_diagonal=(trf_narrowband_remainder(matchatmax_col-lastmatchatmax_col+size,size)==size-1)?
            matches_in_diagonal+1:0;
    }

    loop->r=r; loop->realr=realr; loop->end_of_trace=end_of_trace;
    loop->maxscore=maxscore; loop->matches_in_diagonal=matches_in_diagonal;
    loop->matchatmax_col=matchatmax_col; loop->center=center;
    loop->best_realr=minrealrow; loop->best_col=mincol;
    loop->best_center=mincolbandcenter; loop->best_position=mincolposition;
    loop->rows=rows;
}

static int trf_narrowband_avx2_domain(int wb,int wf)
{
#ifdef TRF_NARROWBAND_DISABLE_AVX2
    (void)wb; (void)wf;
    return 0;
#else
    return 2*wb+7<=MAXBANDWIDTH && 2*wf+7<=MAXBANDWIDTH &&
           trf_narrowband_ramp_domain(wb,wf) && __builtin_cpu_supports("avx2")!=0;
#endif
}

__attribute__((target("sse4.1")))
void trf_narrowband_vector(int start, int size,int bandradius, int bandradiusforward, 
		int option,int tuplesize)

/* start is end of pattern in text */
/* tuplesize is the size of tuple used for this pattern size */

/* started Feb. 7 1997 */
{
	int g;
	register int *pup, *pdiag, *pcurr,pleft;
	int c,realr,end_of_trace,
	    maxscore,maxrow,maxcol,mincol;
	int maxrealrow, minrealrow;
	char currchar;
	const int *similarity_row;
	/* Diag/Up originally copy the preceding S row and S+Delta. Each
	 * diagonal index is consumed once, so the read-only row plus local
	 * diag_score preserves the immediate match observer without stores. */
	const int *previous, *previous_diag, *previous_up;
	int diag_score;
	int w,matches_in_diagonal,matchatmax_col,i,k,maxrowscore,
	    lastmatchatmax_col,match_yes_no;
	/*** 6/9/05 G. Benson ***/ 
	int mincolbandcenter,zeroat,mincolposition;
	/* Feb 16, 2016 Yozen */
	unsigned int r;

	/* G. Benson */
	/* 1/26/10 */
	/* change MAXWRAPLENGTH to MAXWRAPLENGTHCONST so MAXWRAPLENGTH can be used as an int */
	/* int Bandcenter[MAXWRAPLENGTH+1]; */
	// if (Bandcenter == NULL) {
	// 	Bandcenter = calloc(maxwraplength+1, sizeof(*Bandcenter));
	// }
	/* Bandcenter init removed by Yozen on Feb 16, 2016. */

	w=bandradius;
	if(MAXBANDWIDTH<2*w+1)
	{
		trf_message("\nIn narrowbandwrap, MAXBANDWIDTH: %d exceeded by 2*w+1: %d\n",
				MAXBANDWIDTH,2*w+1);
		exit(-1);
	}
	const int ramp_rows=trf_narrowband_ramp_domain(bandradius,bandradiusforward);
	const int packed_limit=trf_narrowband_packed_limit();
	const int avx2_rows=trf_narrowband_avx2_domain(bandradius,bandradiusforward);
	/* circular copy of the pattern for the AVX2 rows */
	unsigned char extended[4000+2*TRF_NARROWBAND_PAD];
	unsigned char *const circular=extended+TRF_NARROWBAND_PAD;
	int center,predicted,hidden,recenter,fastcenter,better,bestcol;

	/* fill EC */
	if(option==WITHCONSENSUS)
		for(g=0;g<size;g++) EC[g]=Consensus.pattern[g];
	else
		for(g=0;g<size;g++) EC[g]=Sequence[start-size+g+1];
	if(avx2_rows)
	{
		memcpy(circular,EC,(size_t)size);
		for(g=0;g<TRF_NARROWBAND_PAD;g++) circular[size+g]=circular[g];
		for(g=1;g<=TRF_NARROWBAND_PAD;g++) circular[-g]=circular[size-g];
	}

	/* The selections below read these before their first real assignment. */
	maxrow=maxcol=mincol=maxrealrow=minrealrow=mincolbandcenter=mincolposition=0;

	/* backward wdp */
	maxscore=0;
	realr=start+1;
	r=maxwraplength;
	Bandcenter[r]=0;
	center=0;
	matches_in_diagonal=0;
	matchatmax_col=-2;

	pcurr=&S[r][0];
	pdiag=&Diag[0];
	pup=&Up[0];




	/* 3/14/05 gary benson -- reverse direction */
	/* change zeroth row values to put in -1000 in unreachable cells
	   and gap penalty in cells beyond start location */


	for (i=0;i<=w;i++)
	{
		*pup=(*pdiag=(*pcurr=0+Delta*(w-i)))+ Delta;
		pup++;
		pdiag++;
		pcurr++;
	}

	for (i=w+1;i<=2*w;i++)
	{
		*pup=(*pdiag=(*pcurr=-1000))+ Delta;
		pup++;
		pdiag++;
		pcurr++;
	}





	/* backwards */

	end_of_trace=FALSE;  
	while ((!end_of_trace) && (realr>1) && (r>0))
	{
		if(avx2_rows)
		{
			/* run ordinary rows; fall through for any other row */
			struct trf_narrowband_run run;
			run.r=r; run.realr=realr; run.end_of_trace=end_of_trace;
			run.maxscore=maxscore; run.matches_in_diagonal=matches_in_diagonal;
			run.matchatmax_col=matchatmax_col; run.center=center;
			run.best_realr=minrealrow; run.best_row=0; run.best_col=mincol;
			run.best_center=mincolbandcenter; run.best_position=mincolposition;
			trf_narrowband_run_backward_avx2(&run,start,size,w,tuplesize,circular,packed_limit);
			r=run.r; realr=run.realr; end_of_trace=run.end_of_trace;
			maxscore=run.maxscore; matches_in_diagonal=run.matches_in_diagonal;
			matchatmax_col=run.matchatmax_col; center=run.center;
			minrealrow=run.best_realr; mincol=run.best_col;
			mincolbandcenter=run.best_center; mincolposition=run.best_position;
			Rows+=run.rows;
			if(!((!end_of_trace) && (realr>1) && (r>0))) break;
		}
		r--;
		realr--;
		Rows++;
		maxrowscore=-1;
		lastmatchatmax_col=matchatmax_col;
		end_of_trace=TRUE;
		currchar=Sequence[realr];
		similarity_row=SM+256*currchar;
		previous=S[r+1];
		pcurr=&S[r][2*w];
		pleft=-1000;

		/* center mirrors Bandcenter[r+1]. The recentered value nearly always
		 * equals the unrecentered prediction; branching on that equality and
		 * continuing from the prediction keeps the next row's addresses
		 * independent of this row's reduction. Values are unchanged. */
		predicted=trf_narrowband_remainder(center-1+size,size);
		recenter=(matches_in_diagonal>=tuplesize);
		if(__builtin_expect(recenter & (matchatmax_col<0),0))
		{
			/* recenter band from the -2 sentinel: keep C's remainder */
			Bandcenter[r]=trf_narrowband_remainder(matchatmax_col-1+size,size);
			TRF_NARROWBAND_BARRIER();
		}
		else
		{
			/* recenter band (a canonical column wraps once), or don't */
			fastcenter=(matchatmax_col==0)?size-1:matchatmax_col-1;
			Bandcenter[r]=recenter?fastcenter:predicted;
		}
		hidden=predicted;
		TRF_NARROWBAND_OPAQUE(hidden);
		if(__builtin_expect(Bandcenter[r]!=hidden,0))
		{
			predicted=Bandcenter[r];
			TRF_NARROWBAND_BARRIER();
		}

		/* change of bandcenter determines which inputs go into which
		   cells */

		k=trf_narrowband_remainder(predicted-center+size,size);
		if(size-k<=k) k=-(size-k);
		c=trf_narrowband_remainder(predicted+w,size);
		center=predicted;
		if(k==-1 && c>=0 && c<size)
        {
            struct trf_narrowband_row q;
            if(ramp_rows && maxscore<=packed_limit)
                trf_narrowband_scan_backward(S[r],previous,EC,size,c,w,
                                realr<=start-size,currchar,&q);
            else
            trf_narrowband_scan(S[r],previous,EC,size,c,w,-1,
                                realr<=start-size,currchar,&q);
            end_of_trace=q.dead;
            /* selections rather than branches: these outcomes follow the
             * sequence, and a late misprediction discards the next rows */
            better=(q.score>=maxscore);
            bestcol=trf_narrowband_column(c,q.last,-1,size);
            maxscore=better?q.score:maxscore;
            minrealrow=better?realr:minrealrow;
            mincol=better?bestcol:mincol;
            mincolbandcenter=better?predicted:mincolbandcenter;
            mincolposition=better?2*w-q.last:mincolposition;
            matchatmax_col=q.dead?matchatmax_col:q.matchcol;
        }
        else if(k<=-1)           /* band shifts left */
		{
			k=-k;
			previous_diag=&previous[2*w-k+1];
			previous_up=&previous[2*w-k];
			for(i=2*w;i>=k;i--)
			{
				diag_score=*previous_diag+(match_yes_no=similarity_row[EC[c]]);
				pleft=(*pcurr=max4(0,diag_score,*previous_up+Delta,pleft))+Delta;
				test_trace_and_backwards_maxscore;
				test_maxrowscore_with_match;
				pcurr--;
				previous_diag--;
				previous_up--;
				c=(c==0)?size-1:c-1;
			}
			i=k-1;
			diag_score=*previous_diag+(match_yes_no=similarity_row[EC[c]]);
			pleft=(*pcurr=max3(0,diag_score,pleft))+Delta;
			test_trace_and_backwards_maxscore;
			test_maxrowscore_with_match;
			pcurr--;
			c=(c==0)?size-1:c-1;
			for(i=k-2;i>=0;i--)
			{
				pleft=(*pcurr=max2(0,pleft))+Delta;
				test_trace_and_backwards_maxscore;
				test_maxrowscore_without_match;
				pcurr--;
				c=(c==0)?size-1:c-1;
			}
		}
		else            /* band shifts right */
		{
			previous_diag=&previous[2*w];
			previous_up=&previous[2*w];
			for(i=2*w;i>=2*w-k+1;i--)
			{
				pleft=(*pcurr=max2(0,pleft))+Delta;
				test_trace_and_backwards_maxscore;
				test_maxrowscore_without_match;
				pcurr--;
				c=(c==0)?size-1:c-1;
			}
			i=2*w-k;
			pleft=(*pcurr=max3(0,*previous_up+Delta,pleft))+Delta;
			test_trace_and_backwards_maxscore;
			test_maxrowscore_without_match;
			pcurr--;
			previous_up--;
			c=(c==0)?size-1:c-1;
			for(i=2*w-k-1;i>=0;i--)
			{
				diag_score=*previous_diag+(match_yes_no=similarity_row[EC[c]]);
				pleft=(*pcurr=max4(0,diag_score,*previous_up+Delta,pleft))+Delta;
				test_trace_and_backwards_maxscore;
				test_maxrowscore_with_match;
				pcurr--;
				previous_diag--;
				previous_up--;
				c=(c==0)?size-1:c-1;
			}
		}
		matches_in_diagonal=(trf_narrowband_remainder(matchatmax_col-lastmatchatmax_col+size,size)==size-1)?
			matches_in_diagonal+1:0;
	}
	/* Restore backward scratch once; forward may initialize a narrower band. */
	pcurr=&S[r][0];
	pdiag=&Diag[0];
	pup=&Up[0];
	for (i=0;i<=2*w;i++)
	{
		*pup=(*pdiag=*pcurr)+Delta;
		pcurr++;
		pdiag++;
		pup++;
	}


	/**********************************************************************/
	/**********************************************************************/
	/**********************************************************************/
	/**********************************************************************/

	/* forward */

	if (ldong==1) /* null statement */;
	else
	{

		r=0;
		realr=minrealrow-1;
		/* mincol matches minrealrow in previous alignment above */

		/*** 6/14/05 G. Benson ***/
		/* modification version g */
		Bandcenter[0]=trf_narrowband_remainder(mincolbandcenter-1+size,size);
		center=Bandcenter[0];
		matches_in_diagonal=0;
		matchatmax_col=-2;
		maxscore=0;
		pup=&Up[0];
		pdiag=&Diag[0];
		pcurr=&S[0][0];




		/* 3/14/05 gary benson -- forward direction */
		/* change zeroth row values to put in -1000 in unreachable cells
		   and gap penalty in cells beyond start location */

		/* initialize top row */

		/*** 6/14/05 G. Benson ***/
		/* set zero at mincol, not at bandcenter, which is now mincolbandcenter, 
		   so that forward alignment starts out the same as backwards alignment ended. */ 
		zeroat=mincolposition-w;



		/* July 9, 20009, G. Benson */
		/* widening the narrowband in the forward direction to catch alignments that don't return to correct start column */
		w = bandradiusforward;

		if(MAXBANDWIDTH<2*w+1)
		{
			trf_message("\nIn narrowbandwrap, MAXBANDWIDTH: %d exceeded by bandradiusforward 2*w+1: %d\n",
					MAXBANDWIDTH,2*w+1);
			exit(-1);
		}

		//	if((zeroat>w)||(zeroat<-w)) fprintf(stderr,"\nfound zeroat out of bounds");
		/*** 6/14/05 G. Benson restore what was in all earlier versions ***/
		for(i=0;i<w+zeroat;i++)
		{
			*pup=(*pdiag=(*pcurr=-1000))+Delta;
			pup++;
			pdiag++;
			pcurr++;
		}

		/*** 6/9/05 G. Benson ***/
		for(i=w+zeroat;i<=2*w;i++)
		{
			/*** 6/9/05 G. Benson ***/
			*pup=(*pdiag=(*pcurr=0+Delta*(i-(w+zeroat))))+Delta;
			pup++;
			pdiag++;
			pcurr++;
		}




		/* compute until end of trace */      
		end_of_trace=FALSE;  
		while ((!end_of_trace) && (realr<Length) && (r<maxwraplength))
		{
			if(avx2_rows)
			{
				/* run ordinary rows; fall through for any other row */
				struct trf_narrowband_run run;
				run.r=r; run.realr=realr; run.end_of_trace=end_of_trace;
				run.maxscore=maxscore; run.matches_in_diagonal=matches_in_diagonal;
				run.matchatmax_col=matchatmax_col; run.center=center;
				run.best_realr=maxrealrow; run.best_row=maxrow; run.best_col=maxcol;
				run.best_center=0; run.best_position=0;
				trf_narrowband_run_forward_avx2(&run,start,size,w,tuplesize,circular,packed_limit);
				r=run.r; realr=run.realr; end_of_trace=run.end_of_trace;
				maxscore=run.maxscore; matches_in_diagonal=run.matches_in_diagonal;
				matchatmax_col=run.matchatmax_col; center=run.center;
				maxrealrow=run.best_realr; maxrow=run.best_row; maxcol=run.best_col;
				Rows+=run.rows;
				if(!((!end_of_trace) && (realr<Length) && (r<maxwraplength))) break;
			}
			r++;
			realr++;
			Rows++;
			end_of_trace=TRUE;
			maxrowscore=-1;
			lastmatchatmax_col=matchatmax_col;
			currchar=Sequence[realr];
			similarity_row=SM+256*currchar;
			previous=S[r-1];
			pcurr=&S[r][0];
			pleft=-1000;      /* don't use pleft for first entry */   

			/* center mirrors Bandcenter[r-1]; see the backward loop. */
			predicted=trf_narrowband_remainder(center+1,size);
			recenter=(matches_in_diagonal>=tuplesize);
			if(__builtin_expect(recenter & (matchatmax_col<0),0))
			{
				/* recenter band from the -2 sentinel: keep C's remainder */
				Bandcenter[r]=trf_narrowband_remainder(matchatmax_col+1,size);
				TRF_NARROWBAND_BARRIER();
			}
			else
			{
				/* recenter band (a canonical column wraps once), or don't */
				fastcenter=(matchatmax_col+1==size)?0:matchatmax_col+1;
				Bandcenter[r]=recenter?fastcenter:predicted;
			}
			hidden=predicted;
			TRF_NARROWBAND_OPAQUE(hidden);
			if(__builtin_expect(Bandcenter[r]!=hidden,0))
			{
				predicted=Bandcenter[r];
				TRF_NARROWBAND_BARRIER();
			}

			/* change of bandcenter determines which inputs go into which
			   cells */

			k=trf_narrowband_remainder(predicted-center+size,size);
			if(size-k<=k) k=-(size-k);
            c=trf_narrowband_remainder(predicted-w+size,size);
			center=predicted;
            if(k==1 && c>=0 && c<size)
            {
                struct trf_narrowband_row q;
                if(ramp_rows && maxscore<=packed_limit)
                    trf_narrowband_scan_forward(S[r],previous,EC,size,c,w,
                                    realr>=start,currchar,&q);
                else
                trf_narrowband_scan(S[r],previous,EC,size,c,w,1,
                                    realr>=start,currchar,&q);
                end_of_trace=q.dead;
                better=(q.score>=maxscore);
                bestcol=trf_narrowband_column(c,q.last,1,size);
                maxscore=better?q.score:maxscore;
                maxrealrow=better?realr:maxrealrow;
                maxrow=better?(int)r:maxrow;
                maxcol=better?bestcol:maxcol;
                matchatmax_col=q.dead?matchatmax_col:q.matchcol;
            }
            else if(k>=1)          /* band shifts right */
			{
				previous_diag=&previous[k-1];
				previous_up=&previous[k];
				c=trf_narrowband_remainder(Bandcenter[r]-w+size,size);
				for(i=0;i<=2*w-k;i++)
				{
					diag_score=*previous_diag+(match_yes_no=similarity_row[EC[c]]);
					pleft=(*pcurr=max4(0,diag_score,*previous_up+Delta,pleft))+Delta;
					test_trace_and_forward_maxscore;
					test_maxrowscore_with_match;
					pcurr++;
					previous_diag++;
					previous_up++;
					c=(c+1==size)?0:c+1;
				}
				i=2*w-k+1;
				diag_score=*previous_diag+(match_yes_no=similarity_row[EC[c]]);
				pleft=(*pcurr=max3(0,diag_score,pleft))+Delta;
				test_trace_and_forward_maxscore;
				test_maxrowscore_with_match;
				pcurr++;
				c=(c+1==size)?0:c+1; 
				for(i=2*w-k+2;i<=2*w;i++)
				{
					pleft=(*pcurr=max2(0,pleft))+Delta;
					test_trace_and_forward_maxscore;
					test_maxrowscore_without_match;
					pcurr++;
					c=(c+1==size)?0:c+1;
				}
			}
			else          /* band shifts left */
			{
				k=-k;
				c=trf_narrowband_remainder(Bandcenter[r]-w+size,size);
				previous_up=&previous[0];
				previous_diag=&previous[0];
				for(i=0;i<=k-1;i++)
				{
					pleft=(*pcurr=max2(0,pleft))+Delta;
					test_trace_and_forward_maxscore;
					test_maxrowscore_without_match;
					pcurr++;
					c=(c+1==size)?0:c+1;
				}
				i=k;
				pleft=(*pcurr=max3(0,*previous_up+Delta,pleft))+Delta;
				test_trace_and_forward_maxscore;
				test_maxrowscore_without_match;
				pcurr++;
				previous_up++;
				c=(c+1==size)?0:c+1;
				for(i=k+1;i<=2*w;i++)
				{
					diag_score=*previous_diag+(match_yes_no=similarity_row[EC[c]]);
					pleft=(*pcurr=max4(0,diag_score,*previous_up+Delta,pleft))+Delta;
					test_trace_and_forward_maxscore;
					test_maxrowscore_with_match;
					pcurr++;
					previous_diag++;
					previous_up++;
					c=(c+1==size)?0:c+1;
				}
			}

			matches_in_diagonal=(trf_narrowband_remainder(matchatmax_col-lastmatchatmax_col+size,size)==1)?
				matches_in_diagonal+1:0;
		}
		/* Restore the caller-visible scratch from the last forward row. */
		pcurr=&S[r][0];
		pdiag=&Diag[0];
		pup=&Up[0];
		for (i=0;i<=2*w;i++)
		{
			*pup=(*pdiag=*pcurr)+Delta;
			pcurr++;
			pdiag++;
			pup++;
		}
		/* store last position scanned with this pattern */
		Wrapend=realr;

		if (Criteria_print) trf_message("   DWPlength:%d",Wrapend-(minrealrow-1)+1);
		/* test for report */
		if (maxscore>=Reportmin)
		{
			Maxrealrow=maxrealrow;
			Maxrow=maxrow;
			Maxcol=maxcol;
			Maxscore=maxscore;

		}

	}

}

#endif
