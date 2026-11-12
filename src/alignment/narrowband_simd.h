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

/* SSE4.1 and AVX2 narrowband kernel and the test for the calls it accepts. */

#ifndef TRF_ALIGNMENT_NARROWBAND_SIMD_H
#define TRF_ALIGNMENT_NARROWBAND_SIMD_H

#include <limits.h>

#if (defined(__i386__) || defined(__x86_64__)) && \
    (defined(__GNUC__) || defined(__clang__)) && INT_MAX == 2147483647 && \
    !defined(TRF_NARROWBAND_DISABLE_SIMD)
#define TRF_NARROWBAND_HAVE_SIMD 1

/* nonzero if trf_narrowband_vector computes this call exactly on this CPU */
int trf_narrowband_vector_domain(int start, int size, int wb, int wf);
void trf_narrowband_vector(int start, int size, int bandradius, int bandradiusforward,
                           int option, int tuplesize);
#else
#define TRF_NARROWBAND_HAVE_SIMD 0
#endif

#endif
