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

/* Narrowband wraparound alignment: selection of the kernel for one call. */

#include "alignment/narrowband.h"
#include "alignment/narrowband_scalar.h"
#include "alignment/narrowband_simd.h"

void narrowbandwrap(int start,int size,int bandradius,int bandradiusforward,
                    int option,int tuplesize)
{
#if TRF_NARROWBAND_HAVE_SIMD
    if(trf_narrowband_vector_domain(start,size,bandradius,bandradiusforward)) {
        trf_narrowband_vector(start,size,bandradius,bandradiusforward,option,tuplesize);
        return;
    }
#endif
    trf_narrowband_scalar(start,size,bandradius,bandradiusforward,option,tuplesize);
}
