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

/* Windows of k-tuple matches for each distance, and the links between active distances. */

#include "core/globals.h"
#include "core/report.h"
#include "detection/distance_list.h"

/****************************************************/

struct distancelist *new_distancelist()
{
	int g,N,K;
	struct distanceentry *ptr;

	struct distancelist *objptr=(struct distancelist *)
		calloc(MAXDISTANCE+1,sizeof(struct distancelist));

	K = Min_Distance_Entries + 1; 
	N = MAXDISTANCE + 1;
	//ptr = _DistanceEntries = malloc( (((N*(N+1) - K*(K+1))/2) + K*(K-1)) * sizeof(struct distanceentry) );
	ptr = _DistanceEntries = malloc( ((K+N)*(N-K+1)/2 + K*(K-1) ) * sizeof(struct distanceentry) );

	for(g=1;g<=MAXDISTANCE;g++) {
		objptr[g].entry = ptr;
		memset(objptr[g].entry, 0, (max(g,Min_Distance_Entries)+1) * sizeof(struct distanceentry));
		ptr += (max(g,Min_Distance_Entries)+1);
	}

	/*

	   for(g=1;g<=MAXDISTANCE;g++)
	   {
	   objptr[g].entry=(struct distanceentry *)
	   calloc(max(g,Min_Distance_Entries)+1,sizeof(struct distanceentry));
	   }
	 */  

	return(objptr);
}

void clear_distancelist(struct distancelist *objptr)
{
	int g;
	for(g=1;g<=MAXDISTANCE;g++)
	{
		objptr[g].lowindex=0;    
		objptr[g].highindex=max(g,Min_Distance_Entries);
		objptr[g].numentries=0;
		objptr[g].nummatches=0;
	}
}
void init_links()
{
	Distance[0].linkup=MAXDISTANCE+1;
}

/* Exact, division-free form of a hot remainder. It falls back to the
 * original expression outside the range in which the shortcut is exact. */
static inline int trf_wrap_index(int value, int modulus)
{
	/* value%modulus for an index that was in range before one increment */
	if(value>=0 && value<modulus) return value;
	if(value==modulus) return 0;
	return value%modulus;
}

void add_tuple_match_to_Distance_entry(int location, int size, int d,
		struct distancelist *objptr)
{
	int *lo, *hi;
	int *z, *m;
	int windowleftend, windowsize;
	struct distancelist *list;

	/* get distance window pointers */
	windowsize=max(d,Min_Distance_Entries)+1; /* this value is used to */
	/* mod the index to the */
	/* entries.  The entries */
	/* run from 0 to */
	/* max(d,Min_Distance_Entries) */ 
	list=&(objptr[d]);
	hi=&(list->highindex);
	lo=&(list->lowindex);
	z=&(list->numentries);
	m=&(list->nummatches);
	if ((*z)!=0) /* there exist previous entries on this list */
	{
		/* first remove trailing entries, i.e., */
		/* those with location<location-d or */
		/* <location-Min_Distance_Window, whichever is further back */
		windowleftend=min(location-d,location-Min_Distance_Window)+1;
		/* note that preceding is based on Min_Distance_Window */
		/* rather than Min_Distance_Entries */

		while((list->entry[*lo].location<windowleftend)&&((*z)>0))
		{
			(*z)--;
			(*m)-=list->entry[*lo].size;
			(*lo)++;
			(*lo)=trf_wrap_index(*lo,windowsize);
		}

	}
	if(((*z)!=0)&&(list->entry[*hi].location==location-1)) 
		/* there are still more entries and this is an adjacent tuple, */
		/* just add on to last entry */
	{
		list->entry[*hi].location++;
		list->entry[*hi].size++;
		(*m)++;
	}
	else          /* need a new entry here */
	{
		(*z)++;
		(*hi)++;
		(*hi)=trf_wrap_index(*hi,windowsize);
		list->entry[*hi].location=location;
		list->entry[*hi].size=size;
		(*m)+=size;
	}
}

void link_Distance_window(int d)
{

	int t,f,h;
	/* get next highest tag */
	t=d/TAGSEP;	/* was (int)ceil(d/TAGSEP): integer division, ceil is the identity */
	if(Tag[t]<d)            /* Tag[t] is the largest index less or */
		/* equal to (t)x(TAGSEP) that is linked */
	{
		f=Tag[t];           /* f is first linked index below d */
		while((t<=Toptag)&&(Tag[t]<d))
		{
			Tag[t]=d;
			t++;
		}
	}
	else if(Tag[t]>d)     /* follow links to insert d */
	{
		f=Tag[t];
		while(f>d) f=Distance[f].linkdown;
		/* f is first linked index below d */
		if (f==d)trf_message("\nTag error following links.  f==d=%d",d);
	}
	else //if(Tag[t]==d) // 'if' Not needed; only remaining case
	{
		trf_message("\nTag error Tag[%d]=%d",t,d);
		exit(-2);
	}

	/* link in d */
	Distance[d].linkdown=f;
	h=Distance[f].linkup;
	Distance[d].linkup=h;
	Distance[f].linkup=d;
	if(h<=MAXDISTANCE) Distance[h].linkdown=d;
	Distance[d].linked=TRUE;
}


void untag_Distance_window(int d, int linkdown)
{

	int t;
	/* get next highest tag */
	t=d/TAGSEP;	/* was (int)ceil(d/TAGSEP): integer division, ceil is the identity */
	if(Tag[t]!=d) return;    /* Tag[t] is the largest index less or */
	/* equal to (t)x(TAGSEP) that is linked */
	else
	{
		while((t<=Toptag)&&(Tag[t]==d)) /* check higer tags and replace d */
			/* with its linkdown */
		{
			Tag[t]=linkdown;
			t++;
		}
	}
}

int no_matches_so_unlink_Distance(int d, int location, 
		struct distancelist *objptr)
{
	int *lo, *hi, g, h;
	int *z, *m;
	int windowleftend, windowsize;
	struct distancelist *list;

	windowsize=max(d,Min_Distance_Entries)+1; /* this value is used to */
	/* mod the index to the */
	/* entries.  The entries */
	/* run from 0 to */
	/* max(d,Min_Distance_Entries) */ 

	list=&(objptr[d]);
	hi=&(list->highindex);
	lo=&(list->lowindex);
	z=&(list->numentries);
	m=&(list->nummatches);
	if ((*z)!=0) /* there exist previous entries on this list */
	{
		/* first remove trailing entries, i.e., */
		/* those with location<location-d or */
		/* <location-Min_Distance_Window, whichever is further back */
		windowleftend=min(location-d,location-Min_Distance_Window)+1;
		/* note that preceding is based on Min_Distance_Window */
		/* rather than Min_Distance_Entries */

		while((list->entry[*lo].location<windowleftend)&&((*z)>0))
		{
			(*z)--;
			(*m)-=list->entry[*lo].size;
			(*lo)++;
			(*lo)=trf_wrap_index(*lo,windowsize);
			if((*lo)>windowsize)trf_message("\n    no_matches; *lo:%d windowsize:%d",*lo,windowsize);
		}
	}
	if((*z)==0) /* no more matches, so unlink */
	{
		/* given a distance d with zero matches, unlink from list */
		g=Distance[d].linkdown;
		h=Distance[d].linkup;

		Distance[g].linkup=h;
		if (h<=MAXDISTANCE) Distance[h].linkdown=g;
		Distance[d].linked=FALSE;

		/* once unlinked, we must also untag if it is a tag */
		untag_Distance_window(d,g);


		return(1);
	}
	return(0);  
}
