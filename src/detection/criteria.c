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

/* Statistical criteria that decide whether a candidate distance is aligned and reported. */

#include "core/globals.h"
#include "core/report.h"
#include "detection/criteria.h"
#include "detection/distance_list.h"
#include "periods/period_ranking.h"

int GetTopPeriods(unsigned char *pattern, int length, int *toparray)
{
    return trf_get_top_periods(pattern, length, toparray, Index,
                               MAXDISTANCE, NUMBER_OF_PERIODS);
}

int multiples_criteria_4(int found_d)
{
	int g,lowerindex,upperindex;
	int topperiods[NUMBER_OF_PERIODS];
	unsigned char* pattern;
	int length;

	lowerindex=AlignPair.indexprime[AlignPair.length];
	upperindex=AlignPair.indexprime[1];
	pattern = Sequence+lowerindex;
	length = upperindex-lowerindex+1;

	/* size one is verified based on composition */
	/* only returns true if > 80% comp for any nucleotide */
	if(found_d==1)
	{
		int comps[4]={0};
		int total=0;
		int maxind=0;
		float percentmatch;

		for(g=0;g<length;g++)
		{
			comps[Index[pattern[g]]]++;
		}
		for(g=0;g<4;g++)
		{
			total+=comps[g];
			if(comps[g]>comps[maxind]) maxind=g;
		}
		percentmatch = comps[maxind]*100.0f/total;
		if(percentmatch>=80.0f) return TRUE; 
		else return FALSE;
	}


	if(GetTopPeriods(pattern, length, topperiods)) /* gettopperiods returns zero on success */
	{
		fprintf(stderr,"\nUnable to allocate counts array in GetTopPeriods()!!!");
		exit(-1);
	}

	/* copy into global for compatibility with best_period_list */
	/*
	   for(g=0;g<NUMBER_OF_PERIODS; g++)
	   {
	   Sortmultiples[g] = topperiods[g];
	   }
	 */

	/* modified 5/23/05 G. Benson */
	for(g=0;g<NUMBER_OF_PERIODS_INTO_SORTMULTIPLES; g++)
		Sortmultiples[g] = topperiods[g];

	/*
	   for(g=0;g<NUMBER_OF_PERIODS; g++)
	   {
	   if(found_d==topperiods[g]) return(TRUE);
	   }
	 */

	/* modified 5/23/05 G. Benson */
	for(g=0;g<NUMBER_OF_PERIODS_TO_TEST; g++)
		if(found_d==topperiods[g]) return(TRUE);

	return(FALSE);  
}


/********************* new_meet_criteria_3 ***************************/

int new_meet_criteria_3(int d, int location, int tuplesize)
{

	struct distancelist *main_d_info, *range_d_info;
	int min_krun_matches,max_first_match_location,main_d_matches,
	    range_d_min_for_waiting_time_test,main_d_first_match_location,
	    range_d_first_match_location,range_d_matches;
	int low_end_of_range, high_end_of_range,lopointer,hipointer,
	    d_range,waiting_time_ok,waiting_time_d;
	int t,m,s,d_still_best;

	/***********************************************************/



	/* collect info about d for criteria tests */
	main_d_info=&(Distance[d]);
	min_krun_matches=main_d_info->k_run_sums_criteria;
	max_first_match_location
		= max(0,location - max(d,Min_Distance_Window)) 
		+ main_d_info->waiting_time_criteria;

	low_end_of_range=main_d_info->lo_d_range;
	high_end_of_range=main_d_info->hi_d_range;

	/* get number of matches in d */
	main_d_matches=main_d_info->nummatches;

	/* set up a test for checking waiting time on range distances */
	/* criteria is that they must have at least 35% of the krun matches */
	range_d_min_for_waiting_time_test=(int)(0.35*min(min_krun_matches,main_d_matches));

	/* calculate the location of the first match for d */
	main_d_first_match_location 
		= main_d_info->entry[main_d_info->lowindex].location /* end of first run of matches */
		- main_d_info->entry[main_d_info->lowindex].size     /* minus size of run */
		+ tuplesize;                                         /* plus tuplesize */


	/* test waiting time criteria for d */
	waiting_time_ok=0;
	waiting_time_d=0;
	if(main_d_first_match_location<=max_first_match_location)
	{
		waiting_time_ok=1;
		waiting_time_d=d;
	}


	/* does d meet all criteria? */
	if((main_d_matches>=min_krun_matches)&&(waiting_time_ok))
	{	
		return(TRUE);
	}

	/* no, so now look in the low range to see if d has most matches */
	/* accumulate matches in low range */
	/* test for waiting time if main d didn't pass test */
	hipointer=d;
	lopointer=d;
	t=main_d_info->linkdown;
	m=main_d_matches;
	d_still_best=TRUE;
	while((t>=low_end_of_range) && (d_still_best))
	{
		range_d_info=&(Distance[t]);
		s=range_d_info->linkdown;
		if(!no_matches_so_unlink_Distance(t,location,Distance))
		{
			if((range_d_matches=range_d_info->nummatches)>main_d_matches) 
				d_still_best=FALSE;
			else
			{
				lopointer=t;
				m+=range_d_info->nummatches;
				if((!waiting_time_ok)&&(range_d_matches>=range_d_min_for_waiting_time_test))
				{
					/* calculate the location of the first match for range d */
					range_d_first_match_location 
						= range_d_info->entry[range_d_info->lowindex].location /* end of first run of matches */
						- range_d_info->entry[range_d_info->lowindex].size     /* minus size of run */
						+ tuplesize;                                           /* plus tuplesize */

					/* do waiting time test; if successful, store t*/
					if(range_d_first_match_location<=max_first_match_location)
					{
						waiting_time_ok=1;
						waiting_time_d=t;
					}
				}
			}
		}
		t=s;
	}

	/* stop if d not best */
	if(!d_still_best)
		return(FALSE);

	/* now check that d is best in upper range */
	t=main_d_info->linkup;
	while((t<=high_end_of_range)&&(d_still_best))
	{
		range_d_info=&(Distance[t]);
		s=range_d_info->linkup;
		if(!no_matches_so_unlink_Distance(t,location,Distance))
		{
			if(range_d_info->nummatches>main_d_matches) 
				d_still_best=FALSE;
		}
		t=s;
	}

	/* stop if d not best */
	if(!d_still_best)
		return(FALSE);

	/* now test lowest range */
	if((m>=min_krun_matches)&&(waiting_time_ok))
	{	
		return(TRUE);
	}

	/* lower range didn't work, now test higher ranges */

	d_range=d-low_end_of_range+1;
	t=main_d_info->linkup;
	while(t<=high_end_of_range)
	{
		range_d_info=&(Distance[t]);
		s=range_d_info->linkup;
		hipointer=t;
		range_d_matches=range_d_info->nummatches;
		m+=range_d_matches;
		while(lopointer<hipointer-d_range+1)
		{
			m-=Distance[lopointer].nummatches;
			if(m<=0){trf_message("\n*** error in meet crieria, m<=0"); exit(-16);}
			/* if lopointer made waiting time okay, clear out waiting time variables */
			if((waiting_time_ok)&&(lopointer==waiting_time_d))
			{
				waiting_time_ok=0;
				waiting_time_d=0;
			}
			lopointer=Distance[lopointer].linkup;
		}
		/* test for waiting time on range d no matter if waiting time is okay or not */
		if(range_d_matches>=range_d_min_for_waiting_time_test)
		{
			/* calculate the location of the first match for range d */
			range_d_first_match_location 
				= range_d_info->entry[range_d_info->lowindex].location /* end of first run of matches */
				- range_d_info->entry[range_d_info->lowindex].size     /* minus size of run */
				+ tuplesize;                                           /* plus tuplesize */

			/* do waiting time test; if successful, store t*/
			if(range_d_first_match_location<=max_first_match_location)
			{
				waiting_time_ok=1;
				waiting_time_d=t;
			}
		}

		/* now test range */
		if((m>=min_krun_matches)&&(waiting_time_ok))
		{	
			return(TRUE);
		}
		t=s;
	}

	return(FALSE);
}
