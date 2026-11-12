#ifndef TRF_PERIODS_PERIOD_RANKING_H
#define TRF_PERIODS_PERIOD_RANKING_H

/* Return the same ranked periods as GetTopPeriods, or 1 on allocation failure.
 * As in the original routine, length is at least 2 and index maps each input
 * byte to 0..3. The caller supplies room for number_of_periods results.
 */
int trf_get_top_periods(const unsigned char *pattern, int length, int *toparray,
                        const int *index, int maxdistance,
                        int number_of_periods);

#endif
