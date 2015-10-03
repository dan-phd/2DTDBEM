// Common functions

#pragma once

#include "Define.h"
#include <omp.h>

//Timing
#ifdef OS_WIN
#include <time.h>
void start_timing(clock_t& t)
{
	t = clock();
}
void finish_timing(clock_t& t)
{
	t = clock() - t;
	printf("\t - complete. \n\nThe elapsed time is %f seconds\n\n",
		((float)t) / CLOCKS_PER_SEC);
}
#else
#include <sys/time.h>	// this isn't in MSVC
void start_timing(struct timeval *start_time)
{
	gettimeofday(start_time, NULL);
}
void finish_timing(struct timeval *start_time)
{
	struct timeval end_time;
	gettimeofday(&end_time, NULL);
	printf("\n\nComplete. The elapsed time is %f seconds\n\n",
		end_time.tv_sec - start_time->tv_sec + (end_time.tv_usec - start_time->tv_usec) / 1e6);
}
#endif


//omp setup
void setup_omp()
{
	int nthreads, tid;
#pragma omp parallel private (tid)
	{
		tid = omp_get_thread_num();
		if (tid == 0)
		{
			nthreads = omp_get_num_threads();
			printf("\nTotal threads = %i\n\n", nthreads);
		}
	}
	omp_set_num_threads(nthreads);
}