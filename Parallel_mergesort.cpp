/*
*  This file is part of Christian's OpenMP software lab 
*
*  Copyright (C) 2016 by Christian Terboven <terboven@itc.rwth-aachen.de>
*  Copyright (C) 2016 by Jonas Hahnfeld <hahnfeld@itc.rwth-aachen.de>
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with this program; if not, write to the Free Software
*  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
*
*/
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <sys/time.h>
#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <cstring>

/**
  * helper routine: check if array is sorted correctly
  */
bool isSorted(int ref[], int data[], const size_t size){
	std::sort(ref, ref + size);
	for (size_t idx = 0; idx < size; ++idx){
		if (ref[idx] != data[idx]) {
			return false;
		}
	}
	return true;
}


/**
  * sequential merge step (straight-forward implementation)
  */
// TODO: cut-off could also apply here (extra parameter?)
// TODO: optional: we can also break merge in two halves
void MsMergeSequential(int *out, int *in, long begin1, long end1, long begin2, long end2, long outBegin) {
	long left = begin1;
	long right = begin2;

	long idx = outBegin;

	while (left < end1 && right < end2) {
		if (in[left] <= in[right]) {
			out[idx] = in[left];
			left++;
		} else {
			out[idx] = in[right];
			right++;
		}
		idx++;
	}

	while (left < end1) {
		out[idx] = in[left];
		left++, idx++;
	}

	while (right < end2) {
		out[idx] = in[right];
		right++, idx++;
	}
}

/**
  * sequential MergeSort
  */
// TODO: remember one additional parameter (depth)
// TODO: recursive calls could be taskyfied
// TODO: task synchronization also is required
void MsSequential(int *array, int *tmp, bool inplace, long begin, long end) {
	if (begin < (end - 1)) {
		const long half = (begin + end) / 2;
		MsSequential(array, tmp, !inplace, begin, half);
		MsSequential(array, tmp, !inplace, half, end);
		if (inplace) {
			MsMergeSequential(array, tmp, begin, half, half, end, begin);
		} else {
			MsMergeSequential(tmp, array, begin, half, half, end, begin);
		}
	} else if (!inplace) {
		tmp[begin] = array[begin];
	}
}

void MsMergeParallel(int *out, int *in, long begin1, long end1, long begin2, long end2, long outBegin) {
    // Base case: if the problem size is smaller than a certain size, switch to sequential merge
    if ((end1 - begin1) + (end2 - begin2) <= 10000) {
        MsMergeSequential(out, in, begin1, end1, begin2, end2, outBegin);
        return;
    }

    long left = begin1;
    long right = begin2;
    long idx = outBegin;

    // Parallelize the merge process
    #pragma omp parallel sections
    {
        #pragma omp section
        {
            while (left < end1 && right < end2) {
                if (in[left] <= in[right]) {
                    out[idx++] = in[left++];
                } else {
                    out[idx++] = in[right++];
                }
            }
        }

        #pragma omp section
        {
            // Merge the remaining elements of the left half (if any)
            while (left < end1) {
                out[idx++] = in[left++];
            }
        }

        #pragma omp section
        {
            // Merge the remaining elements of the right half (if any)
            while (right < end2) {
                out[idx++] = in[right++];
            }
        }
    }
}

void MsParallel(int *array, int *tmp, bool inplace, long begin, long end, int depth, int cutoff) {
    if (begin < (end - 1)) {
        const long half = (begin + end) / 2;

        if ( depth < cutoff) {
            #pragma omp task
            {
                MsParallel(array, tmp, !inplace, begin, half, depth+1, cutoff);
            }
            #pragma omp task
            {
                MsParallel(array, tmp, !inplace, half, end, depth+1, cutoff);
            }
            #pragma omp taskwait            
        } else {
            // Sequential version for smaller problems
            MsSequential(array, tmp, !inplace, begin, half);
            MsSequential(array, tmp, !inplace, half, end);
        }
        if (inplace) {
            MsMergeSequential(array, tmp, begin, half, half, end, begin);
        } else {
            MsMergeSequential(tmp, array, begin, half, half, end, begin);
        }

		// Sequential merge is more useful for our case, as discussed in the report
        //MsMergeParallel(tmp, array, begin, half, half, end, begin);

    } else if (!inplace) {
        tmp[begin] = array[begin];
    }
}

/**
  * Serial MergeSort
  */
// TODO: this function should create the parallel region
// TODO: good point to compute a good depth level (cut-off)
void MsSerial(int *array, int *tmp, const size_t size, int cutoff) {
    
    // omp environment setting before the function call
    omp_set_dynamic(1);
    omp_set_num_threads(omp_get_max_threads());

    #pragma omp parallel
    {   
        #pragma omp single
        {
            MsParallel(array, tmp, true, 0, size, 0, cutoff);
        }
    }
}

/** 
  * @brief program entry point
  */
int main(int argc, char* argv[]) {

    struct timeval t1, t2;
    double etime_seq, etime_par;
    int cutoff = 4;  // cutoff value

    if (argc != 2) {
        printf("Usage: MergeSort.exe <array size> \n");
        printf("\n");
        return EXIT_FAILURE;
    } else {
        const size_t stSize = strtol(argv[1], NULL, 10);
        int *data = (int*) malloc(stSize * sizeof(int));
        int *tmp = (int*) malloc(stSize * sizeof(int));
        int *ref = (int*) malloc(stSize * sizeof(int));

        printf("Cutoff Value: %d\n", cutoff);
        printf("Initialization...\n");

        srand(95);
        for (size_t idx = 0; idx < stSize; ++idx) {
            data[idx] = (int)(stSize * (double(rand()) / RAND_MAX));
        }
        std::copy(data, data + stSize, ref);

        double dSize = (stSize * sizeof(int)) / 1024 / 1024;
        printf("Sorting %zu elements of type int (%f MiB)...\n", stSize, dSize);

        // SEQUENTIAL MERGE-SORT
        gettimeofday(&t1, NULL);
		MsSequential(data, tmp, true, 0, stSize);
        gettimeofday(&t2, NULL);

        etime_seq = (t2.tv_sec - t1.tv_sec) * 1000 + (t2.tv_usec - t1.tv_usec) / 1000;
        etime_seq = etime_seq / 1000;
        printf("Sequential MergeSort done, took %f sec.\nVerification...", etime_seq);

        if (isSorted(ref, data, stSize)) {
            printf(" successful.\n");
        } else {
            printf(" FAILED.\n");
        }

        // get again the first array
        std::copy(ref, ref + stSize, data);

        // PARALLEL MERGE SORT
        gettimeofday(&t1, NULL);
        MsSerial(data, tmp, stSize, cutoff);
        gettimeofday(&t2, NULL);

        etime_par = (t2.tv_sec - t1.tv_sec) * 1000 + (t2.tv_usec - t1.tv_usec) / 1000;
        etime_par = etime_par / 1000;
        printf("Parallel MergeSort done, took %f sec.\nVerification...", etime_par);

        if (isSorted(ref, data, stSize)) {
            printf(" successful.\n");
        } else {
            printf(" FAILED.\n");
        }

        // Speedup value printing
        double speedup = etime_seq / etime_par;
        printf("Speedup: %f\n", speedup);

        // free the allocated memory
        free(data);
        free(tmp);
        free(ref);
    }

    return EXIT_SUCCESS;
}