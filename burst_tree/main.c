// Author: Hang-Hyun Jo (h2jo23@gmail.com)
// Beginning date: 28 Jan 2018
// Description: main c file for burst-tree decomposition and analysis

#include "../header/head.h"

int main(int argv, char *argc[])
{
	time_t t_start, t_end, t_duration;
    long num_event_max, num_event, num_iet, iet_max, ens, binsize, isTreeRead, isTreeRead_ens, ENS;
    long *timings, *iet_sequence, **tree;
    double memory_resol, autocorrel_dts[3];
    char options[10], filename[100], folder[500];
    FILE *time_out;

	strcpy(filename, argc[1]); // file name
	strcpy(folder, argc[2]); // a FULL PATH to the file

    init_genrand(0); // initialize the random number generator

    // record the computation time
	time_out = fopen("time_record.txt", "a"); time(&t_start);
	fprintf(time_out, "\nexe_filename=%s\nstart time=%s", filename, ctime(&t_start));
	printf("exe_filename=%s\nstart time=%s", filename, ctime(&t_start));

    // read timings from the data
    num_event_max = 2000000;
    timings = vector_long(0, num_event_max);
    num_event = read_timings(folder, filename, timings, num_event_max);

    // perform the burst-tree analysis
    num_iet = num_event - 1;
    iet_sequence = vector_long(0, num_iet - 1);
    tree = matrix_long(0, num_iet - 1, 0, 5);

    isTreeRead = 0; // 0: detect the tree from iet sequence
                    // 1: read the tree from the existing file
    get_ready(folder, filename, timings, iet_sequence, tree, num_event, isTreeRead);
    printf("get ready\n");

    // options: 0 (don't) or 1 (do)
    options[0] = '0'; // autocorrelation function
    options[1] = '1'; // iet distribution
    options[2] = '0'; // memory coefficient between iets
    options[3] = '1'; // burst size distribution
    options[4] = '1'; // percolation, including memory coefficient between bursts
    options[5] = '1'; // memory coefficient between sibling bursts
    options[6] = '1'; // asymmetry
    options[7] = '1'; // kernel 2D or diagonal
    options[8] = '0'; // statistics of cases violating "binary" tree

    ens = -1; // no ens = single event sequence analysis
    burst_analysis_option(folder, filename, ens, timings, iet_sequence, tree, num_event, options);
    printf("burst analysis\n");

    // record the computation time
	time(&t_end); t_duration = t_end - t_start;
	fprintf(time_out, "  end time=%s  duration=%ld\n", ctime(&t_end), t_duration);
	printf("  end time=%s  duration=%ld\n", ctime(&t_end), t_duration);
	fclose(time_out);

    // free the memory allocation
    free_vector_long(timings, 0, num_event_max);
    free_vector_long(iet_sequence, 0, num_iet - 1);
    free_matrix_long(tree, 0, num_iet - 1, 0, 5);
}
