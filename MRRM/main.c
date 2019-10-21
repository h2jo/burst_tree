// Author: Hang-Hyun Jo (h2jo23@gmail.com)
// Description: main c file for MRRMs

#include "../header/head.h"

int main(int argv, char *argc[])
{
	time_t t_start, t_end, t_duration;
    long num_event_max, num_event, num_iet, iet_max, ens, binsize, isTreeRead, isTreeRead_ens, ENS;
    long *timings, *iet_sequence, **tree;
    double memory_resol, autocorrel_dts[3];
	char filename[100], folder[500];
    FILE *time_out;

	strcpy(filename, argc[1]); // file name
	strcpy(folder, argc[2]); // a FULL PATH to the file

    init_genrand(0);

    // record the computation time
	time_out = fopen("time_record.txt", "a"); time(&t_start);
	fprintf(time_out, "\nexe_filename=%s\nstart time=%s", filename, ctime(&t_start));
	printf("exe_filename=%s\nstart time=%s", filename, ctime(&t_start));

    // read timings from the data
    num_event_max = 2000000;
    timings = vector_long(0, num_event_max);
    num_event = read_timings(folder, filename, timings, num_event_max);

    // perform the MRRMs 
    num_iet = num_event - 1;
    iet_sequence = vector_long(0, num_iet - 1);
    tree = matrix_long(0, num_iet - 1, 0, 5);

    //iet_max = get_iet_sequence_from_timings(timings, num_event, iet_sequence);

    isTreeRead = 1; // 0: detect the tree from iet sequence
                    // 1: read the tree from the existing file
    get_ready(folder, filename, timings, iet_sequence, tree, num_event, isTreeRead);
    printf("get ready\n");

    ENS = 20; // ensemble size
    isTreeRead_ens = 0; // 0: detect the tree from iet sequence
                        // 1: read the tree from the existing file
    main_RRM(folder, filename, timings, iet_sequence, tree, num_event, ENS, isTreeRead_ens);
    printf("RRM\n");

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
