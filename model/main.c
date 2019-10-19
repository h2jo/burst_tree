// Author: Hang-Hyun Jo (h2jo23@gmail.com)
// Beginning date: 28 Jan 2018
// Description: main c file for kernel-based modeling

#include "../header/head_model.h"

int main(int argv, char *argc[])
{
	time_t t_start, t_end, t_duration;
    long num_event, num_iet, *timings, *iet_sequence, **tree, ENS;
    double model_params[10], alpha, param_a, param_h, param_s;
    char options[10], filename[100], folder[500];
    FILE *time_out;

	strcpy(folder, argc[1]); // a FULL PATH to the file
    num_event = atoi(argc[2]);
    alpha = atof(argc[3]);
    param_a = atof(argc[4]);
    param_h = atof(argc[5]);
    param_s = atof(argc[6]);
    ENS = atoi(argc[7]);
    
    init_genrand((long)(num_event * alpha + ENS + param_a + param_h + param_s));

    sprintf(filename, "n%ld_alpha%.1lf_a%.0lf_h%.0lf_s%.0lf", num_event, alpha, param_a, param_h, param_s); 

    // record the computation time
	time_out = fopen("time_record.txt", "a"); time(&t_start);
	fprintf(time_out, "\nexe_filename=%s\nstart time=%s", filename, ctime(&t_start));
	printf("exe_filename=%s\nstart time=%s", filename, ctime(&t_start));

    // options: 0 (don't) or 1 (do)
    options[0] = '0'; // autocorrelation function
    options[1] = '0'; // iet distribution
    options[2] = '0'; // memory coefficient between iets
    options[3] = '1'; // burst size distribution
    options[4] = '1'; // percolation, including memory coefficient between bursts
    options[5] = '1'; // memory coefficient between sibling bursts
    options[6] = '0'; // asymmetry
    options[7] = '0'; // kernel 2D or diagonal
    options[8] = '0'; // statistics of cases violating "binary" tree

    model_params[0] = alpha; // power-law exponent for IET distribution
    model_params[1] = param_a; // =c1 in Eq. 5
    model_params[2] = param_h; // =c2 in Eq. 5
    model_params[3] = param_s; // =c3 in Eq. 5
    model_params[4] = 0; // 0 for sum model (Eq. 5 in the main text)
                         // 1 for sum and squared model (Eq. S3 in SI)
                         // 2 for squared and sum model (Eq. S4 in SI)

    num_iet = num_event - 1;
    timings = vector_long(0, num_event - 1);
    iet_sequence = vector_long(0, num_iet - 1);
    tree = matrix_long(0, num_iet - 1, 0, 5);

    // generate and analyze ENS event sequences using model kernel
    main_model_ens(tree, iet_sequence, timings, num_event, model_params, folder, filename, options, ENS);

    //burst_analysis_model_option(folder, filename, -1, timings, iet_sequence, tree, num_event, options);
    //printf("burst analysis\n");

    // record the computation time
	time(&t_end); t_duration = t_end - t_start;
	fprintf(time_out, "  end time=%s  duration=%ld\n", ctime(&t_end), t_duration);
	printf("  end time=%s  duration=%ld\n", ctime(&t_end), t_duration);
	fclose(time_out);

    // free the memory allocation
    free_vector_long(timings, 0, num_event - 1);
    free_vector_long(iet_sequence, 0, num_iet - 1);
    free_matrix_long(tree, 0, num_iet - 1, 0, 5);
}
