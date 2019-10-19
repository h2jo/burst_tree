// Author: Hang-Hyun Jo (h2jo23@gmail.com)
// Description: functions for some basic operations

// get distribution
void get_distr_long(long *sequence, long num_item, long *distr, long x_min, long x_max){
    long i, x;

    for(x = x_min; x <= x_max; x ++) distr[x] = 0;

    for(i = 0; i < num_item; i ++){
        x = sequence[i];
        distr[x] ++;
    }
}

void get_distr_double(double *sequence, long num_item, double resolution, long *distr, long x_min, long x_max){ // x is in [x_min, x_max]
    long i, x;

    for(x = x_min; x <= x_max; x ++) distr[x] = 0;

    for(i = 0; i < num_item; i ++){
        x = (long)(sequence[i] / resolution + 0.5);
        if(x < x_min || x > x_max) printf("out of range\n");
        distr[x] ++;
    }
}

long find_max(long *sequence, long num_item){
    long i, x, x_max;

    x_max = -1;
    for(i = 0; i < num_item; i ++){
        if(sequence[i] > x_max)
            x_max = sequence[i];
    }

    return x_max;
}

long find_min(long *sequence, long num_item){
    long i, x, x_min;

    x_min = 1e8;
    for(i = 0; i < num_item; i ++){
        if(sequence[i] < x_min)
            x_min = sequence[i];
    }

    return x_min;
}

// for event sequence analysis
long read_timings(char *folder, char *filename, long *timings, long n_max){
    long i, timing, timing0; 
    char output[500]; 
    FILE *file_in; 

    for(i = 0; i < n_max; i ++) timings[i] = 0;

    sprintf(output, "%stimings_%s.txt", folder, filename);
    file_in = fopen(output, "r"); 
    i = 0;
    timing0 = -1;
    while(fscanf(file_in, "%ld", &timing) && !feof(file_in)){
        if(timing == timing0) continue;
        timings[i] = timing;
        i ++;
        if(i == n_max) printf("data size exceeding the limit in read_timings()\n");
        timing0 = timing;
    }
    fclose(file_in); 

    return i;
}

// exchange one form to another
long get_iet_sequence_from_timings(long *timings, long n, long *iet_sequence){
    long i, timing0, timing, iet, iet_max;

    for(i = 0; i < n - 1; i ++) iet_sequence[i] = 0;

    iet_max = -1;
    timing0 = -1;
    for(i = 0; i < n; i ++){
        timing = timings[i];
        if(timing0 >= 0){
            iet = timing - timing0;
            iet_sequence[i - 1] = iet;
            if(iet > iet_max) iet_max = iet;
        }
        timing0 = timing;
    }

    return iet_max;
}

void get_timing_sequence_from_iet_sequence(long *iet_sequence, long num_iet, long *timing_sequence){
    long i, timing;

    for(i = 0; i <= num_iet; i ++) timing_sequence[i] = 0;

    timing = 0;
    for(i = 0; i < num_iet; i ++){
        timing += iet_sequence[i];
        timing_sequence[i + 1] = timing;
    }
}

void get_iet_sequence_from_tree_direct(long **tree, long *iet_sequence, long num_iet){
    long i, loc;

    for(i = 0; i < num_iet; i ++){
        loc = tree[i][0];
        iet_sequence[loc] = tree[i][1];
    }
}

void get_iet_sequence_from_tree(long **tree, long root_rank, long root_loc, long *iet_sequence, long num_iet){
    long rank_left, rank_right, train0, train1, loc;

    rank_left = tree[root_rank][2];
    rank_right = tree[root_rank][3];

    if(rank_left < 0 && rank_right < 0) return;

    if(rank_left >= 0){
        train1 = tree[rank_left][5]; // train size of right branch of "rank_left"
        loc = root_loc - train1;
        iet_sequence[loc] = tree[rank_left][1];
        get_iet_sequence_from_tree(tree, rank_left, loc, iet_sequence, num_iet);
    }
    if(rank_right >= 0){
        train0 = tree[rank_right][4]; // train size of left branch of "rank_right"
        loc = root_loc + train0;
        iet_sequence[loc] = tree[rank_right][1];
        get_iet_sequence_from_tree(tree, rank_right, loc, iet_sequence, num_iet);
    }
}

// iet distribution
void get_iet_distr(char *folder, char *filename, long *iet_sequence, long num_iet, double *logbin_params){
    long i, iet_min, iet_max, *iet_distr;
    char prefix[500];

    iet_min = 0;
    iet_max = find_max(iet_sequence, num_iet);

    iet_distr = vector_long(iet_min, iet_max);
    get_distr_long(iet_sequence, num_iet, iet_distr, iet_min, iet_max);

    // print
    sprintf(prefix, "iet_distr");
    print_distr_noTotal(folder, filename, prefix, iet_distr, iet_min, iet_max);

    // print logbinned distr
    get_logbin_period_noAvg(folder, filename, prefix, iet_distr, iet_max, logbin_params);

    free_vector_long(iet_distr, iet_min, iet_max);
}

// memory coeff of iet (but also for any sequence)
double calculate_memory_coeff(long *iet_sequence, long num_iet){ 
    long i, iet; 
    double iet_sum, iet2_sum, iet0, iet1, mean1, mean2, std1, std2, ietiet, memory;

    iet_sum = iet2_sum = ietiet = 0;
    for(i = 0; i < num_iet; i ++){
        iet = iet_sequence[i];
        iet_sum += iet;
        iet2_sum += iet * iet;
        if(i > 0) ietiet += iet * iet_sequence[i - 1];
    }
    iet0 = iet_sequence[0];
    iet1 = iet_sequence[num_iet - 1];

    mean1 = (iet_sum - iet1) / (double)(num_iet - 1);
    mean2 = (iet_sum - iet0) / (double)(num_iet - 1);
    std1 = sqrt((iet2_sum - iet1 * iet1) / (double)(num_iet - 1) - mean1 * mean1);
    std2 = sqrt((iet2_sum - iet0 * iet0) / (double)(num_iet - 1) - mean2 * mean2);
    ietiet /= (double)(num_iet - 1); //ietiet /= (double)(n - 2); <- wrong?

    memory = (ietiet - mean1 * mean2) / std1 / std2;

    return memory;
}

// memory coeff of iet (but also for any sequence)
void get_memory_iet(char *folder, char *filename, char *rrm, long *iet_sequence, long num_iet, long option){ 
    long i, iet; 
    double memory;
    char prefix[500], output[500];
    FILE *file_out;

    memory = calculate_memory_coeff(iet_sequence, num_iet);

    sprintf(prefix, "memory_iet");

    if(option == 0){ // original value
        sprintf(output, "%s%s_%s.txt", folder, prefix, filename);
        file_out = fopen(output, "w");
        fprintf(file_out, "%lf\n", memory);
        fclose(file_out);
    }
    else if(option == 1){ // ensemble
        sprintf(output, "%s%s_%s%s_ens.txt", folder, prefix, filename, rrm);
        file_out = fopen(output, "a"); // add not write
        fprintf(file_out, "%lf\n", memory);
        fclose(file_out);
    }
}

// autocorrelation function
void get_autocorrel(char *folder, char *filename, long *timings, long num_event, double *params){
    long i, j, lag, timing, timing_lag, timing1, index, auto_size;
    double **autocorrel, lag_min, lag_max, lag_step, lag_real, T, x_sum, x_avg, x_var, correl;
    char prefix[500];

    lag_min = params[0];
    lag_max = params[1];
    lag_step = params[2];
    auto_size = (long)(log(lag_max / lag_min) / log(lag_step) + 1);

    autocorrel = matrix_double(0, auto_size, 0, 1);
    for(i = 0; i < auto_size; i ++) autocorrel[i][0] = autocorrel[i][1] = 0;

    T = timings[num_event - 1] - timings[0];
    x_sum = num_event;
    x_avg = x_sum / T;
    x_var = x_avg - x_avg * x_avg;

    autocorrel[0][0] = 0;
    autocorrel[0][1] = 1;
    index = 1;
    lag_real = lag_min; 
	while(lag_real <= lag_max){
		lag = (long)(lag_real + 0.5); 
		//printf("lag %ld\n", lag); 

		correl = 0; 
		for(i = 0; i < num_event - 1; i ++){
            timing = timings[i];
            timing_lag = timing + lag;
            j = i + 1;
			while(j < num_event){
				timing1 = timings[j]; 
				if(timing1 == timing_lag){
					correl ++; 
					break; 
				}
				else if(timing1 > timing_lag)
					break; 
				j ++; 
			}
		}
        correl = (correl / T - x_avg * x_avg) / x_var;
		lag_real *= lag_step; 

        autocorrel[index][0] = lag;
        autocorrel[index][1] = correl;
        index ++;
	}

    sprintf(prefix, "autocorrel");
    print_double_matrix(folder, filename, prefix, autocorrel, 2, auto_size);

    free_matrix_double(autocorrel, 0, auto_size, 0, 1);
}

// train size distributions
long get_train_sequence_given_Dt(long *iet_sequence, long num_iet, long *train_sequence, long Dt){
    long i, iet, train, num_train, num_event;

    num_event = num_iet + 1;

    for(i = 0; i < num_event; i ++) train_sequence[i] = 0;

    train = 1;
    num_train = 0;
    for(i = 0; i < num_iet; i ++){
        iet = iet_sequence[i];
        if(iet <= Dt) train ++;
        else{
            train_sequence[num_train] = train;
            num_train ++;
            train = 1;
        }
    }

    return num_train;
}

void get_train_distr(char *folder, char *filename, long *iet_sequence, long num_iet, long *Dts, long num_Dt, double *logbin_params){
    long i, Dt, num_train, num_event, train_min, train_max;
    long *train_sequence, *train_distr;
    char prefix[500];

    num_event = num_iet + 1;
    
    train_sequence = vector_long(0, num_event);
    train_distr = vector_long(0, num_event);

    for(i = 0; i < num_Dt; i ++){
        Dt = Dts[i];
        printf("Dt=%ld\n", Dt);
        num_train = get_train_sequence_given_Dt(iet_sequence, num_iet, train_sequence, Dt);
        //train_min = find_min(train_sequence, num_train);
        train_min = 0;
        train_max = find_max(train_sequence, num_train);
        get_distr_long(train_sequence, num_train, train_distr, train_min, train_max);

        // print
        sprintf(prefix, "train_distr_Dt%ld", Dt);
        print_distr_noTotal(folder, filename, prefix, train_distr, train_min, train_max);
        
        // print logbinned distr
        get_logbin_period_noAvg(folder, filename, prefix, train_distr, train_max, logbin_params);
    }

    free_vector_long(train_sequence, 0, num_event);
    free_vector_long(train_distr, 0, num_event);
}

// get the statistics for the consecutive taus having the same value
void get_nonbinary_distr(char *folder, char *filename_ens, long *iet_sequence, long num_event){
    long i, j, iet, check, count, count_max, iet_max, iet_min, num_iet;
    long **nonbinary_distr, *count_iet;
    char prefix[500];
    
    num_iet = num_event - 1;
    count_iet = vector_long(0, num_iet); 
    nonbinary_distr = matrix_long(0, num_event, 0, 2);

    iet_max = find_max(iet_sequence, num_iet);
    iet_min = find_min(iet_sequence, num_iet);

    for(iet = iet_min; iet <= iet_max; iet ++){
        count_iet[iet] = 0; // initialize as 0
    }

    for(i = 0; i < num_event - 1; i ++){
        nonbinary_distr[i][0] = 0; // count
        nonbinary_distr[i][1] = num_event; // minimum of iet for each count
        nonbinary_distr[i][2] = 0; // maximum of iet for each count
    }

    count_max = 1;
    for(i = 0; i < num_event; i ++){
        iet = iet_sequence[i];
        count_iet[iet] ++; // increase count by 1

        for(j = iet_min; j < iet; j ++){ // for IETs less than iet
            count = count_iet[j];
            if(count >= 2){ // ~more than two consecutive bursts
                nonbinary_distr[count][0] ++;
                if(count > count_max) count_max = count;
                if(j < nonbinary_distr[count][1]) nonbinary_distr[count][1] = j;
                if(j > nonbinary_distr[count][2]) nonbinary_distr[count][2] = j;
            }
            count_iet[j] = 0; // reset the count
        }
    }

    for(j = iet_min; j <= iet_max; j ++){
        count = count_iet[j];
        if(count >= 2){
            nonbinary_distr[count][0] ++;
            if(count > count_max) count_max = count;
            if(j < nonbinary_distr[count][1]) nonbinary_distr[count][1] = j;
            if(j > nonbinary_distr[count][2]) nonbinary_distr[count][2] = j;
        }
    }

    sprintf(prefix, "nonbinary_distr");
    print_distr_matrix(folder, filename_ens, prefix, nonbinary_distr, 0, count_max);

    free_vector_long(count_iet, 0, num_iet);
    free_matrix_long(nonbinary_distr, 0, num_event, 0, 2);
}
