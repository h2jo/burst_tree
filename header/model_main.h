// alpha0 = -1/(alpha-1)
long get_powerlaw_alpha0(long xmin, double alpha0){
    double rand;

    rand = 1. - genrand_real2();
    return (long)(pow(rand, alpha0) * xmin + 0.5);
}

void generate_sorted_iet_sequence(long num_iet, double alpha, long *iet_sequence){
    long i, iet, iet_min, iet_max;
    long *iet_distr;
    double alpha0;

    for(i = 0; i < num_iet; i ++) iet_sequence[i] = 0;

    iet_min = 1;
    iet_max = 1e7;
    alpha0 = -1. / (alpha - 1.);
    for(i = 0; i < num_iet; i ++){
        do{
            iet = get_powerlaw_alpha0(iet_min, alpha0);
        }while(iet > iet_max);
        iet_sequence[i] = iet;
    }
    QuickSort3_descend_long1(iet_sequence, 0, num_iet - 1);
}

void update_tree(long **tree, long rank, long num_iet, long b1, long b2, long b_sum){
    long i, rank0, temp, j0, j1; 

    tree[rank][0] = b_sum;
    tree[rank][4] = b1; // left branch
    tree[rank][5] = b2; // right branch
    
    rank0 = -1;
    if(b1 > 1){
        for(i = rank + 1; i < num_iet; i ++){
            if(tree[i][1] == -1 && tree[i][0] == b1){
                tree[i][1] = 0; // if checked already
                tree[rank][2] = i; // rank of the left branch
                rank0 = i;
                break;
            }
        }
    }
    if(b2 > 1){
        for(i = rank + 1; i < num_iet; i ++){
            if(i == rank0) continue; // avoid double choice
            if(tree[i][1] == -1 && tree[i][0] == b2){
                tree[i][1] = 0; // if checked already
                tree[rank][3] = i; // rank of the right branch
                break;
            }
        }
    }
}

// using general kernel
double get_kernel_value(long b1, long b2, double *model_params){
    double a, h, s, logb1, logb2, temp0, temp1, model;

    if(b1 <= 0 || b2 <= 0) printf("NaN\n");

    a = model_params[1];
    h = model_params[2];
    s = model_params[3];
    model = (long)model_params[4];

    logb1 = log(b1); logb2 = log(b2);
    temp0 = (logb1 - logb2) * (logb1 - logb2);

    if(model == 0){
        temp1 = logb1 + logb2; // logsum
    }
    else if(model == 1){
        temp1 = (logb1 + logb2) * (logb1 + logb2); // logsumsq
    }
    else if(model == 2){ 
        temp1 = logb1 * logb1 + logb2 * logb2; // logsqsum
    }

    // multiply model
    return (1. + a * temp1) * (1. + h * exp(- temp0 / s));
}

void draw_bursts_using_kernel(long *burst_distr, long b_max, double sum, double *model_params, long *b1, long *b2){
    long bi, bj, ni, nj, b_sum;
    double rand, kernel_nn;
    
    rand = genrand_real2() * sum;
    for(bi = 1; bi <= b_max; bi ++){
        ni = burst_distr[bi];
        if(ni == 0) continue;

        if(ni >= 2){
            kernel_nn = get_kernel_value(bi, bi, model_params) * (double)(ni * (ni - 1)) * 0.5;
            if(rand < kernel_nn){
                *b1 = bi;
                *b2 = bi;
                return;
            }
            rand -= kernel_nn;
        }

        for(bj = bi + 1; bj <= b_max; bj ++){
            nj = burst_distr[bj];
            if(nj == 0) continue;

            kernel_nn = get_kernel_value(bi, bj, model_params) * (double)(ni * nj);
            if(rand < kernel_nn){
                *b1 = bi;
                *b2 = bj;
                return;
            }
            rand -= kernel_nn;
        }
    }
    //printf("wrong %ld %ld\n", *b1, *b2);
    return;
}

double get_kernel_sum(long *burst_distr, long b_max, double *model_params){
    long bi, bj, ni, nj;
    double sum, kernel_nn;

    sum = 0;
    for(bi = 1; bi <= b_max; bi ++){
        ni = burst_distr[bi];
        if(ni == 0) continue;

        if(ni >= 2){
            kernel_nn = get_kernel_value(bi, bi, model_params) * (double)(ni * (ni - 1)) * 0.5;
            sum += kernel_nn;
        }

        for(bj = bi + 1; bj <= b_max; bj ++){
            nj = burst_distr[bj];
            if(nj == 0) continue;

            kernel_nn = get_kernel_value(bi, bj, model_params) * (double)(ni * nj);
            sum += kernel_nn;
        }
    }

    return sum;
}

void merge_bursts_kernel(long **tree, long num_iet, long *burst_distr, double *model_params){
    long i, j, k, num_event, rank, b, b_sum, b1, b2, b_left, b_right, b_max;
    double sum;

    num_event = num_iet + 1;
    burst_distr[1] = num_event;
    for(i = 2; i <= num_event; i ++) burst_distr[i] = 0;

    b_max = 1;
    for(rank = num_iet - 1; rank >= 0; rank --){
        sum = get_kernel_sum(burst_distr, b_max, model_params); // = sum_{b,b'} K(b,b')
        draw_bursts_using_kernel(burst_distr, b_max, sum, model_params, &b1, &b2);

        burst_distr[b1] --; 
        burst_distr[b2] --; 
        b_sum = b1 + b2;
        burst_distr[b_sum] ++;
        if(b_sum > b_max) b_max = b_sum;

        // assign left and right
        if(b1 == b2){
            b_left = b_right = b1;
        }
        else if(genrand_real2() < 0.5){ 
            b_left = b1; b_right = b2;
        }
        else{
            b_left = b2; b_right = b1;
        }

        // update tree
        update_tree(tree, rank, num_iet, b_left, b_right, b_sum);
    }
}

void generate_tree(long **tree, long *iet_sequence, long num_event, double *model_params){
    long i, iet, root_rank, root_train, root_loc, num_iet;
    long *burst_distr;
    double alpha;

    num_iet = num_event - 1;
    burst_distr = vector_long(1, num_event);

    for(i = 0; i < num_iet; i ++){ // i can be a rank
        tree[i][0] = 0; // size of the root, and then location
        tree[i][1] = -1; // whether visited, and then iet
        tree[i][2] = -1; // rank of left branch
        tree[i][3] = -1; // rank of right branch
        tree[i][4] = 0; // size of left branch
        tree[i][5] = 0; // size of right branch
    }
    merge_bursts_kernel(tree, num_iet, burst_distr, model_params);

    alpha = model_params[0];
    generate_sorted_iet_sequence(num_iet, alpha, iet_sequence);

    for(i = 0; i < num_iet; i ++){ // i can be a rank
        tree[i][1] = iet_sequence[i]; // assign sorted iet
    }

    tree[0][0] = tree[0][4] - 1;
    update_tree_loc(tree, 0); // update loc

    free_vector_long(burst_distr, 1, num_event);
}

void main_model(long **tree, long *iet_sequence, long *timings, long num_event, double *model_params, char *folder, char *filename){
    long i, num_iet;
    
    num_iet = num_event - 1;

    generate_tree(tree, iet_sequence, num_event, model_params);
    print_tree(folder, filename, tree, num_iet);
    get_iet_sequence_from_tree_direct(tree, iet_sequence, num_iet);
    get_timing_sequence_from_iet_sequence(iet_sequence, num_iet, timings);
}


///////////////////////////////////// ensemble analysis
void summarize_ensemble_option(char *folder, char *filename, char *options, long ENS){
    long i, j, k, ens, num_iet;
    long train_Dts[6];
    char prefix[100];

    // train size distribution
    if(options[3] == '1'){
        train_Dts[0] = 16;
        train_Dts[1] = 64;
        train_Dts[2] = 256;
        train_Dts[3] = 1024;
        for(i = 0; i < 4; i ++){
            sprintf(prefix, "train_distr_Dt%ld", train_Dts[i]);
            printf("%s\n", prefix);
            summarize_distributions(folder, filename, prefix, ENS);
        }
    }

    // USING tree
    // memory coeff of trains, LC, suscept
    if(options[4] == '1'){
        sprintf(prefix, "percol");
        printf("%s\n", prefix);
        summarize_curves(folder, filename, prefix, ENS);
    }

    // memory coeff of b_v and b_w
    if(options[5] == '1'){
        sprintf(prefix, "memory_vw_Dt");
        printf("%s\n", prefix);
        summarize_curves(folder, filename, prefix, ENS);
    }
}

void main_model_ens(long **tree, long *iet_sequence, long *timings, long num_event, double *model_params, char *folder, char *filename, char *options, long ENS){
    long ens, num_iet;
    char filename_ens[500];
    
    num_iet = num_event - 1;

    for(ens = 0; ens < ENS; ens ++){
        sprintf(filename_ens, "%s_ens%ld", filename, ens);
        generate_tree(tree, iet_sequence, num_event, model_params);
        print_tree(folder, filename_ens, tree, num_iet);
        get_iet_sequence_from_tree_direct(tree, iet_sequence, num_iet);
        get_timing_sequence_from_iet_sequence(iet_sequence, num_iet, timings);

        burst_analysis_model_option(folder, filename, ens, timings, iet_sequence, tree, num_event, options);
    }

    summarize_ensemble_option(folder, filename, options, ENS);
}
