void get_train_distr_model(char *folder, char *filename, long *iet_sequence, long num_iet, long *Dts, long num_Dt, double *logbin_params){
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
        get_logbin_period(folder, filename, prefix, train_distr, train_max, logbin_params);
    }

    free_vector_long(train_sequence, 0, num_event);
    free_vector_long(train_distr, 0, num_event);
}

void burst_analysis_model_option(char *folder, char *filename, long ens, long *timings, long *iet_sequence, long **tree, long num_event, char *options){
    long i, j, k, num_iet, train_Dts[50], num_Dt;
    double autocorrel_params[3], logbin_params[3];
    char filename_ens[500];

    num_iet = num_event - 1;
    if(ens >= 0) sprintf(filename_ens, "%s_ens%ld", filename, ens);
    else sprintf(filename_ens, "%s", filename);

    // USING timings
    // autocorrel
    if(options[0] == '1'){
        printf("get autocorrel\n");
        autocorrel_params[0] = 1;
        autocorrel_params[1] = 1e6;
        autocorrel_params[2] = 1.2;
        get_autocorrel(folder, filename_ens, timings, num_event, autocorrel_params);
    }

    // USING iet_sequence
    // iet distribution
    if(options[1] == '1'){
        printf("get iet distr\n");
        logbin_params[0] = 2; // binStart
        logbin_params[1] = 0.2; // binSize
        logbin_params[2] = 10; // xc = crossover of x (?)
        get_iet_distr(folder, filename_ens, iet_sequence, num_iet, logbin_params);
    }

    // memory coeff of iet
    if(options[2] == '1'){
        printf("get memory iet\n");
        get_memory_iet(folder, filename, "", iet_sequence, num_iet, 1); // 1 = ensemble
    }

    // train size distribution
    if(options[3] == '1'){
        printf("get train distr\n");
        /*train_Dts[0] = 8;
        train_Dts[1] = 16;
        train_Dts[2] = 32;
        train_Dts[3] = 64;
        train_Dts[4] = 128;
        train_Dts[5] = 256;
        train_Dts[6] = 512;
        train_Dts[7] = 1024;
        train_Dts[8] = 2048;*/
        train_Dts[0] = 16;
        train_Dts[1] = 64;
        train_Dts[2] = 256;
        train_Dts[3] = 1024;
        logbin_params[0] = 1; // binStart
        logbin_params[1] = 0.2; // binSize
        logbin_params[2] = 10; // xc = crossover of x (?)
        num_Dt = 4;
        get_train_distr_model(folder, filename_ens, iet_sequence, num_iet, train_Dts, num_Dt, logbin_params);
    }

    // USING tree
    // memory coeff of trains, LC, suscept
    if(options[4] == '1'){
        printf("get percolation\n");
        logbin_params[0] = 5; logbin_params[1] = 0.18;
        get_percolation(folder, filename_ens, tree, num_iet, logbin_params);
    }

    // memory coeff of b_v and b_w
    if(options[5] == '1'){
        printf("get memory_vw\n");
        logbin_params[0] = 5; logbin_params[1] = 0.18;
        get_memory_vw_DtLog(folder, filename_ens, tree, num_iet, logbin_params);
    }

    // asymmertry of b_v and b_w
    if(options[6] == '1'){
        printf("get asymmetry\n");
        logbin_params[0] = 5; logbin_params[1] = 0.18;
        get_asymmetry_DtLog(folder, filename_ens, tree, num_iet, logbin_params);
    }

    // kernel diagonal
    if(options[7] == '1'){
        printf("get kernelDiag\n");
        logbin_params[0] = 5; logbin_params[1] = 0.18;
        get_kernelDiagLog(folder, filename_ens, tree, num_iet, logbin_params);
    }
}
