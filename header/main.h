// Author: Hang-Hyun Jo (h2jo23@gmail.com)
// Description: main functions

void get_ready(char *folder, char *filename, long *timings, long *iet_sequence, long **tree, long num_event, long isTreeRead){
    long num_iet;

    /* read timings
    timings = vector_long(0, num_event_max);
    num_event = read_timings(timings, num_event_max, filename, folder);
    */

    num_iet = num_event - 1;

    // get interevent time sequence
    get_iet_sequence_from_timings(timings, num_event, iet_sequence);

    // get tree
    if(isTreeRead == 0){ // detect tree after obtaining inter-event time sequence
        get_tree(iet_sequence, num_iet, tree);
        print_tree(folder, filename, tree, num_iet);
        printf("get and print tree done\n");
    }
    else{ // read tree from the file
        read_tree(folder, filename, tree, num_iet);
        printf("read tree done\n");
    }
}

void burst_analysis_option_rrm(char *folder, char *filename, char *rrm, long ens, long *timings, long *iet_sequence, long **tree, long num_event, char *options){
    long i, j, k, num_iet, train_Dts[6], num_Dt;
    double autocorrel_params[3], logbin_params[3];
    char filename_ens[500];

    num_iet = num_event - 1;
    sprintf(filename_ens, "%s%s_ens%ld", filename, rrm, ens);

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
        if(rrm[4] == '1' || rrm[4] == '2'){
            logbin_params[0] = 5; // binStart
            logbin_params[1] = 0.2; // binSize
            logbin_params[2] = 10; // xc = crossover of x (?)
        }
        else{
            logbin_params[0] = 2; // binStart
            logbin_params[1] = 0.2; // binSize
            logbin_params[2] = 10; // xc = crossover of x (?)
        }
        get_iet_distr(folder, filename_ens, iet_sequence, num_iet, logbin_params);
    }

    // memory coeff of iet
    if(options[2] == '1'){
        printf("get memory iet\n");
        get_memory_iet(folder, filename, rrm, iet_sequence, num_iet, 1); // 1 = ensemble
    }

    // train size distribution
    if(options[3] == '1'){
        printf("get train distr\n");
        train_Dts[0] = 60;    train_Dts[1] = 300;
        train_Dts[2] = 1800;  train_Dts[3] = 3600;
        train_Dts[4] = 10800; train_Dts[5] = 86400;
        if(rrm[4] == '1'){
            logbin_params[0] = 3; // binStart
            logbin_params[1] = 0.2; // binSize
            logbin_params[2] = 10; // xc = crossover of x (?)
        }
        else{
            logbin_params[0] = 1; // binStart
            logbin_params[1] = 0.2; // binSize
            logbin_params[2] = 10; // xc = crossover of x (?)
        }
        num_Dt = 3; // Wikipedia & JUNEC
        //num_Dt = 2; // Twitter
        get_train_distr(folder, filename_ens, iet_sequence, num_iet, train_Dts, num_Dt, logbin_params);
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

void burst_analysis_option(char *folder, char *filename, long ens, long *timings, long *iet_sequence, long **tree, long num_event, char *options){
    long i, j, k, num_iet, train_Dts[6], num_Dt;
    double autocorrel_params[3], logbin_params[3];
    char filename_ens[500];

    num_iet = num_event - 1;
    if(ens >= 0) sprintf(filename_ens, "%s_ens%ld", filename, ens);
    else sprintf(filename_ens, "%s", filename);

    // USING timings
    // autocorrelation function
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
        train_Dts[0] = 60;    train_Dts[1] = 300;
        train_Dts[2] = 1800;  train_Dts[3] = 3600;
        train_Dts[4] = 10800; train_Dts[5] = 86400;
        logbin_params[0] = 1; // binStart
        logbin_params[1] = 0.2; // binSize
        logbin_params[2] = 10; // xc = crossover of x (?)
        num_Dt = 3; // Wikipedia & JUNEC
        //num_Dt = 2; // Twitter
        //num_Dt = 6; // heartbeat
        get_train_distr(folder, filename_ens, iet_sequence, num_iet, train_Dts, num_Dt, logbin_params);
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

    // kernel 2D
    if(options[7] == '1'){
        printf("get kernel2D\n");
        logbin_params[0] = 5; logbin_params[1] = 0.18;
        get_kernel2DLog(folder, filename_ens, tree, num_iet, logbin_params);
    }

    /* kernel diagonal
    if(options[7] == '1'){
        printf("get kernelDiag\n");
        logbin_params[0] = 5; logbin_params[1] = 0.18;
        get_kernelDiagLog(folder, filename_ens, tree, num_iet, logbin_params);
    }*/

    // stats of nonbinary merging
    if(options[8] == '1'){
        get_nonbinary_distr(folder, filename_ens, iet_sequence, num_event);
    }
}

void summarize_ensemble_option(char *folder, char *filename, char *rrm, char *options, long ENS){
    long i, j, k, ens, num_iet;
    long train_Dts[6];
    double resolution;
    char prefix[100], filename_rrm[500];

    sprintf(filename_rrm, "%s%s", filename, rrm); 

    //num_iet = num_event - 1;

    // USING timings
    // autocorrel
    if(options[0] == '1'){
        sprintf(prefix, "autocorrel");
        printf("%s\n", prefix);
        summarize_curves(folder, filename_rrm, prefix, ENS, 0); // 0 for non-log
    }

    // USING iet_sequence
    // iet distribution
    if(options[1] == '1'){
        sprintf(prefix, "iet_distr");
        printf("%s\n", prefix);
        summarize_curves(folder, filename_rrm, prefix, ENS, 1); // 1 for log
    }

    // memory coeff of iet
    if(options[2] == '1'){
        sprintf(prefix, "memory_iet");
        printf("%s\n", prefix);
        //resolution = 0.005; // only for Twitter user1, RRM3
        resolution = 0.001;
        summarize_values(folder, filename, rrm, prefix, ENS, resolution);
    }

    // train size distribution
    if(options[3] == '1'){
        train_Dts[0] = 60;    train_Dts[1] = 300;
        train_Dts[2] = 1800;  train_Dts[3] = 3600;
        train_Dts[4] = 10800; train_Dts[5] = 86400;
        //for(i = 0; i < 3; i ++){
        for(i = 0; i < 2; i ++){ // for RRM
            sprintf(prefix, "train_distr_Dt%ld", train_Dts[i]);
            printf("%s\n", prefix);
            summarize_curves(folder, filename_rrm, prefix, ENS, 1); // 1 for log
        }
    }

    // USING tree
    // memory coeff of trains, LC, suscept
    if(options[4] == '1'){
        sprintf(prefix, "percol");
        printf("%s\n", prefix);
        summarize_curves(folder, filename_rrm, prefix, ENS, 0); // 0 for non-log
    }

    // memory coeff of b_v and b_w
    if(options[5] == '1'){
        sprintf(prefix, "memory_vw_Dt");
        printf("%s\n", prefix);
        summarize_curves(folder, filename_rrm, prefix, ENS, 1); // 1 for log
    }

    // asymmertry of b_v and b_w
    if(options[6] == '1'){
        sprintf(prefix, "asymmetry_Dt");
        printf("%s\n", prefix);
        summarize_curves(folder, filename_rrm, prefix, ENS, 1); // 1 for log
    }

    // kernel diagonal
    if(options[7] == '1'){
        sprintf(prefix, "kernelDiag");
        printf("%s\n", prefix);
        summarize_curves(folder, filename_rrm, prefix, ENS, 1); // 1 for log
    }
}

void RRM_timing(char *folder, char *filename, long *timings, long num_event, long ENS, long isTreeRead_ens){
    long i, ens, num_iet;
    long *timings_shuffle, *iet_sequence, **tree;
    char options[10], filename_ens[500], rrm[500];

    sprintf(rrm, "_rrm1");

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

    num_iet = num_event - 1;

    timings_shuffle = vector_long(0, num_event);
    iet_sequence = vector_long(0, num_iet - 1);
    tree = matrix_long(0, num_iet - 1, 0, 5);

    for(ens = 0; ens < ENS; ens ++){
        printf("ens = %ld\n", ens);

        shuffle_timing(timings, num_event, timings_shuffle);
        get_iet_sequence_from_timings(timings_shuffle, num_event, iet_sequence);

        // get tree
        sprintf(filename_ens, "%s%s_ens%ld", filename, rrm, ens);
        if(isTreeRead_ens == 0){ // detect tree after obtaining inter-event time sequence
            get_tree(iet_sequence, num_iet, tree);
            print_tree(folder, filename_ens, tree, num_iet);
        }
        else{ // read tree from the file
            read_tree(folder, filename_ens, tree, num_iet);
        }
        burst_analysis_option_rrm(folder, filename, rrm, ens, timings_shuffle, iet_sequence, tree, num_event, options);
    }
    summarize_ensemble_option(folder, filename, rrm, options, ENS);

    free_vector_long(timings_shuffle, 0, num_event);
    free_vector_long(iet_sequence, 0, num_iet - 1);
    free_matrix_long(tree, 0, num_iet - 1, 0, 5);
}

void RRM_iet(char *folder, char *filename, long *iet_sequence, long num_event, long ENS, long isTreeRead_ens){
    long i, ens, num_iet;
    long *timings, *iet_sequence_shuffle, **tree;
    char options[10], filename_ens[500], rrm[500];

    sprintf(rrm, "_rrm2");

    // options: 0 (don't) or 1 (do)
    options[0] = '0'; // autocorrelation function
    options[1] = '0'; // iet distribution
    options[2] = '0'; // memory coefficient between iets
    options[3] = '1'; // burst size distribution
    options[4] = '1'; // percolation, including memory coefficient between bursts
    options[5] = '1'; // memory coefficient between sibling bursts
    options[6] = '1'; // asymmetry
    options[7] = '1'; // kernel 2D or diagonal
    options[8] = '0'; // statistics of cases violating "binary" tree

    num_iet = num_event - 1;

    iet_sequence_shuffle = vector_long(0, num_iet - 1);
    timings = vector_long(0, num_event - 1);
    tree = matrix_long(0, num_iet - 1, 0, 5);

    for(ens = 0; ens < ENS; ens ++){
        printf("ens = %ld\n", ens);

        shuffle_iet(iet_sequence, num_event, iet_sequence_shuffle);
        get_timing_sequence_from_iet_sequence(iet_sequence_shuffle, num_iet, timings); 

        // get tree
        sprintf(filename_ens, "%s%s_ens%ld", filename, rrm, ens);
        if(isTreeRead_ens == 0){ // detect tree after obtaining inter-event time sequence
            get_tree(iet_sequence_shuffle, num_iet, tree);
            print_tree(folder, filename_ens, tree, num_iet);
        }
        else{ // read tree from the file
            read_tree(folder, filename_ens, tree, num_iet);
        }
        burst_analysis_option_rrm(folder, filename, rrm, ens, timings, iet_sequence_shuffle, tree, num_event, options);
    }
    summarize_ensemble_option(folder, filename, rrm, options, ENS);

    free_vector_long(timings, 0, num_event);
    free_vector_long(iet_sequence_shuffle, 0, num_iet - 1);
    free_matrix_long(tree, 0, num_iet - 1, 0, 5);
}

void RRM_leftright(char *folder, char *filename, long **tree, long num_event, long ENS, long isTreeRead_ens){
    long i, j, ens, num_iet;
    long *timings, *iet_sequence, **tree_shuffle;
    char options[10], filename_ens[500], rrm[500];

    sprintf(rrm, "_rrm3");

    // options: 0 (don't) or 1 (do)
    options[0] = '0'; // autocorrelation function
    options[1] = '0'; // iet distribution
    options[2] = '0'; // memory coefficient between iets
    options[3] = '0'; // burst size distribution
    options[4] = '1'; // percolation, including memory coefficient between bursts
    options[5] = '1'; // memory coefficient between sibling bursts
    options[6] = '1'; // asymmetry
    options[7] = '1'; // kernel 2D or diagonal
    options[8] = '0'; // statistics of cases violating "binary" tree

    num_iet = num_event - 1;

    timings = vector_long(0, num_event);
    iet_sequence = vector_long(0, num_iet - 1);
    tree_shuffle = matrix_long(0, num_iet - 1, 0, 5);

    for(ens = 0; ens < ENS; ens ++){
        printf("leftright ens = %ld\n", ens);

        sprintf(filename_ens, "%s%s_ens%ld", filename, rrm, ens);
        for(i = 0; i < num_iet; i ++) for(j = 0; j <= 5; j ++)
            tree_shuffle[i][j] = tree[i][j]; 
        
        shuffle_tree_leftright_root(tree_shuffle);
        update_tree_loc(tree_shuffle, 0);
        get_iet_sequence_from_tree_direct(tree_shuffle, iet_sequence, num_iet);
        get_timing_sequence_from_iet_sequence(iet_sequence, num_iet, timings); 

        burst_analysis_option_rrm(folder, filename, rrm, ens, timings, iet_sequence, tree_shuffle, num_event, options);
    }
    summarize_ensemble_option(folder, filename, rrm, options, ENS);

    free_vector_long(timings, 0, num_event);
    free_vector_long(iet_sequence, 0, num_iet - 1);
    free_matrix_long(tree_shuffle, 0, num_iet - 1, 0, 5);
}

void RRM_leftright_less(char *folder, char *filename, long **tree, long num_event, long ENS, long isTreeRead_ens){
    long i, j, ens, num_iet, dt;
    long *timings, *iet_sequence, **tree_shuffle;
    //char options[10] = "10101010";
    char options[10] = "00001110";
    char filename_ens[500], rrm[500];

    sprintf(rrm, "_rrm4");
    dt = 86400;

    // options: 0 (don't) or 1 (do)
    options[0] = '0'; // autocorrelation function
    options[1] = '0'; // iet distribution
    options[2] = '0'; // memory coefficient between iets
    options[3] = '0'; // burst size distribution
    options[4] = '1'; // percolation, including memory coefficient between bursts
    options[5] = '1'; // memory coefficient between sibling bursts
    options[6] = '1'; // asymmetry
    options[7] = '1'; // kernel 2D or diagonal
    options[8] = '0'; // statistics of cases violating "binary" tree

    num_iet = num_event - 1;

    timings = vector_long(0, num_event);
    iet_sequence = vector_long(0, num_iet - 1);
    tree_shuffle = matrix_long(0, num_iet - 1, 0, 5);

    for(ens = 0; ens < ENS; ens ++){
        printf("leftright ens = %ld\n", ens);

        sprintf(filename_ens, "%s%s_ens%ld", filename, rrm, ens);
        for(i = 0; i < num_iet; i ++) for(j = 0; j <= 5; j ++)
            tree_shuffle[i][j] = tree[i][j]; 
        
        shuffle_tree_leftright_root_less(tree_shuffle, dt);
        update_tree_loc(tree_shuffle, 0);
        get_iet_sequence_from_tree_direct(tree_shuffle, iet_sequence, num_iet);
        get_timing_sequence_from_iet_sequence(iet_sequence, num_iet, timings); 

        burst_analysis_option_rrm(folder, filename, rrm, ens, timings, iet_sequence, tree_shuffle, num_event, options);
    }
    summarize_ensemble_option(folder, filename, rrm, options, ENS);

    free_vector_long(timings, 0, num_event);
    free_vector_long(iet_sequence, 0, num_iet - 1);
    free_matrix_long(tree_shuffle, 0, num_iet - 1, 0, 5);
}

void RRM_leftright_more(char *folder, char *filename, long **tree, long num_event, long ENS, long isTreeRead_ens){
    long i, j, ens, num_iet, dt;
    long *timings, *iet_sequence, **tree_shuffle;
    //char options[10] = "10101010";
    char options[10] = "00001110";
    char filename_ens[500], rrm[500];

    sprintf(rrm, "_rrm5");
    dt = 86400;

    // options: 0 (don't) or 1 (do)
    options[0] = '0'; // autocorrelation function
    options[1] = '0'; // iet distribution
    options[2] = '0'; // memory coefficient between iets
    options[3] = '0'; // burst size distribution
    options[4] = '1'; // percolation, including memory coefficient between bursts
    options[5] = '1'; // memory coefficient between sibling bursts
    options[6] = '1'; // asymmetry
    options[7] = '1'; // kernel 2D or diagonal
    options[8] = '0'; // statistics of cases violating "binary" tree

    num_iet = num_event - 1;

    timings = vector_long(0, num_event);
    iet_sequence = vector_long(0, num_iet - 1);
    tree_shuffle = matrix_long(0, num_iet - 1, 0, 5);

    for(ens = 0; ens < ENS; ens ++){
        printf("leftright ens = %ld\n", ens);

        sprintf(filename_ens, "%s%s_ens%ld", filename, rrm, ens);
        for(i = 0; i < num_iet; i ++) for(j = 0; j <= 5; j ++)
            tree_shuffle[i][j] = tree[i][j]; 
        
        shuffle_tree_leftright_root_more(tree_shuffle, dt);
        update_tree_loc(tree_shuffle, 0);
        get_iet_sequence_from_tree_direct(tree_shuffle, iet_sequence, num_iet);
        get_timing_sequence_from_iet_sequence(iet_sequence, num_iet, timings); 

        burst_analysis_option_rrm(folder, filename, rrm, ens, timings, iet_sequence, tree_shuffle, num_event, options);
    }
    summarize_ensemble_option(folder, filename, rrm, options, ENS);

    free_vector_long(timings, 0, num_event);
    free_vector_long(iet_sequence, 0, num_iet - 1);
    free_matrix_long(tree_shuffle, 0, num_iet - 1, 0, 5);
}

void main_RRM(char *folder, char *filename, long *timings, long *iet_sequence, long **tree, long num_event, long ENS, long isTreeRead_ens){
    RRM_timing(folder, filename, timings, num_event, ENS, isTreeRead_ens); printf("RRM1\n");
    RRM_iet(folder, filename, iet_sequence, num_event, ENS, isTreeRead_ens); printf("RRM2\n");
    RRM_leftright(folder, filename, tree, num_event, ENS, isTreeRead_ens); printf("RRM3\n");
    //RRM_leftright_less(folder, filename, tree, num_event, ENS, isTreeRead_ens); printf("RRM4\n");
    //RRM_leftright_more(folder, filename, tree, num_event, ENS, isTreeRead_ens); printf("RRM5\n");
}
