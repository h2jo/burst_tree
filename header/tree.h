// Author: Hang-Hyun Jo (h2jo23@gmail.com)
// Description: functions for detecting and processing tree structure

/************ DETECTING TREE *****************/

void search_branch(long **tree, long root_rank, long list_size){
    long root_loc, bound_left, bound_right, k, k_loc;

    if(root_rank == list_size - 1)
        return;

    root_loc = tree[root_rank][0];
    bound_left = tree[root_rank][4];
    bound_right = tree[root_rank][5];

    // find the "left" branch
    if(root_loc > 0){
        for(k = root_rank + 1; k < list_size; k ++){
            k_loc = tree[k][0];
            if(k_loc < root_loc && k_loc >= bound_left){
                tree[root_rank][2] = k; // rank of left branch
                tree[k][4] = bound_left; // inherit bound_left from the root
                tree[k][5] = root_loc - 1; // new bound_right
                break;
            }
        }
    }

    // find the "right" branch
    if(root_loc < list_size - 1){
        for(k = root_rank + 1; k < list_size; k ++){
            k_loc = tree[k][0];
            if(k_loc > root_loc && k_loc <= bound_right){
                tree[root_rank][3] = k; // rank of right branch
                tree[k][4] = root_loc + 1; // new bound_left
                tree[k][5] = bound_right; // inherit bound_right from the root
                break;
            }
        }
    }
}

void merge_trains(long **tree, long root_rank){
    long rank_left, rank_right, train;

    rank_left = tree[root_rank][2];
    rank_right = tree[root_rank][3];

    if(rank_left == -1){
        if(rank_right == -1){ // no branches
            tree[root_rank][4] = 1;
            tree[root_rank][5] = 1;
        }
        else{ // only right branch
            tree[root_rank][4] = 1;
            tree[root_rank][5] = tree[rank_right][4] + tree[rank_right][5];
        }
    }
    else{
        if(rank_right == -1){  // only left branch
            tree[root_rank][4] = tree[rank_left][4] + tree[rank_left][5];
            tree[root_rank][5] = 1;
        }
        else{ // both branches
            tree[root_rank][4] = tree[rank_left][4] + tree[rank_left][5];
            tree[root_rank][5] = tree[rank_right][4] + tree[rank_right][5];
        }
    }
}

void get_tree(long *iet_sequence, long num_iet, long **tree){
    long i, iet, loc;
    double **iet_sequence_sort;

    iet_sequence_sort = matrix_double(0, num_iet - 1, 0, 1);

    for(i = 0; i < num_iet; i ++){ // rank of iet
        iet_sequence_sort[i][0] = i;
        iet_sequence_sort[i][1] = iet_sequence[i] + genrand_real2() * 0.1;
        // <-- avoid the degeneracy
        tree[i][0] = -1; // loc (=original index)
        tree[i][1] = -1; // iet
        tree[i][2] = -1; // rank of left branch
        tree[i][3] = -1; // rank of right branch
        tree[i][4] = 0; // left bound
        tree[i][5] = num_iet - 1; // right bound
    }

    QuickSort3_descend_double(iet_sequence_sort, 1, 1, 0, num_iet - 1);

    for(i = 0; i < num_iet; i ++){
        tree[i][0] = (long)iet_sequence_sort[i][0]; // loc
        tree[i][1] = (long)iet_sequence_sort[i][1]; // iet
    }

    for(i = 0; i < num_iet; i ++){
        search_branch(tree, i, num_iet);
        if(i % 10000 == 0) printf("searching %ldth node\n", i);
    }
    printf("search done\n");

    for(i = 0; i < num_iet; i ++){ // rank
        tree[i][4] = 0; // train size of left branch of node i
        tree[i][5] = 0; // train size of right branch of node i
    }

    for(i = num_iet - 1; i >= 0; i --) merge_trains(tree, i);
    printf("merge done\n");

    free_matrix_double(iet_sequence_sort, 0, num_iet - 1, 0, 1);
}

void print_tree(char *folder, char *filename, long **tree, long num_iet){
    long i, j;
    char output[500];
    FILE *tree_out;

    sprintf(output, "%stree_%s.txt", folder, filename);
    tree_out = fopen(output, "w");

    for(i = 0; i < num_iet; i ++){
        fprintf(tree_out, "%ld", i);
        for(j = 0; j < 6; j ++){
            fprintf(tree_out, " %ld", tree[i][j]);
        }
        fprintf(tree_out, "\n");
    }
    fclose(tree_out);
}

void read_tree(char *folder, char *filename, long **tree, long num_iet){
    long i, j, x;
    char output[500], temp[500];
    FILE *tree_in;

    sprintf(output, "%stree_%s.txt", folder, filename);
    tree_in = fopen(output, "r");

    for(i = 0; i < num_iet; i ++) for(j = 0; j < 6; j ++) tree[i][j] = 0;

    for(i = 0; i < num_iet; i ++){
        fscanf(tree_in , "%s", temp);
        for(j = 0; j < 6; j ++){
            fscanf(tree_in , "%ld", &x);
            tree[i][j] = x;
        }
    }

    fclose(tree_in);
}

/************ QUANTIFYING TREE: asymmetry index, kernel, etc.  *********/

void get_memory_vw_DtLog(char *folder, char *filename, long **tree, long num_iet, double *logbin_params){ 
    long i, j, *iet_sequence, iet_max, num_event, iet, index, binStart;
    double **memory_vw_Dt, memory_vw, num, bv, bw, bv_avg, bw_avg, bv_std, bw_std, bb_avg, bv_sum, bw_sum, bv2_sum, bw2_sum, bb_sum, iet_bin;
    double binSize, x_start, x_end, x_factor, x_factor_2;
    char output1[500], output2[500], prefix[500];
    FILE *file1_out;

    strcpy(prefix, "memory_vw_Dt");
    sprintf(output1, "%s%s_%sLog.txt", folder, prefix, filename); 
    file1_out = fopen(output1, "w");

    num_event = num_iet + 1;
    iet_max = tree[0][1];
    //iet_max = 0; for(i = 0; i < num_iet; i ++){ iet = tree[i][1]; if(iet > iet_max) iet_max = iet; }

    binStart = logbin_params[0];
    binSize = logbin_params[1];

    memory_vw_Dt = matrix_double(0, iet_max, 0, 5);
    for(i = 0; i <= iet_max; i ++) for(j = 0; j <= 5; j ++) memory_vw_Dt[i][j] = 0;

    for(i = 0; i < num_iet; i ++){
        iet = tree[i][1];
        bv = tree[i][4]; // burst size of the left-child
        bw = tree[i][5]; // burst size of the right-child

        memory_vw_Dt[iet][0] ++; // num of data
        memory_vw_Dt[iet][1] += bv; 
        memory_vw_Dt[iet][2] += bw; 
        memory_vw_Dt[iet][3] += bv * bv; 
        memory_vw_Dt[iet][4] += bw * bw; 
        memory_vw_Dt[iet][5] += bv * bw; 
    }

    index = 0;
    iet = 0;
    while(index < binStart){ // no binning
        num = memory_vw_Dt[iet][0];
        if(num >= 2){ // only when the std exists
            bv_avg = memory_vw_Dt[iet][1] / num;
            bw_avg = memory_vw_Dt[iet][2] / num;
            bv_std = sqrt(memory_vw_Dt[iet][3] / num - bv_avg * bv_avg);
            bw_std = sqrt(memory_vw_Dt[iet][4] / num - bw_avg * bw_avg);
            bb_avg = memory_vw_Dt[iet][5] / num;
            if(bv_std > 0 && bw_std > 0)
                memory_vw = (bb_avg - bv_avg * bw_avg) / bv_std / bw_std;
            else memory_vw = 0;
            //fprintf(file1_out, "%ld %.10lf %ld\n", iet, memory_vw, (long)num);
        }
        else{
            memory_vw = 0;
        }
        fprintf(file1_out, "%ld %.10lf\n", iet, memory_vw); // for RRM
        index ++;
        iet ++;
    }

    // log-binning
    x_factor = exp(binSize);
    x_factor_2 = exp(-binSize * 0.5);
    x_end = iet * x_factor;
    num = bv_sum = bv2_sum = bw_sum = bw2_sum = bb_sum = 0;
    while(iet <= iet_max){
        if(iet < x_end){
            num += memory_vw_Dt[iet][0];
            bv_sum += memory_vw_Dt[iet][1];
            bw_sum += memory_vw_Dt[iet][2];
            bv2_sum += memory_vw_Dt[iet][3];
            bw2_sum += memory_vw_Dt[iet][4];
            bb_sum += memory_vw_Dt[iet][5];
            iet ++;
        }
        else{
            iet_bin = x_end * x_factor_2;
            if(num >= 2){ // only when the std exists
                bv_avg = bv_sum / num;
                bw_avg = bw_sum / num;
                bv_std = sqrt(bv2_sum / num - bv_avg * bv_avg);
                bw_std = sqrt(bw2_sum / num - bw_avg * bw_avg);
                bb_avg = bb_sum / num;
                if(bv_std > 0 && bw_std > 0)
                    memory_vw = (bb_avg - bv_avg * bw_avg) / bv_std / bw_std;
                else memory_vw = 0;
                //fprintf(file1_out, "%.10lf %.10lf\n", iet_bin, memory_vw);
            }
            else{
                memory_vw = 0;
            }
            fprintf(file1_out, "%.10lf %.10lf\n", iet_bin, memory_vw);

            x_end *= x_factor;
            num = bv_sum = bw_sum = bv2_sum = bw2_sum = bb_sum = 0;
        }
    }
    fclose(file1_out);

    free_matrix_double(memory_vw_Dt, 0, iet_max, 0, 5);
}

void get_asymmetry_DtLog(char *folder, char *filename, long **tree, long num_iet, double *logbin_params){ 
    long i, j, *iet_sequence, iet_max, num_event, iet, index, binStart;
    double **asymm_Dt, asymm, num, bv, bw, asymm_avg, asymm_std, asymm_sum, asymm2_sum, iet_bin;
    double binSize, x_start, x_end, x_factor, x_factor_2;
    char output1[500], output2[500], prefix[500];
    FILE *file1_out;

    strcpy(prefix, "asymmetry_Dt");
    sprintf(output1, "%s%s_%sLog.txt", folder, prefix, filename); 
    file1_out = fopen(output1, "w");

    num_event = num_iet + 1;
    iet_max = tree[0][1];
    //iet_max = 0; for(i = 0; i < num_iet; i ++){ iet = tree[i][1]; if(iet > iet_max) iet_max = iet; }

    binStart = logbin_params[0];
    binSize = logbin_params[1];

    asymm_Dt = matrix_double(0, iet_max, 0, 2);
    for(i = 0; i <= iet_max; i ++) for(j = 0; j <= 2; j ++) asymm_Dt[i][j] = 0;

    for(i = 0; i < num_iet; i ++){
        iet = tree[i][1];
        bv = tree[i][4]; // burst size of the left-child
        bw = tree[i][5]; // burst size of the right-child

        asymm = (bw - bv) / (bw + bv); // right burst size - left burst size
        asymm_Dt[iet][0] ++;
        asymm_Dt[iet][1] += asymm;
        asymm_Dt[iet][2] += asymm * asymm;
    }

    index = 0;
    iet = 0;
    while(index < binStart){ // no binning
        num = asymm_Dt[iet][0];
        if(num == 1){ // std = 0 if the number of samples is 1
            asymm_avg = asymm_Dt[iet][1];
            asymm_std = 0;
            //fprintf(file1_out, "%ld %.10lf %.10lf %ld\n", iet, asymm_avg, asymm_std, (long)num);
            index ++;
        }
        else if(num >= 2){ // only when the std exists
            asymm_avg = asymm_Dt[iet][1] / num;
            asymm_std = sqrt(asymm_Dt[iet][2] / num - asymm_avg * asymm_avg);
            //fprintf(file1_out, "%ld %.10lf %.10lf %ld\n", iet, asymm_avg, asymm_std, (long)num);
            index ++;
        }
        else{
            asymm_avg = 0;
        }
        fprintf(file1_out, "%ld %.10lf\n", iet, asymm_avg); // for RRM
        iet ++;
    }

    // log-binning
    x_factor = exp(binSize);
    x_factor_2 = exp(-binSize * 0.5);
    x_end = iet * x_factor;
    num = asymm_sum = asymm2_sum = 0;
    while(iet <= iet_max){
        if(iet < x_end){
            num += asymm_Dt[iet][0];
            asymm_sum += asymm_Dt[iet][1];
            asymm2_sum += asymm_Dt[iet][2];
            iet ++;
        }
        else{
            if(num == 1){ // std = 0 if the number of samples is 1
                asymm_avg = asymm_Dt[iet][1];
                asymm_std = 0;
                //iet_bin = x_end * x_factor_2;
                //fprintf(file1_out, "%ld %.10lf %.10lf %ld\n", iet, asymm_avg, asymm_std, (long)num);
            }
            else if(num >= 2){ // only when the std exists
                asymm_avg = asymm_sum / num;
                asymm_std = sqrt(asymm2_sum / num - asymm_avg * asymm_avg);
                //iet_bin = x_end * x_factor_2;
                //fprintf(file1_out, "%ld %.10lf %.10lf %ld\n", iet, asymm_avg, asymm_std, (long)num);
            }
            else{
                asymm_avg = 0;
            }
            iet_bin = x_end * x_factor_2;
            fprintf(file1_out, "%lf %.10lf\n", iet_bin, asymm_avg); // for RRM
            x_end *= x_factor;
            num = asymm_sum = asymm2_sum = 0;
        }
    }
    fclose(file1_out);

    free_matrix_double(asymm_Dt, 0, iet_max, 0, 2);
}

// fns for percolation
double get_memory_train(long **train_structure, long loc_start, long loc_end, long num_train){
    long train, train_before, loc, k;
    double memory, train_sum, train2_sum, train0, train1, traintrain, mean1, mean2, std1, std2;

    if(num_train <= 2) return -1;

    train_sum = train2_sum = traintrain = 0;
    loc = loc_start;
    while(loc != -1 && loc <= loc_end){
        train = train_structure[loc][0]; // train size
        train_sum += train;
        train2_sum += train * train;
        if(loc > loc_start) traintrain += train * train_before;
        train_before = train;

        loc = train_structure[loc][2];
        //if(loc == loc_end) break;
    }
    train0 = train_structure[loc_start][0];
    train1 = train_structure[loc_end][0];

    mean1 = (train_sum - train1) / (double)(num_train - 1);
    mean2 = (train_sum - train0) / (double)(num_train - 1);
    std1 = sqrt((train2_sum - train1 * train1) / (double)(num_train - 1) - mean1 * mean1);
    std2 = sqrt((train2_sum - train0 * train0) / (double)(num_train - 1) - mean2 * mean2);
    traintrain /= (double)(num_train - 1);

    memory = (traintrain - mean1 * mean2) / std1 / std2;

    return memory;
}

void measure_percol(long **train_structure, long num_event, long loc_start, long loc_end, long num_train, double *memory, long *LC, long *suscept){
    long train, train_before, loc, k, train_max;
    double train_sum, train2_sum, train0, train1, traintrain, mean1, mean2, std1, std2;

    if(num_train <= 2){
        *memory = -1;
        *LC = num_event;
        *suscept = 0;
        return;
    }

    train_sum = train2_sum = traintrain = train_max = 0;
    loc = loc_start;
    while(loc != -1 && loc <= loc_end){
        train = train_structure[loc][0]; // train size
        if(train > train_max) train_max = train;

        train_sum += train;
        train2_sum += train * train;
        if(loc > loc_start) traintrain += train * train_before;
        train_before = train;

        loc = train_structure[loc][2];
    }
    train0 = train_structure[loc_start][0];
    train1 = train_structure[loc_end][0];

    mean1 = (train_sum - train1) / (double)(num_train - 1);
    mean2 = (train_sum - train0) / (double)(num_train - 1);
    std1 = sqrt((train2_sum - train1 * train1) / (double)(num_train - 1) - mean1 * mean1);
    std2 = sqrt((train2_sum - train0 * train0) / (double)(num_train - 1) - mean2 * mean2);
    traintrain /= (double)(num_train - 1);

    *memory = (traintrain - mean1 * mean2) / std1 / std2;
    *LC = train_max;
    *suscept = train2_sum - train_max * train_max;
}

void print_train_structure(long **train_structure, long num_event){
    long i, left, right;
    for(i = 0; i < num_event; i ++){
        left = train_structure[i][1];
        right = train_structure[i][2];
        if(left != -1 || right != -1){
            printf("loc %ld (b=%ld): prev_loc %ld, next_loc %ld\n", i, train_structure[i][0], left, right);
        }
    }
}

void get_percolation_all(char *folder, char *filename, long **tree, long num_iet, double *logbin_params){ 
    long i, k, num_event, iet_buffer, iet_max, iet, bv, bw, rank_left, rank_right;
    long **train_structure, loc_start, loc_end, loc, loc1, loc2, loc_left, loc_right;
    long LC, suscept;
    double memory, **memories;
    char output1[500], output2[500], prefix[500];
    FILE *file1_out;

    strcpy(prefix, "percol_all");
    sprintf(output1, "%s%s_%s.txt", folder, prefix, filename); 
    file1_out = fopen(output1, "w");

    num_event = num_iet + 1;

    train_structure = matrix_long(0, num_event - 1, 0, 2);
    for(i = 0; i < num_event; i ++){
        train_structure[i][0] = 0; // train size
        train_structure[i][1] = -1; // for i's previous train
        train_structure[i][2] = -1; // for i's next train
    }

    iet_max = tree[0][1]; // iet_max
    //memories = matrix_double(0, iet_max, 0, 1);
    //for(i = 0; i <= iet_max; i ++) memories[i][0] = memories[i][1] = 0; 

    // percolation
    iet_buffer = iet_max;
    loc_start = loc_end = tree[0][0]; // loc of the root 
    //k = 0;
    for(i = 0; i < num_iet; i ++){ // rank or node

        loc = tree[i][0];
        iet = tree[i][1];
        rank_left = tree[i][2];
        rank_right = tree[i][3];
        bv = tree[i][4];
        bw = tree[i][5]; 

        //printf("percol %ld\n", iet);
        if(iet != iet_buffer){ // measure M_b, LC, suscept whenever iet changes
            measure_percol(train_structure, num_event, loc_start, loc_end, i + 1, &memory, &LC, &suscept);
            //printf("loc_start %ld, loc_end %ld, number of trains %ld, iet %ld\n", loc_start, loc_end, i + 1, iet);
            //memory = get_memory_train(train_structure, loc_start, loc_end, i + 1);
            //memories[k][0] = iet;
            //memories[k][1] = memory;
            //k ++;
            fprintf(file1_out, "%ld %lf %lf %lf\n", iet, memory, LC / (double)num_event, suscept / (double)num_event);
            iet_buffer = iet;
        }
        
        loc_left = train_structure[loc][1]; // i's previous train
        loc_right = train_structure[loc][2]; // i's next train

        // locations of i's left and right child nodes
        if(rank_left >= 0) loc1 = tree[rank_left][0];
        else loc1 = loc; // leaf
        if(rank_right >= 0) loc2 = tree[rank_right][0];
        else loc2 = loc + 1; // leaf
        
        // construct train "chain" in terms of the location info.:
        // loc_left (or -1) -- loc1 -- loc2 -- loc_right (or -1)
        train_structure[loc1][0] = bv;
        train_structure[loc1][2] = loc2;
        if(loc_left >= 0){
            train_structure[loc_left][2] = loc1;
            train_structure[loc1][1] = loc_left;
        }
        else loc_start = loc1; 
        
        train_structure[loc2][0] = bw;
        train_structure[loc2][1] = loc1; 
        if(loc_right >= 0){
            train_structure[loc_right][1] = loc2;
            train_structure[loc2][2] = loc_right;
        }
        else loc_end = loc2;

        //print_train_structure(train_structure, num_event);
    }

    //for(i = k - 1; i >= 0; i --){
    //    fprintf(file1_out, "%lf %lf\n", memories[i][0], memories[i][1]);
    //}
    fclose(file1_out);

    //Reverse_list_double(memories, 1, 0, k - 1);
    //get_logbin_number_direct(folder, filename, prefix, memories, k, logbin_params);

    free_matrix_long(train_structure, 0, num_event - 1, 0, 2);
    //free_matrix_double(memories, 0, iet_max, 0, 1);
}

void get_percolation(char *folder, char *filename, long **tree, long num_iet, double *logbin_params){ 
    long i, k, num_event, iet_buffer, iet_max, iet, bv, bw, rank_left, rank_right;
    long **train_structure, loc_start, loc_end, loc, loc1, loc2, loc_left, loc_right;
    long LC, suscept;
    double memory, **memories;
    char output1[500], output2[500], prefix[500];
    FILE *file1_out;

    num_event = num_iet + 1;

    train_structure = matrix_long(0, num_event - 1, 0, 2);
    for(i = 0; i < num_event; i ++){
        train_structure[i][0] = 0; // train size
        train_structure[i][1] = -1; // for i's previous train
        train_structure[i][2] = -1; // for i's next train
    }

    iet_max = tree[0][1]; // iet_max
    memories = matrix_double(0, iet_max, 0, 1);
    for(i = 0; i <= iet_max; i ++) memories[i][0] = memories[i][1] = 0; 

    // percolation
    iet_buffer = iet_max;
    loc_start = loc_end = tree[0][0]; // loc of the root 
    k = 0;
    for(i = 0; i < num_iet; i ++){ // rank or node

        loc = tree[i][0];
        iet = tree[i][1];
        rank_left = tree[i][2];
        rank_right = tree[i][3];
        bv = tree[i][4];
        bw = tree[i][5]; 

        //printf("percol %ld\n", iet);
        if(iet != iet_buffer){ // measure M_b, LC, suscept whenever iet changes
            //measure_percol(train_structure, num_event, loc_start, loc_end, i + 1, &memory, &LC, &suscept);
            //printf("loc_start %ld, loc_end %ld, number of trains %ld, iet %ld\n", loc_start, loc_end, i + 1, iet);
            memory = get_memory_train(train_structure, loc_start, loc_end, i + 1);
            memories[k][0] = iet;
            memories[k][1] = memory;
            k ++;
            //fprintf(file1_out, "%ld %lf %lf %lf\n", iet, LC / (double)num_event, suscept / (double)num_event, memory);
            iet_buffer = iet;
        }
        
        loc_left = train_structure[loc][1]; // i's previous train
        loc_right = train_structure[loc][2]; // i's next train

        // locations of i's left and right child nodes
        if(rank_left >= 0) loc1 = tree[rank_left][0];
        else loc1 = loc; // leaf
        if(rank_right >= 0) loc2 = tree[rank_right][0];
        else loc2 = loc + 1; // leaf
        
        // construct train "chain" in terms of the location info.:
        // loc_left (or -1) -- loc1 -- loc2 -- loc_right (or -1)
        train_structure[loc1][0] = bv;
        train_structure[loc1][2] = loc2;
        if(loc_left >= 0){
            train_structure[loc_left][2] = loc1;
            train_structure[loc1][1] = loc_left;
        }
        else loc_start = loc1; 
        
        train_structure[loc2][0] = bw;
        train_structure[loc2][1] = loc1; 
        if(loc_right >= 0){
            train_structure[loc_right][1] = loc2;
            train_structure[loc2][2] = loc_right;
        }
        else loc_end = loc2;

        //print_train_structure(train_structure, num_event);
    }

    strcpy(prefix, "percol");
    sprintf(output1, "%s%s_%s.txt", folder, prefix, filename); 
    file1_out = fopen(output1, "w");
    for(i = k - 1; i >= 0; i --){
        fprintf(file1_out, "%lf %lf\n", memories[i][0], memories[i][1]);
    }
    fclose(file1_out);

    Reverse_list_double(memories, 1, 0, k - 1);

    get_logbin_number_direct(folder, filename, prefix, memories, k, logbin_params);

    free_matrix_long(train_structure, 0, num_event - 1, 0, 2);
    free_matrix_double(memories, 0, iet_max, 0, 1);
}

// estimate K(b,b') for conserving the correlations between b and b'
void get_kernel2D(long **tree, long num_iet, long num_event_temp, double **kernel2D){
    long i, j, k, num_event, train_min, train, train1, train2;
    long *train_distr;
    double n_1, **QQ;

    num_event = num_iet + 1;
    train_distr = vector_long(1, num_event);

    QQ = matrix_double(1, num_event_temp, 1, num_event_temp);
    for(i = 1; i <= num_event_temp; i ++){ // i = 0 is for Q(b), otherwise kernel2D[b][b']
        for(j = 1; j <= num_event_temp; j ++){
            kernel2D[i][j] = 0;
            QQ[i][j] = 0;
        }
    }

    for(i = 1; i < num_event; i ++) train_distr[i] = 0;
    train_distr[num_event] = 1;
    train_min = num_event;

    //printf("kernel 1\n");
    for(i = 0; i < num_iet; i ++){
        train1 = tree[i][4];
        train2 = tree[i][5];
        train = train1 + train2;

        if(train1 < train_min) train_min = train1;
        if(train2 < train_min) train_min = train2;

        train_distr[train] --;
        train_distr[train1] ++;
        train_distr[train2] ++;

        n_1 = 1. / (double)(i + 1); // inverse of the number of bursts
        for(j = train_min; j <= num_event_temp; j ++){
            for(k = train_min; k <= num_event_temp; k ++){
                if(train_distr[j] && train_distr[k]){
                    QQ[j][k] += n_1 * train_distr[j] * n_1 * train_distr[k];
                    //kernel2D[0][j] += n_1 * train_distr[j];
                }
            }
        }

        if(train1 <= num_event_temp && train2 <= num_event_temp)
            kernel2D[train1][train2] ++;
    } 
    //printf("kernel 2\n");
    
    for(i = train_min; i <= num_event_temp; i ++){
        for(j = train_min; j <= num_event_temp; j ++){
            if(kernel2D[i][j] && QQ[i][j])
                kernel2D[i][j] /= QQ[i][j];
            //if(kernel2D[i][j] && kernel2D[0][i] && kernel2D[0][j])
            //   kernel2D[i][j] /= (kernel2D[0][i] * kernel2D[0][j] / num_event_temp);
        }
    }

    free_vector_long(train_distr, 1, num_event);
    free_matrix_double(QQ, 1, num_event_temp, 1, num_event_temp);
}

// for binning
long get_mid_bin(double *mid_bin, double binStart, double binSize, long num_event){
    long train2, index2;
    double x_factor, x_factor_2, y_end; 
    
    x_factor = exp(binSize);
    x_factor_2 = exp(-binSize * 0.5);

    index2 = 0; train2 = 1; // initialize
    while(index2 < binStart){ // no binning
        mid_bin[index2] = train2;
        index2 ++;
        train2 ++;
    } 
    y_end = train2 * x_factor; // bin range starting from train2
    while(train2 <= num_event){ // log binning
        if(train2 < y_end){
            train2 ++;
        }
        else{
            mid_bin[index2] = y_end * x_factor_2;
            y_end *= x_factor;
            index2 ++;
        }
    }

    return index2; // number of bins
}
        
// repeating module for a given index1 and train1
void kernel2D_yLog(double **kernel2DLog, long **nums, double **kernel2D, long num_event, double binStart, double x_factor, long index1, long train1){
    long train2, index2;
    double y_end; 
    
    // initialize
    index2 = 0; train2 = 1;

    // no binning
    while(index2 < binStart){ 
        kernel2DLog[index1][index2] += kernel2D[train1][train2];
        nums[index1][index2] ++;
        index2 ++;
        train2 ++;
    } 
    
    // bin range starting from train2
    y_end = train2 * x_factor;

    // log binning
    while(train2 <= num_event){
        if(train2 < y_end){
            kernel2DLog[index1][index2] += kernel2D[train1][train2];
            nums[index1][index2] ++;
            train2 ++;
        }
        else{
            y_end *= x_factor;
            index2 ++;
        }
    }
}
        
// get the log-binning of K_{b,b'}
void get_kernel2DLog(char *folder, char *filename, long **tree, long num_iet, double *logbin_params){
    long i, j, k, num_event, train_min, train, train1, train2, index1, index2, num_event_temp;
    long **nums, num_bin, num;
    double **kernel2D, **kernel2DLog, *mid_bin, binStart, binSize, x_factor, x_end, avg;
    char output[500], prefix[500];
    FILE *file_out;

    num_event = num_iet + 1;
    num_event_temp = 1e4;
    //num_event_temp = (long)(num_event / 100.);
    kernel2D = matrix_double(0, num_event_temp, 1, num_event_temp);

    get_kernel2D(tree, num_iet, num_event_temp, kernel2D);
    printf("kernel 2d done\n");

    binStart = logbin_params[0];
    binSize = logbin_params[1];
    x_factor = exp(binSize);

    // initialize bins' boundaries and mid values
    mid_bin = vector_double(0, num_event_temp);
    for(i = 0; i <= num_event_temp; i ++) mid_bin[i] = 0;
    num_bin = get_mid_bin(mid_bin, binStart, binSize, num_event_temp);

    // declare arrays for log binning using "num_bin"
    kernel2DLog = matrix_double(0, num_bin, 0, num_bin);
    nums = matrix_long(0, num_bin, 0, num_bin);
    for(i = 0; i <= num_bin; i ++){
        for(j = 0; j <= num_bin; j ++){
            kernel2DLog[i][j] = 0;
            nums[i][j] = 0;
        }
    }

    // initialize
    index1 = 0; train1 = 1;

    // no binning
    while(index1 < binStart){ // no binning
        kernel2D_yLog(kernel2DLog, nums, kernel2D, num_event_temp, binStart, x_factor, index1, train1);
        index1 ++;
        train1 ++;
    }
    
    // bin range starting from train2
    x_end = train1 * x_factor;

    // log binning
    while(train1 <= num_event_temp){
        if(train1 < x_end){ 
            kernel2D_yLog(kernel2DLog, nums, kernel2D, num_event_temp, binStart, x_factor, index1, train1);
            train1 ++;
        }
        else{
            x_end *= x_factor;
            index1 ++;
        }
    }
        
    // print
    strcpy(prefix, "kernel2D");
    sprintf(output, "%s%s_%sLog.txt", folder, prefix, filename);
    file_out = fopen(output, "w");

    for(i = 0; i < num_bin; i ++){
        for(j = 0; j < num_bin; j ++){
            num = nums[i][j];
            if(num){
                avg = kernel2DLog[i][j] / num;
            }
            else avg = 0;

            fprintf(file_out, "%lf %lf %.10lf\n", mid_bin[i], mid_bin[j], avg);
        }
        fprintf(file_out, "\n");
    }
    fclose(file_out);

    free_matrix_double(kernel2D, 0, num_event_temp, 1, num_event_temp);
    free_vector_double(mid_bin, 0, num_event_temp);
    free_matrix_double(kernel2DLog, 0, num_bin, 0, num_bin);
    free_matrix_long(nums, 0, num_bin, 0, num_bin);
}

// get the log-binning of K_{b,b'}
void get_kernelDiagLog(char *folder, char *filename, long **tree, long num_iet, double *logbin_params){
    long i, j, k, num_event, train_min, train, train1, train2, index1, index2, num_event_temp;
    long **nums, num_bin, num;
    double **kernel2D, **kernel2DLog, *mid_bin, binStart, binSize, x_factor, x_end, avg;
    char output[500], prefix[500];
    FILE *file_out;

    num_event = num_iet + 1;
    num_event_temp = 1e4;
    //num_event_temp = (long)(num_event / 100.);
    kernel2D = matrix_double(0, num_event_temp, 1, num_event_temp);

    get_kernel2D(tree, num_iet, num_event_temp, kernel2D);
    printf("kernel 2d done\n");

    binStart = logbin_params[0];
    binSize = logbin_params[1];
    x_factor = exp(binSize);

    // initialize bins' boundaries and mid values
    mid_bin = vector_double(0, num_event_temp);
    for(i = 0; i <= num_event_temp; i ++) mid_bin[i] = 0;
    num_bin = get_mid_bin(mid_bin, binStart, binSize, num_event_temp);

    // declare arrays for log binning using "num_bin"
    kernel2DLog = matrix_double(0, num_bin, 0, num_bin);
    nums = matrix_long(0, num_bin, 0, num_bin);
    for(i = 0; i <= num_bin; i ++){
        for(j = 0; j <= num_bin; j ++){
            kernel2DLog[i][j] = 0;
            nums[i][j] = 0;
        }
    }

    // initialize
    index1 = 0; train1 = 1;

    // no binning
    while(index1 < binStart){ // no binning
        kernel2D_yLog(kernel2DLog, nums, kernel2D, num_event_temp, binStart, x_factor, index1, train1);
        index1 ++;
        train1 ++;
    }
    
    // bin range starting from train2
    x_end = train1 * x_factor;

    // log binning
    while(train1 <= num_event_temp){
        if(train1 < x_end){ 
            kernel2D_yLog(kernel2DLog, nums, kernel2D, num_event_temp, binStart, x_factor, index1, train1);
            train1 ++;
        }
        else{
            x_end *= x_factor;
            index1 ++;
        }
    }
        
    // print
    strcpy(prefix, "kernelDiag");
    sprintf(output, "%s%s_%sLog.txt", folder, prefix, filename);
    file_out = fopen(output, "w");

    for(i = 0; i < num_bin; i ++){
        num = nums[i][i];
        if(num) avg = kernel2DLog[i][i] / num;
        else avg = 0; 
        
        fprintf(file_out, "%lf %.10lf\n", mid_bin[i], avg);
    }
    fclose(file_out);

    free_matrix_double(kernel2D, 0, num_event_temp, 1, num_event_temp);
    free_vector_double(mid_bin, 0, num_event_temp);
    free_matrix_double(kernel2DLog, 0, num_bin, 0, num_bin);
    free_matrix_long(nums, 0, num_bin, 0, num_bin);
}

/* estimate the diagonal of K(b,b') for conserving the correlations between b and b'
void get_kernelDiagLog(char *folder, char *filename, long **tree, long num_iet, double *logbin_params){
    long i, j, k, num_event, train_min, train, train1, train2;
    long *train_distr;
    double n_1, **kernelDiag, *QQ;
    char prefix[500], output[500];
    FILE *file_out;

    num_event = num_iet + 1;
    train_distr = vector_long(1, num_event);

    kernelDiag = matrix_double(0, num_event, 0, 1);
    QQ = vector_double(0, num_event);
    for(i = 0; i <= num_event; i ++) kernelDiag[i][0] = kernelDiag[i][1] = QQ[i] = 0;

    for(i = 1; i < num_event; i ++) train_distr[i] = 0;
    train_distr[num_event] = 1;
    train_min = num_event;

    //printf("kernel 1\n");
    for(i = 0; i < num_iet; i ++){
        train1 = tree[i][4];
        train2 = tree[i][5];
        train = train1 + train2;

        if(train1 < train_min) train_min = train1;
        if(train2 < train_min) train_min = train2;

        train_distr[train] --;
        train_distr[train1] ++;
        train_distr[train2] ++;

        if(train1 == train2){
            n_1 = 1. / (double)(i + 1); // inverse of the number of bursts
            for(j = train_min; j <= num_event; j ++){
                if(train_distr[j]){
                    QQ[j] += n_1 * train_distr[j] * n_1 * train_distr[j];
                }
            } 
            if(train1 <= num_event) kernelDiag[train1][1] ++;
        }
    } 
    //printf("kernel 2\n");
    
    for(i = train_min; i <= num_event; i ++){
        if(kernelDiag[i][1] && QQ[i]) kernelDiag[i][1] /= QQ[i];
    }

    // print
    strcpy(prefix, "kernelDiag");
    sprintf(output, "%s%s_%s.txt", folder, prefix, filename);
    file_out = fopen(output, "w");
    for(i = train_min; i <= num_event; i ++){
        fprintf(file_out, "%ld %lf\n", i, kernelDiag[i][1]);
        kernelDiag[i][0] = i;
    }
    fclose(file_out);

    get_logbin_number_direct(folder, filename, prefix, kernelDiag, num_event, logbin_params);

    free_vector_long(train_distr, 1, num_event);
    free_matrix_double(kernelDiag, 0, num_event, 0, 1);
    free_vector_double(QQ, 0, num_event);
}*/
