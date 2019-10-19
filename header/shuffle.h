// Author: Hang-Hyun Jo (h2jo23@gmail.com)
// Description: functions for shuffling

/*
void generate_Poisson(long *timings, long num_event, long period){
    long i;

    timings[0] = 0;
    timings[num_event - 1] = period;
    period --;
    for(i = 1; i < num_event - 1; i ++){
        timings[i] = genrand_int32() % period + 1;
    }
    QuickSort2_ascend_long1(timings, 0, num_event - 1);
}*/

// shuffle {xs[index0], xs[index0 + 1], ..., xs[index1]}
void shuffle_long1(long *xs, long index0, long index1){
    long i, k, range, temp;

    range = index1 - index0 + 1;
    for(i = index0; i <= index1; i ++){
        k = genrand_int32() % range + index0;
        if(k == i) continue;
        temp = xs[k];
        xs[k] = xs[i];
        xs[i] = temp;
    }
}

void shuffle_timing(long *timings, long num_event, long *timings_shuffle){
    long i, period;

    period = timings[num_event - 1];
    timings_shuffle[0] = 0;
    timings_shuffle[num_event - 1] = period;
    period --;
    for(i = 1; i < num_event - 1; i ++){
        timings_shuffle[i] = genrand_int32() % period + 1;
    }
    QuickSort2_ascend_long1(timings_shuffle, 0, num_event - 1);
}

void shuffle_iet(long *iet_sequence, long num_event, long *iet_sequence_shuffle){
    long i, num_iet;

    num_iet = num_event - 1;

    for(i = 0; i < num_iet; i ++)
        iet_sequence_shuffle[i] = iet_sequence[i];

    shuffle_long1(iet_sequence_shuffle, 0, num_iet - 1);
}

void swap_branch(long **tree, long root_rank, long rank_left, long rank_right){
    long bv, bw, new_rank_left, new_rank_right; 

    new_rank_left = rank_right;
    new_rank_right = rank_left; 
    
    // swap left and right children
    tree[root_rank][2] = new_rank_left; 
    tree[root_rank][3] = new_rank_right; 

    // swap bv and bw
    bv = tree[root_rank][4]; 
    bw = tree[root_rank][5]; 
    tree[root_rank][4] = bw;
    tree[root_rank][5] = bv;
}

void shuffle_tree_leftright(long **tree, long root_rank){
    long rank_left, rank_right, temp;

    rank_left = tree[root_rank][2];
    rank_right = tree[root_rank][3];

    if(rank_left < 0 && rank_right < 0) return;

    if(genrand_real2() < 0.5) swap_branch(tree, root_rank, rank_left, rank_right);

    if(rank_left >= 0) shuffle_tree_leftright(tree, rank_left);
    if(rank_right >= 0) shuffle_tree_leftright(tree, rank_right);
}

void shuffle_tree_leftright_root(long **tree){
    long rank_left, rank_right, temp;

    rank_left = tree[0][2];
    rank_right = tree[0][3];

    if(rank_left < 0 && rank_right < 0) return;

    if(genrand_real2() < 0.5){
        tree[0][0] += (tree[0][5] - tree[0][4]); // special for root
        swap_branch(tree, 0, rank_left, rank_right);
    }

    if(rank_left >= 0) shuffle_tree_leftright(tree, rank_left);
    if(rank_right >= 0) shuffle_tree_leftright(tree, rank_right);
}

// dt is a parameter for the dt-shuffling method (shuffle nodes only with its IET <= dt)
void shuffle_tree_leftright_less(long **tree, long root_rank, long dt){
    long rank_left, rank_right, iet, temp;

    iet = tree[root_rank][1];
    rank_left = tree[root_rank][2];
    rank_right = tree[root_rank][3];

    if(rank_left < 0 && rank_right < 0) return;

    if(iet <= dt && genrand_real2() < 0.5) swap_branch(tree, root_rank, rank_left, rank_right);
    //if(genrand_real2() < 0.5) swap_branch(tree, root_rank, rank_left, rank_right);

    if(rank_left >= 0) shuffle_tree_leftright_less(tree, rank_left, dt);
    if(rank_right >= 0) shuffle_tree_leftright_less(tree, rank_right, dt);
}

void shuffle_tree_leftright_root_less(long **tree, long dt){
    long rank_left, rank_right, iet, temp;

    iet = tree[0][1];
    rank_left = tree[0][2];
    rank_right = tree[0][3];

    if(rank_left < 0 && rank_right < 0) return;

    if(iet <= dt && genrand_real2() < 0.5){
    //if(genrand_real2() < 0.5){
        tree[0][0] += (tree[0][5] - tree[0][4]); // special for root
        swap_branch(tree, 0, rank_left, rank_right);
    }

    if(rank_left >= 0) shuffle_tree_leftright_less(tree, rank_left, dt);
    if(rank_right >= 0) shuffle_tree_leftright_less(tree, rank_right, dt);
}

// dt is a parameter for the dt-shuffling method (shuffle nodes only with its IET > dt)
void shuffle_tree_leftright_more(long **tree, long root_rank, long dt){
    long rank_left, rank_right, iet, temp;

    iet = tree[root_rank][1];
    rank_left = tree[root_rank][2];
    rank_right = tree[root_rank][3];

    if(rank_left < 0 && rank_right < 0) return;

    if(iet > dt && genrand_real2() < 0.5) swap_branch(tree, root_rank, rank_left, rank_right);
    //if(genrand_real2() < 0.5) swap_branch(tree, root_rank, rank_left, rank_right);

    if(rank_left >= 0) shuffle_tree_leftright_more(tree, rank_left, dt);
    if(rank_right >= 0) shuffle_tree_leftright_more(tree, rank_right, dt);
}

void shuffle_tree_leftright_root_more(long **tree, long dt){
    long rank_left, rank_right, iet, temp;

    iet = tree[0][1];
    rank_left = tree[0][2];
    rank_right = tree[0][3];

    if(rank_left < 0 && rank_right < 0) return;

    if(iet > dt && genrand_real2() < 0.5){
    //if(genrand_real2() < 0.5){
        tree[0][0] += (tree[0][5] - tree[0][4]); // special for root
        swap_branch(tree, 0, rank_left, rank_right);
    }

    if(rank_left >= 0) shuffle_tree_leftright_more(tree, rank_left, dt);
    if(rank_right >= 0) shuffle_tree_leftright_more(tree, rank_right, dt);
}

void update_tree_loc(long **tree, long root_rank){
    long rank_left, rank_right, root_loc;

    root_loc = tree[root_rank][0];
    rank_left = tree[root_rank][2];
    rank_right = tree[root_rank][3];

    if(rank_left >= 0){
        tree[rank_left][0] = root_loc - tree[rank_left][5];
        update_tree_loc(tree, rank_left);
    }
    if(rank_right >= 0){
        tree[rank_right][0] = root_loc + tree[rank_right][4];
        update_tree_loc(tree, rank_right);
    }
}
