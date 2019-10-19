void find_percentile(double **curves, double *temp, long index, long ENS, double *x5, double *x50, double *x95){
    long i, ens, index_x;

    for(ens = 0; ens < ENS; ens ++){
        temp[ens] = curves[index][ens + 1];
    }
    QuickSort2_ascend_double1(temp, 0, ENS - 1);

    index_x = (long)((double)ENS * 0.05 + 0.5);
    *x5 = temp[index_x];
    index_x = (long)((double)ENS * 0.5 + 0.5);
    *x50 = temp[index_x];
    index_x = (long)((double)ENS * 0.95 + 0.5);
    *x95 = temp[index_x];
}

void find_avg_std(double **curves, long index, long ENS, double *avg, double *std){
    long i, ens, index_x;
    double x, avg0, std0;

    avg0 = std0 = 0;
    for(ens = 0; ens < ENS; ens ++){
        x = curves[index][ens + 1];
        avg0 += x;
        std0 += x * x;
    }
    avg0 /= (double)ENS;
    std0 = sqrt(std0 / (double)ENS - avg0 * avg0);

    *avg = avg0;
    *std = std0;
}

// summarize curves with avg and std
void summarize_curves(char *folder, char *filename, char *prefix, long ENS){
    long i, j, k, x_min, x_max, ens;
    double x, y, z, **curves, **summary, *temp, x5, x50, x95, avg, std;
    char input[500], output[500]; 
    FILE *file_in, *file_out;

    x_min = 0;
    x_max = 0;
    for(ens = 0; ens < ENS; ens ++){
        sprintf(input, "%s%s_%s_ens%ldLog.txt", folder, prefix, filename, ens);
        k = count_line(input);
        if(k > x_max) x_max = k;
    }
    printf("max line count=%ld\n", k);

    curves = matrix_double(x_min, x_max, 0, ENS);
    //summary = matrix_double(x_min, x_max, 0, 3);
    temp = vector_double(0, ENS);

    for(i = x_min; i <= x_max; i ++){
        for(j = 0; j <= ENS; j ++) curves[i][j] = 0;
        //for(j = 0; j <= 3; j ++) summary[i][j] = 0;
    }

    for(ens = 0; ens < ENS; ens ++){
        sprintf(input, "%s%s_%s_ens%ldLog.txt", folder, prefix, filename, ens);
        file_in = fopen(input, "r");
        k = 0;
        while(fscanf(file_in, "%lf %lf", &x, &y) && !feof(file_in)){
            curves[k][0] = x;
            curves[k][ens + 1] = y;
            k ++;
        }
        fclose(file_in);
    }

    sprintf(output, "%s%s_%s_summary.txt", folder, prefix, filename);
    file_out = fopen(output, "w");
    for(i = 0; i < k; i ++){
        find_avg_std(curves, i, ENS, &avg, &std);
        fprintf(file_out, "%lf %.10lf %.10lf\n", curves[i][0], avg, std);

        //find_percentile(curves, temp, i, ENS, &x5, &x50, &x95);
        //fprintf(file_out, "%lf %.10lf %.10lf %.10lf\n", curves[i][0], x5, x50, x95);
    }
    fclose(file_out);

    free_matrix_double(curves, x_min, x_max, 0, ENS);
    //free_matrix_double(summary, x_min, x_max, 0, 3);
    free_vector_double(temp, 0, ENS);
}

// summarize distributions
void summarize_distributions(char *folder, char *filename, char *prefix, long ENS){
    long i, j, k, x_min, x_max, ens, x, y;
    long *distr;
    double logbin_params[3];
    char input[500], output[500]; 
    FILE *file_in, *file_out;

    x_min = 0;
    x_max = 0;
    for(ens = 0; ens < ENS; ens ++){
        sprintf(input, "%s%s_%s_ens%ld.txt", folder, prefix, filename, ens);
        k = count_line(input);
        if(k > x_max) x_max = k;
    }
    printf("max line count=%ld\n", k);

    distr = vector_long(x_min, x_max);

    for(i = x_min; i <= x_max; i ++){
        distr[i] = 0;
    }

    for(ens = 0; ens < ENS; ens ++){
        sprintf(input, "%s%s_%s_ens%ld.txt", folder, prefix, filename, ens);
        file_in = fopen(input, "r");
        while(fscanf(file_in, "%ld %ld", &x, &y) && !feof(file_in)){
            distr[x] += y;
        }
        fclose(file_in);
    }

    sprintf(output, "%s%s_%s_summary.txt", folder, prefix, filename);
    file_out = fopen(output, "w");
    for(i = x_min; i <= x_max; i ++){
        y = distr[i];
        if(y) fprintf(file_out, "%ld %ld\n", i, y);
    }
    fclose(file_out);

    // print logbinned distr 
    logbin_params[0] = 1; // binStart
    logbin_params[1] = 0.2; // binSize
    logbin_params[2] = 10; // xc = crossover of x (?)
    get_logbin_period(folder, filename, prefix, distr, x_max, logbin_params);

    free_vector_long(distr, x_min, x_max);
}

