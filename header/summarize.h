// Author: Hang-Hyun Jo (h2jo23@gmail.com)
// Description: functions for summarizing analysis results

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

// summarize curves with 5%, 50%, 95% percentile
void summarize_curves(char *folder, char *filename, char *prefix, long ENS, long option){
    long i, j, k, x_min, x_max, ens;
    double x, y, z, **curves, **summary, *temp, x5, x50, x95;
    char input[500], output[500], postfix[100];
    FILE *file_in, *file_out;

    if(option == 0) sprintf(postfix, "");
    else if(option == 1) sprintf(postfix, "Log");

    x_min = 0;
    x_max = 0;
    for(ens = 0; ens < ENS; ens ++){
        sprintf(input, "%s%s_%s_ens%ld%s.txt", folder, prefix, filename, ens, postfix);
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
        sprintf(input, "%s%s_%s_ens%ld%s.txt", folder, prefix, filename, ens, postfix);
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
        find_percentile(curves, temp, i, ENS, &x5, &x50, &x95);
        fprintf(file_out, "%lf %.10lf %.10lf %.10lf\n", curves[i][0], x5, x50, x95);
    }
    fclose(file_out);

    free_matrix_double(curves, x_min, x_max, 0, ENS);
    //free_matrix_double(summary, x_min, x_max, 0, 3);
    free_vector_double(temp, 0, ENS);
}

// summarize values into p-value and distribution
void summarize_values(char *folder, char *filename, char *rrm, char *prefix, long ENS, double resolution){
    long i, j, k, *x_distr, ens, x_min_long, x_max_long;
    double x0, x, p_value, *xs, x_min, x_max;
    char input[500], output[500], filename_distr[500];
    FILE *file_in, *file_out;

    // original value from the original event sequence
    sprintf(input, "%s%s_%s.txt", folder, prefix, filename); // e.g., prefix = "memory_iet"
    file_in = fopen(input, "r");
    fscanf(file_in, "%lf", &x0);
    fclose(file_in);

    // from ensemble
    xs = vector_double(0, ENS - 1);
    for(i = 0; i < ENS; i ++) xs[i] = 0;

    sprintf(input, "%s%s_%s%s_ens.txt", folder, prefix, filename, rrm);
    file_in = fopen(input, "r");

    k = 0;
    x_min = 1e8;
    x_max = -1e8;
    for(ens = 0; ens < ENS; ens ++){
        fscanf(file_in, "%lf", &x);
        xs[ens] = x;
        if(x > x0) k ++;
        if(x < x_min) x_min = x;
        if(x > x_max) x_max = x;
    }
    fclose(file_in);

    // distribution
    x_min_long = (long)(x_min / resolution - 0.5);
    x_max_long = (long)(x_max / resolution + 0.5);
    x_distr = vector_long(x_min_long, x_max_long);
    printf("x_min_long %ld x_max_long %ld\n", x_min_long, x_max_long);

    get_distr_double(xs, ENS, resolution, x_distr, x_min_long, x_max_long);
    sprintf(filename_distr, "%s%s_distr", filename, rrm);
    print_distr(folder, filename_distr, prefix, x_distr, x_min_long, x_max_long);

    // p-value
    p_value = (double)k / (double)ENS;
    if(p_value > 0.5) p_value = 1. - p_value;

    sprintf(output, "%s%s_%s%s_pvalue.txt", folder, prefix, filename, rrm);
    file_out = fopen(output, "w");
    fprintf(file_out, "%lf\n", p_value);
    fclose(file_out);

    free_vector_double(xs, 0, ENS - 1);
    free_vector_long(x_distr, x_min_long, x_max_long);
}
