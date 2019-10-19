// Author: Hang-Hyun Jo (h2jo23@gmail.com)
// Description: functions used to read and print files

long count_line(char *input){
    long i;
    char c;
    FILE *file_in;

    file_in = fopen(input, "r");
    i = 0;
    c = fgetc(file_in);
    while(c != EOF){
        if(c == '\n') i ++;
        c = fgetc(file_in);
    }
    fclose(file_in);

    return i;
}

void read_long_vector(char *folder, char *filename, long *items, long num_item){
    long i;
    long item;
    char output[500];
    FILE *file_in;

    for(i = 0; i < num_item; i ++) items[i] = 0;

    sprintf(output, "%s%s.txt", folder, filename);
    file_in = fopen(output, "r");
    for(i = 0; i < num_item; i ++){
        fscanf(file_in, "%ld", &item);
        items[i] = item;
    }
    fclose(file_in);
}

void print_long_vector(char *folder, char *filename, char *prefix, long *items, long num_item){
    long i;
    char output[500];
    FILE *file_out;

    sprintf(output, "%s%s_%s.txt", folder, prefix, filename); 
    file_out = fopen(output, "w");
    for(i = 0; i < num_item; i ++){
        fprintf(file_out, "%ld\n", items[i]);
    }
    fclose(file_out);
}

/*
void print_double_vector(char *folder, char *filename, char *prefix, double *items, long num_item){
    long i;
    char output[500];
    FILE *file_out;

    sprintf(output, "%s%s_%s.txt", folder, prefix, filename); 
    file_out = fopen(output, "w");
    for(i = 0; i < num_item; i ++){
        fprintf(file_out, "%.10lf\n", items[i]);
    }
    fclose(file_out);
}*/

void print_double_matrix(char *folder, char *filename, char *prefix, double **items, long num_column, long num_row){
    long i, j;
    char output[500];
    FILE *file_out;

    sprintf(output, "%s%s_%s.txt", folder, prefix, filename); 
    file_out = fopen(output, "w");
    for(i = 0; i < num_row; i ++){
        fprintf(file_out, "%ld", (long)items[i][0]);
        for(j = 1; j < num_column; j ++){
            fprintf(file_out, " %.10lf", items[i][j]);
        }
        fprintf(file_out, "\n");
    }
    fclose(file_out);
}

void print_distr(char *folder, char *filename, char *prefix, long *distr, long x_min, long x_max){
    long i, x, total;
    char output[500];
    FILE *file_out;

    sprintf(output, "%s%s_%s.txt", folder, prefix, filename); 
    file_out = fopen(output, "w");

    total = 0; for(x = x_min; x <= x_max; x ++) total += distr[x];

    for(x = x_min; x <= x_max; x ++){
        //if(distr[x]){
            fprintf(file_out, "%ld %ld %ld\n", x, distr[x], total);
        //}
    }
    fclose(file_out);
}

void print_distr_noTotal(char *folder, char *filename, char *prefix, long *distr, long x_min, long x_max){
    long i, x, total;
    char output[500];
    FILE *file_out;

    sprintf(output, "%s%s_%s.txt", folder, prefix, filename); 
    file_out = fopen(output, "w");
    for(x = x_min; x <= x_max; x ++){
        fprintf(file_out, "%ld %ld\n", x, distr[x]);
    }
    fclose(file_out);
}

void print_distr_matrix(char *folder, char *filename, char *prefix, long **distr, long x_min, long x_max){
    long i, x, total, num;
    char output[500];
    FILE *file_out;

    sprintf(output, "%s%s_%s.txt", folder, prefix, filename); 
    file_out = fopen(output, "w");

    total = 0; for(x = x_min; x <= x_max; x ++) total += distr[x][0];

    num = 0;
    for(x = x_min; x <= x_max; x ++){
        if(distr[x][0]){
            fprintf(file_out, "%ld %ld %ld %ld %ld\n", x, distr[x][0], total, distr[x][1], distr[x][2]);
            num += distr[x][0] * x;
        }
    }
    fclose(file_out);
    //printf("num=%ld\n", num);
}
