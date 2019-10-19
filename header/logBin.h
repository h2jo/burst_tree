// Author: Hang-Hyun Jo (h2jo23@gmail.com)
// Description: functions for log-binning

//get_logbin_period_noAvg() calls
//  LogBin_period_original()
//get_logbin_period() calls
//  LogBin_period_mix_nx()
//get_logbin_number_direct() calls
//  LogBin_number()
//get_logbin_number() calls
//  LogBin_number()

void get_logbin_period_noAvg(char *folder, char *filename, char *prefix, long *input_distr, long x_max, double *logbin_params){
	long i, numBin, binStart, index;
	double **distr, y, xc, binSize, avg, num;
	char output[500];
	FILE *file_out; 
    
    binStart = logbin_params[0];
    binSize = logbin_params[1];
    xc = logbin_params[2];

	distr = matrix_double(0, x_max, 0, 1);
    for(i = 0; i <= x_max; i ++) distr[i][0] = distr[i][1] = 0;

    index = 0;
    num = 0;
    for(i = 0; i <= x_max; i ++){
        y = input_distr[i];
        //if(y > 0){
            distr[index][0] = i;
            distr[index][1] = y;
            num += y;
            index ++;
        //}
    }
	for(i = 0; i < index; i ++) distr[i][1] /= num;

	numBin = LogBin_period_original(distr, 0, index - 1, binStart, binSize);
	//numBin = LogBin_period_mix_nx(distr, 0, index - 1, binStart, xc, binSize);

	sprintf(output, "%s%s_%sLog.txt", folder, prefix, filename);
    file_out = fopen(output, "w");
    for(i = 0; i < numBin; i ++){
        fprintf(file_out, "%.15lf %.15lf\n", distr[i][0], distr[i][1]);
    }
    fclose(file_out);

	free_matrix_double(distr, 0, x_max, 0, 1);
}

void get_logbin_period(char *folder, char *filename, char *prefix, long *input_distr, long x_max, double *logbin_params){
	long i, numBin, binStart, index;
	double **distr, y, xc, binSize, avg, num;
	char output[500];
	FILE *file_out;

    binStart = logbin_params[0];
    binSize = logbin_params[1];
    xc = logbin_params[2];

	distr = matrix_double(0, x_max, 0, 1);
    for(i = 0; i <= x_max; i ++) distr[i][0] = distr[i][1] = 0;

    index = 0;
    avg = num = 0;
	for(i = 0; i <= x_max; i ++){
        y = input_distr[i];
        if(y > 0){
            distr[index][0] = i;
            distr[index][1] = y;
            num += y;
            avg += (i * y);
            index ++;
        }
    }
    avg /= num;
	for(i = 0; i < index; i ++) distr[i][1] /= num;

	numBin = LogBin_period_mix_nx(distr, 0, index - 1, binStart, xc, binSize);

	sprintf(output, "%s%s_%sLog.txt", folder, prefix, filename);
    file_out = fopen(output, "w");
    for(i = 0; i < numBin; i ++)
        fprintf(file_out, "%.15lf %.15lf %lf\n", distr[i][0], distr[i][1], avg);
    fclose(file_out);

	//numBin = LogBin_period_mix_nx(distr, 0, i-1, binStart, xc, binSize);
	//numBin = LogBin_period_mix(distr, 0, i-1, binStart, binSize);
	//numBin = LogBin_number(distr, 0, i-1, binStart, binSize);

	free_matrix_double(distr, 0, x_max, 0, 1);
}

void get_logbin_number_direct(char *folder, char *filename, char *prefix, double **input, long x_max, double *logbin_params){
	long i, j, numBin, binStart, index;
	double binSize, avg, num;
	char output[500];
	FILE *file_out;

    binStart = logbin_params[0];
    binSize = logbin_params[1];

	numBin = LogBin_number(input, 0, x_max - 1, binStart, binSize);

	sprintf(output, "%s%s_%sLog.txt", folder, prefix, filename);
    file_out = fopen(output, "w");
    for(i = 0; i < numBin; i ++)
        fprintf(file_out, "%.15lf %.15lf\n", input[i][0], input[i][1]);
    fclose(file_out);
    //printf("logbin ok %s\n", output);
}

void get_logbin_number(char *folder, char *filename, char *prefix, double **input_distr, long x_max, double *logbin_params){
	long i, j, numBin, binStart, index;
	double **distr, x, y, binSize, avg, num;
	char output[500];
	FILE *file_out;

    binStart = logbin_params[0];
    binSize = logbin_params[1];

	distr = matrix_double(0, x_max, 0, 1);
    for(i = 0; i <= x_max; i ++) for(j = 0; j <= 1; j ++) distr[i][j] = 0;

    index = 0;
	for(i = 1; i <= x_max; i ++){
        x = input_distr[i][0];
        y = input_distr[i][1];
        if(x > 0 && y > 0){
            distr[index][0] = i;
            distr[index][1] = x / y;
            index ++;
        }
    }

	numBin = LogBin_number(distr, 0, index - 1, binStart, binSize);

	sprintf(output, "%s%s_%sLog.txt", folder, prefix, filename);
    file_out = fopen(output, "w");
    for(i = 0; i < numBin; i ++)
        fprintf(file_out, "%.15lf %.15lf\n", distr[i][0], distr[i][1]);
    fclose(file_out);
    printf("logbin ok %s\n", output);

	free_matrix_double(distr, 0, x_max, 0, 1);
}
