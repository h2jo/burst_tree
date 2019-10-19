// quick sort

// double,  ascending
int Partition_ascend_double1(double a[], int low, int high); 

void QuickSort_ascend_double1(double a[], int low, int high){
	int pivot; 
	if(high>low){
		pivot=Partition_ascend_double1(a, low, high); 
		QuickSort_ascend_double1(a, low, pivot-1); 
		QuickSort_ascend_double1(a, pivot+1, high); 
	}
}

int Partition_ascend_double1(double a[], int low, int high){
	int pivot=low, left=low, right=high; 
	double pivot_item=a[low], tmp; 
	while(left<right){
		while(a[left]<=pivot_item) left++; 
		while(a[right]>pivot_item) right--; 
		if(left<right){
			tmp=a[left];  a[left]=a[right];  a[right]=tmp; 
		}
	}
	a[low]=a[right]; 
	a[right]=pivot_item; 
	return right; 
}

// 1 column,  ascending,  long
int Partition_ascend_long1(long a[], int low, int high); 

void QuickSort_ascend_long1(long a[], int low, int high){
	int pivot; 
	if(high>low){
		pivot=Partition_ascend_long1(a, low, high); 
		QuickSort_ascend_long1(a, low, pivot-1); 
		QuickSort_ascend_long1(a, pivot+1, high); 
	}
}

int Partition_ascend_long1(long a[], int low, int high){
	int pivot=low, left=low, right=high; 
	long pivot_item=a[low], tmp; 
	while(left<right){
		while(a[left]<=pivot_item) left++; 
		while(a[right]>pivot_item) right--; 
		if(left<right){
			tmp=a[left];  a[left]=a[right];  a[right]=tmp; 
		}
	}
	a[low]=a[right]; 
	a[right]=pivot_item; 
	return right; 
}

// 1 column,  descending,  long
int Partition_descend_long1(long a[], int low, int high); 

void QuickSort_descend_long1(long a[], int low, int high){
	int pivot, i; 
	if(high>low){
		pivot=Partition_descend_long1(a, low, high); 
		QuickSort_descend_long1(a, low, pivot-1); 
		QuickSort_descend_long1(a, pivot+1, high); 
	}
}

int Partition_descend_long1(long a[], int low, int high){
	int pivot=low, left=low, right=high; 
	long pivot_item=a[high], tmp; 
	while(left<right){
		while(a[left]>=pivot_item) left++; 
		while(a[right]<pivot_item) right--; 
		if(left<right){
			tmp=a[left];  a[left]=a[right];  a[right]=tmp; 
		}
	}
	a[high]=a[left]; 
	a[left]=pivot_item; 
	return right; 
}

// double,  ascending,  two columns
int Partition_ascend_double2(double a[], double b[], int low, int high); 

void QuickSort_ascend_double2(double a[], double b[], int low, int high){
	int pivot; 
	if(high>low){
		pivot=Partition_ascend_double2(a, b, low, high); 
		QuickSort_ascend_double2(a, b, low, pivot-1); 
		QuickSort_ascend_double2(a, b, pivot+1, high); 
	}
}

int Partition_ascend_double2(double a[], double b[], int low, int high){
	int pivot=low, left=low, right=high; 
	double pivot_item=a[low], tmp; 
	double b0=b[low]; 
	while(left<right){
		while(a[left]<=pivot_item) left++; 
		while(a[right]>pivot_item) right--; 
		if(left<right){
			tmp=a[left];  a[left]=a[right];  a[right]=tmp; 
			tmp=b[left];  b[left]=b[right];  b[right]=tmp; 
		}
	}
	a[low]=a[right];  b[low]=b[right]; 
	a[right]=pivot_item;  b[right]=b0; 
	return right; 
}

// long,  ascending order,  many columns,  sort by a[][sort_col]
long Partition_ascend_long(long **a, long a_column, long sort_col, long low, long high); 

void QuickSort_ascend_long(long **a, long a_column, long sort_col, long low, long high){
	long pivot; 
	if(high>low){
		pivot=Partition_ascend_long(a, a_column, sort_col, low, high); 
		QuickSort_ascend_long(a, a_column, sort_col, low, pivot-1); 
		QuickSort_ascend_long(a, a_column, sort_col, pivot+1, high); 
	}
}

long Partition_ascend_long(long **a, long a_column, long sort_col, long low, long high){
	long pivot=low, left=low, right=high; 
	long pivot_item=a[low][sort_col], tmp; 
	long *a_tmp, i; 

	a_tmp=vector_long(0, a_column); 
	for(i=0; i<=a_column; i++) a_tmp[i]=a[low][i]; 

	while(left<right){
		while(a[left][sort_col]<=pivot_item) left++; 
		while(a[right][sort_col]>pivot_item) right--; 
		if(left<right)
			for(i=0; i<=a_column; i++) tmp=a[left][i], a[left][i]=a[right][i], a[right][i]=tmp; 
	}
	for(i=0; i<=a_column; i++) a[low][i]=a[right][i], a[right][i]=a_tmp[i]; 

	free_vector_long(a_tmp, 0, a_column); 
	return right; 
}

// long,  descending order,  many columns,  sort by a[][sort_col]
long Partition_descend_long(long **a, long a_column, long sort_col, long low, long high); 

void QuickSort_descend_long(long **a, long a_column, long sort_col, long low, long high){
	long pivot; 
	if(high>low){
		pivot=Partition_descend_long(a, a_column, sort_col, low, high); 
		QuickSort_descend_long(a, a_column, sort_col, low, pivot-1); 
		QuickSort_descend_long(a, a_column, sort_col, pivot+1, high); 
	}
}

long Partition_descend_long(long **a, long a_column, long sort_col, long low, long high){
	long pivot=low, left=low, right=high; 
	long pivot_item=a[high][sort_col], tmp; 
	long *a_tmp, i; 

	a_tmp=vector_long(0, a_column); 
	for(i=0; i<=a_column; i++) a_tmp[i]=a[high][i]; 

	while(left<right){
		while(a[left][sort_col]>=pivot_item) left++; 
		while(a[right][sort_col]<pivot_item) right--; 
		if(left<right)
			for(i=0; i<=a_column; i++) tmp=a[left][i], a[left][i]=a[right][i], a[right][i]=tmp; 
	}
	for(i=0; i<=a_column; i++) a[high][i]=a[left][i], a[left][i]=a_tmp[i]; 

	free_vector_long(a_tmp, 0, a_column); 
	return right; 
}

// for sorting general array: sort a[][], b[][] using b[][column]
// ascending order
long PartitionGeneral_ascend(long **a, double **b, long a_column, long b_column, long column, long low, long high); 

void QuickSortGeneral_ascend(long **a, double **b, long a_column, long b_column, long column, long low, long high){
	long pivot; 
	if(high>low){
		pivot=PartitionGeneral_ascend(a, b, a_column, b_column, column, low, high); 
		QuickSortGeneral_ascend(a, b, a_column, b_column, column, low, pivot-1); 
		QuickSortGeneral_ascend(a, b, a_column, b_column, column, pivot+1, high); 
	}
}

long PartitionGeneral_ascend(long **a, double **b, long a_column, long b_column, long column, long low, long high){
	long i, j, k, pivot=low, left=low, right=high, tmp, *a_tmp; 
	double pivot_item=b[low][column], tmp2, *b_tmp; 

	a_tmp=vector_long(0, a_column-1); 
	b_tmp=vector_double(0, b_column-1); 
	for(i=0; i<a_column; i++) a_tmp[i]=a[low][i]; 
	for(i=0; i<b_column; i++) b_tmp[i]=b[low][i]; 
	
	while(left<right){
		while(b[left][column]<=pivot_item) left++; 
		while(b[right][column]>pivot_item) right--; 
		if(left<right){
			for(i=0; i<a_column; i++)  tmp=a[left][i], a[left][i]=a[right][i], a[right][i]=tmp; 
			for(i=0; i<b_column; i++) tmp2=b[left][i], b[left][i]=b[right][i], b[right][i]=tmp2; 
		}
	}
	for(i=0; i<a_column; i++) a[low][i]=a[right][i], a[right][i]=a_tmp[i]; 
	for(i=0; i<b_column; i++) b[low][i]=b[right][i], b[right][i]=b_tmp[i]; 
	free_vector_long(a_tmp, 0, a_column-1); 
	free_vector_double(b_tmp, 0, b_column-1); 

	return right; 
}

// for sorting general array: sort a[][], b[][] using b[][column]
// descending order
long PartitionGeneral_descend(long **a, double **b, long a_column, long b_column, long column, long low, long high); 

void QuickSortGeneral_descend(long **a, double **b, long a_column, long b_column, long column, long low, long high){
	long pivot; 
	if(high>low){
		pivot=PartitionGeneral_descend(a, b, a_column, b_column, column, low, high); 
		QuickSortGeneral_descend(a, b, a_column, b_column, column, low, pivot-1); 
		QuickSortGeneral_descend(a, b, a_column, b_column, column, pivot+1, high); 
	}
}

long PartitionGeneral_descend(long **a, double **b, long a_column, long b_column, long column, long low, long high){
	long i, j, k, pivot=low, left=low, right=high, tmp, *a_tmp; 
	double pivot_item=b[low][column], tmp2, *b_tmp; 

	a_tmp=vector_long(0, a_column-1); 
	b_tmp=vector_double(0, b_column-1); 
	for(i=0; i<a_column; i++) a_tmp[i]=a[low][i]; 
	for(i=0; i<b_column; i++) b_tmp[i]=b[low][i]; 
	
	while(left<right){
		while(b[left][column]>=pivot_item) left++; 
		while(b[right][column]<pivot_item) right--; 
		if(left<right){
			for(i=0; i<a_column; i++)  tmp=a[left][i], a[left][i]=a[right][i], a[right][i]=tmp; 
			for(i=0; i<b_column; i++) tmp2=b[left][i], b[left][i]=b[right][i], b[right][i]=tmp2; 
		}
	}
	for(i=0; i<a_column; i++) a[low][i]=a[right][i], a[right][i]=a_tmp[i]; 
	for(i=0; i<b_column; i++) b[low][i]=b[right][i], b[right][i]=b_tmp[i]; 
	free_vector_long(a_tmp, 0, a_column-1); 
	free_vector_double(b_tmp, 0, b_column-1); 

	return right; 
}

//////////////////////////////////////// Working functions!

// reverse the list
void Reverse_list(long **a, long a_column, long low, long high){
	long i, left=low, right=high, tmp; 

	while(left<right){
		for(i=0; i<=a_column; i++){
			tmp=a[left][i]; a[left][i]=a[right][i]; a[right][i]=tmp; 
        }
		left++, right--; 
	}
}

void Reverse_list_double(double **a, long a_column, long low, long high){
	long i, left=low, right=high;
    double tmp; 

	while(left<right){
		for(i=0; i<=a_column; i++){
			tmp=a[left][i]; a[left][i]=a[right][i]; a[right][i]=tmp; 
        }
		left++, right--; 
	}
}

void Reverse_long1(long *a, long low, long high){
	long left=low, right=high, tmp; 

	while(left<right){
		tmp=a[left], a[left]=a[right], a[right]=tmp; 
		left++, right--; 
	}
}

void Reverse_double1(double *a, long low, long high){
	long left=low, right=high; 
    double tmp;

	while(left<right){
		tmp=a[left], a[left]=a[right], a[right]=tmp; 
		left++, right--; 
	}
}

void QuickSort2_ascend_long(long **a, long a_column, long sort_col, long low, long high){
	long key, i, j, k, index, temp; 

	if(low<high){
		k=(long)((low+high)/2); 
		for(index=0; index<=a_column; index++)
			temp=a[k][index], a[k][index]=a[low][index], a[low][index]=temp; 
		key=a[low][sort_col]; 
		i=low+1, j=high; 
		while(i<=j){
			while(i<=high && a[i][sort_col]<=key) i++; 
			while(j>=low && a[j][sort_col]>key) j--; 
			if(i<j)
				for(index=0; index<=a_column; index++)
					temp=a[i][index], a[i][index]=a[j][index], a[j][index]=temp; 
		}
		for(index=0; index<=a_column; index++)
			temp=a[j][index], a[j][index]=a[low][index], a[low][index]=temp; 
		QuickSort2_ascend_long(a, a_column, sort_col, low, j-1); 
		QuickSort2_ascend_long(a, a_column, sort_col, j+1, high); 
	}
}

void QuickSort2_ascend_long1(long *a, long low, long high){
	long key, i, j, k, temp; 

	if(low<high){
		k=(long)((low+high)/2); 
		temp=a[k], a[k]=a[low], a[low]=temp; 
		key=a[low]; 
		i=low+1, j=high; 
		while(i<=j){
			while((i<=high) && (a[i]<=key)) i++; 
			while((j>=low) && (a[j]>key)) j--; 
			if(i<j) temp=a[i], a[i]=a[j], a[j]=temp; 
		}
		temp=a[j], a[j]=a[low], a[low]=temp; 
		QuickSort2_ascend_long1(a, low, j-1); 
		QuickSort2_ascend_long1(a, j+1, high); 
	}
}

void QuickSort2_descend_long(long **a, long a_column, long sort_col, long low, long high){
	QuickSort2_ascend_long(a, a_column, sort_col, low, high); 
	Reverse_list(a, a_column, low, high); 
}

void QuickSort2_descend_long1(long *a, long low, long high){
	QuickSort2_ascend_long1(a, low, high); 
	Reverse_long1(a, low, high); 
}

void QuickSort2_ascend_double1(double *a, long low, long high){
	long i, j, k; 
    double temp, key;

	if(low<high){
		k=(long)((low+high)/2); 
		temp=a[k], a[k]=a[low], a[low]=temp; 
		key=a[low]; 
		i=low+1, j=high; 
		while(i<=j){
			while((i<=high) && (a[i]<=key)) i++; 
			while((j>=low) && (a[j]>key)) j--; 
			if(i<j) temp=a[i], a[i]=a[j], a[j]=temp; 
		}
		temp=a[j], a[j]=a[low], a[low]=temp; 
		QuickSort2_ascend_double1(a, low, j-1); 
		QuickSort2_ascend_double1(a, j+1, high); 
	}
}


// print sorted array
void print_sorted(long **a, long a_column, long low, long high){
    long i, j;
    FILE *sort_out;

    sort_out = fopen("sort.txt", "w");

    for(i=low; i<=high; i++){
        for(j=0; j<=a_column; j++) fprintf(sort_out, "%ld ", a[i][j]);
        fprintf(sort_out, "\n");
    }

    fclose(sort_out);
}

// print sorted array
void print_sorted_double(double **a, long a_column, long low, long high){
    long i, j;
    FILE *sort_out;

    sort_out = fopen("sort.txt", "w");

    for(i=low; i<=high; i++){
        for(j=0; j<=a_column; j++) fprintf(sort_out, "%lf ", a[i][j]);
        fprintf(sort_out, "\n");
    }

    fclose(sort_out);
}

// print sorted array
void print_sorted_int(int **a, long a_column, long low, long high){
    long i, j;
    FILE *sort_out;

    sort_out = fopen("sort.txt", "w");

    for(i=low; i<=high; i++){
        for(j=0; j<=a_column; j++) fprintf(sort_out, "%d ", a[i][j]);
        fprintf(sort_out, "\n");
    }

    fclose(sort_out);
}

// added on October 16, 2015
void QuickSort3_ascend_long1(long *a, long first, long last){
	long key, i, j, k, index, temp, pivot; 

	if(first<last){
        pivot = first;
        i = first;
        j = last;

        while(i < j){
            while(a[i] <= a[pivot] && i < last) i ++;
            while(a[j] > a[pivot]) j --;
            if(i < j){
                temp = a[i]; a[i] = a[j]; a[j] = temp;
            }
        }

        temp = a[pivot]; a[pivot] = a[j]; a[j] = temp;
        QuickSort3_ascend_long1(a, first, j - 1);
        QuickSort3_ascend_long1(a, j + 1, last);
    }
}

void QuickSort3_descend_long1(long *a, long first, long last){
	QuickSort3_ascend_long1(a, first, last); 
	Reverse_long1(a, first, last); 
}

void QuickSort3_ascend_double1(double *a, long first, long last){
	long key, i, j, k, index, pivot; 
    double temp;

	if(first<last){
        pivot = first;
        i = first;
        j = last;

        while(i < j){
            while(a[i] <= a[pivot] && i < last) i ++;
            while(a[j] > a[pivot]) j --;
            if(i < j){
                temp = a[i]; a[i] = a[j]; a[j] = temp;
            }
        }

        temp = a[pivot]; a[pivot] = a[j]; a[j] = temp;
        QuickSort3_ascend_double1(a, first, j - 1);
        QuickSort3_ascend_double1(a, j + 1, last);
    }
}

void QuickSort3_ascend_double(double **a, long a_column, long sort_col, long first, long last){
	long key, i, j, k, index, pivot; 
    double temp;

	if(first<last){
        pivot = first;
        i = first;
        j = last;

        while(i < j){
            while(a[i][sort_col] <= a[pivot][sort_col] && i < last) i ++;
            while(a[j][sort_col] > a[pivot][sort_col]) j --;
            if(i < j){
                for(index = 0; index <= a_column; index ++){
                    temp = a[i][index]; a[i][index] = a[j][index]; a[j][index] = temp;
                }
            }
        } 
        
        for(index = 0; index <= a_column; index ++){
            temp = a[pivot][index]; a[pivot][index] = a[j][index]; a[j][index] = temp;
        }
        QuickSort3_ascend_double(a, a_column, sort_col, first, j - 1);
        QuickSort3_ascend_double(a, a_column, sort_col, j + 1, last);
    }
}

void QuickSort3_descend_double(double **a, long a_column, long sort_col, long first, long last){
	QuickSort3_ascend_double(a, a_column, sort_col, first, last); 
	Reverse_list_double(a, a_column, first, last); 
}

void QuickSort3_ascend_long(long **a, long a_column, long sort_col, long first, long last){
	long key, i, j, k, index, temp, pivot; 

    //printf("first=%ld last=%ld\n", first, last);

	if(first<last){
        pivot = first;
        i = first;
        j = last;

        while(i < j){
            while(a[i][sort_col] <= a[pivot][sort_col] && i < last) i ++;
            //while(a[j][sort_col] > a[pivot][sort_col] && j >= first) j --;
            while(a[j][sort_col] > a[pivot][sort_col]) j --;
            if(i < j){
                for(index = 0; index <= a_column; index ++){
                    temp = a[i][index]; a[i][index] = a[j][index]; a[j][index] = temp;
                }
            }
        } 
        
        for(index = 0; index <= a_column; index ++){
            temp = a[pivot][index]; a[pivot][index] = a[j][index]; a[j][index] = temp;
        }
        QuickSort3_ascend_long(a, a_column, sort_col, first, j - 1);
        QuickSort3_ascend_long(a, a_column, sort_col, j + 1, last);
    }
}

void QuickSort3_descend_long(long **a, long a_column, long sort_col, long first, long last){
	QuickSort3_ascend_long(a, a_column, sort_col, first, last); 
	Reverse_list(a, a_column, first, last); 
}

// when the array is declared as INT
void Reverse_list_int(int **a, long a_column, long low, long high){
	long i, left=low, right=high, tmp; 

	while(left<right){
		for(i=0; i<=a_column; i++){
			tmp=a[left][i]; a[left][i]=a[right][i]; a[right][i]=tmp; 
        }
		left++, right--; 
	}
}

void QuickSort3_ascend_int(int **a, long a_column, long sort_col, long first, long last){
	long key, i, j, k, index, temp, pivot; 

	if(first<last){
        pivot = first;
        i = first;
        j = last;
        key = a[pivot][sort_col];

        printf("first %ld value %d, last %ld value %d, key %ld\n", first, a[first][sort_col], last, a[last][sort_col], key);

        while(i < j){
            while(a[i][sort_col] <= key && i < last) i ++;
            while(a[j][sort_col] > key) j --;
            if(i < j){
                for(index = 0; index <= a_column; index ++){
                    temp = a[i][index]; a[i][index] = a[j][index]; a[j][index] = temp;
                }
            }
        } 
        
        for(index = 0; index <= a_column; index ++){
            temp = a[pivot][index]; a[pivot][index] = a[j][index]; a[j][index] = temp;
        }
        QuickSort3_ascend_int(a, a_column, sort_col, first, j - 1);
        QuickSort3_ascend_int(a, a_column, sort_col, j + 1, last);
    }
}

void QuickSort3_descend_int(int **a, long a_column, long sort_col, long first, long last){
	QuickSort3_ascend_int(a, a_column, sort_col, first, last); 
	Reverse_list_int(a, a_column, first, last); 
}

// reverse the list
void Reverse_list_long_double(long **a, double *b, long a_column, long low, long high){
	long i, left=low, right=high, tmp; 
    double tmp2;

	while(left<right){
		for(i=0; i<=a_column; i++){
			tmp=a[left][i]; a[left][i]=a[right][i]; a[right][i]=tmp; 
        }
        tmp2=b[left]; b[left]=b[right]; b[right]=tmp2; 
		left++, right--; 
	}
}

void QuickSort3_ascend_long_double(long **a, double *b, long a_column, long first, long last){
	long key, i, j, k, index, pivot, temp; 
    double temp2;

	if(first<last){
        pivot = first;
        i = first;
        j = last;

        while(i < j){
            while(b[i] <= b[pivot] && i < last) i ++;
            while(b[j] > b[pivot]) j --;
            if(i < j){
                for(index = 0; index <= a_column; index ++){
                    temp = a[i][index]; a[i][index] = a[j][index]; a[j][index] = temp;
                }
                temp2 = b[i]; b[i] = b[j]; b[j] = temp2;
            }
        } 
        
        for(index = 0; index <= a_column; index ++){
            temp = a[pivot][index]; a[pivot][index] = a[j][index]; a[j][index] = temp;
        }
        temp2 = b[pivot]; b[pivot] = b[j]; b[j] = temp2;

        QuickSort3_ascend_long_double(a, b, a_column, first, j - 1);
        QuickSort3_ascend_long_double(a, b, a_column, j + 1, last);
    }
}

void QuickSort3_descend_long_double(long **a, double *b, long a_column, long first, long last){
	QuickSort3_ascend_long_double(a, b, a_column, first, last); 
	Reverse_list_long_double(a, b, a_column, first, last); 
}

