#include <stdlib.h>
#define NR_END 1
#define FREE_ARG char*

long *vector_long(long nl,long nh){
        long *v;
        v=(long *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(long)));
        return v-nl+NR_END;
}

void free_vector_long(long *v,long nl,long nh){
       free((FREE_ARG) (v+nl-NR_END));
}

long **matrix_long(long nrl,long nrh,long ncl,long nch){
        long i,nrow=nrh-nrl+1,ncol=nch-ncl+1;
        long **m;

        m=(long **)malloc((size_t)((nrow+NR_END)*sizeof(long*)));
        m+=NR_END;
        m-=nrl;

        m[nrl]=(long *)malloc((size_t)((nrow*ncol+NR_END)*sizeof(long)));
        m[nrl]+=NR_END;
        m[nrl]-=ncl;

        for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
        return m;
}

void free_matrix_long(long **m,long nrl,long nrh,long ncl,long nch){
        free((FREE_ARG) (m[nrl]+ncl-NR_END));
        free((FREE_ARG) (m+nrl-NR_END));
}

char *vector_char(long nl,long nh){
        char *v;
        v=(char *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(char)));
        return v-nl+NR_END;
}

void free_vector_char(char *v,long nl,long nh){
       free((FREE_ARG) (v+nl-NR_END));
}

char **matrix_char(long nrl,long nrh,long ncl,long nch){
        long i,nrow=nrh-nrl+1,ncol=nch-ncl+1;
        char **m;

        m=(char **)malloc((size_t)((nrow+NR_END)*sizeof(char*)));
        m+=NR_END;
        m-=nrl;

        m[nrl]=(char *)malloc((size_t)((nrow*ncol+NR_END)*sizeof(char)));
        m[nrl]+=NR_END;
        m[nrl]-=ncl;

        for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
        return m;
}

void free_matrix_char(char **m,long nrl,long nrh,long ncl,long nch){
        free((FREE_ARG) (m[nrl]+ncl-NR_END));
        free((FREE_ARG) (m+nrl-NR_END));
}

int *vector_int(long nl,long nh){
        int *v;
        v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
        return v-nl+NR_END;
}

void free_vector_int(int *v,long nl,long nh){
       free((FREE_ARG) (v+nl-NR_END));
}

int **matrix_int(long nrl,long nrh,long ncl,long nch){
        long i,nrow=nrh-nrl+1,ncol=nch-ncl+1;
        int **m;

        m=(int **)malloc((size_t)((nrow+NR_END)*sizeof(int*)));
        m+=NR_END;
        m-=nrl;

        m[nrl]=(int *)malloc((size_t)((nrow*ncol+NR_END)*sizeof(int)));
        m[nrl]+=NR_END;
        m[nrl]-=ncl;

        for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
        return m;
}

void free_matrix_int(int **m,long nrl,long nrh,long ncl,long nch){
        free((FREE_ARG) (m[nrl]+ncl-NR_END));
        free((FREE_ARG) (m+nrl-NR_END));
}

float *vector_float(long nl,long nh){
        float *v;
        v=(float *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(float)));
        return v-nl+NR_END;
}

void free_vector_float(float *v,long nl,long nh){
       free((FREE_ARG) (v+nl-NR_END));
}

float **matrix_float(long nrl,long nrh,long ncl,long nch){
        long i,nrow=nrh-nrl+1,ncol=nch-ncl+1;
        float **m;

        m=(float **)malloc((size_t)((nrow+NR_END)*sizeof(float*)));
        m+=NR_END;
        m-=nrl;

        m[nrl]=(float *)malloc((size_t)((nrow*ncol+NR_END)*sizeof(float)));
        m[nrl]+=NR_END;
        m[nrl]-=ncl;

        for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
        return m;
}

void free_matrix_float(float **m,long nrl,long nrh,long ncl,long nch){
        free((FREE_ARG) (m[nrl]+ncl-NR_END));
        free((FREE_ARG) (m+nrl-NR_END));
}

double *vector_double(long nl,long nh){
        double *v;
        v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
        return v-nl+NR_END;
}

void free_vector_double(double *v,long nl,long nh){
       free((FREE_ARG) (v+nl-NR_END));
}

double **matrix_double(long nrl,long nrh,long ncl,long nch){
        long i,nrow=nrh-nrl+1,ncol=nch-ncl+1;
        double **m;

        m=(double **)malloc((size_t)((nrow+NR_END)*sizeof(double*)));
        m+=NR_END;
        m-=nrl;

        m[nrl]=(double *)malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
        m[nrl]+=NR_END;
        m[nrl]-=ncl;

        for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
        return m;
}

void free_matrix_double(double **m,long nrl,long nrh,long ncl,long nch){
        free((FREE_ARG) (m[nrl]+ncl-NR_END));
        free((FREE_ARG) (m+nrl-NR_END));
}

short *vector_short(long nl,long nh){
        short *v;
        v=(short *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(short)));
        return v-nl+NR_END;
}

void free_vector_short(short *v,long nl,long nh){
       free((FREE_ARG) (v+nl-NR_END));
}

short **matrix_short(long nrl,long nrh,long ncl,long nch){
        long i,nrow=nrh-nrl+1,ncol=nch-ncl+1;
        short **m;

        m=(short **)malloc((size_t)((nrow+NR_END)*sizeof(short*)));
        m+=NR_END;
        m-=nrl;

        m[nrl]=(short *)malloc((size_t)((nrow*ncol+NR_END)*sizeof(short)));
        m[nrl]+=NR_END;
        m[nrl]-=ncl;

        for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
        return m;
}

void free_matrix_short(short **m,long nrl,long nrh,long ncl,long nch){
        free((FREE_ARG) (m[nrl]+ncl-NR_END));
        free((FREE_ARG) (m+nrl-NR_END));
}

unsigned short *vector_unsigned_short(long nl,long nh){
        unsigned short *v;
        v=(unsigned short *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(short)));
        return v-nl+NR_END;
}

void free_vector_unsigned_short(unsigned short *v,long nl,long nh){
       free((FREE_ARG) (v+nl-NR_END));
}

unsigned short **matrix_unsigned_short(long nrl,long nrh,long ncl,long nch){
        long i,nrow=nrh-nrl+1,ncol=nch-ncl+1;
        unsigned short **m;

        m=(unsigned short **)malloc((size_t)((nrow+NR_END)*sizeof(short*)));
        m+=NR_END;
        m-=nrl;

        m[nrl]=(unsigned short *)malloc((size_t)((nrow*ncol+NR_END)*sizeof(short)));
        m[nrl]+=NR_END;
        m[nrl]-=ncl;

        for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
        return m;
}

void free_matrix_unsigned_short(unsigned short **m,long nrl,long nrh,long ncl,long nch){
        free((FREE_ARG) (m[nrl]+ncl-NR_END));
        free((FREE_ARG) (m+nrl-NR_END));
}

unsigned long *vector_unsigned_long(long nl,long nh){
        unsigned long *v;
        v=(unsigned long *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(long)));
        return v-nl+NR_END;
}

void free_vector_unsigned_long(unsigned long *v,long nl,long nh){
       free((FREE_ARG) (v+nl-NR_END));
}

unsigned long **matrix_unsigned_long(long nrl,long nrh,long ncl,long nch){
        long i,nrow=nrh-nrl+1,ncol=nch-ncl+1;
        unsigned long **m;

        m=(unsigned long **)malloc((size_t)((nrow+NR_END)*sizeof(long*)));
        m+=NR_END;
        m-=nrl;

        m[nrl]=(unsigned long *)malloc((size_t)((nrow*ncol+NR_END)*sizeof(long)));
        m[nrl]+=NR_END;
        m[nrl]-=ncl;

        for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
        return m;
}

void free_matrix_unsigned_long(unsigned long **m,long nrl,long nrh,long ncl,long nch){
        free((FREE_ARG) (m[nrl]+ncl-NR_END));
        free((FREE_ARG) (m+nrl-NR_END));
}

// variable size for each row
long **matrix_var_long(long nrl,long nrh,long ncl,long *nch){
//long **matrix_var_long(long nrl,long nrh,long ncl,long nch){
//	long i,nrow=nrh-nrl+1,ncol=nch-ncl+1;
	long i,nrow=nrh-nrl+1,ntotal=0;
	long **m,*nch_tmp;

	m=(long **)malloc((size_t)((nrow+NR_END)*sizeof(long*)));
	m+=NR_END;
	m-=nrl;

	nch_tmp=vector_long(nrl,nrh);
	for(i=nrl;i<=nrh;i++){
		nch_tmp[i]=nch[i]-ncl+1;
		ntotal+=nch_tmp[i];
	}
//	m[nrl]=(long *)malloc((size_t)((nrow*ncol+NR_END)*sizeof(long)));
	m[nrl]=(long *)malloc((size_t)((ntotal+NR_END)*sizeof(long)));
	m[nrl]+=NR_END;
	m[nrl]-=ncl;
	
//	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+nch_tmp[i-1];
	free_vector_long(nch_tmp,nrl,nrh);
	return m;
}

void free_matrix_var_long(long **m,long nrl,long nrh,long ncl,long *nch){
        free((FREE_ARG) (m[nrl]+ncl-NR_END));
        free((FREE_ARG) (m+nrl-NR_END));
}

// variable size for each row, int *nch instead of long *nch
long **matrix_var_long_int(long nrl,long nrh,long ncl,int *nch){
//long **matrix_var_long(long nrl,long nrh,long ncl,long nch){
//	long i,nrow=nrh-nrl+1,ncol=nch-ncl+1;
	long i,nrow=nrh-nrl+1,ntotal=0;
	long **m;
	int *nch_tmp;

	m=(long **)malloc((size_t)((nrow+NR_END)*sizeof(long*)));
	m+=NR_END;
	m-=nrl;

	nch_tmp=vector_int(nrl,nrh);
	for(i=nrl;i<=nrh;i++){
		nch_tmp[i]=nch[i]-ncl+1;
		ntotal+=nch_tmp[i];
	}
//	m[nrl]=(long *)malloc((size_t)((nrow*ncol+NR_END)*sizeof(long)));
	m[nrl]=(long *)malloc((size_t)((ntotal+NR_END)*sizeof(long)));
	m[nrl]+=NR_END;
	m[nrl]-=ncl;
	
//	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+nch_tmp[i-1];
	free_vector_int(nch_tmp,nrl,nrh);
	return m;
}

void free_matrix_var_long_int(long **m,long nrl,long nrh,long ncl,int *nch){
        free((FREE_ARG) (m[nrl]+ncl-NR_END));
        free((FREE_ARG) (m+nrl-NR_END));
}

// variable size for each row, int *nch instead of long *nch
int **matrix_var_int_int(long nrl,long nrh,long ncl,int *nch){
//long **matrix_var_long(long nrl,long nrh,long ncl,long nch){
//	long i,nrow=nrh-nrl+1,ncol=nch-ncl+1;
	long i,nrow=nrh-nrl+1,ntotal=0;
	int **m,*nch_tmp;

	m=(int **)malloc((size_t)((nrow+NR_END)*sizeof(int*)));
	m+=NR_END;
	m-=nrl;

	nch_tmp=vector_int(nrl,nrh);
	for(i=nrl;i<=nrh;i++){
		nch_tmp[i]=nch[i]-ncl+1;
		ntotal+=nch_tmp[i];
	}
//	m[nrl]=(long *)malloc((size_t)((nrow*ncol+NR_END)*sizeof(long)));
	m[nrl]=(int *)malloc((size_t)((ntotal+NR_END)*sizeof(int)));
	m[nrl]+=NR_END;
	m[nrl]-=ncl;
	
//	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+nch_tmp[i-1];
	free_vector_int(nch_tmp,nrl,nrh);
	return m;
}

void free_matrix_var_int_int(int **m,long nrl,long nrh,long ncl,int *nch){
        free((FREE_ARG) (m[nrl]+ncl-NR_END));
        free((FREE_ARG) (m+nrl-NR_END));
}
