/* Log Binning functions by H.-H. Jo
*/

#include <math.h>
#define Nmix 5 // 5

// log binning without error: divide by the length of period
long LogBin_period(double **p_x,long n1,long n2,long n0,double binsize){
	long record_i=n1,read_i=n0,i,count;
	double xj,sum_y,bs_2=binsize*.5,factor=-log(1.-exp(-binsize));
	double tmp;

	if(n0<n1){
		printf("n0 should not be smaller than n1!\n");
		return -1;
	}

	for(i=n1;i<=n2;i++) if(p_x[i][0]>0) p_x[i][0]=log(p_x[i][0]);
	for(i=n1;i<n0;i++){
		tmp=p_x[i][1];
		if(tmp>0) p_x[i][1]=log(tmp),record_i++;
	}

	xj=p_x[n0][0];
	do{
		sum_y=0.;
		count=0;
		xj+=binsize;
		while(read_i<=n2 && p_x[read_i][0]<xj){ // sum_i for xj <= p_x[i][0] < xj+binsize
			sum_y+=p_x[read_i][1];
			count++;
			read_i++;
		}
		if(count==1){ // case of n=1
			i=read_i-1;
			if(p_x[i][1]){
				p_x[record_i][0]=p_x[i][0];
				p_x[record_i][1]=log(p_x[i][1]);
				//p_x[record_i][1]=log(sum_y)-xj+factor;
				record_i++;
			}
		}
		else if(count>1){ // case of n>1
			p_x[record_i][0]=xj-bs_2;
			p_x[record_i][1]=log(sum_y)-xj+factor;
			//p_x[record_i][1]=log(sum_y/(double)count);
			record_i++;
		}
	}while(read_i<=n2);
	record_i--;
	return record_i;
}

// log binning without error: divide by the length of period
long LogBin_period_all(double **p_x,long n1,long n2,long n0,double binsize){
	long record_i=n1,read_i=n0,i,count;
	double xj,sum_y,bs_2=binsize*.5,factor=-log(1.-exp(-binsize));
	double tmp;

	if(n0<n1){
		printf("n0 should not be smaller than n1!\n");
		return -1;
	}

	for(i=n1;i<=n2;i++) if(p_x[i][0]>0) p_x[i][0]=log(p_x[i][0]);
	for(i=n1;i<n0;i++){
		tmp=p_x[i][1];
		if(tmp>0) p_x[i][1]=log(tmp),record_i++;
	}

	xj=p_x[n0][0];
	do{
		sum_y=0.;
		count=0;
		xj+=binsize;
		while(read_i<=n2 && p_x[read_i][0]<xj){ // sum_i for xj <= p_x[i][0] < xj+binsize
			sum_y+=p_x[read_i][1];
			count++;
			read_i++;
		}
		if(count>=1){
			p_x[record_i][0]=xj-bs_2;
			p_x[record_i][1]=log(sum_y)-xj+factor;
			//p_x[record_i][1]=log(sum_y/(double)count);
			record_i++;
		}
	}while(read_i<=n2);
	record_i--;
	return record_i;
}

// log binning without error: divide by the length of period
long LogBin_period_count1except(double **p_x,long n1,long n2,long n0,long n3,double binsize){
	long record_i=n1,read_i=n0,i,count;
	double xj,sum_y,bs_2=binsize*.5,factor=-log(1.-exp(-binsize));
	double tmp;

	if(n0<n1){
		printf("n0 should not be smaller than n1!\n");
		return -1;
	}

	for(i=n1;i<=n2;i++) if(p_x[i][0]>0) p_x[i][0]=log(p_x[i][0]);
	//for(i=n1;i<=n2;i++) p_x[i][0]=log(p_x[i][0]);
	for(i=n1;i<n0;i++){
		tmp=p_x[i][1];
		if(tmp) p_x[i][1]=log(tmp),record_i++;
	}

	xj=p_x[n0][0];
	do{
		sum_y=0.;
		count=0;
		xj+=binsize;
		while(read_i<=n2 && p_x[read_i][0]<xj){ // sum_i for xj <= p_x[i][0] < xj+binsize
			if(p_x[read_i][1]>0){
				sum_y+=p_x[read_i][1];
				count++;
			}
			read_i++;
		}
		if(count==0) continue;
		if(count==1 && xj<log(n3)){ // case of n=1
			i=read_i-1;
			if(p_x[i][1]){
				p_x[record_i][0]=p_x[i][0];
				p_x[record_i][1]=log(p_x[i][1]);
				//p_x[record_i][1]=log(sum_y)-xj+factor;
				record_i++;
			}
		}
		else{ // case of n>1
			p_x[record_i][0]=xj-bs_2;
			p_x[record_i][1]=log(sum_y)-xj+factor;
			//p_x[record_i][1]=log(sum_y/(double)count);
			record_i++;
		}
	}while(read_i<=n2);
	record_i--;
	return record_i;
}

// NEW
// log binning without error: divide by # of points within the period
// n1 < n0 < n2 : [n1,n0] no binning, [n0,n2] log binning
long LogBin_number(double **p_x,long n1,long n2,long n0,double binsize){
	long record_i,read_i,i,j,count;
	double xj,sum_y,tmp,expBin=exp(binsize),expBin2=exp(-binsize*0.5),period0=1./(1.-exp(-binsize));

	if(n0<n1){
		printf("n0 should not be smaller than n1!\n");
		return -1;
	}

	read_i=record_i=n0;
	xj=p_x[n0][0];

	do{
		sum_y=count=0;
		xj*=expBin;
		while(read_i<=n2 && p_x[read_i][0]<xj){ // sum_i for xj <= p_x[0] < xj+binsize
			sum_y+=p_x[read_i][1];
			count++;
			read_i++;
		}
		if(count==0) continue;
		else if(count==1){ // case of n=1
			i=read_i-1;
			p_x[record_i][0]=p_x[i][0];
			p_x[record_i][1]=p_x[i][1];
		}
		else if(count>1){ // case of n>1
			p_x[record_i][0]=xj*expBin2;
			p_x[record_i][1]=sum_y/(double)count;
		}
		record_i++;
	}while(read_i<=n2);
	return record_i;
}

// NEW
// log binning without error: divide by the length of period
// n1 < n0 < n2 : [n1,n0] no binning
//				  [n0,n2] division by # of points if count<Nmix
//						  division by period length if count>=Nmix
long LogBin_period_mix(double **p_x,long n1,long n2,long n0,double binsize){
	long record_i,read_i,i,j,count;
	double xj,sum_y,tmp,expBin=exp(binsize),expBin2=exp(-binsize*0.5),period0=1./(1.-exp(-binsize));

	if(n0<n1){
		printf("n0 should not be smaller than n1!\n");
		return -1;
	}

	read_i=record_i=n0;
	xj=p_x[n0][0];

	do{
		sum_y=count=0;
		xj*=expBin;
		while(read_i<=n2 && p_x[read_i][0]<xj){ // sum_i for xj <= p_x[0] < xj+binsize
			sum_y+=p_x[read_i][1];
			count++;
			read_i++;
		}
		if(count==0) continue;
		else if(count==1){ // case of n=1
			i=read_i-1;
			p_x[record_i][0]=p_x[i][0];
			p_x[record_i][1]=p_x[i][1];
		}
		else if(count>1 && count<Nmix){ // case of n>1
			p_x[record_i][0]=xj*expBin2;
			p_x[record_i][1]=sum_y/(double)count;
		}
		else if(count>=Nmix){ // case of n>1
			p_x[record_i][0]=xj*expBin2;
			p_x[record_i][1]=sum_y/xj*period0;
		}
		record_i++;
	}while(read_i<=n2);
	return record_i;
}

// NEW
// log binning without error: divide by the length of period
// n1 < n0 < n2 : [n1,n0] no binning
//				  [n0,n2] division by # of points if count<Nmix
//						  division by period length if count>=Nmix
long LogBin_period_original(double **p_x,long n1,long n2,long n0,double binsize){
	long record_i,read_i,i,j,count;
	double xj,sum_y,tmp,expBin=exp(binsize),expBin2=exp(-binsize*0.5),period0=1./(1.-exp(-binsize));

	if(n0<n1){
		printf("n0 should not be smaller than n1!\n");
		return -1;
	}

	read_i=record_i=n0;
	xj=p_x[n0][0];

	do{
		sum_y=count=0;
		xj*=expBin;
		while(read_i<=n2 && p_x[read_i][0]<xj){ // sum_i for xj <= p_x[0] < xj+binsize
			sum_y+=p_x[read_i][1];
			count++;
			read_i++;
		}
		if(count==0) continue;
		else if(count==1){ // case of n=1
			i=read_i-1;
			p_x[record_i][0]=p_x[i][0];
			p_x[record_i][1]=p_x[i][1];
		}
		else{ // case of n>1
			p_x[record_i][0]=xj*expBin2;
			p_x[record_i][1]=sum_y/xj*period0;
		}
		record_i++;
	}while(read_i<=n2);
	return record_i;
}

// NEW
// log binning without error: divide by # of points within the period
// log binning without error: divide by the length of period
// n1 < n0 < n2 : [n1,n0] no binning
//				  [n0,n2] division by # of points if count<Nmix and x<xc
//						  division by period length if (count<Nmix and x>xc) or count>=Nmix
long LogBin_period_mix_nx(double **p_x,long n1,long n2,long n0,double xc,double binsize){
	long record_i,read_i,i,j,count;
	double xj,sum_y,tmp,expBin=exp(binsize),expBin2=exp(-binsize*0.5),period0=1./(1.-exp(-binsize));

	if(n0<n1){
		printf("n0 should not be smaller than n1!\n");
		return -1;
	}

	read_i=record_i=n0;
	xj=p_x[n0][0];

	do{
		sum_y=count=0;
		xj*=expBin;
		while(read_i<=n2 && p_x[read_i][0]<xj){ // sum_i for xj <= p_x[0] < xj+binsize
			sum_y+=p_x[read_i][1];
			count++;
			read_i++;
		}
		if(count==0) continue;
		else if(count==1){ // case of n=1
			if(xj<xc){
				i=read_i-1;
				p_x[record_i][0]=p_x[i][0];
				p_x[record_i][1]=p_x[i][1];
			}
			else{
				p_x[record_i][0]=xj*expBin2;
				p_x[record_i][1]=sum_y/xj*period0;
			}
		}
		else if(count>1 && count<Nmix){ // case of n>1
			if(xj<xc){
				p_x[record_i][0]=xj*expBin2;
				p_x[record_i][1]=sum_y/(double)count;
			}
			else{
				p_x[record_i][0]=xj*expBin2;
				p_x[record_i][1]=sum_y/xj*period0;
			}
		}
		else if(count>=Nmix){ // case of n>1
			p_x[record_i][0]=xj*expBin2;
			p_x[record_i][1]=sum_y/xj*period0;
		}
		record_i++;
	}while(read_i<=n2);
	return record_i;
}

// NEW
// linear binning without error: divide by # of points within the period
// n1 < n0 < n2 : [n1,n0] no binning, [n0,n2] log binning
long LinBin_number(double **p_x,long n1,long n2,long n0,double binsize){
	long record_i,read_i,i,j,count;
	double xj,sum_y,bs_2=binsize*.5,factor=-log(1.-exp(-binsize));
	double tmp;

	if(n0<n1){
		printf("n0 should not be smaller than n1!\n");
		return -1;
	}

	read_i=record_i=n0;
	xj=p_x[n0][0];

	do{
		sum_y=count=0;
		xj+=binsize;
		while(read_i<=n2 && p_x[read_i][0]<xj){ // sum_i for xj <= p_x[0] < xj+binsize
			sum_y+=p_x[read_i][1];
			count++;
			read_i++;
		}
		if(count==0) continue;
		else if(count==1){ // case of n=1
			i=read_i-1;
			p_x[record_i][0]=p_x[i][0];
			p_x[record_i][1]=p_x[i][1];
		}
		else if(count>1){ // case of n>1
			p_x[record_i][0]=xj-bs_2;
			p_x[record_i][1]=sum_y/(double)count;
		}
		record_i++;
	}while(read_i<=n2);
	return record_i;
}

// added on 2019.1.3.

// log binning without error: divide by # of points within the period
// n1 < n0 < n2 : [n1,n0] no binning, [n0,n2] linear binning if count<Nmix, log if count>=Nmix
long LogBin_number_mix(double **p_x,long n1,long n2,long n0,double binsize){
	long record_i,read_i,i,j,count;
	double xj,sum_y,bs_2=binsize*.5,factor=-log(1.-exp(-binsize));
	double tmp;

	if(n0<n1){
		printf("n0 should not be smaller than n1!\n");
		return -1;
	}

	i=j=n1;
	do{
		if(p_x[i][0]>0 && (tmp=p_x[i][1])>0){ // skip the negative elements
			p_x[j][0]=log(p_x[i][0]);
			if(i<n0) p_x[j][1]=log(tmp);
			else if(i==n0) p_x[j][1]=tmp,xj=p_x[j][0],read_i=record_i=j;
			else p_x[j][1]=tmp;
			j++;
		}
		if(i==n0 && (p_x[i][0]<=0 || p_x[i][1]<=0)) n0++; // when n0 is skipped
		i++;
	}while(i<=n2);
	n2=j-1;

	do{
		sum_y=0.;
		count=0;
		xj+=binsize;
		while(read_i<=n2 && p_x[read_i][0]<xj){ // sum_i for xj <= p_x[0] < xj+binsize
			sum_y+=p_x[read_i][1];
			count++;
			read_i++;
		}
		if(count==1){ // case of n=1
			i=read_i-1;
			p_x[record_i][0]=p_x[i][0];
			p_x[record_i][1]=log(p_x[i][1]);
			record_i++;
		}
		else if(count>1 && count<Nmix){ // case of n>1
			p_x[record_i][0]=xj-bs_2;
			//p_x[record_i][1]=log(sum_y)-xj+factor;
			p_x[record_i][1]=log(sum_y/(double)count);
			record_i++;
		}
		else if(count>=Nmix){ // case of n>1
			p_x[record_i][0]=xj-bs_2;
			//p_x[record_i][1]=log(sum_y)-xj+factor;
			p_x[record_i][1]=log(sum_y/(double)count);
			record_i++;
		}
	}while(read_i<=n2);
	record_i--;
	return record_i;
}
