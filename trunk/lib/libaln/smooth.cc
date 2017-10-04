/******************************************************************************************
    Copyright (C) 1997-2014 Andrew F. Neuwald, Cold Spring Harbor Laboratory
    and the University of Maryland School of Medicine.

    Permission is hereby granted, free of charge, to any person obtaining a copy of 
    this software and associated documentation files (the "Software"), to deal in the 
    Software without restriction, including without limitation the rights to use, copy, 
    modify, merge, publish, distribute, sublicense, and/or sell copies of the Software,
    and to permit persons to whom the Software is furnished to do so, subject to the 
    following conditions:

    The above copyright notice and this permission notice shall be included in all 
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
    INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR 
    PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE 
    LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT 
    OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR 
    OTHER DEALINGS IN THE SOFTWARE.

    For further information contact:
         Andrew F. Neuwald
         Institute for Genome Sciences and
         Department of Biochemistry & Molecular Biology
         University of Maryland School of Medicine
         801 West Baltimore St.
         BioPark II, Room 617
         Baltimore, MD 21201
         Tel: 410-706-6724; Fax: 410-706-1482; E-mail: aneuwald@som.umaryland.edu
 ******************************************************************************************/

/* spouge.c - john spouge statistics */
#include "smooth.h"

double  *SmoothData(Int4 nl, Int4 nr, Int4 m, Int4 maxLength, Int4 *Data)
// smooth data using Savitzky-Golay smoothing Filter...
/*********************************************************************
  WARNING: ans must have dimensions [1...2*n] and n MUST be an integer
  power of two!
/*********************************************************************/
{
        double  *data,*ans,*respns,total,v,sum,*dist;
        UInt4   np,n;
        Int4    t,g,num,s;

        n = maxLength +1; np = 1;
        /*** set n = to the smallest integer power of two >= maxLength +1. ***/
        do { np *= 2; } while(np < n); n = np;
        np = nl+nr+1;
        NEW(data,n+3,double); NEW(ans,2*n+3,double); NEW(respns,n+3,double);
        NEW(dist,n+3,double);
        data++;
        for(total=0.0,g=0; g<= maxLength; g++){
          data[g] = (double) Data[g];
          total += data[g];
        }
        // fprintf(stderr,"np =%d; maxLength = %d\n",np,n);
        data--;
        savgol(respns, np, nl, nr, 0, m);
        convlv(data, n, respns, np, 1, ans);
        ans++;
        for(num=g=0; g<= maxLength; g++){ dist[g] = ans[g]; } ans--;
        for(g=0; g<= maxLength; g++){ 
		fprintf(stderr,"dist[%d] = %g\n",g,dist[g]); 
	} free(data); free(ans); free(respns);
        return dist;
}

double	*ExpectedGaps(int length_of_seq, int number_of_blocks)
/****************************************
         (Len - k + N - 1)!
  d(k)=  ------------------
         (Len-k)! (N - 1)!
 ****************************************/
{
        Int4	k,m, leftover_length;
        double	*gap_distr,max,sum;

        Int4 tot=number_of_blocks;
        leftover_length=length_of_seq-tot;
        NEW(gap_distr,length_of_seq+2,double);
        for(max = 0.0,k=0;k<=leftover_length;k++){
                gap_distr[k]=lgamma(leftover_length+number_of_blocks-k)-
                        lgamma(leftover_length-k+1)-lgamma(number_of_blocks);
                max = MAXIMUM(double,max,gap_distr[k]);
        }
        for(sum=0.0,k=0;k<=leftover_length;k++){
                gap_distr[k] = exp(gap_distr[k] - max);
// fprintf(stderr,"gap_distr[%d] = %g\n",k,gap_distr[k]);
		sum += gap_distr[k];
        }
        for(k=0;k<=leftover_length;k++){ gap_distr[k] = gap_distr[k]/sum; }
        return gap_distr;
}

double	*SmoothGapsX(Int4 nl, Int4 nr, Int4 m, Int4 maxLength, Int4 *Gaps, 
	double rpts_per_res)
// SmoothGapsX(12,12,4,maxLength,Gaps,rpts_per_res);
/*********************************************************************
  WARNING: ans must have dimensions [1...2*n] and n MUST be an integer
  power of two!
/*********************************************************************/
{
	double	*data,*ans,*respns,total,v,factor=5.0;
	Int4    adj,mode,s;
	double	sum,ave,*dist;
	// double	pseudo=0.2;
	double	pseudo=0.0;
	UInt4	np,n;
	Int4	cutoff,t,g,num,maxgap;
	BooLean	flag;

	n = maxLength +1; np = 1;
	/*** set n = to the smallest integer power of two >= maxLength +1. ***/
	do { np *= 2; } while(np < n); n = np;
	np = nl+nr+1;
	NEW(data,n+3,double); NEW(ans,2*n+3,double); NEW(respns,n+3,double);
	NEW(dist,n+3,double);
	// fprintf(stderr,"\n******************** input gaps ********************\n");

	data++; 
#if 1	// when reach sparce data then use expected gaps...
	Int4 g_cut,Sum=0,w,SumTail=0;
	Int4 window_size=5,*window;
	NEW(window,window_size+3,Int4);
	for(g_cut=maxLength+1, g=0; g<= maxLength; g++){ 
	   if(g > g_cut){
		SumTail+=Gaps[g];	// get counts in tail...
	   } else {
		w = g % window_size; // get current position
		Sum -= window[w]; window[w] = Gaps[g]; Sum += window[w];
		if(g > window_size && Sum < window_size){
			g_cut = g;
		}
	   }
	} free(window);
	double	*exp_dist=ExpectedGaps(maxLength,(Int4)(maxLength*rpts_per_res+1.0));
#endif
	for(total=0.0,g=0; g<= maxLength; g++){
	  if(g <= g_cut){ data[g] = (double) Gaps[g]; }
	  else {	// use expected gaps...
	    data[g] = (double)SumTail*exp_dist[g];
	  }
	  total += data[g];
	} free(exp_dist);
// fprintf(stderr,"g_cut = %d; SumTail = %d; rpts_per_res = %g\n",g_cut,SumTail,rpts_per_res);
	// fprintf(stderr,"np =%d; maxLength = %d\n",np,n);
	data--;
	// savgol(respns, np, nl, nr, 0, 2);
	savgol(respns, np, nl, nr, 0, m);
	convlv(data, n, respns, np, 1, ans);
	ans++;
#if 0
/****** NEW: redo with more stringent settings ****/
	ans--;
	savgol(respns, np, nl, nr, 0, m);
	convlv(data, n, respns, np, 1, ans);
	ans++;
/****** NEW: redo with more stringent settings ****/
#endif
	double min=DBL_MAX;
	for(num=g=0; g<= maxLength; g++){
		dist[g] = (double) ans[g];
		if(dist[g] > 0.0 && min > dist[g]) min=dist[g];
	} ans--;
	for(sum=0.0,g=0; g<= maxLength; g++){
		if(dist[g] <= 0.0) dist[g]=min;
		sum+=dist[g];
	}
	for(g=0; g<= maxLength; g++){ dist[g] = dist[g]/sum; }
	// for(g=0; g<= maxLength; g++){ fprintf(stderr,"dist[%d] = %g\n",g,dist[g]); }
	free(data); free(ans); free(respns); 
	return dist;
}

/*********** WARNING: copyrighted NumRecipes routines *****************/
#include <math.h>
#define NRANSI
// #include "nrutil.h"
#define TINY 1.0e-20;

void ludcmp(double **a, int n, int *indx, double *d)
{
	int i,imax,j,k;
	double big,dum,sum,temp;
	double *vv;

	vv=vector(1,n);
	*d=1.0;
	for (i=1;i<=n;i++) {
		big=0.0;
		for (j=1;j<=n;j++)
			if ((temp=fabs(a[i][j])) > big) big=temp;
		if (big == 0.0) nrerror("Singular matrix in routine ludcmp");
		vv[i]=1.0/big;
	}
	for (j=1;j<=n;j++) {
		for (i=1;i<j;i++) {
			sum=a[i][j];
			for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
		}
		big=0.0;
		for (i=j;i<=n;i++) {
			sum=a[i][j];
			for (k=1;k<j;k++)
				sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
			if ( (dum=vv[i]*fabs(sum)) >= big) {
				big=dum;
				imax=i;
			}
		}
		if (j != imax) {
			for (k=1;k<=n;k++) {
				dum=a[imax][k];
				a[imax][k]=a[j][k];
				a[j][k]=dum;
			}
			*d = -(*d);
			vv[imax]=vv[j];
		}
		indx[j]=imax;
		if (a[j][j] == 0.0) a[j][j]=TINY;
		if (j != n) {
			dum=1.0/(a[j][j]);
			for (i=j+1;i<=n;i++) a[i][j] *= dum;
		}
	}
	free_vector(vv,1,n);
}
#undef TINY
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software =$j7`5|,$. */
void lubksb(double **a, int n, int *indx, double b[])
{
	int i,ii=0,ip,j;
	double sum;

	for (i=1;i<=n;i++) {
		ip=indx[i];
		sum=b[ip];
		b[ip]=b[i];
		if (ii)
			for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
		else if (sum) ii=i;
		b[i]=sum;
	}
	for (i=n;i>=1;i--) {
		sum=b[i];
		for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
		b[i]=sum/a[i][i];
	}
}
/* (C) Copr. 1986-92 Numerical Recipes Software =$j7`5|,$. */
#include <math.h>
#define NRANSI
#include "nrutil.h"

void savgol(double c[], int np, int nl, int nr, int ld, int m)
{
	void lubksb(double **a, int n, int *indx, double b[]);
	void ludcmp(double **a, int n, int *indx, double *d);
	int imj,ipj,j,k,kk,mm,*indx;
	double d,fac,sum,**a,*b;

	if (np < nl+nr+1 || nl < 0 || nr < 0 || ld > m || nl+nr < m)
	nrerror("bad args in savgol");
	indx=ivector(1,m+1);
	a=matrix(1,m+1,1,m+1);
	b=vector(1,m+1);
	for (ipj=0;ipj<=(m << 1);ipj++) {
		sum=(ipj ? 0.0 : 1.0);
		for (k=1;k<=nr;k++) sum += pow((double)k,(double)ipj);
		for (k=1;k<=nl;k++) sum += pow((double)-k,(double)ipj);
		mm=FMIN(ipj,2*m-ipj);
		for (imj = -mm;imj<=mm;imj+=2) a[1+(ipj+imj)/2][1+(ipj-imj)/2]=sum;
	}
	ludcmp(a,m+1,indx,&d);
	for (j=1;j<=m+1;j++) b[j]=0.0;
	b[ld+1]=1.0;
	lubksb(a,m+1,indx,b);
	for (kk=1;kk<=np;kk++) c[kk]=0.0;
	for (k = -nl;k<=nr;k++) {
		sum=b[1];
		fac=1.0;
		for (mm=1;mm<=m;mm++) sum += b[mm+1]*(fac *= k);
		kk=((np-k) % np)+1;
		c[kk]=sum;
	}
	free_vector(b,1,m+1);
	free_matrix(a,1,m+1,1,m+1);
	free_ivector(indx,1,m+1);
}
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software =$j7`5|,$. */
#include <math.h>
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

void four1(double data[], UInt4 nn, int isign)
{
	UInt4 n,mmax,m,j,istep,i;
	double wtemp,wr,wpr,wpi,wi,theta;
	double tempr,tempi;

	n=nn << 1;
	j=1;
	for (i=1;i<n;i+=2) {
		if (j > i) {
			SWAP(data[j],data[i]);
			SWAP(data[j+1],data[i+1]);
		}
		m=n >> 1;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	}
	mmax=2;
	while (n > mmax) {
		istep=mmax << 1;
		theta=isign*(6.28318530717959/mmax);
		wtemp=sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi=sin(theta);
		wr=1.0;
		wi=0.0;
		for (m=1;m<mmax;m+=2) {
			for (i=m;i<=n;i+=istep) {
				j=i+mmax;
				tempr=wr*data[j]-wi*data[j+1];
				tempi=wr*data[j+1]+wi*data[j];
				data[j]=data[i]-tempr;
				data[j+1]=data[i+1]-tempi;
				data[i] += tempr;
				data[i+1] += tempi;
			}
			wr=(wtemp=wr)*wpr-wi*wpi+wr;
			wi=wi*wpr+wtemp*wpi+wi;
		}
		mmax=istep;
	}
}
#undef SWAP
/* (C) Copr. 1986-92 Numerical Recipes Software =$j7`5|,$. */
#include <math.h>

void realft(double data[], UInt4 n, int isign)
{
	void four1(double data[], UInt4 nn, int isign);
	UInt4 i,i1,i2,i3,i4,np3;
	double c1=0.5,c2,h1r,h1i,h2r,h2i;
	double wr,wi,wpr,wpi,wtemp,theta;

	theta=3.141592653589793/(double) (n>>1);
	if (isign == 1) {
		c2 = -0.5;
		four1(data,n>>1,1);
	} else {
		c2=0.5;
		theta = -theta;
	}
	wtemp=sin(0.5*theta);
	wpr = -2.0*wtemp*wtemp;
	wpi=sin(theta);
	wr=1.0+wpr;
	wi=wpi;
	np3=n+3;
	for (i=2;i<=(n>>2);i++) {
		i4=1+(i3=np3-(i2=1+(i1=i+i-1)));
		h1r=c1*(data[i1]+data[i3]);
		h1i=c1*(data[i2]-data[i4]);
		h2r = -c2*(data[i2]+data[i4]);
		h2i=c2*(data[i1]-data[i3]);
		data[i1]=h1r+wr*h2r-wi*h2i;
		data[i2]=h1i+wr*h2i+wi*h2r;
		data[i3]=h1r-wr*h2r+wi*h2i;
		data[i4] = -h1i+wr*h2i+wi*h2r;
		wr=(wtemp=wr)*wpr-wi*wpi+wr;
		wi=wi*wpr+wtemp*wpi+wi;
	}
	if (isign == 1) {
		data[1] = (h1r=data[1])+data[2];
		data[2] = h1r-data[2];
	} else {
		data[1]=c1*((h1r=data[1])+data[2]);
		data[2]=c1*(h1r-data[2]);
		four1(data,n>>1,-1);
	}
}
/* (C) Copr. 1986-92 Numerical Recipes Software =$j7`5|,$. */
void twofft(double data1[], double data2[], double fft1[], double fft2[],
	UInt4 n)
{
	void four1(double data[], UInt4 nn, int isign);
	UInt4 nn3,nn2,jj,j;
	double rep,rem,aip,aim;

	nn3=1+(nn2=2+n+n);
	for (j=1,jj=2;j<=n;j++,jj+=2) {
		fft1[jj-1]=data1[j];
		fft1[jj]=data2[j];
	}
	four1(fft1,n,1);
	fft2[1]=fft1[2];
	fft1[2]=fft2[2]=0.0;
	for (j=3;j<=n+1;j+=2) {
		rep=0.5*(fft1[j]+fft1[nn2-j]);
		rem=0.5*(fft1[j]-fft1[nn2-j]);
		aip=0.5*(fft1[j+1]+fft1[nn3-j]);
		aim=0.5*(fft1[j+1]-fft1[nn3-j]);
		fft1[j]=rep;
		fft1[j+1]=aim;
		fft1[nn2-j]=rep;
		fft1[nn3-j] = -aim;
		fft2[j]=aip;
		fft2[j+1] = -rem;
		fft2[nn2-j]=aip;
		fft2[nn3-j]=rem;
	}
}
/* (C) Copr. 1986-92 Numerical Recipes Software =$j7`5|,$. */
#define NRANSI
#include "nrutil.h"

void convlv(double data[], UInt4 n, double respns[], UInt4 m,
	int isign, double ans[])
{
	void realft(double data[], UInt4 n, int isign);
	void twofft(double data1[], double data2[], double fft1[], double fft2[],
		UInt4 n);
	UInt4 i,no2;
	double dum,mag2,*fft;

	fft=vector(1,n<<1);
	for (i=1;i<=(m-1)/2;i++)
		respns[n+1-i]=respns[m+1-i];
	for (i=(m+3)/2;i<=n-(m-1)/2;i++)
		respns[i]=0.0;
	twofft(data,respns,fft,ans,n);
	no2=n>>1;
	for (i=2;i<=n+2;i+=2) {
		if (isign == 1) {
			ans[i-1]=(fft[i-1]*(dum=ans[i-1])-fft[i]*ans[i])/no2;
			ans[i]=(fft[i]*dum+fft[i-1]*ans[i])/no2;
		} else if (isign == -1) {
			if ((mag2=SQR(ans[i-1])+SQR(ans[i])) == 0.0)
				nrerror("Deconvolving at response zero in convlv");
			ans[i-1]=(fft[i-1]*(dum=ans[i-1])+fft[i]*ans[i])/mag2/no2;
			ans[i]=(fft[i]*dum-fft[i-1]*ans[i])/mag2/no2;
		} else nrerror("No meaning for isign in convlv");
	}
	ans[2]=ans[n+1];
	realft(ans,n,-1);
	free_vector(fft,1,n<<1);
}
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software =$j7`5|,$. */
/*********** WARNING: copyrighted NumRecipes routine *****************/


