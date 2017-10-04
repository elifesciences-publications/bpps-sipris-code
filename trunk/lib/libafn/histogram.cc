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

#include "histogram.h"

h_type	Histogram(const char *id,Int4 start,Int4 end,double inc)
/* bin0[0] = underflow; bin0[1->nbins]= values; bin0[nbins+1]=overflow */
{
	h_type	H;

	NEW(H,1,histogram_type);
	NEW(H->id,(strlen(id)+1),char); strcpy(H->id,id);
	if(end <= start || inc <= 0.0){
		fprintf(stderr,"Histogram input error 1 ('%s')\n",id); 
		fprintf(stderr,"start=%d; end=%d; inc=%.3f\n",start,end,inc); 
		exit(1); 
	}
	H->min = (double) start;
	H->inc = inc;
	H->nbins = (Int4) ceil((double)(end - start)/inc);
	if(H->nbins < 1) {
		fprintf(stderr,"start=%d; end = %d; inc = %g; nbins = %d\n",start,end,inc,H->nbins);
		fprintf(stderr,"Histogram input error 2 ('%s')\n",id); exit(1); 
	}
	H->max = (double) start + (inc * (double)(H->nbins));
	H->maxval = - DBL_MAX; H->minval = DBL_MAX;
	H->n = 0; H->total = 0.0; H->total_sq = 0.0;
	NEW(H->bin0,H->nbins+3,Int4);
	H->bin=(H->bin0+1);
	return(H);
}

void	NilHist(h_type H){ free(H->id); free(H->bin0); free(H); }

Int4	RtnHistArray(double **xval, double **yval, Int4 *Total, h_type H)
/************************************************************************
 return number of bins and set xval and yval arrays 
 ************************************************************************/
{
	Int4	start,total,sum,i,j,n,m, max = 0;
	double	val,*Xval,*Yval;
	Int4	NumBins= H->nbins;
	
	NEW(Xval, NumBins+5,double);
	NEW(Yval, NumBins+5,double);

	total = 0;
	Xval[0] = H->min;
	Yval[0] = H->bin0[0];
	for(i=1;i<=H->nbins+1;i++) {
		total += H->bin0[i];
		Xval[i] = Xval[i-1] + H->inc;
		Yval[i]=H->bin0[i];
	}
	*Total=total;
	*xval=Xval;
	*yval=Yval;
	return NumBins;
}

void	PutHistArray(FILE *fptr,h_type H)
/************************************************************************
	note: Var(X) = E(X*X) - E(X)*E(X).
 ************************************************************************/
{
	Int4	start,total,sum,i,j,n,m, max = 0;
	double	val;
	
	total = 0;
	for(i=0;i<=H->nbins+1;i++) {
		total += H->bin0[i];
		max = H->bin0[i] > max ? H->bin0[i] : max;
	}
	if(total == 0) { fprintf(fptr,"empty.\n"); return; }

	//************* print out X axis... *******************
	sum = 0;
	if((n=H->bin0[0]) > 0){
		sum = n;
		fprintf(fptr,"<%.2f ",H->min);
		start = 1;
		val=H->min;
	} else {
	   for(start=0,val=H->min; H->bin0[start]==0 ;start++){
	      val+=H->inc;
	      if(start >= H->nbins) {
        	fprintf(fptr,"Histogram: variance too small\n"); return; 
	      }
	   }
	   val-=H->inc;
	}
	for(i=start; i<=H->nbins; val+= H->inc, sum+=H->bin0[i],i++) {
	   n = H->bin0[i];
	   if(sum==total) { fprintf(fptr,"%.2f ",val); break; }
	   else fprintf(fptr,"%.2f ",val);
	}
	if((n=H->bin0[H->nbins+1]) > 0){
		fprintf(fptr,">=%.2f ",H->max); 
	}
	fprintf(fptr,"\n"); 

	//************* print out Y axis... ***************
	sum = 0;
	if((n=H->bin0[0]) > 0){
		sum = n;
		fprintf(fptr,"%ld ",n);
		start = 1;
	} else {
	   for(start=0; H->bin0[start]==0 ;start++){
	      if(start >= H->nbins) {
        	fprintf(fptr,"Histogram: variance too small\n"); return; 
	      }
	   } 
	}
	for(i=start; i<=H->nbins; sum+=H->bin0[i],i++) {
	   n = H->bin0[i];
	   if(sum==total) { fprintf(fptr,"%ld ",n); break; }
	   else fprintf(fptr,"%ld ",n);
	}
	if((n=H->bin0[H->nbins+1]) > 0){ fprintf(fptr,"%ld ",n); }
	fprintf(fptr,"\n"); 
}


void    PutHistX(FILE *fp, char *str, double inc, Int4 maxLength, Int4 *dist)
// dist[g] = observed.
{
        h_type  H;

        H=Histogram(str,0,maxLength,inc);
        for(Int4 g=0; g<= maxLength; g++){
            if(dist[g] >= 0) IncdMHist(g, dist[g], H);
        } PutHist(fp,60,H); NilHist(H);
}

#if 0
void    PutHistX(FILE *fp, char *str, double inc, Int4 maxLength, double *prob)
// NOTE: prob[g] sum to 1.0.
{
        h_type  H;

        H=Histogram(str,0,maxLength,inc);
        for(Int4 g=0; g<= maxLength; g++){
                if(prob[g] >= 0.0) IncdMHist(g, 1000.0*prob[g], H);
        } PutHist(fp,60,H); NilHist(H);
}
#endif

Int4	TotalDataHist(h_type H)
{ Int4	i,total=0; for(i=0;i<=H->nbins+1;i++) { total += H->bin0[i]; } return total; }

void	PutHist(FILE *fptr,Int4 line_leng,h_type H)
/************************************************************************
	note: Var(X) = E(X*X) - E(X)*E(X).
 ************************************************************************/
{
	Int4	start,total,sum,i,j,n,m, max = 0,mult_factor;
	double	val;
	
	total = 0;
	for(i=0;i<=H->nbins+1;i++) {
		total += H->bin0[i];
		max = H->bin0[i] > max ? H->bin0[i] : max;
	}
	if(max <= line_leng) mult_factor = 1;
	else mult_factor = (max/line_leng) + 1;
	fprintf(fptr,"Distribution of %s:\n",H->id);
	if(total == 0) { fprintf(fptr,"empty.\n"); return; }

	if(mult_factor==1) fprintf(fptr, " '=' is %d count.\n", mult_factor);
	else fprintf(fptr, "  ( '=' is %d counts. )\n", mult_factor);

	sum = 0;
	if((n=H->bin0[0]) > 0){
		sum = n;
		fprintf(fptr,"\n<%7.2f : %-8d |",H->min,n);
		m = n/mult_factor; for(j=0;j<m;j++) fprintf(fptr,"=");
		start = 1;
		val=H->min;
	} else {
	   for(start=0,val=H->min; H->bin0[start]==0 ;start++){
	      val+=H->inc;
	      if(start >= H->nbins) {
        	fprintf(fptr,"Histogram: variance too small\n"); return; 
	      }
	   }
	   val-=H->inc;
	}

	for(i=start; i<=H->nbins; val+= H->inc, sum+=H->bin0[i],i++) {
	   n = H->bin0[i];
	   if(sum==total) { 
		fprintf(fptr,"\n%8.2f : %-8d |",val,n); 
		m = n/mult_factor; for(j=0;j<m;j++) fprintf(fptr,"=");
		break; 
	   } else fprintf(fptr,"\n%8.2f : %-8d |",val,n);
	   m = n/mult_factor;  for(j=0;j<m;j++) fprintf(fptr,"=");
	}
	if((n=H->bin0[H->nbins+1]) > 0){
		fprintf(fptr,"\n>=%6.2f : %-8d |",H->max,n); 
		m = n/mult_factor; for(j=0;j<m;j++) fprintf(fptr,"=");
	}
	fprintf(fptr,"\n%8s : %-8d\n\n","total",total);
   if(H->n > 1){
	fprintf(fptr,"    mean = %g\n",MeanHist(H));
	fprintf(fptr,"   stdev = %g\n",sqrt(VarianceHist(H)));
	fprintf(fptr,"   range = %g .. %g\n",H->minval,H->maxval);
	fprintf(fptr,"  median = %g\n\n",MedianHist(H));
   }
}

Int4	*RtnHist(Int4 line_leng,h_type H)
/************************************************************************
	note: Var(X) = E(X*X) - E(X)*E(X).
  returns an array (from min to max) of number of '=' to print...
 ************************************************************************/
{
	Int4	start,total,sum,i,j,n,m, max = 0,mult_factor;
	double	val;
	Int4	*rtn;

	total = 0;
	for(i=0;i<=H->nbins+1;i++) {
		total += H->bin0[i];
		max = H->bin0[i] > max ? H->bin0[i] : max;
	}
	if(max <= line_leng) mult_factor = 1;
	else mult_factor = (max/line_leng) + 1;
	if(total == 0) { return 0; }

	sum = 0;
	if((n=H->bin0[0]) > 0){
		sum = n;
		start = 1;
		val=H->min;
	} else {
	   for(start=0,val=H->min; H->bin0[start]==0 ;start++){
	      val+=H->inc;
	      if(start >= H->nbins) { return 0; }
	   }
	   val-=H->inc;
	}

	NEW(rtn,H->nbins+3,Int4);
	Int4	pos=0;
	for(i=start; i<=H->nbins; val+= H->inc, sum+=H->bin0[i],i++) {
	   n = H->bin0[i];
	   pos++;
	   rtn[pos] = n/mult_factor;  
	   if(sum==total) break;
	}
	if((n=H->bin0[H->nbins+1]) > 0){
	   	pos++; rtn[pos] = n/mult_factor; 
	}
	rtn[0]=pos;
	return rtn;
}


double	MeanHist(h_type H) { return H->total/(double)H->n; }

double	MedianHist(h_type H)
{ 
	Int4	i,j;
	double	val,cut;

	cut = (double)H->n/2.0;
	j=H->bin0[0];
	for(val=H->min,i=0; i<=H->nbins; i++,val+=H->inc){
	      j+=H->bin[i];
	      if((double)j >= cut) return val;
	}
	print_error("MedianHist( ) This should not happen");
	return 0;
}

double	VarianceHist(h_type H)
{
	double	m;
	if(H->n <= 1) return 0.0;
	m = H->total/(double)H->n;
	return (((H->total_sq/(double)H->n) - m*m))*H->n/(double)(H->n-1); 
}

void	IncdMHist(double x, Int4 number, h_type H)
{
	Int4	i;
	if(x < H->min) H->bin0[0] += number;
	else if(x < H->max) { 
		i= (Int4)floor((x-H->min)/H->inc); 
		H->bin[i] += number; 
	} else H->bin[H->nbins]+=number;
	H->total += number*x; H->n += number;
	H->total_sq += number*x*x;
	H->minval = MINIMUM(double,x,H->minval);
	H->maxval = MAXIMUM(double,x,H->maxval);
}

void	IncdHist(double x,h_type H)
{
	Int4	i;

	H->total += x; H->n++;
	H->total_sq += x*x;
	H->minval = MINIMUM(double,x,H->minval);
	H->maxval = MAXIMUM(double,x,H->maxval);
	if(x < H->min){ H->bin0[0]++; }	/*** underflow ***/
	else if(x < H->max) { 		/*** min ... < max ***/
		i= (Int4)floor((x-H->min)/H->inc); 
		H->bin[i]++; 
	} else /** x>=H->max **/ { H->bin[H->nbins]++; } 
}

Int4	IncdHistBin(double x,h_type H)
/** Increment histogram and return the number of the bin incremented **/
{
	Int4	i;

	H->total += x; H->n++;
	H->total_sq += x*x;
	H->minval = MINIMUM(double,x,H->minval);
	H->maxval = MAXIMUM(double,x,H->maxval);
	if(x < H->min){ H->bin0[0]++; return 0; }
	else if(x < H->max) { 
		i= (Int4)floor((x-H->min)/H->inc); 
#if 0
		printf("x=%g; H->bin[%d]=%d\n",x,i,H->bin[i]);
#endif
		H->bin[i]++; 
		return (i+1);
	} else { H->bin[H->nbins]++; return (H->nbins+1); }
}

