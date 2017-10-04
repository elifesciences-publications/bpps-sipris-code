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

#include "evalue.h"

evl_typ MkEvalue(Int4 d_max, Int4 d_min, Int4 k_max, double pmin, double pmax,
        Int4 nsing, Int4 mode, Int4 *N, double *freq, Int4 C_min, a_type A)
{
	evl_typ	E;
	Int4	n,c,c2,k,s,s_min,v,d;
	double	p,pmax2,x,x1,x2,y,y1,y2;
	Int4        time1;

	NEW(E,1,evalue_type);
	pmax2 = pmax; pmax += 3; /****/
	E->d_max = d_max; E->d_min = d_min; E->k_max = k_max;
	E->pmin = pmin; E->pmax = pmax;
	E->mode = mode; E->nsingle = nsing;
	E->freq = freq; E->A = A;
	E->N = N; E->pval0 = NULL;
	E->C_min = C_min;
	NEWPP(E->NP,d_max+1,double);
	for(E->K=0.0,d=d_min;d<=d_max;d++){
	    NEWP(E->NP[d],k_max+1,double);
	    for(k=d; k<=k_max; k++){
	        NEW(E->NP[d][k],d_max+1,double);
		s_min = MINIMUM(Int4,nsing,d);
		for(s=s_min; s<=d; s++){
			E->NP[d][k][s]=numpat_eval(d,k,s,E);
			E->K += E->NP[d][k][s];
		}
		/*****
		for(s=s_min; s<=d; s++){
			V = E->NP[d][k][s]/E->NP[d_min][d_min][d_min];
                        fprintf(stderr,"N[%d][%d][%d] = %-10.0f\n",
				d,k,s,E->NP[d][k][s]);
		}
		/*****/
	    }
	}
        if(mode > 0){
          NEW(E->pfreq,mode+2,double);
          for(d=0.0,n=1,c=1; c<=nAlpha(A);c++){
             for(c2=c+1; c2<=nAlpha(A);c2++){
                if(valAlphaP(c,c2,A) <= mode){
                        E->pfreq[n++] = freq[c]+freq[c2];
                        d+=E->pfreq[n-1];
                }
             }
          }
          /**** fprintf(stderr,"ave = %g\n",d/(double)(n-1));/**/
        } else E->pfreq = NULL;
	E->logK = log10(E->K);
	E->pmin0 = E->pmin - E->logK; 
	E->pmax0 = E->pmax - E->logK;
        E->max = (Int4) (E->pmax0*10.0 +1.0);
        E->min = (Int4) (E->pmin0*10.0 - 1.0);
	E->pval0 = E->pval = NULL;
	time1=time(NULL);
	E->adj_pmax = (double) E->max/10.0;
	for(v=E->max; v >=E->min; v--){
		p = FastEvalue(E, E->adj_pmax);
	/******
		printf("p=%.3f; adj_p = %.3f; v=%d; E[v] = %.3f\n",
			p, E->adj_pmax, v,E->pval[v]);
	/******/
		if(p <= pmax2) break;
		else E->adj_pmax -= 0.1;
	}
	if(E->adj_pmax < (double) E->max/10.0){
		x = pmax2;
		y1 = E->adj_pmax;
		y2 = E->adj_pmax+1.0; 
		x1 = FastEvalue(E,y1);
		x2 = FastEvalue(E,y2);
		y = ((y2-y1)/(x2-x1))*(x-x1) + y1;
		if(y > E->pmax0) y = E->pmax0; 
		E->adj_pmax = y;
	/******
		printf("x =%f\n",x);
		printf("y1 = %f\n",y1);
		printf("y2 = %f\n",y2);
		printf("x1 = %f\n",x1);
		printf("x2 = %f\n",x2);
		printf("y =%f\n",y);
		printf("cutoff = %f (%d seconds)\n",
			E->adj_pmax,time(NULL)-time1);
	/******/
	}
	return E;
}

void	PutEvalue(FILE *fptr, evl_typ E)
{
	Int4	n,d,c,c2;
        for(n=1,c=1; c<=nAlpha(E->A);c++){
             for(c2=c+1; c2<=nAlpha(E->A);c2++){
                if(valAlphaP(c,c2,E->A) <= E->mode){
		     n++;
                     fprintf(fptr,"pfreq[%d] (%c%c): %f\n",n-1,
                        AlphaChar(c,E->A), AlphaChar(c2,E->A),
                        E->pfreq[n-1]);
		}
	     }
	}
	fprintf(fptr,"\n");
	FastEvalue(E,E->pmin+(E->pmax-E->pmin)/2.0);
        for(n=E->min; n<=E->max;n++){
			fprintf(fptr,"%d: %g\n", n, E->pval[n]);
	}
	fprintf(fptr,"\n");
	fprintf(fptr,"pmin = %g\n",E->pmin);
	fprintf(fptr,"pmin0 = %g\n",E->pmin0);
	fprintf(fptr,"pmax = %g\n",E->pmax);
	fprintf(fptr,"pmax0 = %g\n",E->pmax0);
	fprintf(fptr,"min = %d\n",E->min);
	fprintf(fptr,"max = %d\n",E->max);
}

void    NilEvalue(evl_typ E)
{
	Int4	d,k;

	for(d=E->d_min;d<=E->d_max;d++){
	    for(k=d; k<=E->k_max; k++) free(E->NP[d][k]);
	    free(E->NP[d]);
	}
	free(E->NP);
	if(E->pval0 != NULL) free(E->pval0);
	if(E->pfreq != NULL) free(E->pfreq);
	free(E);
}

void	eval_error(char *s){fprintf(stderr,"Evalueerate: %s\n",s);exit(1);}

double  numpat_eval(Int4 d, Int4 k, Int4 s, evl_typ E)
/* returns the # of patterns of length k, depth d and s 1-residue sets */
{
        double  N,n1,n2;

        n1 = (double) nAlpha(E->A);
        N = bico((k-2),(d-2));		/* pick positions */
        N *= pow(n1,(double)(s));
	if(d > s){
        	n2 = (double) E->mode;
        	N *= bico(d,s);     	/* pick singles */
        	N *= pow(n2,(double)(d-s));
	} else if(d < s) print_error(" d < s in NumPatEvalue( )");
        return N; 
} 

double  FastEvalue(evl_typ E, double prob)
/***************************************************************
longerpolate E-value for pattern prob from array of E-values.

		*<-- t
	       /|         t-b      r-b
	      / |        ------ = -----  --> t-b = (r-b)/(p-x+1)
	     /  |       x-(x-1)   p-(x-1)
	    /...|<-- r
	   /:   |
	  / :   |        --> r = (t-b)*(p-x+1)+b
	 /  :   |       
	*---+---+<-- b
       x-1  p   x

 **************************************************************/
{
	Int4	x;
	double	t,b,r,p;

	if(E->pval0 == NULL) estprob_eval(10000, E);
	p = prob*10.0;
	x = (Int4)ceil(p);
	if(x <= E->min || x > E->max) return EstEvalue(10000,E,prob);
	else {
		t = E->pval[x];
		b = E->pval[x-1];
		r = (t-b)*(p-(double)x+1)+b;
		return r;
	}
}

Int4	sample_pttrn_type(Int4 *D, Int4 *K, Int4 *S, register evl_typ E)
/** sample a d,k,s-pattern type **/
{
	register double  R,***NP = E->NP; 
	register Int4     d,k,s;

        R = E->K*((double) Random()/(double) RANDOM_MAX);
	for(d=E->d_max; d >= E->d_min; d--){
	    s = MINIMUM(Int4,E->nsingle,d);
	    for( ; s <= d; s++){
		for(k=E->k_max; k >= d; k--){
		   R -= NP[d][k][s];
		   if(R <= 0.0){ *K=k; *D=d; *S=s; return E->N[k]; }
		}
	    }	
	}	
	eval_error("This should not happen!");
}

double	sample_pttrn_freq(register Int4 d, register Int4 s, 
					register evl_typ E)
/*************************************************************
Sample a pattern frequency for input into CumBinomProb function
from d,k,s-pattern type.
 *************************************************************/
{
	register double  p,R,*freq; 
	register Int4     x,nres;

        for(p=1.0; d > 0; d--) {
        	if(d <= s){ nres = nAlpha(E->A); freq = E->freq;}
        	else { nres = E->mode; freq = E->pfreq;}
        	do{             
        		R = (double) Random()/(double) RANDOM_MAX;
        	} while((x = (Int4)(R*(double) nres)) >= nres);
        	p *= freq[x+1];
        }	
	return p;
}

void	estprob_eval(Int4 ntest, evl_typ E)
{
   double  p,prob2,*prob,*prob0,*sum,*sum0,ave;
   Int4     v,N,i,min,max,d,k,s,num;

    min = E->min; max = E->max;
    if((max - min) <= 0) print_error("EstProb: input error");
    num = max - min +1;
    NEW(prob0,num+2,double); NEW(sum0,num+2,double); 
    prob = prob0 - min; sum = sum0 - min; 
    for( v=min; v <= max; v++){
	sum[v] = 0.0;
	prob[v] = pow(10.0,(double)v/10.0);
        prob[v] *= 1.00001;        /* AVOID ROUNDING ERRORS */
    }
    for(i=0; i< ntest; i++){
	N=sample_pttrn_type(&d, &k, &s, E);
	p=sample_pttrn_freq(d,s,E);
        for(v=max, k=E->C_min; k <= N; k++){
                prob2 = CumBinomProb(k,N,p); 
		while(v >= min && prob2 <= prob[v]){
                        sum[v] += prob2; v--;
                }
		if(v < min) break;
        }
    }
    for(v = min; v <= max ; v++){
    	ave = sum[v]/(double) ntest;
    	prob[v] = prob[v] + (E->K-1.0) * ave;
	prob[v] = log10(prob[v]);
    }
    E->pval0 = prob0; E->pval = prob;
    free(sum0);
}

double EstEvalue2(Int4 ntest, evl_typ E,double prob)
{
   double  tot_p,p,prob2,ave,sum;
   Int4     N,i,x,d,k,s,s_min;

    ntest /=100;
    prob = pow(10.0,prob);
    prob *= 1.00001;        /* AVOID ROUNDING ERRORS */
    for(tot_p=0.0,d=E->d_max; d >= E->d_min; d--){
	for(k=E->k_max; k >= d; k--){
	    N = E->N[k];
	    s_min = MINIMUM(Int4,E->nsingle,d);
	    for(s=s_min; s <= d; s++){
    		for(sum=0.0, i=0; i< ntest; i++){
		   p=sample_pttrn_freq(d,s,E);
        	   for(x=E->C_min; x <= N; x++){
                	prob2 = CumBinomProb(x,N,p);
                	if(prob2 <= prob){
                       		sum += prob2;
                       		break;
                	}
        	   }
		}
		ave = sum/(double) ntest;
		tot_p += (E->NP[d][k][s]) * ave;
	    }
	}
    }
    tot_p += prob; 
    return log10(tot_p);
}

double EstEvalue(Int4 ntest, evl_typ E,double prob)
{
   double  p,prob2,ave,sum;
   Int4     N,i,d,k,s;

    prob = pow(10.0,prob);
    prob *= 1.00001;        /* AVOID ROUNDING ERRORS */
    for(sum=0.0, i=0; i< ntest; i++){
	N=sample_pttrn_type(&d, &k, &s, E);
	p=sample_pttrn_freq(d,s,E);
        for(k=E->C_min; k <= N; k++){
                prob2 = CumBinomProb(k,N,p);
                if(prob2 <= prob){
                        sum += prob2;
                        break;
                }
        }
	
    }
    ave = sum/(double) ntest;
    p = prob + (E->K-1.0) * ave;
    return log10(p);
}


