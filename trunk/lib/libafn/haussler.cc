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

#include "haussler.h"
#include "dirichlet.h"

#define MAXDM	200
	
dmp_typ MkDMPriors(BooLean full)
{
	FILE	*fptr;
	dmp_typ	plib;
	char	*dm;
	static char    aa[22] = {"XACDEFGHIKLMNPQRSTVWY"};
	char	i;

	if(full) dm=Dirichlet;
	// if(full) dm=Dirichlet20;
	else dm = Dirichlet9;
	fptr = tmpfile();
	fprintf(fptr,"%s",dm); rewind(fptr);
	plib = read_DMPriors(fptr);
	fclose(fptr);
	for(i=1;i<=20;i++) plib->let2code[aa[i]] = i;
	return plib;
}

void    NilDMPriors(dmp_typ D) { free_DMPriors(D); }

void    CalcDMPriors(double *counts, double *pseudo, dmp_typ D,
	a_type A)
{
	Int4	i,r,i2;
	double	sum,*pseudo_counts,cnts[22];
	static char    aa[22] = {"XACDEFGHIKLMNPQRSTVWY"};

	if(nAlpha(A) != 20) print_error("incompatible alphabet");
	for(i=1; i<=20; i++){
		r = AlphaChar(i,A); i2 = D->let2code[r];
		cnts[i2] = counts[i];
	}
	pseudo_counts = mixture_regularizer(cnts, D);
	for (sum=0,i=1;i<=20;i++){
		i2 = AlphaCode(aa[i],A);
		sum+=pseudo[i2]=cnts[i]+7.5*pseudo_counts[i];
#if 0
		sum+=pseudo[i2]=pseudo_counts[i]; /** NEW method **/
		sum+=pseudo[i2]=cnts[i]+pseudo_counts[i];
#endif
	}
	for(i=1;i<=20;i++) pseudo[i]/=sum;
}

void    CalcDMPriors2(Int4 *counts, double *pseudo, dmp_typ D,
	a_type A)
{
	Int4	i,r,i2;
	double	sum,*pseudo_counts,cnts[22];
	static char    aa[22] = {"ACDEFGHIKLMNPQRSTVWY"};

	if(nAlpha(A) != 20) print_error("incompatible alphabet");
	for(i=1; i<=20; i++){
		r = AlphaChar(i,A); i2 = D->let2code[r];
		cnts[i2] = (double) counts[i];
	}
	pseudo_counts = mixture_regularizer(cnts, D);
	for (sum=0,i=1;i<=20;i++){
		i2 = AlphaCode(aa[i],A);
		sum+=pseudo[i2]=cnts[i]+pseudo_counts[i];
	}
	for(i=1;i<=20;i++) pseudo[i]/=sum;
}

DMPriors *alloc_DMPriors(Int4 l, Int4 Alpha)
/************************************************************************
        L       Number of distributions
        Alpha   Number of alphabet characters
  Allocate space for a prior library with L priors and Alpha characters.
 ************************************************************************/
{
        dmp_typ	temp;
	Int4	i;

        temp = (DMPriors *)malloc(sizeof(DMPriors));
        temp->L = l;
        temp->AlphaChar = Alpha;

        temp->Mix = (double *)malloc(sizeof(double)*l);
        temp->FullUpdate = (Int4 *)malloc(sizeof(Int4)*l);
        temp->QUpdate = (Int4 *)malloc(sizeof(Int4)*l);

        temp->StructID = (char **)malloc(sizeof(char *)*l);
        temp->Comment = (char **)malloc(sizeof(char *)*l);
        temp->Distr = (double **)malloc(sizeof(double *)*l);
        for (i=0; i<l; i++) {
           temp->Distr[i] = (double *)malloc(sizeof(double)*(Alpha+1));
           temp->StructID[i] = (char *)malloc(sizeof(char)*MAXDM);
           temp->Comment[i] = (char *)malloc(sizeof(char)*MAXDM);

        } /* endfor */

        return(temp);
}

void free_DMPriors (DMPriors *lib)
{
  Int4 i;

  if(lib==NULL) return;
  for(i=0; i<lib->L; i++) {
    free(lib->Distr[i]); free(lib->StructID[i]); free(lib->Comment[i]);
  }
  free(lib->Mix); free(lib->FullUpdate); free(lib->QUpdate);
  free(lib->StructID); free(lib->Comment); free(lib->Distr);
  free(lib);
}

void	write_DMPriors(FILE *fp, DMPriors *Library)
/*************************************************************************
   fp      File to write information to Library Information to write 
 *************************************************************************/
{
        Int4 i,j;

        fprintf(fp,"AlphaChar= %d\n",Library->AlphaChar);
        fprintf(fp,"NumDistr= %d\n",Library->L);

        for(i=0; i<Library->L; i++) {
          fprintf(fp,"Number= %d\n",i);
          fprintf(fp,"Mixture= %lf\n", Library->Mix[i]);
          fprintf(fp,"Alpha= ");
          for(j=0; j<=(Library->AlphaChar); j++)
                        fprintf(fp,"%lf ", Library->Distr[i][j] );
          fprintf(fp,"\n");
          fprintf(fp,"FullUpdate= %d\n",Library->FullUpdate[i]);
          fprintf(fp,"QUpdate= %d\n", Library->QUpdate[i]);
          fprintf(fp,"StructID= %s\n", Library->StructID[i]);
          fprintf(fp,"Comment= %s\n", Library->Comment[i]);
        }
}

DMPriors *read_DMPriors(FILE *fp)
/*********************************************************************
   fp      File to write information to Library Information to write 

 This reads prior information from a file longo a DMPriors.
 *********************************************************************/
{
        Int4     i,j,Alpha, l;
        dmp_typ	temp;
        char    input[MAXDM], foo[MAXDM];
        double	x;

        fscanf(fp,"%*s %d\n", &Alpha);
        fscanf(fp,"%*s %d\n", &l);
        temp=alloc_DMPriors(l,Alpha); temp->AlphaChar=Alpha; temp->L=l;
        for(i=0; i<temp->L; i++) {
        	/* Get rid of number= */
        	fscanf(fp,"%*s %*s\n");

        	/* Mixture */
        	fscanf(fp,"%*s");
        	fscanf(fp,"%lf\n", &x);
        	temp->Mix[i] = (double)x;

        	/* Alpha */
        	fscanf(fp,"%*s");
        	for(j=0; j <= temp->AlphaChar; j++) {
               		fscanf(fp,"%lf", &x); temp->Distr[i][j] = x;
        	}

        	/* FullUpdate */
        	fscanf(fp,"%*s");
        	fscanf(fp,"%d\n", &(temp->FullUpdate[i]));

        	/* QUpdate */
        	fscanf(fp,"%*s");
        	fscanf(fp,"%d\n", &(temp->QUpdate[i]));

        	/* StructID */
        	fgets(input, MAXDM, fp);
        	sscanf(input,"%s",foo);
        	strcpy( (temp->StructID[i]), (input + strlen(foo)) );

        	/* Comments */
        	fgets(input, MAXDM, fp);
        	sscanf(input,"%s",foo);
        	strcpy( (temp->Comment[i]), (input + strlen(foo)) );
        }
        return(temp);
}

double	*mixture_regularizer(double *counts, dmp_typ Lib)
/******************************************************************
	counts	 observed residue counts
	Lib	 priors

  (see eqn. (4) in Brown et al. 1993)
  This calculates the regularizer given the observed counts, the 
  prior library, and a weight of the priors
  Old version.
******************************************************************/
{
	double	*n, sum, *prob;
	Int4 	i,j;

	/* prob[0..19] */
	MEW(prob,(Lib->AlphaChar+2),double);
	
	/* Put counts longo array with n[0] = sum n_i */
	MEW(n,(Lib->AlphaChar+2),double);
	for(sum=0.0,i=1; i<=Lib->AlphaChar; i++){ 
		sum+=counts[i]; n[i]=counts[i];
	} n[0]=sum;
	logpajgy(n, Lib, 0, TRUE);  /* Calculate probs */

	/* Calculate new regularizer */
	for(prob[0]=0.0, i=1; i<=Lib->AlphaChar; i++) {
	   for(prob[i]=0.0,j=0; j< Lib->L; j++){
		prob[i] += (exp(logpajgy(n, Lib, j, FALSE)))*
		      		((Lib->Distr[j])[i]); 
	   }
	} free(n);
	return(prob);
}
double	*mixture_regularizer2(double *counts, dmp_typ Lib)
/******************************************************************
  New version (afn implementation).
	counts	 observed residue counts
	Lib	 priors

  (see eqn. (42) in Sjolander et al. 1996)
  This calculates the regularizer given the observed counts, the 
  prior library, and a weight of the priors
******************************************************************/
{
	double	*n, sum, tmp,*prob,max,*logB;
	Int4 	i,j;

	/* prob[0..19] */
	MEW(prob,(Lib->AlphaChar+2),double);
	MEW(logB,(Lib->L+2),double);
	
	/* Put counts longo array with n[0] = sum n_i */
	MEW(n,(Lib->AlphaChar+2),double);
	for(sum=0.0,i=1; i<=Lib->AlphaChar; i++){ 
		sum+=counts[i]; n[i]=counts[i];
	} n[0]=sum;

	for(max=-DBL_MAX,j=0; j< Lib->L; j++){
	    logB[j] = log_beta2_dmp(Lib->Distr[j], n)
                                - log_beta1_dmp(Lib->Distr[j]);
	    max = MAXIMUM(double,max,logB[j]);
	}
	for(j=0; j< Lib->L; j++) logB[j] -= max;
	
	/* Calculate new regularizer */
	for(prob[0]=0.0, i=1; i<=Lib->AlphaChar; i++) {
	   for(prob[i]=0.0,j=0; j< Lib->L; j++){
		tmp = Lib->Mix[j]*exp(logB[j]);
		tmp *= (Lib->Distr[j][i] + n[i]);
		tmp /= (Lib->Distr[j][0] + n[0]);
		prob[i] += tmp;
		prob[0] += tmp;
	   }
	}
	for(i=1; i<= Lib->AlphaChar; i++) prob[i] /= prob[0];
	free(n); free(logB);
	return(prob);
}

double	log_beta1_dmp(double *n)
{
	Int4	i;
	double	v=-lngamma(n[0]);

	for(i=1; i<=AlphLength; i++) v+=lngamma(n[i]); 
	return v;
}

double	log_beta2_dmp(double *a, double *n)
{
	Int4	i;
	double	v=-lngamma(n[0]+a[0]);

	for(i=1; i<=AlphLength; i++) v+=lngamma(n[i]+a[i]); 
	return v;
}

double LogAddLog(double x, double y)
/**************************************************************************
 This function returns log(e^x + e^y) in a nice way. It is used to
  prevent overflow errors when dealing with large x and y.
 **************************************************************************/
{
	if (x>y) { return(x + log(1.0+exp(y-x))); }
	else { return(y + log(1.0+exp(x-y))); }
}

double	logpajgy(double	*n, dmp_typ Library, Int4 j, BooLean calc)
/**************************************************************************
   This function computes log(p(a^j|y)) used in the calculation of theta. 
   It is defined to be

    p(j|n) = q_j x p(n|a^j) divided by sum (q_k x p(n|a^k))
                                        k

  (Using Bayes Rule: see Eqn. (3) in Brown et al., 1993)
  where n = obs_freq.  The j'th prior is examined.
 **************************************************************************/
{
	Int4 i;
	double	tmp;
	static double	logprob[MAXDM], logdenom;  /* Holders for probabilities */

	/* calculate log probs if not already done */
	if(calc){	/* if calc=TRUE calculate all probs for denominator. */
	   tmp=log(Library->Mix[0]) + logpygaj(n,Library->Distr[0]);
	   logdenom = tmp; logprob[0] = tmp;

	   for(i=1; i<Library->L; i++) {	/* Do remaining terms */
		tmp=log(Library->Mix[i]) + logpygaj(n,Library->Distr[i]);
		logdenom = LogAddLog(logdenom, tmp); logprob[i] = tmp;
	   }
	}
	return(logprob[j] - logdenom);
}

double	logpygaj(double *n, double *a)
/**************************************************************************
   This function computes log(p(n|p_j)) used in the calculation of theta. 
   It is defined to be
                                           20
                  / Gamma(n+1)Gamma(a) \  ---  /     Gamma(n_i+a_i)     \
    p(n|a) = log{| ---------------------| | | | ------------------------ |}
                  \    Gamma(n + a)    /  i=1  \ Gamma(n_i+1)Gamma(a_i) /

  (see Eqn. ("3b") in Brown et al., 1993)
  where n = obs_cnts and a[i]=Dirichlet distribution parameters.
 **************************************************************************/
{
	Int4	i;
	double	temp;
	
	temp = lngamma(n[0]+1.0) + lngamma(a[0]); 
	temp -= lngamma(n[0]+a[0]);
	for(i=1; i<=AlphLength; i++){
		temp+=lngamma(n[i]+a[i]); 
		temp -= lngamma(n[i]+1.0) + lngamma(a[i]);
	}
	return(temp);
}

