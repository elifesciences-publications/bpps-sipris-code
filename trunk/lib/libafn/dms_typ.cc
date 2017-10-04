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

#include "dms_typ.h"
#include "dms_data.h"
#include "blosum62.h"

/*	Dirichlet-mixture score code for proteins			*/
/*	Version 1.1							*/
/*	April 11, 2014							*/
/*	Program by Stephen Altschul					*/
/*									*/
/*	Reference:  Sjolander, K., et al. (1996) "Dirichlet mixtures:	*/
/*	a method for improved detection of weak but significant protein	*/
/*	sequence homology." Comput. Appl. Biosci. 12:327-345.		*/

/*	Dirichlet mixture "recode5" from:				*/
/*	http://compbio.soe.ucsc.edu/dirichlets/index.html		*/
/*									*/
/*	Assumed amino acid order: ARNDCQEGHILKMFPSTWYV			*/

/*	Call this program once, to initialize Dirichlet Mixture data, and to	*/
/*	base the amino acid background frequencies on the Dirichlet mixture.	*/

void	dms_typ::Free()
{
	for(Int4 m=1; m <= nDC; ++m) { free(alpha[m]); }
	free(alpha); free(Alpha); free(weight);
	free(countX); free(count); free(back); NilAlpha(AB); 
	free(offset);
}

double	*dms_typ::RtnData1(double *data,a_type xAB)
{
	char	r;
	Int4	i,j;
	double	pseudo=20;
	assert(nAlpha(xAB) == 20); data[0]=1.0;
	for(i=1;i <= nAlpha(xAB); ++i){
		r=AlphaChar(i,xAB); j=AlphaCode(r,AB); data[i]=pseudo*blosum62freq[j];
	} return data;
}

void	dms_typ::init(char mode, a_type ab)
{
	int	r,i,m;
	double	*ptrdata;		/*  Pointer to DM data			*/
	double	data_one[50];
	char	c;
	a_type xAB=MkAlpha("XARNDCQEGHILKMFPSTWYV",PROT_BLOSUM62);  // WARNING!!! PROT_BLOSUM62 is invalid!
	AB=CopyAlphabet(ab);
	switch (mode){
	   case 'O': ptrdata=RtnData1(data_one,xAB); nDC=1; break;	// one component only.
	   case 'T': ptrdata = ucsc_data; nDC=20; break;	// Tiny = twenty = 20.
	   case 'f': ptrdata = data_52; nDC=52; break;		// fifty-two = = 52.
	   case 'F': ptrdata = data_58; nDC=58; break;		// Fifty-eight= 58.
	   case 'S': ptrdata = data_72; nDC=72; break;		// small = seventy-two = 72.
	   case 'M': ptrdata = data_110; nDC=110; break;	// medium = 110.
	   case 'L': ptrdata = data_134; nDC=134; break;	// large = 134.
	   default: ptrdata = data_20; nDC=20; break;		// Altschul's 20-component mixture.
	} NEW(weight,nDC+2,double); NEW(Alpha,nDC+2,double); NEWP(alpha,nDC+2,double);
	assert(nAlpha(AB) == 20); assert(nAlpha(xAB) == 20);
	NEW(count,nAlpha(AB)+2, double);
	NEW(countX,nAlpha(AB)+2, double);
	for(m=1; m <= nDC; ++m) {
		NEW(alpha[m],nAlpha(AB)+2, double);
		weight[m]= *ptrdata++;
		Alpha[m]=0;
		for(i=1;i <= nAlpha(xAB); ++i) Alpha[m]+=alpha[m][i] = *ptrdata++;
	}
	// convert alpha[m][i] over to global alphabet...
	double *tmp; NEW(tmp,nAlpha(AB)+2,double);
	for(m=1; m <= nDC; ++m) {
	   for(r=1; r <= nAlpha(AB); ++r) tmp[r]=0;
	   for(i=1; i <= nAlpha(xAB); ++i){
		c=AlphaChar(i,xAB); r=AlphaCode(c,AB);
		assert(tmp[r]==0); tmp[r]=alpha[m][i];	// i --> r mapping.
	   }
	   for(r=1; r <= nAlpha(AB); ++r){ assert(tmp[r] != 0); alpha[m][r]=tmp[r]; }
	} free(tmp);

	NEW(back,nAlpha(AB)+2, double);
	for(r=1; r <= nAlpha(AB); ++r){
	   for(m=1; m <= nDC; ++m){
		back[r]+=weight[m]*alpha[m][r]/Alpha[m];
// std::cerr << AlphaChar(r,AB) << "." << m << " = " << weight[m]*alpha[m][r]/Alpha[m] << std::endl;
	   }
// std::cerr << AlphaChar(r,AB) << " = " << back[r] << std::endl;
	}
	NEW(offset,nDC+2,double);	// speed up for code...
        for(m=1; m <= nDC; m++){
             offset[m]=log(weight[m])+lgamma(Alpha[m]);
             for(r=1; r <= nAlpha(AB); ++r){ offset[m]-=lgamma(alpha[m][r]); }
        } NilAlpha(xAB);
}


/*	Call this program once for each column, with a vector of independent	*/
/*	counts for the 20 amino acids, in this order:   ARNDCQEGHILKMFPSTWYV	*/
/*	The program returns a doubleing point vector of scores, with the units	*/
/*	of bits.  This vector may scaled as one chooses.			*/

// double	lgamma(),log(),exp();

// void  dms_typ::score(double count[],double score[])
void  dms_typ::score(UInt4 *kount,double *score)
        // double	count[];		/*  Counts of indep. a.a. observations	*/
	// double	score[];		/*  Amino acid scores, in bits		*/
{
	int	i,m;
	double	Count;			/*  Aggregate number of observations	*/
	double	temp;
	double	sum;
	double	q;			/*  Posterior for the amino acid	*/
	double	pr[nDC+2];		/*  Posteriors of the components	*/
	
	//  Calculate posteriors for the Dirichlet components.
	for(Count=0.0,i=1; i <= nAlpha(AB); ++i){
		count[i]=(double) kount[i]/(double) wt_factor; Count+=count[i]; 
// std::cerr << "count(" << AlphaChar(i,AB) << ") = " << count[i] << std::endl;
	}
	for(m=1; m <= nDC; ++m) {
		temp=log(weight[m])+lgamma(Alpha[m])-lgamma(Alpha[m]+Count);
		for(i=1; i <= nAlpha(AB); ++i)
			temp+=lgamma(alpha[m][i]+count[i])-lgamma(alpha[m][i]);
		if(m==1) sum=temp;
		else if (sum>temp) sum+=log(1.0+exp(temp-sum));
		else sum=temp+log(1.0+exp(sum-temp));
		pr[m]=temp;
	}
	for (m=1; m <= nDC; ++m) pr[m]=exp(pr[m]-sum);

	//  Calculate log-odds scores for the amino acids, in bits.
	for (i=1; i <= nAlpha(AB); ++i) {
		temp=count[i];
		for (q=0,m=1; m <= nDC; ++m) q+=pr[m]*(alpha[m][i]+temp)/(Alpha[m]+Count);
		// score[i]=log(q/back[i])/log(2.0);
		score[i]=pernats*log(q/back[i]);
// std::cerr << "score(" << AlphaChar(i,AB) << ") = " << score[i] << std::endl;
	}
// exit(1);
}

/*	BILD information-content code for proteins			*/
/*	Version 1.2							*/
/*	May 27, 2014							*/
/*	Program by Stephen Altschul					*/
/*									*/
/*	Reference:  Altschul, S.F., Wootton, J.C., Zaslavsky, E. &	*/
/*	Yu, Y.-K. (2010) "The construction and use of log-odds		*/
/*	substitution scores for multiple sequence alignment."		*/
/*	PLoS Comput. Biol. 6:e1000852.					*/

/*	Dirichlet mixture data						*/

#if 0
/*	Call this program once, to initialize Dirichlet Mixture data, and to	*/
/*	base the amino acid background frequencies on the Dirichlet mixture.	*/

void	bildinit(back)
	double	back[];			/*  Amino acid background frequencies	*/
{
	int	i,m;
	double	*ptrdata;		/*  Pointer to DM data			*/

	ptrdata = data;
	for (m=0;m< nDC;++m) {
		weight[m]= *ptrdata++;
		Alpha[m]=0;
		for (i=0;i<20;++i) Alpha[m]+=alpha[m][i]= *ptrdata++;
	}

	for (i=0;i<20;++i) for (back[i]=m=0;m<M;++m)
		back[i]+=weight[m]*alpha[m][i]/Alpha[m];
}
#endif


double  dms_typ::bild(UInt4 *kount)
/*	Call this program once for each column, with a vector of independent	*/
/*	counts for the 20 amino acids, in this order:   ARNDCQEGHILKMFPSTWYV	*/
/*      as well as a vector of background frequencies in the same order.  The	*/
/*	program returns a floating point column score, with the unit of nats.	*/
{
	int	i,m;
	double	Count;			/*  Aggregate number of observations	*/
	double	info;			/*  BILD score in nats			*/
	double	logprob;		/*  Log-likelihood of data		*/
	
/*  Calculate aggregate counts for the column					*/

#if 1	// add in background counts for each deleted '-' residue.
	// fprintf(stderr,"Counts[X] = %d\n",kount[0]);
	for(Count=0.0,i=1; i <= nAlpha(AB); ++i){ 
		count[i]=back[i]*((double) kount[0]/(double) wt_factor);
		count[i] += (double) kount[i]/(double) wt_factor; Count+=count[i]; 
	} if(Count<=1.0) return(0.0);
#else	// use this to debug hieraln...
	for(Count=0.0,i=1; i <= nAlpha(AB); ++i){ 
		count[i]=(double) kount[i]/(double) wt_factor; Count+=count[i]; 
	} if(Count<=1.0) return(0.0);
#endif

//  Calculate BILD score in nats						*/

#if 1	// New, faster code...
         for (m=1; m <= nDC; ++m) {
                 logprob=offset[m]-lgamma(Alpha[m]+Count);
                 for(i=1; i <= nAlpha(AB); ++i){
                         logprob+=lgamma(alpha[m][i]+count[i]);
                 }
                 if(m==1) info=logprob;
                 else if(info > logprob) info+=log(1.0+exp(logprob-info));
                 else info=logprob+log(1.0+exp(info-logprob));
         }
#else
	for (m=1; m <= nDC; ++m) {
		logprob=log(weight[m])+lgamma(Alpha[m])-lgamma(Alpha[m]+Count);
		for(i=1; i <= nAlpha(AB); ++i){ 
			logprob+=lgamma(alpha[m][i]+count[i])-lgamma(alpha[m][i]);
		}
		if(m==1) info=logprob;
		else if(info > logprob) info+=log(1.0+exp(logprob-info));
		else info=logprob+log(1.0+exp(info-logprob));
	}
#endif
	for (i=1; i <= nAlpha(AB); ++i) { info-= count[i]*log(back[i]); }
	// return(info*pernats);
	return(info);
}

#if 0	// Stephen's suggested speed up for this code 
Looking over your BILD score code, adapted from code I had sent you, 
I realize it is rather inefficient, requiring twice as many calls to 
lgamma as are necessary.  It can be improved by making the changes suggested below.

Run once, at end of initialization routine, to calculate offset[m];

         for (m=1; m <= nDC; ++m) {
                 offset[m]=log(weight[m])+lgamma(Alpha[m]);
                 for(i=1; i <= nAlpha(AB); ++i){ offset[m]-=lgamma(alpha[m][i]); }
         }

Run subsequently, to calculate BILD scores:

         for (m=1; m <= nDC; ++m) {
                 logprob=offset[m]-lgamma(Alpha[m]+Count);
                 for(i=1; i <= nAlpha(AB); ++i){
                         logprob+=lgamma(alpha[m][i]+count[i]);
                 }
                 if(m==1) info=logprob;
                 else if(info > logprob) info+=log(1.0+exp(logprob-info));
                 else info=logprob+log(1.0+exp(info-logprob));
         }
#endif


double	dms_typ::ComputeLogLike(double *obs,FILE *efp)
#if 0	//****************************************************************
	Note: Have log a and log b and want log(a + b).
	log(a + b) = log a *(1+b/a) = log a + log (1 + b/a) 
		   = log a + log (1 + e(log b - log a)).
#endif	//****************************************************************
{
	Int4	x,b;
        double  logprob,total,tPrb;
	for(total=0.0,b=1;b <= nAlpha(AB); b++) { total += obs[b]; }
        for(tPrb=0.0,x=1;x<=this->nDC; x++){
                 logprob=offset[x]-lgamma(Alpha[x]+total);
                 for(b=1; b <= nAlpha(AB); b++){ logprob += lgamma(alpha[x][b]+obs[b]); }
                 if(x==1) tPrb=logprob;
                 else if(tPrb > logprob) tPrb+=log(1.0+exp(logprob-tPrb));
                 else tPrb=logprob+log(1.0+exp(tPrb-logprob));
                 if(efp) fprintf(efp,"   tPrb[%d] = %.2f; lgprob=%.2f; total=%.2f\n",x,tPrb,logprob,total);
        } if(efp) fprintf(efp,"   tPrb = %.2f\n",tPrb);
#if 0
if(fabs(tPrb) > 9999999){
	for(b=1; b <= nAlpha(AB); b++){
	   fprintf(stderr," obs[%d]=%.2f.\n",b,obs[b]);
	}
	for(x=1;x<=this->nDC; x++){
	   for(b=1; b <= nAlpha(AB); b++){
		fprintf(stderr," alpha[%d][%d]=%.2f.\n",x,b,alpha[x][b]);
	   }
	}
        fprintf(stderr,"   tPrb[%d] = %.2f; lgprob=%.2f; total=%.2f\n",x,tPrb,logprob,total);
	assert(fabs(tPrb) < 9999999);
}
#endif
	return tPrb;
}

double  dms_typ::sdotfs(UInt4 *kount1,UInt4 *kount2)
#if 0
        double	count1[];		/*  Counts of indep. a.a. observations	*/
        double	count2[];		/*  Counts of indep. a.a. observations	*/
	double	back[];			/*  Background amino acid frequencies	*/
	Call this program once for each column, with a vector of independent	
	counts for the 20 amino acids, in this order:   ARNDCQEGHILKMFPSTWYV	
	The program returns a floating point vector of scores, with the units	
	of bits.  This vector may scaled as one chooses.			
#endif
{
	int	i,m;
	double	Count1,Count2;		//  Aggregate number of observations
	double	temp,sum;
	double	q;			//  Posterior for the amino acid
	double	*pr;			//  Posteriors of the components
	double	score1,score2;

	NEW(pr,nDC+3,double);
	for (Count1=Count2=0,i=1; i<=nAlpha(AB); i++) {
		count[i]=(double) kount1[i]/(double) wt_factor; Count1 += count[i];
		countX[i]=(double) kount2[i]/(double) wt_factor; Count2 += countX[i];
	}
	if(Count1 == 0 || Count2 == 0){ free(pr); return 0; }

	// 1. Calculate component posteriors implied by the first column.
	for (m=1; m <= nDC; ++m) {
		temp=log(weight[m])+lgamma(Alpha[m])-lgamma(Alpha[m]+Count1);
		for (i=1; i<= nAlpha(AB); ++i){
			temp+=lgamma(alpha[m][i]+count[i])-lgamma(alpha[m][i]);
		}
		if(m==1) sum=temp;
		else if (sum>temp) sum+=log(1.0+exp(temp-sum));
		else sum=temp+log(1.0+exp(sum-temp));
		pr[m]=temp;
	} for (m=1; m <= nDC; ++m) pr[m]=exp(pr[m]-sum);

	// 2. Calculate weighted log-odds scores.
	for (score1=0,i=1; i <= nAlpha(AB); ++i) {
		temp=count[i];
		for (q=0,m=1; m <= nDC; ++m) q+=pr[m]*(alpha[m][i]+temp)/(Alpha[m]+Count1);
		score1+=log(q/back[i])*countX[i];
	}

	// 3. Calculate component posteriors implied by the second column.
	for (m=1; m <= nDC; ++m) {
		temp=log(weight[m])+lgamma(Alpha[m])-lgamma(Alpha[m]+Count2);
		for (i=1; i<= nAlpha(AB); ++i){
			temp+=lgamma(alpha[m][i]+countX[i])-lgamma(alpha[m][i]);
		} if(m==1) sum=temp;
		else if (sum>temp) sum+=log(1.0+exp(temp-sum));
		else sum=temp+log(1.0+exp(sum-temp));
		pr[m]=temp;
	}
	for (m=1; m <= nDC; ++m) pr[m]=exp(pr[m]-sum);

	// 4. Calculate weighted log-odds scores.
	for (score2=0,i=1; i <= nAlpha(AB); ++i) {
		temp=countX[i];
		for (q=0,m=1; m <= nDC; ++m) q+=pr[m]*(alpha[m][i]+temp)/(Alpha[m]+Count2);
		score2+=log(q/back[i])*count[i];
	} free(pr);
#if 0
	// Return average log-odds score for the two columns, in bits	
	return((score1/Count2+score2/Count1)/log(4.0));
#elif 0
	score1*=log(Count1)/Count2;
	score2*=log(Count2)/Count1;
	return( (score1+score2) / ( (log(Count1)+log(Count2)) * log(4.0) ) );
#else
	// Return average log-odds score for the two columns, in bits	
	return((score1/Count2)+(score2/Count1));
#endif
}

#if 0	// From Stephen:

   The scoring system gives equal weight to the dot product of the "profile scores" of column 1 with the frequencies of column 2, and the dot product of the profile scores of column 2 with the frequencies of column 1.  It does not deal with scores for null characters, so this is a matter for experimentation.

   There are various other ways one can modify this basic scoring system.  If you want to give a greater weight the profile scores from larger columns, as we discussed, you can modify the scoring system, for example, as follows:

Replace

return((score1/Count2+score2/Count1)/log(4.0));

with

score1*=log(Count1)/Count2;
score2*=log(Count2)/Count1;
return( (score1+score2) / ( (log(Count1)+log(Count2)) * log(4.0) ) );

This has the advantage or reducing to a profile-sequence alignment score when one of the columns has only a single observation.  You need to watch out, however, for the case where both columns have only a single observation, which you can deal with by giving equal weight to both scores.

   A separate issue might be that one wants to give lesser weight to the overall score of the alignment of two columns with few observations than to the alignment of two columns with many observations.  I have not looked into this.  There are several obvious things one could try, and you might want to experiment with some of them.

#endif

/*****************************************************************************
   Below is code to calculate a "typicality index", in bits, for a column of 
amino acid observations.  The code assumes a Dirichlet mixture prior and set 
of background frequencies, usually calculated from the Dirichlet mixture.  
I believe you already have code to populate these numbers.  You can of course 
provide these data to the subroutine by passing in variables, rather than by 
using "extern" variables.
   The subroutine requires as input a set of weighted a.a. counts, without any 
added pseudocounts.  It returns an index with the units of bits.  This index 
should be positive if the observed counts are typical of proteins, and negative 
if the counts are atypical, with greater absolute values indicating greater 
significance.
   I haven't tested the code, so I'm not sure whether it contains any bugs, 
or whether it produces sensible results.  I'll be interested to learn what 
you find out.
 *****************************************************************************/
/*	Typicality index code							*/
/*      Version 1.0								*/
/*	August 24, 2015								*/
/*	Program by Stephen Altschul						*/

// double	lgamma(),log(),exp();

double  dms_typ::Typical(UInt4 *kount)
#if 0
         double	*c = letter counts.
	double	weight[];		/*  Weights of Dirichlet components	*/
	double	Alpha[];		/*  Concentration parameters		*/
	double	alpha[][20];		/*  Dirichlet parameters		*/
	double	back[];			/*  Background a.a. frequencies		*/
	int	M;			/*  Number of Dirichlet components	*/
#endif
{
 	double	lq[30];			/*  log adjusted observed a.a. freqs.	*/
 	// double	count;			/*  Number of observations		*/
 	double	logprob;		/*  Log-prob. density of freq. vector	*/
 	double	typ;			/*  Typicality index, in bits		*/
 	int	i,m;
	double  Count;

#if 1	// add in background counts for each deleted '-' residue.
	// fprintf(stderr,"Counts[X] = %d\n",kount[0]);
	for(Count=0.0,i=1; i <= nAlpha(AB); ++i){ 
		count[i]=back[i]*((double) kount[0]/(double) wt_factor);
		count[i] += (double) kount[i]/(double) wt_factor; Count+=count[i]; 
	} 
	// if(Count<=1.0) return(0.0);
#else	// use this to debug hieraln...
	for(Count=0.0,i=1; i <= nAlpha(AB); ++i){ 
		count[i]=(double) kount[i]/(double) wt_factor; Count+=count[i]; 
	} 
	// if(Count<=1.0) return(0.0);
#endif
 	// for (count=0,i=1; i <= nAlpha(AB); i++) count+=c[i]; original code...
/*  Calculate log frequencies for column, using 10 pseudocounts			*/
 	// for (i=1; i <= nAlpha(AB); i++) lq[i]=log((count[i]+10*back[i])/(Count+10));
 	for (i=1; i <= nAlpha(AB); i++) lq[i]=log((count[i]+back[i])/(Count+1));

/*  Calculate typicality index                                                  */

 	for (m=1; m <= nDC; m++) {
 		logprob=log(weight[m])+lgamma(Alpha[m]);
 		for (i=1; i <= nAlpha(AB); i++) logprob+=(alpha[m][i]-1)*lq[i]-lgamma(alpha[m][i]);
 		if (m==1) typ=logprob;
 		else if (logprob < typ) typ+=log(1.0+exp(logprob-typ));
 		else typ=logprob+log(1.0+exp(typ-logprob));
 	} typ -= lgamma(20.0);
 	// was typ -= lgamma(10.0);
 	// for (i=1; i <= nAlpha(AB); i++) typ-= (10*back[i]-1)*lq[i]-lgamma(10*back[i]);
 	typ/=log(2.0);
 	return(typ);
}



