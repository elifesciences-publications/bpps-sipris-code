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

#include "firm.h"

double	erfcc(double x)
{
        double	t,z,ans;

        z=fabs(x);
        t=1.0/(1.0+0.5*z);
        ans=t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+
                t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+
                t*(-0.82215223+t*0.17087277)))))))));
        return x >= 0.0 ? ans : 2.0-ans;
}
/* (C) Copr. 1986-92 Numerical Recipes Software =$j7`5|,$. */

double	lnerfcc(double x)
/** return the natural log of the complementary error function **/
{
        double	t,ans;

	if(x < 0.0) return log(erfcc(x));
        t=1.0/(1.0+0.5*x);
        ans=log(t)
		+(-x*x-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+
                t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+
                t*(-0.82215223+t*0.17087277)))))))));
        return ans;
}

double	FirmRandomWalk(Int4 iterations, char *name, void *X, double target, 
	double (*TotalScore)(void *), double (*Mutate)(double,double,void *,BooLean*))
{
	double	*barriers,p;

	NEW(barriers,1003,double);
	p=RandWalkBarriersFIRM(iterations,barriers,X,target,TotalScore,Mutate);
printf("rough p-value = %g\n",p);
	p=FirmRandomWalkPval(iterations,name,barriers,X,target,TotalScore,Mutate);
	free(barriers);
	return p;
}

double	FirmRandomWalkPval(Int4 iterations, char *name, double *barriers, 
	void *X, double target, double (*TotalScore)(void *),
	double (*Mutate)(double,double,void *,BooLean*))
/**********************************************************************
   iterations = number of iterations for random walk
   name = output file name;
   X = data structure to be mutated
   TotalScore( ) = function to compute the score
   Mutate( ) = function to mutate data structure
***********************************************************************/
{
	Int4	i,k,nbins,*dist,cycle,maxbins,recalc,bin,leftcnt,half;
	Int4	bar,maxcycle = INT4_MAX;
	double	score,log_tail,ave,stdev,d,*score_array;
	double	barrier,oldbarrier,fraction,binsize,max,min,norm;
	double	*hist,a,b,error,var,target_p,norm_p,final_fract;
	BooLean	accept,flag=TRUE,pflag=TRUE,debug=FALSE,pval=TRUE;
	Int4	naccept,stop=20;
	Int4	tbin,lookcnt;	/***** 3/4 above barrier *****/
	double	three4ths;	/***** 3/4 above barrier *****/
	
	if(pval==2) flag = FALSE;	/** pval=2: determine score **/
	if(pval==4) { stop = 0; pval = TRUE; } /** pval=4; stop immediately **/
	target_p = -DBL_MAX;
	bar = 0;

	recalc=10000; maxbins = 100; 
	MEW(dist,iterations+2,Int4); 
	MEW(score_array,iterations+2,double); 
	NEW(hist,maxbins+2,double);
	printf("target score = %g\n",target);
	/*************************************************************
         1a.Move score into the expected range for a random walk without 
		a barrier.
	 *************************************************************/
	for(barrier=-DBL_MAX,i=0;i<iterations;i++) {
	   if(i%recalc==0) score = TotalScore(X);
	   if(debug && i%50==0)
		   fprintf(stderr,"%d: score = %g\n", i,score);
	   score = Mutate(score, barrier,X, &accept); 
	}
	if(debug) printf("%d: score = %g\n",i,score); 
	/*************************************************************
         1b. Determine mean and standard deviation for Normal approximation.
	 *************************************************************/
        score = TotalScore(X); ave=stdev=0;
        barrier = -DBL_MAX;
	k = iterations*20;
	for(i=1;i<=k;i++){
	   if(i%recalc==0){
		score = TotalScore(X);
	   	if(pflag && i) fprintf(stderr,"count=%d\n", i);
	   }
	   score = Mutate(score, barrier,X, &accept); 
	   ave+=score; stdev+=score*score;
	}
	ave=ave/(double)k;
	stdev=sqrt((stdev/(double)k)-ave*ave);
	fprintf(stderr,"ave %g stdev %g\n",ave,stdev);
        score = TotalScore(X); 
     /************************** MAIN LOOP ************************
     2.Begin a process of random walks behind barriers. At each stage 
        the old or previous barrier is "oldbarrier" and "barrier" 
	is the new barrier we wish to use.
     *************************************************************/
     oldbarrier = -DBL_MAX;
     barrier = oldbarrier/1.0001;
     for(log_tail=0,cycle=1; cycle < maxcycle && oldbarrier < barrier; cycle++){
	/*************************************************************
	 3.Iterate using the old barrier till the new barrier is exceeded.
		When the new barrier is exceeded than we will use it
		instead of the old and start our random walk again.
	 *************************************************************/
	for(i=0; score < barrier && i<iterations;i++){
	   score = Mutate(score, oldbarrier,X, &accept); 
	   if(!accept) i--;
	   if(i%recalc==0) {
		score = TotalScore(X);
	   	if(pflag && i) fprintf(stderr,"count=%d\n",i);
	   }
	}
	/*************************************************************
	 4.If we are unable to exceed the new barrier we stop the whole 
		process at this point.
	 *************************************************************/
	if(score<barrier){ 
	   printf("score can't be reached\n");
	   printf("%d: score = %g; barrier = %g; oldbarrier = %g\n",
			i, score, barrier, oldbarrier);
	   break;
	}
	/*************************************************************
	 5.Having gotten above the new barrier we iterate a random walk 
		behind the new barrier.
	 *************************************************************/
	max=barrier; 
	min=DBL_MAX;
	for(naccept=i=0;i<iterations;i++){
	   if(i%recalc==0){
		score = TotalScore(X);
	   	if(pflag && i) fprintf(stderr,"count=%d\n", i);
	   }
	   score = Mutate(score, barrier,X, &accept); 
	   score_array[i]=score;
	   if(score<min)min=score;
	   if(score>max)max=score;
	   if(accept){ naccept++; }
	   if(score < barrier){
		fprintf(stderr,"%d) score = %g; barrier = %g; accept = %d\n", 
			i,score,barrier,accept);
		print_error("execution error in randwalk.c");
	   }
	   dist[i]=0;
	}
	dist[i]=0; /** bin CAN = iterations if score == max **/
	if(min == max) break;	/** if can't get beyond barrier quit **/
	/*************************************************************
	 6.We roughly find the midpoint of the new distribution and use 
		this to set the new barrier value for the next iteration 
		if conditions are right.
	 *************************************************************/
	binsize=(max-min)/(double)iterations;
	for(i=0;i<iterations;i++){
	   bin=(Int4)floor((score_array[i]-min)/binsize); dist[bin]++;
        }
	half=(iterations-1)/2;
	for(leftcnt=bin=0; leftcnt<half; ) leftcnt+=dist[bin++];

	if(flag) { /******* Find 3/4 score for target pvalue *********/
	   three4ths = (3*(iterations-1))/4;
	   for(lookcnt=tbin=0; lookcnt<three4ths; ) lookcnt+=dist[tbin++];
	   three4ths =tbin*binsize+min;
	}
	oldbarrier = barrier;
	bar++; 
	if(leftcnt >= (9*iterations)/10){ 
	   fraction=1.0; 
	   d = (double)(maxbins*iterations)/(double)iterations;
	   nbins=(Int4)floor(d+.5);
	} else { 
#if 0      /** OLD **/
	   max=barrier=bin*binsize+min;
#endif
#if 1	   /** NEW **/
	   max=barrier=barriers[bar];
#endif
	   fraction=((double)(iterations-leftcnt))/(double)iterations;
	   d = (double)(maxbins*leftcnt)/(double)iterations;
	   nbins=(Int4)floor(d+.5);
	   /********* Find and return score for target pvalue **********/
	   if(pval==2 && (log_tail + log(fraction)) < target){
	      for(score=min, leftcnt=bin=0; TRUE; score+=binsize){
	      	 final_fract=((double)(iterations-leftcnt))/(double)iterations;
                 target_p = log_tail + log(final_fract);
		 if(target_p <= target){
     			free(score_array); free(dist); free(hist);
			return score;
		 }
		 leftcnt+=dist[bin++]; 
	      }
	   /***********************************************************/
	   } else if(flag && three4ths > target){
	      /*********************************************************
	      7. Compute target statistics.
	        A. Calculate the target tail probability.
	      *********************************************************/
	      for(score=min, leftcnt=bin=0; score<target; score+=binsize) {
			leftcnt+=dist[bin++];
	      } score = TotalScore(X);
	      final_fract=((double)(iterations-leftcnt))/(double)iterations;
	      if(final_fract == 0){
		print_error("fatal error: target score exceeds limits.");
	      }
              target_p = log_tail + log(final_fract);
	      norm_p = lnerfcc((target-ave)/(stdev*sqrt(2)))-log(2);
	      /*********************************************************
	        B. Calculate the error.
	      *********************************************************/
              /** beta dist. for error **/
              a = b = (double) iterations/2.0;
              var = a*b/((a+b)*(a+b)*(a+b+1));
              error = 2*sqrt((double)(cycle) * var);
#if 0
              /** binomial dist. for error ****/
              error = sqrt((double)(cycle)/(double) iterations);
              printf("error_binomial = %g percent\n",error*100);
#endif
	      if(pval){
		/*********************************************************
		C. If we are computing the p-value then output the results 
			and stop at this point.
		*********************************************************/
		printf("error_beta = %g percent\n",error*100);
		printf("target(obs): S = %g; S0 = %g; p-value = 1.0e%.1f\n",
		target,(target - ave),target_p*log10(exp(1.0)));
		printf("target score is %g sd above the mean\n",
			(target-ave)/stdev);
		printf("normal p-value = %g\n", 
			(norm=erfcc((target-ave)/(stdev*sqrt(2)))/2));
		printf("computed p-value %.2f\n", target_p*log10(exp(1.0)));
		printf("%d iterations per barrier\n", iterations);
		if(norm > 0.0){
		   printf("#%s|I%.1lf|t%.1lf|a%.1lf|d%.1lf|e%.1lf|N%.1lf|\n",
		   name,-target_p,target,ave,stdev,error*100,-log(norm));
		} else {
		   printf("#%s|I%.1lf|t%.1lf|a%.1lf|d%.1lf|e%.1lf|N-inf|\n",
			name,-target_p,target,ave,stdev,error*100);
		}
		break;
	      } else maxcycle = cycle + stop;
	      flag = FALSE;
	   }
	}
	fprintf(stderr,"bar %g fac %g sh %g accept %d/%d\n",
		oldbarrier,log_tail,fraction,naccept,iterations);
	/*************************************************************
	 8.log_tail is updated ready for the next random walk. This value 
              gives the size discount for the histogram values that 
              will be generated at each stage.
	 *************************************************************/
	log_tail +=log(fraction);
     }
     free(score_array); free(dist); free(hist);
     return target_p;
}

double	RandWalkBarriersFIRM(Int4 iterations, double *barriers, 
	void *X, double target, double (*TotalScore)(void *),
	double (*Mutate)(double,double,void *,BooLean*))
/**********************************************************************
   iterations = number of iterations for random walk
   name = output file name;
   X = data structure to be mutated
   TotalScore( ) = function to compute the score
   Mutate( ) = function to mutate data structure
***********************************************************************/
{
	Int4	i,k,nbins,*dist,cycle,maxbins,recalc,bin,leftcnt,half;
	Int4	maxbar,bar,maxcycle = INT4_MAX;
	double	score,log_tail,ave,stdev,d,*score_array;
	double	barrier,oldbarrier,fraction,binsize,max,min;
	double	*hist,a,b,error,var,target_p,norm_p,final_fract;
	BooLean	accept,flag=TRUE,pflag=TRUE,debug=FALSE;
	Int4	naccept;
	Int4	tbin,lookcnt;	/***** 3/4 above barrier *****/
	double	three4ths;	/***** 3/4 above barrier *****/
	
	target_p = -DBL_MAX; 
	bar = 0;		/** initialize to 1st barrier **/
	maxbar = 1000;

	recalc=10000; maxbins = 100;
	MEW(dist,iterations+2,Int4);
	MEW(score_array,iterations+2,double);
	NEW(hist,maxbins+2,double);
	printf("target score = %g\n",target);
	/*************************************************************
         1a.Move score into the expected range for a random walk without 
		a barrier.
	 *************************************************************/
	for(barrier=-DBL_MAX,i=0;i<iterations;i++) {
	   if(i%recalc==0) score = TotalScore(X);
	   if(debug && i%50==0)
		   fprintf(stderr,"%d: score = %g\n", i,score);
	   score = Mutate(score, barrier,X, &accept);
	}
	if(debug) printf("%d: score = %g\n",i,score);
	/*************************************************************
         1b. Determine mean and standard deviation for Normal approximation.
	 *************************************************************/
        score = TotalScore(X); ave=stdev=0;
        barrier = -DBL_MAX;
	k = iterations*20;
	for(i=1;i<=k;i++){
	   if(i%recalc==0){
		score = TotalScore(X);
	   	if(pflag && i) fprintf(stderr,"count=%d\n", i);
	   }
	   score = Mutate(score, barrier,X, &accept); 
	   ave+=score; stdev+=score*score;
	}
	ave=ave/(double)k;
	stdev=sqrt((stdev/(double)k)-ave*ave);
	fprintf(stderr,"ave %g stdev %g\n",ave,stdev);
        score = TotalScore(X); 
     /************************** MAIN LOOP ************************
     2.Begin a process of random walks behind barriers. At each stage 
        the old or previous barrier is "oldbarrier" and "barrier" 
	is the new barrier we wish to use.
     *************************************************************/
     oldbarrier = -DBL_MAX;
     barrier = oldbarrier/1.0001;
     for(log_tail=0,cycle=1; cycle < maxcycle && oldbarrier<barrier; cycle++){
	/*************************************************************
	 3.Iterate using the old barrier till the new barrier is exceeded.
		When the new barrier is exceeded than we will use it
		instead of the old and start our random walk again.
	 *************************************************************/
	for(i=0; score < barrier && i<iterations;i++){
	   score = Mutate(score, oldbarrier,X, &accept); 
	   if(!accept) i--;
	   if(i%recalc==0) {
		score = TotalScore(X);
	   	if(pflag && i) fprintf(stderr,"count=%d\n",i);
	   }
	}
	/*************************************************************
	 4.If we are unable to exceed the new barrier we stop the whole 
		process at this point.
	 *************************************************************/
	if(score < barrier){ 
	   printf("score can't be reached\n");
	   printf("%d: score = %g; barrier = %g; oldbarrier = %g\n",
			i, score, barrier, oldbarrier);
	   break;
	}
	/*************************************************************
	 5.Having gotten above the new barrier we iterate a random walk 
		behind the new barrier.
	 *************************************************************/
	max=barrier; 
	min=DBL_MAX;
	for(naccept=i=0;i<iterations;i++){
	   if(i%recalc==0){
		score = TotalScore(X);
	   	if(pflag && i) fprintf(stderr,"count=%d\n", i);
	   }
	   score = Mutate(score, barrier,X, &accept); 
	   score_array[i]=score;
	   if(score<min)min=score;
	   if(score>max)max=score;
	   if(accept){ naccept++; }
	   if(score < barrier){
		fprintf(stderr,"%d) score = %g; barrier = %g; accept = %d\n", 
			i,score,barrier,accept);
		print_error("execution error in randwalk.c");
	   }
	   dist[i]=0;
	}
	dist[i]=0; /** bin CAN = iterations if score == max **/
	if(min == max) break;	/** if can't get beyond barrier quit **/
	/*************************************************************
	 6.We roughly find the midpoint of the new distribution and use 
		this to set the new barrier value for the next iteration 
		if conditions are right.
	 *************************************************************/
	binsize=(max-min)/(double)iterations;
	for(i=0;i<iterations;i++){
	   bin=(Int4)floor((score_array[i]-min)/binsize); dist[bin]++;
        }
	half=(iterations-1)/2;
	for(leftcnt=bin=0; leftcnt<half; ) leftcnt+=dist[bin++];

	if(flag) { /******* Find 3/4 score for target pvalue *********/
	   three4ths = (3*(iterations-1))/4;
	   for(lookcnt=tbin=0; lookcnt<three4ths; ) lookcnt+=dist[tbin++];
	   three4ths =tbin*binsize+min;
	}
	oldbarrier = barrier;
	if(leftcnt >= (9*iterations)/10){ 
	   fraction=1.0; 
	   d = (double)(maxbins*iterations)/(double)iterations;
	   nbins=(Int4)floor(d+.5);
	} else { 
	   max=barrier=bin*binsize+min;
	   fraction=((double)(iterations-leftcnt))/(double)iterations;
	   d = (double)(maxbins*leftcnt)/(double)iterations;
	   nbins=(Int4)floor(d+.5);
	   if(flag && three4ths > target){
	      /*********************************************************
	      7. Compute target statistics.
	        A. Calculate the target tail probability.
	      *********************************************************/
	      for(score=min, leftcnt=bin=0; score<target; score+=binsize) {
			leftcnt+=dist[bin++];
	      } score = TotalScore(X);
	      final_fract=((double)(iterations-leftcnt))/(double)iterations;
	      if(final_fract == 0){
		print_error("fatal error: target score exceeds limits.");
	      }
              target_p = log_tail + log(final_fract);
	      norm_p = lnerfcc((target-ave)/(stdev*sqrt(2)))-log(2);
	      /*********************************************************
	        B. Calculate the error.
	      *********************************************************/
              /** beta dist. for error **/
              a = b = (double) iterations/2.0;
              var = a*b/((a+b)*(a+b)*(a+b+1));
              error = 2*sqrt((double)(cycle) * var);
	      maxcycle = cycle + 10;
	      flag = FALSE;
	   }
	}
	fprintf(stderr,"bar %g fac %g sh %g accept %d/%d\n",
		oldbarrier,log_tail,fraction,naccept,iterations);
	if(++bar >= 1000) print_error("reset MAX_NO_BARRIER");
	barriers[bar] = barrier; 
	/*************************************************************
	 9.log_tail is updated ready for the next random walk. This value 
              gives the size discount for the histogram values that 
              will be generated at each stage.
	 *************************************************************/
	log_tail +=log(fraction);
     }
     free(score_array); free(dist); free(hist);
     return target_p;
}

