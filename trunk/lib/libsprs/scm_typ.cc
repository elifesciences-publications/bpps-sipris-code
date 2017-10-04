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

#include "scl_typ.h"

#if 0
/*	Subroutine to calculate P-values for			*/
/*	initial match clusters in sequences			*/
/*								*/
/*	Program by Stephen F. Altschul				*/
/*	Version 1.11.3;  January 26, 2017			*/
/* 	int	jeff;		/*   Flat: 0;  Jeffreys': 1	*/
/* 	int	cflag;          /*   p-value correction		*/
/* 	int	pflag;          /*   Logistic p-values		*/
#endif

Int4    scm_typ::pvcalcF(Int4 L,Int4 D,Int4 *pos,double &pval,BooLean jeff,BooLean cflag, BooLean pflag)
{
 	Int4	M,m,n,x,y,minm,maxm,bestm;
 	double	sum,off,score,bestsc,pv;
 	static	double	*lni,*lrat,*offset;
 	static	Int4	oldL=0;
 	static	Int4	oldD=0;

/*	Check that input makes sense				*/
	// fprintf(stderr,"L=%d; D=%d.\n",L,D);
 	if (L<2 || D<1 || D>=L) return(-1);
 	pos[0]=0;
 	for (m=1;m<=D;++m) if (pos[m]<=pos[m-1] || pos[m]>L) return(-1);

/*	Initialize arrays					*/

 	if (L>oldL) {
 		if (oldL) { free(lni); free(lrat); free(offset); }
 		lni= (double *) calloc(L+1,sizeof(double));
 		lrat= (double *) calloc(L,sizeof(double));
 		offset= (double *) calloc(L,sizeof(double));
 		lni[0]=lrat[0]=sum=0;
 		for (x=1;x<L;++x) { lni[x]=x*log(x); sum+=log(x); lrat[x]=lni[x]-sum; }
 		lni[L]=L*log(L);
 	}
 	if (L!=oldL || D!= oldD) for (x=1;x<L;++x) offset[x]=0;
 	oldL=L; oldD=D;

/*	Find best front-weighted cut				*/

 	bestm=0;
 	for (M=1;M<=D;++M) if (M*L > D*(x=pos[M])) {
 		y=L-x;
 		if (offset[x]==0) {
 			sum=0;
 			off=lrat[x]+lrat[y];
 			minm= (y>=D) ? 0 : D-y;
 			maxm= (x<=D) ? x : D;
 			for (m=minm;m<=maxm;++m) {
 				n=D-m;
 				sum+=exp(lrat[m]+lrat[x-m]+lrat[n]+lrat[y-n]-off);
 			}
 			offset[x]=lni[x]+lni[y]+log(sum);
 			if (jeff==0) offset[x]+=log(1.0/(x*x)+1.0/(y*y))/2;
 			offset[y]=offset[x];
 		}
 		n=D-M;
 		score=lni[M]+lni[x-M]+lni[n]+lni[y-n]-offset[x];
 		if (bestm==0 || score>bestsc) { bestm=M; bestsc=score; }
 	}

/*	Return best number of matches and p-value		*/
 	pval=1.0;
 	if (bestm) {
 		bestsc-=lni[D]+lni[L-D]-lni[L];
 		pv = exp(-bestsc)*sqrt(D/3.14159)/2;
 		pv*= (jeff) ? log(1.0242*L) : L-1;
 		if (cflag) pv*=(L-D+1-2.0/(D+1))/(L-1);
 		if (pv<1) pval = pflag ? pv/(1+pv) : pv;
 		else { bestm=0; if (pflag) pval = pv/(1+pv); }

 	} return(bestm);
}

Int4    scm_typ::pvcalcD(Int4 L,Int4 D,Int4 *pos,long double &pval,BooLean jeff,BooLean cflag, BooLean pflag)
// uses long double...
{
 	Int4	M,m,n,x,y,minm,maxm,bestm;
 	long double	sum,off,score,bestsc,pv;
 	static	long double	*lni,*lrat,*offset;
 	static	Int4	oldL=0,oldD=0;

/*	Check that input makes sense				*/
 	if (L<2 || D<1 || D>=L) return(-1);
 	pos[0]=0;
 	for (m=1;m<=D;++m) if (pos[m]<=pos[m-1] || pos[m]>L) return(-1);

/*	Initialize arrays					*/

 	if (L>oldL) {
 		if (oldL) { free(lni); free(lrat); free(offset); }
 		lni= (long double *) calloc(L+1,sizeof(long double));
 		lrat= (long double *) calloc(L,sizeof(long double));
 		offset= (long double *) calloc(L,sizeof(long double));
 		lni[0]=lrat[0]=sum=0;
 		for (x=1;x<L;++x) { lni[x]=x*logl(x); sum+=logl(x); lrat[x]=lni[x]-sum; }
 		lni[L]=L*logl(L);
 	}
 	if (L!=oldL || D!= oldD) for (x=1;x<L;++x) offset[x]=0;
 	oldL=L; oldD=D;

/*	Find best front-weighted cut				*/

 	bestm=0;
 	for (M=1;M<=D;++M) if (M*L > D*(x=pos[M])) {
 		y=L-x;
 		if (offset[x]==0) {
 			sum=0; off=lrat[x]+lrat[y];
 			minm= (y>=D) ? 0 : D-y;
 			maxm= (x<=D) ? x : D;
 			for (m=minm;m<=maxm;++m) {
 				n=D-m; sum+=expl(lrat[m]+lrat[x-m]+lrat[n]+lrat[y-n]-off);
 			}
 			offset[x]=lni[x]+lni[y]+logl(sum);
 			if (jeff==0) offset[x]+=logl(1.0/(x*x)+1.0/(y*y))/2;
 			offset[y]=offset[x];
 		} n=D-M; score=lni[M]+lni[x-M]+lni[n]+lni[y-n]-offset[x];
 		if (bestm==0 || score>bestsc) { bestm=M; bestsc=score; }
 	}

/*	Return best number of matches and p-value		*/
 	pval=1.0;
 	if (bestm) {
 		bestsc-=lni[D]+lni[L-D]-lni[L];
 		pv = expl(-bestsc)*sqrtl(D/3.14159)/2;
 		pv*= (jeff) ? logl(1.0242*L) : L-1;
 		if (cflag) pv*=(L-D+1-2.0/(D+1))/(L-1);
 		if (pv<1) pval = pflag ? pv/(1+pv) : pv;
 		else { bestm=0; if (pflag) pval = pv/(1+pv); }

 	} return(bestm);
}

Int4    scm_typ::pvcalc1(Int4 L,Int4 D,Int4 *pos,double &pval,BooLean jeff,BooLean cflag)
/*      Subroutine to calculate P-values for                    */
/*      initial match clusters in sequences                     */
/*                                                              */
/*      Program by Stephen F. Altschul                          */
/*      Version 1.11.2;  Decemeber 16, 2016                     */
{
        Int4     M,m,n,x,y,minm,maxm,bestm;
        double  sum,off,score,bestsc,pv;
        static  double  *lni,*lrat,*offset;
        static  Int4	oldL=0;
        static  Int4    oldD=0;

/*      Check that input makes sense                            */

        if (L<2 || D<1 || D>=L) return(-1);
        pos[0]=0;
        for (m=1;m<=D;++m) if (pos[m]<=pos[m-1] || pos[m]>L) return(-1);

/*      Initialize arrays                                       */

        if (L>oldL) {
                if (oldL) {
                        free(lni);
                        free(lrat);
                        free(offset);
                }
                lni= (double *) calloc(L+1,sizeof(double));
                lrat= (double *) calloc(L,sizeof(double));
                offset= (double *) calloc(L,sizeof(double));
                lni[0]=lrat[0]=sum=0;
                for (x=1;x<L;++x) {
                        lni[x]=x*log(x);
                        sum+=log(x);
                        lrat[x]=lni[x]-sum;
                }
                lni[L]=L*log(L);
        }
        if (L!=oldL || D!= oldD) for (x=1;x<L;++x) offset[x]=0;
        oldL=L;
        oldD=D;

/*      Find best front-weighted cut                            */

        bestm=0;
        for (M=1;M<=D;++M) if (M*L > D*(x=pos[M])) {
                y=L-x;
                if (offset[x]==0) {
                        sum=0;
                        off=lrat[x]+lrat[y];
                        minm= (y>=D) ? 0 : D-y;
                        maxm= (x<=D) ? x : D;
                        for (m=minm;m<=maxm;++m) {
                                n=D-m;
                                sum+=exp(lrat[m]+lrat[x-m]+lrat[n]+lrat[y-n]-off);
                        }
                        offset[x]=lni[x]+lni[y]+log(sum);
                        if (jeff==0) offset[x]+=log(1.0/(x*x)+1.0/(y*y))/2;
                        offset[y]=offset[x];
                }
                n=D-M;
                score=lni[M]+lni[x-M]+lni[n]+lni[y-n]-offset[x];
                if (bestm==0 || score>bestsc) {
                        bestm=M;
                        bestsc=score;
                }
        }

/*      Return best number of matches and p-value               */

        pval=1.0;
        if (bestm) {
                bestsc-=lni[D]+lni[L-D]-lni[L];
                pv= exp(-bestsc)*sqrt(D/3.14159)/2;
                pv*= (jeff) ? log(1.024*L) : L-1;
                if (cflag) pv*=(L-D+1-2.0/(D+1))/(L-1);
                if (pv < pval) pval=pv;
                else bestm=0;
        }
        return(bestm);
}

Int4    scm_typ::pvcalc3(Int4 L,Int4 D,Int4 *pos,double &pval,BooLean jeff,BooLean cflag)
/*****************************************************************/
/*      Subroutine to calculate P-values for                    */
/*      initial match clusters in sequences                     */
/*                                                              */
/*      Program by Stephen F. Altschul                          */
/*      Version 1.11.1;  Decemeber 16, 2016                     */
/*****************************************************************/
{
        Int4     M,m,n,x,y,minm,maxm,bestm;
        double  sum,off,score,bestsc,pv;
        static  double  *lni,*lrat,*offset;
        static  Int4	oldL=0;
        static  Int4    oldD=0;

/*      Check that input makes sense                            */

        if (L<2 || D<1 || D>=L) return(-1);
        pos[0]=0;
        for (m=1;m<=D;++m) if (pos[m]<=pos[m-1] || pos[m]>L) return(-1);

/*      Initialize arrays                                       */

        if (L>oldL) {
                if (oldL) {
                        free(lni);
                        free(lrat);
                        free(offset);
                }
                lni= (double *) calloc(L+1,sizeof(double));
                lrat= (double *) calloc(L,sizeof(double));
                offset= (double *) calloc(L,sizeof(double));
                lni[0]=lrat[0]=sum=0;
                for (x=1;x<L;++x) {
                        lni[x]=x*log(x);
                        sum+=log(x);
                        lrat[x]=lni[x]-sum;
                }
                lni[L]=L*log(L);
        }
        if (L!=oldL || D!= oldD) for (x=1;x<L;++x) offset[x]=0;
        oldL=L;
        oldD=D;

/*      Find best front-weighted cut                            */
        bestm=0;
        for (M=1;M<=D;++M) if (M*L > D*(x=pos[M])) {
                y=L-x;
                if (offset[x]==0) {
                        sum=0;
                        off=lrat[x]+lrat[y];
                        minm= (y>=D) ? 0 : D-y;
                        maxm= (x<=D) ? x : D;
                        for (m=minm;m<=maxm;++m) {
                                n=D-m;
                                sum+=exp(lrat[m]+lrat[x-m]+lrat[n]+lrat[y-n]-off);
                        }
                        offset[x]=lni[x]+lni[y]+log(sum);
                        if (jeff==0) offset[x]+=log(1.0/(x*x)+1.0/(y*y))/2;
                        offset[y]=offset[x];
                }
                n=D-M;
                score=lni[M]+lni[x-M]+lni[n]+lni[y-n]-offset[x];
                if (bestm==0 || score>bestsc) {
                        bestm=M;
                        bestsc=score;
                }
        }

/*      Return best number of matches and p-value               */
        pval=1.0;
        if (bestm) {
                bestsc-=lni[D]+lni[L-D]-lni[L];
                pv= exp(-bestsc)*sqrt(D/3.14159)/2;
                pv*= (jeff) ? log(1.024*L) : L-1;
                if (cflag && !jeff) pv*=((L-D+1)*(D+1)-2.0)/((L-1)*(D+1));
                if (pv < pval) pval=pv;
                else bestm=0;
        }
        return(bestm);
}

Int4    scm_typ::pvcalc0(Int4 L,Int4 D,Int4 *pos,double &pval,BooLean jeff)
// new routine that compensates for edge effects...
/*	Subroutine to calculate P-values for			*/
/*	initial match clusters in sequences			*/
/*								*/
/*	Program by Stephen F. Altschul				*/
/*	Version 1.10.11;  December 12, 2016			*/
{
 	Int4	M,m,n,x,y,minm,maxm,bestm;
 	double	sum,summ,sumt,temp,off,score,bestsc,pv;
 	static	double	*lni,*lrat,*offset;
 	static	Int4	oldL=0;
 	static	Int4	oldD=0;

/*	Check that input makes sense				*/
 	if (L<2 || D<1 || D>=L) return(-1);
 	pos[0]=0;
 	for (m=1;m<=D;++m) if (pos[m]<=pos[m-1] || pos[m]>L) return(-1);

/*	Initialize arrays					*/
 	if (L>oldL) {
 		if (oldL) { free(lni); free(lrat); free(offset); }
 		lni= (double *) calloc(L+1,sizeof(double));
 		lrat= (double *) calloc(L,sizeof(double));
 		offset= (double *) calloc(L,sizeof(double));
 		lni[0]=lrat[0]=sum=0;
 		for (x=1;x<L;++x) {
 			lni[x]=x*log(x);
 			sum+=log(x);
 			lrat[x]=lni[x]-sum;
 		} lni[L]=L*log(L);
 	}
 	if (L!=oldL || D!= oldD) for (x=1;x<L;++x) offset[x]=0;
 	oldL=L; oldD=D;

/*	Find best front-weighted cut				*/
 	bestm=0;
 	for (M=1;M<=D;++M) if (M*L > D*(x=pos[M])) {
 		y=L-x;
 		if (offset[x]==0) {
 			sumt=summ=sum=0;
 			off=lrat[x]+lrat[y];
 			minm= (y>=D) ? 0 : D-y;
 			maxm= (x<=D) ? x : D;
 			for (m=minm;m<=maxm;++m) {
 				n=D-m;
 				sum+=exp(lrat[m]+lrat[x-m]+lrat[n]+lrat[y-n]-off);
 				sumt+=temp=exp(lrat[m]+lrat[n]);
 				summ+=m*temp;
 			}
 			offset[x]=lni[x]+lni[y]+log(sum);
 			summ/=sumt;
 			if (jeff==0) offset[x]+=log(summ/(x*x)+(D-summ)/(y*y))/2;
 			offset[y]=offset[x];
 		}
 		n=D-M;
 		score=lni[M]+lni[x-M]+lni[n]+lni[y-n]-offset[x];
 		if (bestm==0 || score>bestsc) {
 			bestm=M;
 			bestsc=score;
 		}
 	}

/*	Return best number of matches and p-value		*/

 	pval=1.0;
 	if (bestm) {
 		bestsc-=lni[D]+lni[L-D]-lni[L];
 		pv= exp(-bestsc)/sqrt(2*3.14159);
 		pv*= (jeff) ? sqrt(D/2.0)*log(1.024*L) : L-1;
 		if (pv < pval) pval=pv;
 		else bestm=0;
 	}
 	return(bestm);
}


Int4	scm_typ::pvcalc2(Int4 L,Int4 D,Int4 *pos,double &pval, BooLean jeff)
/*	Subroutine to calculate P-values for			*/
/*	initial match clusters in sequences			*/
/*								*/
/*	Program by Stephen F. Altschul				*/
/*	Version 1.10.8;  September 26, 2016			*/
/********************************************************************
 Concerning the modification of the problem you described to me last week, I think there is 
probably a simple solution.  If you restrict the possible cuts to positions <= C (which 
can be done by skipping the main loop whenever x>C) then the only modification in the 
calculation of P-values for the x>flat-priors case is to replace the factor (L-1) by 
the factor C.  
 For the Jeffreys' prior case, it is slightly more complicated.  You need to replace the 
factor log(1.024*L) by the 0.5 times the sum, for i=1 to C, of sqrt[1/i**2 + 1/(L-i)**2].
(Note that for C==L-1, the formula log(1.024*L) is a very close approximation for this 
sum when L is greater than some very small integer.  You can check this out.)

 ********************************************************************/
{
 	Int4	M,m,n,x,y,minm,maxm,bestm;
 	double	sum,off,score,bestsc,pv;
 	static	double	*lni,*lrat,*offset;
 	static	Int4	oldL=0,oldD=0;

/*	Check that input makes sense				*/
 	if (L<2 || D<1 || D>=L) return(-1);
 	pos[0]=0;
 	for (m=1;m<=D;++m) if (pos[m]<=pos[m-1] || pos[m]>L) return(-1);

/*	Initialize arrays					*/
 	if (L>oldL) {
 		if (oldL) {
 			free(lni);
 			free(lrat);
 			free(offset);
 		}
 		lni= (double *) calloc(L+1,sizeof(double));
 		lrat= (double *) calloc(L,sizeof(double));
 		offset= (double *) calloc(L,sizeof(double));
 		lni[0]=lrat[0]=sum=0;
 		for (x=1;x<L;++x) {
 			lni[x]=x*log(x);
 			sum+=log(x);
 			lrat[x]=lni[x]-sum;
 		}
 		lni[L]=L*log(L);
 	}
 	if (L!=oldL || D!= oldD) for (x=1;x<L;++x) offset[x]=0;
 	oldL=L; oldD=D;

/*	Find best front-weighted cut				*/

 	bestm=0;
 	for (M=1;M<=D;++M){
	   if (M*L > D*(x=pos[M])) {
 		y=L-x;
 		if (offset[x]==0) {
 			sum=0;
 			off=lrat[x]+lrat[y];
 			minm= (y>=D) ? 0 : D-y;
 			maxm= (x<=D) ? x : D;
 			for (m=minm;m<=maxm;++m) {
 				n=D-m;
 				sum+=exp(lrat[m]+lrat[x-m]+lrat[n]+lrat[y-n]-off);
 			}
 			offset[x]=lni[x]+lni[y]+log(sum);
 			if (!jeff) offset[x]+=log(1.0/(x*x)+1.0/(y*y))/2;
 			offset[y]=offset[x];
 		}
 		n=D-M;
 		score=lni[M]+lni[x-M]+lni[n]+lni[y-n]-offset[x];
 		if (bestm==0 || score>bestsc) { bestm=M; bestsc=score; }
 	   }
 	}

/*	Return best number of matches and p-value		*/
 	pval=1.0;
 	if (bestm) {
 		bestsc-=lni[D]+lni[L-D]-lni[L];
 		pv= exp(-bestsc)*sqrt(D/3.14159)/2;
 		pv*= (jeff) ? log(1.024*L) : L-1;
 		if (pv < pval) pval=pv;
 		else bestm=0;
 	}
 	return(bestm);
}

float	**scm_typ::GetHbondDistMtrx( )
//============= Move this routine to libpdb.a at some point: ================
{
	Int4	ai,aj,r,i,j;

	h_type HG=0;
	if(efptr) HG=Histogram("residue H-bond distances",0,100,1.0);
	// Construct minimum Res-Res distance matrix...
	atm_typ	atmI,atmJ;
	res_typ	ResI,ResJ;
	char	cI,cJ;
	float	max,min,dd;	// 350 x 350 = 122,500 x 16 = 1,470 kBytes.
	double	Dd;
	NEWP(DistMtrx, num_resC + 3, float);
	for(i=1; i <= num_resC; i++){
	   NEW(DistMtrx[i], num_resC + 3, float); 
	   for(j=1; j <= num_resC; j++) DistMtrx[i][j] = 999999.0;
	}
	for(i=1; i <= num_resC; i++){
	   ResI=ResALL[i]; 
	   Int4 si=ResidueID(ResI);
	   if(si < resStart || si > resEnd) continue;
	   cI=GetCharResidue(ResI,AB);
	   for(j=1; j <= num_resC; j++){
	     if(i == j) continue;
	     ResJ=ResALL[j]; 
	     Int4 sj=ResidueID(ResJ);
	     if(sj < resStart || sj > resEnd) continue;
	     cJ=GetCharResidue(ResJ,AB); dd=0.0;
	     Int4	h,ii,jj;
		// create arrays of all atoms...
	     atm_typ atmsI[200],atmsJ[200],*H;
	     //=========  Find acceptor atoms... ===========
	     for(ii=0,ai=1; ai <= ResidueAtomNumber(ResI); ai++){
		   atmI=AtomResidue(ai,ResI);
		   if(CarbonAtom(atmI)){
		      if(!(SideAtom(atmI) && (cI=='W' || cI=='Y' || cI == 'F'))) continue;
		   }
		   if(!SideAtom(atmI) && !OxygenAtom(atmI)) continue;
		   ii++; atmsI[ii]=atmI;
	     } ii++; atmsI[ii]=0;
	     //=========  Find donor atoms... =============
	     for(jj=0,aj=1; aj <= ResidueAtomNumber(ResJ); aj++){
		   // atmJ=AtomResidue(aj,ResJ);	// Donor atom
		   // if(CarbonAtom(atmJ)) continue;
		   H=ResidueAtomHydrogens(aj,ResJ);
		   if(H) for(h=1; H[h]; h++){ jj++; atmsJ[jj]=H[h]; }
	     } jj++; atmsJ[jj]=0;
	     //======== Find H-bonds based on D-H...A distances... ==========
	     for(min=999999.9,ai=1; atmsI[ai]; ai++){
		   for(atmI=atmsI[ai],aj=1; (atmJ=atmsJ[aj]) != 0; aj++){
		      // ignore backbone to backbone H-bonds.
		      if(BackboneHbonds){ 
		        if(IsBackboneAtom(atmI) && (IsBackboneAtom(atmJ) && !IsBackboneHydrogenAtom(atmJ))) continue;
		      } else if(IsBackboneAtom(atmI) && IsBackboneHydrogenAtom(atmJ)) continue;
		      dd=DistanceAtoms(atmI,atmJ);
		      if(dd < min) min=dd;
		   }
	     } 
#if 1	// find CH-pi, NH-pi, etc. H-bonds.
	     if(AroPiBonds){
	       double Dd=9999999.9;
	       float dmax=4.5;
	       Int4 hits=0;
	       if(cI=='W' || cI=='Y' || cI == 'F' || cI == 'H'){  // ResJ=Donor; ResI=Acceptor; 
	     	  hits=FindAromaticHbondsPDB(0,ResJ,ResI,C,C,dmax,0,P,Dd);
	       } else if(cJ=='W' || cJ=='Y' || cJ == 'F' || cJ == 'H'){ // ResI=Donor ; ResJ=Acceptor
	     	  hits=FindAromaticHbondsPDB(0,ResI,ResJ,C,C,dmax,0,P,Dd);
	       } Dd = 0.75*Dd;	// Adjust for the fact that CH-pi bonds are further appart.
	       if(hits > 0 && min > (float) Dd) min=(float)Dd;
	     }
	// find CH-pi, NH-pi, etc. H-bonds.
	     if(OtherPiBonds){
	       double Dd=9999999.9;
	       float dmax=4.5;
	       Int4 hits=0;
	       if(cI == 'N' || cI == 'Q' || cI == 'D' || cI == 'E' || cI == 'R'){
	          // ResJ=Donor; ResI=Acceptor; 
	     	  hits=FindOtherPiHbondsPDB(0,ResJ,ResI,C,C,dmax,0,P,Dd);
	       } else if(cJ == 'N' || cJ == 'Q' || cJ == 'D' || cJ == 'E' || cJ == 'R'){
		  // ResI=Donor ; ResJ=Acceptor
	     	  hits=FindOtherPiHbondsPDB(0,ResI,ResJ,C,C,dmax,0,P,Dd);
	       } Dd = 0.75*Dd;
	       if(hits > 0 && min > (float) Dd) min=(float)Dd;
	     }
#endif
#if 0	// DEBUG..
	     if( (cI=='Y' || cJ=='Y') && (ResidueID(ResI)==89 || ResidueID(ResJ)==89)){
		// fprintf(stderr,"%c%d ... %c%d: min=%3.f\n",cI,i,cJ,j,min);
		if(efptr) fprintf(stderr,"%c%d ... %c%d: min=%.3f\n",cI,ResidueID(ResI),cJ,ResidueID(ResJ),min);
	     }
#endif
	     if(min < DistMtrx[i][j]) DistMtrx[i][j]=DistMtrx[j][i]=min; 
	     if(efptr) IncdHist(min,HG);
	   }
	} if(efptr){ PutHist(stderr,60,HG); NilHist(HG); }
	return DistMtrx;
}

float	**scm_typ::GetMtrxDCA(char *dca_file, BooLean **&Used)
//============= Move this routine to libpdb.a at some point: ================
// WARNING: make sure that EVfold corresponds to the residue numbering in pdb files!
{
	Int4	r,I,J,i,j,ii,jj;

	// Construct minimum Res-Res distance matrix...
	h_type HG=0;
	if(efptr) HG=Histogram("residue distances",-1,1,0.05);
        res_typ ResI,ResJ;
	Int4	rI,rJ;
	float	max,min,dd,**DM;	// 350 x 350 = 122,500 x 16 = 1,470 kBytes.
	float	**DstMtx;
	double	Dd;
	Int4	max_res=MaxResPDB(C,P);
	assert(max_res >= resEnd);
	NEWP(Used,num_resC+3, BooLean);
	for(i=1; i <= num_resC; i++){ NEW(Used[i], num_resC + 3, BooLean); }
	NEWP(DstMtx, num_resC + 3, float);
	for(i=1; i <= num_resC; i++){ NEW(DstMtx[i], num_resC + 3, float); }
	NEWP(DM, max_res+ 3, float);
	for(i=1; i <= max_res; i++){ NEW(DM[i], max_res+ 3, float); }
	char	Str[105],cI,cJ,ci,cj;
	e_type	E=GetPDBSeq(C,P),xE=0;
	Int4	os=OffSetSeq(E);
	xE=MkEmptySeq(1,"DCA seq",LenSeq(E) + os + 5);
	FILE *fp=0;
#if 0	// create DCA Seq...move to scl_typ eventually.
	fp=open_file(dca_file,"","r");
	Int4	offset=0,max_ij=0,min_ij=9999;
	unsigned char *seq; 	NEW(seq,1005,unsigned char); // 30-500 is evfold max...
	while(fgets(Str,100,fp) != NULL){
		Int4 dm[12];	// dummy variables.
		if(sscanf(Str,"%d,%d,%lf,%d,%d,%d,%d,%d,%d,%d,%c,%c,",
		   &i,&j,&Dd,&dm[1],&dm[2],&dm[3],&dm[4],&dm[5],&dm[6],&dm[7],&cI,&cJ) != 12) 
			print_error("GetMtrxDCA() input error 1.");
		assert(i < 1000 && j < 1000 && i > 0 && j > 0);
		if(i > max_ij) max_ij=i; if(j > max_ij) max_ij=j; 
		if(i < min_ij) min_ij=i; if(j < min_ij) min_ij=j; 
		seq[i]=AlphaCode(cI,AB); seq[j]=AlphaCode(cJ,AB);
		// EqSeq(i,AlphaCode(cI,AB),dcaE); EqSeq(j,AlphaCode(cJ,AB),dcaE);
	} fclose(fp);
	e_type dcaE=MkSeq("DCA seq",max_ij-min_ij +1,seq + min_ij -1); free(seq);
	SetOffSetSeq(min_ij -1,dcaE);
	PutSeq(stderr,dcaE,AB);
#endif
#if 0
	AlnSeqSW(stderr,11, 1,dcaE, E, AB);
	if(!this->OverlappingSeqs(offset,30,0,dcaE,E)){
		// OverlappingSeq() sets offset > 0 if query starts before pdbseq.
                // It also allows for 'X' residues in pdb files...
                // else it sets offset < 0.
	   fprintf(stderr,"Sequences mismatch\n");
	   PutSeq(stderr,dcaE,AB);
	   PutSeq(stderr,E,AB);
	} else if(offset != 0) fprintf(stderr,"Offset = %d\n",offset);
#endif
	fp=open_file(dca_file,"","r");
	while(fgets(Str,100,fp) != NULL){
		Int4 dm[12];	// dummy variables.
		if(sscanf(Str,"%d,%d,%lf,%d,%d,%d,%d,%d,%d,%d,%c,%c,",
		   &i,&j,&Dd,&dm[1],&dm[2],&dm[3],&dm[4],&dm[5],&dm[6],&dm[7],&cI,&cJ) != 12) 
			print_error("GetMtrxDCA() input error 1.");
		i = i - dcaOS; j = j - dcaOS; 
		ii= i - OffSetSeq(E); ci=ResSeq(ii,E);
		jj= j - OffSetSeq(E); cj=ResSeq(jj,E);
		if(ci == 0 || cj == 0) continue;	// residue not visible in pdb file.
		ii=this->FindIndex(i); jj=this->FindIndex(j);
		if(i < 1 || j < 1 || ii < 1 || jj < 1){
		   fprintf(stderr,"i=%d; j=%d; ii=%d; jj=%d\n",i,j,ii,jj);
		   AlnSeqSW(stderr,11, 1,dcaE, E, AB);
		   PutSeq(stderr,E,AB); PutSeq(stderr,xE,AB); PutSeq(stderr,dcaE,AB);
		   print_error("GetMtrxDCA() input error 2.");
		}
	        ci=GetCharResidue(ResALL[ii],AB); cj=GetCharResidue(ResALL[jj],AB);
		if(ci != cI || cj != cJ){
		   AlnSeqSW(stderr,11, 1,dcaE, E, AB);
		   PutSeq(stderr,E,AB); PutSeq(stderr,xE,AB); PutSeq(stderr,dcaE,AB); 
		   fprintf(stderr,"%c%d == %c%d?\n",ci,ResidueID(ResALL[ii]),cI,i);
		   fprintf(stderr,"%c%d == %c%d?\n",cj,ResidueID(ResALL[jj]),cJ,j);
		   fprintf(stderr,"Make sure that the EVfold query sequence matches the pdb sequence.\n");
		   print_error("GetMtrxDCA() input error 3.");
		} else {
		  // EqSeq(i,AlphaCode(cI,AB),xE); EqSeq(j,AlphaCode(cJ,AB),xE);
		}
		// fprintf(stderr,"%s\n",Str);
           	if(ii < resStart || ii > resEnd) continue;
           	if(jj < resStart || jj > resEnd) continue;
		DM[ii][jj]=DM[jj][ii]=(float) Dd;
		Used[ii][jj]=Used[jj][ii]=TRUE;
	} fclose(fp); // PutSeq(stderr,xE,AB);
	NilSeq(E); NilSeq(xE);
	for(i=1; i <= num_resC; i++){
	   ResI=ResALL[i]; I=ResidueID(ResI);
           if(I < resStart || I > resEnd) continue;
	   for(j=i+1; j <= num_resC; j++){
	     ResJ=ResALL[j]; J=ResidueID(ResJ);
             if(J < resStart || J > resEnd) continue;
	     DstMtx[i][j]=DstMtx[j][i]=-DM[i][j]; 
	     if(efptr) IncdHist(DM[i][j],HG);
	   }
	} if(efptr){ PutHist(stderr,60,HG); NilHist(HG); }
	for(i=1; i <= max_res; i++) free(DM[i]); free(DM);
	return DstMtx;
}

float	**scm_typ::CalcDistMtrx()
//============= Move this routine to libpdb.a at some point: ================
{
	Int4	ai,aj,r,i,j;

	// Construct minimum Res-Res distance matrix...
	h_type HG=0;
	if(efptr) HG=Histogram("residue distances",0,100,1.0);
	atm_typ	atmI,atmJ;
	res_typ	ResI,ResJ;
	char	cI,cJ;
	float	max,min,dd;	// 350 x 350 = 122,500 x 16 = 1,470 kBytes.
	float	**DstMtx;
	double	Dd;
	NEWP(DstMtx, num_resC + 3, float);
	for(i=1; i <= num_resC; i++){ NEW(DstMtx[i], num_resC + 3, float); }
	for(i=1; i <= num_resC; i++){
	   ResI=ResALL[i]; 
	   Int4 site=ResidueID(ResI);
	   if(site < resStart || site > resEnd) continue;
	   cI=GetCharResidue(ResI,AB);
	   for(j=i + 1; j <= num_resC; j++){
	     ResJ=ResALL[j]; site=ResidueID(ResJ);
	     if(site < resStart || site > resEnd) continue;
	     cJ=GetCharResidue(ResJ,AB); dd=0.0; min=999999.0;
	     if(UseHydrogens){
		Int4	h,ii,jj;
		// create arrays of all atoms...
		atm_typ atmsI[200],atmsJ[200],*H;
		for(ii=0,ai=1; ai <= ResidueAtomNumber(ResI); ai++){
		   atmI=AtomResidue(ai,ResI);
		   if(!(SideAtom(atmI) || (cI=='G' && AlphaCarbonAtom(atmI)))) continue;
		   ii++; atmsI[ii]=atmI;
		   H=ResidueAtomHydrogens(ai,ResI);
		   if(H) for(h=1; H[h]; h++){ ii++; atmsI[ii]=H[h]; } 
		} ii++; atmsI[ii]=0;
		for(jj=0,aj=1; aj <= ResidueAtomNumber(ResJ); aj++){
		   atmJ=AtomResidue(aj,ResJ);
		   if(!(SideAtom(atmJ) || (cJ=='G' && AlphaCarbonAtom(atmJ)))) continue;
		   jj++; atmsJ[jj]=atmJ;
		   H=ResidueAtomHydrogens(aj,ResJ);
		   if(H) for(h=1; H[h]; h++){ jj++; atmsJ[jj]=H[h]; }
		} jj++; atmsJ[jj]=0;
		for(ai=1; atmsI[ai]; ai++){
		   for(atmI=atmsI[ai],aj=1; atmsJ[aj]; aj++){
		      dd=DistanceAtoms(atmI,atmsJ[aj]);
		      if(dd < min) min=dd;
		   }
		}
	    } else {
		for(ai=1; ai <= ResidueAtomNumber(ResI); ai++){
		  atmI=AtomResidue(ai,ResI);
		  if(!(SideAtom(atmI) || (cI=='G' && AlphaCarbonAtom(atmI)))) continue;
		  //fprintf(stderr,"%s: \n",ResidueName(ResI)); PutAtom(stderr,atmI);
		  for(aj = 1; aj <= ResidueAtomNumber(ResJ); aj++){
		    atmJ=AtomResidue(aj,ResJ);
		    if(!(SideAtom(atmJ) || (cJ=='G' && AlphaCarbonAtom(atmJ)))) continue;
		    dd=DistanceAtoms(atmI,atmJ);
		    if(dd < min) min=dd;
		  }
		}
	    } DstMtx[i][j]=DstMtx[j][i]=min; if(efptr) IncdHist(min,HG);
	   }
	} if(efptr){ PutHist(stderr,60,HG); NilHist(HG); }
	return DstMtx;
}

float	*scm_typ::CalcBuried( )
{
        Int4 r,i,j,n,num_resA,num_resD,C2=0,NumResC,num_resO,*num_resALL;
        res_typ *ResCD,**ResEvery;
        NEWP(ResEvery,nChainsPDB(P)+3,res_typ);
        NEW(num_resALL,nChainsPDB(P)+3,Int4);
        for(C2=1; C2 <= nChainsPDB(P); C2++){
               ResEvery[C2] = MakeResPDB(C2,&num_resALL[C2],P);
        }
	assert(num_resC==num_resALL[C]); 
	ResCD=ResEvery[C];
	Int4 NumInChain[MAX_NGRPS_PDB];
	assert(nChainsPDB(P) < MAX_NGRPS_PDB); 
	for(C2=1; C2 <= nChainsPDB(P); C2++) NumInChain[C2]=MaxAtomsPDB(C2,P);
        for(C2=1; C2 <= nChainsPDB(P); C2++){
	   char chainB=ChainCharPDB(C2,P);
	   if(strchr(ChainI,chainB) != NULL){
		Int4 CB=GetChainNumberPDB(P,chainB); // C = query residue chain.
		NumInChain[CB]=0;  // compute without this subunit.
	   }
	}
	double buried; NEW(Buried,num_resC+9,float);
	for(n=0,j=1; j <= num_resC; j++){
		double free=CalculateResidueSurface(ResCD[j],nChainsPDB(P),NumInChain,
					AtomsPDB(P));
		double bound=CalculateResidueSurface(ResCD[j],nChainsPDB(P),NumAtomsPDB(P),
					AtomsPDB(P));
		buried = free-bound;
		if(buried >= min_buried){ Buried[j]=buried; n++; }
	}
	fprintf(stderr,"Number buried = %d; min_buried=%.2f\n",n,min_buried);
	// if(n < 1) print_error("Minimum surface area set too high; reset -mb option.");
#if 0
	dh_type dH = dheap(num_resC+3,4);
	for(j=1; j <= num_resC; j++){
		if(Buried[j] >= min_buried) insrtHeap(j,(-(keytyp)Buried[j]),dH);
	}
	a_type AB = AminoAcidAlphabetPDB(P);
	while((buried=(double)-minkeyHeap(dH)) >= 0.0 && (j=delminHeap(dH)) != 0){
		fprintf(stderr,"#%c%d%c.W // %.2f buried\n",GetCharResidue(ResCD[j],AB),
					ResidueID(ResCD[j]),chain,buried);
	} Nildheap(dH);
#endif
        for(C2=1; C2 <= nChainsPDB(P); C2++){
                for(i=1; i <= num_resALL[C2]; i++) NilRes(ResEvery[C2][i]);
                free(ResEvery[C2]);
        } free(ResEvery); free(num_resALL);
	if(n==0){ free(Buried); Buried=0; }
	return Buried;
}

void    scm_typ::Init(scl_typ *S)
{
         // don't copy over these unused parameters;
         time1=0; PutPvals=FALSE; seed=0; vsifp=0; SetI=0;  tab_fp=0;
         KeyAtm=0;KeyAtmDist=0;

         // pass over these parameters;
         AB=S->AB;
         UnionCons=S->UnionCons;
         PutAllSCH=S->PutAllSCH;
         KeyAtmNum=S->KeyAtmNum;
         AroPiBonds=S->AroPiBonds;
         OtherPiBonds=S->OtherPiBonds;
         program=S->program; infile=S->infile;
         dca_vs_pdb=S->dca_vs_pdb;
         dca_file=S->dca_file;
         dcaE=S->dcaE; dcaOS=S->dcaOS;
         OnlyOne=S->OnlyOne; Ignore=S->Ignore; shuffle=S->shuffle;
         KeyRes=S->KeyRes; KeyChain=S->KeyChain; KeyResX=S->KeyResX;
         FindHbonds=S->FindHbonds; PutTheRest=S->PutTheRest;
         ChainI=S->ChainI; FullCluster=S->FullCluster;
         BestOnly=S->BestOnly;
         BackboneHbonds=S->BackboneHbonds; min_buried = S->min_buried;
         MaxPrint=S->MaxPrint; MaxNumPttrns=S->MaxNumPttrns; method=S->method;
         MaxMinAdj=S->MaxMinAdj; EdgeDrop=S->EdgeDrop;
         cutoff=S->cutoff;
         vsi_file=S->vsi_file; Simulate=S->Simulate; UseHydrogens=S->UseHydrogens;
         UseJeffreys=S->UseJeffreys;
         xHG=S->xHG; zHG=S->zHG; xxHG=S->xxHG;
         efptr=S->efptr;

	if(chain) C=GetChainNumberPDB(P,chain); // C = query residue chain.
	if(C==0){ fprintf(stderr,"chain = %c\n",chain); print_error("chain not found!"); }
	if(chain==0) chain='?';
	a_type ab=AminoAcidAlphabetPDB(P); 
	ResALL=MakeResPDB(C,&num_resC,P);
	if(KeyResX){ 	// Use virtual atom at center of residue.
	   // PrintResidueAtoms(stderr,KeyResX);
	   Int4 i,j,ai,aj,h,ii,jj;
	   MaxPrint=1; 
	   double dd,d_min=999999999999.9,D_sum,D_ave,D_min=999999999999.9;
	   if(KeyAtmNum > 0){	// Add an atom identifier later...??
		atm_typ atmX=AtomResidue(KeyAtmNum,KeyResX); assert(atmX != 0);
	   	KeyAtm=CopyCoreAtom(atmX); PutAtom(stderr,KeyAtm);
	   } else KeyAtm=MkAveAtomResidue(KeyResX);
	   NEW(KeyAtmDist,num_resC+5,float);
	   atm_typ atmI,atmJ;
	   res_typ ResI;
	   for(i=1; i <= num_resC; i++){
	      ResI=ResALL[i]; 
	      Int4 site=ResidueID(ResI);
	      // fprintf(stderr,"jj=%d; j=%d; site=%d; start=%d; end=%d\n",jj,j,site,resStart,resEnd);
	      if(site < resStart || site > resEnd) continue;
	      if(!MemberSet(site,resSet)) continue;
	      // Int4 site=ResidueID(ResI); if(site < resStart || site > resEnd) continue;
	      // char cI=GetCharResidue(ResI,ab);
	      for(ii=0,D_sum=0,ai=1; ai <= ResidueAtomNumber(ResI); ai++){
		 atmI=AtomResidue(ai,ResI);
		 dd=DistanceAtoms(atmI,KeyAtm);
		 D_sum += dd; ii++;
		 if(d_min > dd) d_min=dd;
	      } D_ave= D_sum/(double)ii;
	      KeyAtmDist[i]=D_ave;
	      if(D_ave < D_min){ D_min=D_ave; }
	   }
	   // fprintf(stderr,"D_min=%.3f; d_min=%.3f.\n",D_min,d_min);
	   if(d_min > 5.0){ NilAtom(KeyAtm); free(KeyAtmDist); KeyAtmDist=0; KeyAtm=0; }
	} Buried=0; DistMtrx=0;
	if(FindHbonds) DistMtrx=this->GetHbondDistMtrx(); 
	else if(dca_file!=0){
	    BooLean **used=0; DistMtrx=GetMtrxDCA(dca_file,used);
	    for(Int4 j=1; j <= num_resC; j++) free(used[j]); free(used);
	} else if(method < 6) DistMtrx=this->CalcDistMtrx();
	if(method==5){
	  for(Int4 j=1; j <= num_resC; j++){
		CalculateResidueSurface(ResALL[j],nChainsPDB(P),NumAtomsPDB(P),AtomsPDB(P));
	  }
	}
	// fprintf(stderr,"shuffle=%d\n",shuffle);
	if(method >= 6){ Buried=CalcBuried( ); }
	if(Simulate){ if(shuffle) Shuffle(); else Permute(); } // randomize pattern sites.
	SetX=MakeSet(num_resC+5);
}

Int4	scm_typ::Shuffle()
{
	Int4	i,j,r,jj,*PttrnSites;
	for(jj=1; PttrnRes[jj] != 0;jj++){
	      PttrnSites = PttrnRes[jj];
// fprintf(stderr,"DOING SIMULATION\n");
	      dh_type dH0=dheap(num_resC+5,4);
	      for(j=1; j <= num_resC; j++){
		Int4 site=ResidueID(ResALL[j]);
	        // fprintf(stderr,"jj=%d; j=%d; site=%d; start=%d; end=%d\n",jj,j,site,resStart,resEnd);
		if(site < resStart || site > resEnd) continue;
		if(!MemberSet(site,resSet)) continue;
		keytyp kk=(keytyp)Random(); insrtHeap(j,kk,dH0);
	      }
	      if(emptyHeap(dH0)){
		fprintf(stderr,"WARNING: possible residue numbering inconsistency\n");
	      } else for(r=1; PttrnSites[r] > 0; r++){
	 	  j=delminHeap(dH0); 
		  // fprintf(stderr,"j=%d; r=%d\n",j,r);
		  PttrnSites[r] = ResidueID(ResALL[j]);
	      } Nildheap(dH0);
	} return jj-1;
}

Int4	scm_typ::Permute()
// circularly permute pattern sites.
{
	Int4	Ni,i,j,I,J,r,s,x,jj,*PttrnSites,*Site;
	FILE	*efp=0; // efp=stderr; 
	static Int4 os=0;
	const Int4 os_omit=3;

	// 1. Create an array of all sites...
	NEW(Site,num_resC +9,Int4);
	for(Ni=0,j=1; j <= num_resC; j++){
	    Int4 site=ResidueID(ResALL[j]);
	    if(site < resStart || site > resEnd) continue;
	    if(!MemberSet(site,resSet)) continue;
	    Ni++; Site[Ni]=j;
	}

	// 2. For each category move the patterns by a fixed number of positions.
	do { os=((os+1)%Ni); } while(os==(Ni-1) || os < os_omit || (Ni-os) < os_omit);
	// do { os=((os+1)%Ni); } while(os==(Ni-1));
	// do { os=(Random()%Ni); } while(os==(Ni-1));
	if(1 || efp) fprintf(stderr,"os=%d; Ni=%d\n",os,Ni);
	for(jj=1; PttrnRes[jj] != 0;jj++){
	    PttrnSites = PttrnRes[jj];
	    if(efp){
	      for(r=1; (J=PttrnSites[r]) > 0; r++) fprintf(stderr,"%3d ",J); fprintf(stderr,"\n");
	    }
	    for(r=1; (J=PttrnSites[r]) > 0; r++){
		for(i=1; i <= Ni; i++){
		   j=Site[i]; I=ResidueID(ResALL[j]);
		   if(J==I){
			i = ((i+os)%Ni)+1;
			j=Site[i];
			x=ResidueID(ResALL[j]);
			if(efp) fprintf(stderr,"%3d ",x-PttrnSites[r]);
			PttrnSites[r]=x; 
	    		assert(MemberSet(x,resSet)); break;
		   }
		} assert(i <= Ni);
	    } if(efp) fprintf(stderr,"\n");
	    if(efp){
	      for(r=1; (J=PttrnSites[r]) > 0; r++) fprintf(stderr,"%3d ",J); fprintf(stderr,"\n\n");
	    }
	} free(Site);
	return Ni;
}

BooLean	scm_typ::Run(FILE *fptr, FILE *vsifptr, FILE *tfp) 
/************************************************************
Find strings using method 1-6.
 ************************************************************/
{
	if(KeyRes && KeyAtm == 0) return FALSE;	// skip these...
	if(method >= 6 && Buried==0) return FALSE;	// skip these...
#if 1   // print out categories for Aravind: 'R' 'I' 'T' ... root one two...
        char *aStr=0; 
	if(dcaE){
	  NEW(aStr,SetN(SetX)+4,char);
          for(Int4 s=0; s < SetN(SetX); s++) aStr[s]='-';
	}
#endif
	for(Int4 jj=1; PttrnRes[jj] != 0;jj++){
	   Int4	r,i,j,StartRes,ii,*PttrnSites,NumToCheck;
	   char	cI,cJ;
	   if(Ignore == res_class[jj]) continue;
	   if(OnlyOne && OnlyOne != res_class[jj]) continue;
	   PttrnSites = PttrnRes[jj]; ClearSet(SetX);
	   NumToCheck=MINIMUM(Int4,MaxPrint,PttrnSites[0]); assert(NumToCheck <= MaxPrint);
	   Int4 *Map; NEW(Map,num_resC +4, Int4);
	   for(j=1; j <= num_resC; j++){
		for(r=1; PttrnSites[r] > 0; r++){
		    if(ResidueID(ResALL[j])==PttrnSites[r]){ AddSet(j,SetX); Map[j]=r; break; }
		} // if(Map[j] != 0) fprintf(stdout,"site[%d] = %d\n",j,Map[j]);
	   } // PutSet(stdout,SetX);
	   // Hbond networks... need to add simulation option and Map[j].
	   if(isupper(FindHbonds)){	// use all higher level residues for network...
	    for(Int4 jjj=1; jjj < jj; jjj++){
	     Int4 *PttrnSts = PttrnRes[jjj];
	     for(j=1; j <= num_resC; j++){
		for(r=1; PttrnSts[r] > 0; r++){
		    if(ResidueID(ResALL[j])==PttrnSts[r]){ AddSet(j,SetX); break; }
		}
	     }
	    }
	   }
	   Int4	xx,II=0,JJ,JJ_min,nF,nB,s,j_min,Dist,bestm,bestr,D;
	   double N,prob;

if(FullCluster == res_class[jj]) NumToCheck=1; // use with -f= option...print best only..
if(BestOnly) NumToCheck=1; // -first option...print best only..
#if 1	// count number of sites.
	   for(D=0,r=1; PttrnSites[r] > 0; r++) D++;
#endif

#if 0
// fprintf(stdout,"SetX = "); PutSet(stdout,SetX); // NOTE: invisible residues are not counted...
	   sch_typ *sch = new sch_typ(MaxPrint,D,resStart,resEnd,res_class[jj],SetX,ResALL,AB);
#elif 1
	   sch_typ *sch = new sch_typ(NumToCheck,D,resStart,resEnd,res_class[jj],SetX,ResALL,AB);
#else
	   sch_typ *sch = new sch_typ(NumToCheck,MaxNumPttrns,resStart,resEnd,res_class[jj],
					SetX,ResALL,AB);
#endif
	   for(ii=1; ii <= MaxPrint; ii++) {
	     StartRes=PttrnSites[ii];   // with -K= option this is ignored...
	     if(StartRes==0) break;
	     for(II=0,j=1; j <= num_resC; j++) if(ResidueID(ResALL[j]) == StartRes){ II=j; break; }
	     if(II <= 0) continue;	// invisible residues..
	     Int4 MaxStr=0,*String; NEW(String, num_resC+9,Int4); 
	     char *hits; NEW(hits,num_resC+9,char); 
	     s=1; String[s]=II; hits[0]='*';

	     if(method==1){
		prob=NetworkClust(II,StartRes,hits,String,MaxStr,bestm,Dist);
	     } else if(method==2 || (method==4 && res_class[jj] != 'M')){
		prob=CoreClust(II,StartRes,hits,String,MaxStr,bestm,Dist);
	     } else if(method==3 || (method==4 && res_class[jj] == 'M')){
		prob=SphericalClust(II,StartRes,hits,String,MaxStr,bestm,Dist);
	     } else if(method==5){ 	// CalcSurfaceAccess(): use surface to define set.
		prob=SurfaceClust(II,jj, hits,String,MaxStr,bestm,Dist);
	     } else if(method==6){	// Subunit interface-based fixed set.  uses CalcBuried():
		prob=FixedInterface(hits,String,MaxStr,bestm,Dist);
	     } else if(method==7){	// Subunit interface-based fixed set.  uses CalcBuried():
		prob=OptInterfaceClust(hits,String,MaxStr,bestm,Dist);
	     } else if(method==8){	// Subunit interface-based fixed set.  uses CalcBuried():
		prob=SphericalClust(II,StartRes,hits,String,MaxStr,bestm,Dist);
	     } else print_error("illegal method");
	     // nF=CardSet(SetX); nB = nAln - nF; N=(double)(nF+nB); 
	     // if((efptr || FindHbonds) && prob <= cutoff)
	     if(efptr && prob <= cutoff){
	       fprintf(fptr,"%d: D=%d; bestm=%d; maxstr=%d; L=%d\n",ii,Dist,bestm,MaxStr,nAln);
	       fprintf(fptr," str=\"%s\"\n",hits);
	     } free(hits);
	     char rc='X';
	     // if(prob <= cutoff && MaxStr > 0)
	     if(MaxStr > 0)
	     {
		Int4 xxx=sch->Insert(ii,bestm,MaxStr,String,prob); 
		// sch->PutItem(stderr,xxx);
	     } else {
	       if(efptr && MaxStr > 0){
		 rc=GetCharResidue(ResALL[II],AB);
	         fprintf(fptr,"%d: site=%c%d; p = %.3g: \n",ii,rc,StartRes,prob);
	       } fflush(fptr); free(String);
	     }
	     if(zHG){ IncdHist(bestm,zHG); }
	     if(xHG) IncdHist(prob,xHG); // long double dd = -expm1l(-prob); // expm1(x) = exp(x) - 1
	     if(xxHG){ if(prob < 0.025) IncdHist(prob,xxHG); }
	} // end ii loop.
	set_typ tmpSet=sch->Put(fptr,vsifptr,efptr,tfp,FullCluster,aStr,cutoff,PutAllSCH,UnionCons);
	if(tmpSet){
// PutSet(stderr,tmpSet);
	  if(PutTheRest == res_class[jj]){ assert(SetI==0); SetI=tmpSet; }
	  else NilSet(tmpSet); tmpSet=0;
	}
	delete sch; free(Map);
    } // end of jj loop..
#if 1   // print out categories for Aravind: 'R' 'I' 'T' ... root one two...
    if(aStr){
	Int4 Ns=0;
        for(Int4 s=0; aStr[s] != 0; s++) if(aStr[s]!='-') Ns++;
	fprintf(stdout,"\naStr(%d)=\"%s\"\n\n",Ns,aStr); free(aStr); 
    }
#endif
    //========== Output non-intersecting discriminating residues. =========
    if(vsifptr && PutTheRest && SetI && CardSet(SetI) > 0){
	for(Int4 jj=1; PttrnRes[jj] != 0;jj++){
	   if(PutTheRest != res_class[jj]) continue;
	   if(Ignore == res_class[jj]) continue;
	   if(OnlyOne && OnlyOne != res_class[jj]) continue;
	   Int4 *PttrnSites = PttrnRes[jj],nput=0; 
	   for(Int4 r=1; PttrnSites[r] > 0; r++){
	      for(Int4 j=1; j <= num_resC; j++){
		if(!MemberSet(j,SetI)) continue;
		if(ResidueID(ResALL[j])==PttrnSites[r]){
		   char rc = GetCharResidue(ResALL[j],AB);
		   if(nput > 0) fprintf(vsifptr,",");
		   fprintf(vsifptr,"%c%d.L",rc,PttrnSites[r]); nput++;
		   break; 
		}
	      }
	   } fprintf(vsifptr,"\n");
	}
    }
    if(vsifptr && ChainI){
	for(Int4 i=0; ChainI[i] != 0; i++){ fprintf(vsifptr,"\n1-1000%c.L80\n",ChainI[i]); }
    }
    return TRUE;
} 

double	scm_typ::FixedInterface(char *hits,Int4 *String, Int4 &L, Int4 &bestm, Int4 &D)
#if 0	//****************************************************
 Simple Ball-In-Urn model:
	N total balls with N1 red balls and N2 = N-N1 black balls.
	Choose n balls at random.  The probability that the
	group so chosen will contain x or more red balls is
	given by: p=CumHypGeomProb(N1,N2,n,x).
	    red == at interface; black = not at interface.
	out = discriminating;  in = non-discriminating.
#endif	//****************************************************
{
	Int4	i,j,s,ri=0,ro=0,bi=0,bo=0,JJ;
	char	cJ;
	double	prob;
	h_type HG=0; // HG=Histogram("Pattern residue surface access",0,1,0.02);
	dh_type dH = dheap(num_resC+5,4);
	for(j=1; j <= num_resC; j++){
	  Int4 site=ResidueID(ResALL[j]);
	  if(site < resStart || site > resEnd) continue;
	  if(!MemberSet(site,resSet)) continue;		// Note: some residues are invisible...
	  double dd=Buried[j]; 
	  if(dd <= 0.0){	// not at interface.
		if(MemberSet(j,SetX)) bo++; else bi++;
	  } else {		// at interface.
		if(MemberSet(j,SetX)) ro++; else ri++;
		insrtHeap(j,(keytyp)-dd,dH);	// all at interface...
	  } if(HG) IncdHist(dd,HG); 
	} if(HG){ PutHist(stderr,60,HG); NilHist(HG); }
	for(s=0,i=0; !emptyHeap(dH); i++) {	// these are at interface.
		double dd = (double) minkeyHeap(dH);
		assert((j=delminHeap(dH)) != 0);
		if(MemberSet(j,SetX)) hits[s]='*'; else hits[s]=' ';
		if(efptr && hits[s]=='*'){
	           cJ=GetCharResidue(ResALL[j],AB); JJ=ResidueID(ResALL[j]);
		   fprintf(efptr,"%d: %c%d (%.3f)%c\n",i,cJ,JJ,dd,hits[s]);
		} s++; String[s]=j; 
	} Nildheap(dH); nAln=s;
	bestm=ro; D=ro; L=s;
	prob = CumHypGeomProb(ro+ri,bo+bi,ro+bo,ro);
#if 0
	Int4	out=0,in=0;
	out=bo+ro; in=bi+ri;
	fprintf(fptr,"nAln=%d; s=%d; CardSet=%d\n",nAln,s,CardSet(resSet));
	fprintf(fptr,"%d: hits=\"%s\"; prob=%.3g.\n",ii,hits,prob);
	fprintf(fptr," ro=%d; ri=%d; bo=%d; bi=%d; in=%d; out=%d.\n",ro,ri,bo,bi,in,out);
	// PutSet(fptr,resSet);
#endif
	return prob;
}

double  scm_typ::SphericalClust(Int4 II,Int4 StartRes,char *hits,Int4 *String, Int4 &L, Int4 &bestm, Int4 &D)
// Method 3 (and method 4 + 'M').
{
	char    cJ,rc;
	double	prob;
	Int4    i,j,s,JJ,bestr;

	if(KeyAtmDist) { s=1; String[s]=0; hits[0]=' '; II=0; }
	else { s=1; String[s]=II; hits[0]='*'; }
	dh_type dH = dheap(num_resC+5,4);
	for(j=1; j <= num_resC; j++){
	   Int4 site=ResidueID(ResALL[j]);
	   if(site < resStart || site > resEnd) continue;
	   if(!MemberSet(site,resSet)) continue;
	   if(KeyAtmDist) insrtHeap(j,(keytyp)KeyAtmDist[j],dH);
	   else insrtHeap(j,(keytyp)DistMtrx[j][II],dH); 
	}
	for(i=1; !emptyHeap(dH); i++) {
		double dd = (double) minkeyHeap(dH);
		assert((j=delminHeap(dH)) != 0);
		if(j==II) continue;
		if(FindHbonds && dd > 4.0) continue;
		if(MemberSet(j,SetX)) hits[s]='*'; else hits[s]=' ';
		if(efptr && hits[s]=='*'){
	           cJ=GetCharResidue(ResALL[j],AB); JJ=ResidueID(ResALL[j]);
		   fprintf(efptr,"%d: %c%d (%.3f)%c\n",i,cJ,JJ,dd,hits[s]);
		} s++; String[s]=j; 
	} Nildheap(dH); nAln=s-1;	// s++ at the end is not used.
	if(efptr){
		PutSet(efptr,SetX); 
		fprintf(efptr,"nAln=%d; s=%d; CardResSet=%d.\n",nAln,s,CardSet(resSet));
	}

        Int4       *pos; NEW(pos, nAln +6, Int4);
        for(i=1,s=0; i <= nAln; i++) if(hits[i]=='*'){ s++; pos[s]=i; } // else assert(hits[i] == ' ');  
        D=s;
        bestm=pvcalc(nAln,D,pos,prob,UseJeffreys);
        if(bestm == 0){
                if(efptr){ 
			if(II > 0) rc=GetCharResidue(ResALL[II],AB); else rc='X'; 
			fprintf(stderr,"%c%d: (null is best) D=%d; nAln=%d; prob=%3f.\n",rc,StartRes,D,nAln,prob);
		}
                bestr=0;
        } else if(bestm < 0){
                // if(D > 0) fprintf(stderr,"Error: X=%d; D = %d; L=%d; prob=%3f.\n",bestm,D,nAln,prob);
                bestm=0; prob=1.0;
                bestr=pos[1];
                // print_error("pvcalc input error");
        } else  bestr=pos[bestm];
#if 1	// Use for KeyAtm (-K=) option.
	if(hits[0]=='*'){ bestm++; } 
#endif
	bestr++; L=bestr; hits[L] = 0; free(pos); 
	return prob;
}

double  scm_typ::CoreClust(Int4 II,Int4 StartRes, char *hits,Int4 *String, Int4 &L, Int4 &bestm, Int4 &D)
{
	Int4    xx,i,j,s,JJ,bestr,j_min;
	char    cJ,cI;
	double	prob;

	// Create and initialize the residue cluster set with the starting site.
        set_typ SetC=MakeSet(num_resC+5); ClearSet(SetC);
	double	DD_max=0.0;	// max distance used thus far
	Int4	z,i_min,MinAdj=0,NumInClst=1,*NumAdj=0; NEW(NumAdj, num_resC + 3, Int4);
#if 1
	s=1; String[s]=II; hits[0]='*'; AddSet(II,SetC); 
	// Find the Optimum cluster of residues...adding each next closest residue to the set.
	set_typ SetE=CopySet(SetC);	// SetE == excluded sequences from consideration.
	for(j=1; j <= num_resC; j++){
		if((JJ=ResidueID(ResALL[j])) <= 0) AddSet(j,SetE); 	// invisible.
	        else if(JJ > resEnd) AddSet(j,SetE);	// out of range...
	        else if(!MemberSet(JJ,resSet)) AddSet(j,SetE); 	// insert region.
	        // else if(JJ < resStart) AddSet(j,SetE);	// out of range...
	}
#else	// extended interface clustering...Needs much more work...
	Buried=CalcBuried( );
	for(j=1; j <= num_resC; j++){
	  Int4 site=ResidueID(ResALL[j]);
	  if(site < resStart || site > resEnd) continue;
	  if(!MemberSet(site,resSet)) continue;		// Note: some residues are invisible...
	  double dd=Buried[j]; 
	  if(dd <= 0.0){	// not at interface.
		if(MemberSet(j,SetX)) bo++; else bi++;
	  } else {		// at interface.
		if(MemberSet(j,SetX)) ro++; else ri++;
	  } 
	} 
#endif
	for(xx=1; CardSet(SetE) < num_resC && xx <= num_resC; xx++){
	   double DD_min,DD;
	   float  *distmtx=0;
	   Int4 MinEdge=MinAdj;
	   // if(MinEdge > EdgeDrop) MinEdge--;
	   j_min=0; DD_min=9999999.9;
	   for(i=1; i <= num_resC; i++){	// i are members of cluster set.
	      if(!MemberSet(i,SetC)) continue;
	      if(NumAdj[i] < MinEdge) continue;	// use only 'internal' residues...
	      distmtx = DistMtrx[i];
	      for(j=1; j <= num_resC; j++){	// j are not members of cluster set.
	         if(MemberSet(j,SetE)) continue; 
	   	 if(distmtx[j] < DD_min){ j_min=j; i_min=i; DD_min=distmtx[j]; }
	      }
	   } assert(j_min != 0);
	   if(j_min == 0) continue; 	// nothing found???
	   if(FindHbonds && DD_min > 4.0) continue;
	   j=j_min; i=i_min; 
	   if(MemberSet(j,SetX)){
		 hits[s]='*';
	   } else hits[s]=' ';  // hits[s]=String[s+1].
	   if(DD_min > DD_max) DD_max=DD_min;

	   AddSet(j,SetC); AddSet(j,SetE); NumInClst++;
	   for(Int4 y=1; y <= num_resC; y++) NumAdj[y]=0;
	   for(Int4 y=1; y < num_resC; y++){
	      if(!MemberSet(y,SetC)) continue;
	      for(distmtx=DistMtrx[y],z=y+1; z <= num_resC; z++){
	        if(!MemberSet(z,SetC)) continue;
	        if(distmtx[z] <= DD_max){ NumAdj[z]++; NumAdj[y]++;  }
	      }
	   }
	   if(MinAdj < MaxMinAdj){
	     for(z=1; z <= num_resC; z++){ if(MinAdj < NumAdj[z]) MinAdj=NumAdj[z]; }
	   }
	   if(0 || efptr && hits[s]=='*'){
		Int4 rI=ResidueID(ResALL[i]); cI=GetCharResidue(ResALL[i],AB);
		JJ=ResidueID(ResALL[j]); cJ=GetCharResidue(ResALL[j],AB);
		fprintf(stderr,"%d: %c%d (%.3f)%c [%c%d;%d] %.3f A; %d\n",
			xx,cJ,JJ,DD_min,hits[s],cI,rI,NumAdj[i],DD_max,MinAdj);
	   } s++; String[s]=j; 
	} nAln=s-1;
	// fprintf(stderr,"nAln=%d; num_resC=%d; CardSetE=%d.\n", nAln,num_resC,CardSet(SetE));
	NilSet(SetC); NilSet(SetE); free(NumAdj);

        Int4       *pos; NEW(pos, nAln +6, Int4);
        for(i=1,s=0; i <= nAln; i++) { if(hits[i]=='*'){ s++; pos[s]=i; } }
        D=s;
        bestm=pvcalc(nAln,D,pos,prob,UseJeffreys);
        if(bestm == 0){
                char rc=GetCharResidue(ResALL[II],AB);
                if(efptr) fprintf(stderr,"%c%d: (null is best) D=%d; nAln=%d; prob=%3f.\n",
                        rc,StartRes,D,nAln,prob);
                bestr=0;
        } else if(bestm < 0){
                // if(D > 0) fprintf(stderr,"Error: X=%d; D = %d; L=%d; prob=%3f.\n",bestm,D,nAln,prob);
                bestm=0; prob=1.0;
                bestr=pos[1];
                // print_error("pvcalc input error");
        } else  bestr=pos[bestm];
        free(pos); bestm++; bestr++; L=bestr; hits[L] = 0;
        return prob;
}

double  scm_typ::NetworkClust(Int4 II,Int4 StartRes, char *hits,Int4 *String, Int4 &L, Int4 &bestm, Int4 &D)
{
        Int4    xx,i,j,s,JJ,bestr,j_min;
        char    cJ,cI,rc;
        double  prob;
	set_typ	SetC=MakeSet(num_resC+5); ClearSet(SetC); 

	// Create and initialize the residue cluster set with the starting site.
	if(KeyAtmDist){ s=1; String[s]=0; hits[0]=' '; II=0; }
	else { s=1; String[s]=II; hits[0]='*'; AddSet(II,SetC); }

	// Find the Optimum cluster of residues...adding each next closest residue to the set.
	set_typ SetE=CopySet(SetC);	// SetE == excluded sequences from consideration.
	for(j=1; j <= num_resC; j++){
		if((JJ=ResidueID(ResALL[j])) <= 0) AddSet(j,SetE); 	// invisible.  
		else if(JJ > resEnd) AddSet(j,SetE);	// out of range...
	        else if(!MemberSet(JJ,resSet)) AddSet(j,SetE); 	// insert region.
	        // else if(JJ < resStart) AddSet(j,SetE);	// out of range...
	}
	for(xx=1; CardSet(SetE) < num_resC && xx <= num_resC; xx++){
	   double DD_min,DD,dd;
	   j_min=0; DD_min=9999999.9;
	   for(i=1; i <= num_resC; i++){	// i are members of cluster set.
	     if(s==1 && KeyAtmDist){	// for -K=<int> option...
	        for(j=1; j <= num_resC; j++){	// j are not members of cluster set.
	            if(MemberSet(j,SetE)) continue; 
		    if((dd=KeyAtmDist[j]) < DD_min){ j_min=j; DD_min=dd; }
		}
	     } else {
	        if(!MemberSet(i,SetC)) continue;
	        float *distmtx = DistMtrx[i];
	        for(j=1; j <= num_resC; j++){	// j are not members of cluster set.
	            if(MemberSet(j,SetE)) continue; 
	   	    if(distmtx[j] < DD_min){ j_min=j; DD_min=distmtx[j]; }
	        }
	     }
	   } if(j_min == 0) continue; 	// nothing found???
	   j=j_min;
	   if(!(s==1 && KeyAtmDist)){	// allow larger distance at the start in this case.
		if(FindHbonds && DD_min > 4.0) continue;
	   } if(MemberSet(j,SetX)) hits[s]='*'; else hits[s]=' ';
#if 1	// DEBUG...
	   if(s==1 && KeyAtmDist){
		JJ=ResidueID(ResALL[j]); cJ=GetCharResidue(ResALL[j],AB);
		fprintf(stderr,"Start: %d: %c%d (%.3f)%c\n",xx,cJ,JJ,DD_min,hits[xx]);
	   }
#endif
	   AddSet(j,SetC); AddSet(j,SetE);
	   if(efptr && hits[s]=='*'){
		JJ=ResidueID(ResALL[j]); cJ=GetCharResidue(ResALL[j],AB);
		fprintf(stderr,"%d: %c%d (%.3f)%c\n",xx,cJ,JJ,DD_min,hits[xx]);
	   } s++; String[s]=j; 
	} nAln=s-1;
	// fprintf(stderr,"nAln=%d; num_resC=%d; CardSetE=%d.\n", nAln,num_resC,CardSet(SetE));
	NilSet(SetC); NilSet(SetE);

        Int4       *pos; NEW(pos, nAln +6, Int4);
        for(i=1,s=0; i <= nAln; i++) { if(hits[i]=='*'){ s++; pos[s]=i; } }
	D=s;
        bestm=pvcalc(nAln,D,pos,prob,UseJeffreys);
        if(bestm == 0){
                if(II > 0) rc=GetCharResidue(ResALL[II],AB); else rc='X';
                if(efptr) fprintf(stderr,"%c%d: (null is best) D=%d; nAln=%d; prob=%3f.\n",
                        rc,StartRes,D,nAln,prob);
                bestr=0;
        } else if(bestm < 0){	 // D and L too small...
                // if(D > 0) fprintf(stderr,"Error: X=%d; D = %d; L=%d; prob=%3f.\n",bestm,D,nAln,prob);
                bestm=0; prob=1.0;
                bestr=pos[1];
                // print_error("pvcalc input error");
        } else  bestr=pos[bestm];
#if 1	// Use for KeyAtm (-K=) option.
	if(hits[0]=='*'){ bestm++; } 
#endif
	bestr++; L=bestr; hits[L] = 0; free(pos); 
        return prob;
}

double	scm_typ::SurfaceClust(Int4 II,Int4 jj, char *hits,Int4 *String, Int4 &L, Int4 &bestm, Int4 &D)
{
        Int4    i,j,s,JJ,bestr;
        char    cJ,cI;
        double  prob;

	// h_type HG=Histogram("Pattern residue surface access",0,1,0.02);
	dh_type dH = dheap(num_resC+5,4);
	for(j=1; j <= num_resC; j++){
	  Int4 site=ResidueID(ResALL[j]);
	  if(site < resStart || site > resEnd) continue;
	  if(!MemberSet(site,resSet)) continue;
	  if(DistMtrx[j][II] == 0) insrtHeap(j,(keytyp)-2,dH);
	  else {
	    double dd=ResidueSurfaceAccessibility(ResALL[j]);
	    if(res_class[jj]=='M') insrtHeap(j,(keytyp)dd,dH); 
	    else insrtHeap(j,(keytyp)-dd,dH);
	  } // IncdHist(dd,HG); 
        } // PutHist(stderr,60,HG); NilHist(HG);
	for(s=0,i=0; !emptyHeap(dH); i++) {
		double dd = (double) minkeyHeap(dH);
		assert((j=delminHeap(dH)) != 0);
		if(MemberSet(j,SetX)) hits[s]='*'; else hits[s]=' ';
		if(efptr && hits[s]=='*'){
	           cJ=GetCharResidue(ResALL[j],AB); JJ=ResidueID(ResALL[j]);
		   fprintf(efptr,"%d: %c%d (%.3f)%c\n",i,cJ,JJ,dd,hits[s]);
		} s++; String[s]=j; 
	} Nildheap(dH); nAln=s-1;
	// fprintf(efptr,"nAln=%d; s=%d; CardSet=%d\n",nAln,s,CardSet(resSet));
	// fprintf(efptr,"%d: hits=%s\n",ii,hits);

        Int4       *pos; NEW(pos, nAln +6, Int4);
        for(i=1,s=0; i <= nAln; i++) { if(hits[i]=='*'){ s++; pos[s]=i; } }
        D=s;
        bestm=pvcalc(nAln,D,pos,prob,UseJeffreys);
        if(bestm == 0){
                char rc=GetCharResidue(ResALL[II],AB);
                if(efptr) fprintf(stderr,"%c: (null is best) D=%d; nAln=%d; prob=%3f.\n",
                        rc,D,nAln,prob);
                bestr=0;
        } else if(bestm < 0){
                // if(D > 0) fprintf(stderr,"Error: X=%d; D = %d; L=%d; prob=%3f.\n",bestm,D,nAln,prob);
                bestm=0; prob=1.0;
                bestr=pos[1];
                // print_error("pvcalc input error");
        } else  bestr=pos[bestm];
        free(pos); bestm++; bestr++; L=bestr; hits[L] = 0;
        return prob;
}

double  scm_typ::OptInterfaceClust(char *hits,Int4 *String, Int4 &L, Int4 &bestm, Int4 &D)
// Method 7: Ordered interface residue clustering.
{
	Int4    i,j,s,n,JJ,bestr,x;
	double	prob;

	dh_type dH = dheap(num_resC+5,4);
	for(j=1; j <= num_resC; j++){
	   Int4 site=ResidueID(ResALL[j]);
	   if(site < resStart || site > resEnd) continue;
	   if(!MemberSet(site,resSet)) continue;
	   double dd=Buried[j]; 
	   if(dd == 0.0){ 	// put discrinating residues at the end...
		if(MemberSet(j,SetX)) x=num_resC+j; else x=j;
		insrtHeap(j,(keytyp)x,dH);   // use positive numbers.
	   } else insrtHeap(j,(keytyp)-dd,dH); 
	}
	s=0; hits[s]=' '; String[s]=0; 
	for(n=i=0; !emptyHeap(dH); i++) {
		double dd = (double) minkeyHeap(dH);
		assert((j=delminHeap(dH)) != 0);
		if(MemberSet(j,SetX)) hits[s]='*'; else hits[s]=' ';
		if(efptr && hits[s]=='*'){
		   n++;
		   char cJ=GetCharResidue(ResALL[j],AB); JJ=ResidueID(ResALL[j]);
		   // fprintf(efptr,"%d: %c%d (%.3f)%c\n",i,cJ,JJ,dd,hits[s]);
		   fprintf(efptr,"%d(%d;%d): %c%d (%.3f)%c\n",i,j,n,cJ,JJ,dd,hits[s]);
		} s++; String[s]=j; 
	} Nildheap(dH); nAln=s-1;
	// fprintf(efptr,"nAln=%d; s=%d; CardSet=%d\n",nAln,s,CardSet(resSet));

        Int4       *pos; NEW(pos, nAln +6, Int4);
        for(i=s=0; i <= nAln; i++) { if(hits[i]=='*'){ s++; pos[s]=i+1; } }
        D=s;
        bestm=pvcalc(nAln,D,pos,prob,UseJeffreys); nAln--;
        if(bestm == 0){
                if(efptr) fprintf(stderr,"%... (null is best) D=%d; nAln=%d; prob=%3f.\n", D,nAln,prob);
                bestr=0;
        } else if(bestm < 0){
                // if(D > 0) fprintf(stderr,"Error: X=%d; D = %d; L=%d; prob=%3f.\n",bestm,D,nAln,prob);
                bestm=0; prob=1.0;
                bestr=pos[1];
                // print_error("pvcalc input error");
        } else  bestr=pos[bestm];
        free(pos); L=bestr; hits[L] = 0; 
	return prob;
}

