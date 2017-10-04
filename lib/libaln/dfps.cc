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

#include "asset.h" /***** asset search routines ******/

ast_typ	AssetSearch(ast_typ D)
{
	b_type	*IB=D->B,B;
	Int4	j;
	double	p;
	char	c;

   if(D==NULL) asset_error("fatal");
   if(D->P==NULL) InitAsset(D);
   if(FALSE && D->verbose) { PutEvalue(stdout,D->E); fprintf(stdout,"\n"); }
   if(D->pheap == NULL) D->pheap =
		MkPheap(D->d_max,D->d_min,D->k_max,D->hpsz,D->A);
   if(D->mode == 0 || D->s_min>= D->d_max){
	B = Block(BlockN(EBlocksBu(D->eblocks)));
	for(j=0; j < nAlpha(D->A); j++) {
		c=EBlocksRes(j,D->eblocks);
		fprintf(stderr,"%c",AlphaChar(c,D->A));
		if((p= D->tfreq[c]) < 1.0 && p > 0.0) {
		     D->cnts[1]++;
		     IntersectEBlockBu(0,c,B,D->eblocks);
		     CopyBlockL(IB[1],B);
		     setPattern(c,0,D->Q);
		     asset_dfps(D,2,1,CardBlockL(IB[1]),IB,p);
		     setVPattern(0,D->Q);
		}
	} NilBlock(B);
   }else {
        CopyBlockL(IB[0],EBlocksBu(D->eblocks));
        if(D->s_min < 2) asset_error("s_min must be at least 2");
	else asset_dfps12A(D,1,0,CardBlockL(IB[0]),IB,1.0);
   }
   return D;
}

Int4     EvaluateAsset(double *p, Int4 *idex, ast_typ D)
{
        double  p1;
        Int4     c,C=0,i;
        BooLean universe=TRUE;

    if(D->Q == NULL) asset_error("Q == NULL in EvaluateAsset( )");
    else {
        CopyBlockL(D->IB,EBlocksBu(D->eblocks));
        for(*p=1.0,i=0;i<D->k_max;i++) {
           if(NotEqUPattern(i,D->Q) && NotVarPattern(i,D->Q)) {
                *idex = i;
                universe=FALSE;
                ClearBlock(D->UB);
                for(p1=0.0,c=1;c<=nAlpha(D->A);c++){
                        if(MemPattern(i,c,D->Q)){
                                p1+=D->tfreq[c];
                                UnionEBlock(i,c,D->UB,D->eblocks);
                        }
                }
                C=IntersectBlockCF(D->UB,D->IB,D->IB);  /*function */
                (*p)= (*p)*p1;
           }
        }
    }
    if(universe) return CardBlock(EBlocksBu(D->eblocks));
    else return C;
}

void	asset_dfps(ast_typ D,Int4 n, Int4 idex, Int4 CB, b_type *B,double p1)
{
	register b_type	IIB=B[n],IB=B[n-1];
	b_type	EOB=D->EOB[n];
	double	prob,p,E,sd;
	Int4	i,j,c,C,N,end,*L,Ctotal;
	a_type	A=D->A; 

	sd = D->sd;
	end = D->k_max - (D->d_max- n);
	for(i=idex; i < end; i++){
	    CopyBlockL2(EOB, IB); /* */
	    N = EBlocksN_k(i+1,D->eblocks);
	    for(Ctotal=CB,j=0; j < nAlpha(A); j++) {
		c = EBlocksRes(j,D->eblocks);
		C=IntersectEBlockLCXOR(i,c,EOB,IIB,D->eblocks);
#if 0
		/* C=IntersectEBlockCF(i,c,IB,IIB,D->eblocks); /* OFF */
#endif 
		D->cnts[n]++;
		if(C >= D->c_min){
		   p=p1*D->tfreq[c];  
		   if(n >= D->d_min){
		     E=(p*(double)N);
		     if((double)C-E > sd*sqrt(E*(1.0-p))) {
		       if((prob=Log10CBP(C,N,p))!=ILLEGAL){
		       D->ncalc++;
		       if(D->H != NULL) IncdHist(MultEvalue(D->E,prob),D->H);
		       if(prob <= pmax0Evalue(D->E)){ /* top scores saved */
			 if(NumEBlocksPL(&L,IIB,D->eblocks) >= D->n_min){
			    setPattern(c,i,D->Q);
			    InsertPheap(D->Q,D->pheap,
				List2SegsEBlocks(L,D->eblocks),C,prob);
			    setVPattern(i,D->Q);
			 } free(L);
		       }
		     }
		   }
		 }
		 if(n < D->d_max) {
			setPattern(c,i,D->Q);
		   	asset_dfps(D,n+1,i+1,C,B,p);
			setVPattern(i,D->Q);
		  } 
		}
		if((Ctotal -= C) < D->c_min) break;
	   }
	}
}

void	asset_dfps12A(ast_typ D,Int4 n, Int4 idex, Int4 CB, b_type *B,double p1)
{
	register b_type	IIB=B[n],IB=B[n-1];
	b_type	EOB=D->EOB[n];
	double	prob,p,E,sd;
	Int4	i,j,c,C,N,end,*L;
	a_type	A=D->A; 
	Int4	Ctotal;

	sd = D->sd;
	if(n == 1) end = 1; 
	else if(n == 2 && D->s_min > 2) end = D->k_max -1;
	else end = D->k_max;
#if 0
	/** else end = MINIMUM(Int4, D->k_max, D->k_max - (D->s_min - n));/***/
#endif
	for(i=idex; i < end; i++){
	    CopyBlockL2(EOB, IB);
	    N = EBlocksN_k(i+1,D->eblocks);
	    for(Ctotal=CB,j=0; j < nAlpha(A); j++) {
		c = EBlocksRes(j,D->eblocks);
	        if(n==1) fprintf(stderr,"%c",AlphaChar(c,D->A));
		C=IntersectEBlockLCXOR(i,c,EOB,IIB,D->eblocks);
		D->cnts[n]++;
		if(C >= D->c_min){ 
		  p=p1*D->tfreq[c];  
		  if(n >= D->d_min) {
	           E=(p*(double)N);
		   if((double)C-E > sd*sqrt(E*(1.0-p)) &&
			(prob=Log10CBP(C,N,p))!=ILLEGAL){
		      D->ncalc++;
		      if(D->H != NULL) IncdHist(MultEvalue(D->E,prob),D->H);
                      if(prob <= pmax0Evalue(D->E)){ /* top scores saved */
			if(NumEBlocksPL(&L,IIB,D->eblocks) >= D->n_min){
                            setPattern(c,i,D->Q);
                            InsertPheap(D->Q,D->pheap,
				List2SegsEBlocks(L,D->eblocks),C,prob);
                      	    setVPattern(i,D->Q);
                        } free(L);
		      }
		   } 
		  }
		  if(n < D->d_max) {
			setPattern(c,i,D->Q);
			if(n < D->s_min) asset_dfps12A(D,n+1,i+1,C,B,p);
			else asset_dfps12T(D,n+1,i,C,B,p);
			setVPattern(i,D->Q);
		  } 
		}
		if((Ctotal -= C) < D->c_min) break;
	   }
	}
}

void	asset_dfps12T(ast_typ D,Int4 n,Int4 i,Int4 CB, b_type *B,double p1)
/*********************************************************************
	       +---------+	  0	(k_max = 11; i=7)
     ..........G...G..A...	offset	(w = 2*11-1=21; w-k_max = 10)
              +---------+	  1
             +---------+	  2
            +---------+		  3   (stop = MINIMUM(k_max-i,w-k_max) = 4) 

*********************************************************************/
{
	Int4	offset,stop;

	D->end = i+1; 
	asset_dfps12B(D,n,1,CB,B,p1);
	stop = (D->k_max-i);
	for(offset = 1; offset < stop; offset++){
		ShiftRPattern(D->Q); 
		ShiftREBLocks(D->eblocks);;
		D->end++; 
		asset_dfps12B(D,n,0,CB,B,p1);
	}
	for(offset = 1; offset < stop; offset++){
		ShiftLPattern(D->Q);
		ShiftLEBLocks(D->eblocks);;
		D->end--; 
	}
}

void	asset_dfps12B(ast_typ D,Int4 n,Int4 i,Int4 CB, b_type *B,double p1)
{
	b_type	IB,IIB,***b2,**b;
	b_type  EOB=D->EOB[n];
	double	prob,p,p2,E,*freq,sd;
	Int4	c,c2,C,N,j,k,c_min,end,*L;
	Int4	c_end=EndFlagEBlocks(D->eblocks);
	Int4	Ctotal, C2,*Card=D->Card[n];
	char	*rel;
	a_type	A=D->A; 

	IB = B[n-1]; IIB = B[n]; sd = D->sd;
	freq = D->tfreq; c_min = D->c_min;
	if(i==0) end = 1;
	else end = D->k_max;
#if 0
	/** else end = D->k_max - (D->d_max - n); /***/
#endif
	while(i<end) {
	  if(VarPattern(i,D->Q)) {
	   CopyBlockL2(EOB, IB);
	   k = MAXIMUM(Int4,i+1,D->end);
	   N = EBlocksN_k(k,D->eblocks);
	   for(Ctotal=CB,j=0; j < nAlpha(A); j++) {
	      c=EBlocksRes(j,D->eblocks);
	      if(PairedAlpha(c,A) || i >= D->end){
		Card[c]=C=IntersectEBlockLCXOR(i,c,EOB,IIB,D->eblocks);
		D->cnts[n]++;
		if(i >= D->end && C>=c_min) {
	           p2=freq[c]; p=p1*p2;  E=(p*(double)N);
	  	   if(n >= D->d_min && (double)C-E > sd*sqrt(E*(1.0-p)) &&
			(prob=Log10CBP(C,N,p))!=ILLEGAL){
		      D->ncalc++;
		      if(D->H != NULL) IncdHist(MultEvalue(D->E,prob),D->H);
                      if(prob <= pmax0Evalue(D->E)){ /* top scores saved */
			if(NumEBlocksPL(&L,IIB,D->eblocks) >= D->n_min){
                            setPattern(c,i,D->Q);
                            InsertPheap(D->Q,D->pheap,
				List2SegsEBlocks(L,D->eblocks),C,prob);
                      	    setVPattern(i,D->Q);
                        } free(L);
		      }
		   } 
		   if(n < D->d_max) {
			setPattern(c,i,D->Q);
			asset_dfps12B(D,n+1,i+1,C,B,p);
			setVPattern(i,D->Q);
		   } 
		}
	        /**** search for 2 residues *****/
	 	for(rel=RelatedEBlocks(c,D->eblocks);(c2=(*rel))!=c_end;rel++){
		   if(Card[c] && Card[c2] 	/* note: sets disjoint */
			&& (C2=Card[c]+Card[c2]) >= c_min){
	        	p2 = freq[c] + freq[c2];
			p=p1*p2;  E=(p*(double)N);
	  	        if(n >= D->d_min && ((double)C2-E)>sd*sqrt(E*(1.0-p)) 
				&& (prob=Log10CBP(C2,N,p))!=ILLEGAL){
		      	    D->ncalc++;
		      	    if(D->H != NULL) IncdHist(MultEvalue(D->E,prob),D->H);
                            if(prob <= pmax0Evalue(D->E)){
			      IntersectEBlockF2(i,c,c2,IB,IIB,D->eblocks);
			      if(NumEBlocksPL(&L,IIB,D->eblocks) >= D->n_min){
                  	          setPattern(c,i,D->Q);
			          AddPattern(c2,i,D->Q);
                            	  InsertPheap(D->Q,D->pheap,
					List2SegsEBlocks(L,D->eblocks),C2,prob);
                            	  setVPattern(i,D->Q);
			      } free(L);
			    }
			} 
			if(n < D->d_max) {
			    setPattern(c,i,D->Q);
			    AddPattern(c2,i,D->Q);
			    IntersectEBlockF2(i,c,c2,IB,IIB,D->eblocks);
			    D->cnts[n]++;
			    asset_dfps12B(D,n+1,i+1,C2,B,p);
			    setVPattern(i,D->Q);
			}
		   }  
		}
		if((Ctotal -= C) < 1) break;
	      }
	   }
	  } i++;
	}
}

