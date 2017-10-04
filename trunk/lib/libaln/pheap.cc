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

#include "pheap.h"

ph_type	MkPheap(Int4 d_max, Int4 d_min, Int4 k_max, Int4 hpsz, a_type A)
{
	ph_type H;
	NEW(H,1,pheap_type); 
	NEW(H->pattern,hpsz+1,ptn_typ); 
	NEW(H->C, hpsz+1, Int4);
	NEWP(H->S, hpsz+1, s_type);
	H->mheap = Mheap(hpsz,4);
	H->A = A;
	H->d_max = d_max; H->d_min = d_min; H->k_max = k_max;
	return H;
}

ph_type	NilPheap(ph_type H)
{
	Int4 i;

	if(H==NULL) return (ph_type) NULL;
	for(i=0;i<=SizeMheap(H->mheap); i++){
		if(H->pattern[i] != NULL){
			NilPattern(H->pattern[i]);
		}
		if(H->S[i] != NULL) NilSegmentList(H->S[i]);
	}
	free(H->S); free(H->pattern); free(H->C);
	NilMheap(H->mheap);
	free(H); 
	return (ph_type) NULL;
}

ptn_typ	DelMinPheap(s_type **S, ph_type H)
{
	Int4 i;
	ptn_typ G;
	
	if((i=DelMinMheap(H->mheap))!=NULL){
		G=H->pattern[i];
		H->pattern[i]=NULL;
		if(S==NULL) free(H->S[i]); else *S = H->S[i]; 
		H->S[i] = NULL;
		return G;
	} else return (ptn_typ) NULL;
}

ptn_typ	DelMaxPheap(s_type **S, ph_type H)
{
	Int4 i;
	ptn_typ G;
	
	if((i=DelMaxMheap(H->mheap))!=NULL){
		G=H->pattern[i];
		H->pattern[i]=NULL;
		if(S==NULL) free(H->S[i]); else *S = H->S[i]; 
		H->S[i] = NULL;
		return G;
	} else return (ptn_typ) NULL;
}

/********************* storage of significant motifs ********************/
Int4	InsertPheap(ptn_typ G, ph_type H, s_type *S, Int4 C, double prob)
{	
	Int4 i;

	// if(QuickPurgePheap(G, C, H)) { NilSegmentList(S); return NULL; }
	i = InsertMheap(prob,H->mheap);
	if(i != NULL){
		if(H->pattern[i] != NULL) NilPattern(H->pattern[i]);
		H->pattern[i]= CopyPattern(G);
		if(H->S[i] != NULL) NilSegmentList(H->S[i]);
		H->S[i] = S; H->C[i] = C;
	} else NilSegmentList(S);
	return i;
}

BooLean	QuickPurgePheap(ptn_typ Q, Int4 C, ph_type H)
/* eliminate safe submotifs on the fly. */
{
	Int4	j,end;
	Int4	*card = H->C;
	s_type	**S = H->S;
	ptn_typ	*pattern=H->pattern;

	end= SizeMheap(H->mheap);
	for(j=1;j <= end;j++){
	   if(pattern[j] != NULL && C == card[j]){
		if(SubPattern(pattern[j],Q)){   	/* Qx > Q */
			NilPattern(pattern[j]); pattern[j]=NULL;
			NilSegmentList(S[j]); S[j]=NULL;
			RmMheap(j,H->mheap);
		} else if(SubPattern(Q,pattern[j])){	/* Q > Qx */
			return TRUE;
		}
	   }
	}
	return FALSE;
}

Int4	PurgePheap(ph_type H)
/* eliminates all supermotifs - returns number of patterns removed */
{
	Int4	n,i,j,end;
	Int4	*card = H->C;
	s_type	**S = H->S;
	ptn_typ	G1,G0,*array=H->pattern;

	end= SizeMheap(H->mheap);
	for(n=0,i=1;i <= end;i++){
	  if((G1=array[i]) !=NULL) {
	    for(j=i+1;j <= end;j++){
		if((G0=array[j]) != NULL && card[i] == card[j]){
			if(SubPattern(G0,G1)){   		/* G0 > G1 */
				/*******
				PutPattern(stderr,G0,H->A);
				fprintf(stderr," > \n\t");
				PutPattern(stderr,G1,H->A);
				fprintf(stderr,"\n\n");
				/*******/
				NilPattern(G0); array[j]=NULL;
				NilSegmentList(S[j]); S[j]=NULL;
				RmMheap(j,H->mheap); n++; 
			} else if(SubPattern(G1,G0)){	/* G1 > G0 */
				NilPattern(G1); array[i]=NULL;
				NilSegmentList(S[i]); S[i]=NULL;
				RmMheap(i,H->mheap); n++; 
				break;
			}
		} 
	    }
	  }
	}
	return n;
}

Int4	PutPheap(FILE *fptr, ph_type H)
{
	Int4	i,end=SizeMheap(H->mheap);
	ptn_typ	G;
	UInt4 total=0,n,j,I;

	for(i=1;i <= end;i++){
	   if((G=H->pattern[i])!=NULL) {
		PutPattern(fptr,G,H->A); 
		BubbleSortSegments(H->S[i]);
		for(n=0,I=0,j=1; j<=H->C[i]; j++){
			if(I!=SegmentI(H->S[i][j])){ n++; I=SegmentI(H->S[i][j]); }
		}
		fprintf(fptr," (%d hits in %d seqs.)\n",H->C[i],n);
		total+=H->C[i];
	   }
	}
	fprintf(fptr,"\t%d total hits\n",total);
	return end;
}

ptn_typ	*CombinePheap(ptn_typ **motif, s_type ***slist, Int4 *ngroup, 
	Int4 *npat, Int4 N, ph_type H)
{
	Int4	N1,N2,n,x,g,i,j,k,k1,k2;
	Int4	os,group,numG,ngrps,C,sC,*Card;
	s_type	*L,*LC,*L0,*L1,**SGL,**SList,**xlist;
	ptn_typ	G1,G2,*GA,sG,*SGA,*xmotif;
	double	p,K;

  // if(H->d_max < 4) print_error("depth must be greater than 3");
  if(H->d_max < 3) print_error("depth must be greater than 2");
  if(!EmptyPheap(H)) {
	NEW(GA,ItemsInPheap(H)+2,ptn_typ);
	NEW(Card,ItemsInPheap(H)+2,Int4);
	NEWP(SList,ItemsInPheap(H)+2,s_type);
	for(numG=C=0; !EmptyPheap(H); numG++) {
		Card[numG] = MinCardPheap(H);
		C += Card[numG];
		G1=DelMinPheap(&L,H);
		SList[numG] = L; GA[numG] = G1;
	}
	NEW(LC,C+2,s_type); NEW(L0,C+2,s_type); NEW(L1,C+2,s_type); 
	NEWP(SGL,numG,s_type); MEW(SGA,numG,ptn_typ); 
	for(i=0;i<numG;i++){ SGA[i]=NULL; SGL[i]=NULL; }
	fprintf(stderr,"Combining primary patterns... ");
/** TEST **/
	ngrps = CombineOSPHeap(GA,SGA,SGL,numG,N, SList, C); 
/** TEST **
	if(H->d_min < H->d_max){
		ngrps = Combine1PHeap(GA,SGA,SGL,numG,N,SList,C,H); 
	} else ngrps = CombineOSPHeap(GA,SGA,SGL,numG,N, SList, C); 
/** TEST **/
	fprintf(stderr,"(%d groups)\n",ngrps);
	group = ngrps;
      /*** USE HYPERGEOMETRIC PROBABILITY TO COMBINE SUPERPATTERNS ***/
     if(group > 1){
/*******
	fprintf(stderr,"Estimate number of comparisons... ");
	for(K=0.0,j=0; group > 1 && j < ngrps; j++){ 
	   if(SGA[j] !=NULL) {
	      k = kPattern(SGA[j]);
	      for(i=j+1; i < ngrps; i++){	  
		if(SGA[i] != NULL){ 
			K+= (double)(k + (Int4) kPattern(SGA[i])); 
		}
	      }
	   }
	} 
	K = K/2.0;
	fprintf(stderr,"(%g comparisons)\n",K);
/*** DEBUG ***/ K = 1000; /****/
	fprintf(stderr,"Combine using hypergeometric probability... ");/***/
	for(j=0; group > 1 && j < ngrps; j++){
	 if((G1=SGA[j]) !=NULL){
	   k1 = LengthPattern(G1); k = j; 
	   for(i=0; i < ngrps; i++){	  
	      if((G2=SGA[i]) != NULL && G1 != G2){
	        k2 = LengthPattern(G2);
		if(IntersectOSegments(&x,3,&os,k1,k2,SGL[k],SGL[i])==1 && x>1){
		     if(os > 0) {
		   	N1 = CopySegments(L0,SGL[k]);
		        n = OffsetSegments(-os,SGL[i],L1);
		     } else {
		   	n = CopySegments(L1,SGL[i]);
			N1 = OffsetSegments(os,SGL[k],L0);
		     }
		     N2 = N - N1; 
		     if((p=CumHypGeomProb(N1,N2,n,x))<0.001/K || x==N1 || x==n){
		        if(os > 0) { sG = MergePatterns(G2,G1,-os); } 
			else { sG = MergePatterns(G1,G2,os); }
			sC = UnionSegments(LC,L0,L1);
			/******
		   	fprintf(stderr,"k2=%d; k1=%d; os=%d; x=%d\n",k2,k1,os,x);
	   		PutPattern(stderr,G1,H->A); 
			fprintf(stderr,"\n\t+ ");
	   		PutPattern(stderr,G2,H->A); 
			fprintf(stderr,"\n\t\t  = ");
	   		PutPattern(stderr,sG,H->A); 
			fprintf(stderr,"\n CumHypGeom(%d,%d,%d,%d) = %g\n",
				N1,N2,n,x,p);
			/*** DEBUG ***/
			NilSegmentList(SGL[k]); NilSegmentList(SGL[i]); 
			SGL[k] = NULL; SGL[i] = NULL; 
			NilPattern(G1); NilPattern(G2);
			SGA[k] = NULL; SGA[i] = NULL; 
			if(i > k) i = k;
			NEW(SGL[i],sC+2,s_type); CopySegments(SGL[i],LC);
			SGA[i]=G1=sG; k=i; 
	   		k1 = LengthPattern(G1);
			group--; j=(-1); i=ngrps; break;    /** start over **/
		     }
		   }
		}
	   }
	 }
	}
	fprintf(stderr,"(%d groups)\n",group); /*****/
     } /*** OUTPUT PATTERNS AND SUPERPATTERNS ***/
     for(group=i=0; i < ngrps; i++) if(SGA[i] != NULL) group++;
     NEW(xmotif,group+2,ptn_typ); NEWP(xlist,group+2,s_type);
     for(group=i=0; i < ngrps; i++){
	if((sG=SGA[i]) != NULL){
	   	xmotif[group] = SGA[i]; SGA[i] = NULL;
		xlist[group] = SGL[i];  SGL[i] = NULL;
		group++;
	}
     }
     for(g=0;g<numG;g++) { NilSegmentList(SList[g]); }
     free(SGA); free(SList); free(Card);
     free(LC); free(L0); free(L1); free(SGL); 
     *motif = xmotif; *slist = xlist; *ngroup = group; *npat = numG;
     return GA;
  } else return NULL;
}

Int4	CombineOSPHeap(ptn_typ *GA, ptn_typ *SGA, s_type **SGL, Int4 ngrps, 
	Int4 N, s_type **List, Int4 TotCard)
{
	Int4	add,out,*idex,*odex,*temp,group,os,i,sC;
	Int4	n,N1,N2,x,k1,k2;
	double	p;
	s_type	*L,*L0,*L1,*sL,*L2;
	ptn_typ	G,G1,G2,sG,fsG;

	if(ngrps == 1) { 
		SGA[0] = CopyPattern(GA[0]);
		for(sC=0; List[0][sC]!=NULL; sC++);
	 	NEW(SGL[0],sC+2,s_type); CopySegments(SGL[0],List[0]);
		return 1; 
	}
	MEW(idex,ngrps+1,Int4); MEW(odex,ngrps+1,Int4);
	NEW(sL,TotCard+2,s_type); NEW(L0,TotCard+2,s_type); 
	NEW(L1,TotCard+2,s_type); NEW(L2,TotCard+2,s_type); 
	for(i=0;i<ngrps;i++){ SGA[i]=NULL; idex[i] = i; } idex[i] = -1;
	for(group=0; (*idex) != -1; group++){
	 G1=GA[(*idex)];
	 k1 = LengthPattern(G1);
	 sG=CopyPattern(G1); sC = CopySegments(sL,List[(*idex)]); 
         if(idex[1]!=(-1)){
	      for(i=1,add=1; idex[i]!=(-1) && add!=0; i=0){
	       for(out=add=0; idex[i]!=(-1); i++){ 
	   	G2=GA[idex[i]]; 
	        k2 = LengthPattern(G2);
		if(IntersectOSegments(&x,3,&os,k1,k2,sL,List[idex[i]]) == 1 && x>1){
		     if(os > 0) {
		   	N1 = CopySegments(L0,sL);
		        n = OffsetSegments(-os,List[idex[i]],L2);
		     } else {
		   	n = CopySegments(L2,List[idex[i]]);
			N1 = OffsetSegments(os,sL,L0);
		     }
		     N2 = N - N1; 
		     if((p=CumHypGeomProb(N1,N2,n,x)) < 0.000001 || x==N1 || x==n){
		        if(os > 0) { fsG = MergePatterns(G2,sG,-os); } 
			else { fsG = MergePatterns(sG,G2,os); }
// fprintf(stderr,"\n CumHypGeom(%d,%d,%d,%d) = %g\n", N1,N2,n,x,p);
			sC = UnionSegments(L1,L0,L2);
			L = sL, sL = L1;  L1 = L;
			NilPattern(sG); sG = fsG; add++;
	 		k1 = LengthPattern(sG);
		   } else { odex[out++]=idex[i]; }
		} else odex[out++]=idex[i]; 
	       }
	       odex[out]=(-1);
	       temp=idex; idex=odex; odex=temp;
	      }
	 } else idex[0]=(-1);
	 SGA[group] = sG;
	 NEW(SGL[group],sC+2,s_type); CopySegments(SGL[group],sL);
	}
	free(odex); free(idex); 
	free(L0); free(L1); free(L2); free(sL); 
	return group;
}

Int4	Combine1PHeap(ptn_typ *GA, ptn_typ *SGA, s_type **SGL, Int4 ngrps, 
	Int4 N, s_type **List, Int4 TotCard, ph_type H)
{
	Int4	add,out,*idex,*odex,*temp,group,os,i,sC;
	Int4	n,N1,N2,x;
	double	p;
	s_type	*L,*L0,*L1,*sL,*L2;
	ptn_typ	G1,G2,sG,fsG;
	a_type	A=H->A;
	Int4	z,k1,k2,o;

	if(ngrps == 1) { 
		SGA[0] = CopyPattern(GA[0]);
		for(sC=0; List[0][sC]!=NULL; sC++);
	 	NEW(SGL[0],sC+2,s_type); CopySegments(SGL[0],List[0]);
		return 1; 
	}
	MEW(idex,ngrps+1,Int4); MEW(odex,ngrps+1,Int4);
	NEW(sL,TotCard+2,s_type); NEW(L0,TotCard+2,s_type); 
	NEW(L1,TotCard+2,s_type); NEW(L2,TotCard+2,s_type); 
	for(i=0;i<ngrps;i++){ SGA[i]=NULL; idex[i] = i; } idex[i] = -1;
	for(group=0; (*idex) != -1; group++){
	 G1=GA[(*idex)];
	 sG=CopyPattern(G1); sC = CopySegments(sL,List[(*idex)]); 
	 k1 = LengthPattern(sG);
         if(idex[1]!=(-1)){
	      for(i=1,add=1; idex[i]!=(-1) && add!=0; i=0){
	       for(out=add=0; idex[i]!=(-1); i++){ 
	   	G2=GA[idex[i]]; 
	        k2 = LengthPattern(G2);
		if(IntersectOSegments(&z,3,&o,k1,k2,sL,List[idex[i]]) == 1 
		  && z>1 &&
		  (fsG=CombinePatterns(sG,G2,&os,H->d_min,0,100))!=NULL){
		   if(os > 0) {
		   	N1 = CopySegments(L0,sL);
			n = OffsetSegments(-os,List[idex[i]],L2);
		   } else {
		   	n = CopySegments(L2,List[idex[i]]);
			N1 = OffsetSegments(os,sL,L0);
		   }
		   N2 = N - N1; x = IntersectSegments(L1,L0,L2);
		   if(x > 1 && (p=CumHypGeomProb(N1,N2,n,x)) <= 0.000001){
			sC = UnionSegments(L1,L0,L2);
			L = sL, sL = L1;  L1 = L;
			NilPattern(sG); sG = fsG; add++;
	 		k1 = LengthPattern(sG);
		   } else { odex[out++]=idex[i]; NilPattern(fsG); }
		} else odex[out++]=idex[i];
	       }
	       odex[out]=(-1);
	       temp=idex; idex=odex; odex=temp;
	      }
	 } else idex[0]=(-1);
	 SGA[group] = sG;
	 NEW(SGL[group],sC+2,s_type); CopySegments(SGL[group],sL);
	}
	free(odex); free(idex); 
	free(L0); free(L1); free(L2); free(sL); 
	return group;
}

