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

#include "eblocks.h"

#if 0
Int4     NonOverlapEBlocks(Int4 **L, Int4 minseq, Int4 length, b_type B,
	ebs_typ D, Int4 *N)
/* if n >= minseq; returns the number of non-overlapping segments */
/* else returns -1; note: length = LengthPattern(G) */
{
        s_type  S1,S2;
        Int4     C,II,I,i,n,*L0,o1,o2,os;

        L0=ListBlockL(B);
	/******/
	if((os=(D->w - D->k_max)-D->offset)!=0) OffsetIArray2(os,L0);
	/******/
        S2=NULL; II=-1;
        for(C=i=n=0; L0[i] != -1; i++) {
           S1 = D->segment[L0[i]];
           I=SegmentI(S1);
           o1 = SegmentStart(S1);
           if(I != II) n++;
           else if(S2 != NULL) {        /* S2 != NULL && I == II */
                while((abs(o1 - o2) < length)){
                        if(L0[i+1] == -1) { C--; break; }
                        i++; S1 = D->segment[L0[i]];
                        I=SegmentI(S1);
                        o1 = SegmentStart(S1);
                        if(I != II) { n++; break; }
                }
           }
           C++; S2 = S1; o2 = o1; II=I;
        }
	*L = L0;
	*N = n;
        return C;
}
#endif

Int4     NumEBlocksP(b_type B,ebs_typ D)
/* returns the number of sequence entities for the segments of block B */
{ Int4     n,*L; n= NumEBlocksPL(&L, B,D); free(L); return n; }

Int4     NumEBlocksPL(Int4 **List, b_type B, ebs_typ D)
/* returns the number of sequence entities for the segments of block B */
{
	Int4     II,I,i,n,*L,os;

	L=ListBlockL(B);/**/
	if((os=(D->w - D->k_max)-D->offset)!=0) OffsetIArray2(os,L);
	for(II = -1,i=n=0; L[i] != -1; i++) {
		assert(L[i] <= D->nseg);
		I=SegmentI(D->segment[L[i]]);
		if(I != II) {n++; II=I; }
	}
	*List = L;
	return n;
}

s_type  *List2SegsEBlocks(Int4 *L, ebs_typ D)
{
        Int4  i,C;
        s_type  *seg;

        for(C=0; L[C]!=-1; C++) ;
        NEW(seg,C+2,s_type);
        for(i=0; L[i]!=-1; i++){ seg[i] = D->segment[L[i]]; }
        seg[i] = NULL;
        return seg;
}

char    SegValXEBlocks(Int4 i, s_type S, ebs_typ D)
{
        Int4 x,I;
        e_type E;
        x = i + SegmentStart(S);
        I=SegmentI(S); E=SeqSetE(I,D->P);
        if(x > 0 && x <= (Int4) LenSeq(E)){ return XSeq(x,E); }
        else return UndefAlpha(D->A);
}

void    eblocks_error(char *s) { fprintf(stderr,"EBlocks: %s\n",s); exit(1); }

void	ShiftREBLocks(ebs_typ D)
{
	D->b--; D->b2--; D->offset++; 
	if(D->offset >= (D->w-D->k_max))
		eblocks_error("offset limit exceeded");
}

void	ShiftLEBLocks(ebs_typ D)
{
	D->b++; D->b2++; D->offset--;
	if(D->offset < 0) eblocks_error("offset limit exceeded");
}

ebs_typ MkEBlocks(ss_type P, Int4 k_max)
/**************** Allocate fully overlapping segments *************
sequence: MAVSTDEDSL			lenseq = 10	k_max = 5
      xxxxMAVSTDEDSLxxxxxx		width = 2 x k_max - 1 = 9
	  v   		       :...|.:...| start: (2-k_max..lenseq)
      xxxxMAVST 		xxxxMAVST   -3
       xxxMAVSTD		xxxMAVSTD   -2
        xxMAVSTDE		xxMAVSTDE   -1
         xMAVSTDED		xMAVSTDED    0
          MAVSTDEDS		MAVSTDEDS    1
           AVSTDEDSL		AVSTDEDSL    2
            VSTDEDSLx		VSTDEDSLx    3
             STDEDSLxx		STDEDSLxx    4
              TDEDSLxxx		TDEDSLxxx    5
               DEDSLxxxx	DEDSLxxxx    6
                EDSLxxxxx	EDSLxxxxx    7
                 DSLxxxxxx	DSLxxxxxx    8
                  SLxxxxxxx	SLxxxxxxx    9
                   Lxxxxxxxx	Lxxxxxxxx   10
		       ^
   Number overlaping seg = 14 = 5 -1 + len = (5-1) + 10.
   
******************************************************************/
{
	e_type	E;
	Int4	i,j,I,k,c,sdex,start;
	double	*tfreq;
	s_type	S;
	a_type	A=SeqSetA(P);
	dh_type H;
	ebs_typ D;

	NEW(D,1,eblocks_type);
	D->P = P;
	D->A = A;
	D->mode = 0;
	D->k_max = k_max; 
	D->w = 2*D->k_max-1;
	D->rel = NULL;
	D->res = NULL;
	D->offset = 0;
	D->eb2 = D->b2 = NULL;
	
	/** order residues **/
	H = dheap((nAlpha(A)+1),3);
	tfreq = tFreqSeqSet(P);
	for(c=1; c <= nAlpha(A); c++){
		insrtHeap(c,-(keytyp) tfreq[c],H);
	}
	NEW(D->res, nAlpha(A)+1,char);
	for(i=0;(c=delminHeap(H)) != NULL; i++){D->res[i]=(char)(c);}
	Nildheap(H);

	/**** fully overlapping segments ***/
	for(D->nseg=0,I=1 ;I <= NSeqsSeqSet(D->P);I++) {
		E=SeqSetE(I,D->P);
		D->nseg += k_max - 1 + LenSeq(E);  
	}
	MEW(D->segment,D->nseg+2,s_type);
	// NEW(D->segment,D->nseg+k_max+2,s_type); // doesn't help...
	for(j=0,I=1;I <= NSeqsSeqSet(D->P);I++) {
		E=SeqSetE(I,D->P);
		for(start=2-k_max; start <= (Int4)LenSeq(E); start++) {
		    D->segment[j]= Segment(I,start);
		    j++;  
		}
	}
	D->Bu=Block(D->nseg); D->dummy = Block(D->nseg); 
	FillBlock(D->Bu);
	/* D->N[k] == # of segments of length k in population */
	if(D->N==NULL) NEW(D->N,k_max+1,Int4); 
	for(k=1; k<=k_max; k++){
	    for(D->N[k]=0, I=1 ;I <= NSeqsSeqSet(D->P);I++) {
		E=SeqSetE(I,D->P);
		D->N[k] += MAXIMUM(Int4,0, LenSeq(E) - k + 1);  
	    }
	}
	/* create k_max x nAlpha elementary sets representing the segments */
        NEWP(D->eb,(D->w+1),b_type); D->b = D->eb + (D->w - k_max);
        for(i=0;i<D->w;i++) {
        	NEW(D->eb[i],(nAlpha(A)+2),b_type); /* eb[i][0] for dummy */
	}
        for(i=(D->w - k_max); i < D->w; i++) {
		D->eb[i][0] = D->dummy; /* for dummy residues */
        	for(c=1;c<=nAlpha(A);c++) D->eb[i][c] = Block(D->nseg);
	}
	for(sdex=0; sdex < D->nseg; sdex++) {
	    S=D->segment[sdex];
            for(i=(D->w - k_max); i < D->w; i++) {
	      c=SegValXEBlocks(i,S,D);
	      if(c != UndefAlpha(A)) AddBlock(sdex, D->eb[i][c]);
	    }
	}                               
	D->eb2 = NULL; 
	D->b2 = NULL; 
	return D;
}

void	NilEBlocks(ebs_typ D)
{
	Int4	i,c;

	if(D->b2 != NULL) NilEBlocks2(D);
	free(D->segment);
	free(D->N);
	free(D->res);
	if(D->Bu != NULL) NilBlock(D->Bu);
	if(D->dummy != NULL) NilBlock(D->dummy);
	if(D->eb != NULL) {
            for(i=0;i<D->w;i++) {
        	for(c=1;c<=nAlpha(D->A);c++) 
			if(D->eb[i][c]!=NULL) NilBlock(D->eb[i][c]);
		free(D->eb[i]);
	    }
	    free(D->eb); D->eb = NULL; D->b = NULL;
	}
	free(D);
}

void	AddEBlocks2(ebs_typ D, Int4 mode)
/* make two residue eblocks */
{
	b_type	***b2; 
	a_type	A=D->A; 
	s_type	S;
	Int4	i,j,k,sdex;
	char	c,c2,**rel;

	D->mode = mode;
/*******/
        for(i=0; i < (D->w - D->k_max); i++) {
		D->eb[i][0] = D->dummy; 
        	for(c=1;c<=nAlpha(A);c++) D->eb[i][c] = Block(D->nseg);
	}
	for(sdex=0; sdex < D->nseg; sdex++) {
	    S=D->segment[sdex];
            for(i=0; i < (D->w - D->k_max); i++) {
	      c=SegValXEBlocks(i,S,D);
	      if(c != UndefAlpha(A)) AddBlock(sdex, D->eb[i][c]);
	    }
	}                               
/*******/
	NEWP(rel,nAlpha(A)+2,char);
	for(j=0; j < nAlpha(A); j++) {
		c=D->res[j];
		NEW(rel[c],nAlpha(A)+2,char);
		for(i=0,k=j-1; k >= 0; k--) {
			c2=D->res[k];
			if(valAlphaP(c,c2,A) <= mode){ rel[c][i++]=c2; }
		}
		rel[c][i] = -1;
	}
	D->rel = rel;
	NEWPP(b2,D->w+1,b_type);
	for(i=0;i<D->w;i++){
		NEWP(b2[i],nAlpha(A)+2,b_type);
		for(j=0; j < nAlpha(A); j++) {
			c=D->res[j];
			NEW(b2[i][c],nAlpha(A)+2,b_type);
			b2[i][c][c]=D->eb[i][c]; 
			for(k=j-1; k >= 0; k--) {
				c2=D->res[k];
		      		if(valAlphaP(c,c2,A) <= mode){
					b2[i][c][c2]=Block(BlockN(D->Bu)); 
			    		UnionBlock3(D->eb[i][c],
						    D->eb[i][c2],b2[i][c][c2]);
				} else b2[i][c][c2]=NULL;
			}
		}
	}
	D->eb2=b2; D->b2 = b2 + (D->w - D->k_max);
}

void	NilEBlocks2(ebs_typ D)
{
	b_type	***b2=D->eb2; 
	a_type	A=D->A; 
	Int4	i,j,k;
	char	c,c2;

	/*******/
        for(i=0; i < (D->w - D->k_max); i++) {
		D->eb[i][0] = NULL;	/* for dummy residues */
        	for(c=1;c<=nAlpha(A);c++) {
			NilBlock(D->eb[i][c]); D->eb[i][c] = NULL;
		}
	}
	/*******/
	for(i=0;i<D->w;i++){
		for(j=0; j < nAlpha(A); j++) {
	 	  c=D->res[j];
		  for(k=j-1; k >= 0; k--) {
			c2=D->res[k];
			if( b2[i][c][c2]!=NULL) {
				NilBlock(b2[i][c][c2]);
				b2[i][c][c2] = NULL;
			}
		  }
		  free(b2[i][c]);
		}
		free(b2[i]);
	}
	free(b2);
	D->eb2 = NULL; D->b2 = NULL;
	for(c=1;c<=nAlpha(A);c++) free(D->rel[c]);
	free(D->rel);
	D->rel = NULL;
}

