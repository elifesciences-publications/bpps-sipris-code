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

#include "pattern.h"

Int4	LengthPattern(ptn_typ G)
{ Int4	i,l=0; for(i=0;i<(Int4)G->k;i++) if(G->m[i]!=0) l=i+1; return l; } 

ptn_typ	Pattern(Int4 k,Int4 nlet)
/* Create a pattern of length k for an alphabet of n letters */
{
	ptn_typ G;
	if(k > MAX_PATTERN_LENGTH) pattern_error("pattern too long");
	NEW(G,1,pattern_type);
	NEW(G->m,k,sst_typ); 
	G->k = (unsigned char) k; G->nlet = (unsigned char) nlet; G->i=0;
	return G;
}

ptn_typ	CopyPattern(ptn_typ H)		/* Create a copy of the pattern H */
{
	Int4	i,k=(Int4)H->k; ptn_typ G;

	if(H == NULL) pattern_error("Attempt to copy null pattern");
	NEW(G,1,pattern_type);
	NEW(G->m,k,sst_typ);
	for(i=0;i<k;i++) G->m[i] = H->m[i];
	G->i = H->i; G->k = H->k; G->nlet = H->nlet; 
	return G;
}

Int4	DepthPattern(ptn_typ G)
/* Return the search depth needed to for the pattern G */
/* WARNING: Assumes that "." == 0 */
{ Int4	d,i,k=G->k; for(d=i=0;i<k;i++) if(G->m[i]) d++; return d; }

void	ShiftLPattern(ptn_typ P)
/* shift the pattern 1 space to the left: ...PPY. -> ..PPY.. */
{
	sst_typ *m; 
	Int4	i,k=(Int4)P->k;
	for(m=P->m,i=0; i<k-1;i++){ m[i] = m[i+1]; }
	m[k-1] = 0;
	if((Int4) P->i > 0) P->i--;
}

void	ShiftRPattern(ptn_typ P)
/* shift the pattern 1 space to the right:  .PPY... -> ..PPY..  */
{
	sst_typ *m; 
	Int4	i,k=(Int4)P->k;
	for(m=P->m,i=k-1; i>0;i--){ m[i] = m[i-1]; }
	m[0] = 0;
	if((Int4)P->i < k-2) P->i++;
}

ptn_typ	NilPattern(ptn_typ G)			/* Destroy the pattern G */
{ if(G!=NULL) {free(G->m); free(G);} return (ptn_typ) NULL; }

ptn_typ	PutPattern(FILE *fptr,ptn_typ G, a_type A)
/* Put the elements of each set for the 1-pattern pattern G for alphabet A */
{
	Int4	i,j,end; 
	UInt4	s,u=UPattern(G);
	char 	s2[33],ptr;
	char	c;
	/*char 	*str,pnt[]="<-  ",blk[]="    "; /*******/

	if(G == NULL) { fprintf(stderr,"\tillegal\n"); return NULL; }
	for(end=i=0; i < (Int4) G->k; i++){
	    if((s=G->m[i])!=u && s!=0) end=i;
	}
	for(i=0; i <= end; i++){
	    if((s=G->m[i])==u || s==0) fprintf(fptr,".");
	    else {
		for(ptr=0,c=1; c <= (Int4)nAlpha(A);c++){
			if(MemPattern(i,c,G)) s2[ptr++] = c;
		}
		if(ptr > 1){
			fprintf(fptr,"[");
			for(j=0; j < ptr; j++)
				fprintf(fptr,"%c",AlphaChar(s2[j],A));
			fprintf(fptr,"]");
		} else if(ptr ==1  ) fprintf(fptr,"%c",AlphaChar(s2[0],A));
		else pattern_error("Illegal pattern configuration");
	    }
	}
	/*******
        for(s2[32]=0,i=0;i<32;i++)
           {if((j=31-i)<=nAlpha(A)) s2[i]=AlphaChar(j,A);else s2[i]=' ';}
        fprintf(fptr,"POS     %s  CARD\n",s2);
	for(j=0;(Int4)j<(Int4)G->k;j++) {
	   if(j!=G->i) str=blk; else str=pnt; 
	   for(i=0;i<32;i++)
		if(G->m[j] & (1<<i)) s2[31-i]='1';
		else s2[31-i]='0';
	   for(k=1,n=0;(Int4)k<=(Int4)G->nlet;k++) if(MemSset(k,G->m[j])) n++;
	   fprintf(fptr,"%2d%s '%s'  %d\n", j,str,s2,n);
	}
	/*******/
	return G;
}

/*****
    from os = minid - k2 to k1 - minid ... 
	            v
	      D....NPPY			os = s1 - k2 = 7 - 13 = -6
	I..NPPYV....G	 
	      ^
    to
	D....NPPY			os = k1 - s2 = 9 - 3 = 6
	      I..NPPYV....G	 

    case 1: offset >= 0; 
	k1=9; k2=13; j=abs(*os).

	start1 = *os; start2=0 

	D....NPPY			for i= 0..k1-1	 -> m[i] = m1[i]
	  I..NPPYV....G	 *os = +2	for i= j..k2+j-1 -> m[i] |= m2[i-j]
	---------------
	D.I..NPPYV....G   k = max(k2+j,k1) = 13 + 2 = 15

    case 2: offset < 0.

	start2 = abs(*os); start1=0 

	  I..NPPYV....G	 		for i= j..k1+j-1 -> m[i] = m1[i-j]
	D....NPPY    	 *os = -2	for i= 0..k2-1 	 -> m[i] |= m2[i]
	---------------
	D.I..NPPYV....G   k = max(k1+j,k2) = 13 + 2 = 15

******/
ptn_typ	CombinePatterns(ptn_typ G1,ptn_typ G2,Int4 *os, Int4 minid, Int4 max,
	Int4 wt)
/* if superpattern G1 can be combined with pattern G2 then 
   create and return pattern G which is G1 and G2 */
{
	/*** static Int4 i=0,z=0,p=0,w=2;/****/
	Int4	score;
	ptn_typ	G=NULL;

	score = AlignScorePatterns(G1,G2, os, minid, max, wt);
	/*** if(score > 0) return MergePatterns(G1, G2, *os); /*** OLD ***/
/***
	i++; if(score > w*wt) p++;
	else z++;
	if(i%500 == 499) fprintf(stderr,"score > %d? %d: %d\n ",w,p,z); /****/
	if(score >= wt*minid) G = MergePatterns(G1, G2, *os);
	return G;
}

Int4     AlignScorePatterns(ptn_typ G1,ptn_typ G2, Int4 *os, Int4 minid, 
        Int4 max, Int4 wt)
{
	sst_typ	*m1,*m2;
	Int4	i,j,k1,k2,conflict,s1,s2;
	Int4	nident,start1,start2,score,tmp,offset,pident;
	Int4	spid, New,snew;

	/*** SHORTEN PATTERNS AS MUCH AS POSSIBLE **/
	for(m1=G1->m,k1=(Int4)G1->k; k1 > 0 && m1[k1-1] == 0; ) k1--;
	for(m2=G2->m,k2=(Int4)G2->k; k2 > 0 && m2[k2-1] == 0; ) k2--;
	for(s1=i=0; i < minid && s1 < k1; s1++) if(m1[s1]!=0)i++;
	for(s2=i=0; i < minid && s2 < k2; s2++) if(m2[s2]!=0)i++;
	score = 0; spid=0; snew=0;
	for(offset=s1-k2; (offset)<=(k1-s2); offset++) {
	    if(offset >= 0) { start2 = 0; start1 = offset; }
	    else { start2 = abs(offset); start1 = 0; }
	    conflict=0; pident=0; New=0;
	    for(nident=0,i=start1,j=start2; j<k2 && i<k1; i++,j++) {
		if(m1[i] != 0){
		   if(m2[j] != 0){
			if(m1[i] == m2[j]) nident++;
			else if(m1[i] & m2[j])pident++;
			else conflict++;
		   } 
		} if(m2[j] != 0) New++;
	    }
	    if(conflict <= max){
		tmp = (wt*nident + pident - wt*conflict);
		/** if(tmp >= wt*minid){ ... } /**OLD**/
		if(score < tmp){
		   snew = New; spid = pident;
		   score = tmp; *os = offset;
		}else if(score == tmp){
		   if(spid > pident){	/* prefer full ident */
			snew = New; spid = pident; *os = offset;
		   } else if(spid == pident){
		   	if(New < snew){	/* & fewer New */
				snew = New; *os = offset;
			} else if(New == snew){
		     	   if(abs(*os) > abs(offset)){ /* & short offsets */
				*os = offset;
			   } /* if |os| < |offset| ... */
			} /* if New > snew  ... */
		   } /* if spid < pident then don't change */
		}
	    }
	}
	return score;
}

ptn_typ	MergePatterns(ptn_typ G1, ptn_typ G2, Int4 os)
/* set G = G1 U G2 */
{
	sst_typ		*m1,*m2,*m;
	Int4		i,j,k1,k2,k,end;
	ptn_typ		G;

	for(m1=G1->m,k1=(Int4)G1->k; k1 > 0 && m1[k1-1] == 0; ) k1--;
	for(m2=G2->m,k2=(Int4)G2->k; k2 > 0 && m2[k2-1] == 0; ) k2--;
	j=abs(os);
	if(os >= 0){			/* case 1: */
	   k = MAXIMUM(Int4,(j+k2),k1);
	   G = Pattern(k,G1->nlet); m = G->m;
	   for(m1=G1->m, i=0; i < k1; i++,m1++) m[i] = *m1;
	   m2=G2->m; 		/* see below */
	   end = k2+j;
	   for(i=j; i< end; i++,m2++) m[i] |= *m2;
	} else {			/* case 2: */
	   k = MAXIMUM(Int4,(j+k1),k2);
	   G = Pattern(k,G1->nlet); m = G->m;
	   m1=G1->m; 		/* m1[0]=m1[i-j] where i=j */
	   for(end=k1+j,i=j; i < end; i++,m1++) m[i] = *m1;
	   for(m2=G2->m,i=0; i < k2; i++,m2++) m[i] |= *m2;
	}
	return G;
}

BooLean	CombinablePatterns(ptn_typ G1,ptn_typ G2,Int4 *os)
/* tests whether G1 is a superpattern of (and thus longer) than G2 */
{
	sst_typ	*m1,*m2;
	Int4	i,j,k1=G1->k,k2=G2->k,len;
	BooLean	conflict;
	Int4	nident,pident,start1,start2,score,offset;
	Int4	tmp,spid;

	m2=G2->m; m1=G1->m;
	score = 0; spid=0; *os=k1;
	if(k1 < (len=LengthPattern(G2))) return FALSE;
	for(offset=0; (offset)<=(k1-3); offset++) {
	    start2 = 0; start1 = offset;
	    conflict = FALSE;
	    for(nident=pident=0,i=start1,j=start2; j<k2 && i<k1; i++,j++) {
		if(m1[i] != 0){
			if(m1[i] == m2[j]) nident++;
			else if(m2[j]==(m1[i]&m2[j])) pident++;
			else if(m2[j]!=0){ conflict=TRUE; break; }
		} else if(m2[j] != 0) { conflict = TRUE; break; }
	    }
	    if((tmp=pident + 2*nident) >= 3 && !conflict){
		if(score < tmp){
			spid = pident; score = tmp; *os = offset;
		}else if(score == tmp){
			if(spid > pident){	/* prefer full ident */
				spid = pident; *os = offset;
			} else if(spid == pident){
			   if(abs(*os) > abs(offset)){ /* & short offsets */
				*os = offset;
			   } 
			} /* if spid < pident then don't change */
		   }
	    }
	}
	if(k1 < (len + *os)) return FALSE;
	if(score > 0) return TRUE;
	else return FALSE;
}

Int4	ResPatterns(Int4 i,char *L, ptn_typ G, a_type A)
{
	Int4	j,n;
	for(n=0,j=1; j<=nAlpha(A); j++)if(G->m[i] & (1<<j)) L[n++] = j;
	L[n] = -1;
	return n;
}

BooLean	IdenticalPatterns(ptn_typ G1,ptn_typ G2)
{
	sst_typ	*m1,*m2;
	Int4	i,k;
	
	if(G1->k != G2->k) return FALSE;
	k = (Int4) G1->k; m1=G1->m,m2=G2->m;
	for(i=0; i<k; i++,m1++,m2++) { if(*m1 != *m2) return FALSE; }
	return TRUE;
}

/************************************************
    G0    >    G1 	-> return T else return F
-----------------------------------
 ..PP{YF} >  .NPPY 	-> return T 
 ..DPPY   <  ..PPY 	-> return F
 .{ND}PPY >  .DPPY 	-> return T  

i.e., if all sets in G1 are subsets of G0 return TRUE
	m1 00100000000000000000
		  &
	m0 00100010000010000000

	=  00100000000000000000 == m1 if m1 is a subset of m0
 note: assumes '.' == 0 == {A,C,...,Y}

	yyyyy...  .yyyyyy.. ..yyyyy.   ...yyyyy  k0 = 5
	xxxxxxxx  xxxxxxxx  xxxxxxxx   xxxxxxxx  k1 = 8
offset: 0          1          2           3

********************************************************************/
BooLean	SubPattern(ptn_typ G0,ptn_typ G1)
/* if G1 is a subpattern of G0 then return true */
{
	register sst_typ	*m1=G1->m,*m0=G0->m;
	register Int4		i,j,k1,k0,offset;

	/*** SHORTEN PATTERNS AS MUCH AS POSSIBLE **/
	for(k1=G1->k; m1[k1-1] == 0; ) {
	   k1--; if(k1==0) pattern_error("zero length pattern?!?");
	}
	for(k0=G0->k; m0[k0-1] == 0; ) {
	   k0--; if(k0==0) pattern_error("zero length pattern?!?");
	}
	if(k1 >= k0){
	  for(offset=0; offset <= (k1-k0); offset++) {
	    for(i=0, j=offset; i<k0; i++,j++) {
		if(m0[i] != 0){		/* i.e., G0(i) != {A,C,...,Y} */
		   if(m1[j]==0) break;  /* G0(i) != '.' && G1(j) == '.' */
                   else if(m1[j]==(m1[j]&m0[i])) {
			if(i==(k0-1)) return TRUE;
		   } else break;
		}
	    }
	  }
	  return FALSE;
	} else return FALSE;
}

ptn_typ MultiVar2Pattern(Int4 n, ptn_typ G)
/* Add n variable sets to pattern G */
{Int4 i;for(i=0;i<n;i++) Var2Pattern(G); return G;}

ptn_typ Var2Pattern(ptn_typ G) 
/* Add a variable set '?' = {{x} element A} to the next position of G */
{       
        if(G==NULL) return NULL; 
        if((G->i++) < G->k) return G; else return NilPattern(G); 
}

ptn_typ MultiUSet2Pattern(Int4 n, ptn_typ G)
/* Add n universal sets to pattern G */
{Int4 i;for(i=0;i<n;i++)USet2Pattern(G);return G;}

ptn_typ USet2Pattern(ptn_typ G) 
/* Add a universal set to next position in pattern G */
{       
        if(G==NULL) return NULL; 
        if(G->i < G->k) { G->m[G->i++]=UPattern(G); return G; }
        else return NilPattern(G); 
}

Int4     LengFormatPattern(ptn_typ G)
/* return the length of a pattern if blanks are inserted for formating */
{
        Int4     j,k;
        BooLean up;
 
	if(G == NULL) return 0;
        for(j=k=0,up=FALSE; j < (Int4) G->k;j++,k++) {
           if(NotEqUPattern(j,G) && NotEmptyPattern(j,G)) {
                if(!up){ k++; up=TRUE;}
           } else if(up){ k++; up=FALSE; }
        } 
        return k; 
} 

Int4	GetPattern(char **pattern, ptn_typ G, a_type A)
/* return a 2 dimensional character array for the 1-pattern pattern G */
{
	Int4	i,c,d,tmax,k;
	UInt4	s,u=UPattern(G);

	if(G == NULL) return 0;
	k = LengthPattern(G);
	for(c=0; c<(Int4)nAlpha(A);c++)
		for(i=0;i < k;i++) pattern[c][i]=' ';
	for(i=0,tmax=1;(Int4)i< k;i++) {
	    if((s=G->m[i])==u || s==0) pattern[0][i]='.';
	    else {
		for(d=0,c=1;c<=nAlpha(A);c++)
			if(s & (1<<c)) pattern[d++][i]=AlphaChar(c,A); 
		tmax = MAXIMUM(Int4,tmax,d);
	    }
	}
	for(c=0;c<tmax;c++) pattern[c][i]='\0';
	return tmax;
}

void	pattern_error(char *s){fprintf(stderr,"Pattern: %s\n",s);exit(1);}


