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

#include "cmsa.h"

/*********************************************************************

	  1      2            3        4
       0-===-----===---------===------===---5   aln1
         ...     ....      ......     ....
       ./   \.../    \..../      \.../    \..   
       0-----===------===----===---===------5   aln2
              1        2      3     4

       01   12  12   23   23     34  34   45	cross1 = 1,2,2,3,3,4,4,5
       --   --  --   --   --     --  --   --	crossover positions 
       01   01  12   12   23     34  45   45	cross2 = 1,1,2,2,3,4,5,5
						     c = 1,2,3,4,5,6,7,8

 aln1=[0, 1,  0,  2,   0,     3,    0, 4,   5]
 aln2=[0, 0,  1,  0,   2,     3,    4, 0,   5]

for(c=1; c<=numcross; c++)...
 case A:
  if(cross1[c] < cross1[c+1] && cross2[c] < cross2[c+1])
     then blocks cross1[c]  & cross2[c] overlap 
	-->	aln1[b] = cross1[c];
		aln2[b++] = cross2[c];

 case B:
  if(cross1[c] == cross1[c+1] && cross2[c] < cross2[c+1])
     then block cross1[c] > cross2[c] 
	-->	aln1[b] = 0;
		aln2[b++] = cross2[c];

 case C:
  if(cross1[c] < cross1[c+1] && cross2[c] == cross2[c+1]) 
     then block cross1[c] > cross2[c] 
	-->	aln1[b] = cross1[c];
		aln2[b++] = 0;

 case D:
  if(cross1[c] == cross1[c+1] && cross2[c] == cross2[c+1]) 
	SHOULD NOT HAPPEN: print_error();

 But may have:  
                2      3  4  5      6
		==----***===###----==
                   X            X
                ==---***===###----==
                3     4  5  6      7
                 23/34         56/67
 so will get:
    
 aln1=[0, ...,  3,  6, ...]
 aln2=[0, ...,  4,  7, ...]

   implies that must go with either 345 from aln1 or 456 from aln2.

/**********************************************************************/
static Int4	gaps_blks4cross_cmsa(Int4 numcross, Int4 *aln1, Int4 *aln2,
	Int4 *cross1, Int4 *cross2)
{
	Int4	b,c,t1,t1n,t2,t2n;

	for(b=0,c=1; c<numcross; c++){
	  t1=cross1[c]; t1n=cross1[c+1];
	  t2=cross2[c]; t2n=cross2[c+1];
  	  if(t1 < t1n && t2 < t2n){
		b++; aln1[b] = t1; aln2[b] = t2;
  	  } else if(t1 == t1n && t2 < t2n){
		b++; aln1[b] = 0; aln2[b] = t2;
  	  } else if(t1 < t1n && t2 == t2n){
		b++; aln1[b] = t1; aln2[b] = 0;
	  } else print_error("gaps_blks4cross_cmsa( ) input error");

	}
	return b;
}

/**********************************************************************
  Find all crossover points between alignment gaps in ma1 & ma2.

A crossover occurs between t1 and t2 if for ith aligned sequence:

	(s2[i] < e1[i] && s1[i] < e2[i])

              s1          e1
	=======------. .--==================   aln1
                      X
	===========--' `------==============   aln2
                  s2          e2

 *********************************************************************

 repeats: -==---===--			aln1	(This needs to be able to
          -==---===--(-==---===--)		 deal with repeats.)
         (-==---===--)

          -==---===--			aln2
          -==---===--(-==---===--)
         (-==---===--)

 Only need crossovers BETWEEN blocks (not at ends).
 at most ntyp1 + ntyp2 types;
/**********************************************************************/
static Int4	GCrossPtsCMSA(cma_typ ma1, cma_typ ma2,Int4 *aln1,Int4 *aln2)
{
    Int4	t,t1,t2,s1,s2,e1,e2,c,s,s1o,s2o,i,pt,ncross;
    Int4	n,n1,n2,N,ntyp1,ntyp2,len1,len2;
    Int4	**start1,**start2,**end1,**end2,**tmp;
    Int4	**Tstart1,**Tstart2,**Tend1,**Tend2;
    Int4	cross1[60],cross2[60];
    ss_type	data1 = DataCMSA(ma1);
    ss_type	data2 = DataCMSA(ma2);
    st_type	S1,S2;
    double	bad,total;
    // double	cutoff=0.000010;
    double	cutoff=0.10;
    Int4	minoverlap=10;
    gss_typ	*gss1,*gss2;

    gss1=gssCMSA(ma1); gss2=gssCMSA(ma2);
    S1=SitesCMSA(ma1); ntyp1 = nTypeSites(S1);
    S2=SitesCMSA(ma2); ntyp2 = nTypeSites(S2);
    N = NSeqsSeqSet(data1);
    NEWP(start1,N+2,Int4); NEWP(start2,N+2,Int4);
    NEWP(end1,N+2,Int4);  NEWP(end2,N+2,Int4);
    for(n=1; n<=N; n++){
	NEW(start1[n],55,Int4); NEW(start2[n],55,Int4);
	NEW(end1[n],55,Int4);  NEW(end2[n],55,Int4);
    }
    for(n=1; n<=N; n++){	// Delineate all gaps between sites.
	n1=PosSites(n, end1[n], S1);
	n2=PosSites(n, end2[n], S2); // ASSUME NO REPEATS.
	assert(n1 == ntyp1); assert(n2 == ntyp2);
	end1[n][n1+1]=SqLenSeqSet(n,data1)+1; 
	end2[n][n2+1]=SqLenSeqSet(n,data2)+1;	
	start1[n][1] = 0; start2[n][1] = 0;
	for(i=1,t=2; i<=n1; i++,t++){
		start1[n][t] = end1[n][t-1] + SiteLen(t-1,S1) - 1;
	}
	for(i=1,t=2; i<=n2; i++,t++){
		start2[n][t] = end2[n][t-1] + SiteLen(t-1,S2) - 1;
	}
    }
// NEW PART...  // Convert sites from FakeSeq to RealSeq.
      NEWP(Tstart1,N+2,Int4); NEWP(Tstart2,N+2,Int4);
      NEWP(Tend1,N+2,Int4);  NEWP(Tend2,N+2,Int4);
      for(n=1; n<=N; n++){
	NEW(Tstart1[n],55,Int4); NEW(Tstart2[n],55,Int4);
	NEW(Tend1[n],55,Int4);  NEW(Tend2[n],55,Int4);

	Tstart1[n][1] = gss1->TruePos(n,1)-1;
	Tstart2[n][1] = gss2->TruePos(n,1)-1;
	Tend1[n][1]=gss1->TruePos(n,end1[n][1]);
	Tend2[n][1]=gss2->TruePos(n,end2[n][1]);

        for(t1=2; t1<=ntyp1; t1++){
	   Tstart1[n][t1]=gss1->TruePos(n,start1[n][t1]);
	   Tend1[n][t1]=gss1->TruePos(n,end1[n][t1]);
	}
	Tstart1[n][t1]=gss1->TruePos(n,start1[n][t1]);
	Tend1[n][t1]=gss1->TruePos(n,end1[n][t1]-1)+1;

        for(t2=2; t2<=ntyp2; t2++){
	   Tstart2[n][t2]=gss2->TruePos(n,start2[n][t2]);
	   Tend2[n][t2]=gss2->TruePos(n,end2[n][t2]);
	}
	Tstart2[n][t2]=gss2->TruePos(n,start2[n][t2]);
	Tend2[n][t2]=gss2->TruePos(n,end2[n][t2]-1)+1;
      }

    // Now check for all crossover points.
    for(ncross=0,t1=1; t1<=ntyp1+1; t1++){
      for(t2=1; t2<=ntyp2+1; t2++){
	total=(double)N;
    	for(bad=0.0, n=1; n<=N; n++){
	   s1=Tstart1[n][t1]; e1=Tend1[n][t1]; 
	   s2=Tstart2[n][t2]; e2=Tend2[n][t2];
	   if(!(s2<e1 && s1<e2)){ 
/*********************************************************************
	(s2[i] < e1[i] && s1[i] < e2[i])   minoverlap==3
                 s1                             e1
   ================--===========   ===========--===========-----
   ===========--===========-----   =================--==============
     e2-s1=2    e2   (overlap=3)    e1-s2=3       s2  (overlap = 4)
 *********************************************************************/
#if 0		// QUANTIFY HOW BAD FOR GAP PENALTIES!?
		if((e2-s1) < minoverlap && (e1-s2) < minoverlap){
		   bad+=1.0; // printf("...fail\n\n");**/ break;
		   if((bad/total) > cutoff) break;  // if > 10 % conflicts
		} else { bad+=total; break; }
#endif
#if 1
		bad+=1.0; // printf("...fail\n\n");**/ break;
		if((bad/total) > cutoff) break;  // if > 10 % conflicts
#endif
	   }
#if 0
	   else {
	   	printf("1. %2d(%d): %d..%d\n",n,t1,s1,e1);
	   	printf("2. %2d(%d): %d..%d\n",n,t2,s2,e2);
	   	printf("\tcrosspt = %2d\n",MINIMUM(Int4,e1,e2));
	   } 	// DEBUG: potential crossover point.
#endif
	} // printf("\n"); 
	if((bad/total) <= cutoff){ // then store block type following crossover.
	    ncross++; cross1[ncross]=t1; cross2[ncross]=t2;
	}
      }
    }
#if 1 // DEBUG 
    fprintf(stderr,"%d crossover points...\n",ncross);
    fprintf(stderr,"cross1: ");
    for(n=1; n<=ncross; n++) fprintf(stderr,"%2d ",cross1[n]);
    fprintf(stderr,"\ncross2: ");
    for(n=1; n<=ncross; n++) fprintf(stderr,"%2d ",cross2[n]);
    fprintf(stderr,"\n\n");
#endif

    for(n=1; n<=N; n++){
	free(start1[n]); free(start2[n]); free(end1[n]); free(end2[n]);
	free(Tstart1[n]); free(Tstart2[n]); free(Tend1[n]); free(Tend2[n]);
    }
    free(start1); free(start2); free(end1); free(end2);
    free(Tstart1); free(Tstart2); free(Tend1); free(Tend2);
    return gaps_blks4cross_cmsa(ncross, aln1, aln2, cross1, cross2);
}

/************************************************************************
 Find all possible recombination events using a depth-first-search.  

aln1:  0  1  2  3  4 
aln2:  1  2  3  0  5 

/************************************************************************/
static void	gaps_dfs_config_cmsa(Int4 *aln1, Int4 *aln2, Int4 **config, 
	Int4 depth, Int4 *nc, Int4 *rc, Int4 len, Int4 nrc)
{
	Int4	m,n,i,t;
	
	if(depth == 0){
#if 1	// NEW (HOPEFULLY FIXED) CODE.
	   n = *nc; n++;	// increment number of configurations.
	   if(nrc == 0) { config[n][0] = 0; return; }
	   // if(nrc == 0) { config[n][0] = 0; *nc=n; return; }
	   Int4 j;
	   if(rc[1] < 0){
		t=abs(rc[1]); config[n][0] = t-1;
		for(j=1; j < t; j++) config[n][j]=-j;
	   } else {
		t=rc[1]; config[n][0] = t-1;
		for(j=1; j < t; j++) config[n][j]=j;
	   } config[n][0] += nrc;
	   for(i=1; i<= nrc; j++,i++){ config[n][j] = rc[i]; }
	   *nc = n;
#endif
#if 0	// OLD CODE WITH BUG (IN UNGAPPED VERSION TOO!!!)
		n = *nc; n++;	/** number of configurations **/
		config[n][0] = nrc;
		for(i=1; i<= nrc; i++){ config[n][i] = rc[i]; }
		*nc = n;
#endif
	} else { // aln1 case.
	   m=nrc; m++; n=len; n++;
	   if((rc[m]=aln1[n]) !=0){	/** rc > 0 -> aln1 **/
		t=aln1[n];	/** current type of block **/
		for(i=1; aln1[n+i] < t; ) i++;
		i=aln1[n+i];	/** i = next type of block **/
		for(t++; t < i; t++) {	/** record the block **/
			m++; rc[m] = t;
		}
	   } else m--;	/** don't include zero's in configuration **/
	   gaps_dfs_config_cmsa(aln1, aln2, config, depth-1, nc,rc,len+1,m);
	   // aln2 case.
	   m=nrc; m++; n=len; n++;
	   if((rc[m]=-aln2[n]) !=0){	/** rc < 0 -> aln2 **/
		t=aln2[n];	/** current type of block **/
		for(i=1; aln2[n+i] < t; ) i++;
		i=aln2[n+i];	/** i = next type of block **/
		for(t++; t < i; t++) {	/** record the block **/
			m++; rc[m] = -t;
		}
	   } else m--;
	   gaps_dfs_config_cmsa(aln1, aln2, config, depth-1, nc,rc,len+1,m);
	}
}

cma_typ	GRecombineCMSA(cma_typ ma1, cma_typ ma2)
{
	Int4	ntyp1,ntyp2,num,i,j,nc,n,*bconfig;
	Int4	*aln1,*aln2,**config,rc[100];
	double	map,bmap;
	double	*Map,sum,rand;
	cma_typ ma;
	ss_type	data;
	fm_type	*model,M,M1,M2;
	BooLean	okay;

	// if(DataCMSA(ma1) != DataCMSA(ma2)) return NULL; // distinct gapped SqSets
	// if(TrueDataCMSA(ma1) != TrueDataCMSA(ma2)) return NULL; // distinct SqSets
	ntyp1 = nBlksCMSA(ma1); ntyp2 = nBlksCMSA(ma2);
	NEW(aln1,ntyp1+ntyp2+5,Int4); NEW(aln2,ntyp1+ntyp2+5,Int4);
	num = GCrossPtsCMSA(ma1, ma2,aln1,aln2);
#if 1   //************* DEBUG *******************
        printf("ncross = %d\n\naln1: ",num);
        for(i=1;i <= num; i++){ printf("%2d ",aln1[i]); } printf("\naln2: ");
        for(i=1;i <= num; i++){ printf("%2d ",aln2[i]); } printf("\n\n");
#endif

	nc = (Int4) floor(pow(2.,num)+0.1);
	NEWP(config,nc+2,Int4);
	for(i=1; i <= nc; i++) NEW(config[i],62,Int4);

	rc[0]=0; n=0; aln1[0]=0; aln2[0]=0; 
	aln1[num+1]=ntyp1+1; aln2[num+1]=ntyp2+1;
	gaps_dfs_config_cmsa(aln1, aln2, config, num, &n,rc,0,0);

#if 1
        fprintf(stderr,"%d recombinants: \n",n);
#endif
#if 0
for(j=1; j <= nc; j++){
bconfig=config[j];
fprintf(stderr,"config = [ %d ",bconfig[0]);
for(i=1; i<=bconfig[0]; i++){ fprintf(stderr,"%d ",bconfig[i]); }
fprintf(stderr,"]\n");
}
#endif
	for(bmap=0,i=1; i<=n; i++){
	   map=GRecombinantMapCMSA(ma1, ma2, config[i]);
	   if(map > bmap){ 
		bmap=map; bconfig=config[i]; 
	   } /** pick a parent if given a choice **/
	   else if(map == bmap){  
    	        for(okay=FALSE,j=1; j<=config[i][0]; j++){ 
			if(config[i][j] != j){ okay=TRUE; break; }
		}
		if(okay){
    	           for(okay=FALSE,j=1; j<=config[i][0]; j++){ 
			if(config[i][j] != -j){ okay=TRUE; break; }
		   }
		}
		if(!okay) bconfig=config[i];
	   }
	}
	if(bmap > 0){ /** see whether either parent is optimal **/
    	    for(okay=FALSE,i=1; i<=bconfig[0]; i++){ 
#if 1
fprintf(stderr,"A: bconfig[%d] = %d\n",bconfig[i],i);
#endif
		if(bconfig[i] != i){ okay=TRUE; break; }
	    } // config = [ n, 1, 2, 3, 4, ..., n]
	    if(okay){
    	      for(okay=FALSE,i=1; i<=bconfig[0]; i++){ 
#if 1
fprintf(stderr,"B: bconfig[%d] = %d\n",bconfig[i],-i);
#endif
		if(bconfig[i] != -i){ okay=TRUE; break; }
	      }
	    } // config = [ n,-1,-2,-3,-4, ...,-n] 
	} else okay=FALSE;
#if 0
if(okay){
	fprintf(stderr,"bconfig = [ %d ",bconfig[0]);
	for(i=1; i<=bconfig[0]; i++){ fprintf(stderr,"%d ",bconfig[i]); }
	fprintf(stderr,"]\n");
}
#endif
	if(okay) ma=RecombineGapsCMSA(ma1,ma2,bconfig);
	else ma = NULL;
	for(i=1; i <= nc; i++) free(config[i]);
	free(config); free(aln1); free(aln2);
	return ma;
}

double	GRecombinantMapCMSA(cma_typ ma1, cma_typ ma2, Int4 *config)
/*********************************************************************
  Compute the posterior probability for recombinant ma1 & ma2.

	  1      2            3        4
	-===-----===---------===------===---	aln1
        ....     ....      ......     ......
            \.../    \..../      \.../     	(crossover points)
	-----===------===----===---===------ 	aln2
             -1       -2     -3    -4

  config = [ n, 1,-1, 2,-2, 3,-4, 4] 	n = 7 = number of recombinant blocks.

 **********************************************************************/
{
    Int4        N,n,len,i,j,t,T,s,ncol,off_col,on_col,nblks;
    Int4 	totsite,b,*counts,nonsite[30],*site,**observed;
    double      MAP,weight,total,Ps[30],npseudo,*freq,v;
    ss_type	data = TrueDataCMSA(ma1);
    ss_type	data1 = DataCMSA(ma1);
    ss_type	data2 = DataCMSA(ma2);
    st_type	S1,S2;
    fm_type	M,*model1,*model2;
    a_type	A;
    double      indel_penalty=0.0;

    indel_penalty=IndelPenaltySeqSet(data1);
    // indel_penalty=IndelPenaltySeqSet(data2);
    // for(i=1; i<=config[0]; i++){ printf(",%d",config[i]); } printf("]\n");

    model1=ModelsCMSA(ma1); model2=ModelsCMSA(ma2);
    S1=SitesCMSA(ma1); S2=SitesCMSA(ma2);
    nblks=config[0];
    A = SeqSetA(data); counts = CntsSeqSet(data);
    freq=tFreqSeqSet(data); N = NSeqsSeqSet(data);
    for(len=0, i=1; i <= nblks; i++) {
		if((t=config[i]) < 0){ t = abs(t); M = model2[t]; }
		else M = model1[t];
		n = MaxLenFModel(M);
		if(n > len) len = n;
    }
    NEWP(observed,len+1,Int4);
    for(weight=0.0, s=1; s<=N; s++){
	for(T=0, i=1; i <= nblks; i++){
		if((t=config[i]) < 0){ t = abs(t); M = model2[t]; }
		else M = model1[t];
		T -= LenFModel(M);
	}
	T += SqLenSeqSet(s,data) + nblks;
	weight += lnbico(T, nblks);
    }
    npseudo = (double)N*PSEUDO_CMSA;
    if(npseudo > 20.0) npseudo = 20.0;
    for(b=1;b<=nAlpha(A);b++){nonsite[b]=counts[b];Ps[b]=npseudo*freq[b];}
    for(on_col=off_col=0,MAP=0.0,j=1; j <= nblks; j++) {
	if((t=config[j]) < 0){ t = abs(t); M = model2[t]; }
	else M = model1[t];
	len = ObservedFModel(observed, M);
        ncol = nColsFModel(M); on_col += ncol;
	off_col += LenFModel(M) - ncol;
	totsite = TotSitesFModel(M);
	/*************** NEW: column width weight ******************/
	MAP -= lnbico(len - 2, ncol - 2);
	/***********************************************************/
	/***** MAP *******/
	for(i=1; i <= len; i++) {
	  site = observed[i];
	  if(site != NULL){
	    for(b=1;b<= nAlpha(A); b++) {
	      if(Ps[b] > 0){
		MAP += lngamma((double)site[b] + Ps[b]); 
		nonsite[b] -= site[b];
	      }
	    }
            MAP -= lngamma((double)(totsite-site[0]) + npseudo);
	    // MAP -= lngamma((double)totsite + npseudo);
	  } 
	}
    }
    v = lngamma(npseudo);
    /***** NONSITES MAP *******/
    for(total = 0.0, b=1;b<= nAlpha(A); b++) {
	   if(Ps[b] > 0){
		MAP += lngamma((double) nonsite[b] + Ps[b]);
		total += (double) nonsite[b];
        	v -= lngamma(Ps[b]);   /** term pulling down MAP **/
	   }
    }
    MAP += v *(double) on_col; 
    MAP -= lngamma((double) total + npseudo);
    /***** subtract NULL MAP ******/
    for(total = 0.0, b=1;b<= nAlpha(A); b++) {
	   if(counts[b] > 0.0){
		MAP -= lngamma(counts[b] + Ps[b]);
		total += counts[b];
	   }
    }
    MAP += lngamma((double)total + npseudo);
    /**************** weight for column transfers *************/
    if(nblks > 1) {
	weight +=  lnbico(on_col-nblks-1,nblks-1); 
	if(off_col > 0) weight += lnbico(off_col+nblks-1,nblks-1);
    }
    /**************** return MAP *************/
    free(observed);
    // return (MAP - weight);
    return (MAP - weight - indel_penalty);
}

cma_typ  RecombineGapsCMSA(cma_typ cma1, cma_typ cma2, Int4 *config)
// sample a recombinant alignment from cma1 and cma2.
/*****************************************************************************
 Set *this == to the recombinant of gsq1 & gsq2 at crossover points xoverpt1
  and xoverpt2.  Set pos ==

          1      2            3        4
        -===-----===---------===------===---    aln1
        ....     ....      ......     ......
            \.../    \..../      \.../          (crossover points)
        -----===------===----===--===-------    aln2
             -1       -2     -3    -4

  Int4 *config = [ n, 1,-1, 2,-2, 3,-4, 4] n = 7 = number of recombinant blocks.

          1      2            3        4
        -===-----===---------===------===---    aln1
        .............      ......     ......
                     \..../      \.../          (crossover points)
        -----===------===----===--===-------    aln2
             -1       -2     -3    -4

  Int4 *config = [ 5, 2,-2, 3,-4, 4] 
 *****************************************************************************/
{
	Int4	m,n,i,t,start,trace_length,indels,gapopen,gapextend;
	Int4    score;
	Int4	*newpos,b,*blklen,nblks;
	smx_typ	*smx;
	a_type	A=AlphabetCMSA(cma1);
	e_type	E;
	double	pernats;
	char	*operation;
	cma_typ	cma;
    	gss_typ	*gss,*gss1,*gss2;

    assert(NumSeqsCMSA(cma1)==NumSeqsCMSA(cma2));
    gss1=gssCMSA(cma1); gss2=gssCMSA(cma2); 
    gapopen=gss1->GapOpen(); gapextend=gss1->GapExtend(); pernats=gss1->PerNats();

    nblks=config[0];
    blklen = new Int4 [nblks+2];
    for(b=1; b <=nblks; b++){
    	if((t=config[b]) < 0) blklen[b]= LengthCMSA(-t,cma2);
	else blklen[b]= LengthCMSA(t,cma1);
    }
    cma=EmptyCMSA(nblks,blklen,TrueDataCMSA(cma1),gapopen,gapextend,pernats,
			gss1->LeftFlank(),gss1->RightFlank());
#if 0
PutSeqSetEs(stderr,TrueDataCMSA(cma));
// PutCMSA(stderr,cma);  // NEED TO DEBUG THIS FOR EMPTY SEQS???
if(TRUE) exit(1);
#endif
    NEW(newpos,nblks+2,Int4);  // for new sites.
    // 2. Obtain scoring matrix from rest of alignment.
    NEW(smx,nblks + 2,smx_typ);
    for(n=0,b=1; b<= nblks; b++){
	if((t=config[b]) < 0){		// second alignment...
		smx[b] = GetSmatrixFModel(pernats,ModelCMSA(-t,cma2));
	} else {
		smx[b] = GetSmatrixFModel(pernats,ModelCMSA(t,cma1));
	} n+=LenSMatrix(smx[b]);
    }

    gss=gssCMSA(cma);
    for(Int4 s=1;s<=NumSeqsCMSA(cma);s++){
	// ========== 4. Sample a gapped alignment for sequence s. ===========
        E = gss->TrueSeq(s); 	
	operation=GapAlnTraceSMatrix(gapopen,gapextend,LenSeq(E),
			SeqPtr(E),nBlksCMSA(cma),smx, NULL,&start);
	trace_length=strlen(operation);
#if 0
std::cerr << operation;
        PutSeqInfo2(stderr,E);
        PutGappedSeqAlnSMatrix(stderr, operation, OffSetSeq(E)+start-1, 
                LenSeq(E)-start+1, SeqPtr(E)+start-1, nBlksCMSA(cma), smx);
#endif
	// ========== 5. Create a gapped sequence. ===========
	gsq_typ	*gsq; gsq = new gsq_typ[1];
        gsq->initialize(gss->LeftFlank(),gss->RightFlank(),
				operation,trace_length,start,E,newpos);
	// ========== 7. Add sites to gsq. =========
	ReplaceCMSA(s,gsq,cma); // replace sequence s in CMSA & fmodel.
	for(t=1; t<=nBlksCMSA(cma); t++) AddSiteCMSA(t,s,newpos[t], cma);
	free(operation);
   }
   // ========== 6. Deallocate memory. ===========
   for(m=1; m<=nBlksCMSA(cma); m++) NilSMatrix(smx[m]); 
   free(smx); free(newpos); delete [] blklen;
   return cma;
}

