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
#include "sma.h"

/******************** Genetic Operations ****************************/
cma_typ	FuseBlksCMSA(Int4 x, Int4 maxlen, cma_typ L)
/***********************************************************************
  take two blocks and fuse them into one block.
 ***********************************************************************/
{
    cma_typ	ma;
    Int4	j,t2,t,N,T,T2;
    st_type	S,S2=NULL;
    fm_type	*model;
    BooLean	**null;

    S=SitesCMSA(L); T = nBlksCMSA(L); T2 = T-1;
    if(T < 2 || x < 1 || x >= T) return NULL;
    N = NumSeqsCMSA(L);
#if 0 // new code started on 5/13/04
    S2 = FuseElementsSites2(x, maxlen, S);
#else
    S2 = FuseElementsSites(x, maxlen, S);
#endif
    if(S2 == NULL) return NULL;
    model=ModelsCMSA(L); 
    NEWP(null,T2+2,BooLean); 
    for(t2=t=1,j=1; j<= T2; j++,t2++,t++){ 
	if(j==x){ t++; continue; } // use all columns in split block 
	NEW(null[t2],LenFModel(model[t])+5,BooLean);
	NullSitesFModel(null[t2], model[t]);
    }
    ma=MakeCMSA(S2,null);
    for(t=1; t<=T2; t++) if(null[t]!=NULL)free(null[t]);
    free(null);
#if 0
    PutSites(stderr,x,S,NULL, NULL);
    PutSites(stderr,x+1,S,NULL, NULL);
    PutSites(stderr,x,S2,NULL, NULL);
#endif
    return ma;
}

cma_typ	SplitBlkCMSA(Int4 x, Int4 left_leng, Int4 minlen, cma_typ L)
/***********************************************************************
 Splits a block into left_leng making sure that right_leng >= minlen.
 ***********************************************************************/
{
    cma_typ	ma;
    Int4	j,t2,t,T,T2;
    st_type	S,S2=NULL;
    fm_type	*model;
    BooLean	**null;

    S=SitesCMSA(L); T = nTypeSites(S); T2=T+1;
    if(x < 1 || x > T) return NULL;
    S2 = SplitElementSites(x, left_leng, minlen, S);
    // rest same as below... eventually merge...
    if(S2 == NULL) return NULL;
    model=ModelsCMSA(L); 
    NEWP(null,T2+2,BooLean); 
    if(x==0) t2=2; else t2=1;
    for(t=1,j=1; j<= T; j++,t2++,t++){ 
	if(j==x){ t2++; continue; } // use all columns in split block 
	NEW(null[t2],LenFModel(model[t])+5,BooLean);
	NullSitesFModel(null[t2], model[t]);
    }
    ma=MakeCMSA(S2,null);
    for(t=1; t<=T2; t++) if(null[t]!=NULL)free(null[t]);
    free(null);
    return ma;
}

cma_typ	SplitBlkCMSA(Int4 x, Int4 minlen, cma_typ L)
/***********************************************************************
 Splits a block about evenly in half...
 ***********************************************************************/
{
    cma_typ	ma;
    Int4	j,t2,t,T,T2;
    st_type	S,S2=NULL;
    fm_type	*model;
    BooLean	**null;

    S=SitesCMSA(L); T = nTypeSites(S); T2=T+1;
    if(x < 1 || x > T) return NULL;
    S2 = SplitElementSites(x, minlen, S);
    if(S2 == NULL) return NULL;
    model=ModelsCMSA(L); 
    NEWP(null,T2+2,BooLean); 
    if(x==0) t2=2; else t2=1;
    for(t=1,j=1; j<= T; j++,t2++,t++){ 
	if(j==x){ t2++; continue; } // use all columns in split block 
	NEW(null[t2],LenFModel(model[t])+5,BooLean);
	NullSitesFModel(null[t2], model[t]);
    }
    ma=MakeCMSA(S2,null);
    for(t=1; t<=T2; t++) if(null[t]!=NULL)free(null[t]);
    free(null);
    return ma;
}

cma_typ	AddBlkCMSA(Int4 x, Int4 lenx, cma_typ L)
/***********************************************************************
  Create and return L2: an exact copy of L but with an additional block 
  randomly inserted between motif x and motif x+1 in each sequence. 
  If x < 1 then new motif is added before the first site.
  If x >= ntyp then new motif is added after the last site.
  if there is not room to add an element of length lenx or if
  for some other reason the operation can't be done then NULL is returned.
 ***********************************************************************/
{
    cma_typ	ma;
    Int4	j,t2,t,N,T,T2,min_free=10,min_squeeze;
    st_type	S,S2=NULL;
    fm_type	*model;
    BooLean	**null;

    S=SitesCMSA(L); T = nBlksCMSA(L); T2=T+1; N=NumSeqsCMSA(L);
    min_squeeze = (Int4) ((double)N*0.20);
    S2 = AddElementSites(lenx, x, min_free, min_squeeze, S);
    if(S2 == NULL) return NULL;
    model=ModelsCMSA(L); 
    NEWP(null,T2+2,BooLean); 
    if(x==0) t2=2; else t2=1;
    for(t=1,j=1; j<= T; j++,t2++,t++){ 
	NEW(null[t2],LenFModel(model[t])+5,BooLean);
	NullSitesFModel(null[t2], model[t]);
	if(j==x) t2++; // use contiguious columns in new block
    }
    ma=MakeCMSA(S2,null);
    for(t=1; t<=T2; t++) if(null[t]!=NULL)free(null[t]);
    free(null);
    return ma;
}

Int4	blks4cross_cmsa(Int4 numcross, Int4 *aln1, Int4 *aln2,
	Int4 *cross1, Int4 *cross2)
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
	  } else print_error("blks4cross_cmsa( ) input error");

	}
	return b;
}

/**********************************************************************
  Find all crossover points between alignment gaps in ma1 & ma2.

  A crossover occurs between t1 and t2 if for all i=sequences in
  the alignment:

	(s2[i] < e1[i] && s1[i] < e2[i])

              s1          e1
	=======-----------==================   aln1
                      X
	===========-----------==============   aln2
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
Int4	CrossPtsCMSA(cma_typ ma1, cma_typ ma2,Int4 *aln1,Int4 *aln2)
{
    Int4	t,t1,t2,s1,s2,e1,e2,i,ncross;
    Int4	n,n1,n2,N,ntyp1,ntyp2;
    Int4	**start1,**start2,**end1,**end2;
    Int4	cross1[60],cross2[60];
    ss_type	data = DataCMSA(ma1);
    st_type	S1,S2;
    BooLean	okay;

    S1=SitesCMSA(ma1); ntyp1 = nTypeSites(S1);
    S2=SitesCMSA(ma2); ntyp2 = nTypeSites(S2);
    N = NSeqsSeqSet(data);
    NEWP(start1,N+2,Int4); NEWP(start2,N+2,Int4);
    NEWP(end1,N+2,Int4);  NEWP(end2,N+2,Int4);
    for(n=1; n<=N; n++){
	NEW(start1[n],55,Int4); NEW(start2[n],55,Int4);
	NEW(end1[n],55,Int4);  NEW(end2[n],55,Int4);
    }
    for(n=1; n<=N; n++){	// Delineate all gaps between sites.
	n1=PosSites(n, end1[n], S1);
	n2=PosSites(n, end2[n], S2);
				/** ASSUME NO REPEATS for NOW **/
	if(n1 != ntyp1) print_error("CrossPtsCMSA() input error");
	if(n2 != ntyp2) print_error("CrossPtsCMSA() input error");
	end1[n][n1+1]=SqLenSeqSet(n,data)+1; 
	end2[n][n2+1]=SqLenSeqSet(n,data)+1;	
	start1[n][1] = 0; start2[n][1] = 0;
	for(i=1,t=2; i<=n1; i++,t++){
		start1[n][t] = end1[n][t-1] + SiteLen(t-1,S1) - 1;
	}
	for(i=1,t=2; i<=n2; i++,t++){
		start2[n][t] = end2[n][t-1] + SiteLen(t-1,S2) - 1;
	}
    }
    // Now check for all crossover points.
    for(ncross=1,t1=1; t1<=ntyp1+1; t1++){
      for(t2=1; t2<=ntyp2+1; t2++){
    	for(okay=TRUE, n=1; n<=N; n++){
	   s1=start1[n][t1]; e1=end1[n][t1]; 
	   s2=start2[n][t2]; e2=end2[n][t2];
	   if(!(s2<e1 && s1<e2)){ okay=FALSE; /*printf("...fail\n\n");**/ break;}
#if 0
	   else {
	   	printf("1. %2d(%d): %d..%d\n",n,t1,s1,e1);
	   	printf("2. %2d(%d): %d..%d\n",n,t2,s2,e2);
	   	printf("\tcrosspt = %2d\n",MINIMUM(Int4,e1,e2));
	   } 	/** DEBUG: potential crossover point **/
#endif
	}
	// printf("\n"); 
	if(okay){ // then store block type following crossover.
	    cross1[ncross]=t1;	cross2[ncross]=t2; ncross++;
	}
      }
    }
#if 0
    /**** DEBUG ****/
    fprintf(stderr,"cross1: ");
    for(n=1; n< ncross; n++) fprintf(stderr,"%2d ",cross1[n]);
    fprintf(stderr,"\ncross2: ");
    for(n=1; n< ncross; n++) fprintf(stderr,"%2d ",cross2[n]);
    fprintf(stderr,"\n\n");
    /**** DEBUG ****/
#endif

    for(n=1; n<=N; n++){
	free(start1[n]); free(start2[n]); free(end1[n]); free(end2[n]);
    }
    free(start1); free(start2); free(end1); free(end2);
    return blks4cross_cmsa(ncross-1, aln1, aln2, cross1, cross2);
}

void	dfs_config_cmsa(Int4 *aln1, Int4 *aln2, Int4 **config, Int4 depth,
	Int4 *nc, Int4 *rc, Int4 len, Int4 nrc)
/************************************************************************
 Find all possible recombination events using a depth-first-search.  

aln1:  0  1  2  3  4 
aln2:  1  2  3  0  5 

/************************************************************************/
{
	Int4	m,n,i,t;
	
	if(depth == 0){
		n = *nc; n++;	/** number of configurations **/
		config[n][0] = nrc;
		for(i=1; i<= nrc; i++){ config[n][i] = rc[i]; }
		*nc = n;
	} else {
	   /*** aln1 case **/
	   m=nrc; m++; n=len; n++;
	   if((rc[m]=aln1[n]) !=0){	/** rc > 0 -> aln1 **/
		t=aln1[n];	/** current type of block **/
		for(i=1; aln1[n+i] < t; ) i++;
		i=aln1[n+i];	/** i = next type of block **/
		for(t++; t < i; t++) {	/** record the block **/
			m++; rc[m] = t;
		}
	   } else m--;	/** don't include zero's in configuration **/
	   dfs_config_cmsa(aln1, aln2, config, depth-1, nc,rc,len+1,m);

	   /*** aln2 case **/
	   m=nrc; m++; n=len; n++;
	   if((rc[m]=-aln2[n]) !=0){	/** rc < 0 -> aln2 **/
		t=aln2[n];	/** current type of block **/
		for(i=1; aln2[n+i] < t; ) i++;
		i=aln2[n+i];	/** i = next type of block **/
		for(t++; t < i; t++) {	/** record the block **/
			m++; rc[m] = -t;
		}
	   } else m--;
	   dfs_config_cmsa(aln1, aln2, config, depth-1, nc,rc,len+1,m);
	}
}

cma_typ	RecombineCMSA(cma_typ ma1, cma_typ ma2)
{
	Int4	ntyp1,ntyp2,num,i,j,nc,n,*bconfig;
	Int4	*aln1,*aln2,**config,rc[100];
	double	map,bmap;
	cma_typ ma;
	st_type	S1,S2;
	BooLean	okay;

	if(DataCMSA(ma1) != DataCMSA(ma2)){
		return GRecombineCMSA(ma1,ma2);
		//  return NULL; // distinct gapped SqSets
	}
	// NOTE: This ^^^ needs to be fixed for recombine_msa to work!!
	S1 = SitesCMSA(ma1); ntyp1 = nTypeSites(S1); 
	S2 = SitesCMSA(ma2); ntyp2 = nTypeSites(S2); 
	NEW(aln1,ntyp1+ntyp2+5,Int4);
	NEW(aln2,ntyp1+ntyp2+5,Int4);
	num = CrossPtsCMSA(ma1, ma2,aln1,aln2);
#if 0
	/************* DEBUG *******************/
        printf("ncross = %d\n",num);
        printf("\n\naln1: ");
        for(i=1;i <= num; i++){ printf("%2d ",aln1[i]); }
	printf("\naln2: ");
        for(i=1;i <= num; i++){ printf("%2d ",aln2[i]); }
	printf("\n\n");
	/************* DEBUG *******************/
        PutTypeSites(stdout,S1);
        PutTypeSites(stdout,S2);
	/************* DEBUG *******************/
#endif
	nc = (Int4) floor(pow(2.,num)+0.1);
	NEWP(config,nc+2,Int4);
	for(i=1; i <= nc; i++) NEW(config[i],62,Int4);

	rc[0]=0; n=0; aln1[0]=0; aln2[0]=0; 
	aln1[num+1]=ntyp1+1; aln2[num+1]=ntyp2+1;
	dfs_config_cmsa(aln1, aln2, config, num, &n,rc,0,0);

#if 0
        fprintf(stderr,"%d recombinants: \n",n);
#endif
	for(bmap=0,i=1; i<=n; i++){
	   map=RecombinantMapCMSA(ma1, ma2, config[i]);
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
#if 0
	/***** NEW: replace above (sampling) *****/
	NEW(Map, n+2, double);
	for(bmap=0,i=1; i<=n; i++){
		Map[i]=RecombinantMapCMSA(ma1, ma2, config[i]);
		if(Map[i] > bmap){ bmap = Map[i]; }
		}
	}
	for(sum=0.0,i=1; i<=n; i++) sum+=pow(2,Map[i]/bmap);
	do{
	   rand=sum*(double)Random()/(double)RANDOM_MAX;
	} while(rand >= sum);
	for(i=1; i<=n; i++){
	   rand-=pow(2,Map[i]/bmap);
	   if(rand <= 0.0){ bmap=Map[i]; bconfig=config[i]; break; }
	} /***** NEW *****/
#endif
	if(bmap > 0){ /** see whether either parent is optimal **/
    	    for(okay=FALSE,i=1; i<=bconfig[0]; i++){ 
#if 0
fprintf(stderr,"A: bconfig[%d] = %d\n",bconfig[i],i);
#endif
		if(bconfig[i] != i){ okay=TRUE; break; }
	    } /**** config = [ n, 1, 2, 3, 4, ..., n] ****/
	    if(okay){
    	      for(okay=FALSE,i=1; i<=bconfig[0]; i++){ 
#if 0
fprintf(stderr,"B: bconfig[%d] = %d\n",bconfig[i],-i);
#endif
		if(bconfig[i] != -i){ okay=TRUE; break; }
	      }
	    } /**** config = [ n,-1,-2,-3,-4, ...,-n] ****/
	} else okay=FALSE;
#if 0
fprintf(stderr,"bconfig = [ %d ",bconfig[0]);
for(i=1; i<=bconfig[0]; i++){ fprintf(stderr,"%d ",bconfig[i]); }
fprintf(stderr,"]\n");
#endif
	if(okay) ma = MkRecombinantCMSA(ma1, ma2, bconfig);
	else ma = NULL;
	for(i=1; i <= nc; i++) free(config[i]);
	free(config); free(aln1); free(aln2);
	return ma;
}

cma_typ	MkRecombinantCMSA(cma_typ ma1, cma_typ ma2, Int4 *config)
/*********************************************************************
    recombine ma1 & ma2 using the configuration config.

 *********************************************************************/
{
	Int4	N,j,t,n,nblks,tmp1[30],tmp2[30];
	Int4	s,*len_elem,n1,n2;
	fm_type	M,*model1,*model2,*model;
	st_type	S,S1,S2;
	BooLean	**null;
	cma_typ	ma;

	model1=ModelsCMSA(ma1); model2=ModelsCMSA(ma2);
	S1=SitesCMSA(ma1); S2=SitesCMSA(ma2);
	nblks=config[0]; 
	NEWP(null,nblks+2,BooLean); NEW(len_elem,nblks+2,Int4);
	NEW(model,nblks+2,fm_type); 
	for(j=1; j<= nblks; j++){ 
		if((t=config[j]) < 0){ M = model2[-t]; }
		else { M = model1[t]; }
		model[j] = M;
		NEW(null[j],LenFModel(model[j])+5,BooLean);
	    	NullSitesFModel(null[j], model[j]);
		len_elem[j] = LenFModel(model[j]);
	}
// std::cerr << "debug8\n\n";
	S = CreateNewSites(nblks,len_elem, S1);
	N = NSeqsSites(S);
        for(n=1; n <= N; n++){
                n1 = PosSites(n, tmp1, S1); 
                n2 = PosSites(n, tmp2, S2); 
		for(j=1; j <=nblks; j++){
		   if((t=config[j]) < 0){ s=tmp2[-t]; }
		   else { s=tmp1[t]; }
		   AddSite(j,n,s,S);
		}
	}
// std::cerr << "debug9\n\n";
	ma = MakeCMSA(S, null);
        for(t=1; t <= nblks; t++) free(null[t]); free(null);
	free(model); free(len_elem); 
        return ma;
}

cma_typ	IntersectionCMSA(cma_typ ma1, cma_typ ma2)
/*************************************************************
 returns the 'intersection' of ma1 and ma2 
 This is defined as an alignment where:
   1. if > 50% of the segments in ma1  
  also recombines sequences between the intersection of ma1 and ma2.

 Do this quick and dirty for now; will check & improve later if it works.
 WARNING: 
 ***************************************************************/
{
	cma_typ ma1n,ma2n;
	char	name[15]="merge_junk",nameX[15]="merge_junk.msa";
	char	*name1,*name2;
	Int4	*site_pos;
	sma_typ	sma1,sma2,sma1n,sma2n;
	FILE	*fp=NULL;

fprintf(stderr,"WARNING: need to modify DiffSMA for this to work!\n");
fprintf(stderr,"Only works now if the same number of blocks in both alignments!\n");
	PutAlnCMSA(name, ma1,NULL); sma1 = ReadSMA(nameX);
	PutAlnCMSA(name, ma2,NULL); sma2 = ReadSMA(nameX);

	if(!DiffSMA(fp,sma1,sma2)) return NULL; 
	sma1n = ReadSMA("junk_diffmsa.msa");

	if(!DiffSMA(fp,sma2,sma1)) return NULL; 
	sma2n = ReadSMA("junk_diffmsa.msa");

	name1=NameCMSA(ma1); name2=NameCMSA(ma2);
	ma1n = SMA2CMSA(name1, sma1n);
	ma2n = SMA2CMSA(name2, sma2n);

{	/** Merge Intersecting alignments **/
	Int4	ntyp,n,t,s2;
	ss_type	data=DataCMSA(ma1);
	st_type	S1,S2;
	Int4	N=NumSeqsCMSA(ma1);

	ntyp = nBlksCMSA(ma1n);
	S1=SitesCMSA(ma1n);
	S2=SitesCMSA(ma2n);
#if 0
        PutTypeSites(stdout,S1);
        PutTypeSites(stdout,S2);
#endif
	NEW(site_pos,ntyp+5,Int4);
	for(n =1; n <= N; n++){
// ADD routine to select recombinants here
	  Int4	num,i,aln[2][60];
	  Int4	use;

	  num = SeqCrossPtsCMSA(n, S1, S2, aln[0], aln[1]);
    	  PosSites(n, site_pos, S1);
#if 0
          fprintf(stderr,"seq %d: ncross = %d\n",n,num);
          fprintf(stderr,"\n\naln0: ");
          for(i=1;i <= num; i++){ fprintf(stderr,"%2d ",aln[0][i]); }
	  fprintf(stderr,"\naln1: ");
          for(i=1;i <= num; i++){ fprintf(stderr,"%2d ",aln[1][i]); }
	  fprintf(stderr,"\n\n");
#endif

	  aln[0][num+1] = aln[1][num+1] = 0;
	  for(t =1; t <= ntyp; t++) RmSiteCMSA(t,n,site_pos[t],ma1n);
	  use = random_integer(2);
	  for(i=1,t =1; t <= ntyp; i++,t++){
		while(aln[0][i] != aln[1][i]){
		   if(aln[use][i] != 0){
			if(use == 1){
			   s2 = SitePos(t,n,1,S2); AddSiteCMSA(t,n,s2,ma1n); t++;
			} else { AddSiteCMSA(t,n,site_pos[t],ma1n); t++; }
		   } i++;
		   if(t > ntyp) break;
		} 
		while(aln[use][i] > t){  // Fill in blocks between cross points
			if(use == 1){
			   s2 = SitePos(t,n,1,S2); AddSiteCMSA(t,n,s2,ma1n); t++;
			} else { AddSiteCMSA(t,n,site_pos[t],ma1n); t++; }
			if(t > ntyp) break;
		}
		if(t > ntyp) break;
		if(aln[0][i] == 0){  // then both must be zero
		   while(t <= ntyp) {
			if(use == 1){
			   s2 = SitePos(t,n,1,S2); AddSiteCMSA(t,n,s2,ma1n); t++;
			} else { AddSiteCMSA(t,n,site_pos[t],ma1n); t++; }
		   }
		   break;
		}
		use = random_integer(2);
	        if(use == 1){  // use site from second alignment
			 s2 = SitePos(t,n,1,S2); AddSiteCMSA(t,n,s2,ma1n);
	   	} else AddSiteCMSA(t,n,site_pos[t],ma1n);
		if(t == ntyp) break;
		while(aln[use][i] != 0 && (aln[use][i] +1) < aln[use][i+1]){
		// deal with case where aln[use] goes from t to t+2 or more 
		     if(aln[use][i] == ntyp) break;
		     i++; t++;
		     if(use == 1){
			s2 = SitePos(t,n,1,S2); AddSiteCMSA(t,n,s2,ma1n); 
		     } else AddSiteCMSA(t,n,site_pos[t],ma1n); 
		}
	  }
	}
#if 0
        for(n =1; n <= N; n++){		// old routine...
           p1=p2=0.0;
           for(t =1; t <= ntyp; t++){
                p1 += probSMA(n,t,sma1n);
                p2 += probSMA(n,t,sma2n);
           }
           if(p1 < p2){
             for(t =1; t <= ntyp; t++){
                s1 = SitePos(t,n,1,S1);
                RmSiteCMSA(t,n,s1,ma1n);
             }
             for(t =1; t <= ntyp; t++){
                s2 = SitePos(t,n,1,S2);
                AddSiteCMSA(t,n,s2,ma1n);
             }
           } /** else leave ma1 as it is **/
        }
#endif
	free(site_pos);
	NilCMSA(ma2n); 
}
        NilSMA(sma1); NilSMA(sma2);
        NilSMA(sma1n); NilSMA(sma2n);
	return ma1n;
}

Int4	SeqCrossPtsCMSA(Int4 n, st_type S1, st_type S2, Int4 *aln1, Int4 *aln2)
/**********************************************************************
  Find all crossover points between alignment gaps in ma1 & ma2.

  A crossover occurs between t1 and t2 if for all i=sequences in
  the alignment:

	(s2[i] < e1[i] && s1[i] < e2[i])

              s1          e1
	=======-----------==================   aln1
                      X
	===========-----------==============   aln2
                  s2          e2

 Only need crossovers BETWEEN blocks (not at ends).
/**********************************************************************/
{
    Int4	t,t1,t2,s1,s2,e1,e2,i,ncross;
    Int4	n1,n2,N,ntyp;
    Int4	*start1,*start2,*end1,*end2,*cross1,*cross2;
    ss_type	data = SitesSeqSet(S1);
    BooLean	okay;

    ntyp = nTypeSites(S1); 
    if(ntyp != nTypeSites(S2)) print_error("SeqCrossPtsCMSA() input error");
    N = NSeqsSeqSet(data); i = NSeqsSeqSet(SitesSeqSet(S2));
    if(i != N) print_error("SeqCrossPtsCMSA() input error");
    NEW(start1,ntyp+5,Int4); NEW(start2,ntyp+5,Int4);
    NEW(end1,ntyp+5,Int4);  NEW(end2,ntyp+5,Int4);
    NEW(cross1,2*ntyp+5,Int4);  NEW(cross2,2*ntyp+5,Int4);

    /*** Delineate all gaps between sites ***/
    n1=PosSites(n, end1, S1); n2=PosSites(n, end2, S2);
    if(n1 != ntyp || n2 != ntyp) print_error("CrossPtsCMSA() input error");
    end1[n1+1]=end2[n2+1]=SqLenSeqSet(n,data)+1;	
    start1[1] = start2[1] = 0;
    for(i=1,t=2; i<=n1; i++,t++){ start1[t]=end1[t-1] + SiteLen(t-1,S1) - 1; }
    for(i=1,t=2; i<=n2; i++,t++){ start2[t]=end2[t-1] + SiteLen(t-1,S2) - 1; }

    /*** Now check for crossover points ***/
    for(ncross=1,t1=1; t1<=ntyp+1; t1++){
      s1=start1[t1]; e1=end1[t1]; 
      for(t2=1; t2<=ntyp+1; t2++){
    	okay=TRUE;
	s2=start2[t2]; e2=end2[t2];
	if(!(s2<e1 && s1<e2)) okay=FALSE; 
#if 0
	else {
	   	printf("1. %2d(%d): %d..%d\n",n,t1,s1,e1);
	   	printf("2. %2d(%d): %d..%d\n",n,t2,s2,e2);
	   	printf("\tcrosspt = %2d\n",MINIMUM(Int4,e1,e2));
	}
	printf("\n"); 
#endif
	if(okay){ /** store block type following crossover **/
	    cross1[ncross]=t1;	cross2[ncross]=t2; ncross++;
	}
      }
    }
#if 0
    fprintf(stderr,"cross1: ");
    for(i=1; i< ncross; i++) fprintf(stderr,"%2d ",cross1[i]);
    fprintf(stderr,"\ncross2: ");
    for(i=1; i< ncross; i++) fprintf(stderr,"%2d ",cross2[i]);
    fprintf(stderr,"\n\n");
#endif
    free(start1); free(start2); free(end1); free(end2);
    ncross--;
    {		// for explanation & comments see blks4cross_cmsa( );
	Int4	b,c,t1n,t2n;

	for(b=0,c=1; c<ncross; c++){
	  t1=cross1[c]; t1n=cross1[c+1]; t2=cross2[c]; t2n=cross2[c+1];
  	  if(t1 < t1n && t2 < t2n){ b++; aln1[b] = t1; aln2[b] = t2; }
	  else if(t1 == t1n && t2 < t2n){ b++; aln1[b] = 0; aln2[b] = t2; }
	  else if(t1 < t1n && t2 == t2n){ b++; aln1[b] = t1; aln2[b] = 0; }
	  else print_error("blks4cross_cmsa( ) input error");
	}
	free(cross1); free(cross2);
	return b;
    }
}


