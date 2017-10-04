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

void    SetPseudoToMapCMSA(cma_typ cmsa)
// Set pseudocounts to be the same as for CMSA map.
{
	cmsa->mdl->SetPseudo(NSeqsSeqSet(DataCMSA(cmsa)),PSEUDO_CMSA);
#if 0
        Int4    t,N,ntyps;
        N=NSeqsSeqSet(DataCMSA(cmsa));
        ntyps=nBlksCMSA(cmsa);
        for(t=1; t<=ntyps; t++) {
           SetPseudoFModel((double)N*PSEUDO_CMSA,ModelCMSA(t,cmsa));
        }
#endif
}

void    calc_lngamma_cmsa(cma_typ C)
// static variables are guaranteed to be initialized to zero.
// May want to use this here??? (Instead of cmsa array).
{
        Int4    i,b,N;
        a_type  A;
        double  Ps[30],npseudo,*freq;

	N=NumSeqsCMSA(C); A = AlphabetCMSA(C);
        if(C->lngamma == NULL){
           NEWP(C->lngamma, N+2, double);
           for(i=0; i<=N; i++) NEW(C->lngamma[i], nAlpha(A)+2, double);
        }
        ss_type data = TrueDataCMSA(C); freq=tFreqSeqSet(data);
        npseudo = (double)N*PSEUDO_CMSA;
	if(npseudo > 20.0) npseudo = 20.0;
        for(b=0; b<=nAlpha(A); b++){ Ps[b]=npseudo*freq[b]; }
        for(i=0; i<=N; i++){
           for(b=0; b<=nAlpha(A); b++){
             if(freq[b] > 0.){
                C->lngamma[i][b]=lngamma((double)i + Ps[b]);
             } else C->lngamma[i][b] = 0.0;
           }
        }
}

#include "blosum62.h"
// #include "dirichlet.h"

double  DirichletRelMap(cma_typ L) { return DirichletRelMap('b',L); }

double  DirichletRelMap(char method, cma_typ L)
/**********************************************************************
 return the relative map for aligment given the number of blocks.
 note: this version is independent of the model pseudocounts.
 (it always uses PSEUDO_CMSA*N where N is the number of sequences.
 **********************************************************************/
{
    Int4        N,n,l,len,i,j,t,T,s,ncol,off_col,on_col;
    Int4 	totsite,b,x,*counts,nonsite[30],*site,**observed;
    double      MAP,weight,total,Ps[30],npseudo,*freq,v;
    ss_type	data = DataCMSA(L);
     Int4	nblks=nTypeSites(SitesCMSA(L));
    fm_type	*model=ModelsCMSA(L);
    a_type	A;
    double	b62sites[30],b62nonsites[30];
    double	wt,frq,Nps,indel_penalty=0.0;

    indel_penalty=IndelPenaltySeqSet(data);

    A = SeqSetA(data); counts = CntsSeqSet(data);
    freq=tFreqSeqSet(data); N = NSeqsSeqSet(data);
    for(len=0, t=1; t <= nblks; t++) {
		n = MaxLenFModel(model[t]);
		if(n > len) len = n;
    } NEWP(observed,len+1,Int4);
    /******** compute weight for number of blocks ********/
    for(weight=0.0, s=1; s<=N; s++){
	for(T=0, t=1; t <= nblks; t++) T -= LenFModel(model[t]);
	T += SqLenSeqSet(s,data) + nblks; weight += lnbico(T, nblks);
    }
    /******** end compute weight for number of blocks ********/
    /******** end NEW compute weight for number of blocks ********/
    npseudo = (double)N*PSEUDO_CMSA;
    if(npseudo > 20.0) npseudo = 20.0;
    for(b=0;b<=nAlpha(A);b++){ Ps[b]=npseudo*freq[b]; nonsite[b]=counts[b]; }
    for(off_col=on_col=0,MAP=0.0,t=1; t <= nblks; t++) {
	len = ObservedFModel(observed, model[t]);
        ncol = nColsFModel(model[t]); on_col += ncol;
	off_col += LenFModel(model[t]) - ncol;
	totsite = TotSitesFModel(model[t]);
	/*************** NEW: column width weight ******************/
	MAP -= lnbico(len - 2, ncol - 2);
	/***********************************************************/
	/***** MAP *******/
	for(i=1; i <= len; i++){
	  site = observed[i];
	  if(site != NULL){
	    for(total=0.0,Ps[0]=0.0,b=1;b<= nAlpha(A); b++){ Ps[b]=0.0; total+=site[b]; }
	    if(method == 'b'){
	      for(b=1;b<= nAlpha(A); b++) {
		frq=npseudo*((double)site[b]/(double)total);
	        for(x=1;x<= nAlpha(A); x++){ Ps[x]+=frq*blosum62P[b][x]; }
	      }
	    } else {    // Dirichlet mixture priors...
	      double Prb[50],tPrb=0.0,p,mix,con,v,a;
	      // for(x=1;x<=30; x++){ Prb[x]=0.0; }
	      for(x=1;x<=30; x++){
		mix=Dirichlet30Mix[x]; con=Dirichlet30Conc[x];
		v = lgamma(con) + lgamma(total+1.0) - lgamma(total+con);
	        for(b=1;b <= nAlpha(A); b++) {
		    a=con*Dirichlet30P[x][b]; 
		    v+=lgamma((double)site[b]+a) - lgamma((double)site[b] + 1) - lgamma(a);
		} Prb[x]=mix*exp(v); tPrb+=Prb[x];
	      }
	      for(x=1;x<=30; x++){ Prb[x]=Prb[x]/tPrb; }
	      // fprintf(stderr,"%d: ",i);
	      for(b=1;b <= nAlpha(A); b++) {
	         for(x=1;x<=30; x++){ con=Dirichlet30Conc[x]; Ps[b]+=Prb[x]*con*Dirichlet30P[x][b]; }
	         // fprintf(stderr,"%c%.2f ",AlphaChar(b,A),Ps[b]);
	      } // fprintf(stderr,"\n");
	    }
    	    for(Nps=0.0,b=1;b<= nAlpha(A); b++) Nps+=Ps[b];
	    double d;
	    fprintf(stderr,"%d_%d:",t,i); 
	    for(d=0,x=1;x<= nAlpha(A); x++) {
		d+=Ps[x];
		fprintf(stderr," %c=%.2f",AlphaChar(x,A),Ps[x]); 
		if(x==nAlpha(A)) fprintf(stderr," (total=%.2f)\n",d);
		else if(x%10==0) fprintf(stderr,"\n   ");
	    } if(i==len) fprintf(stderr,"\n");
	    for(b=1;b<= nAlpha(A); b++) {
// fprintf(stderr,"%d_%d.rolling map=%.3f; Ps[%d]=%g; site=%d\n",t,i,MAP,b,Ps[b],site[b]);
		MAP += lngamma(site[b] + Ps[b]); nonsite[b] -= site[b];
	    }
	    // MAP -= lngamma((double)totsite + npseudo);
            // MAP -= lngamma((double)(totsite-site[0]) + npseudo);
            MAP -= lngamma((double)(totsite-site[0]) + Nps);
	  }
// fprintf(stderr,"%d_%d.rolling map=%.3f\n",t,i,MAP);
	}
    }
    /***** NONSITES MAP *******/
    for(Ps[0]=0.0,total=0.0,b=1;b<= nAlpha(A); b++) { Ps[b]=0.0; total += nonsite[b]; }
    if(method == 'b'){
      v = lngamma(npseudo);
      for(b=1;b<= nAlpha(A); b++){
	for(x=1;x<= nAlpha(A); x++) {
	    Ps[x]+=npseudo*((double)nonsite[b]/total)*blosum62P[b][x];
	} 
      }
    } else  {	// Dirichlet mixture priors...
	double p,mix,con,v,a,Prb[50],tPrb=0.0;
	// for(x=1;x<=30; x++){ Prb[x]=0.0; }
	for(x=1;x<=30; x++){
		mix=Dirichlet30Mix[x]; con=Dirichlet30Conc[x];
		v = lgamma(con) + lgamma(total+1.0) - lgamma(total+con);
	        for(b=1;b <= nAlpha(A); b++) {
		    a=con*Dirichlet30P[x][b]; 
		    v+=lgamma((double)nonsite[b]+a) - lgamma((double)nonsite[b] + 1) - lgamma(a);
		} Prb[x]=mix*exp(v); tPrb+=Prb[x];
	} 
	for(x=1;x<=30; x++){ Prb[x]=Prb[x]/tPrb; }
	// fprintf(stderr,"%d: ",i);
	for(b=1;b <= nAlpha(A); b++) {
	        for(x=1;x<=30; x++){ con=Dirichlet30Conc[x]; Ps[b]+=Prb[x]*con*Dirichlet30P[x][b]; }
	        // fprintf(stderr,"%c%.2f ",AlphaChar(b,A),Ps[b]);
	} // fprintf(stderr,"\n");
    }
    for(Nps=0.0,b=1;b<= nAlpha(A); b++) Nps+=Ps[b];
    for(total=0.0,b=1;b<= nAlpha(A); b++) {
	MAP += lngamma((double) nonsite[b] + Ps[b]);
	total += (double) nonsite[b];
        v -= lngamma(Ps[b]);   /** term pulling down MAP **/
    }
    MAP += v *(double) on_col; 
    // MAP -= lngamma((double) total + npseudo);
    MAP -= lngamma((double) total + Nps);

    /***** subtract NULL MAP ******/
    for(b=0;b<= nAlpha(A); b++) Ps[b]=0.0;
    if(method == 'b'){
      for(b=1;b<= nAlpha(A); b++) {
	for(x=1;x<= nAlpha(A); x++) {
	    Ps[x]+=npseudo*freq[b]*blosum62P[b][x];
	}
      }
    } else { 	// Dirichlet mixture priors...
        for(total=0.0,Ps[0]=0.0,b=1;b<= nAlpha(A); b++){ Ps[b]=0.0; total+=counts[b]; }
	double p,mix,con,v,a,Prb[50],tPrb=0.0;
	for(x=1;x<=30; x++){
		mix=Dirichlet30Mix[x]; con=Dirichlet30Conc[x];
		v = lgamma(con) + lgamma(total+1.0) - lgamma(total+con);
	        for(b=1;b <= nAlpha(A); b++) {
		    a=con*Dirichlet30P[x][b]; 
		    v+=lgamma((double)counts[b]+a) - lgamma((double)counts[b] + 1) - lgamma(a);
		} 
		Prb[x]=mix*exp(v); tPrb+=Prb[x];
	} 
	for(x=1;x<=30; x++){ Prb[x]=Prb[x]/tPrb; }
	// fprintf(stderr,"%d: ",i);
	for(b=1;b <= nAlpha(A); b++) {
	        for(x=1;x<=30; x++){ con=Dirichlet30Conc[x]; Ps[b]+=Prb[x]*con*Dirichlet30P[x][b]; }
	        // fprintf(stderr,"%c%.2f ",AlphaChar(b,A),Ps[b]);
	} // fprintf(stderr,"\n");
    }
    for(Nps=0.0,b=1;b<= nAlpha(A); b++) Nps+=Ps[b];
    for(total=0.0, b=1;b<= nAlpha(A); b++) {
	   if(counts[b] > 0.0){
		MAP -= lngamma(counts[b] + Ps[b]);
		total += counts[b];
	   }
    }
    // MAP += lngamma((double)total + npseudo);
    MAP += lngamma((double)total + Nps);
    /**************** weight for column transfers *************/
    if(nblks > 1){
	weight +=  lnbico(on_col-nblks-1,nblks-1);
	if(off_col > 0) weight += lnbico(off_col+nblks-1,nblks-1);
    }
    /**************** return MAP *************/
    free(observed);
    // return (MAP - weight);
    return (MAP - weight - indel_penalty);
}

double  RelMapCMSA2(cma_typ cma)
{
	Int4	*SeqLen = new Int4 [NumSeqsCMSA(cma) +2];
	// ss_type	data = TrueDataCMSA(cma);
	ss_type	data = DataCMSA(cma);

	for(Int4 s=1; s<=NumSeqsCMSA(cma); s++) SeqLen[s]=SqLenSeqSet(s,data);
	double lpr=cma->mdl->LogProbRatio(NumSeqsCMSA(cma),SeqLen,
				IndelPenaltySeqSet(DataCMSA(cma)));
	delete [ ] SeqLen;
	return lpr;
}

double  RelMapCMSA(cma_typ cma)
{
    ss_type	data = DataCMSA(cma);
    double	indel_penalty=0.0;
    double map = UnGappedRelMapCMSA(cma);
    indel_penalty=IndelPenaltySeqSet(data);
    return (map - indel_penalty);
}

double  UnGappedRelMapCMSA(cma_typ cma)
/**********************************************************************
 return the relative map for aligment given the number of blocks.
 note: this version is independent of the model pseudocounts.
 (it always uses PSEUDO_CMSA*N where N is the number of sequences.
 **********************************************************************/
{
    Int4        N,n,len,i,j,t,T,s,ncol,off_col,on_col;
    Int4 	totsite,b,*counts,nonsite[30],*site,**observed;
    double      MAP,weight,total,Ps[30],npseudo,*freq,v;
    ss_type	data = DataCMSA(cma);
    fm_type	*model = ModelsCMSA(cma);    
    Int4	nblks=nBlksCMSA(cma);
    a_type	A;
#if 0 // TEST NEW BLOCK WEIGHTS.
    st_type	S;
    double	weight2;
#endif
    double	indel_penalty=0.0;
    ss_type	truedata = TrueDataCMSA(cma);

    // indel_penalty=IndelPenaltySeqSet(data);

    A = SeqSetA(data); N = NSeqsSeqSet(data);
    if(cma->FullSeq){	// then use full_counts and freq;
	counts = CntsSeqSet(cma->FullSeq); freq=tFreqSeqSet(cma->FullSeq);
    } else { // OLD: prior to full_counts.
       // counts = CntsSeqSet(data); freq=tFreqSeqSet(data); 
       counts = CntsSeqSet(truedata); freq=tFreqSeqSet(truedata); 
    	// freq=tFreqSeqSet(data); 
    }
    for(len=0, t=1; t <= nblks; t++) {
	n = MaxLenFModel(model[t]); if(n > len) len = n;
    }
    NEWP(observed,len+1,Int4);
    /******** compute weight for number of blocks ********/
#if 1	// test domain sampling 
    st_type	S=SitesCMSA(cma);
#endif
    if(cma->FullSeq){ 	// NEW: full_counts.
#if 0	//OLD
      assert(nBlksCMSA(cma) == 1); // use for one block only for now.
      // determine total number of ways to distribute sites ignoring gaps.
      for(weight=0.0, s=1; s<=NSeqsSeqSet(cma->FullSeq); s++){
	if(cma->FullRpts[s] > 0){
	  n = cma->FullRpts[s];
	  T = SqLenSeqSet(s,cma->FullSeq) + n;
	  T -= n*LenFModel(model[1]);
	  if(T < n) T = n;  // allows for gaps...
	  weight += lnbico(T,n);  // number of positions...
	}
      }
#endif
#if 1 	// NEW
      for(weight=0.0, s=1; s<=NSeqsSeqSet(cma->FullSeq); s++){
	if(cma->FullRpts[s] > 0){
	  n = cma->FullRpts[s]*nBlksCMSA(cma);
	  T = SqLenSeqSet(s,cma->FullSeq) + n;
	  T -= n*TotalLenCMSA(cma);
	  if(T < n) T = n;  // allows for gaps...
	  weight += lnbico(T,n);  // number of positions...
	}
      }
#endif
    } else {
      for(weight=0.0, s=1; s<=N; s++){
#if 1	// test domain sampling 
	if(nBlksCMSA(cma) ==1 && nSites(1,s,S) == 0) continue;
#endif
	T = SqLenSeqSet(s,data) + nblks;
	for(t=1; t <= nblks; t++) T -= LenFModel(model[t]);
	if(T < nblks) T = nblks;  // allows for gaps...
	weight += lnbico(T, nblks);
      }
    }
    /******** end compute weight for number of blocks ********/
#if 0
	/** number of places to put the alignment aInt4 the seq. **/
	/** times number of ways to put the blocks within the alignment **/
    /*********** TEST NEW compute weight for number of blocks *********
    S = SitesCMSA(cma);
    for(weight2=0.0, s=1; s<=N; s++){
	i = SitePos(1,s,1,S);
	j = SitePos(nblks,s,nblks,S) + LenFModel(model[nblks]);
	T = SqLenSeqSet(s,data) - (j - i) + 1;
	weight2 += log((double) T);
	if(nblks > 2){
	  for(t=1; t <= nblks; t++) T -= LenFModel(model[t]);
	  b = nblks - 2;
	  weight2 += lnbico(T+b, b);
	}
    }
    /******** end NEW compute weight for number of blocks ********/
#endif
    npseudo = (double)N*PSEUDO_CMSA;
    if(npseudo > 20.0) npseudo = 20.0;
    for(b=1;b<=nAlpha(A);b++){nonsite[b]=counts[b];Ps[b]=npseudo*freq[b];}
    for(off_col=on_col=0,MAP=0.0,t=1; t <= nblks; t++) {
	len = ObservedFModel(observed, model[t]);
        ncol = nColsFModel(model[t]); on_col += ncol;
	off_col += LenFModel(model[t]) - ncol;
	totsite = TotSitesFModel(model[t]);
	/*************** NEW: column width weight ******************/
	MAP -= lnbico(len - 2, ncol - 2);
	/***********************************************************/
	/***** MAP *******/
	for(i=1; i <= len; i++) {
	  site = observed[i];
	  if(site != NULL){
	    for(b=1;b<= nAlpha(A); b++) {
		MAP += cma->lngamma[site[b]][b]; 
		nonsite[b] -= site[b];
#if 0		/***** OLD *****/
	      if(Ps[b] > 0){
		MAP += lngamma((double)site[b] + Ps[b]); 
		nonsite[b] -= site[b];
	      } /***** OLD *****/
#endif
	    }
	    MAP -= lngamma((double)(totsite-site[0]) + npseudo);
	    // MAP -= lngamma((double)totsite + npseudo);
	  } 
	}
    }
    v = lngamma(npseudo);
    /***** NONSITES MAP *******/
if(cma->FullSeq){	// check for problems due to flanking regions...
    for(b=1;b<= nAlpha(A); b++) assert(nonsite[b] >= 0);
}
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
    if(nblks > 1){
	weight +=  lnbico(on_col-nblks-1,nblks-1);
	if(off_col > 0) weight += lnbico(off_col+nblks-1,nblks-1);
    }
    /**************** return MAP *************/
    free(observed);
    return (MAP - weight);
    // return (MAP - weight - indel_penalty);
}

double  single_rel_map_cmsa(cma_typ L, Int4 t)
/**********************************************************************
 return the relative map for aligment assuming a single block.
 note: this version is independent of the model pseudocounts.
 (it always uses PSEUDO_CMSA*N where N is the number of sequences.
 WARNING: does not adjust for indel penalties!!!!
 **********************************************************************/
{
    Int4        N,len,i,T,s,ncol;
    Int4 	totsite,b,*counts,nonsite[30],*site,**observed;
    double      MAP,weight,total,Ps[30],npseudo,*freq,v;
    ss_type	data = DataCMSA(L);
    fm_type	model;
    a_type	A;

    model=ModelCMSA(t,L);
    A = SeqSetA(data); counts = CntsSeqSet(data);
    freq=tFreqSeqSet(data); N = NSeqsSeqSet(data);
    len = MaxLenFModel(model);
    NEWP(observed,len+1,Int4);
    for(weight=0.0, s=1; s<=N; s++){
	T = SqLenSeqSet(s,data) + 1;
	T -= LenFModel(model);
	weight += lnbico(T, 1);
    }
    npseudo = (double)N*PSEUDO_CMSA;
    if(npseudo > 20.0) npseudo = 20.0;
    for(b=1;b<=nAlpha(A);b++){nonsite[b]=counts[b];Ps[b]=npseudo*freq[b];}
    len = ObservedFModel(observed, model);
    ncol = nColsFModel(model); 
    totsite = TotSitesFModel(model);
    MAP = -lnbico(len - 2, ncol - 2); 	/** column width weight **/
    for(i=1; i <= len; i++) {
	site = observed[i];
	if(site != NULL){
	  for(b=1;b<= nAlpha(A); b++) {
		MAP += L->lngamma[site[b]][b]; 
		nonsite[b] -= site[b];
#if 0		/***** OLD *****/
	      if(Ps[b] > 0){
		MAP += lngamma((double)site[b] + Ps[b]); 
		nonsite[b] -= site[b];
	      }
#endif 
	  }
	  // MAP -= lngamma((double)totsite + npseudo);
	  MAP -= lngamma((double)(totsite-site[0]) + npseudo);
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
    MAP += v *(double) ncol; 
    MAP -= lngamma((double) total + npseudo);
    /***** subtract NULL MAP ******/
    for(total = 0.0, b=1;b<= nAlpha(A); b++) {
	   if(counts[b] > 0.0){
		MAP -= lngamma(counts[b] + Ps[b]);
		total += counts[b];
	   }
    }
    MAP += lngamma((double)total + npseudo);
    free(observed);
    return (MAP - weight);
}

double  RelMapMinusCMSA(Int4 tdel, cma_typ L)
/**********************************************************************
 Return the relative map for aligment but with the t_rm block removed.
 note: this version is independent of the model pseudocounts.
 (it always uses PSEUDO_CMSA*N where N is the number of sequences.
 NOTE: this could be used for general purpose by setting tdel=0;
 **********************************************************************/
{
    Int4        N,n,len,i,t,T,s,ncol,off_col,on_col,nblks;
    Int4 	totsite,b,*counts,nonsite[30],*site,**observed;
    double      MAP,weight,total,Ps[30],npseudo,*freq,v;
    fm_type	*model=ModelsCMSA(L);
    ss_type	data = DataCMSA(L);
    st_type	S=SitesCMSA(L);
    a_type	A;
    double	indel_penalty=0.0;

    indel_penalty=IndelPenaltySeqSet(data);

    nblks=nTypeSites(S);
    A = SeqSetA(data); counts = CntsSeqSet(data);
    freq=tFreqSeqSet(data); N = NSeqsSeqSet(data);
    for(len=0, t=1; t <= nblks; t++) {
	if(t != tdel){
		n = MaxLenFModel(model[t]);
		if(n > len) len = n;
	}
    }
    NEWP(observed,len+1,Int4);
    for(weight=0.0, s=1; s<=N; s++){
	for(T=0, t=1; t <= nblks; t++){
	   if(t != tdel) T -= LenFModel(model[t]);
		
	}
	T += SqLenSeqSet(s,data) + nblks -1;
	weight += lnbico(T, (nblks-1));
    }
    npseudo = (double)N*PSEUDO_CMSA;
    if(npseudo > 20.0) npseudo = 20.0;
    for(b=1;b<=nAlpha(A);b++){nonsite[b]=counts[b];Ps[b]=npseudo*freq[b];}
    for(off_col=on_col=0,MAP=0.0,t=1; t <= nblks; t++) {
      if(t != tdel){
	len = ObservedFModel(observed, model[t]);
        ncol = nColsFModel(model[t]); on_col += ncol;
	off_col += LenFModel(model[t]) - ncol;
	totsite = TotSitesFModel(model[t]);
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
    if(nblks > 2){
	weight +=  lnbico(on_col-nblks-2,nblks-2);
	if(off_col > 0) weight += lnbico(off_col+nblks-2,nblks-2);
    }
    /**************** return MAP *************/
    free(observed);
    // return (MAP - weight);
    return (MAP - weight - indel_penalty);
}

/***************************** TEST ********************************/
double  PutRelMapCMSA(FILE *fp, cma_typ L)
{
    st_type	S = SitesCMSA(L);
    Int4	nblks = nTypeSites(S);
    fm_type 	*model = ModelsCMSA(L);
    Int4        N,n,len,i,k,t,T,s,ncol,off_col,on_col;
    Int4 	totsite,b,*counts,nonsite[30],*site,**observed;
    double      d,MAP,weight,total,Ps[30],npseudo,*freq,v,map,map2;
    ss_type	data = DataCMSA(L);
    a_type	A;
    double	indel_penalty=0.0;

    indel_penalty=IndelPenaltySeqSet(data);

    A = SeqSetA(data); counts = CntsSeqSet(data);
    freq=tFreqSeqSet(data); N = NSeqsSeqSet(data);
    for(len=0, t=1; t <= nblks; t++) {
		n = MaxLenFModel(model[t]);
		if(n > len) len = n;
    }
    NEWP(observed,len+1,Int4);
    /**************** MAP ******************/
    npseudo = (double)N*PSEUDO_CMSA;
    if(npseudo > 20.0) npseudo = 20.0;
    for(b=1;b<=nAlpha(A);b++){nonsite[b]=counts[b];Ps[b]=npseudo*freq[b];}
MAP = -indel_penalty;
// MAP=0.0;
    for(off_col=on_col=0,t=1; t <= nblks; t++) {
	len = ObservedFModel(observed, model[t]);
        ncol = nColsFModel(model[t]); on_col += ncol;
	off_col += LenFModel(model[t]) - ncol;
	totsite = TotSitesFModel(model[t]);
	/*************** NEW: column width weight ******************/
	d = -lnbico(len - 2, ncol - 2);
	fprintf(fp,"model %d fragmentation penalty: %g\n",t,d);
	MAP += d;
	/***********************************************************/
	/***** MAP *******/
	for(i=1; i <= len; i++) {
	  site = observed[i];
	  if(site != NULL){
	    for(d=0.,b=1;b<= nAlpha(A); b++) {
	      if(Ps[b] > 0){
		d += lngamma((double)site[b] + Ps[b]); 
        	d -= lngamma(Ps[b]);   /** term pulling down MAP **/
		nonsite[b] -= site[b];
	      }
	    }
    	    d += lngamma(npseudo);
	    // d -= lngamma((double)totsite + npseudo);
            d -= lngamma((double)(totsite-site[0]) + npseudo);
	    fprintf(fp,"model %d column %2d: %g\n",t,i,d);
	    MAP += d;
	  } 
	}
    }
    fprintf(fp,"Raw MAP: %g\n",MAP); 
    /***** NONSITES MAP *******/
    for(d=0.,total=0.0, b=1;b<= nAlpha(A); b++) {
	   if(Ps[b] > 0){
		d += lngamma((double) nonsite[b] + Ps[b]);
        	d -= lngamma(Ps[b]);   /** term pulling down MAP **/
		total += (double) nonsite[b];
	   }
    }
    d += lngamma(npseudo);
    d -= lngamma((double) total + npseudo);
    fprintf(fp,"nonsites: %g\n",d); 
    MAP += d;
    fprintf(fp,"Running MAP: %g\n",MAP); 
    /***** subtract NULL MAP ******/
    for(v=0.,total=0.0, b=1;b<= nAlpha(A); b++) {
	   if(counts[b] > 0.0){
		v -= lngamma(counts[b] + Ps[b]);
        	v += lngamma(Ps[b]);   /** term pulling down MAP **/
		total += counts[b];
	   }
    }
    v -= lngamma(npseudo);
    v += lngamma((double)total + npseudo);
    fprintf(fp,"(null map deficit per column: %g)(%d columns)\n",
		-(d+v)/(double)on_col,on_col); 
    fprintf(fp,"null map adjustment: %g\n",v); 
    MAP += v;
    fprintf(fp,"Running MAP: %g\n",MAP); 
    /**************** weight for blocks *************/
    for(weight=0.0, s=1; s<=N; s++){
	for(T=0, t=1; t <= nblks; t++) T -= LenFModel(model[t]);
	T += SqLenSeqSet(s,data) + nblks;
	weight -= lnbico(T, nblks);
    }
    fprintf(fp,"block configuration penalty = %g\n",weight);
    /**************** weight for column transfers *************/
    /** at least two columns required for each block **/
    d = 0.;
    if(nblks > 1){
	d -=  lnbico(on_col-nblks-1,nblks-1);
	if(off_col > 0) d -= lnbico(off_col+nblks-1,nblks-1);
    }
    fprintf(fp,"weight for column transfers: %g\n",d);
    weight += d;
    /**************** return MAP *************/
    MAP += weight; free(observed);
    fprintf(fp,"Net MAP = %g\n",MAP);

    k = nTypeSites(S);
    for(t=1; t <= k; t++){
	map = single_rel_map_cmsa(L, t);
	map2=RelMapMinusCMSA(t, L);
	fprintf(fp,"\n(%d)\tsingle block NetMAP = %.2f\n", t,map);
	fprintf(fp,"   \tNetMAP with block deleted = %.2f\n", map2);
	map2=FieldRelMapCMSA(L, t);
	fprintf(fp,"   \tField NetMAP = %.2f\n", map2);

    }
    return MAP;
}

double  FieldRelMapCMSA(cma_typ L, Int4 t)
/**********************************************************************
 Return the relative map for aligment assuming that it was found
 within its field.
 Note: this version is independent of the model pseudocounts.
 (it always uses PSEUDO_CMSA*N where N is the number of sequences.
 WARNING: does not adjust for indel penalties!!!!
 **********************************************************************/
{
    Int4        N,len,i,T,s,ncol,*fieldlen;
    Int4 	totsite,b,*counts,nonsite[30],*site,**observed;
    double      MAP,weight,total,Ps[30],npseudo,*freq,v;
    ss_type	data = DataCMSA(L);
    fm_type	model;
    a_type	A;

    A = SeqSetA(data); N = NSeqsSeqSet(data);
    model = ModelCMSA(t,L);
    NEW(fieldlen,N+2,Int4);
    NEW(freq,nAlpha(A)+2,double);
    counts = CntsFieldCMSA(L, t, fieldlen);
    for(total=0.0,i=0; i<= nAlpha(A); i++){ total += (double)counts[i]; }
    for(i=0; i<= nAlpha(A); i++){ freq[i] = (double)counts[i]/total; }

    len = MaxLenFModel(model);
    NEWP(observed,len+1,Int4);
    for(weight=0.0, s=1; s<=N; s++){
	T = fieldlen[s] - LenFModel(model) + 1;
	weight += lnbico(T, 1);
    }
    npseudo = (double)N*PSEUDO_CMSA;
    if(npseudo > 20.0) npseudo = 20.0;
    for(b=1;b<=nAlpha(A);b++){nonsite[b]=counts[b];Ps[b]=npseudo*freq[b];}
    len = ObservedFModel(observed, model);
    ncol = nColsFModel(model); 
    totsite = TotSitesFModel(model);
    MAP = -lnbico(len - 2, ncol - 2); 	/** column width weight **/
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
    v = lngamma(npseudo);
    /***** NONSITES MAP *******/
    for(total = 0.0, b=1;b<= nAlpha(A); b++) {
	   if(Ps[b] > 0){
		MAP += lngamma((double) nonsite[b] + Ps[b]);
		total += (double) nonsite[b];
        	v -= lngamma(Ps[b]);   /** term pulling down MAP **/
	   }
    }
    MAP += v *(double) ncol; 
    MAP -= lngamma((double) total + npseudo);
    /***** subtract NULL MAP ******/
    for(total = 0.0, b=1;b<= nAlpha(A); b++) {
	   if(counts[b] > 0.0){
		MAP -= lngamma(counts[b] + Ps[b]);
		total += counts[b];
	   }
    }
    MAP += lngamma((double)total + npseudo);
    free(observed); free(fieldlen); free(freq); free(counts);
    return (MAP - weight);
}

double	RecombinantMapCMSA(cma_typ ma1, cma_typ ma2, Int4 *config)
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
    ss_type	data = DataCMSA(ma1);
    st_type	S1,S2;
    fm_type	M,*model1,*model2;
    a_type	A;
    double      indel_penalty=0.0;

    indel_penalty=IndelPenaltySeqSet(data);

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
	    // MAP -= lngamma((double)totsite + npseudo);
            MAP -= lngamma((double)(totsite-site[0]) + npseudo);
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

double  RelMapMinusSeqCMSA(Int4 sq, cma_typ cma)
/**********************************************************************
 return the relative map for aligment without sequence sq.
 **********************************************************************/
{
    Int4        N,n,len,i,j,t,T,s,ncol,off_col,on_col;
    Int4 	totsite,b,*counts,nonsite[30],*site,**observed;
    double      MAP,weight,total,Ps[30],npseudo,*freq,v;
    Int4	nblks=nBlksCMSA(cma);
    a_type	A;
    double	indel_penalty=0.0;
    ss_type	truedata = TrueDataCMSA(cma);

    A = SeqSetA(truedata); N = NSeqsSeqSet(truedata);
    counts = CntsSeqSet(truedata); freq=tFreqSeqSet(truedata);
    
    // 1. Remove sequence sq from alignment.
    assert(sq > 0 && sq <= N);
    Int4 *oldpos; NEW(oldpos,nBlksCMSA(cma)+2,Int4);  // save old sites.
    for(t=1; t<=nBlksCMSA(cma); t++){
            PosSiteCMSA(t,sq,cma->pos,cma); oldpos[t]=cma->pos[1];
    } VacateSitesCMSA(sq,cma);

    fm_type	*model = ModelsCMSA(cma);    
    ss_type	data = DataCMSA(cma);
    indel_penalty=IndelPenaltySeqSet(data);
#if 1
    gss_typ& gss=*gssCMSA(cma);
    indel_penalty-=gss.GapOpen()*gss.NumOpen(sq);
    indel_penalty-=gss.GapExtend()*gss.NumExtend(sq);
#endif
    for(len=0, t=1; t <= nblks; t++) {
	n = MaxLenFModel(model[t]); if(n > len) len = n;
    }
    NEWP(observed,len+1,Int4);
    /******** compute weight for number of blocks ********/
    for(weight=0.0, s=1; s<=N; s++){
	if(s != sq){
	  T = SqLenSeqSet(s,data) + nblks;
	  for(t=1; t <= nblks; t++) T -= LenFModel(model[t]);
	  weight += lnbico(T, nblks);
	}
    }
    /******** end compute weight for number of blocks ********/
    npseudo = (double)N*PSEUDO_CMSA;
    if(npseudo > 20.0) npseudo = 20.0;
    for(b=1;b<=nAlpha(A);b++){nonsite[b]=counts[b];Ps[b]=npseudo*freq[b];}
    for(off_col=on_col=0,MAP=0.0,t=1; t <= nblks; t++) {
	len = ObservedFModel(observed, model[t]);
        ncol = nColsFModel(model[t]); on_col += ncol;
	off_col += LenFModel(model[t]) - ncol;
	totsite = TotSitesFModel(model[t]);
	MAP -= lnbico(len - 2, ncol - 2); // NEW: column width weight 
	for(i=1; i <= len; i++) {
	  site = observed[i];
	  if(site != NULL){
	    for(b=1;b<= nAlpha(A); b++) {
		MAP += cma->lngamma[site[b]][b]; 
		nonsite[b] -= site[b];
	    }
	    MAP -= lngamma((double)(totsite-site[0]) + npseudo);
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
    if(nblks > 1){	// weight for column transfers 
	weight +=  lnbico(on_col-nblks-1,nblks-1);
	if(off_col > 0) weight += lnbico(off_col+nblks-1,nblks-1);
    } free(observed);
    MAP = (MAP - weight - indel_penalty);
    for(t=1; t <= nBlksCMSA(cma); t++) AddSiteCMSA(t,sq,oldpos[t],cma);
    free(oldpos); 
    return MAP;
}

