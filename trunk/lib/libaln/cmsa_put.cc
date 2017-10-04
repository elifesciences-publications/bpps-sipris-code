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

void	PutRptPartionCMSA(FILE *fp, cma_typ cma)
// Put the number of repeats on either side of the midpoint partition.
{
	gss_typ	*gss=gssCMSA(cma);
	e_type	E;
	Int4	J,m,x,*s,pos[5],start,end,repeats;
	Int4	number=gss->NumSeq();
	char	id0[202],id[202];

	if(nBlksCMSA(cma) != 1) print_error("PutRptSpacingCMSA( ) input error");
	Int4	Sq,*NumFake,*NumReal;
	NEW(NumFake,number+3,Int4);
	NEW(NumReal,number+3,Int4);
    // 1. Loop over all sequences in aligment.
	Sq=end=0;
    	for(id0[0]=0,repeats=0,J=1; J <= number; J++) {
	   E = gss->TrueSeq(J);
	   Int4 FullLen = OffSetSeq(E) + LenSeq(E) + CtermExtendSeq(E);
	   Int4 MidLen = (FullLen/2) - (LengthCMSA(1,cma)/2);
	   if(MidLen < 1) MidLen=1;
	   StrSeqID(id,200,E);
	   if(strncmp(id,id0,200) != 0){ // New sequence...
		repeats=0; Sq++; strcpy(id0,id);
	   }
	   PosSiteCMSA(1,J,pos,cma);
	   start = gss->TrueSite(J,pos[1]); // offset is added here...
	   end = pos[1]+LengthCMSA(1,cma)-1; end = gss->TrueSite(J,end);
	   if(start >= MidLen) NumFake[Sq]++; else NumReal[Sq]++;
	   repeats++;
	}
	fprintf(fp," Seq.: Real Fake\n");
	Int4 f=0,r=0;
	for(Int4 i=1; i <= Sq; i++){
	   fprintf(fp," %4d: %4d %4d\n",i,NumReal[i],NumFake[i]);
	   f += NumFake[i]; r += NumReal[i];
	}
	fprintf(fp,"Total: %4d %4d\n",r,f);
	free(NumFake); free(NumReal);
}

void	PutBlockSpacingCMSA(FILE *fp, Int4 block, cma_typ cma)
/*************** Codes used: ************************
  Eventually use true seqs to get accurate estimates...
  ***************************************************/
{
	h_type	HG;
        Int4   n,res;
	char	str[100];
	mh_type	mH=0;

        if(block < 0 || block > nBlksCMSA(cma)) return;
	// eventually use heap to get spacings first...
	sprintf(str,"Spacing after block %d",block);
	HG = Histogram(str,0,100,2.0);
	// 1. Loop over all sequences in aligment.
	Int4	N  = NumSeqsCMSA(cma);
	if(N >= 20){ mH=Mheap(N/4,3); }
        for(res=0,n=1; n <= NumSeqsCMSA(cma); n++){
           res = GapBetweenSites(n,block,SitesCMSA(cma));
	   if(mH) InsertMheap(-res,mH);
	   IncdHist(res,HG);
        } PutHist(stdout,60,HG); NilHist(HG);
	if(mH) {
	  sprintf(str,"Top 25%% largest gaps");
	  HG = Histogram(str,0,100,2.0);
	  while(!EmptyMheap(mH)){
		double d = MinKeyMheap(mH);
		IncdHist(-d,HG);
		DelMinMheap(mH);
	  } NilMheap(mH); PutHist(stdout,60,HG); NilHist(HG);
	}
}

#if 0	// EXPERIMENTAL: want to create cobbled sequences....
void	PutRptSpacingsCMSA(FILE *fp, cma_typ cma)
/*************** Codes used: ************************
Put cobbled sequences with rpts number of repeats and with
spacers taken from the actual sequences...
  ***************************************************/
{
	gss_typ	*gss=gssCMSA(cma);
	e_type	E;
	Int4	J,m,x,*s,pos[5],start,end,repeats;
	Int4	i,j,number=gss->NumSeq();
	char	id0[202],id[202];
	a_type	A=AlphabetCMSA(cma);
	h_type	HG;

	if(nBlksCMSA(cma) != 1) print_error("PutRptSpacingsCMSA( ) input error");
	HG = Histogram("repeat spacing",-100,250,5);
    // 1. Loop over all sequences in aligment.
	end = 0;
    	for(id0[0]=0,repeats=0,J=1; J <= number; J++) {
	   E = gss->TrueSeq(J);
	   PosSiteCMSA(1,J,pos,cma);
	   start = gss->TrueSite(J,pos[1]); // offset is added here...
	   StrSeqID(id,200,E);
	   if(strncmp(id,id0,200) != 0){ // New sequence...
		if(J > 1){ printf("\n"); repeats=0; } strcpy(id0,id);
		PutSeq(stdout,E,A);
		printf(">fake sequence (%s)\n",id);
		end=OffSetSeq(E);
	   }
	   for(i=end+1; 0 && i < start; i++){
		printf("%c",AlphaChar(ResSeq(i-OffSetSeq(E),E),A));
	   } // printf("\n");
	   if(repeats!=0) IncdHist(start-end, HG);
	   if(repeats!=0){
	      for(i=start; i <= end; i++){
		printf("%c",AlphaChar(ResSeq(i-OffSetSeq(E),E),A));
	      } if(i > start) printf("\n");
	   }
	   end = pos[1]+LengthCMSA(1,cma)-1;
	   end = gss->TrueSite(J,end);
	   for(i=start; i <= end; i++){
		printf("%c",AlphaChar(ResSeq(i-OffSetSeq(E),E),A));
	   } if(i > start) printf("\n");
	   repeats++;
	}
	PutHist(stdout,60,HG); NilHist(HG);
}
#endif

void	PutRptSpacingCMSA(FILE *fp, cma_typ cma)
/*************** Codes used: ************************
  ***************************************************/
{
	gss_typ	*gss=gssCMSA(cma);
	e_type	E;
	Int4	J,m,x,*s,pos[5],start,end,repeats;
	Int4	number=gss->NumSeq();
	char	id0[202],id[202];
	h_type	HG;

	if(nBlksCMSA(cma) != 1) print_error("PutRptSpacingCMSA( ) input error");
	HG = Histogram("repeat spacing",-100,250,1);
    // 1. Loop over all sequences in aligment.
	end = 0;
    	for(id0[0]=0,repeats=0,J=1; J <= number; J++) {
	   E = gss->TrueSeq(J);
	   StrSeqID(id,200,E);
	   if(strncmp(id,id0,200) != 0){ // New sequence...
		if(J > 1){ repeats=0; } strcpy(id0,id);
	   }
	   PosSiteCMSA(1,J,pos,cma);
	   start = gss->TrueSite(J,pos[1]); // offset is added here...
	   if(repeats!=0) IncdHist(start-end, HG);
	   end = pos[1]+LengthCMSA(1,cma)-1;
	   end = gss->TrueSite(J,end);
	   repeats++;
	}
	PutHist(stdout,60,HG); NilHist(HG);
}

void    PutModelsCMSA(FILE *fp, cma_typ cma)
{
    Int4 k=nTypeSites(SitesCMSA(cma));
    for(Int4 t=1; t <= k; t++){
        fprintf(fp,"================ Motif %c ================\n",
                        ('A' + (char) t - 1));
        PutFModel(fp,ModelCMSA(t,cma));
    } 
}

void	WriteMtfCMSA(char *name, cma_typ cma, gd_type G)
/*********************************************************************
 Open files with name and extensions and write important cma information. 
 WARNING: this clobbers any files with the same names.
 *********************************************************************/
{
    FILE	*fp;

    fp = open_file(name,".mdl","w");
    PutModelsCMSA(fp,cma); fclose(fp);
    PutAlnCMSA(name, cma,G); 
    WriteConsensusSeqCMSA(name, cma);
}

void	PutInterBlockCMSA(FILE *fptr, Int4 t1, Int4 t2, cma_typ L)
// put a fasta formated sequence of the region between blocks t1 & t2
{ 
    ss_type		data=DataCMSA(L);
    st_type		sites = SitesCMSA(L);
    Int4		n,end,N,ntyps;
    Int4 		start,length;
    e_type		E;
    
    N = NSeqsSeqSet(data); ntyps=nTypeSites(sites);
    if(t1 >= t2 || t2 > ntyps || t1 < 1) print_error("PutInterBlockCMSA( ) input error");
    for(n = 1; n <= N; n++) {
        E=SeqSetE(n,data);
        length = SiteLen(t1,sites);
        start= SitePos(t1,n,1,sites);
        start = start + length - 1;
	end = SitePos(t2,n,1,sites);
	if(start < end) PutSubSeq(fptr, start, end, E, SeqSetA(data));
    }
}

void	PutStuffedSeqSetCMSA(FILE *fptr, Int4 mingap, cma_typ L)
// put a fasta formated file for each sequence filling in the gaps
// between blocks with x's so that at least mingap residues are
// between each motif...
{ 
    ss_type		data=DataCMSA(L);
    st_type		sites = SitesCMSA(L);
    Int4		n,end,N,ntyps,t,s,start,length;
    unsigned char	*seq;
    e_type		E;
    a_type		A;
    
    if(mingap < 1) print_error("PutStuffedSeqSetCMSA( ) input error");
    N = NSeqsSeqSet(data); ntyps=nTypeSites(sites); A = SeqSetA(data);
    for(n = 1; n <= N; n++) {
      E=SeqSetE(n,data);
      fprintf(fptr,">"); PutSeqInfo(fptr,E);
      seq = SeqPtr(E); 
      start = 1; end = SitePos(1,n,1,sites);
      for(t=1; t <= ntyps; t++){
	// print gap...
	s=end-start;
	while(s < mingap){ fprintf(fptr,"x"); s++; }
	while(start < end){ 
	   fprintf(fptr,"%c",AlphaChar(seq[start],A)); start++; 
	}

	// print matching region...
        length = SiteLen(t,sites);
        end = start + length;
	while(start < end) {
	   fprintf(fptr,"%c",AlphaChar(seq[start],A)); start++;
	}

	if(t < ntyps) {
		start = end;
		end = SitePos(t+1,n,1,sites);
	} else if(t == ntyps){
		start = end;
		end = LenSeq(E) + 1;
		s=end-start;
		while(start < end){ 
		   fprintf(fptr,"%c",AlphaChar(seq[start],A)); start++; 
		}
		while(s < mingap){ fprintf(fptr,"x"); s++; }
		fprintf(fptr,"\n\n");
	}
      }
    }
}


void	WritePhylipCMSA(char *name, cma_typ L)
/*************************************************************************
 output a *.phy (phylogenetic analysis) file of the sites 
 WARNING: assumes that "InitMAPMSA(L); " has been invoked!!!
 WARNING: will be right if not COLINEAR!!
 *************************************************************************/
{
	st_type	S=SitesCMSA(L);
	FILE	*fptr;
/**** NEW: ignore non-informative positions ****/
        Int4    t,ntyps=nTypeSites(S);
	fm_type	M;
	BooLean	**null;
/**** NEW: ignore non-informative positions ****/

	fptr = open_file(name,".phy","w");
/**** NEW: ignore non-informative positions ****/
	NEWP(null, ntyps +2, BooLean);
        for(t=1; t <= ntyps; t++) {
	  M=ModelCMSA(t,L);
	  NEW(null[t], LenFModel(M) +2, BooLean);
          NullSitesFModel(null[t], M);
	}
        PutPhylipSites(fptr,S,null);
        for(t=1; t <= ntyps; t++) free(null[t]);
	free(null);
/**** NEW: ignore non-informative positions ****/
	fclose(fptr);
}

void	WriteConsensusSeqCMSA(char *name, cma_typ L)
/*************************************************************************
 output a *.csq file of a consensus sequence for the alignment
 WARNING: assumes that "InitMAPMSA(L); " has been invoked!!!
 WARNING: assumes alignment is COLINEAR!!
 *************************************************************************/
{
	
	ss_type	data=DataCMSA(L);
	st_type	S=SitesCMSA(L);
	Int4	n,s,t,N=NSeqsSeqSet(data),ntyps=nTypeSites(S);
	Int4	*gaplen,r,i,j,*cnts,start,end,d;
	Int4	seqlen;
	unsigned char	*seq;
	double	*freq = tFreqSeqSet(data),rand,likelihood,best;
	fm_type	M;
	a_type	A = SeqSetA(data);
	FILE	*fptr;

    fptr = open_file(name,".csq","w");
    NEW(gaplen, ntyps+3, Int4);
    /*** 1. find average gap lengths ***/
    for(n = 1; n <= N; n++) {
        for(end = 0, t=0; t <= ntyps; t++) {
	  if(t == ntyps) start = SqLenSeqSet(n,data)+1;
	  else start = SitePos(t+1,n,1,S);
	  gaplen[t] += start-(end+1);
	  if(t < ntyps) { end = EndSitePos(t+1,n,1,S); }
	}
    }
    for(seqlen=0,t=0; t <= ntyps; t++) { 
	if(t > 0){ M=ModelCMSA(t,L); seqlen+=LenFModel(M); }
	gaplen[t] = 1 + gaplen[t]/N;  seqlen+=gaplen[t];
    }
    NEW(seq,seqlen+3,unsigned char);
    fprintf(fptr,">maxseq "); 
    /*** 2. generate consensus sequence ***/
    for(i=0, t = 0; t <= ntyps; t++) {
	/** 2.a. fill in motif region ***/
	if(t > 0){
	  M = ModelCMSA(t,L);
	  start = i+1;
	  for(d=0; d < LenFModel(M); d++){
		cnts = GetSiteFreq(S, t, d);
#if 0
		for(best=0.0,j=1; j <= nAlpha(A); j++){
		   likelihood = (double) cnts[j]/(double) N;
		   likelihood /= freq[j];
		   if(best < likelihood){ best=likelihood; r=j; }
		}
#else		// highest frequency concensus
		for(best=0.0,j=1; j <= nAlpha(A); j++){
		   if(best < (double)cnts[j]){ best=(double)cnts[j]; r=j; }
		}
#endif
		/*** fprintf(stderr,"best = %g; r = %c\n",
			best, AlphaChar(r,A)); /*****/
		free(cnts);
	   	i++; seq[i] = r;
	  }
	  end = i;
	  if(t < ntyps) fprintf(fptr,"%d-%d;",start,end);
	  else fprintf(fptr,"%d-%d.",start,end);
	}
	/** 2.b. fill in gap regions ***/
	for(s=1; s <= gaplen[t]; s++){
           rand = (double) Random()/(double) RANDOM_MAX;
           for(r = 0; (rand-=freq[r]) > 0.0; r++){ if(r == nAlpha(A)) break; }
	   i++; seq[i] = r;
	}
    }
    if(seqlen != i) print_error("Consensus length error");
    fprintf(fptr,"\n");
    for(j=1; j <=i; j++){
	r = seq[j];
	fprintf(fptr,"%c",AlphaChar(r,A));
	if(j % 60 == 0) fprintf(fptr,"\n");
    }
    fprintf(fptr,"\n");
    free(seq); free(gaplen);
    fclose(fptr);
}

void    PutSitesCMSA(FILE *fptr,Int4 t,double **site_prob, cma_typ L)
{
        NullSitesFModel(L->null,ModelCMSA(t,L));
        PutSites(fptr,t,L->sites,site_prob,L->null);
        PutFModel(fptr,ModelCMSA(t,L));
}

void    PutMSA(FILE *fptr, cma_typ L)
/*************************************************************************
 output the sites
 WARNING: you may want to invoke "InitMAPMSA(L); " prior to output!!!
 *************************************************************************/
{
        Int4    n,s,t,N=NumSeqsCMSA(L),end;
        fm_type *finalmodel=NULL;
        st_type mapsites=NULL;
        double  **pos_prob;
        char    str[50];
        unsigned char   *seq;
	mdl_typ *mdl=mdlCMSA(L);

    finalmodel = ModelsCMSA(L);
    mapsites = L->sites;
    for(t=1; t <= nTypeSites(mapsites); t++) {
        sprintf(str,"MOTIF %c",t + 'A' -1);
        fprintf(fptr,"\n\n%*s%*s%s",
                23,"",(SiteLen(t,mapsites)-7)/2,"",str);
        /****** compute probabilites *****/
        NEWP(pos_prob,N+2,double);
        for(n = 1; n <= N; n++) {
                end = LenSeqCMSA(n,L);
                NEW(pos_prob[n],end+5,double);
                end -= LenFModel(finalmodel[t]) + 1;
                seq = SeqPtrCMSA(n,L);
                for(pos_prob[n][0]=0.0, s= 1; s<= end; s++){
                   if(TypeSite(n,s,mapsites) == t){
			mdl->Remove(t,seq,s);
                        pos_prob[n][s] = (double)
                           LikelihoodFModel(seq, s, finalmodel[t]);
			mdl->Add(t,seq,s);
                   } else pos_prob[n][s] = (double)
                           LikelihoodFModel(seq, s, finalmodel[t]);
                   pos_prob[n][0] += pos_prob[n][s];
                }
        }
        for(n = 1; n <= N; n++) {
          end = LenSeqCMSA(n,L) - LenFModel(finalmodel[t]) + 1;
          for(s = 1; s <= end; s++) {
            if(pos_prob[n][s]>0.0) pos_prob[n][s]=log10(pos_prob[n][s]);
          }
        }
        NullSitesFModel(L->null, finalmodel[t]);
        PutSites(fptr,t,mapsites,pos_prob,L->null);
        PutFModel(fptr,finalmodel[t]);
        for(n = 1; n <= N; n++) free(pos_prob[n]);
        free(pos_prob);
    }
}

void    PrintCMSA(FILE *fptr, cma_typ L)
/*************************************************************************
 output the alignment (WARNING: you probably want to invoke InitMAPCMSA(L)!!!)
 *************************************************************************/
{
    Int4        t,k;
    double      map,map2;

    k = nBlksCMSA(L);
    PutMSA(fptr, L); 
    for(t=1; t <= k; t++){
        map = single_rel_map_cmsa(L, t);
        map2=RelMapMinusCMSA(t, L);
        fprintf(stderr,"\n(%d)\tsingle block NetMAP = %.2f\n",
                        t,map);
        fprintf(stderr,"   \tNetMAP with block deleted = %.2f\n", map2);
    }
}

void	PutSeqAlnCMSA(FILE *fptr, cma_typ L)
{ 
    ss_type		data=DataCMSA(L);
    st_type		sites = SitesCMSA(L);
    fm_type		*model=ModelsCMSA(L);
    Int4		n,s,t,end,N,ntyps;
    Int4 		i,e,start,length,offset,*maxgap,gap;
    unsigned char	*seq;
    BooLean		*null = nullCMSA(L);
    e_type		E;
    char		r,c;
    Int4		*pos;
    
    N = NSeqsSeqSet(data); 
    ntyps=nTypeSites(sites);
    for(n=0,t=1; t <= ntyps; t++) n+=nColsFModel(model[t]);
    NEW(maxgap,ntyps+3,Int4);
    for(t=1; t < ntyps; t++) {
	maxgap[t]=0;
        for(n = 1; n <= N; n++) {
	  e = SitePos(t,n,1,sites) + SiteLen(t,sites);
	  gap = SitePos(t+1,n,1,sites) - e;
	  if(gap > maxgap[t]) maxgap[t]=gap;
	}
    }
    fprintf(fptr,"     ");
    for(t=1; t <= ntyps; t++) {
        if(null != NULL){
	       for(s=1; s <= SiteLen(t,sites); s++){
                      if(null[s]==TRUE)fprintf(fptr,".");
                      else if(null[s]==FALSE)fprintf(fptr,"*");
                      else fprintf(fptr,"^");
               }
        } else for(s=1; s <= SiteLen(t,sites); s++) fprintf(fptr,"*");
	for(gap=1; gap <= maxgap[t]; gap++,i++) fprintf(fptr," ");
    }
    fprintf(fptr,"\n");
    for(n = 1; n <= N; n++) {
        E=SeqSetE(n,data);
        offset = OffSetSeq(E);
        fprintf(fptr,"%4d ",SitePos(1,n,1,sites) + offset);
    	for(t=1; t <= ntyps; t++) {
           length = SiteLen(t,sites);
           start= SitePos(t,n,1,sites);
           e = start + length - 1;
           end = e;
           for(i=start; i <= end; i++){
                   if(i < 1 || i > (Int4) LenSeq(E)) fprintf(fptr," ");
                   else{
                        r = ResSeq(i,E);
                        c = AlphaChar(r,SeqSetA(data));
                        fprintf(fptr,"%c", c);
                   }
           }
	   if(t < ntyps){
	     for(gap=1; gap <= maxgap[t]; gap++,i++){
		if(i < SitePos(t+1,n,1,sites)){
			r = ResSeq(i,E);
			c = AlphaChar(r,SeqSetA(data)); c = tolower(c);
                        fprintf(fptr,"%c", c);
		} else fprintf(fptr,".");
	     }
	   }
        }
        fprintf(fptr," %-4d ",end+offset);
        PutShortSeqID(fptr,E);
        fprintf(fptr,"\n");
    }
    free(maxgap);
    fprintf(fptr,"\n");
}

void	PutAlnCMSA(FILE *fptr, cma_typ cma)
{
    gss_typ	*gss=gssCMSA(cma);
    if(gss->Gapped()) PutGappedAlnCMSA(fptr,cma,NULL);
    else PrintAlnCMSA(fptr,NameSeqSet(DataCMSA(cma)),cma,NULL);
}

void	PutAlnCMSA(char *name, cma_typ cma,gd_type G,Int4 Rpts)
{ 
    FILE	*fptr;
    gss_typ	*gss=gssCMSA(cma);
    
    assert(Rpts > 0);
    fptr = open_file(name,".msa","w");
    if(gss->Gapped()) PutGappedAlnCMSA(fptr,cma,G,Rpts);
    else {
	if(Rpts > 1) print_error("Ungapped repeats option not implemented");
	PrintAlnCMSA(fptr,name,cma,G);
    }
    fclose(fptr);
    fptr = open_file(name,".cma","w");
    PutCMSA(fptr,cma); fclose(fptr);
}

void	PutRepSetCMSA(FILE *fp, Int4 percent_ident,Int4 *Nset,cma_typ cma)
{ PutRepSetCMSA(0, fp, percent_ident,Nset,cma); }

void	PutRepSetCMSA(FILE *fp_err, FILE *fp, Int4 percent_ident,Int4 *Nset,cma_typ cma)
{
	BooLean	*list=RtnRepSetCMSA(fp_err,percent_ident,Nset,cma);
	PutSelectCMSA(fp,list,cma); free(list); 
}

BooLean	*RtnRepSetCMSA(FILE *fp_err, Int4 percent_ident,Int4 *Nset,cma_typ cma)
// return a representative set of sequences from cma.
{
	ds_type sets;
	Int4	b,i,j,k,N=NumSeqsCMSA(cma);
	Int4	score,nblk=nBlksCMSA(cma);
	Int4	s,si,sj,pos[3],seti,setj;
	ss_type	data=DataCMSA(cma);
        a_type  A=AlphabetCMSA(cma);
	unsigned char	*isq,*jsq;
	Int4	start=1;

	if(percent_ident < 0){ 	// then keep first sequence in cma output file...
		start=2; percent_ident = -percent_ident;
	}
	assert(percent_ident > 0 && percent_ident <= 100);
	// 1. Cluster sequences into sets at the percent identity cutoff.
	Int4 total = TotalLenCMSA(cma);
	sets = DSets(N);
	for(i = start; i < N; i++){
	  isq = SeqPtrCMSA(i,cma);
	  seti=findDSets(i,sets);
	  if(fp_err && i % 1000 == 0) fprintf(stderr,"\r%.1f",100.0*((double)i/(double)N));
	  for(j=i+1; j <= N; j++){
	     setj=findDSets(j,sets);
	     if(seti != setj){
	      jsq = SeqPtrCMSA(j,cma);
	      for(score=0,b=1; b <= nblk; b++){
		PosSiteCMSA(b,i,pos,cma); si=pos[1];
		PosSiteCMSA(b,j,pos,cma); sj=pos[1];
		for(s=1; s <= LengthCMSA(b,cma); s++,si++,sj++){
			if(isq[si] == jsq[sj]) score++;
		}
	      } // score = (Int4) floor(((double)score*100.0/(double)total) +0.5);
	      score = (Int4) floor(((double)score*100.0/(double)total));
	      if(score >= percent_ident) seti=linkDSets(seti,setj,sets);
	     }
	  }
	}
	// 2. Within each set pick the sequence with the highest profile score.
	double	iprob,jprob;
	Int4	best;
	BooLean	*skipseq,*bestskipseq;
	NEW(skipseq, NumSeqsCMSA(cma)+2, BooLean);
        for(i=1; i <= N; i++) skipseq[i] = TRUE;
        for(s=0,i=1; i <= N; i++){
	  seti=findDSets(i,sets);
	  if(i==seti){	// is this the canonical sequence?
	     best=i;
	     iprob=GetTotalProbCMSA(i,cma);
	     for(j=1; j <= N; j++){
		if(j!=i){
		  setj=findDSets(j,sets);
		  if(seti == setj){
			jprob=GetTotalProbCMSA(j,cma);
			if(jprob > iprob) best = j; 
		  }
		}
	     } skipseq[best] = FALSE; s++;
	  }
	} *Nset=s;
	// 3. output the cma skipping all but the representive sequences
	NilDSets(sets);
	return skipseq;
}

void	PutGoodAlnCMSA(char *name, cma_typ L,double cutoff, gd_type G)
{ 
    gss_typ	*gss=gssCMSA(L);
    
    FILE *fptr = open_file(name,".msa","w");
    if(gss->Gapped()) PutGappedAlnCMSA(fptr,L,G);
    else PrintAlnCMSA(fptr,name,L,G);
    fclose(fptr); fptr = open_file(name,".cma","w");
    PutGoodCMSA(fptr,cutoff,L); fclose(fptr);
}


void	PrintAlnCMSA(FILE *fptr, char *name, cma_typ L, gd_type G)
{
    ss_type		data=DataCMSA(L);
    st_type		sites = SitesCMSA(L);
    fm_type		*model=ModelsCMSA(L);
    Int4		n,s,t,end,N,ntyps;
    Int4 		i,e,start,length,offset;
    double		**pos_prob,pernats=1000.0;
    unsigned char	*seq;
    BooLean		*null = nullCMSA(L);
    e_type		E;
    char		r,c;
    smx_typ		*smx;
    Int4		*pos;
    mdl_typ		*mdl=mdlCMSA(L);

    fprintf(fptr,"//\nID   XXX\nAC   A00000\nDE   %s\n",name);
    N = NSeqsSeqSet(data); 
    ntyps=nTypeSites(sites);
    if(G != NULL){
	E = GuideE(G);
	NEW(smx, ntyps+2, smx_typ);
	NEW(pos, LenSeq(E)+2, Int4);
	for(t=1; t <= ntyps; t++) {
		smx[t] = GetSmatrixFModel(pernats,model[t]);
	}
	AlnSeqSMatrix(LenSeq(E), SeqPtr(E), ntyps, smx, pos);
    }
    for(n=0,t=1; t <= ntyps; t++) n+=nColsFModel(model[t]);
    fprintf(fptr,"LPR  %.2f.\n", RelMapCMSA(L));
    fprintf(fptr,"NU   %d blocks; %d columns.\n",ntyps,n);
    fprintf(fptr,"SQ   %d sequences.\n",N);
    for(n = 1; n <= N; n++) {
        fprintf(fptr,"     ");
        PutSeqInfo2(fptr,SeqSetE(n,data));
    }   
    for(t=1; t <= ntyps; t++) {
// PutFModel(stderr, model[t]); // testing...
    	fprintf(fptr,"BLK  %d: %d residues (%d columns).\n",
		t,LenFModel(model[t]),nColsFModel(model[t]));
	fprintf(fptr,"MAP  %.1f (%.1f only; %.1f w/o).\n",
		FieldRelMapCMSA(L, t),
		single_rel_map_cmsa(L, t),
		RelMapMinusCMSA(t, L));
        NEWP(pos_prob,N+2,double);
        for(n = 1; n <= N; n++) {
                end = SqLenSeqSet(n,data) + 1;
                NEW(pos_prob[n],end+5,double);
                end += 1 - LenFModel(model[t]);
                seq = SeqSeqSet(n,data);
                for(s= 1; s<= end; s++){
                   if(TypeSite(n,s,sites) == t){
			mdl->Remove(t,seq,s);
                        pos_prob[n][s] = (double)
                           LikelihoodFModel(seq, s, model[t]);
			mdl->Add(t,seq,s);
                   } else pos_prob[n][s] = (double)
                           LikelihoodFModel(seq, s, model[t]);
                   if(pos_prob[n][s]>0.0) pos_prob[n][s]=log10(pos_prob[n][s]);
                }
        }
        NullSitesFModel(null, model[t]);
	if(G != NULL){
	   for(s=pos[t],n=1; n <= LenFModel(model[t]); n++,s++){
		if(GuideIgnore(s,G)) null[n] = 2;
	   }
	}
        fprintf(fptr,"AL   ");
        if(null != NULL){
                for(s=1; s <= SiteLen(t,sites); s++){
                        if(null[s]==TRUE)fprintf(fptr,".");
                        else if(null[s]==FALSE)fprintf(fptr,"*");
                        else fprintf(fptr,"^");
                }
        } else for(s=1; s <= SiteLen(t,sites); s++) fprintf(fptr,"*");
        fprintf(fptr,"\n");
        length = SiteLen(t,sites);
        for(n=1; n<= N; n++){
           E=SeqSetE(n,data);
           offset = OffSetSeq(E);
           for(s=1; s<= nSites(t,n,sites); s++){
                start= SitePos(t,n,s,sites);
                e = start + length - 1;
                end = e;
#if 0
/** try this?? Now used in Scanheap **/
	PutSeqRegionFormatMSA(fptr,start, length, E, pos_prob[n][start], SeqSetA(data));
#endif
#if 1
                fprintf(fptr,"     ");
                for(i=start; i <= end; i++){
                   if(i < 1 || i > (Int4) LenSeq(E)) fprintf(fptr," ");
                   else{
                        r = ResSeq(i,E);
                        c = AlphaChar(r,SeqSetA(data));
                        fprintf(fptr,"%c", c);
                   }
                }
                fprintf(fptr," %d-%d ",start+offset,end+offset);
                PutShortSeqID(fptr,E);
                fprintf(fptr," %.2f",pos_prob[n][start]);
                fprintf(fptr,"\n");
#endif
           }
        }
        for(n = 1; n <= N; n++) free(pos_prob[n]);
        free(pos_prob);
    }
    if(G != NULL){
	for(t=1; t <= ntyps; t++) { NilSMatrix(smx[t]); }
	free(pos); free(smx);
    }
}

/************************* GAPPED SEQ_SET ****************************/
Int4    cmsa_fix_alignment(Int4 lenM, Int4 N,char **alignment, BooLean *null)
{
        char	**temp,r;
        Int4    alnlen,i,j,J,k,*index;
        BooLean aln;

        NEW(index,N+3,Int4);
        ///// 1. Find alignment length
        for(alnlen=j=0; j < lenM; alnlen++){ //for each alignment position:
           for(aln=TRUE,i=1; i<=N; i++){     // Is there a insertion?
                if(islower(alignment[i][index[i]])){ aln=FALSE; break; }
           }  // If not then move to the next column for all sequences.
           if(aln) for(j++,i=1; i<= N; i++) index[i]++;
           else for(i=1; i<= N; i++){  
                if(islower(alignment[i][index[i]])) index[i]++;
           } // else move indices for sequences with gaps only.
        } /// 2. Allocate memory for alignment length.
        NEWP(temp,N+3,char);
        for(J=0; J <= N; J++){ NEW(temp[J],alnlen+3,char); index[J]=0; }
        ///// 3. Find configuration and fix alignment.
        for(k=0,j=1; j<= lenM; k++){
           for(aln=TRUE,i=1; i<= N; i++){ // does this column have an insert?
                if(islower(alignment[i][index[i]])){ aln=FALSE; break; }
           }
           if(aln){     // no inserts in this region
	     if(null[j]) temp[0][k] = '.'; else temp[0][k] = '*';
             for(j++,i=1; i<= N; i++){  // mover index for all sequences
                temp[i][k]=alignment[i][index[i]]; index[i]++;
             }
           } else {     // else add inserts and '.' characters.
             temp[0][k] = '_';   // -> delete these columns for goscan
             for(i=1; i<= N; i++){
                r=alignment[i][index[i]];
                if(islower(r)) { temp[i][k] = r; index[i]++; }
                else temp[i][k] = '.';
             }
           }
        }
        for(i=1; i<= N; i++){ free(alignment[i]); alignment[i]=temp[i]; }
        alignment[0] = temp[0]; free(temp); free(index);
        return alnlen;
}

void	PutGappedAlnCMSA(FILE *fp,cma_typ cma,gd_type G,Int4 Rpts)
/// WARNING: assumes that seqset for cma is derived from gss_typ.
// eventually change so that uses gss_typ& gss = gssCMSA(cma);
{
	gss_typ	*gss = gssCMSA(cma);
	Int4	lenM,LenStr,AlignLen,J,m,x,*s;
	char	**alignment;
	st_type	S=SitesCMSA(cma);
	ss_type	data=DataCMSA(cma);
	a_type	A=SeqSetA(data);
	e_type	sgE,E;
	fm_type	*M=ModelsCMSA(cma);
	Int4	site,number=gss->NumSeq(),buffer_size;
	char	*ps;
	BooLean	**null;
// Guide type
    smx_typ	*smx;
    Int4	*pos,i,j;

    if(G != NULL){
        double      pernats=1000.0;
        e_type	gE = GuideE(G);
	NEW(smx, nBlksCMSA(cma)+2, smx_typ);
	NEW(pos, LenSeq(gE)+2, Int4);
	for(m=1; m <= nBlksCMSA(cma); m++) {
		smx[m] = GetSmatrixFModel(pernats,M[m]);
	}
	AlnSeqSMatrix(LenSeq(gE),SeqPtr(gE),nBlksCMSA(cma),smx,pos);
    }
// end guide type...
	// if(data != DataCMSA(cma)) print_error("PutGMSA( ) input error");
	fprintf(fp,"//\nID   XXX\nAC   A00000\nDE   %s\n","testfile");
	fprintf(fp,"LPR  %.2f.\n", RelMapCMSA(cma));
	for(x=0,m=1; m <= nBlksCMSA(cma); m++) x += nColsFModel(M[m]);
	fprintf(fp,"NU   %d blocks; %d columns.\n",Rpts*nBlksCMSA(cma),Rpts*x);
	fprintf(fp,"SQ   %d sequences.\n",number);
    	for(J=1; J <= number; J++) {
		E = gss->TrueSeq(J);
        	fprintf(fp,"     "); PutSeqInfo2(fp,E);
	}
	// PrintAlnCMSA(fp, "test file", cma,NULL);
	buffer_size = gss->MaxTrueSqLen();
	for(J=1; J <= number; J++) buffer_size+=gss->NumIns(J)+gss->NumDel(J);
	s = new Int4[number+3];
	null = NullCMSA(cma);
    for(Int4 r=1; r <= Rpts; r++){
	for(m=1; m <= nBlksCMSA(cma); m++){

	   //// 1. Get locations of sites in gapped sequences from cmsa.
	   lenM = SiteLen(m,S);
	   NEWP(alignment,number+3,char);
	   for(J=1; J <= number; J++){ 
		sgE = gss->FakeSeq(J);
		s[J] =SitePos(m,J,1,S);
// fprintf(stderr,"s[%d] = %d\n",J,s[J]);
		NEW(alignment[J],buffer_size+3,char); 
	   	//// 2. Get corresponding gapped alignable segments.
		LenStr=gss->Region(J,alignment[J],s[J],lenM);
	   }
	   //// 3. Align gapped sequences.
	   AlignLen=cmsa_fix_alignment(lenM,number,alignment,null[m]);
// guide type...
	if(G != NULL){
	   Int4 end=pos[m]+LenFModel(M[m]);
	   for(i=pos[m],j=0; j < AlignLen; j++){
		if(alignment[0][j] != '_'){
		   if(GuideIgnore(i,G) && alignment[0][j] == '*') alignment[0][j] = '^';
		   assert(i < end); i++; 
		}
	   }
	}

	   //// 4. Print gapped alignment.
	   fprintf(fp,"BLK  %d: %d residues (%d with gaps; %d columns).\n",
                m+(nBlksCMSA(cma)*(r-1)),LenFModel(M[m]),AlignLen,nColsFModel(M[m]));
	   fprintf(fp,"MAP  %.1f (%.1f only; %.1f w/o).\n",
                FieldRelMapCMSA(cma, m), single_rel_map_cmsa(cma, m),
                RelMapMinusCMSA(m, cma));
	   fprintf(fp,"AL   ");
	   fprintf(fp,"%s\n",alignment[0]);
	   for(J=1; J <= number; J++){ 
		sgE = gss->FakeSeq(J);
		site = gss->TrueSite(J, s[J]);
		ps=alignment[J];
		fprintf(fp,"     %s ",alignment[J]); 
		fprintf(fp," %d-",site);
		for(ps++; *ps != 0; ) { if(isalpha(*ps)) site++; ps++; }
		fprintf(fp,"%d ",site);
		PutShortSeqID(fp,sgE);
		// fprintf(fp," %.2f\n",GetProbCMSA(m,J,cma));
		fprintf(fp," %.2f\n",GetGappedProbCMSA(m,J,cma));
	   } fflush(fp);
	   for(J=0; J <= number; J++){ free(alignment[J]); } free(alignment);
	   // alignment[0] allocated above....
	}
    } fprintf(fp,"\n");
	for(Int4 t=1; t <= nBlksCMSA(cma); t++) free(null[t]); free(null);
	delete []s;
    if(G != NULL){
	for(m=1; m <= nBlksCMSA(cma); m++) { NilSMatrix(smx[m]); }
	free(pos); free(smx);
    }
}

void	PutAsDomTypCMSA(FILE *fp,char *name,cma_typ cma)
/// WARNING: assumes that seqset for ma is derived from gss_typ.
// eventually change so that uses gss_typ& gss = gssCMSA(cma);
{
	st_type	S=SitesCMSA(cma);
	e_type	E;
	gss_typ	*gss=gssCMSA(cma);
	char    id0[202],id[202];

	assert(gss);

	Int4	lenM,J,m,s,number=gss->NumSeq(),start,end;

	m = nBlksCMSA(cma); 
	assert(m == 1);
	lenM = SiteLen(m,S);
	id0[0]=0;
    	for(J=1; J <= number; J++) {
		E = gss->TrueSeq(J);
		StrSeqID(id,200,E);
	        if(strncmp(id,id0,200) != 0){ // New sequence; output line.
		  if(J>1) fprintf(fp,"};\n");
		  fprintf(fp,"->(%d)",J); PutShortSeqID(fp,E);
		  fprintf(fp,"={"); strcpy(id0,id); 
		}
		s = SitePos(1,J,1,S);  // gives position in fakeE.
		start = gss->TrueSite(J,s);
		s = SitePos(m,J,1,S)+lenM-1;  // gives position in fakeE.
		end = gss->TrueSite(J,s);
		fprintf(fp,"[%d-%d->%s:%d-%d(0);%.2f]",
			start,end,name,1,lenM,GetGappedProbCMSA(m,J,cma));
	} fprintf(fp,"};\n");
}

