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

#if 1	// copied from cma_gmb.cc; still a copy there as well.
static Int4 RtnMisMatch(char sq1[8], char sq2[8])
// Inner loop; most time intensive...
{
	register Int4   n=0;
	if(sq1[0] != sq2[0]) n++;
	if(sq1[1] != sq2[1]) n++;
	if(sq1[2] != sq2[2]) n++;
	if(sq1[3] != sq2[3]) n++;
	if(sq1[4] != sq2[4]) n++;
	if(sq1[5] != sq2[5]) n++;
	if(sq1[6] != sq2[6]) n++;
	if(sq1[7] != sq2[7]) n++;
	return n;
}
#endif

#if 1	// moved from cma_gmb.cc

static Int4 *RtnWorstColumnsCMSA(set_typ Set, cma_typ cma)
// The maximum = ln 20 = 2.9957 (or ln 21 = 3.044).
// find worst positions...
{
        Int4    i,j,t,k,sq,r,n;
        double  p,sum,entropy,cnt[50],tcnt;
        a_type AB = AlphabetCMSA(cma);

	dh_type dH=dheap(TotalLenCMSA(cma)+3,4);
	Int4	*list; NEW(list,TotalLenCMSA(cma) +3,Int4);

     for(sum=0.0,n=0,t=1; t <= nBlksCMSA(cma); t++){
        for(k=1; k <= LengthCMSA(t,cma); k++){
          n++;
          for(tcnt=0.0,r=0; r <= nAlpha(AB); r++) cnt[r]=0.0;
          for(sq=1; sq <= NumSeqsCMSA(cma); sq++){
	      if(Set && !MemberSet(sq,Set)) continue;
              r=ResidueCMSA(t,sq,k,cma);
              if(r){ cnt[r]+=1.0; tcnt +=1.0; }
          }
          for(entropy=0.0,r=0; r <= nAlpha(AB); r++){
                p = cnt[r]/tcnt; 
                if(p > 0.0) entropy += -p * log(p);
          } sum += entropy;
	  insrtHeap(n,(keytyp)-entropy,dH);
        }
     }
     for(i=1; !emptyHeap(dH); i++){
	// entropy=-minkeyHeap(dH);
	n=delminHeap(dH);  list[i]=n; 
	// fprintf(stderr,"%d.%d %.3f\n",i,n,entropy);
     } Nildheap(dH); return list;
}

static UInt8 **RtnBlkAln(Int4 &NumBlks,Int4 &MaxMisMatch, Int4 percent_ident,set_typ InSet,cma_typ cma)
{
	Int4	b,i,j,k,N=NumSeqsCMSA(cma),score,nblk=nBlksCMSA(cma);
	Int4	s,sb,bl,si,sj,i8,pos[3],seti,setj;
	unsigned char	r;
        a_type  AB=AlphabetCMSA(cma);
	Int4	start=1;
	FILE	*efp=0; // efp=stderr;

	if(percent_ident < 0){ 	// then keep first sequence in cma output file...
		start=2; percent_ident = -percent_ident;
	} assert(percent_ident > 0 && percent_ident <= 100);

	// [0,1,2,3, 4,5,6,7] = 8 residues at a time.
	assert(sizeof(UInt4) != sizeof(UInt8)); // make sure this is a 64 bit machine...
	assert(sizeof(UInt8)==8); // make sure this is a 64 bit machine...
	char	*bsq; 
	Int4	col,TotLen=TotalLenCMSA(cma);
	double	d=(double) (TotLen*((double)(100-percent_ident)/100.0));
	double	D=(double) TotLen/8.0;
	MaxMisMatch= (Int4) floor(d); NumBlks=(Int4) ceil(D); 

	// 0. Store each fake sequence in 8 byte arrays.
	UInt8	**SQ; NEWP(SQ,N+3,UInt8);
#if 0	// Original code.
	for(i = start; i <= N; i++){
  	  if(!MemberSet(i,InSet)) continue;
	  NEW(SQ[i],(TotLen/7)+2,UInt8);
#if 0	  // purify uninitialized read error due to short seq?
	  unsigned char *isq = SeqPtrCMSA(i,cma);
	  for(bl=sb=0,b=1; b <= nblk; b++){
	    PosSiteCMSA(b,i,pos,cma); si=pos[1];
	    for(s=1; s <= LengthCMSA(b,cma); s++,si++,sj++){
		if(sb%8==0){ bsq=(char *) &SQ[i][bl]; bl++; sb=0; } 
		bsq[sb] = AlphaChar(isq[si],AB); sb++;
	    }
	  } // fprintf(stderr,"%d: %s\n",i,SQ[i]);
#else
	  for(bl=sb=0,b=1; b <= nblk; b++){
	    for(s=1; s <= LengthCMSA(b,cma); s++){
		if(sb%8==0){ bsq=(char *) &SQ[i][bl]; bl++; sb=0; } 
		r=ResidueCMSA(b,i,s,cma);
		bsq[sb] = AlphaChar(r,AB); sb++;
	    }
	  } // fprintf(stderr,"%d: %s\n",i,SQ[i]);
#endif
	}
#else	//====== NEW faster method.. ===========
	// sort the coluumns from least conserved first to hit mismatches sooner!
	assert(nblk == 1);
	Int4 *list=RtnWorstColumnsCMSA(InSet,cma);
	for(i=start; i <= N; i++){
  	  if(!MemberSet(i,InSet)) continue;
	  NEW(SQ[i],(TotLen/7)+2,UInt8);
	  for(bl=sb=0,j=1; j <= LengthCMSA(1,cma); j++){
		s=list[j];
		if(sb%8==0){ bsq=(char *) &SQ[i][bl]; bl++; sb=0; } 
		r=ResidueCMSA(1,i,s,cma); bsq[sb] = AlphaChar(r,AB); sb++;
	  } // fprintf(stderr,"%d: %s\n",i,SQ[i]);
	} free(list);
#endif
	if(efp) fprintf(efp,"NumBlks = %d; TotLen= %d; MaxMisMatch=%d\n",NumBlks,TotLen,MaxMisMatch);
	return SQ;
}

ds_type GetRepSetsStatic(Int4 start,Int4 N,Int4 NumBlks,Int4 MaxMisMatch,
				set_typ InSet, Int4 *edges,UInt8 **SQ)
// Don't make this static as the optimization option eliminates this and makes it slower!!!
{
	ds_type sets = DSets(N);
	register Int4	i,j,seti,setj,score,bl,sb;
	for(i = start; i < N; i++){
  	  if(!MemberSet(i,InSet)) continue;
	  seti=0;
	  // if(fp && i % 1000 == 0) fprintf(fp,"\r%.1f",100.0*((double)i/(double)N));
	  for(j=i+1; j <= N; j++){
  	     if(!MemberSet(j,InSet)) continue;
	     if(seti == 0) seti=findDSets(i,sets);
	     setj=findDSets(j,sets);
	     if(seti != setj){
	      UInt8 *SQ_i=SQ[i],*SQ_j=SQ[j];
	      for(score=0,bl=0,sb=0; bl <= NumBlks; bl++){
		if(SQ_i[bl] != SQ_j[bl]){
		   score += RtnMisMatch((char *)&SQ_i[bl],(char *)&SQ_j[bl]);
		   if(score > MaxMisMatch) break;
		}
	      }
	      if(score <= MaxMisMatch) { 
		   edges[i]++; edges[j]++;
		   seti=linkDSets(seti,setj,sets);
	      }
	     }
	  }
	} return sets;
}

set_typ	RtnFastRepSetCMSA(FILE *fp, Int4 percent_ident,set_typ InSet,cma_typ cma)
// return a representative set of sequences from cma.
{
	ds_type sets;
	Int4	b,i,j,k,N=NumSeqsCMSA(cma),score,nblk=nBlksCMSA(cma);
	Int4	s,sb,bl,si,sj,i8,pos[3],seti,setj;
        a_type  AB=AlphabetCMSA(cma);
	unsigned char	r;
	Int4	start=1;

	Int4	NumBlks=0,MaxMisMatch=0; 
	UInt8	**SQ=RtnBlkAln(NumBlks,MaxMisMatch,percent_ident,InSet,cma);
	if(fp) fprintf(fp,"NumBlks=%d; MaxMisMatch=%d\n",NumBlks,MaxMisMatch);
	if(percent_ident < 0) start=2; 
	assert(percent_ident > 0 && percent_ident <= 100);

	//============ 1. Cluster sequences into sets at the percent identity cutoff. ============
	
	Int4	*edges; NEW(edges,N+3,Int4);	// number of edges out of each node.
	h_type	HG=0,HG2=0,HG3=0;
	if(fp) HG=Histogram("number of edges",0,25,1.0);

#if 1	// new routine.
	sets=GetRepSetsStatic(start,N,NumBlks,MaxMisMatch,InSet,edges,SQ);
#else
	sets = DSets(N);
	for(i = start; i < N; i++){
  	  if(!MemberSet(i,InSet)) continue;
	  seti=0;
	  // if(fp && i % 1000 == 0) fprintf(fp,"\r%.1f",100.0*((double)i/(double)N));
	  for(j=i+1; j <= N; j++){
  	     if(!MemberSet(j,InSet)) continue;
	     if(seti == 0) seti=findDSets(i,sets);
	     setj=findDSets(j,sets);
	     if(seti != setj){
	      UInt8 *SQ_i=SQ[i],*SQ_j=SQ[j];
	      for(score=0,bl=0,sb=0; bl <= NumBlks; bl++){
		if(SQ_i[bl] != SQ_j[bl]){
		   score += RtnMisMatch((char *)&SQ_i[bl],(char *)&SQ_j[bl]);
		   if(score > MaxMisMatch) break;
		}
	      }
	      if(score <= MaxMisMatch) { 
		   edges[i]++; edges[j]++;
		   seti=linkDSets(seti,setj,sets);
	      }
	     }
	  }
	}
#endif
	if(fp){
	  for(i=1; i <= N; i++){
	    IncdHist((double)edges[i],HG);
	    // if(edges[i] > 70) PutSeq(stderr,TrueSeqCMSA(i,cma),AB);
	  } PutHist(fp,50,HG); NilHist(HG);
	}
	for(i = start; i <= N; i++) if(SQ[i]) free(SQ[i]);  free(SQ);

	// 2. Within each set pick the sequence with the highest score.
	double	bestprob,jprob,var;
	Int4	best;
	set_typ Set=0,SubSet=MakeSet(SetN(InSet)); ClearSet(SubSet);
	if(fp){
	  HG=Histogram("number of edge difference",0,25,1.0);
	  // Int4 inc=(Int4) ceil((double) CardSet(InSet)/100.0); 
	  HG2=Histogram("seqs set sizes",0,30,1);
	  HG3=Histogram("number of seqs in each set",0,100,1);
	  Set=MakeSet(N+3); // ClearSet(Set);
	}
        for(s=0,i=1; i <= N; i++){
  	  if(!MemberSet(i,InSet)) continue;
	  seti=findDSets(i ,sets);
	  if(i != seti) continue;	// skip non-canonical sequences
	  best=i;
	  // bestprob=GetTotalProbCMSA(i,cma);
	  bestprob=(double)edges[i];
	  if(fp){ ClearSet(Set); AddSet(i,Set); }
	  for(j=1; j <= N; j++){
  	        if(!MemberSet(j,InSet)) continue;
		if(j!=i){
		  setj=findDSets(j,sets);
		  if(seti == setj){
		        if(fp) AddSet(j,Set);
			// jprob=GetTotalProbCMSA(j,cma);
			jprob=(double) edges[j];
			if(jprob > bestprob){ best = j; bestprob=jprob; }
		  }
		}
	  } AddSet(best,SubSet); 
	  if(fp){
	     Int4 card=CardSet(Set);
	     if(card > 2){ s++; IncdMHist((double)s,card, HG3); }
	     IncdHist(card,HG2);
	     for(j=1; j <= N; j++){
		if(!MemberSet(j,Set)) continue;
		for(Int4 jj=1; jj <= N; jj++){
		  if(!MemberSet(jj,Set) || j == jj) continue;
		  IncdHist(abs(edges[j]-edges[jj]),HG);
  		}
	     }
	  }
	}
	if(fp){
	   PutHist(fp,50,HG); NilHist(HG); PutHist(fp,50,HG2); NilHist(HG2);
	   PutHist(fp,50,HG3); NilHist(HG3); NilSet(Set);
	} NilDSets(sets); free(edges);
	return SubSet;
}

#endif


cma_typ RmWrinklesCMSA(cma_typ cma)
// remove insertions next to deletions.
{
   assert(nBlksCMSA(cma)==1);
   a_type AB=AlphabetCMSA(cma);
   gss_typ *gss=gssCMSA(cma);
   Int4 len[5]; len[1]=LengthCMSA(1,cma);
   cma_typ rcma=EmptyCMSA(1,len,TrueDataCMSA(cma),gss->GapOpen(),gss->GapExtend(),
                        PerNatsCMSA(cma),0,0);
   Int4 i,pos[9];  pos[2]=0;
   gss_typ *rgss=gssCMSA(rcma);
   for(i=1; i<= NumSeqsCMSA(cma); i++){
        gsq_typ *gsq=gsqCMSA(i,cma);
        // gsq->Put(stderr,AB);
        // gsq->Put_cma_format(stderr,i,nBlksCMSA(cma),sites,LengthsCMSA(cma),AB);
        Int4    *sites=GetPosSitesCMSA(i,cma);
        Int4    X;
        gsq_typ *gsq0 = gsq->IronOut(nBlksCMSA(cma),sites,LengthsCMSA(cma),X);
        // gsq_typ *gsq0 = gsq->ToOneBlk(nBlksCMSA(cma),sites,LengthsCMSA(cma),X);
        // gsq0->Put_cma_format(stderr,i,1,pos,LengthsCMSA(rcma),AB);
        ReplaceCMSA(i,gsq0,rcma); // replace sequence s in CMSA & fmodel.
        AddSiteCMSA(1,i,X,rcma);
        free(sites);
   } return rcma;
}

cma_typ RmOverHangsCMSA(cma_typ cma)
// Rmove sequence regions beyond the ends of the alignment (as for cma2fa routine).
// this is needed for mapgaps.
// WARNING: returns a cma_typ with new data; need to delete using TotalNilCMSA();
{
   	a_type	AB=AlphabetCMSA(cma);
	Int4	sq,blk,nblks=nBlksCMSA(cma),*len=LengthsCMSA(cma),N=NumSeqsCMSA(cma);
	Int4	**pos,*sites; NEWP(pos,N+2,Int4);
	gsq_typ **gsq = new gsq_typ*[N+2];
	for(sq=1; sq <= N; sq++){
	   gsq_typ *gsq0=gsqCMSA(sq,cma);
	   sites=GetPosSitesCMSA(sq,cma);
	   Int4 start=TruePosCMSA(sq,1,1,cma);
	   Int4 end=TruePosCMSA(sq,nblks,len[nblks],cma);
	   gsq[sq]=gsq0->RmOverHangs(nblks,sites,len,start,end);
	   pos[sq]=sites;  // set by gsq_typ routine...
	   // if(sq==164) { gsq[sq]->Put(stderr,AB); gsq[sq]->Put(stderr,60,AB); }
        }
	gss_typ *gss = new gss_typ[1];
	gss->FromArray(N,gsq,gss->GapOpen(),gss->GapExtend(),PerNatsCMSA(cma),0,0,NameCMSA(cma),AB);
	cma_typ rcma=MakeCMSA(MakeSites(nblks,len,*gss),NullCMSA(cma));
	rcma->Level=cma->Level;
	for(sq=1; sq <= N; sq++){
                for(blk=1; blk <= nblks; blk++) AddSiteCMSA(blk,sq,pos[sq][blk],rcma);
                free(pos[sq]);
        } free(pos);
	return rcma;
}

cma_typ AddConsensusCMSA(cma_typ cma)
{
        // cma_typ      CMA[4]; CMA[1]=cma;
        Int4    Number,n;
        a_type  A=AlphabetCMSA(cma);
        if(nBlksCMSA(cma) > 1) print_error("AddConsensusCMSA( ) input error");

        // create concensus cma
        FILE *fp=tmpfile();
        PutConsensusCMSA(fp,cma);
        PutCMSA(fp,cma);
        // PutConsensusCMSA(stderr,cma);
        rewind(fp);
        cma_typ *IN_CMA=MultiReadCMSA(fp,&Number,A);
        fclose(fp);

        // merge csq.cma and original cma.
        fp=tmpfile();
        PutMergedCMSA(fp,Number,IN_CMA);
        rewind(fp);
        cma_typ CMA=ReadCMSA(fp,A);
        fclose(fp);
        for(n = 1; n<= Number; n++) TotalNilCMSA(IN_CMA[n]);
        free(IN_CMA);
        return CMA;
}

cma_typ SimSeqToCMSA(e_type *EList, Int4 N, a_type A)
// convert simulated sequences to a cma file.
{
        Int4    len[3],n,Len;
        st_type S;
        BooLean *null[3];
        cma_typ cmsa;

	Len=LenSeq(EList[1]);
        for(n=1; n<=N; n++) assert(Len == LenSeq(EList[n]));
	ss_type data=Array2SeqSet(EList, N, "simulated",A);
        len[1]=Len;
        S = MkSites(1,len,data,10,5,5.0,0,0);
        NEW(null[1],Len+2,BooLean);
        for(n=1; n<=Len; n++) null[1][n] = FALSE;
        cmsa=MakeCMSA(S, null);
        for(n=1; n<=N; n++) AddSite(1,n,1,S);  // all seq same length.
	free(null[1]);
        return cmsa;
}

cma_typ	ShuffleSeqCMSA(double fraction, cma_typ cma)
// Shuffle all of the sequences for cma: can use this instead of a std background 
// for CHAIN analysis; need to remove and add sequences???  CAREFUL!
{
	assert(nBlksCMSA(cma) == 1);
	ss_type data = TrueDataCMSA(cma);
	Int4	m,n,N=NumSeqsCMSA(cma);
	unsigned char	r,*seq;
	Int4	Len=LengthCMSA(1,cma),col;
	e_type	*SeqE,sE,rE;
	double	*freq=tFreqSeqSet(data);

        assert(nBlksCMSA(cma) == 1);
	NEW(SeqE,N+3,e_type);
        for(n=1; n <= N; n++) {
		Int4 col;
		sE=EmptySeq(n, Len);
	        for(m=0,col=1; col <=Len; col++){
			r=ResidueCMSA(1, n, col,cma); EqSeq(col,r,sE);
        	}
		if(fraction <= 0.0){ SeqE[n]=sE; }
		else {
		  rE=RandomSeq(Len, n, freq,SeqSetA(data));
		  if(fraction >= 1.0){ SeqE[n]=rE; NilSeq(sE); }
		  else {
		    dh_type dH=dheap(Len+2,4);
		    for(col=1; col <= Len; col++){
			r=ResSeq(col,sE);
			if(r) insrtHeap(col,(keytyp) Random(),dH);
		    }
		    m=0; seq=SeqPtr(sE);
		    while((col=delminHeap(dH)) != 0){
			m++; seq[col] = ResSeq(col,rE);
			if(((double)m/(double)Len) >= fraction) break;
		    } NilSeq(rE); Nildheap(dH);
		    SeqE[n]=sE;
		  }
		}
	} return SimSeqToCMSA(SeqE, N,SeqSetA(data));
}

void	ReplaceCMSA(Int4 s, gsq_typ *gsq, cma_typ cma)
// Replace sequence n with gsq putting site[1->ntyp] into alignment
{ assert(SwapGsqCMSA(s,gsq,cma,FALSE) == 0); }

gsq_typ	*SwapGsqCMSA(Int4 s, gsq_typ *gsq, cma_typ cma,BooLean swap)
// Replace sequence n with gsq putting site[1->ntyp] into alignment
{
	Int4	t,max=0,oldmax;
	e_type	E;
	gsq_typ	*rtn_gsq=0;
	
	if(gsq == NULL) E = SeqSetE(s,TrueDataCMSA(cma)); else E = gsq->FakeSeq();
#if 0	// old length...
	oldmax=MaxTrueSeqCMSA(cma);
#else
	oldmax=cma->maxlen[1];
#endif
	if(LenSeq(E) > oldmax){
	   BooLean *bestnull;
	   max=LenSeq(E); 
	   free(cma->pos); free(cma->null);
           NEW(cma->pos,(2*max)+2,Int4); NEW(cma->null, (2*max)+2,BooLean);
           for(t=1; t <= nBlksCMSA(cma); t++){
                NEW(bestnull, (2*max)+2, BooLean);
                for(Int4 i = 0; i <= oldmax; i++) bestnull[i] = cma->bestnull[t][i];
		free(cma->bestnull[t]); cma->bestnull[t]=bestnull;
#if 0
fprintf(stderr,"Sq %d: max=%d; cma->maxlen[t] = %d; LengthCMSA(t,cma)=%d\n",
					s,max,cma->maxlen[t],LengthCMSA(t,cma));
#endif
                cma->maxlen[t] = max;
           }
	}
	VacateSitesCMSA(s,cma);
// std::cerr << "replacing sequence s in data\n";
	if(swap){
	  rtn_gsq=SwapSeqSites(s,gsq,SitesCMSA(cma)); // removes all sites in s.
	} else {
	  ReplaceSeqSites(s,gsq,SitesCMSA(cma)); // removes all sites in s.
	  rtn_gsq=0;
	}
// std::cerr << "replacing sequence information in model\n";
	if(max > 0){
	   for(t=nBlksCMSA(cma); t>0; t--) EnlargeFModel(max,ModelCMSA(t,cma));
	} return rtn_gsq;
}

void    RmSiteCMSA(Int4 t,Int4 n,Int4 s, cma_typ L)
/** Remove a site from the multiple alignment **/
{
	VacateSite(t,n,s,L->sites);
	mdl_typ *mdl=mdlCMSA(L);
	mdl->Remove(t,SeqPtrCMSA(n,L),s);
}

void    AddSiteCMSA(Int4 t,Int4 n,Int4 s, cma_typ L)
/** if possible add a site to the multiple alignment **/
{
	AddSite(t,n,s,L->sites);
	mdl_typ *mdl=mdlCMSA(L);
	mdl->Add(t,SeqPtrCMSA(n,L),s);
}

void    VacateSitesCMSA(Int4 n,cma_typ cmsa)
// WARNING: this has not yet been tested...
{
        st_type S=SitesCMSA(cmsa);
        Int4    t,k,s,ntyp=nBlksCMSA(cmsa);
	mdl_typ *mdl=mdlCMSA(cmsa);

        for(t=1; t <= ntyp; t++){
           for(k=1; k <= nSites(t,n,S); k++){
                s=SitePos(t,n,k,S);
		mdl->Remove(t,SeqPtrCMSA(n,cmsa),s);
#if 0
		if(!mdl->Remove(t,SeqPtrCMSA(n,cmsa),s)){

std::cerr << "sequence = "; std::cerr << n; std::cerr << "; site = "; std::cerr << s;
std::cerr << "; block = "; std::cerr << t; std::cerr << "\n\n";
			gssCMSA(cmsa)->Put(stderr,n);
			PutSites(stderr,t,S,NULL, NULL); exit(1);
		}
#endif
                VacateSite(t,n,s,S);
           }
        }
}

Int4	AddColumnMSA(Int4 t, Int4 pos, cma_typ L)
/***************************************************************************
  WARNING: FModel and Sites have two different conventions for designating 
   sites!  GetSiteFreq numbers from 
		-inf...-1 for the left flank.
		0...len-1 for the model.
		len...+inf for right flank
	FModel numbers from
		-inf...0 for the left flank.
		1...len for the model 
		and len+1...+inf for right flank.
 ***************************************************************************/
{
	Int4	len,j,k,d,*sitefreq;
	st_type	S=L->sites;
	fm_type	M=ModelCMSA(t,L);
	mdl_typ *mdl=mdlCMSA(L);
	
	len = LenFModel(M);
	if(pos ==  1) print_error("error in AddColumnMSA()");
	sitefreq = GetSiteFreq(S,t,pos-1); /** 0..len-1 = in block **/
	// NullSitesFModel(L->null, M); PutSites(stderr,t,S,NULL,L->null);
	d=mdl->AddColumn(t,sitefreq, pos);
	if(d > 0) print_error("error in AddColumnMSA()");
	if(LenFModel(M) > len){		/** model has grown **/
		k = LenFModel(M) - len;
		// fprintf(stderr,"d = %d; k = %d; pos = %d\n",d,k,pos);
		if(d != 0) ShiftSitesM(S, t, d);
		for(j = 1; j <= k; j++) { GrowSites(t,S); }
	}
	return d;
}

BooLean	InsertColCMSA(Int4 t, BooLean right, cma_typ cma)
/***************************************************************************
  Inserts a column at one of the ends of a block whether or not there 
  is an opening by adding gaps into 'tight' sequences.
  returns TRUE if it was necessary to insert gaps in sequences.
 ***************************************************************************/
{
   Int4    n,len,j,k,d,m,*sitefreq,*location,*oldsite,num_blkd;
   char	   *blocked;
   BooLean gapped=FALSE;

   location = new Int4[NumSeqsCMSA(cma) + 2];
   blocked = new char[NumSeqsCMSA(cma) + 2];
   num_blkd=GetEdgeBlocksSite(SitesCMSA(cma),t,right,location,blocked);
   if(num_blkd > 0){
        oldsite = new Int4[nBlksCMSA(cma)+2];
        // NullSitesFModel(cma->null,ModelCMSA(t,cma)); PutSites(stderr,t,S,NULL,cma->null);
   	// fprintf(stderr,"DEBUG 3x (%d seqs)\n",NumSeqsCMSA(cma));
	for(n=1; n <= NumSeqsCMSA(cma); n++){
	   if(blocked[n]=='T'){
		for(m=1; m <= nBlksCMSA(cma); m++){
                    assert(nSites(m,n,SitesCMSA(cma))==1);
		    // for now requires that sequence n has one site of each type.
                    oldsite[m]=SitePos(m,n,1,SitesCMSA(cma));
		    RmSiteCMSA(m,n,oldsite[m], cma);
        	}
		InsertGapSites(n,location[n],1,SitesCMSA(cma)); // calls VacateSites( ).
		for(m=1; m <= nBlksCMSA(cma); m++){
#if 0	// DEBUG...
   	          fprintf(stderr,"DEBUG sq=%d; blk=%d(%d)\n",n,m,right);
   	          fprintf(stderr," oldsite[%d] = %d; location[%d]=%d\n",m,oldsite[m],n,location[n]);
#endif
		  if(right){	// MmmmmMmm --> Mmmmm-Mmmm
                   if(oldsite[m] <= location[n]) AddSiteCMSA(m,n,oldsite[m],cma);
                   else AddSiteCMSA(m,n,oldsite[m]+1,cma);
		  } else {	// mmMmmmiii --> mm-Mmmmiii
                   if(oldsite[m] < location[n]+1) AddSiteCMSA(m,n,oldsite[m],cma);
                   else AddSiteCMSA(m,n,oldsite[m]+1,cma);
		  }
		}
	   }
	} delete []oldsite; gapped=TRUE;
    }
// gssCMSA(cma)->Put(stderr);
// PutCMSA(stderr,cma);
    if(right) d=AddColumnMSA(t,LengthCMSA(t,cma)+1,cma);
    else d=AddColumnMSA(t,0,cma);
    // NullSitesFModel(cma->null,ModelCMSA(t,cma)); PutSites(stderr,t,S,NULL,cma->null);
    delete [] location; delete [] blocked;
    return gapped;
}

void	ExtendFakeToRealCMSA(Int4 sq,cma_typ cma)
{
#if 1	// subcma operations; check for vacant site ...afn: 11-26-2014.
        if(!SeqIsInCMSA(sq,cma)) return;
#else
        assert(SeqIsInCMSA(sq,cma));
#endif
	a_type	AB=AlphabetCMSA(cma);
	gsq_typ *gsq=gsqCMSA(sq,cma);
	Int4	b,*sites=GetPosSitesCMSA(sq,cma);

	char *Op=gsq->Operation(nBlksCMSA(cma),sites,LengthsCMSA(cma));
        gsq_typ *gsq0 = new gsq_typ[1]; 
        gsq0->initialize(Op,gsq->TrueSeq(),sites); free(Op);
        ReplaceCMSA(sq,gsq0,cma); // replace sequence s in CMSA & fmodel.
        for(b=1; b<=nBlksCMSA(cma); b++) AddSiteCMSA(b,sq,sites[b],cma);
	// delete gsq;
	free(sites);
}

void	ExtendFakeToRealCMSA(cma_typ cma)
{
   a_type  AB=AlphabetCMSA(cma);
   gss_typ *gss=gssCMSA(cma); gss->SetLeftFlank(500); gss->SetRightFlank(500);
   for(Int4 i=1; i<= NumSeqsCMSA(cma); i++){
#if 1	// subcma operations; check for vacant site ...afn: 11-26-2014.
        if(!SeqIsInCMSA(i,cma)) continue;
#else
        assert(SeqIsInCMSA(i,cma));
#endif
	gsq_typ *gsq=gsqCMSA(i,cma);
	Int4	*sites=GetPosSitesCMSA(i,cma);
	char *Op=gsq->Operation(nBlksCMSA(cma),sites,LengthsCMSA(cma));
        gsq_typ *gsq0 = new gsq_typ[1]; 
        gsq0->initialize(Op,gsq->TrueSeq(),sites); free(Op);
        ReplaceCMSA(i,gsq0,cma); // replace sequence s in CMSA & fmodel.
        for(Int4 t=1; t<=nBlksCMSA(cma); t++) AddSiteCMSA(t,i,sites[t],cma);
	// delete gsq;
	free(sites);
   }
}


Int4	RmColumnMSA(Int4 t, Int4 lemon, cma_typ L)
{
	Int4	len,j,k,d;
	st_type	S=L->sites;
	mdl_typ *mdl=mdlCMSA(L);
	
	len = LenFModel(ModelCMSA(t,L));
	if(lemon < 1 || lemon > len) print_error("error in RmColumnMSA()");
	d=mdl->RmColumn(t,lemon);
	if(d < 0) print_error("error in RmColumnMSA()");
	// fprintf(stderr," (new=%d; old=%d; w=%d)", i,lemon,len); 
	if(LenFModel(ModelCMSA(t,L)) < len){	/* model shrunk */
		k = len - LenFModel(ModelCMSA(t,L));
		for(j = 1; j <= k; j++) { ShrinkSites(t,S); }
		if(d != 0) ShiftSitesM(S, t, d);
	} 
	return d;
}

cma_typ RmBlkCMSA(Int4 t, cma_typ L)
{
        st_type S= SitesCMSA(L); 
        Int4    ntyp= nTypeSites(S);
        if(t < 1 || t > ntyp) return NULL;
        else return delete_blk_cmsa(L, t);
}

BooLean	TrimCMSA(Int4 blk, unsigned short RmLeft, unsigned short RmRight, cma_typ cma)
{
	Int4 i;
        if(blk < 1 || blk > nBlksCMSA(cma)) return FALSE;
	if((RmLeft+RmRight) >= LengthCMSA(blk,cma)) return FALSE;
        for(i=1; i<=RmRight; i++) RmColumnMSA(blk,LengthCMSA(blk,cma),cma);
        for(i=1; i<=RmLeft; i++) RmColumnMSA(blk, 1, cma);
	return TRUE;
}

cma_typ	TrimCMSA(float info_cut, Int4 *TrimLimit, Int4 *RmLeft, Int4 *RmRight, 
	cma_typ cmsa)
// Trim cma at ends of each block.  Put the number of columns trimmed 
// on each side of each block in RmRight[blk] & RmLeft[blk].
{
	Int4	i,s,blk;
	Int4    *len,t,n,N,x;
	float	*info;
	cma_typ	cma2=0,cma=0;
	a_type	A=AlphabetCMSA(cmsa);

	cma = CopyCMSA(cmsa);
	for(blk = nBlksCMSA(cma); blk > 0; blk--){
	   // 1. first add in all null columns.
	   for(i =1; i <= LengthCMSA(blk,cma); i++){
		if(NullSiteCMSA(blk,i,cma)) { // then add column.
		   AddColumnMSA(blk, i, cma);
		}
	   }
	   fprintf(stderr,"Motif %3d:\n",blk);
	   fm_type fm = ModelCMSA(blk,cma);
	   Int4 lenM = LenFModel(fm);
	   info = InfoFModel(fm);
	   for(i=1; i<=lenM; i++){
		if(info[i] >= 0.0){
		  fprintf(stderr,"%3d: %.2f | ",i,info[i]);
		  for(float flt=0.033; flt<=info[i]; flt+=0.033){
			if(flt > 2.0){ fprintf(stderr,"*"); break; }
			else fprintf(stderr,"=");
		  } fprintf(stderr,"\n");
		} else fprintf(stderr,"%3d: - -- |\n",i);
	   } fprintf(stderr,"\n");

	   // Determine number of columns to trim from either side.
	   Int4 r=lenM,l=1;
	   Int4 max = lenM - TrimLimit[blk];
	   RmRight[blk] = RmLeft[blk] = 0;
	   for(i = 0; i < max; i++){
		if(info[r] < info[l]){ 
		    if(info[r] <= info_cut) RmRight[blk]++; else break; r--; 
		} else { 
		    if(info[l] <= info_cut) RmLeft[blk]++; else break; l++;
		}
	   }

	   // first remove from right end.
	   if((RmRight[blk] + RmLeft[blk]) >= lenM){
		if(nBlksCMSA(cma) <= 1){ NilCMSA(cma); return NULL; }
		cma2 = RmBlkCMSA(blk, cma); NilCMSA(cma); cma=cma2;
		continue;
	   } else if(RmRight[blk] > 0){
		Int4 lemon;
		for(i=1; i<=RmRight[blk]; i++){
		   if(info[lenM-i+1] < 0.0) continue; // null site
                   lemon = LengthCMSA(blk,cma);
                   if(lemon <= 3){
			if(nBlksCMSA(cma) <= 1){ NilCMSA(cma); return NULL; }
			cma2 = RmBlkCMSA(blk,cma); NilCMSA(cma); 
			cma=cma2; continue;
		   } else RmColumnMSA(blk, lemon, cma);
                }
	   } 
	   // next remove from left end.
	   if(RmLeft[blk] > 0){
		for(i=1; i<=RmLeft[blk]; i++){
		   if(info[i] < 0.0) continue; // null site
		   if(LengthCMSA(blk,cma) <= 3){
			if(nBlksCMSA(cma) <= 1){ NilCMSA(cma); return NULL; }
			cma2 = RmBlkCMSA(blk,cma); NilCMSA(cma); 
			cma=cma2; continue;
		   } else RmColumnMSA(blk, 1, cma);
		}
	   } 
	} free(info); // not owned by fmodel.
	return cma;
}

cma_typ delete_blk_cmsa(cma_typ ma1, Int4 tdel)
/*********************************************************************
    Delete block tdel from ma1 and return as new cmsa.

 *********************************************************************/
{
        Int4    N,i,t,t1,n,nblks,tmp[30],ntyp1,*len_elem;
        fm_type *model;
        st_type S,S1;
        ss_type data;
        BooLean **null;
        cma_typ ma2;

        S1 = SitesCMSA(ma1); ntyp1 = nTypeSites(S1);
        if(tdel < 1 || tdel > ntyp1) print_error("delete_blk_cmsa() error");
        nblks=ntyp1-1;

        NEW(model,nblks+2,fm_type);
        NEW(len_elem,nblks+2,Int4);
        for(t=1,t1=1; t1<= ntyp1; t1++){
            if(t1 != tdel) model[t++] = ModelCMSA(t1,ma1);
        }

        for(t=1 ; t<=nblks; t++){ len_elem[t]=LenFModel(model[t]); }
        data = SitesSeqSet(S1); N = NSeqsSeqSet(data);
        NEWP(null,nblks+2,BooLean);
        for(t=1; t <= nblks; t++){
            NEW(null[t],MaxSeqSeqSet(data)+2,BooLean);
            NullSitesFModel(null[t], model[t]);
        }
        S = CreateNewSites(nblks,len_elem,S1);
        for(t1=t=1; t1 <= ntyp1; t1++){
          if(t1 != tdel){
            for(n=1; n <= N; n++){
              PosTSites(t1, n, tmp, S1);
              for(i=1; i<=nSites(t1,n,S1); i++) AddSite(t,n,tmp[i],S);
            }
            t++;
          }
        }
        ma2 = MakeCMSA(S, null);
	if(ma1->FullSeq) CopyFullCountsCMSA(ma2,ma1);
	if(ma1->Domains) CopyDomainsCMSA(ma2,ma1);
        for(t=1; t <= nblks; t++) free(null[t]); free(null);
        free(model); free(len_elem);
        return ma2;
}

#if 0	// moved to cmsa_io.cc
// After testing, can delete this
void	PutFastaAlnCMSA(FILE *fp,cma_typ cma)
// output a single block cmsa file in fasta format...
/// WARNING: assumes that seqset for ma is derived from gss_typ.
// eventually change so that uses gss_typ& gss = gssCMSA(cma);
{
	Int4	lenM,LenStr,AlignLen,J,x,*s;
	char	**alignment;
	st_type	S=SitesCMSA(cma);
	ss_type	data=DataCMSA(cma);
	a_type	A=SeqSetA(data);
	e_type	sgE,E;
	fm_type	*M=ModelsCMSA(cma);
	Int4	site,number,buffer_size;
	char	*ps;
	BooLean	**null;
	gss_typ     *gss=gssCMSA(cma);

	print_error("PutFastaAlnCMSA( ) not working properly...exiting.");
	assert(gss->Gapped());
	number=gss->NumSeq();
	assert(nBlksCMSA(cma) == 1);
	x = nColsFModel(M[1]);
	buffer_size = gss->MaxTrueSqLen();
	for(J=1; J <= number; J++) buffer_size+=gss->NumIns(J)+gss->NumDel(J);
	s = new Int4[number+3];
	null = NullCMSA(cma);

	//// 1. Get locations of sites in gapped sequences from cmsa.
	lenM = SiteLen(1,S);
	NEWP(alignment,number+3,char);
	for(J=1; J <= number; J++){ 
		sgE = gss->FakeSeq(J);
		s[J] =SitePos(1,J,1,S);
// fprintf(stderr,"s[%d] = %d\n",J,s[J]);
		NEW(alignment[J],buffer_size+3,char); 
	   	//// 2. Get corresponding gapped alignable segments.
		LenStr=gss->Region(J,alignment[J],s[J],lenM);
	}
	//// 3. Align gapped sequences.
	AlignLen=cmsa_fix_alignment(lenM,number,alignment,null[1]);

	//// 4. Print gapped alignment.
	for(J=1; J <= number; J++){ 
		sgE = gss->FakeSeq(J);
		site = gss->TrueSite(J, s[J]);
		ps=alignment[J];
		E = gss->TrueSeq(J);
        	fprintf(fp,">"); PutSeqInfo2(fp,E);
		for(Int4 i=0; alignment[J][i]; i++){
		   char r=alignment[J][i];
		   if(r == '.') r='-';
		   fprintf(fp,"%c",r);
		}
		fprintf(fp,"\n\n");
		for(ps++; *ps != 0; ) { if(isalpha(*ps)) site++; ps++; }
	} fflush(fp);
	for(J=0; J <= number; J++){ free(alignment[J]); } free(alignment);
	// alignment[0] allocated above....
	free(null[1]); free(null);
	delete []s;
}
#endif

char	*AddInsertToOperationArray(Int4 start, Int4 end, char *operation)
// ========== Add an insertion to an operational array. ===========
{
	
	char state,*new_operation=0;
	Int4 j,o,no,column;
	Int4 trace_length=strlen(operation);
	NEW(new_operation,trace_length+5,char);
	new_operation[0]='E'; no=1;
	for(o=j=1,column=1,state='E'; operation[o] != 'E'; o++){
          switch(operation[o]){
            case 'M': 
            case 'm': 
		if(column >= start && column <=end){
			   new_operation[no]='I'; 
		} else new_operation[no]=operation[o];
		no++; j++; column++; break;
            case 'D': 
            case 'd': // deletion in sequence relative to profile.
		if(column >= start && column <=end){
			   // do nothing in new_operation; 
		} else { new_operation[no]=operation[o]; no++; }
                column++; break;
            case 'i': // insert is between profile blocks;
		new_operation[no]=operation[o]; no++; 
		j++; break;
            case 'I': // Insert ('-') within a profile block; delete from seq.
		new_operation[no]=operation[o]; no++; 
		j++; break;
            default:
            // fprintf(stderr,"operations = %s\n",operation);
            print_error("operation( ): input error"); break;
          }  state=operation[o];
       	}
	new_operation[no]='E'; no++; new_operation[no]=0;
	// printf("new operation: %s\n",new_operation);
	return new_operation;
}

st_type	FuseElementsSites3(Int4 x, Int4 maxlen, st_type S,a_type A)
/*****************************************************************
 NEW: 5/1/04 - Andy Neuwald
 if element x and x+1 are <= maxlen then fuse them into one 
 element; create and return the resultant object S2.
 WARNING: NOT TESTED.
 Use only if number of blocks == 2 (for now).
 *****************************************************************/
{
    Int4        t,n,s,s2,N,T,T2,*len_elem,len,len2,newlen;
    st_type	S2;

    // print_error("FuseElementsSites3 not yet implemented");
    // 0. Test for input errors
 if(x==1 && nTypeSites(S) == 2){
    len = SiteLen(x,S); len2 = SiteLen(x+1,S);
    NEW(len_elem, 3, Int4);
    len_elem[1]=len+len2; 
    S2=MakeSites(1,len_elem,S->gss); // preserves S->gss in S2 with indels
 } else { assert(x==1 && nTypeSites(S) == 2); }
    free(len_elem);

    // 2. Apply HideInsert to gsq in gss structure...
    // Find start and end for gapping later...
    char *new_operation;
    Int4 pos[4],newpos[4];
    N=NSeqsSites(S2);
    for(Int4 sq=1; sq <= N; sq++){
	Int4 start=EndSitePos(x,sq,1,S) + 1;
        Int4 end=SitePos(x+1,sq,1,S) - 1;
	assert(PosTSites(1,sq,pos,S));
	if(start <= end){  // there is an insert here...
	   char *operation=S2->gss.Operation(sq);
	   Int4 trace_length=strlen(operation);
	   new_operation=AddInsertToOperationArray(start,end,operation);
	   // trace_length=strlen(new_operation);

	   e_type E = S2->gss.TrueSeq(sq);
	   Int4 strt=S2->gss.TruePos(sq,pos[1]);
	   gsq_typ *gsq; gsq = new gsq_typ[1];
	   assert(S2->gss.LeftFlank() == 0 && S2->gss.RightFlank() == 0);
	   gsq->initialize(S2->gss.LeftFlank(),S2->gss.RightFlank(),
                new_operation,trace_length,strt,E,newpos);
	   fprintf(stderr,"pos=%d; newpos=%d\n",strt,newpos[1]);
	   gsq->Put(stdout,A);
	   ReplaceSeqSites(sq,gsq,S2); // removes all sites in S.
	} AddSite(1,sq,pos[1],S2);
    } return S2;
}

cma_typ	SimpleFuseBlksCMSA(Int4 x, Int4 maxlen, cma_typ L)
// Fuse block x and x+1 by deleting block x+1 and lengthing block x
// without modifying underlying gapped sequences...
{
    cma_typ     ma;
    Int4        sq,j,t2,t,N,T,T2;
    st_type     S,S2=NULL;
    fm_type     *model;
    BooLean     **null;
    Int4	*site_len;

    S=SitesCMSA(L); T = nBlksCMSA(L); T2 = T-1;
    if(T < 2 || x < 1 || x >= T) return NULL;
    // NEW(site_pos,T+5,Int4);
    N = NumSeqsCMSA(L);
    S2 = FuseElementsSites3(x, maxlen, S,AlphabetCMSA(L));
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

char    *AddInsertToOperationArray2(Int4 start, Int4 end, char *operation)
// ========== Add an insertion to an operational array. ===========
{

        char state,*new_operation=0;
        Int4 o,no,column;

        Int4 trace_length=strlen(operation);
        NEW(new_operation,trace_length+5,char);
        new_operation[0]='E'; no=1;
	for(no=o=1; isalpha(operation[o]); o++){
		char c = operation[o];
		if(c == 'M' || c == 'D'){ break; }
		else { new_operation[no] = c; no++; }
	}
	if(!isalpha(operation[o])){
	    fprintf(stderr,"operation[%d]='%c'=%d\n",o,operation[o],operation[o]);
	    fprintf(stderr,"operation=%s\n",operation);
	    fprintf(stderr,"start=%d; end=%d\n",start,end);
	    fprintf(stderr,"new operation=%s\n",new_operation);
	    print_error("AddInsertToOperationArray2() input error");
	}
        // Int4 InsSize=end-start+1;
        // NEW(new_operation,trace_length+InsSize+5,char);
        for(column=1,state='E'; operation[o] != 'E'; o++){
          switch(operation[o]){
            case 'M':
            case 'm':
                if(column >= start && column <=end){
                           new_operation[no]='I';
                } else new_operation[no]=operation[o];
                no++; column++; break;
            case 'D':
            case 'd': // deletion in sequence relative to profile.
                if(column >= start && column <=end){
                           // do nothing in new_operation;
                } else { new_operation[no]=operation[o]; no++; }
                column++; break;
            case 'i': // insert is between profile blocks;
                new_operation[no]=operation[o]; no++;
                break;
            case 'I': // Insert ('-') within a profile block; delete from seq.
                new_operation[no]=operation[o]; no++;
                break;
            default:
            // fprintf(stderr,"operations = %s\n",operation);
            print_error("operation( ): input error"); break;
          }  state=operation[o];
        }
        new_operation[no]='E'; no++; new_operation[no]=0;
        // printf("new operation: %s\n",new_operation);
        return new_operation;
}

Int4	IronOutOperation(char *operation)
/*************************************************************************************
 Just switch ddddIIII to IIIIdddd to fix program...so can work with cma_gblastpgp...

    ddddIIIIm
 d: 01233333
 i: 00001234

  Simply move insertions to other side of deletions; don't do this 
 *************************************************************************************/
{
	Int4 length=strlen(operation);
	// char *new_operation;
	char state,last=' ';
	Int4 i,j,strt_i=0,strt_d=0,end;
	Int4 num_i=0,num_d=0,num_fix=0;

	assert(operation[0] == 'E');
	// NEW(new_operation,length+3,char);
	if(operation[1] == 'D'){ for(i=2; i < length; i++) if(operation[i] != 'd') break; }
	else i=1;
	last=operation[i-1];
	for( ; i < length; i++){
		state = operation[i];
		if(state != last){	// then 
		   if(state=='E'){
		     // if "..dE" or "..mE" then do nothing
		     // "..IE"
		     if(num_i > 0) print_error("IronOutOperation() input error 0");
		   } else if(state=='I' && last == 'd'){	// then start collect iron out data.
			// ...dI...
			num_i=1; strt_i=i;
		   } else if(last == 'I' && num_d > 0){	// then iron out...
			num_fix++;
#if 1
			// ..dI..Im or ..dI..Id
			end = strt_d + num_i;
			for(j=strt_d; j <  end; j++) operation[j]='I';
			end += num_d;
			for( ; j <  end; j++) operation[j]='d';
		        num_i=num_d=0; strt_i=0; strt_d=0;
			if(state=='d'){
				strt_d=i; num_d=1;
			}
#else		// change inserts into matches and remove deletions...
#endif
#if 0	// dddDIII
		   } else if(state == 'D' && operation[i+1] == 'I'){
			   operation[i]='M'; i++;
			   operation[i]='d'; last='M'; state='d'; continue;
#endif
		   } else {	// state 
		     num_i=num_d=0; strt_i=strt_d=0;
		     switch(state){
			case 'I': strt_i=i; num_i=1; break;
			case 'm': case 'M': break;
			case 'D': {
			    fprintf(stderr,"i=%d\n",i);
			    fprintf(stderr,"operation = %s\n",operation);
			    print_error("IronOutOperation() input error 1"); 
			  }break;
			case 'd': strt_d=i; num_d=1; break;
		     }
		   }
		} else {
		   switch(state){
			case 'I': num_i++; break;
			case 'm': case 'M': break;
			case 'd': num_d++; break;
			default: print_error("IronOutOperation() input error 2"); break;
		   }
		} last=state;
	} return num_fix;
}


cma_typ	ColumnsToInsertCMSA2(cma_typ cma,Int4 start_ins, Int4 end_ins)
// Convert a region within a block into an insert (lower case in cma file).
{
#if 0	// replaced with newer routine (delete below after testing...)
	cma_typ rcma=ConvertColsToInsertsCMSA(cma, 1, start_ins, end_ins); 
	NilCMSA(cma); *cma=*rcma; return rcma;
#else
	Int4	j,i,N=NumSeqsCMSA(cma);
	a_type	A=AlphabetCMSA(cma);
	char	**Operation;
	Int4	trace_length,*TruePos;
	Int4	newpos[3];

	assert(nBlksCMSA(cma)==1);
	assert(LengthCMSA(1,cma) > end_ins && start_ins <= end_ins);
	gss_typ *gss=gssCMSA(cma);
	NEWP(Operation,NumSeqsCMSA(cma)+5,char);
	NEW(TruePos,NumSeqsCMSA(cma)+5,Int4);
	for(Int4 sq=1; sq <= NumSeqsCMSA(cma); sq++){
	        // ========== X. Get operational array for sequence. ===========
#if 0
		gsq_typ *gsq=gss->GetGSQ(sq);
		char *operation=gsq->Operation3(0,A);
		printf("old operation(%d): %s\n",sq,operation);
		if(sq==23 || sq==21) gsq->Put(stderr,A);
#else
		char *operation=gss->Operation(sq);
#endif
		// WARNING: This maps fake to real seq. NOT seq to alignment model!!!

		// printf("old operation(%d): %s\n",sq,operation);
		TruePos[sq]=TruePosCMSA(sq,1,cma);
		if(TruePos[sq] > 1){	// Fixes problem with N-terminal extension...
			// gsq_typ *gsq0=gss->GetGSQ(sq); gsq0->Put(stderr,A); 
			operation[1] = tolower(operation[1]);
			operation[TruePos[sq]] = toupper(operation[TruePos[sq]]);
			// printf("old operation(%d): %s\n",sq,operation);
			TruePos[sq] =1;
		} 
		char *new_operation=AddInsertToOperationArray2(start_ins,end_ins,operation);
#if 0
		fprintf(stderr,"%d: operation=%s\n",sq,operation);
		fprintf(stderr,"%d: operation=%s\n",sq,new_operation);
#elif 0
		// printf("%d: operation=%s\n",sq,operation);
		Int4	wrinkles=IronOutOperation(new_operation);
		if(wrinkles > 0){
			fprintf(stderr,"%d: %d wrinkles\n",sq,wrinkles);
			fprintf(stderr,"old operation(%d): %s\n",sq,operation);
			fprintf(stderr,"new operation(%d): %s\n",sq,new_operation);
		}
#endif
		Operation[sq]=new_operation;
		VacateSitesCMSA(sq,cma);
		// free(operation);
	}

	// ========== X. modify operational array to add insert. ===========
	for(i=start_ins; i<=end_ins; i++){ 	// shorten model...
		Int4 lemon = LengthCMSA(1,cma); RmColumnMSA(1,lemon,cma);
	}

	// ========== X. Create a gapped sequence. ===========
	for(Int4 sq=1; sq <= NumSeqsCMSA(cma); sq++){
		gsq_typ *gsq0=gss->GetGSQ(sq); 
		// VacateSitesCMSA(sq,cma);
		Int4 start=TruePos[sq];
		trace_length=strlen(Operation[sq]);
		e_type E = gss->TrueSeq(sq);
		gsq_typ *gsq; gsq = new gsq_typ[1];
		// assert(gss->LeftFlank() == 0 && gss->RightFlank() == 0);
#if 1	// Bug fix for deletions at ends...
		if(start > 0 && Operation[sq][1]=='D') start++;
#endif
		// gss->OverHangN(sq) + gss->OverHangC(sq);
		gsq->initialize(gsq0->OverHangN(),gsq0->OverHangC(),
                  Operation[sq],trace_length,start,E,newpos);
		// fprintf(stderr,"pos=%d; newpos=%d\n",TruePos[sq],newpos[1]);
		ReplaceCMSA(sq,gsq,cma); // replace sequence s in CMSA & fmodel.
		AddSiteCMSA(1,sq,newpos[1], cma);
		free(Operation[sq]);
	} free(Operation);
	return cma;
#endif
}

BooLean	ColumnsToInsertCMSA(cma_typ cma,Int4 start_ins, Int4 end_ins)
// Convert a region within a block into an insert (lower case in cma file).
{
#if 0	// replaced with newer routine (delete below after testing...)
	cma_typ rcma=ConvertColsToInsertsCMSA(cma, 1, start_ins, end_ins); 
	NilCMSA(cma); *cma=*rcma; return TRUE;
#else
	Int4	j,i;
	a_type	A=AlphabetCMSA(cma);
	char	**Operation;
	Int4	trace_length,*TruePos;
	Int4	pos[3],newpos[3];

	assert(nBlksCMSA(cma)==1);
	assert(LengthCMSA(1,cma) > end_ins && start_ins <= end_ins);
	gss_typ *gss=gssCMSA(cma);
	NEWP(Operation,NumSeqsCMSA(cma)+5,char);
	NEW(TruePos,NumSeqsCMSA(cma)+5,Int4);
	for(Int4 sq=1; sq <= NumSeqsCMSA(cma); sq++){
	        // ========== X. Get operational array for sequence. ===========
		char *operation=gss->Operation(sq);
		// printf("old operation(%d): %s\n",sq,operation);
		// PosSiteCMSA(1,sq,pos,cma); 
		// fprintf(stdout,"%d: pos=%d\n",sq,pos[1]);
		TruePos[sq]=TruePosCMSA(sq,1,cma);
#if 0	// fix bug...
		char *new_operation=operation;
#else	// regular routine
		char *new_operation=AddInsertToOperationArray(start_ins,end_ins,operation);
#endif
		// printf("new operation(%d): %s\n",sq,new_operation);
		Operation[sq]=new_operation;
	}

	// ========== X. modify operational array to add insert. ===========
	for(i=start_ins; i<=end_ins; i++){
		Int4 lemon = LengthCMSA(1,cma);
		RmColumnMSA(1,lemon,cma);
	}

	// ========== X. Create a gapped sequence. ===========
	for(Int4 sq=1; sq <= NumSeqsCMSA(cma); sq++){
		VacateSitesCMSA(sq,cma);
		Int4 start=TruePos[sq];
		// make sure that start is correct for gsq->initialize();
		trace_length=strlen(Operation[sq]);
		e_type E = gss->TrueSeq(sq);
		gsq_typ *gsq; gsq = new gsq_typ[1];
		assert(gss->LeftFlank() == 0 && gss->RightFlank() == 0);
#if 1	// Bug fix for deletions at ends...
		if(start > 0 && Operation[sq][1]=='D') start++;
#endif
		// WARNING: need to fix this as flanking sequence on the left 
		// "{(AFGV)..." is causing fatal problems....
#if 0	// Diagnostic to find out how to deal with overhangs...
		fprintf(stderr,"operation: %s\n",Operation[sq]);
#endif
#if 1	// Bug fix...
		gsq->initialize(LenSeq(E),LenSeq(E),Operation[sq],trace_length,start,E,newpos);
#else
		gsq->initialize(gss->LeftFlank(),gss->RightFlank(),
                  Operation[sq],trace_length,start,E,newpos);
#endif
		// fprintf(stderr,"pos=%d; newpos=%d\n",TruePos[sq],newpos[1]);
		// gsq->Put(stdout,A);
		ReplaceCMSA(sq,gsq,cma); // replace sequence sq in CMSA & fmodel.
		AddSiteCMSA(1,sq,newpos[1],cma); // afn: 2-21-08
		// AddSiteCMSA(1,sq,start,cma);	// afn: 2-15-08
	}
        // sprintf(str,"%s.insrt",argv[1]);
        // WriteMtfCMSA(str, cma, NULL);
	return TRUE;
#endif
}

#if 0	// afn: 1/26/08  moved back to tweakcma
Int4	PutClusterOfCMSA(char *name, Int4 percent_ident,Int4 min_size,
	BooLean	IncludeFirst,cma_typ cma)
// print out separate cmas for each cluster from the input cma.
{
	ds_type sets;
	Int4	b,i,j,k,N=NumSeqsCMSA(cma);
	Int4	score,nblk=nBlksCMSA(cma);
	Int4	s,si,sj,pos[3],seti,setj;
	ss_type	data=DataCMSA(cma);
        a_type  A=AlphabetCMSA(cma);
	unsigned char	*isq,*jsq;
#if 0
	Int4	NumPhyla=0;
	Int4    *phyla=GetPhylaSeqSet(stderr, &NumPhyla, data);
#endif

	assert(percent_ident > 0 && percent_ident <= 100);
	// 1. Cluster sequences into sets at the percent identity cutoff.
	Int4 total = TotalLenCMSA(cma);
	sets = DSets(N);
	for(i = 1; i < N; i++){
	  isq = SeqPtrCMSA(i,cma);
	  seti=findDSets(i,sets);
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
	      } score = (Int4) floor(((double)score*100.0/(double)total) +0.5);
	      if(score >= percent_ident) seti=linkDSets(seti,setj,sets);
	     }
	  }
	}
	// 2. Within each set pick the sequence with the highest profile score.
	BooLean	**skipseq;
	Int4	*size;
	NEWP(skipseq, NumSeqsCMSA(cma)+2, BooLean);
	NEW(size, NumSeqsCMSA(cma)+2, Int4);
        for(i=1; i <= N; i++){
	  seti=findDSets(i,sets);
	  if(skipseq[seti]==0) {
		NEW(skipseq[seti],NumSeqsCMSA(cma)+2, BooLean); 
        	for(j=1; j <= N; j++) skipseq[seti][j] = TRUE;
	  } skipseq[seti][i]=FALSE; size[seti]++;
	} 
#if 1	// get probabilities...
	double *prob=0,best_prob;
	Int4	best_sq;
	NEW(prob,N+2,double); NEW(skipseq[0],N+2, BooLean); 
        for(i=1; i <= N; i++){
		skipseq[0][i]=TRUE; 
		prob[i]=GetTotalGappedProbCMSA(i,cma);
	}
        for(i=1; i <= N; i++){
	  if(skipseq[i]==0) continue; // no set i; this sequence in another set.
	  best_prob=-DBL_MAX;
	  for(best_sq=0,j=1; j <= N; j++){
		if(skipseq[i][j] == FALSE){	// indicates sequence j is in set i.
		  if(prob[j] >= best_prob){
			best_prob=prob[j]; best_sq=j;
		  }
		}
	  }
	  assert(skipseq[0][best_sq]==TRUE); // must be true at this point.
	  skipseq[0][best_sq]=FALSE;	// retain best sequence in set i.
	} free(prob);
#endif
	// 3. output the cma files;
	char str[100];
	sprintf(str,"_0.cma");
        FILE *fp = open_file(name,str,"w");
	sprintf(str,"cluster0"); RenameCMSA(str,cma);
	PutSelectCMSA(fp,skipseq[0],cma); fclose(fp);
	free(skipseq[0]);
	for(s=0,i=1; i <= N; i++){
	   if(skipseq[i] && size[i] >= min_size){
		if(IncludeFirst) skipseq[i][1]=FALSE; // always include the concensus seq.
		s++;
		fprintf(stderr,"set %d(%d): size = %d\n",s,i,size[i]);
		sprintf(str,"_%d.cma",s);
           	fp = open_file(name,str,"w");
		sprintf(str,"cluster%d",s); RenameCMSA(str,cma);
	   	PutSelectCMSA(fp,skipseq[i],cma); 
		fclose(fp); free(skipseq[i]); 
	   }
	} NilDSets(sets); free(skipseq);
	return s;
}
#endif

double	JunLiuHMM_PenaltyCMA(FILE *fp, cma_typ cma)
// Test Jun Liu's statistical model for GISMO...
{
	gss_typ *gss=gssCMSA(cma);
	assert(gss->Gapped());
	double penalty=0.0;
	UInt4 Nmm,Nmd,Nmi,Nm;
	UInt4 Nii,Ndd,Nim,Ndm;
	UInt4 Nid,Ndi,Nsd,Nsm;
	Int4 ae,be,ao,bo,n0,n1,n2;

	Nmm=Nmi=Nmd=Nii=Nim=Ndd=Ndm=Nid=Ndi=Nsd=Nsm=0;
	Int4 *pos;
	Int4	*lens = LengthsCMSA(cma);
	NEW(pos,nBlksCMSA(cma)+3,Int4);
	for(Int4 sq=1; sq <= NumSeqsCMSA(cma); sq++){
	   	assert(PosSites(sq,pos,SitesCMSA(cma)) == nBlksCMSA(cma));
		gsq_typ *gsq=gss->GetGSQ(sq);
		gsq->FindIndels(nBlksCMSA(cma),pos,lens,Nmm,Nmi,Nmd,Nii,Nim,Ndd,Ndm,Ndi,Nid,Nsd,Nsm);
	} free(pos);

	n0=10000;
	ao=10;	// one insertion in 1000 residues.
	bo=10;	// one deletion every 1000 residues 
	n1=n2=100;
	ae=20;	// L1 = 20/100 --> 5 = average insert length
	be=50;	// L2 = 50/100 --> 2 = average deletion length

	Nm=Nmm+Nmi+Nmd;

	penalty=0.0;
	// numerator...
	penalty += lngamma((double)(Nmi+ao));
	penalty += lngamma((double)(Nmd+bo));
	penalty += lngamma((double)(Nmm+n0-ao-bo));
	penalty += lngamma((double)(n0));
	// denominator...
	penalty -= lngamma((double)(Nm+n0));
	penalty -= lngamma((double)(ao));
	penalty -= lngamma((double)(bo));
	penalty -= lngamma((double)(n0-ao-bo));

	// numerator...
	penalty += lngamma((double)(Nii+ae));
	penalty += lngamma((double)(Nim+n1-ae));
	penalty += lngamma((double)(n1));
	// denominator...
	penalty -= lngamma((double)(Nim+Nii+n1));
	penalty -= lngamma((double)(ae));
	penalty -= lngamma((double)(n1-ae));

	// numerator...
	penalty += lngamma((double)(Ndd+be));
	penalty += lngamma((double)(Ndm+n2-be));
	penalty += lngamma((double)(n2));
	// denominator...
	penalty -= lngamma((double)(Ndd+Ndm+n2));
	penalty -= lngamma((double)(be));
	penalty -= lngamma((double)(n2-be));
	if(fp){
	 double Bo,Be,Ao,Ae;
	 static Int4 Seed=0;
	 if(Seed==0) {
           Seed=-Random();  // need to initialize with a negative number.
           fprintf(stderr,"INITIALIZING SEED (%d)\n",Seed);
         }
	for(Int4 iter=1; 0 && iter < 10; iter++){
	 Bo=betadev(Nmd+Nid+1,Nmm+Nmi+Nii+Nim+1,&Seed);
	 Be=betadev(Ndd,Ndm,&Seed);
	 Ao=(1-Bo)*betadev(Nmi,Nmm,&Seed);
	 Ae=(1-Bo)*betadev(Nii,Nim,&Seed);
	 fprintf(stdout,"Ao = %g; Ae = %g; Bo = %g; Be = %g\n",Ao,Ae,Bo,Be);
	}

	 fprintf(stdout,"Nmm = %d; Nmi = %d; Nmd = %d; Nm=%d; Nii = %d; Nim = %d\n",
			Nmm,Nmi,Nmd,Nm,Nii,Nim);
	 fprintf(stdout,"Ndd = %d; Ndm = %d; Ndi = %d; Nid = %d\n",Ndd,Ndm,Ndi,Nid);
	 fprintf(stderr,"penalty=%g; Nmm = %d\n",penalty,Nmm);
	 fprintf(stderr,"current indel_penalty = %g\n",-IndelPenaltySeqSet(DataCMSA(cma)));
	}
	return penalty;
}

void    PutTableOfCMSA(FILE *fptr, cma_typ cma)
{
    char	r,c,str[200],str2[10];
    a_type	A=AlphabetCMSA(cma);
    BooLean	*use;
    Int4	*pos,t,i,j,n,s,length,ntyp,N;
    char	*gnull,**gseq;
    unsigned char	**seq;

    FILE    *fp=tmpfile();
    PutAlnCMSA(fp,cma); rewind(fp);
    sma_typ MA=ReadSMA(fp); fclose(fp);
    ntyp = ntypSMA(MA);
    N=nseqSMA(MA);
    assert(N <= 100); // don't go beyond 50 sequences...
    NEW(pos,nseqSMA(MA)+3,Int4);
    NEWP(gseq,nseqSMA(MA)+3,char);
    NEWP(seq,nseqSMA(MA)+3,unsigned char);
    // assert(ntyp == 1); t=1;
    for(t=1; t <= ntyp; t++){
        length = lengthSMA(t,MA);
        gnull = gnullSMA(t,MA);
	for(n=1; n<= nseqSMA(MA); n++){
	   if(t==1){
		sprintf(str,"%s",seqidSMA(n,MA));
		strncpy(str2,str,5); str2[5]=0;
		fprintf(fptr,"%-5s ",str2);
		// fprintf(fptr,"%-5s %4d ",str2,startSMA(1,n,MA));
	   } else {
		fprintf(fptr,"mtf%d ",t);
	   }
		pos[n] = startSMA(t,n,MA);
	        gseq[n] = gseqSMA(t,n,MA);
		if(gseq[n][0] == '-') pos[n]++;
		seq[n] = seqSMA(t,n,MA);
	} fprintf(fptr,"\n");
	for(i=j=0; j < lengthSMA(t,MA); i++){
	  if(gnull[i] == '_'){  // deleted in concensus sequence...
	    if(ntyp == 1){
	      for(N=0,n=1; n<= nseqSMA(MA); n++) if(isalpha(gseq[n][i])) N++;
	    } else N=1000;
	    {	// don't print a delete row unless at least 2 seqs with residues.
	      for(n=1; n<= nseqSMA(MA); n++){
		c = gseq[n][i];
		if(isalpha(c)){ if(N > 2) fprintf(fptr,"%c%-4d ",c,pos[n]); pos[n]++; }
		else if(N > 2) fprintf(fptr,"%c     ",c);
	      } if(N > 2) fprintf(fptr,"\n");
	    }
	  } else {
	    j++; 
	    for(n=1; n<= nseqSMA(MA); n++){
		r = seq[n][j]; 
		if(r==0 && gseq[n][i] == '-'){ c = '-'; fprintf(fptr,"%c     ",c); }
		else { c=AlphaChar(r,A); fprintf(fptr,"%c%-4d ",c,pos[n]); pos[n]++; }
	    } fprintf(fptr,"\n");
	  }
	    // fprintf(fptr," %d", endSMA(t,n,MA));
	} fprintf(fptr,"\n");
    } NilSMA(MA); free(pos); free(gseq); free(seq);
}

void	ShuffleColumnsCMA(cma_typ cma)
// Shuffle all of the columns between the sequences for cma
// This routine is for testing Jun's BPPS procedure.
// We would expect the BPPS procedure not to work for these
// shuffled cma's.
{
	assert(nBlksCMSA(cma) == 1);
	unsigned char	*res,r;
	gss_typ *gss=gssCMSA(cma); assert(gss != 0);
	Int4	i,m,n,N=NumSeqsCMSA(cma),*Pos,number=gss->NumSeq();
	e_type	fE,tE;

	assert(number==N); 
	// 0. Go column by column through the alignment.
	for(Int4 col=1; col <=LengthCMSA(1,cma); col++){
		NEW(res,N+3,unsigned char); NEW(Pos,N+3,Int4);
		// 1. store residues and positions for each sequence.
		dh_type dH=dheap(number+2,4);
        	for(m=0,n=1; n <= number; n++) {
			// Pos[n]=TruePosCMSA(n,col,cma);
			// r=ResidueCMSA(1,n,Pos[n],cma);
		   r=ResidueCMSA(1,n,col,cma);	// fake seq?
		   if(r){
			m++; res[m]=r; 
		        insrtHeap(n,(keytyp) Random(),dH);
		   }
		} assert(ItemsInHeap(dH) == m);
		// 2. randomly sort residues on a heap.
		for(m=0; (n=delminHeap(dH)) != 0; ){
		    Int4 site=FakeToRealCMA(n,col,cma);
		    gsq_typ *gsq=gsqCMSA(n,cma);
		    m++; EqSeq(site,res[m],gsq->TrueSeq( )); 
		} Nildheap(dH);
		free(res); free(Pos);
	} 
}

cma_typ	InsertColumnsCMSA(cma_typ cma, Int4 Blk, Int4 start_ins, Int4 length)
// Insert new columns after start_ins (lower to upper case in cma file).
{
	a_type  AB=AlphabetCMSA(cma);
	Int4	blk,*len,numfix,Length;
	assert(Blk > 0 && Blk <= nBlksCMSA(cma));
	if(length < 0){ 
		assert(start_ins <= LengthCMSA(Blk,cma) && start_ins > 1); Length=-length;
	} else {
		assert(start_ins < LengthCMSA(Blk,cma) && start_ins >= 1); Length=length;
	}
	gss_typ *gss=gssCMSA(cma);
	NEW(len,nBlksCMSA(cma) +3, Int4);
	for(blk=1; blk <= nBlksCMSA(cma); blk++) len[blk]=LengthCMSA(blk,cma);
	len[Blk] = len[Blk]+Length;
   	cma_typ rcma=EmptyCMSA(nBlksCMSA(cma),len,TrueDataCMSA(cma),gss->GapOpen(),
		gss->GapExtend(),PerNatsCMSA(cma),0,0);
	for(Int4 sq=1; sq <= NumSeqsCMSA(cma); sq++){
	        // ========== X. Get operational array for sequence. ===========
		gsq_typ *gsq=gsqCMSA(sq,cma);
#if 1	// subcma operations; check for vacant site ...afn: 11-26-2014.
        	if(!SeqIsInCMSA(sq,cma)){
			if(gsq==0){
			   gsq_typ *gsq0 = new gsq_typ[1];
                	   gsq0->initialize(TrueSeqCMSA(sq,cma));    // copies positions to sites array.
			   ReplaceCMSA(sq,gsq0,rcma);
			} continue;
		}
#else
        	assert(SeqIsInCMSA(sq,cma));
#endif
		Int4	*sites=GetPosSitesCMSA(sq,cma);
#if 0	// debug
		gsq->Put(stderr,AB);
		gsq->Put_cma_format(stderr,sq,nBlksCMSA(cma),sites,LengthsCMSA(cma),AB);
		char *operation=gss->Operation(sq);
        	fprintf(stderr,"operation(%d): %s\n",sq,operation);
#endif
		gsq_typ	*gsq0=gsq->InsertColumns(Blk, start_ins, length,
					nBlksCMSA(cma),LengthsCMSA(cma),sites);
		ReplaceCMSA(sq,gsq0,rcma); // replace sequence s in CMSA & fmodel.
		for(blk=1; blk <= nBlksCMSA(rcma); blk++) AddSiteCMSA(blk,sq,sites[blk],rcma);
		free(sites); 
#if 0	// debug
		sites=GetPosSitesCMSA(sq,rcma);
		gsq0->Put_cma_format(stderr,sq,nBlksCMSA(rcma),sites,LengthsCMSA(rcma),AB);
		free(sites); 
#endif
	} free(len);
        // sprintf(str,"%s.insrt",argv[1]); WriteMtfCMSA(str, cma, NULL);
	return rcma;
}

cma_typ InsertColumnsCMA(Int4 start, Int4 N, cma_typ cma)
// add n new columns just to the right of start ...
{
        Int4 min_len=1;
	BooLean	InsertAtEnd=FALSE;
        assert(nBlksCMSA(cma) == 1);
        // e_type qE=SeqSetE(1,TrueDataCMSA(cma));
        // start = start - OffSetSeq(qE);

        // 1. Split block right after start_ins
	cma_typ cma0=0,cma1=0;
#if 1
	if(start == LengthCMSA(1,cma)){	// insert at end of alignment...
        	cma0=CopyCMSA(cma);
		InsertAtEnd=TRUE;
	} else {
        	cma0=SplitBlkCMSA(1,start,min_len,cma);
	}
#else
        cma_typ cma0=SplitBlkCMSA(1,start,min_len,cma);
#endif
        if(cma0 == 0){
          FILE *fp = open_file("debug_error",".cma","w"); PutCMSA(fp,cma); fclose(fp);
          assert(cma0 != 0);
          print_error("Error in InsertColumnsCMA()");
        }
        // WriteCMSA("debug0.cma", cma0);

        // 2. Append n columns to first block right after start of insert
        BooLean right=TRUE;
        for(Int4 col=1; col<= N; col++){
                if(InsertColCMSA(1,right,cma0) == FALSE)
                        fprintf(stderr,"InsertColCMSA() inserted gaps\n");
        }
        // WriteCMSA("debug1.cma", cma0);

        // 3. Fuse blocks again.
	if(InsertAtEnd){
		cma1=CopyCMSA(cma0);
	} else {
        	cma1=SimpleFuseBlksCMSA(1, 10000, cma0);
	}
        if(cma1 == 0){
          FILE *fp = open_file("debug_error",".cma","w"); PutCMSA(fp,cma0); fclose(fp);
	  print_error("Fuse error in InsertColumnsCMA()");
	}
        // WriteCMSA("debug2.cma", cma1);
        // NilCMSA(cma);
        NilCMSA(cma0);
        return cma1;
}

static Int4 MAX_ALIGN_LENGTH=1000000;

cma_typ	RemoveOverhangsCMSA(cma_typ in_cma, BooLean AddX2Ends)
{
	FILE *ofp, *ifp;
	a_type AB=AlphabetCMSA(in_cma);
	ifp = tmpfile(); PutFastaCMSA(ifp,in_cma); rewind(ifp);

	unsigned long    num=0;
	char	*Seq; NEW(Seq,MAX_ALIGN_LENGTH +3, char);
	char	*MasterSeq; NEW(MasterSeq,MAX_ALIGN_LENGTH + 3,char);
	char	*str; NEW(str,MAX_ALIGN_LENGTH+3,char);
	unsigned long Begin=ULONG_MAX,End=0,Length;
	e_type	DummySeq=0;
	long	SeqOffSet=0,SeqExtend=0;
	char	Kingdom,Phylum[100];
        Int4	i,arg,strleng;

   ofp=tmpfile(); // ofp=stderr;
   while(fgets(str,(MAX_ALIGN_LENGTH),ifp)){
	strleng=strlen(str)-1;
	if(strleng > MAX_ALIGN_LENGTH) print_error("RemoveOverhangsCMSA(): sequence too long!!!");
	if(str[0] == '>'){
		Int4 i,os,x;
		char *str1,*str2,*str3,*str4;
		unsigned char seq[5]={1,1,1,1,1,};
		Int4 len=5;
		num++;
		if(num==1){
		   fprintf(ofp,"[0_(1)=%s(%d){go=10000,gx=2000,pn=1000.0,lf=0,rf=0}:\n",
			NameCMSA(in_cma),NumSeqsCMSA(in_cma));
		}
		SeqOffSet=0; SeqExtend=0;
		e_type  E=MkSeq(str, len, seq);
		DummySeq=E;
		Kingdom=KingdomSeq(E);
		if(PhylumSeq(E)) strncpy(Phylum,PhylumSeq(E),50);
		else strncpy(Phylum,"unknown",50);
		SeqOffSet=OffSetSeq(E);
		// ExtendSeq = CtermExtendSeq(E);
	} else if(isspace(str[0])){
		// ignore;
	} else {	// A-Z a-z '-'.
		unsigned long i,j,match,insert,deletion;
		char c;
		if(num < 1) { fprintf(stderr,"Input error 1\n"); exit(1); }
		if(AddX2Ends){
			c=str[0];
			if(c =='-' || c =='.') str[0]='X';
			i=strleng-1; c=str[i];
			if(c =='-' || c =='.') str[i]='X';
		}
		if(num == 1){	// then get Master state information
		   for(Length=i=0; i < strleng; i++){
			if(isalpha(str[i])){
				Length++;
				if(Begin==ULONG_MAX) Begin=i;
				End=i;
				MasterSeq[i]='m'; // match state.
			} else if(str[i] == '.' || str[i] == '-'){
				MasterSeq[i]='d'; // deletion state.
			} else { fprintf(stderr,"Input error 2 ('%c')\n",str[i]); exit(1); }
		   }
		   fprintf(ofp,"(%d)",Length);
		   for(i=1; i <= Length; i++) fprintf(ofp,"*"); 
		   fprintf(ofp,"\n"); 
		} 
		if(num > 1){
		  for(j=0; str[j] && j < Begin; j++) if(isalpha(str[j])) SeqOffSet++;
		}
		for(match=insert=deletion=0, j=0,i=Begin; i <= End; i++){
		   c=str[i];
		   if(isalpha(c)){
			if(MasterSeq[i]=='m'){		// a match in subject
			   Seq[j]=toupper(c); j++;
			   match++;
			} else if(MasterSeq[i]=='d'){  // an insert in subject
			   Seq[j]=tolower(c); j++;
			   insert++;
			} else { fprintf(stderr,"Input error 3\n"); exit(1); }
		   } else if(c == '.' || c == '-'){
			if(MasterSeq[i]=='m'){		// a deletion in subject
			   Seq[j]= '-';	j++;
			   deletion++;
			} else if(MasterSeq[i]=='d'){  // ignore (collapse) these positions 
			   ;
			}
		   } Seq[j]=0;
		}
		fprintf(ofp,"$%d=%d(%d):\n",num,match+insert,match+deletion);
		SetOffSetSeq(SeqOffSet,DummySeq);
		PutSeqInfo(ofp,DummySeq);
		fprintf(ofp,"{()%s()}*\n\n",Seq);
		NilSeq(DummySeq); DummySeq=0;
	}
   } fclose(ifp);
   fprintf(ofp,"\n_0].\n\n"); rewind(ofp);
   cma_typ  out_cma=ReadCMSA(ofp,AB); fclose(ofp);
   free(Seq); free(MasterSeq); free(str); 
   return out_cma;
}

cma_typ OneSeqToCMSA(e_type Seq, char *name, a_type AB)
{
	e_type *EList; NEW(EList,3,e_type); EList[1]=CopySeq(Seq);
	cma_typ tmp=SimSeqToCMSA(EList, 1, AB);
	ReNameCMSA(name,tmp); // free(EList); This is freed by SimSeqToCMSA();
	return tmp;
}

Int4	NumDeletionsCMSA(Int4 blk, Int4 pos, cma_typ cma)
{
	Int4	d,m,sq;
	assert(blk > 0 && blk <= nBlksCMSA(cma));
	assert(pos > 0 && pos <= LengthCMSA(blk,cma));
	for(d=0,m=0,sq=1; sq <= NumSeqsCMSA(cma); sq++){
		if(!SeqIsInCMSA(sq,cma)) continue;
		if(IsDeletedCMSA(blk,sq,pos,cma)) d++;
		else m++;
	} return d;
}

double	FractDeletionsCMSA(Int4 blk, Int4 pos, cma_typ cma)
{ return (double) NumDeletionsCMSA(blk,pos,cma)/(double)NumSeqsCMSA(cma); }

cma_typ ConvertColsToInsertsCMSA(cma_typ cma, Int4 Blk, Int4 start_ins, Int4 end_ins)
// Convert a region within a block into an insert (lower case in cma file).
{
        Int4    blk,*len,numfix;
        assert(Blk > 0 && Blk <= nBlksCMSA(cma));
        assert(LengthCMSA(Blk,cma) > end_ins);
        assert(start_ins <= end_ins);
        assert(start_ins > 1);
	cma_typ	rcma=0;
#if 1	// this is not working for deletions on the ends...in fa2cma.l
        gss_typ *gss=gssCMSA(cma);
        NEW(len,nBlksCMSA(cma) +3, Int4);
        for(blk=1; blk <= nBlksCMSA(cma); blk++) len[blk]=LengthCMSA(blk,cma);
        len[Blk] = len[Blk]-(end_ins - start_ins + 1);
        rcma=EmptyCMSA(nBlksCMSA(cma),len,TrueDataCMSA(cma),gss->GapOpen(),
                gss->GapExtend(),PerNatsCMSA(cma),0,0);
        for(Int4 sq=1; sq <= NumSeqsCMSA(cma); sq++){
                gsq_typ *gsq0,*gsq=gsqCMSA(sq,cma);
                if(!SeqIsInCMSA(sq,cma)){
                    if(gsq==0){
                           gsq0 = new gsq_typ[1];
                           gsq0->initialize(TrueSeqCMSA(sq,cma));    // copies positions to sites array.
                           ReplaceCMSA(sq,gsq0,rcma);
                    } continue;
                }
                // ========== X. Get operational array for sequence. ===========
                // char *operation=gss->Operation(sq);
                Int4    *sites=GetPosSitesCMSA(sq,cma);
#if 0	// DEBUG...
if(sq==43){
gsq->Put(stderr,AlphabetCMSA(cma));
gsq->Put_cma_format(stderr,sq,nBlksCMSA(cma),sites,LengthsCMSA(cma),AlphabetCMSA(cma));
}
#endif
                gsq0=gsq->ConvertColsToInsert(Blk, start_ins, end_ins,
                                        nBlksCMSA(cma),LengthsCMSA(cma),sites);
#if 0	// DEBUG...
if(sq==43){
gsq0->Put(stderr,AlphabetCMSA(cma));
gsq0->Put_cma_format(stderr,sq,nBlksCMSA(cma),sites,LengthsCMSA(cma),AlphabetCMSA(cma));
}
#endif
                ReplaceCMSA(sq,gsq0,rcma); // replace sequence s in CMSA & fmodel.
                for(blk=1; blk <= nBlksCMSA(rcma); blk++) AddSiteCMSA(blk,sq,sites[blk],rcma);
                free(sites);
        } free(len);
        // sprintf(str,"%s.insrt",argv[1]); WriteMtfCMSA(str, cma, NULL);
#else	// Use this instead...
	// for(Int4 s=end_ins; s >= start_ins; s--) ColumnsToInsertCMSA2(cma,s,s);
	ColumnsToInsertCMSA2(cma,start_ins,end_ins);
	FILE *tfp=tmpfile(); PutCMSA(tfp,cma); rewind(tfp); rcma=ReadCMSA(tfp,AlphabetCMSA(cma));
	fclose(tfp); 
#endif
        return rcma;
}

