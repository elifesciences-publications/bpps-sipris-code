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

#include "editcma.h"

ema_typ CMAtoEMA(cma_typ cma)
{
	Int4	start,end,gap0,gap1;
	Int4	N,fsq,sq,sst;
	ema_typ	ema=0;

	gss_typ *gss=gssCMSA(cma);
	ss_type	fulldata = FullSeqCMSA(cma);

	if(!(nBlksCMSA(cma) == 1 && gss && fulldata)) print_error("CMAtoEMA( ) error");
	N = NumSeqsCMSA(cma); 
	for(sq=1; sq<= N; sq++) {
		fsq = SubToFullCMSA(sq,cma);
		sst = RepeatsInfoCMSA(&start,&end,&gap0,&gap1,sq,cma);
		ema = AppendEMA(ema,MakeEMA(fsq,sq,start,end));
	} return ema;
}

cma_typ EMAtoCMA(ema_typ ema, char *psm_arg, cma_typ oldcma)
{
	a_type	A=AlphabetCMSA(oldcma);
	cma_typ	cma=0;
	gss_typ *gss=gssCMSA(oldcma);
	ss_type	fulldata = FullSeqCMSA(oldcma);

	if(!(nBlksCMSA(oldcma) == 1 && gss && fulldata)) print_error("CMAtoEMA( ) error");

// cma=MakeCMSA(ListE,NumHits,operation,start,nblks,LengthsCMSA(cma),
//                 gapo,gapx,5,left,right,NameCMSA(cma),AlphabetCMSA(cma),FullE,FullR);
	e_type	*ListE,fullE;
	ema_typ	tmp;
	Int4	fsq,start,end,s,M=LengthEMA(ema);
	ss_type	truedata = TrueDataCMSA(oldcma);
	NEW(ListE, M+3, e_type);
	unsigned short *FullR;
	NEW(FullR,NSeqsSeqSet(fulldata)+3,unsigned short);
        for(tmp=ema,s=1; tmp; tmp=tmp->next,s++){
		fsq = FullSeqEMA(tmp); FullR[fsq]++;
		fullE = SeqSetE(fsq,fulldata);
		start = StartEMA(tmp);
		start = MAXIMUM(Int4,start - gss->LeftFlank(),1);
		end = EndEMA(tmp);
		end = MINIMUM(Int4,end + gss->RightFlank(),LenSeq(fullE));
// fprintf(stderr,"fsq %d: %d..%d\n",fsq,start,end);
		ListE[s] = MkSubSeq(start,end,fullE);
		if(s > M) print_error("EMAtoCMA( ) runtime error");
	}
        ss_type data=Array2SeqSet(ListE,M,NameCMSA(oldcma),AlphabetCMSA(oldcma));
        cma=EmptyCMSA(nBlksCMSA(oldcma),LengthsCMSA(oldcma),data,gss->GapOpen(),
		gss->GapExtend(), gss->PerNats(),gss->LeftFlank(),gss->RightFlank());
        gss=gssCMSA(cma);
        Int4    *newpos;
        NEW(newpos,nBlksCMSA(cma)+2,Int4);
	psm_typ psm(1,psm_arg,oldcma); 
        for(s=1; s <= M; s++){
          e_type E = gss->TrueSeq(s);
          Int4 oper_len,score;
	  char *operation;
          operation=psm.Align(E,1,&score,&start,&oper_len);
	  oper_len = strlen(operation);
// std::cerr << s; std::cerr << ": "; std::cerr << operation; std::cerr << std::endl;
// fprintf(stderr,"fsq = %d start = %d, end = %d\n",FullSeqEMA(tmp),start,end);
          gsq_typ *gsq; gsq = new gsq_typ[1];
          gsq->initialize(gss->LeftFlank(),gss->RightFlank(),operation,
		oper_len,start,E,newpos);
          ReplaceCMSA(s,gsq,cma);
          for(Int4 t=1; t<=nBlksCMSA(cma); t++) AddSiteCMSA(t,s,newpos[t], cma);
	  free(operation);
        } free(newpos);
	AddFullCountsCMSA(CopySeqSet(fulldata),FullR,cma);
	if(DomainsCMSA(oldcma)) CopyDomainsCMSA(cma,oldcma);
	return cma;
}

cma_typ AddRelatedCMSA(double cluster_cut, double add_cut, char *psm_arg, cma_typ cma)
// Add conserved elements that are pairwise related to other 
// elements in cma.
/*****************************************************************
    --====-----====------====---     --====-----====------====---    
               :  :              -->   
    --====---------------====---     --====-----====------====---
             ('missing')
 *****************************************************************/
// This currently only works for alignments with a single block.
// The input cma_typ cma must contain the full sequences.
// Within a sequence: * --> * --> * --> *  where wt = distance between...?
//                    |           |     |
// Within a sequence: * --------> * --> *  where wt = distance between...?
{ 
	double	CutNeighEval=0.001,prob_cut=0.05;
        double  Ethresh=1,Ecutoff=0.001;
	Int4	s,n,N,fullN;
	ss_type	data,truedata,fulldata;
	a_type	A=AlphabetCMSA(cma);
	Int4	score,s1,s2,q_start,s_start;
	Int4	NumFSets,*FullSet,*FullCard;
	Int4	nhits,qsq,ssq,fsq,sq;
	gss_typ *gss=gssCMSA(cma);
	e_type	subE,fullE,*EList;
	double	evalue;
	Int4	fqsq,fssq;
        Int4    T=11,gap_open=11,gap_extend=1;
        double  x_parameter=15.0;
        UInt4   hpsz=5;
        cma_typ cma2=0;


	Ecutoff=add_cut;
	fulldata = FullSeqCMSA(cma);
	if(fulldata == 0){ print_error("AddRelatedCMSA( ) fulldata == 0 error"); }
	if(!(nBlksCMSA(cma) == 1 && gss && fulldata))
                print_error("AddRelatedCMSA( ) input error");
	fullN = NSeqsSeqSet(fulldata);
	data = DataCMSA(cma); N = NSeqsSeqSet(data); assert(N > 1);
	truedata = TrueDataCMSA(cma); 
	Int4 *FirstRpt=FirstRptFullCMSA(cma);
	unsigned short *FullRpts=FullRptsCMSA(cma);

	ema_typ ema = CMAtoEMA(cma);

	ds_type fullsets = ClusterGPSI('E', cluster_cut, fulldata);
	FullSet = AssignDSets(fullsets, &FullCard, &NumFSets);
	wdg_typ	*WDG;
	Int4	*NFullRpt;
	NEW(WDG, NumFSets+2, wdg_typ);
	NEW(NFullRpt, NumFSets+2, Int4);
	for(fsq=1; fsq <= fullN; fsq++){
	    s = FullSet[fsq];
	    NFullRpt[s] += FullRptsCMSA(fsq, cma);
	}
	for(s=1; s <= NumFSets; s++){ n = NFullRpt[s]; WDG[s]=MkWdgraph(N+2,n*n); }
	st_type	sites=SitesCMSA(cma);
	ds_type	sets = DSets(N);
	h_type HG=Histogram("number of hits for each repeat",0,N,1);
	for(qsq=1; qsq<= N; qsq++) {
		fqsq = SubToFullCMSA(qsq,cma);
		e_type	fqE = SeqSetE(fqsq,fulldata);
		sbp_typ sbp=GBLAST_ScoreBlkNew(A,SeqPtr(fqE)+1,LenSeq(fqE),gap_open,
                	   gap_extend, TotalSeqSet(fulldata), NSeqsSeqSet(fulldata));
		SetStdStatsSBP(sbp);    // use standard statistics
		double x_parameter_final=10.0;
		gab_typ gap_align= GABNew(fqE,gap_open,gap_extend,SMatrixSBP(sbp),
			  -(GBLAST_SCORE_MIN),GapXDropoffSBP(x_parameter_final,sbp),0); 

		s1 = findDSets(qsq,sets);
		e_type qE = SeqSetE(qsq,data);
		double	prob;
		smx_typ smx = SMXforSeqCMSA(1,qsq,cma);
		for(nhits=0,ssq=1; ssq <= N; ssq++) {
		   fssq = SubToFullCMSA(ssq,cma);
		   if(ssq == qsq || FullSet[fqsq] != FullSet[fssq]) continue;
		   // otherwise sequences beInt4 in the same fullset...
		   wdg_typ wdg = WDG[FullSet[fqsq]];
		   e_type sE = SeqSetE(ssq,data);
		   Int4 sts = SitePos(1,ssq,1,sites);
		   score = ScoreSMatrix(SeqPtr(sE),sts, smx);
		   if(score > 0){
			prob = SMatrixProb(score, smx);
			if(prob <= prob_cut){
				Int4 word_score,ssq_start;
				e_type ssqE = MaxWordAndSubSeqCMSA(1,qsq,ssq,
					&word_score,&q_start,&s_start,&ssq_start,cma);
				e_type fsE = SeqSetE(fssq,fulldata);
				SetSubjectGAB(q_start,s_start-ssq_start+1,ssqE, gap_align);
				PerformGappedAlignment(gap_align);
				evalue= GappedScoreToEvalueSBP(gap_align->score,sbp);
				if(evalue <= CutNeighEval){
		   		   s2 = findDSets(ssq,sets);
		   		   if(s1 != s2){ s1 = linkDSets(s1,s2,sets); }
	   			   JoinWdgraph(qsq,ssq,(Int4)(-100.*log10(evalue)),wdg);
				   nhits++;
				} NilSeq(ssqE);
			}
		   }
		} IncdHist(nhits,HG);
		NilSMatrix(smx); GBLAST_ScoreBlkDestruct(sbp); GABDelete(gap_align);
	} // PutDSets(stdout,sets);
	PutHist(stderr,60,HG); NilHist(HG);
	Int4 NumSets=0,*Card;
	Int4 *SubSet=AssignDSets(sets, &Card, &NumSets);
	// store information about matches in related sequences...
	BooLean	**FullMatch;	// FullMatch[sq][fsq] == full match in fsq to sq.
	NEWP(FullMatch,N+3, BooLean);
	for(qsq=1; qsq<= N; qsq++) NEW(FullMatch[qsq],fullN+3, BooLean);
	// find full seqeunces matches for each rpt...
	for(qsq=1; qsq<= N; qsq++) {
	    fqsq = SubToFullCMSA(qsq,cma);
	    s = FullSet[fqsq];
	    for(ssq=1; ssq <= N; ssq++){
		fssq = SubToFullCMSA(ssq,cma);
		if(fqsq != fssq && s == FullSet[fssq]){
		    if(SubSet[qsq] == SubSet[ssq]){
			 // i.e., both in same fullset and Subset...
			FullMatch[qsq][fssq]=FullMatch[ssq][fqsq]=TRUE;
		    }
		}
	    }
	}
	HG=Histogram("SW alignment scores",0,1000,5.0);
	for(qsq=1; qsq<= N; qsq++) {
	    subE = SeqSetE(qsq,truedata);
	    fqsq = SubToFullCMSA(qsq,cma); s = FullSet[fqsq];
	    Int4 qst,qstrt,qend,qap0,qap1,gap0,gap1,r,qos,qes;
	    qst = RepeatsInfoCMSA(&qstrt,&qend,&qap0,&qap1,qsq,cma);

	    for(fsq = 1; fsq <= fullN; fsq++){
	      // if full sequence related but no match...
	      if(fqsq != fsq && s == FullSet[fsq] && !FullMatch[qsq][fsq]){
		fullE = SeqSetE(fsq,fulldata);
		NEW(EList,4,e_type); 
		EList[1] = CopySeq(fullE); EqSeqI(1,EList[1]);
		ss_type tmp_data = Array2SeqSet(EList, 1, "tmp_fafile",A);

		FILE *fp=0;
		gpsi_type gpsi(subE,tmp_data,Ethresh,Ecutoff,x_parameter,hpsz,1,fp);
		gpsi.KeepQuite( ); // don't want output...
		sap_typ sap = gpsi.Search(T,' ');
		if(sap){	// found a potential match...
		    // Assume NextGSAP(sap) == 0 for now...
		    // find region of match...
		    score = ScoreGSAP(sap); evalue = EvalueGSAP(sap);
		    Int4 QS,QE,SS,SE,SSA,SEA,len_ssq,before,after;
		    len_ssq = InfoGSAP(&QS,&QE,&SS,&SE,sap);
// look through other repeats info for additional hits.
		    IncdHist(score,HG);
		    // now look for a segment in this region...
		    Int4 sst,strt,end,flank=3;
		    qos = qstrt - QS; qes = qend - QE;
		    SSA = SS + qos; SEA = SE + qes;
	   	    BooLean flag=TRUE;
		    before=after=0;
		    for(ssq = FirstRpt[fsq], r=1; r <= FullRpts[fsq]; r++, ssq++){
			sst = RepeatsInfoCMSA(&strt,&end,&gap0,&gap1,ssq,cma);
			// does alignment significantly overlap with any repeats?
		        if((SSA+flank) <= end && (SEA-flank) >= strt){
			  if(FALSE) fprintf(stdout,"  %d (%d): (%d) %d-%d (%d) (%d-%d)\n",
				ssq,r,gap0,strt,end,gap1,SSA,SEA);
			  flag=FALSE;
			} else if(end < (SSA+flank)) {   // then repeat is before
			   before = ssq;
			} else if(strt > (SEA-flank)) {  // & (SSA+flank)  > end
			   if(after==0) after = ssq;
			} else print_error("This should not happen!!");
		    }
		    if(flag){	// Found a new motif... Create a new sequence...
			fprintf(stdout,"  ********** Potential New Site!!! **********\n");
		        // gpsi.ShowAlign(stdout,60,1);
		        // gpsi.PutSubSeqs(stderr,5,5);  // left & right flank = 5.
			PutSeqID(stdout,subE);
		        fprintf(stdout," (%d): %d-%d --> %d-%d (score: %d)(E=%.2g) hits ",
				RptNumCMSA(qsq,cma),QS,QE,SS,SE,score,evalue); 
				// RptNumCMSA(qsq,cma),fqsq,QS,QE,SS,SE,score,evalue); 
			PutSeqID(stdout,fullE); fprintf(stdout,"\n");
			fprintf(stdout," motif: (%d) %d-%d (%d) --> %d-%d between %d(%d) & %d(%d)\n",
				qap0,qstrt,qend,qap1,SSA,SEA,
				RptNumCMSA(before,cma),before,RptNumCMSA(after,cma),after);

			ema = AppendEMA(ema,MakeEMA(fsq,before,SSA,SEA,score,evalue));
		    }
		} NilSeqSet(tmp_data);
	      }
	    }
	}
#if 0
	cma2=MakeCMSA(ListE,NumHits,Operation,Start,nblks,lengths,
                        gapo,gapx,5,left,right,name,A,FullE,FullR);
#endif
	PutHist(stderr,60,HG); NilHist(HG);
	for(s = 1; s <= NumFSets; s++){ 
	   printf("set %3d:\n",s);
	   for(fsq = 1; fsq <= fullN; fsq++){
		if(s == FullSet[fsq]){
		   fullE = SeqSetE(fsq,fulldata);
		   printf("%3d: ",SeqI(fullE));
		   PutSeqInfo(stdout,fullE);
		   printf("{ ");
		   for(sq =1; sq <= N; sq++){
			if(fsq == SubToFullCMSA(sq,cma)){ printf("%d ",sq); }
		   } printf("}\n");
		}
	   }
	   printf("\n");
	   PutWdgraph(stdout, WDG[s]); 
	}
	for(qsq=1; qsq<= N; qsq++) free(FullMatch[qsq]); free(FullMatch);
	for(s = 1; s <= NumFSets; s++) NilWdgraph(WDG[s]); 
	free(WDG); NilDSets(sets); NilDSets(fullsets); free(NFullRpt);
	free(FullSet); free(SubSet); free(Card);  free(FullCard);
	PutEMA(stderr,ema);
	cma2 = EMAtoCMA(ema,psm_arg,cma); NilEMA(ema);
        return cma2;
}

// TEST NEW PROCEDURES: 
char	SampGSeqTempCMSA(Int4 s, BooLean *skipseq, double temperature, 
	UInt4 min_rpt, cma_typ *cma)
// sample a gapped alignment for sequence s using alignment cma model.
// need a way to compute map for this to see if should take gaps??
{ double OldMap;
    return SampleGSeqGibbsCMSA(s,skipseq,cma,&OldMap,temperature,1.0,min_rpt,0); }

char	SampGSeqTempCMSA(Int4 s, BooLean *skipseq, double temperature, 
	double *OldMap, UInt4 min_rpt, cma_typ *cma)
// sample a gapped alignment for sequence s using alignment cma model.
// need a way to compute map for this to see if should take gaps??
{ return SampleGSeqGibbsCMSA(s,skipseq,cma,OldMap,temperature,1.0,min_rpt,0); }

char	SampGSeqTempCMSA(Int4 s, BooLean *skipseq, double temperature, 
	double prob_accept, double *OldMap, UInt4 min_rpt,
	cma_typ *cma,psm_typ *psm)
// sample a gapped alignment for sequence s using alignment cma model.
// need a way to compute map for this to see if should take gaps??
{ return SampleGSeqGibbsCMSA(s,skipseq,cma,OldMap,temperature,1.0,min_rpt,psm); }

char	SampGSeqTempCMSA(Int4 s, BooLean *skipseq, double temperature, 
	double prob_accept, double *OldMap, UInt4 min_rpt, hma_typ *hma)
// sample a gapped alignment for sequence s using alignment cma model.
// need a way to compute map for this to see if should take gaps??
{ return hma->SampleGSeqGibbs(s,skipseq,OldMap,temperature,1.0,min_rpt); }

char	comp_map_cmsa3(double rmap, double nmap, double omap, double temperature)
/** short routine to be used for sampling three different maps **/
{
	double	ratio,maxmap,minmap,midmap;
	char	max,min,mid;

	if(rmap > nmap){
	   if(rmap > omap){ 	// i.e., rmap > nmap & rmap > omap.
	     maxmap=rmap; max = 'r'; 
	     if(omap > nmap){			// r > o > n.
		midmap = omap; minmap=nmap; 
		mid = 'o'; min = 'n'; 
	     } else { 				// r > n > o.
		midmap = nmap; minmap=omap; 
		mid = 'n'; min = 'o'; 
	     }
	   } else {				// o > r > n.
		maxmap=omap; midmap=rmap; minmap=nmap; 
		max = 'o'; mid = 'r'; min = 'n';
	   }
	} else {		
	   if(nmap > omap){	// i.e., nmap > rmap & nmap > omap.
		maxmap=nmap; max = 'n'; 
		if(omap > rmap){		// n > o > r.
		   midmap = omap; minmap=rmap; 
		   mid = 'o'; min = 'r'; 
		} else {			// n > r > o.
		   midmap = rmap; minmap=omap;
		   mid='r'; min = 'o'; 
		}
	   } else {				// o > n > r.
		maxmap=omap; midmap= nmap; minmap=rmap; 
		max = 'o'; mid = 'n'; min = 'r';
	   }
	}
	if(temperature <= 50.0) return max;
	else temperature = 300./temperature;
#if 0
 fprintf(stderr,"rmap = %g; nmap = %g; omap = %g\n",rmap,nmap,omap);
 fprintf(stderr,"maxmap = %g; midmap = %g; minmap = %g\n", maxmap,midmap,minmap);
 fprintf(stderr,"max = %c; mid = %c; min = %c\n",max,mid,min);
#endif
	if((maxmap - midmap) > 100) return max;
	else { maxmap -= minmap; midmap -= minmap; }
	if(maxmap > 100) { 		// i.e., discard minmap.
	   maxmap -= midmap;
	   if(maxmap > 100) return max;
	   else {
		ratio = exp(maxmap);	/** ratio of new to old **/
		if(temperature != 1.0) ratio = pow(ratio,temperature);
		ratio = ((double)Random()/(double)RANDOM_MAX)*(ratio+1);
		if(ratio >= 1.) return max; else return mid;
	   }
	} else {
	   double total,rand;
#if 0
fprintf(stderr,"rel: maxmap = %g; midmap = %g; minmap = %g\n", maxmap,midmap,minmap);
	   maxmap = exp(maxmap); midmap = exp(midmap); minmap = 1.0;
fprintf(stderr,"exp: maxmap = %g; midmap = %g; minmap = %g\n", maxmap,midmap,minmap);
#endif
	   if(temperature != 1.0){
		 maxmap = pow(maxmap,temperature);
		 midmap = pow(midmap,temperature);
	   }
	   total = maxmap + midmap + 1.0;
	   rand = total*((double)Random()/(double)RANDOM_MAX);
	   // total = |------max-+----|---mid--|-|
	   // rand  = |----------^
#if 0
fprintf(stderr,"total = %g; rand = %g; maxmap = %g; midmap = %g; minmap = %g\n",
		total,rand,maxmap,midmap,minmap);
 fprintf(stderr,"max = %c; mid = %c; min = %c\n",max,mid,min);
#endif
	   if((rand -= maxmap) <= 0.0) return max; 
	   else if((rand - midmap) <= 0.0) return mid;
	   else return min;
	}
}

double  GetLPRGappedSqCMSA(Int4 s, cma_typ cma, psm_typ *psm)
// get the probability ratio of sequence s being in the alignment
// WARNING: assumes gap penalties are meaningful!!!
{
        Int4    t,*oldpos;
        double  oldmap,newmap;

	if(psm) oldmap=psm->RelMap(cma); 
        else oldmap=RelMapCMSA(cma);
        // 1. Remove sequence s from alignment.
        NEW(oldpos,nBlksCMSA(cma)+2,Int4);  // save old sites.
        for(t=1; t<=nBlksCMSA(cma); t++){
                PosSiteCMSA(t, s, cma->pos, cma); oldpos[t]=cma->pos[1];
        } VacateSitesCMSA(s,cma);
	if(psm) newmap=psm->RelMap(cma); 
        else newmap=RelMapCMSA(cma);
        for(t=1 ; t <= nBlksCMSA(cma); t++) AddSiteCMSA(t,s,oldpos[t], cma);
        free(oldpos);
        return ((oldmap-newmap)*0.43429448); // convert to log10
}

char	SampleGSeqGibbsCMSA(Int4 s, BooLean *skipseq, cma_typ *oldcma,
	double *OldMap, double Temperature, double prob_accept,
	UInt4 min_rpt, psm_typ *psm)
// For sampling gapped sequences into and out of cmsa.
{
	Int4	m,n,t,start,trace_length,gapopen,gapextend;
	Int4    score;
	Int4	*newpos;
	smx_typ	*smx;
	a_type	A=AlphabetCMSA(*oldcma);
	e_type	oldE,E;
	gss_typ	*gss=gssCMSA(*oldcma);
	double	pernats;
	char	*operation;
	Int4	wt;
	cma_typ	cma;
	double	oldmap,newmap,rm_map;
	double	inc=10.0;
	BooLean	CanDelete=FALSE;

	// assert(prob_accept > 0.0);
	cma=CopyCMSA(*oldcma);
	if(psm) oldmap=psm->RelMap(*oldcma); 
	else oldmap=RelMapCMSA(*oldcma);
	*OldMap = oldmap;
	gapopen=gss->GapOpen(); gapextend=gss->GapExtend(); pernats=gss->PerNats();

	if(Temperature >=300.0) wt=1; // i.e., Natural temperature...
	else if(Temperature > 75.0){
	   wt = 1 + (Int4) floor((300.0+inc-Temperature)/inc);
	   // for inc = 10.0 --> wt for T=290-300 is 2;  280-290 is 3; 270-280 is 4; etc...
	   // down to 75-80 = 24.
	}  else wt=0;  // Use maximum likelihood.

	if((cma)->FullSeq){
	    if((cma)->FullRpts[(cma)->SubToFull[s]] > min_rpt) CanDelete=TRUE;
	    else CanDelete=FALSE;
	} else CanDelete=TRUE;
	// 1. Remove sequence s from alignment.
	NEW(newpos,nBlksCMSA(cma)+2,Int4);  // for new sites.
	if(!skipseq[s] && CanDelete){
	  VacateSitesCMSA(s,cma); 
	  if((cma)->FullSeq){ (cma)->FullRpts[(cma)->SubToFull[s]]--; }
	  if(psm) rm_map=psm->RelMap(s,cma); 
	  else rm_map = RelMapCMSA(cma) + gss->IndelPenalty(s);
#if 1	// NEW prob_accept
	  rm_map = rm_map - prob_accept; // scale up rm_map to favor removal by prior prob.
#endif
	  if((cma)->FullSeq){ (cma)->FullRpts[(cma)->SubToFull[s]]++; }
// if(oldmap < rm_map){ fprintf(stderr,"#### oldmap = %g; rm_map = %g\n",oldmap,rm_map); }
	} else {
		// rm_map = -200.0; 
		rm_map = (double)LONG_MIN; 
	}

        E = gss->TrueSeq(s); 	
        oldE = gss->FakeSeq(s); 	
	// PutSeq(stderr,E,AlphabetCMSA(cma)); fflush(stderr);

	// ===== 2. Obtain scoring matrix from rest of alignment.

	// ========== 4. Sample a gapped alignment for sequence s. ===========
	if(psm){
	  // assert(nBlksCMSA(cma) == 1);
	  if(wt == 0) operation=psm->Align(E,1,&score,&start,&trace_length);
	  else operation=psm->SampleAlign(E,1,&score,&start,&trace_length,wt);
	} else {
	 NEW(smx,nBlksCMSA(cma) + 2,smx_typ);
	 for(n=0,m=1; m<= nBlksCMSA(cma); m++){
	   smx[m]=SampleDirichletSmatrixFModel(pernats,wt,ModelCMSA(m,cma));
           // if(wt == 0) smx[m]=GetSmatrixFModel(pernats,ModelCMSA(m,cma));
	   // else smx[m]=SampleSmatrixFModel(pernats,wt,ModelCMSA(m,cma));
	   n+=LenSMatrix(smx[m]);
	   // PutSMatrix(stderr,smx[m]);
	 }
	 operation=GapAlnTraceSMatrix(gapopen,gapextend,LenSeq(E),
			SeqPtr(E),nBlksCMSA(cma),smx, NULL,&start);
	 trace_length=strlen(operation);
	 for(m=1; m<= nBlksCMSA(cma); m++) NilSMatrix(smx[m]); free(smx);
	}
	// ========== 5. Create a gapped sequence. ===========
	gsq_typ	*gsq; gsq = new gsq_typ[1];
        gsq->initialize(gss->LeftFlank(),gss->RightFlank(),
				operation,trace_length,start,E,newpos);
	// ========== 6. Deallocate memory. ===========
        free(operation); 

	// ========== 7. If sequence changed then replace s with gsq. =========
	if(!skipseq[s] && gss->Identical(s,*gsq)){ 
	  free(newpos); delete []gsq; 
	  // newmap = -200.0; 
	  newmap = (double)LONG_MIN; 
	} else {
	  ReplaceCMSA(s,gsq,cma); // replace sequence s in CMSA & fmodel.
	  for(t=1; t<=nBlksCMSA(cma); t++) AddSiteCMSA(t,s,newpos[t], cma);
	  free(newpos);
	  if(psm) newmap=psm->RelMap(cma); 
	  else newmap= RelMapCMSA(cma);
	}
	// Now Sample one of these alignments.
	char option;
// fprintf(stderr,"#### oldmap = %g; newmap = %g; rm_map = %g\n",oldmap,newmap,rm_map); 
// oldmap=0; *OldMap = oldmap;
	if(Temperature<=10.0) option = comp_map_cmsa3(rm_map,newmap,oldmap,0.0);
	else option = comp_map_cmsa3(rm_map,newmap,oldmap,Temperature);
	switch(option){
	  case 'r': {   // remove sequence from alignment.
		skipseq[s]=TRUE;   
		ReplaceCMSA(s,NULL,*oldcma);
		if((*oldcma)->FullSeq){
		  (*oldcma)->FullRpts[(*oldcma)->SubToFull[s]]--;
		}
		gss_typ     *gss0=gssCMSA(cma);
		gss0->~gss_typ(); NilCMSA(cma); 
	    } break;
	  case 'o': {  // retain old sequence.
		gss_typ     *gss0=gssCMSA(cma);
		gss0->~gss_typ(); NilCMSA(cma); 
	    } break;
	  case 'n': {  // convert to new sequence;
		if(skipseq[s]) {
		    skipseq[s]=FALSE;   
		    option = 'a';
		    if(cma->FullSeq){
		      cma->FullRpts[cma->SubToFull[s]]++;
		    }
		}
		gss_typ     *gss0=gssCMSA(*oldcma);
		gss0->~gss_typ(); NilCMSA(*oldcma); 
		*oldcma=cma; 
	    } break;
	  default:  print_error("This should not happen"); break;
	} return option;
}


