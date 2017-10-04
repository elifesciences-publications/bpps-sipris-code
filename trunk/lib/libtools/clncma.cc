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

#include "clncma.h"

#define USAGE_CLEANCMA     "USAGE: cleancma cmafafile [options]\n\
   options:\n\
     -e<real>      - set e-value for pairwise comparisons to <real> (default: 0.05)\n\
     -C<real>      - gapped log10 likelihood cutoff for inclusion in alignment (3.0)\n\
     -c<real>      - adjusted gapped log10 likelihood cutoff for inclusion (3.0)\n\
     -L            - look at histogram of number of hits for each subsequence\n\
     -L<seqid>     - look at subsequences of <seqid> from the alignment\n\
     -l<real>      - look cutoff gapped log10(likelihood)(default: 3.0)\n\
     -N<real>      - Neighbor E-value cutoff to inherit status (default: 0.001)\n\
     -p<int>       - xpurge MA sequences at cutoff <int>\n\
     -S            - shuffle (or reverse) sequences to perform Monte Carlo\n\
     -v            - verbose mode\n\
     -z            - dummy\n\n"

int	CleanCMSA(Int4 argc,char *argv[], a_type A)
{ 
	Int4	arg,i,cutoff=-999;
	Int4    N,fullN,left=0,right=0,TrimMax=3;
	char	str[300],mode='c',*rm=NULL,*look=0;
	double	look_cut=2.0,prob_cut=0.05;
	double	Cut_gP=3.0,Cut_agP=3.0,CutNeighEval=0.001;
	cma_typ	cma=0;
	BooLean	shuffle=FALSE,verbose=FALSE;
	Int4	*value;
	ss_type	data;
	FILE	*fp;
	Int4	min_rpt=0,min_spacing=0;

	if(argc < 2) print_error(USAGE_CLEANCMA);
	for(arg = 2; arg < argc; arg++){
	   if(argv[arg][0] != '-') print_error(USAGE_CLEANCMA);
	   switch(argv[arg][1]) {
	     case 'e': prob_cut = RealOption(argv[arg],'e',0.0,10000.0,USAGE_CLEANCMA);
                        break;
	     case 'L': if(argv[arg][2]!=0) look=argv[arg]+2; mode = 'L'; break;
	     case 'l': look_cut=RealOption(argv[arg],'l',-1000.0,1000.0,USAGE_CLEANCMA); 
			break;
	     case 'N': CutNeighEval=RealOption(argv[arg],'N',0.0,1000.0,USAGE_CLEANCMA); 
			break;
	     case 'C': Cut_gP=RealOption(argv[arg],'C',-1000.0,1000.0,USAGE_CLEANCMA); 
			break;
	     case 'c': Cut_agP=RealOption(argv[arg],'c',-1000.0,1000.0,USAGE_CLEANCMA); 
			break;
	     case 'p': cutoff=IntOption(argv[arg],'p',2,5000,USAGE_CLEANCMA); 
			break;
             case 'S': shuffle = TRUE; break;
	     case 'v': verbose = TRUE; break;
	     case 'z': break;
	     default: print_error(USAGE_CLEANCMA);
	   }
	}
	sprintf(str,"%s.cma",argv[1]);
	cma=ReadCMSA2(str,A);
	gss_typ     *gss=gssCMSA(cma);
	ss_type	fulldata = FullSeqCMSA(cma);

	assert(nBlksCMSA(cma) == 1);
	assert(gss);
	assert(fulldata);

	data = DataCMSA(cma);
	N = NSeqsSeqSet(data);
	fullN = NSeqsSeqSet(fulldata);
	assert(N > 1);
	if(look==0){  // then output histogram of all hits
	   Int4		qst,start,score,num_set;
	   double	*Prob,*gProb;
	   st_type	sites=SitesCMSA(cma);
	   ds_type	sets = DSets(N);
	   Int4		s1,s2,end,*gap0,*gap1;
	   Int4		q_start,s_start;
	   char		Str[108];
	   BooLean	*skip;

// calculate Gap function....
#if 1	// LOOK AT FULLDATA gaps
	Int4	maxLength=MaxSeqSeqSet(fulldata);
	Int4	*Gaps;
	NEW(Gaps, maxLength+3,Int4);
#endif
	   // PutDSets(stdout,sets);
	   h_type HG=Histogram("number of hits for each repeat",0,N,1);
	   Int4 nhits,lenM= LengthCMSA(1,cma);
	   Int4 *Start,*End;
	   NEW(Prob,N+3,double); NEW(gProb,N+3,double); 
	   NEW(gap0,N+3,Int4); NEW(gap1,N+3,Int4); 
	   NEW(Start,N+3,Int4); NEW(End,N+3,Int4); 
	   Int4	qsq,ssq,fsq=0,lastfsq=0;
	   Int4 *netlen,*num_rpts;
	   NEW(netlen, fullN+3,Int4);
	   NEW(num_rpts, fullN+3,Int4);
	   Int4 *lenRpt; NEW(lenRpt, N+3,Int4);
	   Int4 g0,g1;
	   for(i=1; i<= N; i++) {
		gProb[i] = GetGappedProbCMSA(1,i,cma);
		Prob[i] = GetProbCMSA(1,i,cma);
		qst = RepeatsInfoCMSA(&start,&end,&g0,&g1,i,cma);
		lenRpt[i] = end - start + 1;
		gap0[i] = g0; gap1[i] = g1;
		if(g0 < 0) g0=0; if(g1 < 0) g1=0;
		Start[i] = start; End[i] = end;
		fsq=SubToFullCMSA(i,cma);
		if(lastfsq != fsq){ Gaps[g0]++;	lastfsq=fsq; netlen[fsq] += g0; }
		Gaps[g1]++; netlen[fsq] += g1; num_rpts[fsq]++;
		
	   }
	   if(verbose) PutHistX(stderr,"Observed gaps",2.0,maxLength, Gaps);
	   Int4 TotNetLen=0;
	   for(i=1; i<= fullN; i++) TotNetLen+=netlen[i];
	   double rpts_per_res = (double)N/(double)TotNetLen;
	   // double *dist=SmoothGapsX(12, 12, 4, maxLength, Gaps, rpts_per_res);
	   double *dist=SmoothGapsX(9, 9, 4, maxLength, Gaps, rpts_per_res);
	   // There is a problem with this PutHistX( ) routine...(when compiling for 64 bit).
	   // if(verbose) PutHistX(stderr, "Smoothed observed gaps", 2.0, maxLength, dist);
	   double *Like;
	   NEW(Like,N+3,double);
	   for(i=1; i<= N; i++) {
		fsq=SubToFullCMSA(i,cma);
	        double  *exp_dist = ExpectedGaps(netlen[fsq],num_rpts[fsq]);
	        double  *expNot_dist = ExpectedGaps(netlen[fsq]+lenRpt[i],num_rpts[fsq]-1);
		g0=gap0[i]; g1=gap1[i]; 
		if(g0 < 0) g0 = 0; if(g1 < 0) g1 = 0;
		double like0 = log10(dist[g0]/exp_dist[g0]);
		double like1 = log10(dist[g1]/exp_dist[g1]);
		Int4 z = g0+g1+lenRpt[i];
		double	likeNot = log10(dist[z]/expNot_dist[z]);
		// double	likeNot1 = log10(expNot_dist[z]/(exp_dist[g0]+exp_dist[g1]));
		double	likeNot2 = log10((dist[g0]*dist[g1])/(exp_dist[g0]*exp_dist[g1]));
		likeNot2 = likeNot2 - likeNot;
#if 0
	fprintf(stderr,"%d(%d): like0(%d) = %g; like1(%d) = %g; likeNot(%d) = %g; Like=%g\n",
			netlen[fsq],num_rpts[fsq],g0,like0,g1,like1,z,likeNot,likeNot2);
#endif
		if(verbose) fprintf(stderr,
			"%d rpt(%d)[%d..%d] = %g.. %g (net = %d; rpts = %d)\n",
			fsq,i,gap0[i],gap1[i],like0,like1,netlen[fsq],num_rpts[fsq]);
		Like[i] = (like0 + like1)/2.0;
		free(exp_dist);
		free(expNot_dist);
	   } free(dist);
// End calculate Gap function....
	   NEW(skip,N+3,BooLean);
	   for(i=1; i<= N; i++) skip[i]=TRUE;
	   for(i=1; i<= N; i++) {
		if(gProb[i] >= Cut_gP) skip[i] = FALSE;
		double agP = gProb[i] + Like[i];
		if(agP >= Cut_agP) skip[i] = FALSE;
	   }
// if(TRUE) exit(1);
	   Int4 Edges,nEdges=0;
	   if(N > 5000) print_error(USAGE_CLEANCMA);  // don't let edges get out of hand.
	   if(N <= 1000) Edges = N*N;
	   else Edges = 1000000;  // allocates > 20 bytes per edge ==> about 20 MB max.
	   wdg_typ	wdg=MkWdgraph(N+2,N*N); // no edges to itself.
	   for(qsq=1; qsq<= N; qsq++) {
		// gapxdrop routines...
		// qst = RepeatsInfoCMSA(&start,&end,&gap0,&gap1,qsq,cma);
		double	evalue;
		Int4	gap_open=10,gap_extend=1;
		e_type	fqE = SeqSetE(SubToFullCMSA(qsq,cma),fulldata);
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
		   if(ssq == qsq) continue;
		   e_type sE = SeqSetE(ssq,data);
		   Int4 sts = SitePos(1,ssq,1,sites);
		   score = ScoreSMatrix(SeqPtr(sE),sts, smx);
		   if(score > 0){
			prob = SMatrixProb(score, smx);
			if(prob <= prob_cut){
				//  gapxdrop routines...
				Int4 word_score,ssq_start;
				e_type ssqE = MaxWordAndSubSeqCMSA(1,qsq,ssq,
					&word_score,&q_start,&s_start,&ssq_start,cma);
				e_type fsE = SeqSetE(SubToFullCMSA(ssq,cma),fulldata);
				SetSubjectGAB(q_start,s_start-ssq_start+1,ssqE, gap_align);
				PerformGappedAlignment(gap_align);
				evalue= GappedScoreToEvalueSBP(gap_align->score,sbp);
				// if(evalue <= prob_cut){
				if(evalue <= CutNeighEval){
		   		   s2 = findDSets(ssq,sets);
		   		   if(s1 != s2){ s1 = linkDSets(s1,s2,sets); }
				   nEdges++; 
				   if(nEdges > Edges) print_error("out of edges");
	   			   JoinWdgraph(qsq,ssq,(Int4)(-100.*log10(evalue)),wdg);
				   if(evalue <= CutNeighEval){
				     if(skip[qsq] != skip[ssq]){
					if(skip[qsq]) skip[qsq] = FALSE; 
					else skip[ssq] = FALSE; 
				     }
				   }
				   nhits++;
				}
				NilSeq(ssqE);
			}
		   }
		} IncdHist(nhits,HG);
		NilSMatrix(smx); GBLAST_ScoreBlkDestruct(sbp); GABDelete(gap_align);
	   } // PutDSets(stdout,sets);
	   PutHist(stdout,60,HG); NilHist(HG);
	   num_set=0;
	   char	gap1str[10],gap0str[10];
	   BooLean flag;
	   for(qsq=1; qsq<= N; qsq++) {
	     if((s1=findDSets(qsq,sets)) == qsq){
	       	num_set++;
		sprintf(Str,"set(%d)",num_set);
		flag=FALSE;
		if(skip[qsq]){
	          for(ssq=1; ssq <= N; ssq++){
		    s2 = findDSets(ssq,sets);
	            if(s2 == s1){
			if(!skip[ssq]) flag=TRUE;
		    }
		  } if(flag) skip[qsq]=FALSE;
		} else flag=TRUE;
		for(ssq=1; ssq <= N; ssq++){
		  s2 = findDSets(ssq,sets);
	          if(s2 == s1){
		     if(flag && skip[ssq]) skip[ssq]=FALSE;
		     StrSeqID(str,15,SeqSetE(ssq,data));
		     sprintf(gap0str,"(%d)",gap0[ssq]);
		     sprintf(gap1str,"(%d)",gap1[ssq]);
		     if(verbose) fprintf(stderr,
"%-9s: %4d %-15s %5s %4d-%-4d %-5s gp: %.2f p: %.2f g_Like: %.2f => %.2f\n",
			   Str,ssq,str,gap0str,Start[ssq],End[ssq],gap1str,
			   gProb[ssq],Prob[ssq],Like[ssq],gProb[ssq]+Like[ssq]);
		     sprintf(Str,"      ");
	          }
	        }
	     }
	   }
	   NilDSets(sets); free(Prob); free(gProb); 
	   free(gap0); free(gap1); free(Start); free(End); 
	   // PutWdgraph(stdout, wdg);
	   NilWdgraph(wdg);
	   fp = open_file(argv[1],".cln.cma","w");
           PutSelectCMSA(fp,skip,cma);
           fclose(fp);
	   free(skip); free(Like);  free(lenRpt);
	   free(Gaps); free(netlen); free(num_rpts);
	} else {
		char str2[108];
		Int4 v,n2,end,end0,qst,j,start,score;
		Int4	total_internal=0;
		Int4 gap0,gap1;
		double	g_like,ug_like;
		st_type	sites=SitesCMSA(cma);
		NEW(value,N+3,Int4);
		// prob_cut /= (double)N;
		Int4	sum,nhits,lenM= LengthCMSA(1,cma);
		Int4	q_start,s_start;
		Int4	diag=0;
		hsp_typ hsp;
	        for(end=v=n2=0,i=1; i<= N; i++) {
		   e_type qE = SeqSetE(i,data);
		   StrSeqID(str2,100,qE);
		   if(strncmp(str2,look,100) == 0){
			if(n2 == 0) {
			  if(verbose) fprintf(stderr,"%s = seq %d\n",look,i);
			  end0 = 0; 
			} else { end0 = end; }
			n2++;
			g_like = GetGappedProbCMSA(1,i,cma);
			ug_like = GetProbCMSA(1,i,cma);
			qst = RepeatsInfoCMSA(&start,&end,&gap0,&gap1,i,cma);
			if(n2 == 1 && gap0 > 0) total_internal+=gap0;
			if(gap1 > 0) total_internal+=gap1;
			smx_typ smx = SMXforSeqCMSA(1,i,cma);

			// gapxdrop routines...
			double	evalue;
			Int4	gap_open=10,gap_extend=1;
			e_type fqE = SeqSetE(SubToFullCMSA(i,cma),fulldata);
			sbp_typ sbp=GBLAST_ScoreBlkNew(A,SeqPtr(fqE)+1,LenSeq(fqE),gap_open,
                		gap_extend, TotalSeqSet(fulldata), NSeqsSeqSet(fulldata));
			SetStdStatsSBP(sbp);    // use standard statistics
			double x_parameter_final=10.0;
			gab_typ gap_align= GABNew(fqE,gap_open,gap_extend,SMatrixSBP(sbp),
			     -(GBLAST_SCORE_MIN),GapXDropoffSBP(x_parameter_final,sbp),0); 

// h_type HG=Histogram("number of hits for each repeat",-500,500,5.0);
			double j_like,max_like=-99999.0,max_evalue = 999999999.0;
			for(sum=nhits=0,j=1; j <= N; j++) {
			   if(j == i) continue;
			   e_type sE = SeqSetE(j,data);
			   if(verbose) StrSeqID(str2,100,sE);
			   Int4 sst = SitePos(1,j,1,sites);
			   score = ScoreSMatrix(SeqPtr(sE),sst, smx);
// IncdHist(score,HG);
			   if(score > 0){
			     double	prob = SMatrixProb(score, smx);
			     if(prob <= prob_cut){
				//  gapxdrop routines...
				Int4 word_score,ssq_start;
				e_type ssqE = MaxWordAndSubSeqCMSA(1,i,j,
					&word_score,&q_start,&s_start,&ssq_start,cma);
				e_type fsE = SeqSetE(SubToFullCMSA(j,cma),fulldata);
				SetSubjectGAB(q_start,s_start-ssq_start+1,ssqE, gap_align);
				PerformGappedAlignment(gap_align);
				evalue= GappedScoreToEvalueSBP(gap_align->score,sbp);
				if(evalue <= prob_cut){
			  	   j_like = GetGappedProbCMSA(1,j,cma);
				   if(max_evalue > evalue){
					 max_evalue = evalue; max_like = j_like;
				   }
			           sum+=gap_align->score;
			           //sum+=score;
				   nhits++;
	// if(verbose && strncmp(str2,look,100) == 0){
	if(verbose){ // Do full alignment
		unsigned char *fqsp=SeqPtr(fqE),*fssp=SeqPtr(fsE);
		fprintf(stderr,"\n%c%c%c fqE (%3d-%-3d)\n",AlphaChar(fqsp[q_start],A),
                        AlphaChar(fqsp[q_start+1],A), AlphaChar(fqsp[q_start+2],A),
                        q_start,q_start+2);
		fprintf(stderr,"%c%c%c fsE (%3d-%-3d) score = %d\n",
			AlphaChar(fssp[s_start],A),AlphaChar(fssp[s_start+1],A),
			AlphaChar(fssp[s_start+2],A), s_start,s_start+2,word_score);
		// hhp = MkHSPHeap(nhits+2);
		brh_typ results;
		PerformGappedAlignmentWithTraceback(gap_align);
		hsp=MakeHSPTyp(gap_align,sbp);
		// results = MakeGBLASTResultHitlist(n+1,fsE);
		results = MakeGBLASTResultHitlist(2,ssqE);
                AddHspBRH(hsp,results);
  		sap_typ sap;
		// Int4 n=PurgeHSPHeap(hhp);  ... use to eliminate weird alignments.
        	sap = ExtractAlnBRH(results,fqE);
		PutGSeqAlignList(stderr, sap, 60, A);
		fprintf(stderr,"Likelihood(%d) = %.2f\n",j,j_like);
		results=NilBRH(results);
		FreeGSeqAlignList(sap);
	}
				} else {
				}
				NilSeq(ssqE);
			     }
			   }
			}
			if(verbose && nhits) fprintf(stderr,"\n");
			GBLAST_ScoreBlkDestruct(sbp); GABDelete(gap_align);

// PutHist(stdout,60,HG); NilHist(HG);
			// sum = (Int4)(0.5+((double)sum/(double)(N-1)));
			if(nhits > 0) sum = (Int4)(0.5+((double)sum/(double)nhits));
			else sum=0;
			double  **colfreq = ColResFreqsCMSA(1, cma);
			Int4 expscore = (Int4) ExpScoreSMatrix(colfreq,smx);
			for(j=1; j <= lenM; j++) free(colfreq[j]); free(colfreq);
			double adj_like = ug_like + 0.30103*(2*expscore);
			// sum -= ExpScoreSMatrix(smx);
			NilSMatrix(smx);

	// OUTPUT results...
			fprintf(stdout,"%3d: %s ",n2,look);
			// if(n2 > 1) fprintf(stdout," (%3d(%3d)) ",gap0,gap1);
			if(n2 > 1) fprintf(stdout," (%3d) ",gap0);
			else fprintf(stdout,"       ");
			char c1,c2,c4;
			if(g_like < look_cut && nhits < 1 && (gap0 > 20 && gap1 > 20)) {
				value[v] = n2; v++; 
			} 
			char strikes=0;
			if(g_like < look_cut){ c1='*'; strikes +=2; }
			else if(g_like < (look_cut + 2.0)){ c1='?'; strikes ++; } else c1=' ';
			fprintf(stdout,"%c",c1); 
			if(gap0 > 20 && gap1 > 20){ c2='*'; strikes +=2; }
			else if(gap0 > 20 || gap1 > 20){ c2='?'; strikes ++; } else c2=' ';
			fprintf(stdout,"%c",c2); 
			// if(nhits < 1){ c3='*'; strikes +=2; }
			// else if(nhits < 2){ c3='?'; strikes ++; } else c3=' ';
			// fprintf(stdout,"%c",c3); 
			if(adj_like < -2.0){ c4='*'; strikes +=2; }
			else if(adj_like < 0.0){ c4='?'; strikes ++; } else c4=' ';
			fprintf(stdout,"%c%d ",c4,strikes); 
			fprintf(stdout,"%4d-%-4d ", start,end);
			// fprintf(stdout,"{|%d(%d)|} ",OffSetSeq(qE),CtermExtendSeq(qE));
			fprintf(stdout,"[%.2f] (raw: %.2f; adj: %.2f)",g_like,ug_like,adj_like);
			if(nhits > 1) fprintf(stdout," %d hits [%.2f](%.2f) ",
					nhits,max_like,-log10(max_evalue));
			else if(nhits > 0) fprintf(stdout," 1 hit [%.2f](%.2f) ",
					max_like,-log10(max_evalue));
			else fprintf(stdout," no hits ");
			fprintf(stdout,"(%d[%d])\n",sum,expscore);
		   } 
		} 
		if(v){
		  fprintf(stdout,"-x%s:",look);
		  for(i=0; i < v-1; i++) fprintf(stdout,"%d,",value[i]);
		  fprintf(stdout,"%d\n",value[i]); 
	   } free(value);
	   if(total_internal) fprintf(stdout," total_internal = %d; repeats = %d\n",
			total_internal,n2);
	} 
        if(cma) TotalNilCMSA(cma);
        return 0;
}

