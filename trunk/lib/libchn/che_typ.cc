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

#include "che_typ.h"
#include "blosum62.h"

//******************** Start of che_typ routines **********************

che_typ::che_typ(char *sst_str,sst_typ **sst_in,chn_typ *chain,double A0, double B0, 
	BooLean verbose,char mode,double **Rho,double PriorRi,BooLean UseGlobalSqWts)
{
	OwnChn=FALSE; chn = chain; Sigma=0.8; 
	Init(sst_str,sst_in,A0,B0,verbose,mode,Rho,PriorRi,0,UseGlobalSqWts); 
}

che_typ::che_typ(char *sst_str, chn_typ *chain,double A0, double B0, 
	BooLean verbose,char mode,double **Rho,double PriorRi,BooLean UseGlobalSqWts)
//*************** mcBPPS constructor *********************
{ OwnChn=FALSE; chn = chain; Sigma=0.8; 
Init(sst_str,0,A0,B0,verbose,mode,Rho,PriorRi,0,UseGlobalSqWts); }

che_typ::che_typ(char *sst_str, chn_typ *chain,BooLean verbose,char mode, char Type, 
		BooLean UseGlobalSqWts,double PriorRi)
//*************** omcBPPS constructor *********************
{ OwnChn=FALSE; chn = chain; Sigma=0.8; 
Init(sst_str,0,0,0,verbose,mode,0,PriorRi,Type,UseGlobalSqWts); }

void	che_typ::Init(char *sst_str,sst_typ **sst_in,double A0,double B0,
	BooLean verbose,char Mode,double **Rho,double PriorRi, char Type ,BooLean UseGlobalSqWts)
{
	Int4	i,j,s,r;

	// use average sequence weights so that MAP contribution more consistent with Urn model
	// afn: 6-27-08
	if(UseGlobalSqWts) GlobalSqWts=TRUE; else GlobalSqWts=FALSE;
	Verbose=verbose;
	debug_fp=0; // debug_fp=open_file("debug",".err","w");
	MinNumColumns=5;
	MaxNumColumns=1000;
	FIXED_PARTITIONS=FALSE;
	unsigned char del_as_random=0;
	if(islower(Mode)){ del_as_random=1; Mode=toupper(Mode); }
        num_sets=2;	// change the number of sets (partitions) here....Other changes also necessary!
	assert(chn->GetNumAnalysis() > 0);
	AB=chn->GetAlphabet();

	// 1. Make basic assignments from chn object.
	IN_CMA = chn->GetIN_CMSA( );  // get all input cma files.
	qcma=chn->GetIN_CMSA(1); assert(nBlksCMSA(qcma) == 1);
	mcma=chn->GetIN_CMSA(2); assert(nBlksCMSA(mcma) == 1);
	assert(LengthCMSA(1,qcma) == LengthCMSA(1,mcma));
        Int4 start = TruePosCMSA(1,1,qcma);      // ignores offset...
        Int4 end = start+LengthCMSA(1,qcma)-1;   // will jump over gaps...
        Query=MkSubSeq(start,end,SeqSetE(1,DataCMSA(qcma)));  // Fake sequence...
	assert(LengthCMSA(1,qcma) == LenSeq(Query));
	Contrast=(Int4) ceil((double)LenSeq(Query)/10.0);
	Contrast=MAXIMUM(Int4, Contrast,5);	// default contrast setting.
	// assert(IdentSeqs(Query,chn->QuerySeq()));
	NumAnalysis = chn->GetNumAnalysis();
	assert(NumAnalysis==2);

        NEW(qsst,LenSeq(Query)+2,sst_typ);
        if(sst_in){
		for(i=1;i <= LenSeq(Query); i++){ qsst[i]=sst_in[i][0]; }
	} else if(sst_str){
	  Int4		nval,*pos=0;
	  sst_typ	*sst=0;
          NEW(sst,strlen(sst_str)+2,sst_typ);
          NEW(pos,strlen(sst_str)+3,Int4);
          nval=ParseResidueSets(sst_str,pos,sst,AB,"che_typ Init() input error 1");
          // nval=ParseResidueSets(sst_str,pos,sst,AB,"che_typ Init() input error 1");
          for(i=1;i <= nval; i++){
		s = RealToFakeCMSA(1,pos[i],qcma); // s == actual position in array.
		// s == fake seq position for residue s in sq (with offset).
		if(s <= 0){
		   fprintf(stderr,"sst_str=%s\n",sst_str);
		   fprintf(stderr,"Position %d(%d):\n",pos[i],s);
		   print_error("che_typ input pattern option error");
		}
                if(s > LenSeq(Query)){ print_error("che_typ Init() input option error 2"); }
                r=ResSeq(s,Query);
                // fprintf(stdout,"%c%d(%d)\n",AlphaChar(r,AB),pos[i],s);
		if(!MemSset(r,sst[i])){
			fprintf(stderr,
			  "**************** FATAL ERROR! **************\n");
                	fprintf(stderr,"%c%d(%d)\n",AlphaChar(r,AB),pos[i],s);
			fprintf(stderr,"sst[%d]: \"",pos[i]);
			PutSST(stderr,sst[i],AB);
			fprintf(stderr,"\"\n");
			PutSeq(stderr,Query,AB);
			gsq_typ *gsq=gsqCMSA(1,qcma);
			gsq->Put(stderr,AB);
			gsq->Put(stderr,60,AB);
			fprintf(stderr,"%d(%d): %c\n",pos[i],s,AlphaChar(r,AB));
			fprintf(stderr,"Query offset = %d\n",OffSetSeq(Query));
			fprintf(stderr,"%s\n",sst_str);
			print_error("che_typ: pattern (sst) input error\nFor mcBPPS try -nocsq option.\n\n");
		}
		else qsst[s]=sst[i];
          } free(sst); free(pos); // exit(1);
        } else {
		for(i=1;i <= LenSeq(Query); i++){ r=ResSeq(i,Query); qsst[i]=SsetLet(r); }
	}

	// 3. Set objects for new cma files.
	if(Verbose){
	    // char str[15]; StrSeqID(str, 12,Query);
	    // fprintf(stderr,"Set up objects for new cma files (%s).\n",str);
	    fprintf(stderr,"Set up objects for new cma files (%s).\n",NameCMSA(qcma));
	}
	num_seqs=NumSeqsCMSA(mcma);
	gold_set=MakeSet(num_seqs+1);
	ClearSet(gold_set); 
	NEW(column,LenSeq(Query)+3,char);
#if 1	// eventually move below...
	NEW(card_set,num_sets + 3,UInt4);
	NEW(Card_Set,num_sets + 3,UInt4);
#endif
        NEWP(ResWt,5,UInt4 *); 
	for(s=1; s<=LenSeq(Query); s++) column[s]=0;
	set_cols[0]=LenSeq(Query);
	for(s=1; s<=num_sets; s++){ 
		set_cols[s]=0;
		card_set[s]=0;
		map_set[s]=0;
		info_set[s]=0;
		seq_set[s]=MakeSet(num_seqs+1); ClearSet(seq_set[s]);
	}

#if 1	// Added for Recursively nested multiple constraint sampler (afn: 11/12/09).
	RmSet = MakeSet(num_seqs+1); ClearSet(RmSet);
#endif

	GetIntegerWts( );	// get AveSqIWt...
	Int4 set;
	// Fix problem with sequence ordering...
        // Find gold_set (qcma) within mcma with purging.

#if 1	// eventually turn this off and do this below...
        for(set=1; set <= num_sets; set++){
          NEWP(ResWt[set],LenSeq(Query)+3,UInt4);
          for(s=1; s <= LenSeq(Query); s++){
                NEW(ResWt[set][s],nAlpha(AB)+3,UInt4);
          }
        }
#endif
	// Fix bug that adds first seq in mcma by default...
	// fix by treating first seq just like the others in display set...
	// This was due to chn_aln adding a concensus seq. to the main cma.
	UInt4 wt;
	unsigned char *seq;
        e_type  qE,mE;
        for(i=1; i<=NumSeqsCMSA(mcma); i++){
           mE=TrueSeqCMSA(i,mcma);
           for(j=1; j<=NumSeqsCMSA(qcma); j++){
                qE=TrueSeqCMSA(j,qcma);
                if(IsSubSeq(qE,mE)) { 
		   AddSet(i,gold_set); AddSet(i,seq_set[1]);
#if 0	// eventually turn this off
		   card_set[1]++;
		   Card_Set[1]+= AveSqIWt[i];	// i is index for mcma...
		   // WARNING: OVERCOUNTS SEQUENCES WHEN SUBSEQS IN DISPLAY SET!!!
		   // NEED TO FIX THIS AT SOME POINT...SQ_WTS MAY BE WRONG!!!
		   // NEW: integer weights.
		   seq=GetAlnResInSiteCMSA(1,i,mcma);
		   for(s=1; s <= LenSeq(Query); s++){
			if(GlobalSqWts) wt=AveSqIWt[i]; else wt = SqWt[s][i];
			ResWt[1][s][seq[s]] += wt;
		   }
#endif
		   break; 
		} 
           } 
        }

	// NEW: sort sequences by pairwise pseudoscore to chn concensus seq...
	Int4 score;
        dh_type dH = dheap(NumSeqsCMSA(mcma)+3,4);
        for(i=1; i <= NumSeqsCMSA(mcma); i++){
	      if(Mode == 'R'){		// 'R' = initialize partition randomly
		score=Random();
	      } else if(Mode == 'O'){	// 'O' = initialize partition by input order
            	score = -i;
	      } else if(Mode == 'S'){	// 'S' = initialize partition by score to query.
            	score = PseudoAlnScoreCMSA(i,1,mcma);
	      } else print_error("che_typ input error: mode != 'R','S' or 'O'");
              // fprintf(stderr,"%d vs %d: %d\n",i,j,score);
              insrtHeap(i,-score,dH);
              // fprintf(stderr,"%d\n",i);
	}
	double partition_size=((double)num_seqs/(double)num_sets);
	// Partition sequences into distinct sets...
	// e.g., partition_size=1000/3 = 333.333;
        for(i=1; !emptyHeap(dH) && (j=delminHeap(dH)); i++){
	    if(!MemberSet(j,gold_set) && !MemberSet(j,seq_set[1]) && 
			!MemberSet(j,seq_set[2]))
		{
		// e.g., set = ceil(300/333.33)== 1);
		set=(Int4)ceil((double)i/partition_size);
		if(set > num_sets) set = num_sets;
	        AddSet(j,seq_set[set]); 
#if 0	// Eventually turn this off...
		card_set[set]++;
		Card_Set[set]+=AveSqIWt[j];	// integer weights.
		seq=GetAlnResInSiteCMSA(1,j,mcma);
		for(s=1; s <= LenSeq(Query); s++){
			if(GlobalSqWts) wt=AveSqIWt[j]; else wt = SqWt[s][j];
			ResWt[set][s][seq[s]] += wt;
			// SumWts[1][s] += wt;
		}
#endif
	    } // else Card_Set[1] += AveSqIWt[j];  // done above...
	} Nildheap(dH);

#if 1	// Debug inconsistency between card_set/Card_Set and seq_set;
        for(s=1; s <= num_sets; s++){
	   card_set[s] = CardSet(seq_set[s]); Card_Set[s]=0;
	}
        for(set=1; set <= num_sets; set++){
          for(s=1; s <= LenSeq(Query); s++){ free(ResWt[set][s]); } free(ResWt[set]);
          NEWP(ResWt[set],LenSeq(Query)+3,UInt4);
          for(s=1; s <= LenSeq(Query); s++){ NEW(ResWt[set][s],nAlpha(AB)+3,UInt4); }
        }
        for(i=1; i <= NumSeqsCMSA(mcma); i++){
	  seq=GetAlnResInSiteCMSA(1,i,mcma);
          for(s=1; s <= num_sets; s++){
	     if(MemberSet(i,seq_set[s])){
		Card_Set[s] += AveSqIWt[i];;
		for(r=1; r <= LenSeq(Query); r++){
			if(GlobalSqWts) wt=AveSqIWt[i]; else wt = SqWt[r][i];
			ResWt[s][r][seq[r]] += wt;
			// SumWts[1][r] += wt;
		}
	     }
	  }
	}
	NEWP(TotalResWt,num_sets + 2,UInt4);
        for(s=1; s <= num_sets; s++){
            NEW(TotalResWt[s],LenSeq(Query)+3,UInt4);
	    for(i=1; i <= LenSeq(Query); i++){
		TotalResWt[s][i] = 0;
	    	// for(r=startAB; r <= nAlpha(AB); r++) TotalResWt[s][i] += ResWt[s][i][r];
	    	for(r=0; r <= nAlpha(AB); r++) TotalResWt[s][i] += ResWt[s][i][r];
	    }
	}
#endif

   if(Type != 0) { // 
#if 0	// old common default sequence priors..
	pps = new bpps_typ(qsst,SeqPtr(Query),LenSeq(Query),Type, AB,UseGlobalSqWts);
#elif 1
	pps = new bpps_typ(qsst,SeqPtr(Query),LenSeq(Query),Type, AB,UseGlobalSqWts,PriorRi);
#else
	pps = new bpps_typ(qsst,SeqPtr(Query),LenSeq(Query),Type, AB,UseGlobalSqWts,0.5);
#endif
   } else {
	// WtFactor=100;	// Sequence weighting 0..100 --> 0.0 ... 1.0.
	// WARNING: This WtFactor is set independently within swt_typ: need to fix!!!
	// NEW(beta,nAlpha(AB)+3,double);
	// for(r=startAB; r <= nAlpha(AB); r++) beta[r]=WtFactor; // default = non-informed prior.
	// for(r=1;r<=nAlpha(AB);r++) beta[r]=WtFactor*(double)nAlpha(AB)*blosum62freq[r];
	pps = new bpps_typ(qsst,SeqPtr(Query),LenSeq(Query),A0,B0,AB,Rho,PriorRi,UseGlobalSqWts);
   }
	WtFactor = pps->RtnWtFactor(); 
#if 0
	pps->TreatDeletionAsRandom();	// turn this on always...
#else
	if(del_as_random) pps->TreatDeletionAsRandom();
#endif
	if(Verbose) fprintf(stderr,"Set up randomization procedures.\n");
	sqH = dheap(num_seqs+2,4);
	mvH = dheap((num_sets*2)+2,4);
	posH = dheap(LenSeq(Query)+2,4);
	// NullMap=pps->NullMap(ResWt[2],ResWt[1]);
	if(Verbose) fprintf(stderr,"LPR = %g\n",pps->LPR(ResWt[2],ResWt[1]));
}

void	che_typ::Free()
{
	Int4	s,set;
#if 0	// diagnostics: need to pass in mcma file!!!
        if(chn->OutPutAln()) chn->PutAln( );
        else {
          chn->PutHierarchicalAln( );
          chn->PutRasmol( );
        }
#endif
	// free(beta);
	for(s=1; s <= LenSeq(Query); s++) free(SqWt[s]); free(SqWt);
	for(set=1; set <= num_sets; set++){
	   for(s=1; s <= LenSeq(Query); s++){ free(ResWt[set][s]); }
	   free(ResWt[set]); free(TotalResWt[set]);
	   // free(SumWts[set]);
	   // If use ReadCMSA then sequence data sets need to be deallocated.
	   // If use CopyCMSA then sequence data sets beInt4 to chn_typ.
	   NilSet(seq_set[set]);
	} NilSet(gold_set); free(column);
	free(ResWt); free(TotalResWt);
	Nildheap(sqH); Nildheap(mvH); Nildheap(posH);
	free(qsst);
	free(card_set);
	free(Card_Set); free(AveSqIWt);
	if(RmSet) NilSet(RmSet);
	if(Query) NilSeq(Query);	// when creating a subseq...
	if(OwnChn) delete chn;
	delete pps;
}

UInt4   *che_typ::FGResWtsBILD(Int4 s, UInt4 &totWt)
{
        Int4 r,i,nAB=nAlpha(AB); assert(s > 0 && s <= LenSeq(Query));
        UInt4 *Rtn; NEW(Rtn,nAlpha(AB) + 3, UInt4);
        for(totWt=0,r=0; r <= nAlpha(AB); r++){ Rtn[r]=ResWt[1][s][r]; totWt += Rtn[r]; }
        UInt4 diff=TotalResWt[1][s]-totWt;
	UInt4 inc=diff/nAB;
	if(inc > 0) for(r=1; r <= nAB; r++){ Rtn[r] += inc; }
	UInt4 rmn=diff%nAB;
	if(rmn > 0) for(i=1; i <= rmn; i++){ r=(i%nAB)+1; Rtn[r]++; }
        return Rtn;
}

UInt4   *che_typ::BGResWtsBILD(Int4 s, UInt4 &totWt) 
{
        Int4 r,i,nAB=nAlpha(AB); assert(s > 0 && s <= LenSeq(Query));
        UInt4 *Rtn; NEW(Rtn,nAlpha(AB) + 3, UInt4);
        for(totWt=0,r=0; r <= nAlpha(AB); r++){ Rtn[r]=ResWt[2][s][r]; totWt += Rtn[r]; }
        UInt4 diff=TotalResWt[2][s]-totWt;
	UInt4 inc=diff/nAB;
	if(inc > 0) for(r=1; r <= nAB; r++){ Rtn[r] += inc; }
	UInt4 rmn=diff%nAB;
	if(rmn > 0) for(i=1; i <= rmn; i++){ r=(i%nAB)+1; Rtn[r]++; }
        return Rtn;
}

void    che_typ::PutSeqs()
{
	char	str[100];
	for(Int4 s=1; s <=num_sets; s++){
		sprintf(str,"_new%d.seq",s);
		FILE *fp=open_file(FileName(),str,"w");
 		PutSeqsInSetCMSA(fp,seq_set[s],mcma); 
		fclose(fp);
	}
}

void    che_typ::PutInfo( )
{
	FILE *fptr=open_file(FileName(),"_new.info","w");
	PutInfo(fptr); fclose(fptr);
}

void    che_typ::PutInfo(FILE *fp) { PutInfo(fp,-1); }

void    che_typ::PutInfo(FILE *fp,Int4 id)
{
	Int4	s;
#if 1
	char	*str=strrchr(FileName(),'/');
	if(str==NULL) str=FileName();
	else str++;
	
	if(id < 0) fprintf(fp,"%s: (",str);
	else fprintf(fp,"%s_new%d: (",str,id);
#else
	if(id < 0) fprintf(fp,"%s: (",FileName());
	else fprintf(fp,"%s_new%d: (",FileName(),id);
#endif
   	for(s=1; s <= num_sets; s++) fprintf(fp,"%d ",card_set[s]);
	fprintf(fp," seqs)(%g; %.1f nps; %.1f npws)(%d cols)\n",
		Map,Map/(double)card_set[1],this->RtnNatsPerWtSq(),pps->NumColumns());
}

double  che_typ::RtnAveWtSeqs( )
{
	assert(GlobalSqWts==TRUE);
	double	D=0.0;
	for(Int4 s=1; s <= LenSeq(Query); s++) D += TotalResWt[1][s];
	return D/(double)LenSeq(Query);     // Assumes using AveResWt.
}

double  che_typ::RtnNatsPerWtSq( )
// Caution: this assumes that LPR has been calculated and that using AveResWt.
{ return (pps->RtnWtFactor()*Map)/this->RtnAveWtSeqs(); }

double	che_typ::PutInfoShort(FILE *fp)
{
	fprintf(fp,"(");
   	for(Int4 s=1; s <= num_sets; s++) fprintf(fp,"%d ",card_set[s]);
	char c=' '; if(Map <= 0.0 || pps->NumColumns() ==0) c = 'x';
	fprintf(fp," seqs)(%g; %.1f nps; %.1f npws)(%d cols) %c\n",
		Map,Map/(double)card_set[1],this->RtnNatsPerWtSq(),pps->NumColumns(),c);
	return Map;
}

void    che_typ::PutPattern(FILE *fp) {	PutPattern(fp,0); }

void    che_typ::PutPattern(FILE *fp,Int4 *fake2real)
{
	Int4 j;

	fprintf(fp,"Pattern sets:\n");
	for(j=1; j <= LenSeq(Query); j++){
	   if(pps->RtnSST(j)){
	     PutSST(fp,pps->RtnSST(j),AB); 
	     if(fake2real) fprintf(fp,"%d,",fake2real[j]);
	     else fprintf(fp,"%d,",j);
	   }
	} fprintf(fp,"\n");
	
}

void    che_typ::PutVSI_Info(FILE *fp,Int4 *fake2real)
{
	Int4 j;

	fprintf(fp,"# Residues found for vsi file:\n");
	for(j=1; j <= LenSeq(Query); j++){
	   if(pps->RtnSST(j)){
	     fprintf(fp,"%c",AlphaChar(j,AB));
	     if(fake2real) fprintf(fp,"%d.Y\n",fake2real[j]);
	     else fprintf(fp,"%d.Y\n",j);
	   }
	} fprintf(fp,"\n");
	
}

void    che_typ::ContribSeqMAP(Int4 id,char SFBG,double *min_contrib,double min_nats,double min_freq)
// output contribution of each sequence to the MAP
// Also finds intermediate sequences and discards them.
{
	char	str[200];
	FILE	*fptr=stderr,*sfp=NULL;
	Int4	sq,sq_set,num_sq_removed;
	Int4	num_outliers,num_marginal;
	double	map,map0,delta_map;

	assert(num_sets==2);
	sprintf(str,"_new%d.contrib",id);
	fptr=open_file(FileName(),str,"w");
	Int4	mv; 
	map0=pps->LPR(ResWt[2],ResWt[1]);	// 
	dh_type dH = dheap(num_seqs+3,4);
	num_sq_removed=0;
	for(sq_set=1; sq_set <=2; sq_set++){
	  h_type HG=Histogram("sequence contributions to the map",-10,100,0.25);
	  fprintf(fptr,"Partition %d:\n",sq_set);
	  num_outliers=0; num_marginal=0;
	  for(sq=1; sq <= num_seqs; sq++){

	   if(MemberSet(sq,RmSet)) continue; // Ignore removed sequences (afn: 11/12/09).

	   if(!MemberSet(sq,seq_set[sq_set])) continue;
           if(MemberSet(sq,seq_set[2])) mv=-1;	// move down one.
	   else mv=1;	// move up one.
	   assert(SwapSeq(sq,mv));
	   map=pps->LPR(ResWt[2],ResWt[1]); 
	   delta_map=map0 - map;
	// WARNING: this is not working when ModelMode=='B'!! Need to fix this.
	// TEMP fix == prohibit removing marginal sequences in 'B' mode.
	   IncdHist(delta_map,HG);
	   insrtHeap(sq,(keytyp)delta_map,dH);
#if 0
	   e_type E=TrueSeqCMSA(sq,mcma);
	   StrSeqID(str, 60, E);
	   fprintf(fptr,"%d (%s): delta map=%.3f\n",sq,str,delta_map);
#endif
	   SwapSeq(sq,-mv); // Switch back to optimum configuration.
	  } fprintf(fptr,"\n");
	  double mean = MeanHist(HG);
	  double variance = VarianceHist(HG);
	  double stdev=sqrt(variance);
	  double outlier=mean+4.0*stdev;
	  double inlier=mean-stdev;
	  if(inlier < 0.25) inlier=0.25;
	  PutHist(fptr,60,HG); NilHist(HG);
#if 1	// use min_contrib as fraction of seqs to ignore (up to 50%)
	// this is to eliminate pseudogenes
	// ignore the first 
	  assert(min_contrib[sq_set] <= 0.9);
	  Int4 Ignored=(Int4) floor(min_contrib[sq_set]*(double) ItemsInHeap(dH));
	  if(Ignored >= ItemsInHeap(dH)){
		fprintf(stderr,"min_contrib[sq_set] = %f; sq_set=%d\n",
			min_contrib[sq_set],sq_set);
		fprintf(stderr,"WARNING: Ignored = %d >= ItemsInHeap(dH) = %d\n",
			Ignored,ItemsInHeap(dH));
	  	// assert(Ignored < ItemsInHeap(dH));
	  }
	  if(Ignored > 0) {
	   sprintf(str,"_ignore%d.%dsq",id,sq_set);
	   sfp=open_file(FileName(),str,"w");
	  }
#endif
	  Int4 nsq_seen=0;
	  while(!emptyHeap(dH)){
		delta_map=(double)minkeyHeap(dH);
		sq=delminHeap(dH);
#if 1		// remove most ambigous sequences to ignore
		nsq_seen++;
                if(nsq_seen <= Ignored && !MemberSet(sq,gold_set)){
			if(num_marginal < 1) fprintf(fptr,"\n==== Marginal: ====\n");
			num_marginal++;
	   		e_type E=TrueSeqCMSA(sq,mcma);
	   		StrSeqID(str,60,E);
	   		fprintf(fptr,"%d (%s): delta map=%.3f\n",sq,str,delta_map);
			DeleteSet(sq,seq_set[sq_set]); num_sq_removed++;
			if(sfp != NULL) PutSeq(sfp,E,AB);
		}
	        if(delta_map > outlier){
			if(num_outliers < 1) fprintf(fptr,"\n==== Outliers: ====\n");
			num_outliers++;
			e_type E=TrueSeqCMSA(sq,mcma);
	   		StrSeqID(str,60,E);
	   		fprintf(fptr,"%d (%s): delta map=%.3f\n",sq,str,delta_map);
		}
#else
		// remove intermediate sequences from the partition unless in Display set.
		if(delta_map < min_contrib || delta_map < inlier){ 
			if(num_marginal < 1) fprintf(fptr,"\n==== Marginal: ====\n");
			num_marginal++;
	   		e_type E=TrueSeqCMSA(sq,mcma);
	   		StrSeqID(str,60,E);
	   		fprintf(fptr,"%d (%s): delta map=%.3f\n",sq,str,delta_map);
			// remove from both partitions
                	if(delta_map < min_contrib && !MemberSet(sq,gold_set)){
				DeleteSet(sq,seq_set[sq_set]); num_sq_removed++;
			}
		} else if(delta_map > outlier){
			if(num_outliers < 1) fprintf(fptr,"\n==== Outliers: ====\n");
			num_outliers++;
			e_type E=TrueSeqCMSA(sq,mcma);
	   		StrSeqID(str,60,E);
	   		fprintf(fptr,"%d (%s): delta map=%.3f\n",sq,str,delta_map);
		}
#endif
	  }
	} Nildheap(dH);
	if(sfp != NULL) fclose(sfp); 
	fprintf(fptr,"\n"); fclose(fptr);
	if(num_sq_removed > 0){
	   sprintf(str,"_new%d.chn",id);
	   fptr=open_file(FileName(),str,"w");
	   Put(fptr,SFBG,min_nats,min_freq); fclose(fptr);
	}
}

void    che_typ::PutChn(const char *str, char SFBG,double min_nats,double min_freq)
{
	FILE *fptr=open_file(FileName(),str,"w");
	Put(fptr,SFBG,min_nats,min_freq); fclose(fptr);
}

void    che_typ::PutSubLPRs(char *filename)
{
	FILE *fptr=open_file(filename,"","w");
	PutSubLPRs(fptr); fclose(fptr);
}

void    che_typ::PutSubLPRs(FILE *fp)
{
	Int4 *fake2real;

	NEW(fake2real,LenSeq(Query)+4,Int4);
	for(Int4 s=1; s <= LenSeq(Query); s++){
		fake2real[s]=FakeToRealCMA(1,s,qcma); 
	}
#if 0	// compute values for arrays...fix this later.
	pps->SubLPR(ResWt[2],ResWt[1]);
#endif
	assert(num_sets==2);
	fprintf(fp,"Foreground: %d; ",card_set[1]);
	fprintf(fp,"Background: %d; ",card_set[2]);
	fprintf(fp,"LPR: %g (%d cols)\n\n\n",Map,pps->NumColumns());
	pps->PutSubLPR(fp,fake2real,ResWt[2],ResWt[1]);
	free(fake2real);
}

void    che_typ::WriteSubLPR(FILE *fptr)
{
	Int4 *fake2real;

	NEW(fake2real,LenSeq(Query)+4,Int4);
	for(Int4 s=1; s <= LenSeq(Query); s++){ fake2real[s]=FakeToRealCMA(1,s,qcma); }
	pps->WriteSubLPR(fptr,fake2real,ResWt[2],ResWt[1]);
	// Note: this won't work for multinomial version of pps; only for bpps version.
	free(fake2real);
}

void    che_typ::PutBestPatterns(FILE *fptr,BooLean UseRealPos)
{
	Int4 *fake2real=0;
	if(UseRealPos){
	   NEW(fake2real,LenSeq(Query)+4,Int4);
           for(Int4 s=1; s <= LenSeq(Query); s++){ fake2real[s]=FakeToRealCMA(1,s,qcma); }
	}
	pps->PutPattern(fptr,fake2real,ResWt[2],ResWt[1]);
	if(fake2real) free(fake2real);
}

void    che_typ::PutFgBgSeqs(Int4 id)
{
	char	str[200];

	for(Int4 s=1; s <=num_sets; s++){
		sprintf(str,"_new%d.%dsq",id,s);
		FILE *fptr=open_file(FileName(),str,"w");
 		PutSeqsInSetCMSA(fptr,seq_set[s],mcma); 
		fclose(fptr);
	}
}

void    che_typ::PutPttrnFile(FILE *fptr,Int4 id)
{
	Int4 *fake2real; NEW(fake2real,LenSeq(Query)+4,Int4);
	for(Int4 s=1; s <= LenSeq(Query); s++){ fake2real[s]=FakeToRealCMA(1,s,qcma); }
	pps->PutSubLPR(fptr,fake2real,ResWt[2],ResWt[1]);
	free(fake2real);
}

void    che_typ::PutAll(Int4 id,char SFBG,double min_nats,double min_freq)
{
	char	str[200];

	sprintf(str,"_new%d.info",id);
	FILE *fptr=open_file(FileName(),str,"w");
	// fprintf(fptr,"(%d)",id);
        PutInfo(fptr,id); fclose(fptr);

	sprintf(str,"_new%d.pttrn",id);
	fptr=open_file(FileName(),str,"w");
#if 0	// replace with this after testing...
void    che_typ::PutSubLPR(FILE *fptr);	// need to make this routine...
#elif 0
	Int4 *fake2real;
	NEW(fake2real,LenSeq(Query)+4,Int4);
	for(Int4 s=1; s <= LenSeq(Query); s++){ fake2real[s]=FakeToRealCMA(1,s,qcma); }
	pps->PutSubLPR(fptr,fake2real,ResWt[2],ResWt[1]);
	// PutPattern(fptr,fake2real); 
	free(fake2real);
#else
	PutPttrnFile(fptr,id);
#endif
	fclose(fptr);

	sprintf(str,"_new%d.chn",id);
	PutChn(str,SFBG,min_nats,min_freq);

#if 1	// Create omcBPPS files.
	sprintf(str,"_new%d.sma",id);
	fptr=open_file(FileName(),str,"w");
	cma_typ xcma=MakeConsensusCMSA(mcma); RenameCMSA("Set1",xcma);
	PutCMSA(fptr,xcma); RenameCMSA("Set2",qcma); PutCMSA(fptr,qcma); 
	fclose(fptr); TotalNilCMSA(xcma);

	sprintf(str,"_new%d.hpt",id);
	fptr=open_file(FileName(),str,"w");
	Int4 i,M,n,N=1+(NumSeqsCMSA(mcma)/3);
	fprintf(fptr,"HyperParTition:\n!!\n+- 1.Set1?\n++ 2.Set2!\n-o 3.Random=%d.\n\n",N);
	fclose(fptr);
	
	set_typ sets[20]; sets[0]=0; 
	M=SetN(seq_set[1]) + N;
	sets[1]=MakeSet(M); ClearSet(sets[1]); sets[2]=CopySet(sets[1]); sets[3]=CopySet(sets[1]); 
	for(n=1; n <= NumSeqsCMSA(mcma); n++){
		if(MemberSet(n,seq_set[2])) AddSet(n,sets[1]);
		else if(MemberSet(n,seq_set[1])) AddSet(n,sets[2]);
	} 
	for(i=1; i <= N; n++,i++) AddSet(n,sets[3]); 

	sprintf(str,"_new%d.sets",id);
	fptr=open_file(FileName(),str,"w");
	WriteSets(fptr,3,sets);
	fclose(fptr); NilSet(sets[1]); NilSet(sets[2]); NilSet(sets[3]);
#endif
#if 0	// shut this off...
	for(Int4 s=1; s <=num_sets; s++){
		sprintf(str,"_new%d.%dsq",id,s);
		fptr=open_file(FileName(),str,"w");
 		PutSeqsInSetCMSA(fptr,seq_set[s],mcma); 
		fclose(fptr);
	}
#endif
}

void    che_typ::PutDS(FILE *fp,double min_nats,double min_freq)
// Output display set...
{
	PutCMSA(fp,IN_CMA[1]);
}

void    che_typ::PutDSubSets(FILE *fp,double min_nats,double min_freq)
// Output display subset...
{
	for(Int4 i=3; i <= chn->GetNumInCMSA(); i++){
		cma_typ cma=chn->GetIN_CMSA(i);
		PutCMSA(fp,cma);
	}
}

void    che_typ::PutFG(FILE *fp,double min_nats,double min_freq)
{
	sst_typ *mq_sst = pps->ModifySST(Contrast,min_nats,min_freq);
	PutInSetCMSA(fp,seq_set[1],mq_sst,mcma);
	free(mq_sst);
}

void    che_typ::PutBG(FILE *fp,double min_nats,double min_freq)
{
	// PutInSetCMSA(fp,seq_set[2],mcma);
	set_typ Set=MakeSet(num_seqs+1); ClearSet(Set);
	UnionSet(Set,seq_set[2]);
	// PutSelectCMSA(0,fp,0,"BackGround",Set,mcma);
	PutSelectCMSA(fp,0,"BackGround",0,Set,mcma);
	NilSet(Set);
}

void    che_typ::PutMainSet(FILE *fp,double min_nats,double min_freq)
{
	set_typ Set=MakeSet(num_seqs+1); FillSet(Set);
	//              FG,BG,
	// PutSelectCMSA(0,fp,0,"MainSet",Set,mcma);
	PutSelectCMSA(fp,0,"MainSet",0,Set,mcma);
}

void    che_typ::Put(FILE *fp,char SFBG,double min_nats,double min_freq)
// output current contrast heirarchical alignment.
{
	Int4	i;

	set_typ	Set=MakeSet(num_seqs+1); ClearSet(Set); 
	
	PutCMSA(fp,IN_CMA[1]);
	// New routine return only positions in the pattern
	// with >= min_nats (default = 5.0) nats to eliminate 
	// pattern positions with insignificant contributions.
	sst_typ *mq_sst = pps->ModifySST(Contrast,min_nats,min_freq);
	PutInSetCMSA(fp,seq_set[1],mq_sst,mcma);
	free(mq_sst);
	// PutInSetCMSA(fp,seq_set[1],mcma);
	// fprintf(stderr,"GetNumInCMSA = %d\n",chn->GetNumInCMSA());
	for(i=3; i <= chn->GetNumInCMSA(); i++){
		cma_typ cma=chn->GetIN_CMSA(i);
		PutCMSA(fp,cma);
	}
	if(SFBG == 'A') PutSelectCMSA(0,fp,0,"BackGround",Set,mcma); 
	// background for Family.
	// NEW: Cumulative sequence alignments as background.
	if(SFBG != 'B') UnionSet(Set,seq_set[1]);
	for(i=2; i <= num_sets; i++){
	   if(SFBG == 'B') PutInSetCMSA(fp,seq_set[i],mcma); 
	   else { UnionSet(Set,seq_set[i]); PutInSetCMSA(fp,Set,mcma); }
	   if(SFBG == 'A'  && i < num_sets) PutSelectCMSA(0,fp,0,"BackGround",Set,mcma);
	} NilSet(Set);
}

double	che_typ::Sample(BooLean	RmCols,BooLean Lock, sst_typ **res_sst)
// 3. Perform Gibbs sampling on patterns and sequence set assignments.
{
	double	map=-9999.0,*sub_map;
	Int4	j;
	BooLean	RmAllNeg=TRUE;
#if 0	// KEEP FOR DEBUGGING...
	if(res_sst){
	   for(j=1; j <= LenSeq(Query); j++){
		fprintf(stdout,"%d: ",j);
		for(Int4 s=1; res_sst[j][s]; s++){
			fprintf(stdout,"{");
		 	PutSST(stdout,res_sst[j][s],AB);
			fprintf(stdout,"} ");
		} fprintf(stdout,"\n");
	   } // exit(1);
	}
#endif
	// fprintf(stderr,"LPR = %g\n",map=pps->LPR(ResWt[2],ResWt[1])); 
	// sub_map=pps->SubLPR(ResWt[2],ResWt[1]);
	sub_map=SampleSeqs(map);
	map=sub_map[0];
	if(!Lock && !RmCols){
	   for(j=1; j <= LenSeq(Query); j++){
	          if(AddBestColumn(map,j,res_sst[j])){ } 
		  map=sub_map[0];
	   }
	   RmCols=TRUE;
	}
	if(Lock){
		fprintf(stderr,"LPR = %g\n",map=pps->LPR(ResWt[2],ResWt[1]));
		sub_map=SampleSeqsOnce(map);
		map=sub_map[0];
	} else {
	  while(RmCols && RmWorstNegColumn(RmAllNeg,sub_map)){
	     // fprintf(stderr,"LPR = %g\n",map=pps->LPR(ResWt[2],ResWt[1])); 
	     sub_map=SampleSeqsOnce(map);
	     map=sub_map[0];
	     BooLean	ColumnAdded=FALSE;
    	     if(res_sst) // Added to fix chn_ebc core dump!
	        for(j=1; j <= LenSeq(Query); j++){
	          if(AddBestColumn(map,j,res_sst[j])){
			ColumnAdded=TRUE;
			map=sub_map[0];
		  } else map=sub_map[0];
		}
		if(ColumnAdded){
			sub_map=SampleSeqsOnce(map);
		  	map=sub_map[0];
		}
	  }
	}
	if(res_sst){
	   BooLean again=TRUE;
	   Int4		iters=0;
	   double	old_map,Difference;
	   while(again){
	    iters++;
	    again=FALSE;
	    for(j=1; j <= LenSeq(Query); j++) insrtHeap(j,(keytyp)Random(),posH);
	    while((j=delminHeap(posH)) != 0){
	     if((iters >= 2 && res_sst[j]) || !qsst[j]) {
		// Getting stuck in this loop for model = 'U'
	       if(!Lock) while(AddBestColumn(map,j,res_sst[j])){
		if(Verbose) pps->PutSubLPR(stdout,ResWt[2],ResWt[1]);
		again=TRUE;
		old_map=map;
	        if(Verbose){
		  fprintf(stderr,"OldMap = %g\n",old_map);
		  fprintf(stderr,"LPR = %g\n",map=pps->LPR(ResWt[2],ResWt[1]));
		  fprintf(stderr,"Difference = %g\n",map-old_map);
	 	}
	        sub_map=SampleSeqsOnce(map);
		// sub_map=SampleSeqs(map);
		if(map < sub_map[0]) map=sub_map[0];	// possible rounding errors...
		// map=sub_map[0];
		// Remove columns that are now negative.
		while(RmCols && RmWorstNegColumn(RmAllNeg,sub_map)){
		  if(Verbose) pps->PutSubLPR(stdout,ResWt[2],ResWt[1]);
		  if(map < sub_map[0]) map=sub_map[0]; // possible rounding errors...
		  fprintf(stderr,"LPR = %g\n",map=pps->LPR(ResWt[2],ResWt[1]));
	          sub_map=SampleSeqsOnce(map);
		  // sub_map=SampleSeqs(map);
		  if(map < sub_map[0]) map=sub_map[0];	// possible rounding errors...
		  // map=sub_map[0];
	        }
	       }
	     }	// stuck in this loop for urn model...
    	    }
	    sub_map=SampleSeqs(map);
	    if(map < sub_map[0]){
		again=TRUE;
		map=sub_map[0];	// possible rounding errors...
	    }
    	   }
    	} Map=map;
	return Map;
}

BooLean	che_typ::RmColButKeepLPR(Int4 j)
// routine that doesn't recalculate LPR.
{
	if(pps->NumColumns() <= MinNumColumns) return FALSE;
	assert(j > 0 && j <= LenSeq(Query));
	pps->RmColumn(j); 
	return TRUE;
}

double	*che_typ::RemoveColumn(Int4 j)
{
	if(pps->NumColumns() <= MinNumColumns) return NULL;
	assert(j > 0 && j <= LenSeq(Query));
#if 1
	if(Verbose){
	   sst_typ  old_sst = pps->RtnSST(j);
	   fprintf(stderr,"Removing col: ");
	   PutSST(stderr,old_sst,AB); fprintf(stderr,"%d\n",j);
	}
#endif
	pps->RmColumn(j); 
	return pps->SubLPR(ResWt[2],ResWt[1]);
}

BooLean	che_typ::SamplePattern(double best_map,Int4 j,sst_typ *sstJ,double temperature)
{
	double	d,map,*submap;
	Int4	s;
	sst_typ	best_sst=0;
	if(pps->NumColumns() >= MaxNumColumns) return FALSE;
	// assert(qsst[j] == 0);
	sst_typ	old_sst = pps->RtnSST(j);
	assert(j > 0 && j <= pps->LenPattern( ));
	Int4 num_pttrns=0;
	for(s=1; sstJ[s]; s++) num_pttrns++;
	BooLean	*Skip; NEW(Skip,num_pttrns + 4, BooLean);
	double	*deltaLPR; NEW(deltaLPR, num_pttrns + 3,double);
	for(s=1; sstJ[s]; s++){
		pps->AddColumn(j,sstJ[s]);
		d=pps->QuickLPRsubJ(j,ResWt[2],ResWt[1]);
		// if(j==2) { SSetOkay(j,sstJ[s],stderr); fprintf(stderr,"    %.3f nats\n",d); }
		if(d <= 0.0){ Skip[s]=TRUE; pps->RmColumn(j); continue; }
		submap = pps->SubLPR(ResWt[2],ResWt[1]);
		map = submap[0];
#if 0
		deltaLPR[s]=map;
		// deltaLPR[s]=submap[j];
		// if(map > best_map)  // Insist on positive sub_map at column position...
		if(submap[j] > 0.0 && map > best_map) // Insist on positive sub_map at column position...
		{ best_map=map; best_sst=sstJ[s]; }
#else
		deltaLPR[s]=submap[j];
		if(submap[j] > best_map) // Insist on positive sub_map at column position...
		{ best_map=submap[j]; best_sst=sstJ[s]; }
#endif
	}
	if(best_sst==0){ 
		if(old_sst == 0) pps->RmColumn(j); 
		else pps->AddColumn(j,old_sst);
		submap = pps->SubLPR(ResWt[2],ResWt[1]); free(Skip); 
    	  	free(deltaLPR); return FALSE; 
	} else {
	  //**************** Sample a pattern ****************
	  if(temperature > 0.0){		// else use best_sst;
	    double sum=0.0;
	    assert(temperature <= 300.0);
	    for(s=1; sstJ[s]; s++){
		if(Skip[s]) continue;
		double d=deltaLPR[s]-best_map; 
		assert(d <= 0.0);
		if(d < -100){ d = 0.0; deltaLPR[s]=0.0; }
		else {
		  d = exp(d);
		  if(temperature == 300.0){ deltaLPR[s]=d; sum += d; }
		  else {
		    double y=(300.0/temperature);
		    assert(!(d==0 && y <= 0)); assert(d >= 0);
		    if(d == 0) deltaLPR[s]=0;
		    else {
			 d = pow(d,y);
			 // assert(isfinite(d)); this is failing!!!
			 if(!isfinite(d)) assert(!"isfinite failed");
			 deltaLPR[s]=d; sum += d;
		    }
		  }
		}
	    }
	    double rand = sum * ((double) Random()/(double) RANDOM_MAX);        // random number...
	    for(sum=0.0, s=1; sstJ[s]; s++){
		 if(Skip[s]) continue;
		 sum += deltaLPR[s];
		 if(sum >= rand){ best_sst=sstJ[s]; best_map=deltaLPR[s]; break; }
	    } assert(sum >= rand);
    	  } free(deltaLPR); 
	  //**************** End of pattern sampling ****************
	  pps->AddColumn(j,best_sst); 
	  submap = pps->SubLPR(ResWt[2],ResWt[1]);
	  if(Verbose){
	   fprintf(stderr,"Redefined pattern set: '");
	   PutSST(stderr,old_sst,AB); fprintf(stderr,"%d' -> '",j);
	   PutSST(stderr,best_sst,AB); fprintf(stderr,"%d' (%.3f)\n",j,submap[j]); 
	  } free(Skip);
	  return TRUE; 
	}
}

sst_typ	che_typ::ForceBestColumn(Int4 j,sst_typ *sstJ)
{
	double	map,*submap,best_map=-9999999999999999999.9;
	sst_typ	best_sst=0;
	// pps->RmColumn(j); 
	sst_typ	old_sst = pps->RtnSST(j);
	for(Int4 s=1; sstJ[s]; s++){
		pps->AddColumn(j,sstJ[s]);
		submap = pps->SubLPR(ResWt[2],ResWt[1]);
		map = submap[j];
		if(map > best_map){
#if 0
fprintf(stderr,"ForceBestColumn map=%.2f (from %.2f)(%d,%d)\n",map,best_map,pps->NumColumns(),NumSignificantColumns( ));
#endif
			best_map=map; best_sst=sstJ[s]; 
		}
	}
	assert(best_sst != 0); 
	if(Verbose){
	   fprintf(stderr,"Redefined pattern set: '");
	   PutSST(stderr,old_sst,AB); fprintf(stderr,"%d' -> '",j);
	   PutSST(stderr,best_sst,AB); fprintf(stderr,"%d'\n",j); 
	} pps->AddColumn(j,best_sst); 
	submap = pps->SubLPR(ResWt[2],ResWt[1]);
	return best_sst; 
}

BooLean	che_typ::AddBestColumn(double best_map,Int4 j,sst_typ *sstJ)
{
	double	map,*submap;
	sst_typ	best_sst=0;
	if(pps->NumColumns() >= MaxNumColumns) return FALSE;
	// assert(qsst[j] == 0);
	sst_typ	old_sst = pps->RtnSST(j);
	for(Int4 s=1; sstJ[s]; s++){
		pps->AddColumn(j,sstJ[s]);
		submap = pps->SubLPR(ResWt[2],ResWt[1]);
		map = submap[0];
		// if(map > best_map)
		if(submap[j] > 0.0 && map > best_map) {
			// Insist on positive sub_map at column position...
#if 0
fprintf(stderr,"AddBestColumn map=%.2f (from %.2f)(%d,%d)\n",map,best_map,pps->NumColumns(),NumSignificantColumns( ));
#endif
			best_map=map; best_sst=sstJ[s]; 
		}
	}
	if(best_sst==0){ 
		if(old_sst == 0) pps->RmColumn(j); 
		else pps->AddColumn(j,old_sst);
		submap = pps->SubLPR(ResWt[2],ResWt[1]);
		return FALSE; 
	} else {
	  if(Verbose){
	   fprintf(stderr,"Redefined pattern set: '");
	   PutSST(stderr,old_sst,AB); fprintf(stderr,"%d' -> '",j);
	   PutSST(stderr,best_sst,AB); fprintf(stderr,"%d'\n",j); 
	  }
	   pps->AddColumn(j,best_sst); 
	   submap = pps->SubLPR(ResWt[2],ResWt[1]);
	   return TRUE; 
	}
}

Int4	che_typ::ForceRmWorstColumn( )
{
	Int4	j,worst_j=0;
	double	worst=0.0,old_sub_map;
	double  *sub_map=SubMap( );
	
	// sub_map=pps->SubLPR(ResWt[2],ResWt[1]);	// same as when sub_map passed in...
	for(j=1; j <= LenSeq(Query); j++){
		if(sub_map[j] < worst){ worst_j = j; worst=sub_map[j]; }
	}
	if(Verbose) fprintf(stderr,"Columns: %d ",pps->NumColumns());
	if(worst < 0.0) pps->RmColumn(worst_j);  
	if(Verbose) fprintf(stderr,"-> %d\n",pps->NumColumns());
	sub_map = pps->SubLPR(ResWt[2],ResWt[1]);
	return worst_j; 
}


BooLean	che_typ::RmWorstNegColumn(BooLean RmAllNeg,double *sub_map)
{
	Int4	j,worst_j=0;
	double	worst=0.0,old_map,old_sub_map;
	BooLean	success=FALSE;
	sst_typ	sstJ;
	
	// sub_map=pps->SubLPR(ResWt[2],ResWt[1]);	// same as when sub_map passed in...
	old_map=sub_map[0];
#if 1	// Fix tendency to go down to zero columns (afn: 12/6/05)
	if(pps->NumColumns() <= MinNumColumns) return FALSE;
#endif
#if 0	// DEBUG
	double *Old_sub_map; NEW(Old_sub_map,LenSeq(Query)+3,double);
	for(j=1; j <= LenSeq(Query); j++) Old_sub_map[j]=sub_map[j];
#endif
	if(RmAllNeg){ // TEST: remove all negative columns at this point.
	  for(j=1; j <= LenSeq(Query); j++){
	   if(qsst[j] && sub_map[j] <= 0.0){
	     sstJ=pps->RtnSST(j);
	     if(Verbose) fprintf(stderr,"Columns: %d ",pps->NumColumns());
	     old_sub_map=sub_map[j];
	     pps->RmColumn(j);  
	     sub_map = pps->SubLPR(ResWt[2],ResWt[1]);
	     if(sub_map[0] < old_map){	// reject all of these regardless
		// alpha gets reset with RmColumn() which changes subLPR at each position.
#if 0		// DEBUG...
	        if(Verbose) fprintf(stderr,
			"(j=%d new map=%.2f; old_map=%.2f; submap=%.2f; full_map=%.2f; alpha=%.2f)\n",
			j,sub_map[0],old_map,old_sub_map,pps->LPR(ResWt[2],ResWt[1]),pps->Alpha());
		double map=sub_map[0];
		for(Int4 jj=1; jj <= LenSeq(Query); jj++){
			if(Old_sub_map[jj] != sub_map[jj]){
				fprintf(stderr,"%d: Old_sub_map=%.2f; new_sub_map=%.2f\n",
					jj,Old_sub_map[jj],sub_map[jj]);
			}
			map-=sub_map[jj];
		}
		fprintf(stderr,"\n**** map = %.2f; map - Sum(sub_map[j]) = %.2f; alpha=%.2f ****\n",
			sub_map[0],map,pps->Alpha());
#endif
		
		pps->AddColumn(j,sstJ);
	        sub_map = pps->SubLPR(ResWt[2],ResWt[1]);
#if 0
	        if(Verbose) fprintf(stderr,"\n(j=%d new_old_map=%.2f; new old submap=%.2f; full_map=%.2f)",
			j,sub_map[0],sub_map[j],pps->LPR(ResWt[2],ResWt[1]));
#endif
	     } else {
// fprintf(stderr,"RmWorstColumn map=%.2f (from %.2f)\n",sub_map[0],old_map);
		old_map = sub_map[0];
	     	success=TRUE;
	     }
#if 1	// Fix tendency to go down to zero columns (afn: 12/6/05)
	     if(pps->NumColumns() <= MinNumColumns) return success;
#endif
	     if(Verbose) fprintf(stderr,"-> %d\n",pps->NumColumns());
	   }
	  } return success;
	} else {
	 worst = sub_map[1]; worst_j=1;
	 for(j=1; j <= LenSeq(Query); j++){
		if(sub_map[j] < worst){ worst_j = j; worst=sub_map[j]; }
	 }
	 if(worst_j > 0){
	   if(Verbose) fprintf(stderr,"Columns: %d ",pps->NumColumns());
	   pps->RmColumn(worst_j);  
	   if(Verbose) fprintf(stderr,"-> %d\n",pps->NumColumns());
	   sub_map = pps->SubLPR(ResWt[2],ResWt[1]);
	   return TRUE; 
	 } else return FALSE;
	}
}

double	*che_typ::SampleSeqs(double map)
{
	double	*submap,map0=map;

	if(FIXED_PARTITIONS){ return pps->SubLPR(ResWt[2],ResWt[1]); }
	do {
	    map=map0;
	    submap=SampleSeqsOnce(map);
	    map0=submap[0];
	} while(map0 > map);
	return submap;
}

double	*che_typ::SampleSeqsOnce(double map0)
{
	double	map,last_map;
	Int4	sq,num_sampled;

        last_map=map0;
// fprintf(stderr,"DEBUG S\n");
	if(FIXED_PARTITIONS){ return pps->SubLPR(ResWt[2],ResWt[1]); }
	Int4	mv,max_move = num_sets-1; 
	for(sq=1; sq <= num_seqs; sq++){

	   if(MemberSet(sq,RmSet)) continue; // Ignore removed sequences (afn: 11/12/09).

	   insrtHeap(sq,(keytyp)Random(),sqH);
	}
	num_sampled=0;
// fprintf(stderr,"DEBUG S1\n");
	while((sq=delminHeap(sqH)) != 0){
	 if(MemberSet(sq,gold_set)) continue;
	 else if(MemberSet(sq,RmSet)) continue; // This shouldn't happen (afn: 11/12/09).
	 else num_sampled++;

	 // if(num_sampled > 10) break;

// fprintf(stderr,"DEBUG S2\n");
	 for(mv=1; mv <= max_move; mv++) {
	   Int4 mv0;
	   if(random_integer(2)==0) mv0=mv; else mv0=-mv;
	   if(!SwapSeq(sq,mv0)){
		mv0=-mv0;
	    	if(!SwapSeq(sq,mv0)) break; // impossible to make move this large.
	   }
	   map=pps->LPR(ResWt[2],ResWt[1]);	// Otherwise swap worked...
#if 0
	  unsigned char *seq=GetAlnResInSiteCMSA(1,sq,mcma);
	  // double probRatio = pps->ProbRatio(seq,SqWt,ResWt[2],ResWt[1]);
  if(map > 0.0){
	  if(map > map0 && probRatio > 1.0 ) {
	  	fprintf(stderr,"Improved map = %g vs %g; ratio = %g\n",map,map0,probRatio);
	  } else if(map <= map0 && probRatio <= 1.0){
	  	fprintf(stderr,"worse map = %g vs %g; ratio = %g\n",map,map0,probRatio);
	  } else 
	  	fprintf(stderr,"Conflict!: map = %g vs %g; ratio = %g\n",map,map0,probRatio);
  }
#endif
// fprintf(stderr,"map = %g\n",map);
// fprintf(stderr,"DEBUG S3\n");
	   if(map < map0) SwapSeq(sq,-mv0); // Switch back if map not improved.
	   else { 
// fprintf(stderr,"SampleSeqsOnce map=%.2f (from %.2f)\n",map,map0);
		if(Verbose) fprintf(stderr,"LPR= %g (",map);
		for(Int4 s=1; s <= num_sets; s++){
			map_set[s]=info_set[s];
			if(Verbose) fprintf(stderr,"%d ",card_set[s]);
		}
		if(Verbose) fprintf(stderr," seqs; move = %d; seq = %d; alpha = %g; cols: %d(%d))\n",
			mv0,sq,pps->Alpha(),pps->NumColumns(),NumSignificantColumns( ));
		// PutSqDefLineCMSA(stderr,sq,mcma);
		Map=map0=map;
		// while(delminHeap(mvH) != 0) ;
		break;
	   }
	 }
	} return pps->SubLPR(ResWt[2],ResWt[1]);
}

//********************************* mcBPPS ****************************************//
// Added for multiple category BPPS sampler (afn: 11/12/09).

#if 0
BooLean che_typ::RestoreSeq(Int4 sq)
// put arbitrarily into background or foreground.
{
        Int4    s,set;
	UInt4	wt;

	if(!MemberSet(sq,RmSet)) return FALSE; // already removed.
        for(set =1 ; set <= num_sets; set++){
          if(MemberSet(sq,seq_set[set])) print_error("RestoreSeq( ) input error");;
          // if(MemberSet(sq,seq_set[set])) return FALSE;
	}
	set = 1 + random_integer(num_sets);
        DeleteSet(sq,RmSet); AddSet(sq,seq_set[set]);
	unsigned char *seq=GetAlnResInSiteCMSA(1,sq,mcma);
	for(s=1; s <= LenSeq(Query); s++){
		if(GlobalSqWts) wt=AveSqIWt[sq]; else wt = SqWt[s][sq];
		ResWt[set][s][seq[s]] += wt;
		TotalResWt[set][s] += wt;	// added 10-31-2013; afn
	} 
	card_set[set]++; Card_Set[set]+=AveSqIWt[sq];
	return TRUE;
}
#endif

BooLean che_typ::RestoreSeq(char NewPartition, Int4 sq)
// put arbitrarily into background or foreground.
// used by cmcBPPS procedure.
{
        Int4    s,nset,set;
	UInt4	wt;

	if(!MemberSet(sq,RmSet)) return FALSE; // sequence was not removed.
        for(set=1 ; set <= num_sets; set++){
          if(MemberSet(sq,seq_set[set])) print_error("RestoreSeq( ) error: already assigned to a partition!");
          // if(MemberSet(sq,seq_set[set])) return FALSE;
	}

	if(NewPartition == '-'){ nset=2; } else if(NewPartition == '+'){ nset=1; }	
	else print_error("RestoreSeq( ) input error");

        DeleteSet(sq,RmSet); AddSet(sq,seq_set[nset]);
	unsigned char *seq=GetAlnResInSiteCMSA(1,sq,mcma);
	for(s=1; s <= LenSeq(Query); s++){
		if(GlobalSqWts) wt=AveSqIWt[sq]; else wt = SqWt[s][sq];
		ResWt[nset][s][seq[s]] += wt;
		TotalResWt[nset][s] += wt;	// added 10-31-2013; afn
	} card_set[nset]++; Card_Set[nset]+=AveSqIWt[sq];
	return TRUE;
}

BooLean che_typ::RemoveSeq(Int4 sq)
{
        Int4    s,set;
	UInt4	wt;

	if(MemberSet(sq,RmSet)) return FALSE; // already removed.
        for(set =1 ; set <= num_sets; set++){
          if(MemberSet(sq,seq_set[set])){
		unsigned char *seq=GetAlnResInSiteCMSA(1,sq,mcma);
		for(s=1; s <= LenSeq(Query); s++){
			if(GlobalSqWts) wt=AveSqIWt[sq]; else wt = SqWt[s][sq];
			ResWt[set][s][seq[s]] -= wt;
			TotalResWt[set][s] -= wt;	// added 10-31-2013; afn
		} 
                DeleteSet(sq,seq_set[set]); AddSet(sq,RmSet);
		card_set[set]--; Card_Set[set] -= AveSqIWt[sq];	
                return TRUE;	// check to see if in both sets?!?
	  }
        } assert(!"RemoveSeq( ) fatal error");
	return FALSE;
}


char	che_typ::RtnPartition(Int4 sq)
{
	if(MemberFG(sq)) return '+'; if(MemberBG(sq)) return '-'; if(RemovedSeq(sq)) return 'o'; 
	print_error("che_typ::RtnPartition(): this should not happen"); 
}

char	che_typ::SetPartition(char target_p, Int4 sq)
// mv sequence from current partition to Target partition
// return previous partition: FG = '+'; BG = '-'; RM = 'o'.
// For use by cmcBPPS procedure.
{
        Int4    s,oset,nset;
	UInt4	wt;
	char	rtn;

	switch (target_p){
	   case '+':	// move to FG partition.
	     {
		if(MemberFG(sq)) return '+';  // already in FG;
		else if(MemberBG(sq)){ rtn='-'; oset=2; nset=1; }
		else if(RemovedSeq(sq)){ assert(RestoreSeq('+', sq)); return 'o'; } 
		else print_error("SetPartition( ) input error 1");
	     } break;
	   case '-':	// move to BG partition.
	     {
		if(MemberBG(sq)) return '-';  // already in BG;
		else if(MemberFG(sq)){ rtn='+'; oset=1; nset=2; }
		else if(RemovedSeq(sq)){ assert(RestoreSeq('-', sq)); return 'o'; }
		else print_error("SetPartition( ) input error 2");
	     } break;
	   case 'o':	// Remove from alignment.
	     {
		if(RemovedSeq(sq)) return 'o'; 
		else {
		   if(MemberBG(sq)) rtn = '-'; else if(MemberFG(sq)) rtn = '+'; 
		   else print_error("SetPartition( ) input error 0");
		   RemoveSeq(sq);  return rtn;
		}
	     } break;
	   default: print_error("SetPartition( ) input error 3");
	}
        AddSet(sq,seq_set[nset]);
	card_set[nset]++;
	Card_Set[nset]+=AveSqIWt[sq];
	unsigned char *seq=GetAlnResInSiteCMSA(1,sq,mcma);
	for(s=1; s <= LenSeq(Query); s++){
		if(GlobalSqWts) wt=AveSqIWt[sq]; else wt = SqWt[s][sq];
		ResWt[oset][s][seq[s]] -= wt;
		TotalResWt[oset][s] -= wt;	// added 10-31-2013; afn
		// SumWts[oset][s] -= wt;
		ResWt[nset][s][seq[s]] += wt;
		TotalResWt[nset][s] += wt;	// added 10-31-2013; afn
		// SumWts[nset][s] += wt;
	} 
        DeleteSet(sq,seq_set[oset]);
	card_set[oset]--; Card_Set[oset] -= AveSqIWt[sq];	// rounding errors?
	return rtn;
}

char	che_typ::ChngPartition(Int4 sq)
// mv sequence from current partition to alternative partition
// return new partition: FG = '+'; BG = '-'; RM = 'o'.
// Return 0 if sequence has been removed from Main set.
// For use by cmcBPPS procedure.
{
        Int4    s,oset,nset;
	UInt4	wt;
	char	rtn;

	if(MemberSet(sq,RmSet)) return 'o'; // already removed.
	if(MemberFG(sq)){ rtn='-'; oset=1; nset=2; }
	else if(MemberBG(sq)){ rtn='+'; oset=2; nset=1; }	
	else print_error("MoveSeq( ) input error");
        AddSet(sq,seq_set[nset]);
	card_set[nset]++;
	Card_Set[nset]+=AveSqIWt[sq];
	unsigned char *seq=GetAlnResInSiteCMSA(1,sq,mcma);
	for(s=1; s <= LenSeq(Query); s++){
		if(GlobalSqWts) wt=AveSqIWt[sq]; else wt = SqWt[s][sq];
		ResWt[oset][s][seq[s]] -= wt;
		TotalResWt[oset][s] -= wt;	// added 10-31-2013; afn
		// SumWts[oset][s] -= wt;
		ResWt[nset][s][seq[s]] += wt;
		TotalResWt[nset][s] += wt;	// added 10-31-2013; afn
		// SumWts[nset][s] += wt;
	} 
        DeleteSet(sq,seq_set[oset]);
	card_set[oset]--;
	Card_Set[oset] -= AveSqIWt[sq];	// rounding errors?
	return rtn;
}

//********************************* mcBPPS ****************************************//

BooLean che_typ::SwapSeq(Int4 sq, Int4 up)
#if 0   //*****************************************************************
        Moves a sequence up or down one set in the hierarchy...
#endif  //*****************************************************************
{
        Int4    s,oset,nset;
	UInt4 omax,nmax,wt;

	if(MemberSet(sq,RmSet)) return FALSE; // don't allow swapping of these sequences.

        if(up==0) return FALSE;
        for(oset =1 ; oset <= num_sets; oset++){
          if(MemberSet(sq,seq_set[oset])){
                nset = oset + up;
                if(nset > num_sets || nset < 1) return FALSE;
                AddSet(sq,seq_set[nset]);
		card_set[nset]++;
		Card_Set[nset]+=AveSqIWt[sq];
		unsigned char *seq=GetAlnResInSiteCMSA(1,sq,mcma);
		omax=nmax=0;
		for(s=1; s <= LenSeq(Query); s++){
			if(GlobalSqWts) wt=AveSqIWt[sq]; else wt = SqWt[s][sq];
			ResWt[oset][s][seq[s]] -= wt;
			TotalResWt[oset][s] -= wt;	// added 10-31-2013; afn
			// SumWts[oset][s] -= wt;
			ResWt[nset][s][seq[s]] += wt;
			TotalResWt[nset][s] += wt;	// added 10-31-2013; afn
			// SumWts[nset][s] += wt;
		} 
                DeleteSet(sq,seq_set[oset]);
		card_set[oset]--;
		Card_Set[oset] -= AveSqIWt[sq];	// rounding errors?
                return TRUE;
          }
        } assert(!"SwapSeq( ) fatal error");
        return FALSE;
}

Int4    *che_typ::PutBestPttrn(Int4 &N)
{
        Int4 *fake2real;

	fprintf(stderr,"DEBUG: N=%d\n",N);
        NEW(fake2real,LenSeq(Query)+4,Int4);
        for(Int4 s=1; s <= LenSeq(Query); s++){
                fake2real[s]=FakeToRealCMA(1,s,qcma);
        }
        Int4 *Rtn=pps->nBestPttrn(N,fake2real,ResWt[2],ResWt[1]);
        free(fake2real);
        return Rtn;
}

