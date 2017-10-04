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

#if !defined(_CHE_TYP_)
#define _CHE_TYP_

#include "cmsa.h"
#include "dheap.h"
#include "random.h"
#include "chn_typ.h"
#include "set_typ.h"
#include "bpps_typ.h"
#include "pps_typ.h"
#include "residues.h"
#include "histogram.h"

#if 0	//**********************************************************************
	Optimally partitions the aligned input sequences into telescoping
	sets that are hierarchically arranged from the most specifically
	conserved patterns to the least conserved patterns.

	Two types of residues: 
	  '1' == identical to query at corresponding position (and not 'X').
	  '0' == different from the query at corresponding position.


	Think about more efficient (bitwise?) representation and routines.
#endif	//**********************************************************************

class che_typ {         // Comparative Heirarchical Alignment type
public:
                che_typ( ){ assert(!"Illegal constructor"); }
		che_typ(char *sst_str,sst_typ **sst_in,chn_typ *chain,double,double,
			BooLean verbose,char mode,double **Rho,double PriorRi, 
			BooLean UseGlobalSqWts);

		// following constructor is used by the cmcBPPS sampler
		che_typ(char *sst_str, chn_typ *chain,double,double,
			BooLean verbose,char mode,double **Rho,double PriorRi, 
			BooLean UseGlobalSqWts);

		// New constructor for omcBPPS sampler.
		che_typ(char *sst_str, chn_typ *chain, BooLean verbose,char mode,
				char Type, BooLean UseGlobalSqWts,double PriorRi=1.0);

                ~che_typ( ){ Free( ); }
	e_type	KeyE( ){ return Query; }
	//********************* omc_resume routines ***********************
	double  FakeSubLLR(UInt4 wt, char mode, Int4 sq, cma_typ xcma);
	double  FakeSubLLR(UInt4 wt, char mode, unsigned char *seq, Int4 Len);
	//********************* omc_resume routines ***********************
#if 1	// Evolving CSQ.
	e_type	EmitRandomSeq(){	// Update the csq used to restrain the sst.
		  Int4	s,r,N,sim_r,Sum,Rnd;
		  e_type E=CopySeq(Query);
		  for(s=1; s <= LenSeq(Query); s++){
		     N = TotalResWt[1][s]; assert(N < 900000000);  // less than 900 million!
		     Rnd=random_integer(N)+1; sim_r = -1;
	// fprintf(stderr,"%d: N = %d;",s,N);
		     for(Sum=0,r=0; r <= nAlpha(AB); r++){
			Sum += ResWt[1][s][r];		// 1 == Foreground.
	// fprintf(stderr,"  %c = %d sum = %d; Rnd=%d\n",AlphaChar(r,AB),ResWt[1][s][r],Sum,Rnd);
			if(Sum >= Rnd){ sim_r=r; break; }
		     } assert(sim_r >= 0);
		     EqSeq(s,sim_r,E); 
		  } return E;
		}
	BooLean	UpdateConsensus(){	// Update the csq used to restrain the sst.
		  Int4	s,r,max_r=0;	// WARNING: this can get pattern and csq out of syn in mcs_typ!!!
		  for(s=1; s <= LenSeq(Query); s++){
		     UInt4  X,MaxCnt=0;
		     max_r=ResSeq(s,Query);
		     for(r=1; r <= nAlpha(AB); r++){
			X = ResWt[1][s][r];	
			if(X > MaxCnt){ max_r=r; MaxCnt=X; }
		     } EqSeq(s,max_r,Query); 
		     if(XedSeq(Query)) EqXSeq(s,max_r,Query); // 
		  }
		}
	void	SetConsensus(e_type csqE){	// Set the csq used to restrain the sst.
		  Int4	s,r;
		  assert(LenSeq(Query) == LenSeq(csqE));
		  for(s=1; s <= LenSeq(Query); s++){
		     r=ResSeq(s,csqE); EqSeq(s,r,Query); 
		     if(XedSeq(Query)) EqXSeq(s,r,Query); // 
		  }
		}
#endif
	void	PutAll(Int4,char,double,double);
	void	PutParametersBPPS(FILE *fp){ pps->PutParameters(fp); }

	void	PutFgBgSeqs(Int4 id);
	double	PutBoltzmannLike(FILE *fptr,double *frq, Int4 c, Int4 n);  // in omc_debug.cc
	void    PutCDTreePttrnFile(FILE *fptr,char *id){
			pps->PutCDTreeSubLPR(fptr,id,ResWt[2],ResWt[1]);
		}
	void	PutPttrnFile(FILE *fp,Int4 id);

	void	PutAll(Int4 id,char SFBG,double min_nats)
		{ PutAll(id,SFBG,min_nats,0.5); }
	void    WriteSubLPR(FILE *fptr);
	// void	Put( ) { Put(FALSE); }
	// void	Put(BooLean);
	void    PutChn(const char *str,char SFBG,double min_nats){
			PutChn(str,SFBG,min_nats,0.5);
		}
	void    PutChn(const char *str,char SFBG,double min_nats,double min_freq);
	void	Put(FILE *,char ,double,double);
	void	Put(FILE *fp,char SFBG,double min_nats)
		{ Put(fp,SFBG,min_nats,0.5); }
//********************* New routines for mpps ****************************
        void    PutDS(FILE *fp,double min_nats,double min_freq);
        void    PutDSubSets(FILE *fp,double min_nats,double min_freq);
        void    PutFG(FILE *fp,double min_nats,double min_freq);
        void    PutBG(FILE *fp,double min_nats,double min_freq);
        void    PutMainSet(FILE *fp,double min_nats,double min_freq);
        bpps_typ *BPPS() { return pps; }
        rst_typ *BPPS_RST() { return pps->RtnRST(); }
	double	*RemoveColumn(Int4 j); // 
	BooLean RmColButKeepLPR(Int4 j);	// removes column but doesn't recalculate LPR.
	double	*SubMap( ){ double *d=pps->SubLPR(ResWt[2],ResWt[1]); Map=d[0]; return d; }
	double	CalcLLR( ){ double *d=pps->SubLPR(ResWt[2],ResWt[1]); Map=d[0]; return Map; }
	double	CalcRawLLR( ){ double *d=pps->SubRawLPR(ResWt[2],ResWt[1]); return d[0]; }
	double	RawLLR( ){ double d=pps->RawLPR(ResWt[2],ResWt[1]); return d; }
	void	UpdateParametersBPPS( ){ pps->UpdateParameters(ResWt[2],ResWt[1]); }
	BooLean AddBestColumn(double best_map,Int4 j,sst_typ *sstJ);
	sst_typ	ForceBestColumn(Int4 j,sst_typ *sstJ);
	BooLean SamplePattern(double best_map,Int4 j,sst_typ *sstJ,double temperature);
	BooLean AddColumn(Int4 j,sst_typ sstJ){
			if(pps->NumColumns() >= MaxNumColumns) return FALSE;
			else pps->AddColumn(j,sstJ); return TRUE;
		}
	BooLean	PttrnPos(Int4 j){ 
		   if(pps->RtnSST(j) == 0) return FALSE;
		   else return TRUE;
		}
//********************* End new routines for mpps ****************************
	void	PutSeqs( );
	void	PutPattern(FILE *fp);
	void    PutBestPatterns(FILE *fptr,BooLean UseRealPos);
	void	PutSubLPRs(char *);
	void	PutSubLPRs(FILE *fp);
	void	PutPattern(FILE *fp,Int4 *fake2real);
	a_type	Alphabet( ){ return AB; }
	void	PutVSI_Info(FILE *fp,Int4 *fake2real);
	void	PutInfo( );
	void	PutInfo(FILE *);
	double	PutInfoShort(FILE *);
	double	RtnNatsPerWtSq( );	// Assumes LPR calculated and using AveSqWts
	double	RtnAveWtSeqs();
	void	PutInfo(FILE *,Int4);
	char	*FileName( ){ return chn->FileName( ); }
		// Sample columns, residue sets and sequences to optimize map.
	double	Sample(BooLean RmCols,sst_typ **res_sst){
			return Sample(RmCols,FALSE, res_sst);
		}
	double	Sample(BooLean RmCols,BooLean Lock, sst_typ **res_sst);
	UInt4	NumColumns( ){ return pps->NumColumns( ); }
	Int4	LengthPattern( ){ return pps->LenPattern( ); }
	UInt4	NumSignificantColumns( ){
			UInt4	C=0,j;
			double *sub_map=pps->SubLPR(ResWt[2],ResWt[1]);
			for(j=1; j <= LenSeq(Query); j++){
			  if(qsst[j] && sub_map[j] > 1.0) C++;
			} return C; 
		}
	double 	ProbCols( ){ return pps->ProbCols( ); }
#if 1	// return true if both che share the same partitions.
	BooLean TheSameSets(che_typ *che){
		  if(num_sets != che->num_sets) return FALSE;
		  for(Int4 s=1; s <= num_sets; s++){
		     if(CardInterSetINotJ(seq_set[s],che->seq_set[s]) != 0) return FALSE;
		  } return TRUE;
		}
	BooLean TheSameWeights(che_typ *che){
		   Int4 n=0; assert(this->mcma == che->mcma);
		   for(Int4 i=1; i <= NumSeqsCMSA(mcma); i++){
			if(this->AveSqIWt[i] != che->AveSqIWt[i]){
			   fprintf(stderr,"%d: AveSqIWt1 = %d; AveSqIWt2 = %d\n",
				i,this->AveSqIWt[i],che->AveSqIWt[i]); n++;
			} if(n > 10) return FALSE;
		   } if(n > 0) return FALSE; else return TRUE;
		}
	BooLean TheSamePattern(che_typ *che){
			if(NumColumns( ) != che->NumColumns( )) return FALSE;
			return pps->TheSamePattern(che->pps); 
		}
	BooLean TheSameCardSet(che_typ *che){
		   assert(num_sets == che->num_sets);
		   for(Int4 s=1; s <= num_sets; s++){
			if(this->card_set[s] != che->card_set[s]){
			   fprintf(stderr,"%d: card_set1=%d; card_set2=%d\n",
					s,this->card_set[s],che->card_set[s]);
			   return FALSE;
			}
			if(this->Card_Set[s] != che->Card_Set[s]){
			  Int4 total1=0,total2=0;
		   	  for(Int4 x=1; x <= num_sets; x++){
			   fprintf(stderr,"%d: Card_Set1=%d; Card_Set2=%d\n",
					x,this->Card_Set[x],che->Card_Set[x]);
			   total1 += this->Card_Set[x]; total2 += che->Card_Set[x];
			  } fprintf(stderr,"total: Card_Set1: %d; Card_Set2: %d\n", total1,total2);
			  return FALSE;
			}
		   } return TRUE;
		}
	void	PutCardSets(FILE *fp){
			Int4 total=0;
		   	for(Int4 x=1; x <= num_sets; x++){
			   fprintf(fp,"%d: Card_Set=%d\n",x,this->Card_Set[x]);
			   total += this->Card_Set[x]; 
			} fprintf(fp,"total: %d\n", total);
		}
#endif
	BooLean TheSame(che_typ *che){
			if(!pps->TheSamePattern(che->pps)) return FALSE;
			if(num_sets != che->num_sets) return FALSE;
			for(Int4 s=1; s <= num_sets; s++){
			  if(CardInterSetINotJ(seq_set[s],
				che->seq_set[s]) != 0) return FALSE;
			  if(this->card_set[s] != che->card_set[s]) return FALSE;
			}
			if(!this->TheSameCardSet(che)) return FALSE;
			return TRUE;
		}
	void	Forbid(sst_typ **sst){ 
		   sst_typ *Qst=pps->RtnSST(); 
		   for(Int4 i=1; i<= LenSeq(Query); i++){
			// sst[i][0] = Qst[i]; 
			if(!Qst[i]) for(Int4 j=0; sst[i][j]; j++) sst[i][j]=0;
#if 1			// allow only inputed sets...
			if(Qst[i] != 0){ sst[i][0]=Qst[i]; sst[i][1]=0; }
#endif
		   }
		}
	void	Forbid(char *str,sst_typ **sst){ 
		   Int4 *values,len; 
		   e_type qE=chn->TrueFirstSeq();
		   len = LenSeq(qE) + OffSetSeq(qE);
		   // len = strlen(str); old...core dumps when too many positions...
		   NEW(values,len+2,Int4);
		   Int4 n=ParseIntegers(str, values, "che_typ::Forbid: input error");
		   if(n > len) print_error("input values for -I option are out of range");
		   sst_typ *Qst=pps->RtnSST(); 
		   for(Int4 i=1; i<= n; i++){
			Int4 j=values[i] - OffSetSeq(Query);
		   	// WARNING: this could create a problem if the Query has gaps...
			if(j > 0 && j <= LenSeq(Query)){
			    if(!Qst[j]) for(Int4 k=0; sst[j][k]; k++) sst[j][k]=0;
			}
		   } free(values);
		}
	void    ContribSeqMAP(Int4 ,char ,double *,double ,double );
	void    SetContrast(Int4 N){ assert(N >= 1); Contrast=N; }
	void	FixedPartition( ) { FIXED_PARTITIONS=TRUE; }
	void	BeQuiet( ) { Verbose=FALSE; }
	void	SpeakUp( ) { Verbose=TRUE; }
	void	SetMinNumColumns(Int4 mc){ MinNumColumns = mc; }
	void	SetMaxNumColumns(Int4 mc){ MaxNumColumns = mc; }
	Int4	RtnMinNumColumns(){ return MinNumColumns; }
	Int4	RtnMaxNumColumns(){ return MaxNumColumns; }
	Int4    *PutBestPttrn(Int4 &N);
	double	*SampleSeqs(double);	// Sample sequences to optimize map.
	double  *SampleSeqsOnce(double);
#if 1	// Added for cmcBPPS sampler (afn: 11/24/09).
	// double	*SampleSequence(Int4 sq,double map0);	// Sample one sequence.
	char    ChngPartition(Int4 sq);
	char    RtnPartition(Int4 sq);
	char    SetPartition(char target_p, Int4 sq);
	BooLean	MemberGold(Int4 sq){ if(MemberSet(sq,gold_set)) return TRUE; else return FALSE; }
	BooLean	AddToGold(Int4 sq){
		if(MemberSet(sq,gold_set)) return FALSE;
		AddSet(sq,gold_set); return TRUE; 
	}
#endif
#if 1	// Added for Recursively nested multiple constraint sampler (afn: 11/12/09).
	BooLean	IsPatternPos(Int4 j){ 
			if(pps->RtnSST(j) == 0) return FALSE; else return TRUE;
		}
	// BooLean	RestoreSeq(Int4 sq);
	BooLean	RestoreSeq(char new_partition,Int4 sq);
	BooLean RemoveSeq(Int4 sq);
	BooLean	MemberBG(Int4 sq){ if(MemberSet(sq,seq_set[2])) return TRUE; else return FALSE; }
	BooLean	MemberFG(Int4 sq){ if(MemberSet(sq,seq_set[1])) return TRUE; else return FALSE; }
	BooLean	MemberNeitherFG_BG(Int4 sq){ return !(MemberBG(sq) || MemberFG(sq)); }
	BooLean	MemberEitherFG_BG(Int4 sq){ return (MemberBG(sq) || MemberFG(sq)); }
	BooLean	RemovedSeq(Int4 sq){ if(MemberSet(sq,RmSet)) return TRUE; else return FALSE; }
	set_typ	CreateBG_Set( ){
        		set_typ Set=MakeSet(num_seqs+1); ClearSet(Set);
        		UnionSet(Set,seq_set[2]);
        		// PutSelectCMSA(fp,0,"BackGround",0,Set,mcma);
        		return Set;
		}
	UInt4	RtnCardBG( ){ return card_set[2]; }
	UInt4	RtnCardFG( ){ return card_set[1]; }
	set_typ	RtnBG_Set( ){ return seq_set[2]; }
	set_typ	RtnFG_Set( ){ return seq_set[1]; }
	set_typ	RtnRmSet( ){ return RmSet; }
	set_typ	RtnGoldSet( ){ return gold_set; }
	Int4	UpdateMainSet(che_typ *che2){
#if 1
			print_error("UpdateMainSet() routine turned off");
#else
			set_typ	BG=che2->RtnBG_Set( );
			set_typ	FG=che2->RtnFG_Set( );;
			set_typ	Rm=che2->RtnRmSet( );;
			Int4	n=0,sq;

			if(SetN(FG) != SetN(seq_set[1])) print_error("Set inconsistency in che_typ");
			if(SetN(BG) != SetN(seq_set[2])) print_error("Set inconsistency in che_typ");
			if(SetN(Rm) != SetN(RmSet)) print_error("Set inconsistency in che_typ");

			for(sq=1; sq <= num_seqs; sq++){
				if(MemberSet(sq,BG)) if(RemoveSeq(sq)) n--;
				if(MemberSet(sq,Rm)) if(RemoveSeq(sq)) n--;
			}

			for(sq=1; sq <= num_seqs; sq++){
				if(MemberSet(sq,FG)) if(RestoreSeq(sq)) n++;
			}
			return n; 
#endif
		}
	cma_typ	RtnMainSet() { return mcma; }
	cma_typ	RtnQuerySet() { return qcma; }
	BooLean	RemoveWorstColumn( ){ double  *sub_map=SubMap( ); return RmWorstNegColumn(FALSE,sub_map); } 
	Int4    ForceRmWorstColumn( );
	BooLean	RmWorstNegColumns(BooLean RmAllNeg, sst_typ **res_sst){
		   double  map=-9999.0,*sub_map;
		   BooLean rtn=FALSE;
		   // sub_map=SampleSeqsOnce(map); map=sub_map[0];
		   sub_map=SubMap( ); map=sub_map[0];
		   while(RmWorstNegColumn(RmAllNeg,sub_map)){
		        rtn=TRUE;
			// sub_map=SampleSeqsOnce(map); map=sub_map[0];
		        sub_map=SubMap( ); map=sub_map[0];
		   } return rtn;
		}
#endif
	UInt4	GetResWtFG(Int4 s,Int4 r) {
			assert(s > 0 && s <= LenSeq(Query) && r >= 0 && r <= nAlpha(AB));
			return ResWt[1][s][r]; 
		}
	UInt4	GetResWtBG(Int4 s,Int4 r) { 
			assert(s > 0 && s <= LenSeq(Query) && r >= 0 && r <= nAlpha(AB));
			return ResWt[2][s][r]; 
		}
	UInt4	*GetTotalResWtFG(){ return TotalResWt[1]; }
	UInt4	*GetTotalResWtBG(){ return TotalResWt[2]; }
	UInt4   GetWtFactor(){ return WtFactor; }
	//===== New routines for tri-partitioned model. =====
	UInt4	**GetResWtsFG() { return ResWt[1]; }
	UInt4	**GetResWtsBG() { return ResWt[2]; }
	// ===== return weighted counts with background counts for deletions.
	UInt4	*FGResWtsBILD(Int4 s, UInt4 &totWt);
	UInt4	*BGResWtsBILD(Int4 s, UInt4 &totWt);
	//===== End routines for tri-partitioned model. =====
	Int4	Contrast;	// == number of positions to show (default= 1/10 positions.);
	void	PutTailRTF(FILE *fp){ chn->PutHierarchicalAlnTail(fp); }
	UInt4	RtnAveSqIWt(Int4 sq){ 
			assert(sq > 0 && sq <= num_seqs);
			return AveSqIWt[sq];
		}
#if 0
	double	RtnAveSqIWt(Int4 sq){ 
			assert(sq > 0 && sq <= num_seqs);
			return (double) ( (double) AveSqIWt[sq]/ (double) WtFactor);
		}
#endif
private:
	//**************** che_junk.cc *********************
	BooLean SampleAnyPattern(Int4 j,sst_typ **sst, BooLean **skip, unsigned char *NumResSets,
                        double prior_rho, double temperature);
	BooLean ResetColumnJ(Int4 j, Int4 r, double prior_rho, Int4 s, sst_typ *sstJ);
	BooLean SSetOkay(Int4 s, sst_typ sst,FILE *fp=0);
	//**************** che_junk.cc *********************
	BooLean	Verbose,FIXED_PARTITIONS;
	set_typ	RmSet; // Added for Recursively nested multiple constraint sampler (afn: 11/12/09).
	void	Init(char *,sst_typ **,double,double, BooLean verbose, char mode,double **Rho,
			double PriorRi,char Type, BooLean UseGlobalSqWts);
	BooLean	RmWorstNegColumn(BooLean, double *); // 
	Int4	MinNumColumns;	// don't let columns go below this point.
	Int4	MaxNumColumns;	// don't let columns go above this point.
	void	GetIntegerWts(){ AveSqIWt=chn->GetIntWeights(1,&SqWt); }
	BooLean	SwapSeq(Int4 ,Int4);
	void	Free();		// free memory.
// NOTE: ADD SET ZERO TO REPRESENT QUERY SUBFAMILY!!???
// Sample sequences into zero set and these will be eliminated from 
// other alignments???

	double	Sigma;		// stringency of KL divergence...
	bpps_typ *pps;		// binary partition-pattern selection type.
	// pps_typ *pps;		// partition-pattern selection type.
	sst_typ	*qsst;		// query-based residue sets.
	chn_typ	*chn;
	BooLean	OwnChn;
	double	NullMap;	// For Jun's method.
	a_type	AB;		// alphabet.
	e_type	Query;		// query sequence
	cma_typ	*IN_CMA;
	Int4	NumAnalysis;
	cma_typ	qcma;		// query cma file
	cma_typ	mcma;		// main cma file
	unsigned short num_sets;
	Int4	num_seqs;	// total number of sequences.
	set_typ	seq_set[9];
	double	info_set[9];	// information contributed by each set.
	set_typ	gold_set;	// keep these sequences in the query set!
	dh_type	posH,sqH,mvH;	// heaps for random ordering of data.
	double	Map,map_set[9];
	UInt4   WtFactor;
	UInt4	***ResWt;  // ResWt[set][s][r];  // 100 counts for residue weight.
	UInt4	**TotalResWt;	// TotalResWt[set][s]; sums over all residue types.
	// UInt4 *SumWts[9];  // SumWts[set][s] sum of weights at position s;
	unsigned char **SqWt;	   // SqWt[s][sq];  // integerize sequence weights
	Int4	bottom_set;	   // bottom telescoping set for optimization.
	char	*column;
	UInt4	set_cols[9];
	UInt4	*card_set;
	BooLean	GlobalSqWts;
	UInt4	*AveSqIWt;	// averaged weight for each sequence.
	UInt4	*Card_Set;	// weight number of seq. in each set.
	FILE	*debug_fp;
#if 0
	static const Int4 startAB=0; //  treat deletions ('-'='x') as non-functional residues.
#else
	static const Int4 startAB=1;	// treat seqs with deletions ('-'='x') as non-existent
#endif
};

#endif

