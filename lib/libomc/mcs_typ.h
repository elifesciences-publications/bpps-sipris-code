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

#if !defined (_MCS_TYP_)
#define _MCS_TYP_

#include "probability.h"
#include "cmsa.h"
#include "chn_typ.h"
#include "rst_typ.h"
#include "che_typ.h"
#include "hpt_typ.h"
#include "pch_typ.h"
#include "mheap.h"
#include "swt_typ.h"
#include "swaln.h"
#include "tax_typ.h"
#include "lpr_typ.h"	// DEBUG...
#include "dms_typ.h"

#define MCS_USAGE_START "FATAL: parameter setting syntax error within <infile_prefix>.hpt:\n\
   USAGE: GroupName [options}\n\
     -A<real>:<real> - alpha hyperparameters A0:B0 (default: A0=B0=1.0)\n\
     -concise      - Don't create extra output files (for public version)\n\
     -col=<int>:<int>  - Specify the min and max number of columns allowed\n\
     -global       - Use global sequence weighting\n\
     -N=<int>      - Maximum number of significant pattern positions to highlight\n\
                   - (This sets the contrast for the alignment.)\n\
     -P=<str>      - seed pattern string\n\
     -Ri=<real>    - Set prior probability that a row (seq) is in the foreground (default: 0.5)\n\
     -rho=<real>   - set prior probability (rho) that a column is a pattern position (default: 0.5)\n\
                        input is log(rho); e.g., rho=100 --> rho=e^-100 = column 'penalty' of -100 nats\n\
   \n\n"

class mcs_typ {         // multiple category sampler (partitioning with pattern selection)
public:
	mcs_typ( ){ PrintUsage(stderr); }
	mcs_typ(cma_typ,cma_typ,hsw_typ,                                Int4, char *argv[]);
	mcs_typ(cma_typ,cma_typ,hsw_typ,Int4,set_typ*,                  Int4, char *argv[]);
	mcs_typ(cma_typ,cma_typ,hsw_typ,              hpt_typ*,cma_typ*,Int4,char *argv[]);
	mcs_typ(cma_typ in_cma, cma_typ in_mcma, hsw_typ hsw, Int4 NumInSets, set_typ *InSet,
        	hpt_typ *in_hpt, cma_typ *in_sma, Int4 argc, char *argv[]);
	// void	DisOwnInSMA(){ own_in_sma=TRUE; }
	Int4	GetSetSize(){ return SetN(Labeled); }
	Int4	NumRejected(){
			return CardSet(GrpSet[Hpt->NumSets()])-CardSet(RandomSet);
		}
	Int4	SetCard(Int4 s){
			assert(s > 0 && s <= Hpt->NumSets());
			if(GrpSet[s]) return CardSet(GrpSet[s]);
			else return 0;
		}
	cma_typ	MkMainFileCMSA(cma_typ cma, Int4 num_random);
	void	MoveSeqs(Int4, Int4);
	~mcs_typ( ){ Free( ); }
	//**************** rtf routines ***********************
	Int4	font_size;	// default 6 points.
	void	SetFontSize(Int4 f){ assert(f >= 4 && f <= 25); font_size=f; }
	char	page_format;	// l,p,L,P; default 'P'.
	void	SetPageFormat(char x){ assert(strchr(" lLpP",x) != NULL); page_format=x; }
	//**************** omc_resume routines ***********************
	Int4	AssignFakeSeq(Int4 sq, cma_typ xcma);
	Int4    AssignFakeSeq(unsigned char *seq, Int4 Len);
	//**************** omc_resume routines ***********************
	//*************************** checkpoint routines ***************************

	void	PutPttrnLLR(FILE *fp, Int4 n) 
		   { assert(n > 0 && n <= Hpt->NumBPPS()); che[n]->PutPttrnFile(fp,n); }
	void	PutSubLPRs(FILE *fp, Int4 n) 
		   { assert(n > 0 && n <= Hpt->NumBPPS()); che[n]->PutSubLPRs(fp); }
	void	PutSubLPRs(FILE *fp) {
		  for(Int4 n=1; n<= Hpt->NumBPPS(); n++){ che[n]->PutSubLPRs(fp); }
		}
	//**************** lpr routines ***********************
	double  GetPriorRi(Int4 n, char mode='C'){
                    double  sq_prior=0.0;
		    if(mode=='C'){ // weight sequence priors by # children for each node.
                        Int4    nn = Hpt->NumChildren(n);
                        if(n==1) sq_prior=0.5/(1.0 + (double)nn);
                        else sq_prior=0.999/(1.0 + (double)nn);
		    } else if(mode=='E'){ // equal weights on all...
                        sq_prior=1.0/(double)Hpt->NumSets();
		    } else if(mode=='R'){ // equal weights on all except for root node.
                        if(n==1) sq_prior=0.5; else sq_prior=0.5/(double)Hpt->NumSets();
		    } return sq_prior;
                }
	void	PutPttrnLLRs(FILE *fp){
		  for(Int4 n=1; n<= Hpt->NumBPPS(); n++)
			{ che[n]->PutCDTreePttrnFile(fp,Hpt->ElmntSetName(n)); }
		}

	//**************** mcs_rtrn.cc ***********************
	e_type	*RtnKeySeqs( );
	e_type	RtnKeySeq(Int4 n);

	sst_typ *RtnCopyOfWorkingSST(Int4 n);
	sst_typ *RtnCopyOfBestSST(Int4 n);
	sst_typ	**RtnCopyOfSSTs( );
	sst_typ **RtnCopyOfWorkingSSTs( );
	sst_typ **RtnCopyOfBestSSTs( );

	char	**RtnCopyOfPttrns( );
	sst_typ	*RtnCopyOfSST(Int4 n);
	char	*RtnCopyOfPttrn(Int4 i);
	e_type	RtnCopyOfKeySeq(Int4 i,char *name){
		   assert(i > 0 && i < Hpt->NumSets());
		   char tmp_str[100]; sprintf(tmp_str,"%s consensus",name);
		   e_type cE=CopySeq(che[i]->KeyE( ));
		   ChangeInfoSeq(tmp_str,cE); return cE;
		}
	double  **RtnCopyOfLPRs( );	// to be called after sampling is completed.
	cma_typ RtnCsqAsCMSA(Int4 n, char *name,sst_typ *xsst=0);
	cma_typ RtnCsqSstAsCMSA(Int4 n, char *name,sst_typ *xsst=0);
	cma_typ RtnSeqAsCMSA(Int4 n, char *name,sst_typ *xsst=0);
	cma_typ RtnBstAsCMSA(Int4 n, char *name,sst_typ *xsst=0);
	double	RtnNatsPerWtSeq(Int4 i){ 
		   assert(i > 0 && i <= Hpt->NumBPPS());
		   return che[i]->RtnNatsPerWtSq();
		}
	Int4	NumColumns(Int4 i){
		   assert(i > 0 && i <= Hpt->NumBPPS());
		   return che[i]->NumColumns();
		}
	Int4	MinNumColumns(Int4 i){
		   assert(i > 0 && i <= Hpt->NumBPPS());
		   return che[i]->RtnMinNumColumns();
		}
	void	SetMaxNumColumns(Int4 i,Int4 mc){
		   assert(i > 0 && i <= Hpt->NumBPPS());
		   assert(mc >= che[i]->RtnMinNumColumns());
		   che[i]->SetMaxNumColumns(mc);
		}
	e_type	EmitRandomSeq(Int4 i){
		   assert(i > 0 && i <= Hpt->NumBPPS());
		   return che[i]->EmitRandomSeq();
		}
	Int4	ComputeNumRandom(Int4 N){ return 1+(N/3); }
	// need information on how to compute set size available to the outside. 
	void	SetNthSeqForDisplay(Int4 n){ assert(n > 0); NthSeqForDisplay=n; }
	void	SetNumHighlighted(Int4 n){ assert(n > 0); NumHighlighted=n; }
	void    VerboseOff(){ verbose=FALSE; }
        void    VerboseOn(){ verbose=TRUE; }
private:
	BooLean	verbose;
	Int4	NthSeqForDisplay;
	Int4	NumHighlighted;
	BooLean	del_as_random;
	Int4	MinimumSetSize;
	double	MinimumNatsPerWtSeq;
	double	MinimumLLR;
	set_typ	Skip;
public:
	double	PutBoltzmannLike(FILE *fp,Int4 n,Int4 c);	// in omc_debug.cc
	//**************** mcs_rtrn.cc ***********************
	Int4	RtnNumElmntSetCMA( ){ return Hpt->NumSets(); }
	Int4	RtnLengthMainCMSA( ){ return LengthCMSA(1,MainCMA); }
	Int4	RtnNumCategories( ){ return Hpt->NumBPPS(); }
	set_typ	*CopyOfSeqSets(){ return CopyOfSeqSets_Private(); }
	set_typ	*CopyOfBestTreeSets(){ return CopyOfBestTreeSets_Private(); }
	set_typ	*CopyOfPartitionSets(){ return CopyOfPartitionSets_Private(); }
	bpps_typ *RtnBPPS(Int4 n){ assert(n > 0 && n <= Hpt->NumBPPS()); return che[n]->BPPS(); }
	che_typ *RtnChe(Int4 n){ assert(n > 0 && n <= Hpt->NumBPPS()); return che[n]; }
	che_typ **RtnChe( ){ return che; }
	e_type	RtnQryBPPS(Int4 n)
		   { assert(n > 0 && n <= Hpt->NumBPPS()); bpps_typ *bpps=che[n]->BPPS();
			return bpps->RtnQuerySeq(); }
	//----------------------------------------------------

	//**************** mcs_cdd.cc ***********************
	void    PutCDD(FILE *ofp);
	Int4    FirstInSet(set_typ St,cma_typ cma);
	BooLean IsSeqEST_ENV(e_type E);
	Int4    *RtnBestSeqs(set_typ *&GoodSet);
	set_typ RtnPdbSet(Int4 g);
	Int4    NoIndelsSet(set_typ NoIndels);
	//**************** mcs_cdd.cc ***********************
	
	//**************** mcs_put.cc ***********************
	BooLean PutSeqCntrb(FILE *fp,FILE *hfp,FILE *sfp);
	BooLean PutMapContributions(FILE *fp,Int4 g);
	double  NodeDiversity(FILE *fp, Int4 i, set_typ Set);
	void    PutSubHierarchyCntrb(FILE *fp);
	double  Put( ){ Put(TRUE); }
	double  Put(BooLean x, char **names=0, BooLean do_put_cdd=FALSE,FILE *mmafp=0,FILE *ptrnfp=0);
	void    PutRTF(BooLean updateCSQ,BooLean SaveChnFiles, Int4 KeyStart=0,
							Int4 KeyEnd=0,char *KeySetName=0, double **value=0);
	void    PutRTF(BooLean updateCSQ){ PutRTF(updateCSQ,PutIntermediateFiles); }
	double  **GetResEvals(Int4 n);
	void    PutLpr(char *outfile);
	void    PutLpr(FILE *fp);
	void    PutPttrns(FILE *fp);
	void	SetUserDisplaySet(){ user_display_set=TRUE; }
	Int4	PutMajorNodesMMA(FILE *mmafp);
private:
	BooLean	user_display_set;
	BooLean PutEnhancedDisplaySet(FILE *sfp);
public:
	//**************** mcs_put.cc ***********************

	//**************** mcs_typ.cc ***********************
	BooLean	MvUpContribution(FILE *fp, Int4 g,Int4 n, double &Rtn);
	BooLean MvDownContribution(FILE *fp, Int4 node,Int4 sibling,double &Rtn);
	BooLean QuickMvDownContrib(FILE *fp, Int4 node,Int4 sibling,double &Rtn);
	BooLean	PutMapContributions(FILE *fp,lpr_typ *xlpr=0);
	cma_typ	*RtnLeafNodeCMSA(Int4 &Number);	// return the 
	BooLean RtnContribLLR(Int4 row, Int4 col, double &subLpr);
	double	RestoreBest();
	double	RevertToBest();
	void	StoreBest();
	double  CalcTotalLPR(FILE *fp,BooLean StoreBestOK);
	//**************** mcs_typ.cc ***********************
	BooLean	SaveBest;
	double  *RtnSubLPR( ){ CalcTotalLPR(0,FALSE); return Map; }
	double  RtnSubLPR(Int4 n)
		    { assert(n > 0 && n <= Hpt->NumBPPS()); return che[n]->CalcLLR( ); }
	double  CalcTotalLPR(FILE *fp){ return CalcTotalLPR(fp,SaveBest); }
	double  CalcTotalLPR( ){ return CalcTotalLPR(0,SaveBest); }
	double  CalcNoFailLPR( ){ return CalcNoFailLPR(0); }
	double  CalcNoFailLPR(FILE *fp){ 
			double rtn=CalcTotalLPR(fp,FALSE); 
			for(Int4 n=1; n<= Hpt->NumBPPS(); n++) if(Map[n] <= 0.0) return 0.0;
			return rtn;
		}
	double	RtnBestLPR(){ return BestLPR; }
	//-----------------------------------------------------

	//**************** lpr routines ***********************

	//**************** return results routines ***********************
	double	Temperature( ){ return temperature; }
	Int4	ModStart,TargetMod;
	void	ResetTemperature(double T){ assert(T >= 0.0); temperature=T; }
	Int4	GetNumRandom(){ return NumRandom; }
	hpt_typ	*GetHpt( ){ return Hpt; }
	//**************** return results routines ***********************
	//**************** output routines ***********************
	void    PutAllSubLPRs(FILE *fp){
			for(Int4 i=1; i<=Hpt->NumBPPS(); i++){ assert(che[i]); che[i]->PutSubLPRs(fp);  }
		} 
	void    ReSetSubTree(set_typ subtree, Int4 c){ ReSetSubTree(subtree,c,this->Hpt); }
	void    ReSetSubTree(set_typ subtree, Int4 c, hpt_typ *hpt){
                    ClearSet(subtree);
                    for(Int4 row=1; row < hpt->NumSets(); row++)
                        { if(hpt->Cell(row,c) == '+') AddSet(row,subtree); }
                }
	FILE	*GetOutFP( ){ return outfp; }
	void	PutHpt(FILE *fp){ Hpt->Put(fp); }
	void	PutDisplayCMA(FILE *fp){ for(Int4 n=1; n <= NumDisplayCMA; n++) PutCMSA(fp,DisplayCMA[n]); }
	void    PutCheckPoint( ){
			FILE *cfp=open_file(infile,".hpt","w"); PutHpt(cfp);  fclose(cfp);
			cfp=open_file(infile,".sma","w"); PutDisplayCMA(cfp);  fclose(cfp);
		}
	void    PutHyperPartition( ){
			if(outfp == 0){ outfp = open_file(infile,".out","w"); }
			PutHyperPartition(outfp); fflush(outfp);
			// PutHyperPartition(stderr); fflush(stderr);
		}
	void    PutHyperPartition(FILE *fp);
	//**************** output routines ***********************

	//**************** mcs_debug.cc ***********************
	void	PutSetRelations(FILE *fp);
	void    CheckPttrnCsqMatch(char *msg);
	BooLean ConsistencyCheck();
	BooLean CheckInputSets();
	BooLean	ChecksOut();  // check to see whether anything worthwhile was found.
	BooLean	PutMapContributions2(FILE *fp,lpr_typ *xlpr=0);
	void    PutPttrnVsConsSeq(FILE *fp,char *msg);
	void    PutPttrnVsConsSeq(FILE *fp, Int4 i);
	//**************** mcs_debug.cc ***********************

	Int4	SetSize,SizeMain,SizeTrueMain,LengthMain;

	//************************** msc_pttrn.cc ****************************
	Int4	LoadUpColumns( );
	double  SampleColumns(BooLean UseNegCol, char mode=' ');
	double  SampleColumns(){ return SampleColumns(FALSE); }
	Int4    SampledColumn;		// SampleHpt on this column or none if == 0.
	sst_typ *PruneSST(Int4 n, Int4 k, sst_typ *isst);
	BooLean IsConflict(Int4 k, Int4 better, Int4 worse);
	BooLean IsConflict0(Int4 k, Int4 better, Int4 worse);	// old version when modifying.
	Int4    LoadUpBestColumns(Int4 NumToKeep);
	//************************** msc_pttrn.cc ****************************
	void	RmWorstColumn(Int4 n, Int4 NumToKeep)
	    { while(che[n]->NumColumns( ) > NumToKeep){ if(che[n]->ForceRmWorstColumn( ) == 0) break; } }

	//**************** sampling routines (msc_sample.cc) ***********************
	BooLean	Sample( ) { return Sample(IterStart_DF,IterEvolve_DF,NumRounds_DF); }
	BooLean	Sample(Int4 x,Int4 y, Int4 z){ return Sample(x,y,z,ColSampleStart_DF); }
	BooLean	Sample(Int4 x,Int4 y, Int4 z, Int4 c){ return Sample(x,y,z,c,ModStart); }
	BooLean	Sample(Int4 ,Int4, Int4, Int4, Int4);
	double	SampleSeq(Int4 sq,double lpr);
	BooLean ResurrectRejectSeq(Int4 sq);
	set_typ	SampledSet;		// Sample sequences in these groups (or in all sets if == 0).
	//**************** sampling routines ***********************

	//************************** msc_updown.cc ****************************
	BooLean MvUpDown(Int4 grandpa, Int4 parent,Int4 child, set_typ subtree, char State);
	BooLean MoveUp(Int4 grandpa, Int4 parent,Int4 child, set_typ subtree)
		{ return MvUpDown(grandpa, parent,child,subtree,'-'); }
	BooLean MoveDown(Int4 grandpa, Int4 parent,Int4 child, set_typ subtree)
		{ return MvUpDown(grandpa, parent,child,subtree,'+'); }
	//************************** msc_updown.cc ****************************

	//**************** test routines (msc_test.cc) ***********************
	void    PutContinueFile(FILE *fp,BooLean Label); // not currently used.
	Int4	RemoveSimilarSets( ); // Just for testing right now.
	// BooLean	SampleHpt(){ wdg_typ X=0; return SampleHpt(0,0,X); }
	// BooLean	SampleHpt(FILE *fp){ wdg_typ X=0; return SampleHpt(fp,0,X); }
	BooLean	SampleHpt(FILE *fp,Int4 Root, wdg_typ &Tree);
	BooLean SampleHpt(FILE *fp,Int4 SampledCol);
	//**************** test routines ***********************
	double  *ContributionsToLLR(BooLean global=FALSE);

	BooLean	SaveSets;
	BooLean	IsTreeHpt; // Does the input FD-table correspond to a tree? e.g., for pmcBPPS program.
	BooLean	SpeakUp(){ Verbose=TRUE; }
	void	RenameInfile(char *str){ if(infile) free(infile); infile=AllocString(str); }
	void	DoNotEvolve(){ Evolve=FALSE; }
	void	DoEvolve(){ Evolve=TRUE; }
	void    UpdateDisplaySeqs( );
	void	RenameDisplayCMA(Int4 n,char *name)
		    { assert(n > 0 && n <= NumDisplayCMA); RenameCMSA(name,DisplayCMA[n]); }
	double	MinNatsPerWtSq( ){ return this->MinimumNatsPerWtSeq; }
	BooLean	OkayNatsPerWtSq(Int4 i){
			if(this->RtnNatsPerWtSeq(i) >= this->MinimumNatsPerWtSeq) return TRUE;
			else return FALSE;
		}
	Int4	NumRawFailedNodes(Int4 MinSize,double minLLR) {
		   Int4 n,i; 
		   CalcTotalLPR(0,FALSE);
		   for(n=0,i=2; i<= Hpt->NumBPPS(); i++){ 
			if(Skip && MemberSet(Hpt->ItoSetID(i),Skip)) continue;
			else if(Hpt->TypeOfSet(i) == '!' && this->SetCard(i) < MinSize) n++;
			else if(Map[i] < minLLR) n++;
			else if(che[i]->RtnNatsPerWtSq() < MinimumNatsPerWtSeq) n++;
		   } return n; 
		}
	Int4	NumFailedNodes( )
		{ Int4 n,i,id;
		  CalcTotalLPR(0,FALSE);
		  for(n=0,i=1; i<= Hpt->NumBPPS(); i++){
			if(Skip && MemberSet(Hpt->ItoSetID(i),Skip)) continue;
			else if(IsFailedBPPS[i]) n++;
			else if(Hpt->TypeOfSet(i) == '!' && this->SetCard(i) < MinimumSetSize) n++;
			else if(che[i]->RtnNatsPerWtSq() < MinimumNatsPerWtSeq) n++;
			else if(Map[i] < MinimumLLR) n++;
			else if(che[i]->NumColumns( ) < che[i]->RtnMinNumColumns( )) n++;
		  } return n; 
		}
	BooLean	SetStringency(double MinLLR, Int4 MinSetSize, set_typ skip=0){
			MinimumSetSize=MinSetSize; MinimumLLR=MinLLR; Skip=skip; }
	Int4	NumberFailedNodes( ){
		    Int4 n=0,i;
		    for(i=1; i<= Hpt->NumBPPS(); i++){
			double d=che[i]->CalcLLR();
			if(d <= 0.0 || che[i]->NumColumns( ) < 1) n++;
		    }  return n;
		}
	BooLean	NoFailureMode;
	BooLean	IsBestRestored( ){ return DidRestoreBest; }
	// void	SortHpt(){ assert(IsTreeHpt); hpt_typ *hpt=Hpt->Sort( ); delete Hpt; Hpt=hpt; }
	// Need to eliminate copied Hpt info below first...
	Int4	FileID;
	double	MaximumTemperature(){ return MaxTemperature; }
	Int4	RtnDefaultMaxCol(){ return DefaultMaxCol; }
	Int4	RtnDefaultMinCol(){ return DefaultMinCol; }
#if 0
	Int4	MaxNumColumns(Int4 x)
		   { assert(x > 0 && x <= Hpt->NumBPPS()); return MaxNumCol[x]; }
#endif
	Int4	PrintToggle;
private:
	Int4	IterStart_DF,IterEvolve_DF,NumRounds_DF,ColSampleStart_DF;
	BooLean	PartitionRandomly;
	BooLean	Evolve;
	void	UpdateCSQ( ){ UpdateCSQ(0); }
	void	UpdateCSQ(FILE *fp);
	//**************** initialization routines ***********************
	void    RandomlyPartitionSets( );
	void	PartitionBySeedAlnCMSA(Nlm_Int4, colinearmaln_type**, char*);
	void    PartitionByInputSetCMSA( );
        void    ComputeSeedCMAs();
        Int4    MinSeed2CsqScore[MAX_NUM_ELMENTARY_SETS];
        Int4    Seq2SeedCsqScore(Int4 sq,Int4 set);
        Int4    *SortByScoreCMSA(FILE *fp, char mode, Int4 &first_best, double cut, 
			cma_typ cma, Int4 set);
	char    *FindSeedPattern(Int4);

	//**************** sampling routines (msc_sample.cc) ***********************
	void	TransferAllSeqs(Int4 from, Int4 to);
	void	ReSetRelations() { return ReSetRelations(0); }
	void	ReSetRelations(FILE *fp);
	set_typ	*CopyOfPartitionSets_Private();
	set_typ	*CopyOfSeqSets_Private();
	set_typ	*RtnSeqSets();
	set_typ *RtnSubTreeSeqSet( );
	set_typ RtnSubTreeSeqSet(Int4);
	set_typ *CopyOfBestTreeSets_Private( );
	Int4    FindSetForRemoval(Int4 failed_column);
	Int4    RmUnfruitfulSets( );
	BooLean	Unlabeled;

	//********************* mcs_typ.cc ***********************
public:
	void	CopySubTreeSeqSet(Int4 node, set_typ xSet);
	void    CopyBkGrndSeqSet(Int4 node, set_typ rtnSet);
	void    UnLabelAllSeqs( );
	BooLean	IsInSet(Int4 sq, Int4 st){
		   assert(st > 0 && st <= Hpt->NumSets()); 
		   assert(sq > 0 && sq <= NumSeqsCMSA(TrueMainCMA));
		   return MemberSet(sq,GrpSet[st]);
		}
	char	RtnRelateFG(Int4 i, Int4 j){
		   assert(i > 0 && i <= Hpt->NumSets()); assert(i > 0 && i <= Hpt->NumBPPS());
		   assert(j > 0 && j <= Hpt->NumSets()); assert(j > 0 && j <= Hpt->NumBPPS());
		   return RelateFGs[i][j];
		}
	char	RtnRelateBG(Int4 i, Int4 j){
		   assert(i > 0 && i <= Hpt->NumSets()); assert(i > 0 && i <= Hpt->NumBPPS());
		   assert(j > 0 && j <= Hpt->NumSets()); assert(j > 0 && j <= Hpt->NumBPPS());
		   return RelateBGs[i][j];
		}
private:
	//********************* mcs_typ.cc ***********************

	//****************** Initialization routines: ******************
	void    Init(Int4 argc, char *argv[]);
	void	InitDefaults();
	void	InitFlags();
	void	InitAsNull();

	//**************** mcs_arg.cc ***********************
	void    PrintUsage(FILE *fp);
	void	ReadMainArg(Int4 argc, char *argv[]);
	//**************** mcs_arg.cc ***********************

	Int4    ReadArgFile();
	Int4    GetElmntSets( );
	Int4	ReadSeedPttrns( );

	//**************** mcs_init.cc ***********************
	Int4    SetUpNthSrch(Int4 n, Int4 argc,char *argv[]);
	Int4	SemiConvergedState;		// decide when to give up on failed sets.
	//****************** Initialization routines: ******************

	void    Free();         // free memory.
	Int4	TotalColumns( ){
		   Int4 n,nCol=0;
		   for(n=1; n<= Hpt->NumBPPS(); n++) nCol += che[n]->NumColumns( ); return nCol; }

	//****************** input parameters: ******************
	char	set_mode;	// mode for rst_typ; 'L' by default.
	BooLean	PrintEachRTF;

	//************* hpt_typ routines ********************
	hpt_typ	*Hpt;	// hyperpartition type. (to replace arg_typ).
	char	*sst_str[MAX_NUM_ELMENTARY_SETS];
	Int4	MaxNumberBPPS;
	BooLean	IsFailedSet[MAX_NUM_ELMENTARY_SETS];
	BooLean	IsFailedBPPS[MAX_NUM_ELMENTARY_SETS];
	Int4	NumRandom;
	char    **GetSetRelations(char *Title,Int4 *nGrpsX, Int4 **GrpsX,set_typ **RtnSetX);
	char	**RelateFGs;	// subset[n1][n2]='<',superset='>',Intersect='+',Disjoint='0'.
	char	**RelateBGs;	// subset[n1][n2]='<',superset='>',Intersect='+',Disjoint='0'.
	set_typ	RandomSet;
	void	MakeRandomSet();
	//************* hpt_typ routines ********************

	//*************** Saving the Optimum Sets found *******************
	set_typ InitSet[MAX_NUM_ELMENTARY_SETS];

	hpt_typ	*BestHpt;
	e_type	BestCsq[MAX_NUM_ELMENTARY_SETS];	// == Best che[n]->KeyE();
	BooLean IsFailedBestSet[MAX_NUM_ELMENTARY_SETS];
        BooLean IsFailedBestBPPS[MAX_NUM_ELMENTARY_SETS];
	set_typ BestSet[MAX_NUM_ELMENTARY_SETS];
	sst_typ	*best_sst[MAX_NUM_ELMENTARY_SETS];
	double	BestLPR;
	BooLean	DidRestoreBest;
	// hyp_typ BestHpt;

	set_typ	GrpSet[MAX_NUM_ELMENTARY_SETS];
	double  *SubLPR[MAX_NUM_ELMENTARY_SETS];
	double	TotalLPR;
	double	Map[MAX_NUM_ELMENTARY_SETS];
	BooLean CheckValue(double x);
	set_typ *SetFG;
	set_typ *SetBG;
	set_typ	Labeled;		// these are labeled as belonging to a fixed set.
	//*************** Saving the Optimum found *******************

	//==================== InputSets =========================
	Int4	num_passed_in_sets;
	set_typ	*passed_in_sets;
	//==================== InputSets =========================

	//==================== Temporary arrays =========================
	Int4	*WorstToBest[MAX_NUM_ELMENTARY_SETS];
	Int4	Index1stBest[MAX_NUM_ELMENTARY_SETS];
	//==================== Temporary arrays =========================

	//==================== MSAs =========================
	// BooLean	own_in_sma;
	Int4	num_passed_in_sma;
	cma_typ	*passed_in_sma;
	Int4	NumDisplayCMA;
	cma_typ	*DisplayCMA;
        cma_typ SeedCMA[MAX_NUM_ELMENTARY_SETS];	// Display set + consensus seq.
	cma_typ	*QryCMAs;	// used to pass query sequences to chn_typ.
	cma_typ	passed_in_cma,TrueMainCMA;	// TrueMainCMA always passed in..
	cma_typ	passed_in_mcma,MainCMA;		// MainCMA always passed in..	
	cma_typ	dummyCMA;
	hsw_typ	passed_in_hsw;	// with random sequences added.
	//==================== MSAs =========================

	//****************** passed in data: ******************
	hpt_typ	*passed_in_hpt;
	char	*program_name;
	char	*infile;
	a_type  AB;	// always passed in with in_cma...
	//****************** passed in data: ******************

	cma_typ RtnSeqOrCsqAsCMSA(char mode, Int4 n, char *name,sst_typ *xsst=0);

	void	GetChnFiles();
	chn_typ	**chn;
	che_typ	**che;	// == old che_typ modified for cdh_typ.
	sst_typ ***sst;
	char    *SFBG;
	BooLean	PutIntermediateFiles;

	//=============== parameter settings. ===================
	double	GlobalA0,GlobalB0,GlobalRi;
	double	MiscGlobalA0,MiscGlobalB0,MiscGlobalRi;
	double	RootGlobalA0,RootGlobalB0,RootGlobalRi;
	Int4	DefaultMaxCol,DefaultMinCol; // ,DefaultContrast;
	Int4	GlobalN;
	double	Global_rho,MinNats;
	Int4	SeedPttrnLen;
	double	temperature,MinTemperature,MaxTemperature;
	UInt4	ppb_increase;	// limit for

	Int4	*MaxNumCol;
	BooLean	NoSeeds;	// TRUE -> call FindSeedPatterns heuristic.
	BooLean	NoCSQ;		// TRUE -> don't add a consensus sequence to seed alignments.
	BooLean	Verbose;
	UInt8	Iteration;
	Int4	NumCalls;	// how many times was Sampling called.
	BooLean	StrictIndepend;
	//=============== parameter settings. ===================
	
	//=============== FILE pointers. ===================
	FILE	*cfp;		// convergence file pointer (iteration, LPR, temperature).
	FILE	*ifp;		// sampling iteration file pointer with cardinality of sets.
	FILE	*efp;		// stderr file pointer.
	FILE	*outfp;		// infile.out file pointer.

	//********************* mcs_junk.cc ***********************
	double	SampleLeafParent(Int4 n);
	void	AlignToKeyE(Int4 n, e_type E){ AlnSeqSW(11, 1, che[n]->KeyE( ),E,AB); }
public:
	void	SetTheCSQ(FILE *fp,Int4 n, e_type csqE=0);
	//********************* mcs_junk.cc ***********************
	//----------------------------------------------------
};

Int4    RunmcsBPPS(Int4 argc,char *argv[]);

#endif

