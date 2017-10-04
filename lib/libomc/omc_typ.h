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

#if !defined (OMC_TYP)
#define OMC_TYP

#include "stdinc.h"
#include "afnio.h"
#include "tree2hpt.h"
#include "random.h"
#include "dheap.h"
#include "mheap.h"
#include "set_typ.h"
#include "dsets.h"
#include "probability.h"
#include "cmsa.h"
#include "chn_typ.h"
#include "lpr_typ.h"
#include "sqd_typ.h"
#include "che_typ.h"
#include "hpt_typ.h"
#include "pch_typ.h"
#include "mcs_typ.h"
#include "mad_typ.h"
#include "hsi_typ.h"
#include "stack.h"
#include "rsq_typ.h"
#include "bpcp_typ.h"
#include <iostream>

class omc_typ  {	// optimizing multiple category BPPS sampler type.
public:
                omc_typ(char sc,Int4 argc, char *argv[],FILE *mfp=0,FILE *hfp=0,FILE *stfp=0,FILE *smfp=0)
		{ mmafp=mfp; hptfp=hfp; setfp=stfp; smafp=smfp; 
			SecretCode=sc; program_name=AllocString(argv[0]);
			if(mfp!=0) assert(hfp!=0 && stfp!=0 && smfp!=0); Init(argc,argv); }
                omc_typ(Int4 argc, char *argv[],FILE *mfp=0,FILE *hfp=0,FILE *stfp=0,FILE *smfp=0)
		{ mmafp=mfp; hptfp=hfp; setfp=stfp; smafp=smfp; 
			SecretCode=0; program_name=AllocString(argv[0]);
			if(mfp!=0) assert(hfp!=0 && stfp!=0 && smfp!=0); Init(argc,argv); }
                ~omc_typ(){ free(program_name); Free(); }
	Int4	Run(FILE *fpmma=0, FILE *fphpt=0, FILE *ptrnfp=0);		// omc_run.cc
	Int4	Put(BooLean SkipCheckPoint=FALSE,FILE *fpmma=0,FILE *fphpt=0,FILE *ptrnfp=0);	// omc_typ.cc
	void    PutCheckPoint(char *name,BooLean change=TRUE); // in omc_typ.cc
	void	PutToggleCheckPoint( ){
			char	c='j';
			static UInt4 calls=0;
			calls++; if(calls % 2 == 0) c='k';
			if(verbose) sprintf(str,"%s_%c",infile,c); 
			else sprintf(str,"%s",infile);
			PutCheckPoint(str); 
		}
	void	PrintTime(FILE *fp)
		{ fprintf(stderr,"\ttime(omc%c): %d seconds (%0.2f minutes)\n",
                        InitMode,time(NULL)-time1,(float)(time(NULL)-time1)/60.0); }
	//=================== rtf options ==========================
	Int4    font_size;
	char    page_format;
	//=================== hieraln.cc ==========================
	void    PrintOptions(FILE *fpmma=0, FILE *fphpt=0,FILE *ptrnfp=0);
        set_typ *RtnSets( ){
			Int4 n;
			assert(InitMode == 's');
			FILE *fp = open_file(infile,".sets","r");
        		set_typ *sets=ReadSets(fp,n); fclose(fp);
			return sets; 
		}
        hpt_typ *RtnHpt( ){ return Hpt; }
        cma_typ RtnTrueCMA(){ return TrueMainCMA; }
	char	RtnInitMode(){ return InitMode; }
	a_type	RtnAlphabet(){ return AB; }
	void	FindKeyPositions( );
	void	FindKeyPositions2( );
	// UInt4	**EachNodeWtCnts( );
	double	*FindCrossConserved(Int4 key, FILE *ofp=stdout);
	void	PrintRTF(BooLean updateCSQ,BooLean SaveChnFiles=FALSE,Int4 KeyStart=0,
			Int4 KeyEnd=0,char *KeySetName=0);
	void	CompareBPPS_BILD( );
	Int4	*Correspond(omc_typ *that);
	void	VerboseOff(){ verbose=FALSE; }
	void	VerboseOn(){ verbose=TRUE; }
private:
	//=================== omc_init.cc ==========================
#if 1	// <infile>.chk
	void	FromChkPtInit(cma_typ *chk_sma);
	void	FromSetsInit(cma_typ *chk_sma=0, set_typ *chk_sets=0);
#endif
#if 1	// moving nodes to (as a child of) a sibling.
	Int4	MvNode,ToNode;
	Int4	DelNode;
#endif
	FILE	*mmafp,*hptfp,*setfp,*smafp;
	BooLean	verbose;
	void	Init(Int4 argc, char *argv[]);
	void    InitAsNull();
	void    InitDefaults();
	void	MkFileIdHeap( );
	void	FromFileInit();
	void    AbInitioInit( );
	Int4	SetDefaultArguments(Int4 ppb=1000000);
	void	SetStringency(Int4 x);
	void	PutStringency( );
	cma_typ	*OrderNodeSMAs(cma_typ *iSMA,BooLean renamed=TRUE);
	cma_typ	*RenameHptSMA(BooLean rename=TRUE, cma_typ *chk_sma=0);	// renames sets using identifiers.
	void	Free(); 
	char	*program_name;
	void	PrintError(char *);
	char	SecretCode;
	//=================== omc_init.cc ==========================
	BooLean	del_as_random;
	BooLean	use_usr_sma;
	BooLean	OutPutRTF;
	Int4	NthSeqForDisplay;
	Int4	NumHighlighted;
	Int4	MaxDepthHierarchy;
	UInt4	RandomSeed;
	Int4	DefaultMaxCol;
	Int4	stringency;
	BooLean	StrictIndepend;
	char	InitMode;
	double	*aafreq;
	struct	mp_type {	// Parameters for mad_typ (used by AddLeaf() )
		mp_type(double p1,double p2, double p3,double p4,Int4 p5,char p6)
		  { MinKeyFrq=p1; MaxGapFrq=p2; pcut=p3; Exact_pcut=p4; MinClique=p5; sets_mode=p6; }
		~mp_type( ){  }
		double  MinKeyFrq,MaxGapFrq;
		double	pcut;		// the maximum CumHypGeomProb() cutoff for merging patterns.
		double	Exact_pcut; 	// Exact test probability cutoff for identifying a significant pattern.
                Int4    MinClique;	// the minimum size of a clique of patterns to merge (?)
		char	sets_mode;
	};
	mp_type	*MP;

	//=================== omc_resume.cc ==========================
public:
	set_typ	*Resume(cma_typ xcma, Int4 *Mapping=0);
	Int4	CalcSetSize(Int4 &nrand, cma_typ xcma)
		{
		   nrand=1+(NumSeqsCMSA(xcma)/3);
		   Int4 set_size=NumSeqsCMSA(xcma) + nrand+1;
		   return set_size;	// MakeSet(set_size);
		}
	void	PutHyperPartition(FILE *fp){ mcs->PutHyperPartition(fp); }
	sst_typ	*RtnCopyOfSST(Int4 n){ return mcs->RtnCopyOfSST(n); }
private:
	//=================== omc_simulate.cc ==========================
	void    PutSimulatedAln( );
	void    PutAbInitioSimulatedAln( );

	//=================== omc_run.cc ==========================
	double  SproutLeaves(FILE *fp,double Temp, Int4 iter);
	double	AddInternalNodes(FILE *fp,double Temp, Int4 iter);
	double	TrimNodes(FILE *fp,double Temp, Int4 iter);
	double  GrowLeaves(FILE *fp,double Temp, Int4 iter, char Action='A');
	double  RaiseBranches(FILE *fp,double Temp, Int4 iter);
	double  LowerBranches(FILE *fp,double Temp, Int4 iter);
	Int4	ResurrectRejects(FILE *fp, double Temp=300.0,Int4 iter=1);

	Int4	DeleteWorstNode(double Temp);
	Int4	PruneTree(FILE *fp,double Temp, Int4 iter,Int4 Lid);
	void	PrintOperationHeader(FILE *fp,char *msg, Int4 iter);
	void	MixItUp(char mode, Int4 iter);
	double  RePartitionSeqs(FILE *fp, Int4 iter);
	double  RandomizeSeqAssign(FILE *fp, Int4 F, Int4 iter);
	void	StoreRestoreBest(mcs_typ *xmcs);

	double  RmInternalNodes(FILE *fp,double Temp, Int4 iter);
	Int4    DeleteInternalNodes(FILE *fp,double Temp);

	Int4    RunIter(Int4 iter,Int4 jter,double Temp,double temp);
	void    Sample(double T,Int4 p1,Int4 p2,Int4 p3,Int4 p4)
		  { Sample(T,p1,p2,p3,p4,mcs,0.03,0); }
	void    Sample(double T,Int4 p1,Int4 p2,Int4 p3,Int4 p4,mcs_typ *&xmcs,
							double min_gain=0.03,set_typ st=0);
	//----------------- working on --------------------
	double	MergeSiblings(FILE *fp,Int4 iter, double Temp=300.0);
	//=================== omc_run.cc ==========================
	char	PrintMode;	// used to print output only.
	char	*KeySetName;
	char	*NameSeedAln;
	char	**GetNames();
	Int4	KeyStart,KeyEnd;

	//=================== omc_typ.cc ==========================
	double	CalcLLR(mcs_typ *&xmcs);
	BooLean SampleMCS(mcs_typ *x_mcs, Int4 ID);
	mcs_typ *CreateMCS(hpt_typ *hpt,set_typ *set, cma_typ *in_sma, double Temp,char mode,
						Int4 pID, Int4 cID=0);
	void	RestoreFinalBest( );
	double	RestoreVeryBest( );
	Int4	ToggleOutFileName(BooLean NoSubID=FALSE);
	void	DeleteMCS(mcs_typ *x_mcs);
	Int4    DeletableNodes(mcs_typ *x_mcs,double minLLR, Int4 minCard);
	void	PutAsBestHpt(mcs_typ *x_mcs);
	//=================== omc_typ.cc ==========================
	mcs_typ	*best_mcs;
	double	Best_LPR;
	double	RevertToBest(){ double d=mcs->RevertToBest(); Hpt=mcs->GetHpt( ); return d; }
	void	RestoreBest(){ mcs->RestoreBest(); Hpt=mcs->GetHpt( ); }
	//----------------------------------------------------------------

	//=================== omc_up.cc ==========================
	Int4    MoveNodesUp(double Temp, double minLLR=0.0);
	mcs_typ *MoveUp(Int4 node, double MinDeltaLLR);
	double	OnTheFlyMvUp(Int4 i, Int4 *Parent);
	mcs_typ *FlattenHiearchy();
	//=================== omc_up.cc ==========================

	//=================== omc_down.cc ==========================
	Int4    MoveNodesDown(double Temp, double minLLR=0.0);
	mcs_typ *MoveDown(Int4 node, Int4 sibling, double MinDeltaLLR);
	double	OnTheFlyMvDown(Int4 i, Int4 n, Int4 *Parent);
	//=================== omc_down.cc ==========================

	//=================== omc_fuse.cc ==========================
	Int4	FuseNodes(double minLLR=0.0, double Temp=300.0);
	mcs_typ	*Fuse(Int4 node, Int4 sibling, double MinDeltaLLR);
	//=================== omc_fuse.cc ==========================

	//=================== omc_operate.cc ==========================
	Int4	AddLeaves(double Temp=300.0, char Action='A');
	Int4	DeleteNodes(double Temp=300.0){ return DeleteNodes(Temp,MinimumLLR,MinimumSetSize); }
	Int4	DeleteNodes(double Temp,double minLLR){ return DeleteNodes(Temp,minLLR,MinimumSetSize); }
	Int4	DeleteNodes(double Temp,double minLLR,Int4 minCard); 

	mcs_typ *InsertNode(Int4 parent,Int4 ID,hsi_typ *hsi, BooLean sample=FALSE);
	mcs_typ	*DeleteNode(Int4 nodeID,double T=300);
	mcs_typ	*AddLeaf(Int4 pID, Int4 ID,double T=300);

	mcs_typ	*Duplicate(mcs_typ *x_mcs){
			mcs_typ *tmp_mcs = mcs; mcs=x_mcs; Hpt=mcs->GetHpt();
			mcs_typ *rtn_mcs=this->RtnCopy( );
			rtn_mcs->NoFailureMode=FALSE; rtn_mcs->StoreBest(); rtn_mcs->RestoreBest();
			mcs=tmp_mcs; Hpt=mcs->GetHpt(); return rtn_mcs;
		 }
	mcs_typ	*Optimize( );	// Adds optimum sequences to displayed alignment.
	mcs_typ *RtnCopy(char mode='C');
	mcs_typ *Randomize(Int4 F){ return this->Operate('R',F,0); }
	mcs_typ *RePartition( ){ return this->Operate('c',0,0); }
	mcs_typ *BstRePartition( ){ return this->Operate('b',0,0); }
	mcs_typ *RandSqRePartition( ){ return this->Operate('r',0,0); }
	mcs_typ	*Operate(char action,Int4 pID,Int4 ID,double Temp=300, hsi_typ *hsi=0);
	e_type  SubGrpCsq(set_typ set);
	//=================== omc_operate.cc ==========================
	Int4	MinimumSetSize,MinimumSplitSize;
	double	MinimumLLR,MinimumPreLLR,MinimumNatsPerSeq;
	struct Info {
		hpt_typ *hpt;
		set_typ *set;
		cma_typ *sma;
#if 1
		sst_typ **xsst;
		char    **pttrn;
#endif
		double	T;
		Int4	pID;
		void	Free(),Init(Int4 N);
		void	PrintSetSizes(FILE *fp){
			   if(set==0) return;
			   for(Int4 i=1; set[i]; i++) fprintf(fp,"%d: %d seqs\n",i,CardSet(set[i]));
			}
		void	Move(Int4 pI,Int4 I) // // set[pI] = set[pI] U set[I]; set[I]=null.
			  { UnionSet(this->set[pI],this->set[I]); ClearSet(this->set[I]); }
	};
	//----------------------------------------------------------------

	//=================== omc_addleaf.cc ==========================
	Int4    AddFocusedLeaf(double Temp);
	mcs_typ *AddFocusedLeaf(Int4 pID, Int4 ID,double T);
	BooLean	GrowFocusedLeaf(FILE *fp,double Temp, Int4 iter);
	BooLean ReSetUpFocusedSrch(Int4 NewI);

	e_type  *FindDivergentSeqs(Info *info, Int4 pI, Int4 I,sst_typ *osst,FILE *efp=0);
	BooLean AddNewLeaf(Info *info, set_typ *set, Int4 pID,Int4 ID);
	//=================== omc_addleaf.cc ==========================

	//=================== omc_insert.cc ==========================
	Int4	InsertInternalNodes(FILE *fp,double Temp); 
	Int4	SampleInternalNode(Int4 parent, hsi_typ *root,double Temp=300.0);
	hsi_typ	*OnTheFlyLPR_dfs(Int4 , Int4, Int4);
	sst_typ *ComputeHybridLLR(Int4 parent, double &Lpr);
	hsi_typ *CheckForBottom(Int4 depth, Int4 maxdepth,Int4 parent);
	void	ReSetIDSubTree(set_typ subtree, Int4 c);
	void	ReSetSubTree(set_typ subtree, Int4 c, hpt_typ *hpt);
	mcs_typ *OptimizeInsertedNode(hsi_typ *hsi);
	Int4    ForcedInternalNode(Int4 pID, hsi_typ *subroot);
	Int4	ForcedP,*ForcedNode,NumForced;
	//=================== omc_insert.cc ==========================
	void	ReSetSubTree(set_typ subtree, Int4 c){ ReSetSubTree(subtree,c,this->Hpt); }
	lpr_typ	*lpr;
	h_type	dfsHG,preHG;
	set_typ	ChildSetFG,ChildSetBG,ChildSetU,ChildSetBoth;
	set_typ SetFG,SetBG,TmpSet,TmpFG,TmpBG;	// FG and BG sets for each node (compute solely based on tree). 
	set_typ	*Set;		// array of sequence sets corresponding to nodes.
	sst_typ *GetOptPttrnLPR(FILE *f,set_typ S1, set_typ S2,double &L,Int4 x, unsigned char *&csq,char typ)
		{ return lpr->GetOptPttrnLPR(f,S1,S2,FALSE,L,x,csq,typ); }
	sst_typ *GetOptPttrnLPR(FILE *f,set_typ S1, set_typ S2,double &L,Int4 x,char typ,e_type qE=0,double pRi=0)
		{ return lpr->GetOptPttrnLPR(f,S1,S2,FALSE,L,x,typ,qE,pRi); }
	double  CalcSetvsPttrnLPR(FILE *fp,set_typ setFG, set_typ setBG,sst_typ *qsst,char typ,e_type qE=0)
		{ return lpr->CalcSetvsPttrnLPR(fp,setFG, setBG,qsst,FALSE,typ,qE); }
	double	GetPriorRi(Int4 n, char mode='C'){
		    double  sq_prior=0;
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
	//----------------------------------------------------------------

	//=================== omc_tweak.cc ==========================
	Int4    AdjustSetIDs();
	set_typ MkSubTreeSet(Info *info,Int4 pI);
	e_type  FindSeedSeq(Info *info, Int4 pI, Int4 I,sst_typ *osst,FILE *efp=0);
	void    PartitionNode(set_typ *set, Int4 pI, Int4 I, e_type bstE, e_type csqE);
	void    PatternBasedPartition(set_typ *set, sst_typ *sst,Int4 pI, Int4 I, double alpha);
	//=================== omc_tweak.cc ==========================
	Int4	PrintPttrn(FILE *fp,sst_typ *osst){
		    Int4 n=0,k;
		    for(k=1; k <= mcs->RtnLengthMainCMSA(); k++){
			if(osst[k]){ n++; if(fp){ PutSST(fp,osst[k],AB); fprintf(fp,"%d ",k); } }
		    } if(fp) fprintf(fp,"\n"); return n;
		}
private:
	//=================== omc_hihmm.cc ==========================
	Int4    TestHiHMM( );
	Int4    TestHiHMM2( );
	smx_typ	*ProfileHiHMM(Int4 nn);
	//=================== omc_hihmm.cc ==========================
	//=================== omc_cmh.cc ==========================
	class cmh_typ {	// cross-conserved min-max heap.
public:
		cmh_typ(){ print_error("cmh_typ input error"); }
		cmh_typ(Int4 s){
		    size=s; mH=Mheap(size,4); 
		    NEW(FG,size+3,set_typ); NEW(BG,size+3,set_typ);
		    NEWP(SST,size+3,sst_typ); NEW(Num,size+3,Int4);
		}
		Int4 Insert(set_typ fg, set_typ bg, keytyp k, sst_typ *isst)
		{	// absorbs isst but copies fg and bg sets.
			Int4 x;
			if(FullMheap(mH)){
			   if(k <= MinKeyMheap(mH)) return 0;
			   x=DelMinMheap(mH); NilSet(FG[x]); NilSet(BG[x]); free(SST[x]);
			   FG[x]=BG[x]=0; SST[x]=0;
			} x=InsertMheap(k,mH);
			Num[x]=1; FG[x]=CopySet(fg); BG[x]=CopySet(bg); SST[x]=isst;
			return x;
		}
		Int4 DeleteMax(set_typ &fg, set_typ &bg, keytyp &k, sst_typ *&isst){
			k=MaxKeyMheap(mH);
			Int4 rtn,x=DelMaxMheap(mH);
			fg=FG[x]; FG[x]=0; bg=BG[x]; BG[x]=0; isst=SST[x]; SST[x]=0;
			rtn=Num[x]; Num[x]=0; return rtn;
		}
		BooLean	Empty(){ return EmptyMheap(mH); }
		BooLean	NumHits(Int4 x){ if(x > 0 && x <= size) return Num[x]; else return 0; }
		BooLean	Full(){ return FullMheap(mH); }
		sst_typ	*RtnSST(Int4 x){ if(x < 1 || x > size) return 0; else return SST[x]; }
		Int4	Size(){ return size; }
		Int4	InHeap(set_typ fg, set_typ bg){
			// If fg and bg sets are in the heap then return the id; else return 0.
			Int4 i,f=CardSet(fg),b=CardSet(bg);
			for(i=1; i <= size; i++){
			   if(FG[i] == 0) continue;
			   if(CardUnionSet(fg,FG[i]) == f && CardInterSet(fg,FG[i]) == f &&
				CardUnionSet(bg,BG[i])== b && CardInterSet(bg,BG[i]) == b){
				Num[i]++; return i;
			   }
			} return 0;
		}
		~cmh_typ(){
			while(!EmptyMheap(mH)){
			   Int4 x=DelMinMheap(mH); NilSet(FG[x]); NilSet(BG[x]); free(SST[x]);
			} free(FG); free(BG); free(SST); free(Num); NilMheap(mH);
		}
private:
		mh_type mH;
		Int4	size,*Num;
		set_typ	*FG,*BG;
		sst_typ	**SST;
	};
	Int4	XC_line,XC_size;
	BooLean	XCuseX,XCuseL;
	char	XC_mode;
	double	XC_wt;
	double  CrossScoreBILD(Int4 key,Int4 site,set_typ FGSet,set_typ BGSet,dms_typ *dms,FILE *ofp=stdout);
	double  ComputeCrossScore(Int4 key,Int4 site,set_typ FGSet,set_typ BGSet,cmh_typ *cmh);
	Int4	PutCrossConserved(Int4 key,cmh_typ *cmh);
	//=================== omc_cmh.cc ==========================
	//=================== omc_debug.cc ==========================
	BooLean CheckSSTvsCSQ(mcs_typ *xmcs);
	void    CompareInput(FILE *fp,set_typ *set,hpt_typ *hpt,cma_typ *SMA);
	BooLean CheckID4NewMCS(mcs_typ *xmcs, Int4 ID);
	BooLean	CheckBestCopy(mcs_typ *xmcs);
	BooLean	TheSame(mcs_typ *xmcs);
	Int4	RecordChange(char *title);
	Int4	TestSubRoutine(char mode);
	void	PutPatternOverlap(FILE *fp);
	void    PutLineage(FILE *fp,Int4 node);
	BooLean OverlappingLineages(Int4 i, Int4 j);
	void    RandomStart(Int4 NumRandNodes);
	void    FlattenHierarchy( );
	//=================== omc_debug.cc ==========================
	BooLean	flatten_hpt;
	Int4	random_start;
	char	test_mode;
	h_type	goodHG,badHG,testHG,triesHG;
	double	AddLeafLLR;
	FILE	*logfp;
	Int4	status,t_stat;
	time_t	log_time;
	Int4	LastNumNodes;
	// Int4	NumNewNodes( ){ Int4 n=Hpt->NumBPPS()-LastNumNodes; LastNumNodes=Hpt->NumBPPS(); return n; }

	//========================= global objects ====================================
	cma_typ TrueMainCMA;	// main sequence alignment without Random seqs.
	cma_typ MainCMA;	// main alignment with Random set included.
	Int4	NumRandom;
	Int4	SizeMainCMA,SizeTrueMainCMA;
	set_typ	*DisplaySet;	// seed sequence display sets (==0 for internal nodes).
	hsw_typ MainHSW;	// sequence weights with Random seqs.
	hsw_typ hsw;		// 
	BooLean	own_hsw;	// Should this be freed or does swt free it?
	swt_typ *swt;		// sequence weights.
	char	*infile;
	char	*chk_prefix;	// checkpoint source filename.
	a_type	AB;

	//******************** Sampling over a single subtree *********************
	set_typ	stable;		// set of identifiers for nodes that remain fixed.
	Int4	FocusNode;	// node and subtree to be sampled (relabeled).
	Int4	FocusSeq;	// sequence to be sampled from FocusNode.
	BooLean	SetUpFocusedSrch();	// 

	//******************** file management *********************
	Int4	Argc;
	char	**Argv;
	char	str[500];
	dh_type	FileIdHeap;
	char	*outfilename;
	mcs_typ	*mcs;		// conserved domain hierarchy:
	hpt_typ	*Hpt;		// FD-table 
	Int4	MaxNumNodes;	// Maximum number of rows and of columns in the FD-table.
	Int4	MaxNumNodesPlus;	// MaxNumNodes + 2 to accomodate internal nodes at end.
	set_typ	RandomSet;
	BooLean	Evolve;		// Should consensus sequences be allowed to evolve from the start?
	//******************** mcs_typ file management *********************
	Int4	Iteration;
	time_t	time1;
	BooLean	SaveSets;

	//=================== omc_junk.cc ==========================
	//                Currently unused routines.
	double	RaiseMultiBranches(FILE *fp,double Temp, Int4 iter);
	Int4    MultiMoveNodesUp(double Temp, set_typ pre_skip=0, double minLLR=0.0);
	mcs_typ *MultiMoveUp(set_typ BadKids,double MinDeltaLLR);
	e_type	SetToBest(set_typ set, Int4 n, mcs_typ *xmcs);
	void	SortHpt();
	Int4	MoveUpAndDown(Int4 SampledColumn, double MinDeltaLLR=20.0);
	hsi_typ	*ArrangeAsTree0(hsi_typ **hsi,Int4 num,set_typ *Set);
	Int4	ResurrectRejects0(double Temp=300.0);
	hsi_typ *OptimizeInsertion(set_typ ChildFG, set_typ ChildBG, hsi_typ *hsi);
	hsi_typ *Check4Bottom(Int4 depth, Int4 maxdepth);
	//=================== omc_junk.cc ==========================
};

#endif

