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

#if !defined (_CHN_TYP_)
#define	_CHN_TYP_

#include "gpsi_typ.h"
#include "rtf_typ.h"
#include "rst_typ.h"
#include "gpsi2cma.h"
#include "new-dsc.h"
#include "txs_typ.h"
#include "pdb.h"
#include "swt_typ.h"
#include "swaln.h"
#include "mheap.h"
#include "sset.h"
#include "pcr_typ.h"
#include "table.h"
#include "prtn_model.h"
#include "blosum62.h"

// #include "bpps_typ.h"
// #include "pps_typ.h"

#if 0	//********************************************************************
	  GaI vs Std                                        all_proteins (std)
	  GaI vs QueryMainSet (else GaI vs Not_GaI)              |
	                 (GaI distinct from Ga)                  |
                                                             (GTPases)
	  Ga vs Not_Ga (Ga distinct features)                /   |    \
                                                          (Ga) (Ras) (Efs)
	  Ras vs Not_Ras (Ras distinct features)         /  | \
                                                       GaI GaS GaQ
	  EFs vs Not_EFs (EF distinct features)

	  Not_GTPases vs Std (GTPase distinct features)

	Store the following *.cma files in *.chn:
		1. Query_subfamily		{ GaI (txs) }		
		2. alignment consensus seqs	{ for alignments below }
		3. txs subalignments
		4. Query Famliy1		{ Ga family }		
		5. Not 4.			{ Not Ga family }
		6. First major cluster		{ Ras
		7. Not 6.			{ Not Ras-like }
		8. Second major cluster		{ EFs
		9. Not 8.			{ Not EFs }
		2. Superfamily 			{ GTPases (all) }

alignment consensus seqs:		cmafile_code
>QuerySubFamily 			 'Q'
 	:	:
>Family1 				 'F'
 	:	:
>Family2 				 'F'
 	:	:
>BGFamily2 				 'B'
 	:	:
 	:	:
>Cluster1 				 'C'
	:	:
>BGCluster1 				 'B'
	:	:
>Cluster2 				 'C'
 	:	:
>BGCluster2 				 'B'
 	:	:
 	:	:
>SuperFamily  				 'S'
 	:	:
#endif	//********************************************************************

#define MAX_ALN_CHN_TYP	100
#define MAX_ALN_CHN_TYP_PLUS	(MAX_ALN_CHN_TYP + 2)
#define MAX_PDB_CHN_TYP	20
#define MAX_CHN_PATTERNS 100

class chn_typ { 	// PSI-BLAST Compare class
public:
		chn_typ( ){ assert(!"Illegal constructor"); }
		chn_typ(Int4 argc,char *argv[]) {
			   cmc_rtf_mode=FALSE; cmc_res_vals=0;
			   cmc_input=FALSE;
			   Init(argc,argv,0,0,101.0,0,0); 
			   if(!isprint(output_name[0])){
				fprintf(stderr,"Possible library corruption problem\n");
				// fprintf(stderr,"**** output_name(1)=%s ****\n",output_name);
			   	assert(isprint(output_name[0])); assert(!OutputAln);
			   }
			}
		// New constructer for calling from omc_typ or cmc_typ. 
		chn_typ(Int4 argc,char *argv[],Int4 NumCMA,cma_typ *CMA_ARRAY,hsw_typ *HSW,
				char **in_status,double default_percent_cut,
				double ***bpps_res_evals){
			cmc_rtf_mode=TRUE; cmc_rtf_status=in_status;
			cmc_input=FALSE;
			cmc_res_vals=bpps_res_evals;

			Init(argc,argv,NumCMA,CMA_ARRAY,default_percent_cut,0,HSW); 
			assert(isprint(output_name[0])); assert(!OutputAln);
			}
#if 0		// Can delete this one later after switching to above.
		chn_typ(Int4 argc,char *argv[],Int4 NumCMA,cma_typ *CMA_ARRAY,hsw_typ *HSW,
				char **in_status,double default_percent_cut){
			cmc_rtf_mode=TRUE; cmc_rtf_status=in_status;
			cmc_input=FALSE;
			cmc_res_vals=0;
			Init(argc,argv,NumCMA,CMA_ARRAY,default_percent_cut,0,HSW); 
			assert(isprint(output_name[0])); assert(!OutputAln);
			}
#endif
		chn_typ(Int4 argc,char *argv[],Int4 NumCMA,cma_typ *CMA_ARRAY,
				double default_percent_cut,hsw_typ *HSW){
			cmc_rtf_mode=FALSE; cmc_res_vals=0;
			cmc_input=TRUE;
			Init(argc,argv,NumCMA,CMA_ARRAY,default_percent_cut,0,HSW); 
			assert(isprint(output_name[0])); assert(!OutputAln);
			}
		chn_typ(Int4 argc,char *argv[],Int4 NumCMA,cma_typ *CMA_ARRAY,
				double default_percent_cut,char *NameSWT){
			cmc_rtf_mode=FALSE; cmc_res_vals=0;
			cmc_input=TRUE;
			Init(argc,argv,NumCMA,CMA_ARRAY,default_percent_cut,NameSWT,0); 
			assert(isprint(output_name[0])); assert(!OutputAln);
			}
		chn_typ(Int4 argc,char *argv[],Int4 NumCMA,cma_typ *CMA_ARRAY,
				double default_percent_cut){
			cmc_rtf_mode=FALSE; cmc_res_vals=0;
			cmc_input=TRUE;
			Init(argc,argv,NumCMA,CMA_ARRAY,default_percent_cut,0,0); 
			assert(isprint(output_name[0])); assert(!OutputAln);
			}
		chn_typ(Int4 argc,char *argv[],double default_percent_cut) {
			cmc_rtf_mode=FALSE; cmc_res_vals=0;
			cmc_input=FALSE;
			Init(argc,argv,0,0,default_percent_cut,0,0); 
			assert(isprint(output_name[0])); assert(!OutputAln);
			}
		~chn_typ( ){ Free(); }
	void    PutHierarchicalAlignment(){
			// fprintf(stderr,"**** output_name(1)=%s ****\n",output_name);
			PutHierarchicalAln(output_name); 
		}
	void    PutHierarchicalAln(char *filename);
	void    PutHierarchicalAlns(FILE *fp){ PutCHA(fp,ColorString); }
	void    PutHierarchAligns(FILE *fp){ PutCHA(fp,0); }
	void    PrintFileName(FILE *fp){ fprintf(fp,"{\\b\\f2\\fs%d\\cf1 %s }{\\f3 \n\\par }\n",
				fontsize+4,output_name); }
	void    PutPageBreak(FILE *fp){ fprintf(fp,"\n\\par \\page \n"); }
	void    PutHierarchicalAlnHead(FILE *fp){ PutCHA_HEADER(fp); }
	void    PutHierarchicalAlnTail(FILE *fp){ PutCHA_TAIL(fp); }
	BooLean	OutPutAln(){
		  // assert(isprint(output_name[0])); assert(!OutputAln);
		  if(OutputAln) return TRUE; else return FALSE; 
		}
	void    PutAln(char *filename,char mode);
	void    PutAln(){ PutAln(output_name,' '); }
	void    PutAln(char *filename){ PutAln(filename,' '); }
	// char	*FileName( ){ return AllocString(output_name); }
	char	*FileName( ){ return output_name; }
	Int4	GetNumInCMSA(){ return Number; }
	Int4	GetNumAnalysis(){ return NumAnalysis; }
	a_type	GetAlphabet( ) { return AB; }
	e_type	TrueFirstSeq( ) { 
			return (TrueSeqCMSA(1,IN_CMA[1]));
			// return CopySeq(TrueSeqCMSA(1,IN_CMA[1]));
		}
	double	**GetNullFreq(Int4 x)
		{ x = NumberMCMA-x+1; 
		  if(x > 0 && x <= NumberMCMA) return NullFreq[x]; else return 0; }
	hsw_typ	RtnHSW(Int4 x){
		  if(x > 0 && x <= NumAnalysis){ assert(swt[x]); return swt[x]->RtnHSW(); }
		  else return 0;
		}
	BooLean	FWriteSWT(Int4 x,char *FileName){
		  if(x > 0 && x <= NumAnalysis){ swt[x]->FWrite(FileName); return TRUE; } else return FALSE;
		}
	UInt4   GetWtFactor(){ return swt[1]->WtFactor(); }  // all should use the same WtFactor...
	UInt4   GetWtFactor(Int4 x){ 
		  if(x > 0 && x <= NumAnalysis){ return swt[x]->WtFactor(); } else return 0;
		}
	UInt4	*GetIntWeights(Int4 x,unsigned char ***RtnSqWt){
		  if(x > 0 && x <= NumAnalysis){
		    return swt[x]->GetIntegerWts(RtnSqWt);
		  } else return 0;
		}
	double	**GetWeights(Int4 x)
		{ // x = NumAnalysis-x+1;  // not reversed???
		  // if(x > 0 && x <= NumberMCMA) return swt[x]->Weights(); else return 0; }
		  if(x > 0 && x <= NumAnalysis) return swt[x]->Weights(); else return 0; }
	double	**GetWtFreqs(Int4 x)
		{ x = NumAnalysis-x+1; 
		  // if(x > 0 && x <= NumberMCMA) return WtFrq[x]; else return 0; }
		  if(x > 0 && x <= NumAnalysis) return Freq[x]; else return 0; }
	cma_typ	GetMCMA(Int4 x)
		{ if(x > 0 && x <= NumAnalysis) return MCMA[x]; else return 0; }
	cma_typ	*GetIN_CMSA( ) { return IN_CMA; }
	e_type	QuerySeq( ) { return qE; }
	e_type	KeySeq( ) { return keyE; }
	cma_typ	GetIN_CMSA(Int4 x)
		{ if(x > 0 && x <= Number) return IN_CMA[x]; else return 0; }
	rtf_typ	*GetRtfQ(Int4 x) {
		x = NumAnalysis-x+1; if(x > 0 && x < MAX_ALN_CHN_TYP_PLUS) return rtfQ[x];
		else return 0; }
	rtf_typ	*GetRtf(Int4 x) {
		x = NumAnalysis-x+1; if(x > 0 && x < MAX_ALN_CHN_TYP_PLUS) return rtf[x];
		else return 0; }
	pcr_typ	*get_pcr( ) { pcr_typ *tmp=pcr; pcr = NULL; return tmp; }
	void	PutVSI( ) { pcr->print(stdout,0); }
	void	PutVSI(FILE *fp,Int4 offset) {
			pcr->print(fp,offset);
		}
	void	PutVSI(FILE *fp) { pcr->print(fp,0); }
	BooLean	OutputVSI(){ return output_vsi; }
	BooLean	OwnsCMAs() { return OwnCMAs; }
	void	AddXconserved(double *XC){ if(Xconserved) free(Xconserved); Xconserved=XC; }
private:
	double	*Xconserved;	// values for cross conservation... Xconserved[0]=maximum.
	void	Init(Int4 ,char *argv[],Int4 ,cma_typ *,double,char *,hsw_typ*);
	void	InitAsNull();
	char	ModeLPR;
	double	ExpPatterns;

	void    CreateRTF( );	// in chn_rtf.cc
	BooLean	RTFsCreated;
#if 0	// changed for borrowing the input CMA files
	void    ReadChnAnalMultAln(Int4 argc, char *argv[],const char *USAGE); // in chn_read.cc
#else
	void	ReadChnAnalMultAln(Int4 argc, char *argv[],Int4 NumCMA,cma_typ *CMA_ARRAY,const char *USAGE);
	BooLean	OwnCMAs;
#endif
	// **************** End additional routines: *******************

	void    put_marg_prob(FILE *fptr, Int4 gstart, Int4 start, Int4 gend,
        		Int4 color,Int4 gapsize,char *gnull,Int4 Analysis);
	void    PutCHA(FILE *fptr, const char *ColorFont);
	void    PutCHA_RTNS(FILE *fptr,unsigned short numRtns);
	void    PutCHA_TAIL(FILE *fptr);
	void    PutCHA_HEADER(FILE *fptr);
	char    FractionToCharPBC(double fract);
	void    Free( );                // free memory...
	cma_typ *CMA,*IN_CMA;
	pcr_typ	*pcr;
	a_type	AB,ownAB;
	e_type	qE,keyE,full_csq;
	double	Alpha;
	Int4	A0,B0;
	double	*Alphas;	// for multiple input files...
	Int4	*A0s,*B0s;
	char	*SetMode;

	// Alignment parameter options:
	Int4	char_per_line,fontsize;
        double  cbp_cut,percent_cut,mp_cutoff,fsa_cutoff;
        BooLean verbose,NoWeights,use_ocma_pseudo,use_sfbg,extra_files;
        char    PageSetUp;
        double  FractHighlight;  
        double  FractHighlightTOP;
	Int4	InputFractionAsNumber;
	BooLean	InputFraction,show_marg_prob;
        BooLean OutputAln;	// output simple alignment for a cma file.
	BooLean	compute_mpwf,InvertOrder,HideIndelInfo;
        float   mpwf_scale;
        Int4    minbars;
        Int4    Begin, End;	// regions in query to align...
	BooLean	UseColPos;
	BooLean	SelectColPos;
	Int4	Seq4numbering;	// which sequence to use for numbering?
	Int4	maxlen_gnull;
#if 1	// InsertLen = 0...INT4_MAX
	Int4	*gnull_insrt_len;	// length of insertion...
#endif

	// PSI-BC input file information:
        Int4    Number,NumSeqAln,NumberMCMA;

	// PDB input file information:
	// char	checkin[200]
	char	temp_name[300];
	// char	output_name[300];
	char	*output_name;

	// Data storage for analysis.
	sma_typ	MA;
	double	*MargProb[MAX_ALN_CHN_TYP_PLUS],**NullFreq[MAX_ALN_CHN_TYP_PLUS];
	double	*FractSeqAln[MAX_ALN_CHN_TYP_PLUS],**Obs[MAX_ALN_CHN_TYP_PLUS];
	rtf_typ	*rtf[MAX_ALN_CHN_TYP_PLUS],*rtfQ[MAX_ALN_CHN_TYP_PLUS];
	char   	*SecondStrct[MAX_ALN_CHN_TYP_PLUS];
	unsigned char merged_rasmol[MAX_ALN_CHN_TYP_PLUS];
	cma_typ	MCMA[MAX_ALN_CHN_TYP_PLUS],*TAX_CMA;
	char	ColorCode[MAX_ALN_CHN_TYP_PLUS];
	char	*ColorString;
	char	AlignCode[MAX_ALN_CHN_TYP_PLUS];
	double	**Freq[MAX_ALN_CHN_TYP_PLUS];
	Int4	NumSqMain[MAX_ALN_CHN_TYP_PLUS];
	Int4	NumAnalysis,BackGrnd[MAX_ALN_CHN_TYP_PLUS];
	BooLean Hist[MAX_ALN_CHN_TYP_PLUS],SuperAln[MAX_ALN_CHN_TYP_PLUS];
	BooLean	PatternOnlyMode,StdAlignmentOnly,cha_as_input,cmc_input,cmc_rtf_mode;
	double	***cmc_res_vals;

	// Input pattern info...
	BooLean InvertResSet[MAX_CHN_PATTERNS];
	char    Residue[MAX_CHN_PATTERNS];
	Int4	NumResidues;
	Int4	MaxMisMatches,MaxConcensusLines;
	char    residue_str[MAX_CHN_PATTERNS][30];
	sst_typ	Residues[MAX_CHN_PATTERNS];
	Int4    Position[MAX_CHN_PATTERNS],Position0[MAX_CHN_PATTERNS];

	char    *Kingdom;	// Kingdoms for display sequences...
	char    **Phylum;	// Phylum for display sequences...
	char	**status,**cmc_rtf_status;
	BooLean	output_vsi;
	Int4	insert,extend;

	// NEW private data to clean up code...(need to test these...)
	BooLean Show2ndary,RmFamilyAln;
	void	GetArg(Int4 argc,char *argv[]);

	// for rtf_typ
	double	MinResFreqCutoff;

	// BG_CMA file structures:	This is a block-based motif file...
	cma_typ	BG_CMA;
	BooLean	use_bg_cma;
	Int4	StartBG_CMA,EndBG_CMA;
	double	**NullFreqBG,**ObsBG,*FractSeqAlnBG,**WtFrqBG;
        double  FractHighlightBG;
	char	sets_mode,root_sets_mode;
	void	GetSwtBG(Int4 anal);

	// csp_typ with Integer Weights!  // moved here from above...
	void    CheckWeightedMatrix(e_type qE, double  **WtFreq, a_type A);
	double	*GetWtNumSq(double *RtnWtNSq,double **Obsvd, Int4 length);

	swt_typ **swt;
	double	**WtFrq[MAX_ALN_CHN_TYP_PLUS];
	UInt4	**ObsIWt[MAX_ALN_CHN_TYP_PLUS];	// IntegerWeighted counts...
	double	WtNumSq[MAX_ALN_CHN_TYP_PLUS];
	double	*WtNumSeq[MAX_ALN_CHN_TYP_PLUS];

	// void	GetSWT(char *);
	void	GetSWT(hsw_typ*);
	void	WriteMainSWT();
	void	GetMergedIntegerSWT(cma_typ qcma, cma_typ fg_cma, cma_typ bg_cma);
	unsigned char **MainSqWt;
	cma_typ	MainCMA;
	UInt4	*SqIWtMain;
	UInt4	TotalFG_N,TotalBG_M;
	unsigned char StartAlpha;
	BooLean	GlobalSqWts;
#if 0
	void	GetIntWeights(Int4 x,);	// for analysis x..
#endif
};

#endif

