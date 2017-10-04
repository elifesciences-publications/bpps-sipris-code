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

#if !defined (_PSC_TYP_)
#define _PSC_TYP_

#if 0
#include "chn_aln.h"
#include "chn_pps.h"
#include "chn_typ.h"
#include "chn_vsi.h"
#include "tax_typ.h"
#endif
#include "swaln.h"
#include "adh_typ.h"
#include "histogram.h"
#include "dheap.h"
#include "random.h"
#include "residues.h"
#include "blosum62.h"
#include "cmsa.h"
#include "pdb.h"
#include "atom.h"
#include "sequence.h"
#include "mheap.h"
#include "rai_typ.h"
#include "rih_typ.h"
#include "clique.h"
#include "hpt_typ.h"
#include "set_typ.h"
#include "pat_typ.h"
#include "sset.h"
#include "mps_typ.h"
#include "esc_typ.h"
#include "swaln.h"

// #define MAX_NUM_PDB_INPUT_FILES 1000
// #define	MAX_NUMBER_INTERNAL_RPTS 100

#define	USE_COL_PAIR_FILTER 1

// call this from hsc_typ passing in IN_CMA and mpdb objects.
// pass this into psc_typ from hsc_typ instead of mpdb + IN_CMA both...
class p2c_typ {		// pdb to cma mapping type.
public:
		p2c_typ(){ print_error("this constructor disallowed"); }
		p2c_typ(Int4 argc, char *argv[],mps_typ *in_mpdb, cma_typ *in_cma,
			   Int4 number,char *out_file,esc_typ *in_esc,hpt_typ *hpt0,
			   pat_typ *ptrn0, char *trace_colors,char *side_colors,
			   set_typ **cps=0){
#if USE_COL_PAIR_FILTER	// For sas aln to aln compare only!!
			ColPairSet=cps; MapPdb2cmaSq = 0;
#else
			ColPairSet=0; MapPdb2cmaSq = 0;
#endif
			TraceColors=trace_colors;
			SideColors=side_colors;
			hpt=hpt0; ptrn=ptrn0; // ptrn->Put(stdout);
			mpdb=in_mpdb; IN_CMA=in_cma; Number=number;
			AB=AlphabetCMSA(IN_CMA[1]);
			esc=in_esc; OutFile=out_file; 
			Diagnostic=0;
			Init();
			GetArg(argc,argv);
			FullSeq=esc->FullSeq;
			MapSqAlnToStruct( );
		}
		~p2c_typ() {
			if(MapPdb2cmaSq){
			  for(Int4 I=1; I <= mpdb->NumPDB; I++){
			  	if(MapPdb2cmaSq[I]) free(MapPdb2cmaSq[I]);
			  } free(MapPdb2cmaSq);
			}
			Free(); 
		}
	void    GetStructAlnScores(FILE *scrfp, Int4 Row,BooLean PerformControl)
			{ this->RoundTwo(scrfp, Row,PerformControl); return; }
	void	PrintVSI_Files(FILE *vsifp=0);
	int	PrintKLST_Files(char call=0);
	mps_typ	*mpdb;
	esc_typ	*esc;
	char	*OutFile;
	set_typ	*RelevantSet;	    // A element of RelevantSet[S] = pdb seq matching a rpt in cma file A.
	Int4	***RelevantSeq;	    // RelevantSeq[I][C][R] == A set.
	Int4	****Col2pdbSeq;  // Col2pdbSeq[I][C][R][column] = pdb residue position i 
	Int4	**NumRpts;		// NumRpts[I][C] = number of repeats....
	Int4	NumCol( ){ return LengthCMSA(1,IN_CMA[1]); }

	e_type	*FullSeq;	// FullSeq[S] == full sequence for equivalence set S.
	Int4	***Col2FullSeq;	// Col2FullSeq[S][R][col] = seq, repeat, column.
	Int4	**RptCategory;	// RptCategory[S][R] = A
	Int4	*NumFullRpts;	// NumFullRpts[S] 

	double	LastARSD;
	double	FinalARSD;
	Int4	FinalPnts;
	Int4	FinalCols;
	// void	AssignSets(set_typ **sets){ ColPairSet=sets; }

	set_typ	**ColPairSet;
	Int4	**MapPdb2cmaSq;	// map pdb full seq to corresponding CMA seq.

	Int4	LastNumData,LastTotalData;

	// search routines
	void    FindBestPatternContacts2();
	void	ShowTheBest(){ FindBest=TRUE; }
	hpt_typ *hpt;
#if 1	// For using top scores only...assessing structural scores...
	//========================= p2c_CaScrs.cc =====================================
	adv_typ **RunAndRtnADV(FILE *scrfp,Int4 &NumADV);
	void    PutAdvList(FILE *fp, Int4 NumADV, adv_typ **ADV,Int4 NumShown=0,char *msg=0);
	adh_typ *RoundOne(FILE *scrfp, Int4 Row, Int4 OmitCol=-1, Int4 OmitSeq=-1);
	void    RoundTwo(FILE *scrfp, Int4 Row,BooLean PerformControl);
	adv_typ **GetAdvList(Int4 &Num,adh_typ *adh);
	char	*Diagnostic;	// String for diagnosing problems.
	//========================= p2c_CaScrs.cc =====================================
	void	SetTarget(Int4 t){ assert(t >= 5); Target=t; }
	void	SetDataPoints(Int4 t){ assert(t > 9); MaxDataPoints=t; }
#endif
	Int4	NumPDB( ){ return mpdb->NumPDB; }
	e_type	RtnSqPDB(Int4 i){ assert(i > 0 && i <= mpdb->NumPDB); return mpdb->pdbSeq[i][1]; }
	char	*RtnPDB_ID(Int4 i){ assert(i > 0 && i <= mpdb->NumPDB); return mpdb->pdb_file[i]; }
private:
	Int4	Target;	// # seqs to include in histogram.
	Int4	MaxDataPoints;
	char	*TraceColors,*SideColors;
	void	Free(){
		   Int4 S,R,I,C;
		   if(Diagnostic) free(Diagnostic);
		   for(S=1; S <=esc->NumPDB_Sets; S++){
			for(R=1; R <= NumFullRpts[S]; R++){ free(Col2FullSeq[S][R]); }
			NilSet(RelevantSet[S]); free(RptCategory[S]); free(Col2FullSeq[S]);
		   } free(RelevantSet); free(RptCategory); free(Col2FullSeq); free(NumFullRpts);
		   for(I=1; I <=mpdb->NumPDB; I++){
			for(C=1; C <=nChainsPDB(mpdb->pdb[I]); C++){
			   for(R=1; R <= NumRpts[I][C]; R++){
			      if(Col2pdbSeq[I][C][R]) free(Col2pdbSeq[I][C][R]);
			   } if(RelevantSeq[I][C]) free(RelevantSeq[I][C]);
			   if(Col2pdbSeq[I][C]) free(Col2pdbSeq[I][C]);
			} free(RelevantSeq[I]); free(Col2pdbSeq[I]); free(NumRpts[I]);
		   } free(RelevantSeq); free(Col2pdbSeq); free(NumRpts);
		}
	//========================= p2c_debug.cc =====================================
	void    DebugMapSqAln2Strct(FILE *fp,Int4 S,Int4 C,Int4 R,Int4 I,Int4 *TmpColToSeq,
					e_type pdbS, e_type pdbIC, e_type pdbE);
	//========================= p2c_debug.cc =====================================
	//========================= p2c_junk.cc =====================================
	//========================= p2c_junk.cc =====================================
	//========================= p2c_typ.cc =====================================
	Int4	FindHBonds(FILE *fp,res_typ Res,float dmax,res_typ Res2,
                	char color, unsigned short file_id);

	// output routines
	Int4    FindVDWsContacts(FILE *fptr,Int4 file_id,Int4 ResI,Int4 ResJ,float HA_dmax,float dmax,
                Int4 C,char chain,char color,res_typ **ResALL_I,Int4 *num_resALL_I,pdb_typ P);
	Int4    FindResResHBonds(FILE *fptr, Int4 file_id, Int4 ResI,Int4 ResJ,
		float HA_dmax,float dmax,Int4 C, char chain, char color,
		res_typ **ResALL_I, Int4 *num_resALL_I, pdb_typ P);

	// search for the best pairwise interactions
	Int4    ResToHeteroHBonds(FILE *fp,res_typ Res,float dmax,char color,
			unsigned short file_id, pdb_typ P);
	Int4    HeteroToResHBonds(FILE *fp,res_typ ResA,float dmax, char color,
                        unsigned short file_id, res_typ **ResALL_I,Int4 *num_resALL_I, pdb_typ P);

	pat_typ *ptrn;
	char    **AdjacentHeteroMolecule(Int4 RR, Int4 CC, Int4 II, double maxdist);
	BooLean	AddColToSeq(Int4 S, Int4 A, Int4 *ColToSeq);
	BooLean	AddColToSeq(Int4 i, Int4 S, Int4 A, Int4 *ColToSeq);
	a_type	AB;
	void	Init();
	void	GetArg(Int4, char**);
	void	MapSqAlnToStruct( );
	cma_typ	*IN_CMA;	// multiple 
	cma_typ cma;
	Int4	Number;

	//================ for psc_score.cc ================
	float	HA_dmax, dmax;
	// FILE    *ofp;
        BooLean FindBest;
	BooLean	UseBeta;
        Int4    K1,K2;
        double	MaxMeanDist;
	double	MinCaDist,MaxCaDist;
        Int4    MaxSqDist;
        Int4    MinDistInSeq;
        Int4    MinDist;      // what about inserts in some proteins?
        Int4    begin,end;
        Int4    Begin,End;
        Int4    KeyCol;
        double  MinVar;
        double  bin_size;
        Int4    D;              // Separation in sequence.
};

class psc_typ {         // protein structural comparison type
public:
		psc_typ(){ print_error("this constructor disallowed"); }
		psc_typ(p2c_typ *p2c, Int4 cma_file_id, char *cmaInFileName,
						pat_typ *in_ptrn,esc_typ *esc){
			initialize( );
			cmaFileID=cma_file_id; 
			ptrn=in_ptrn; ptrn0=0;

			mpdb=p2c->mpdb; mpdb0=0;
			RelevantSeq=p2c->RelevantSeq;
			RelevantSet=p2c->RelevantSet;
			NumRpts=p2c->NumRpts;
			Col2pdbSeq=p2c->Col2pdbSeq;
			num_cols = p2c->NumCol( );

			FullSeq=p2c->FullSeq;
			Col2FullSeq=p2c->Col2FullSeq;
			RptCategory=p2c->RptCategory;
			NumFullRpts=p2c->NumFullRpts;

			p2c0=p2c;

			cmafilename = AllocString(cmaInFileName);

			// CreateDsets( );
			PDB_SetI=esc->PDB_SetI;
			PDB_SetC=esc->PDB_SetC;
			NumPDB_Set=esc->NumPDB_Set;
			NumPDB_Sets=esc->NumPDB_Sets;
        		FullSeq=esc->FullSeq;	// FullSeq[S].
			FindMatches( );
		}
		~psc_typ( ){ Free(); }
	void    FindBestPatternContacts();
	void    FindBestPatternContacts2(char mode=0,FILE *logfp=0); // new version...
	p2c_typ	*p2c0;
	void    VerboseOff(){ verbose=FALSE; }
        void    VerboseOn(){ verbose=TRUE; }
private:
	BooLean	verbose;
	char	*GetRowName(Int4 x)
		{
			Int4    col,row,*P;
			assert(x > 0 && x <= ptrn->NumPttrns);
			if(p2c0->hpt->IsTree(P) == FALSE){ free(P); return 0; }
			free(P); col=ptrn->PttrnCategory[x];
                	row= ptrn->NumPttrns - ptrn->PttrnCategory[col] + 1; // reversed order for *.pttrns file.
                	assert(row > 0 && row <= p2c0->hpt->NumBPPS());
                	return p2c0->hpt->ElmntSetName(row);
		}
	e_type	*FullSeq;	// FullSeq[S] == full sequence for equivalence set S.
	Int4	***Col2FullSeq;	// Col2FullSeq[S][R][col]
	Int4	**RptCategory;	// RptCategory[S][R] = A
	Int4	*NumFullRpts;	// NumFullRpts[S] 

	void	FindMatches( );
	Int4	num_cols;

	Int4	NumWeakHBonds;
	Int4	***WeakHBonds;	// WeakHBonds[I][C][R]
	Int4	NumModHBonds;
	Int4	***ModHBonds;	// ModHBonds[I][C][R]
	Int4	***PatternResHbonds;
	// PatternResHbonds[I][C][R] = 54. (54 pattern residue interactions for this chain rpt.)

	// clustering of pdb files into equivalence sets.
	// void	CreateDsets();	// Disjoint sets.
	Int4	**PDB_SetI;	// PDB_SetI[S][i] = I; 
	Int4	**PDB_SetC;	// PDB_SetC[S][i] = C; 
	Int4	*NumPDB_Set;	// NumPDB_Set[S] = number in each set.
	Int4	NumPDB_Sets;	// NumPDB_Sets = number of sets.

	// MapSqAlnToStruct( );
	// void	MapSqAlnToStruct( );
	set_typ	*RelevantSet;	    // A element of RelevantSet[S] = pdb seq matching a rpt in cma file A.
	Int4	***RelevantSeq;	    // RelevantSeq[I][C][R] == A set.
	Int4	****Col2pdbSeq;  // Col2pdbSeq[I][C][R][column] = pdb residue position i 
	BooLean	****matches; 	// if(matches[I][C][R][j]==A) then pdb[I][C] matches a set A pattern.
	Int4	**NumRpts;		// NumRpts[I][C] = number of repeats....
	Int4	cmaFileID;
	void	initialize( );
	Int4    FindHBonds(FILE *fp,res_typ Res,Int4 C,Int4 C2,float dmax,Int4 Res2,
                		char color, unsigned short file_id, pdb_typ P);
	Int4    FindHBonds(FILE *fp,res_typ Res,float dmax,res_typ Res2,
                char color, unsigned short file_id);
	void	GetArg(Int4 argc, char *argv[]);
	void	Free();

	// output routines
	Int4    FindVDWsContacts(FILE *fptr,Int4 file_id,Int4 ResI,Int4 ResJ,float HA_dmax,float dmax,
                Int4 C,char chain,char color,res_typ **ResALL_I,Int4 *num_resALL_I,pdb_typ P);
	Int4    FindResResHBonds(FILE *fptr, Int4 file_id, Int4 ResI,Int4 ResJ,
		float HA_dmax,float dmax,Int4 C, char chain, char color,
		res_typ **ResALL_I, Int4 *num_resALL_I, pdb_typ P);

	// global variables
	a_type  AB,nAB;
	time_t  time1;

	//================ for psc_score.cc ================
	float	HA_dmax, dmax;
	

	// pdb coordinates
	mps_typ *mpdb,*mpdb0;
        Int4    NumPDB_Chains;

	// subgroup patterns corresponding to the input pdb sequence alignment.
	pat_typ *ptrn0,*ptrn;
	char	*pttrn_filename;

	// search for the best pairwise interactions
	char    **AdjacentHeteroMolecule(Int4 RR, Int4 CC, Int4 II, double maxdist);
	Int4    ResToHeteroHBonds(FILE *fp,res_typ Res,float dmax,char color,
			unsigned short file_id, pdb_typ P);
	Int4    HeteroToResHBonds(FILE *fp,res_typ ResA,float dmax, char color,
                        unsigned short file_id, res_typ **ResALL_I,Int4 *num_resALL_I, pdb_typ P);

	// cma multiple alignment file
	char	*cmafilename;
};

#endif

