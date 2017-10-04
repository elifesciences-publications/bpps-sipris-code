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

#if !defined (_MGS_TYP_)
#define _MGS_TYP_

#include "hat_typ.h"
#include "dft_typ.h"
#include "stdinc.h"
#include "stack.h"
// #include "stk_typ.h"
#include "fdt_typ.h"

class mgs_typ {         // MAPGAPS search type 
public:
                mgs_typ( ){ assert(!"Illegal constructor"); }
                mgs_typ(int argc,char *argv[],const char *usage) { Init(); GetArg(argc,argv,usage); }
                mgs_typ(int argc,char *argv[]){ Init(); GetArg(argc,argv); }
	// BooLean	Compare(Int4 Ni,char **iR,char **iName,Int4 No, char **oR,char **oName,char *TypeOfSet);
	BooLean	Compare(Int4 Ni,char **iR,char **iName,char *TypeOfSetI,
						Int4 No, char **oR,char **oName,char *TypeOfSetO);
	BooLean	CompareCMSAs(Int4 Ni,char **iR,char **iName,char *TypeOfSetI,
						Int4 No, char **oR,char **oName,char *TypeOfSetO);
	BooLean	Run(BooLean ps=TRUE, FILE *ofp=0){
        	 SetUpSearch();
		 if(Search( )) {
        		ProcessHits( );
			// fprintf(stderr,"DEBUG3\n");
        		ConvertViaTemplate( ); // problem here with deletions on ends of template seq.
			// fprintf(stderr,"DEBUG4\n");
        		if(OutputMAPS) PutMAPInputFile( );
        		PutSeqs(ps); PutAln( ); PutSummary( ); PutMainAln(ofp);
        		CleanUpSearch();
			return TRUE;
		 } return FALSE;
		}
                ~mgs_typ( ){ Free(); }
private:
	void    StoreSummary(Int4 R);
	void    StoreSummaryCMSAs(Int4 R);
	Int4    Intersection(ss_type P1, ss_type P2);
	BooLean CheckInput(Int4 N_cma, Int4 N_rows, cma_typ *CMSA,char **Name_row);
	fdt_typ *FDT;
	dft_typ	*dft;
	Int4	*Level;
        ss_type Data;
	BooLean	OwnData;
	// void	SetUpSearch(){ SetUpSearch(0); }
	void	SetUpSearch(ss_type in_data=0);
	// BooLean	Search( ) { return Search(stderr); }
	BooLean	Search(FILE *f=stderr);
	BooLean	TreeBasedSearch(Int4 *Level);
	Int4    *RtnTree( );
	void	ProcessHits( );
	void	PutSeqs(BooLean ps=TRUE);
	void	ConvertViaTemplate( );
	void	PutAln( );
	void	PutSummary( );
	void	PutMainAln(FILE *ofp=0);
	void    PutMAPInputFile( );
	void    CleanUpSearch();
	void	Init();
	void	GetArg(int argc, char *argv[]);
	void	GetArg(int argc, char *argv[],const char *usage);
	void	Free();

	// ************** files: ****************
        char    *Checkpoint,*infile,*database;
        char    *FalseCheckPoint;

	// ************** objects & arrays: ****************
	a_type	AB;
        ss_type FalseData;
        e_type  *QueryE;
	cma_typ TemplateCMA,*SRCH_CMA,*ALN_CMA;
        sap_typ *SAP;
        BooLean	*IGNORE;
        char    str[900];
	set_typ	UnionSetH,FalseSetH,TrueSetH;

	// ************** psi-blast parameters: ****************
        Int4    T,gap_open,gap_extend;
        Int4    left_flank,right_flank;
        double  gap_trigger_bits;
        double  x_parameter;
        Int4    width,printopt;
        UInt4   hpsz;
        double  Ethresh,Ecutoff;
        char    mask_dbs;
        Int4    dbslen;
        double  misalncut;
        BooLean see_smatrix,seg_seq;

	// ************** database parameters: ****************
        UInt4   max_length_input_seqs;

	// ************** MAPGAPS parameters: ****************
	BooLean	use_patch_sap,found_false,found_seqs;
        BooLean output_by_class,out_all_class;
        BooLean print_names_only,find_repeats;
	BooLean	DoGlobalAlign;
        Int4    num_cma_files,overlap_cutoff;
	char	prune_mode;
        Int4    num_false_pos,first_false_pos;
        Int4    Begin;        // Begin with template (II=0) or with first profile (II=1).
        BooLean OutputMAPS;
	Int4	SeqsTrue,SeqsFalse;
	Int4	HSPsTrue,HSPsFalse;
        BooLean use_weak_cutoff; // find weakly conserved domains as profile candidates.
        double  fractionMS,WeakEcutoff;
	double	fract_IDs;	// fraction of residue identities for weak hits.
	Int4	MinWeakAlnLen;
	BooLean	Use2LevelTree;

	// ************** VSI parameters: ****************
        Int4    NumRegions;
        Int4    *Start,*End,test;
        char    *colors;
};

#endif

