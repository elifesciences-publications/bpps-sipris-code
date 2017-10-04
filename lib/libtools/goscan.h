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

/* goscan.h - codes and constants for ordered scan program. */
#if !defined (GOSCAN)
#define GOSCAN
#include "stdinc.h"
#include "afnio.h"
#include "scanheap.h"
#include "probability.h"
#include "histogram.h"
#include "sequence.h"
#include "smatrix.h"
#include "wmodel.h"
#include "pseg.h"
#include "guide.h"
#include "residues.h"
// #include "spouge.h"
#include "sma.h"
#include "cmsa.h"
#include "prtn_model.h"
#include "psm_typ.h"
#include "pah_typ.h"
#include "cls_typ.h"
#include "HMM_typ.h"

typedef struct {
	char	*snfile;
	Int4	num_models;	// number of protein models
	ptm_typ *pm;		// array of protein models
	ptm_typ	PM;		// first protein model == pm[1]**/
	BooLean	shuffle;
	double	maxEval,log10ecut;
	double	singleEval;	/** max single block Eval **/
	double	repeatEval;	/** max repeat Eval **/
	FILE	*dfp;		/** new database file pointer **/
	FILE	*recomb;	/** file pointer for recombinants **/
	psm_typ	*pssm;		// complex PSSM protein profile.
	cma_typ	pssm_cma;	// complex PSSM input alignment.
	HMM_typ	*hmm;		// hidden Markov model pointer.
} goscan_type;
typedef goscan_type *gsn_typ;

/********************************* private ********************************/
BooLean contain_segs_goscan(Int4 *s, Int4 *e, Int4 R, Int4 *gs,
        Int4 *ge, Int4 gR,double cutoff);

Int4    gap_seq_goscan(char *operation, Int4 start0, gsn_typ F, Int4 *left, Int4 *right);
Int4    gap_seq_goscan(char *operation,Int4 start0,gsn_typ F,Int4 *left,
        Int4 *right, char **SubOper);

/********************************* PUBLIC ********************************/
gsn_typ MakeGOScan(char *snfile, a_type A, Int4 maxrpts, char method, char mode, 
	double ecut,double Ecut, Int4 maxLength, BooLean weight,
	float minmap,double pseudo,double *freq);
gsn_typ MakeMGOScan(char *snfile, a_type A, Int4 maxrpts, char method, char mode, 
	double ecut,double Ecut, Int4 maxLength, BooLean weight,
	float minmap,double pseudo,double *freq);
snh_typ GOScanScan(FILE *fptr,Int4 number,UInt8 total,unsigned short *nsize,gsn_typ F);
snh_typ GOScanScan(FILE *fptr,Int4 number,UInt8 total,unsigned short *nsize,
		Int4 min_rpt, gsn_typ F);
snh_typ GOScanScan(FILE *fptr,Int4 number,UInt8 total,unsigned short *nsize,
		Int4 min_rpt, Int4 MaxBadBlks, gsn_typ F);
snh_typ *MGOScanScan(FILE *fptr,Int4 number,UInt8 total,unsigned short *nsize,gsn_typ F);
void	PutSmxGOScan(FILE *fptr,gsn_typ F);
void	GOScanSWScan(FILE *fptr,Int4 number,UInt8 total,unsigned short *nsize,
        gsn_typ F);
void    GOScanSWScan(FILE *fptr,Int4 number,UInt8 total,unsigned short *nsize,
        Int4 open, Int4 extend,gsn_typ F);
void    GOScanSWScan(FILE *fptr,Int4 number,UInt8 total,unsigned short *nsize,
        Int4 open, Int4 extend,gsn_typ F,char mode);
void    NilGOScan(gsn_typ F);
FILE    *OpenDatabaseGOScan(gsn_typ F);
void	SetEvalGOScan(double ecut, double Ecut,gsn_typ F);
void	NoMaskGOScan(gsn_typ F);
void    MaskNonGlobularGOScan(gsn_typ F);

e_type  **GapSeqGOScanScan(FILE *ofp, FILE *fptr,Int4 a, Int4 b, Int4 left, Int4 right,
	Int4 min_rpt, Int4 number,UInt8 total, unsigned short *nsize, gsn_typ F);
e_type  **GapSeqGOScanScan(FILE *fptr,Int4 a, Int4 b, Int4 left, Int4 right,
	Int4 min_rpt, Int4 number,UInt8 total, unsigned short *nsize, gsn_typ F);
Int4    GapSeqGOScanScan(FILE *fptr,Int4 a, Int4 b, Int4 gapo, Int4 gapx,
	Int4 left, Int4 right, Int4 min_rpt, Int4 number,UInt8 total, 
	unsigned short *nsize, char ***Operation, e_type **RtnE, e_type **FullE, 
	unsigned short **FullR, Int4 **Start, char Mode, gsn_typ F, 
	UInt4 minlen, UInt4 maxlen);
Int4    GapSeqGOScanScan(FILE *fptr,Int4 a, Int4 b, Int4 gapo, Int4 gapx,
	Int4 left, Int4 right, Int4 min_rpt, Int4 number,UInt8 total, 
	unsigned short *nsize, char ***Operation, e_type **RtnE, e_type **FullE, 
	unsigned short **FullR, Int4 **Start, char Mode, gsn_typ F);
cma_typ	GOScanToCMSA(FILE *fptr,Int4 a, Int4 b, Int4 gapo, Int4 gapx,
	Int4 left, Int4 right, Int4 min_rpt, Int4 number,UInt8 total, 
	unsigned short *nsize, char *name, gsn_typ F);
cma_typ GOScanToCMSA(Int4 argc,char *argv[], a_type A);
char    *DescriptionGOScan(Int4 i, gsn_typ F);

/********************************* MACROS ********************************/
#define ShuffleSegsGOScan(F)	((F)->shuffle = TRUE)
#define SetMethodGOScan(x,F)	((F)->method = (char)(x))
#define SetRptEvalGOScan(x,F)	((F)->repeatEval= (double)(x))
#define PrtnModelGOScan(F)	((F)->PM)
#define NumModelsGOScan(F)	((F)->num_models)

#endif

