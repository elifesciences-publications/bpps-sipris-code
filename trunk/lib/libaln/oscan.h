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

/* oscan.h - codes and constants for ordered scan program. */
#if !defined (OSCAN)
#define OSCAN
#include "stdinc.h"
#include "afnio.h"
#include "scanheap.h"
#include "probability.h"
#include "histogram.h"
#include "sequence.h"
#include "smatrix.h"
#include "wmodel.h"
#include "mheap.h"
#include "sites.h"
#include "pseg.h"
#include "seqset.h"
#include "guide.h"
#include "residues.h"
// #include "spouge.h"
#include "sma.h"

typedef struct {
	a_type  A;
	BooLean	gaps;		/** look for gaps in segments **/
	FILE	*fp_gaps;	/** save sequences with gaps. **/
/****  pdb homologs ****/
	char	*pdbfa;		/** pdb fasta sequence file **/
	BooLean	weights;	/** weight input sequences **/
	BooLean	segmask;	/** check for compositional bias? **/
/****  SPOUGE GAPS ****/
	char	mode;		/** mode for spouge functions **/
	BooLean gapfunct;	/** use a gap function for search **/
	Int4	**scoreGap;	/** gap lengths **/
	Int4	maxLength;	/** maximum sequence length **/
	Int4	totGaps;	/** total number of gaps **/
/****  SPOUGE GAPS ****/
	char	*snfile;
	BooLean	shuffle;
	Int4	N;		/** number of motif models **/
	Int4	maxrpts;	/** maximum number of repeats **/
	Int4	permute;	/** circularly permute motifs **/
	double	total;
	double	*freq;
	double	maxEval,log10ecut;
	double	singleEval;	/** max single block Eval **/
	wm_type  *M;		/** weighted model **/
	wm_type  *M2;		/** alternative model (ignoring '^') **/
	char	method;		/** method used for model: g,c,r **/
	FILE	*dfp;		/** new database file pointer **/
	FILE	*recomb;	/** file pointer for recombinants **/
} oscan_type;
typedef oscan_type *osn_typ;

/********************************* private ********************************/
#define MAXSCN_BLOCK_SIZE  5000
#define MAX_BLOCK_LENGTH 200
#define MAX_NUM_MODELS	100

void	ReadOScan(osn_typ F, double pseudo,char *snfile, Int4 *counts);
BooLean ReadGapsOSCAN(FILE *file, osn_typ F);
BooLean ReadOScanCMSA(char *msafile, float minmap, double pseudo,Int4 *counts, 
        osn_typ F);
osn_typ make_oscan(double pseudo, char *snfile, Int4 *counts, a_type A, 
	Int4 maxrpts, char method, double ecut,double Ecut, char mode,
	Int4 maxLength, float minmap, BooLean weight);
/********************************* PUBLIC ********************************/
osn_typ MakeOScan(double pseudo, char *snfile, Int4 *counts, a_type A,
        Int4 maxrpts, char method, char mode, double ecut,double Ecut,
	Int4 maxLength, float minmap, BooLean weight);
osn_typ MkOScan(double pseudo, char *snfile, Int4 *counts, a_type A, 
	Int4 maxrpts, char method, double ecut,double Ecut);
snh_typ OScanScan(FILE *fptr,Int4 number, unsigned short *nsize, osn_typ F);
snh_typ OScanScan1(FILE *fptr,Int4 number, unsigned short *nsize, osn_typ F,
	BooLean combine);
void    NilOScan(osn_typ F);
FILE    *OpenDatabaseOScan(osn_typ F);
FILE    *OpenPermuteOScan(osn_typ F);
void	PermuteOScan(Int4 permute, osn_typ F);
void	GapsOScan(osn_typ F);
void    PutSMXOScan(FILE *fptr,osn_typ F);
void	PDBhomOScan(char *infile, osn_typ F);
void	SetEvalOScan(double ecut, double Ecut,osn_typ F);

/********************************* MACROS ********************************/
#define ShuffleSegsOScan(F)	((F)->shuffle = TRUE)
#define SetMethodOScan(x,F)	((F)->method = (char)(x))
#define NoMaskOScan(F)		((F)->segmask= FALSE)
#define NumModelsOScan(F)    	((F)->N)
#define SMatrixOScan(m,F)	(((m) <= (F)->N && (m) > 0)?\
					GetSMatrixWModel((F)->M[(m)]): NULL)

#endif

