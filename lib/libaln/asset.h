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

#if !defined (ASSET)
#define	ASSET
#include "stdinc.h"
#include "afnio.h"
#include "histogram.h"
#include "pheap.h"
#include "pattern.h"
#include "block.h"
#include "eblocks.h"
#include "align.h"
#include "seqset.h"
#include "segment.h"
#include "probability.h"
#include "residues.h"
#include "evalue.h"
#include "random.h"
/*************************** ADT ASSET ***************************
	
**********************************************************************/

/*************************** ASSET type **************************/
typedef struct {
	/***** POPULATION *****/
	char		*infile;	/* input file name */
	ss_type		P;		/* aligned segment population */
	double		*tfreq;		/* total residue frequencies */
	BooLean		xnu;		/* remove low complexity regions */
	BooLean		shuffle;	/* shuffle sequence set */
	BooLean		query;		/* query sequence mode - see usage */
	Int4		seed;		/* random seed */
	a_type		A;		/* alphabet */
        /**** BLOCKS  ****/
	ebs_typ		eblocks;	/* elementary blocks */
	b_type		*EOB;		/* exclusive or blocks with lists */
	b_type		IB;		/* intersection block with list */
	b_type		UB;		/* union block */
	b_type		*B;		/* working blocks B[d] */
	/***** PATTERNS *****/
	ph_type		pheap;		/* min-max heap for patterns */
	ptn_typ		Q;		/* primary pattern */
	Int4		hpsz;		/* motif heap size */
	char		**motif;	/* buffer for printing motif */
	/***** search parameters *****/
	Int4		k_max;		/* segment size */
	Int4		c_min;		/* minimum block size */
	Int4		d_max;		/* maximum search depth */
	Int4		d_min;		/* minimum search depth */
	Int4		s_min;		/* minimum # 1-residue positions */
	Int4		n_min;		/* minimum # matching sequences */
	Int4		mode;		/* search mode (0-30?) */
	Int4		end;
	Int4		maxpat;		/* maximum pattern printed */
	Int4		minscore;	/* minimum score for segments */
	BooLean		verbose;	/* TRUE == print out lots of info. */
	Int4		percent;	/* sequences w/ motif for scan file */
	/***** STORAGE *****/
	Int4		*C;		/* intersect cardinality C[d][c] */
	Int4		**Card;		/* intersect cardinality C[d][c] */
	Int4		*cnts;		/* # of intersections at depth d*/
	/***** E-values *****/
	evl_typ		E;
	double		min_prob;	/* log10 cutoff for calc. prob. */
	double		pout;		/* output significance level */
	/***** STATS *****/
	Int4		ncalc;		/* number of probabilities calculated */
	double		sd;		/* standard deviations */
	h_type		H;		/* histogram */
} asset_type;
typedef asset_type *ast_typ;

/******************************* PRIVATE *****************************/
ast_typ InitAsset(ast_typ D);
BooLean PutAsset(FILE *fptr,ast_typ D, Int4 offset, Int4 length);
Int4	EvaluateAsset(double *p, Int4 *idex, ast_typ D);
void	asset_dfps(ast_typ D,Int4 n,Int4 idex, Int4 CB, b_type *B,double p1);
void    asset_dfps12B(ast_typ D,Int4 n,Int4 i,Int4 CB, b_type *B,double p1);
void    asset_dfps12A(ast_typ D,Int4 n, Int4 idex, Int4 CB, b_type *B,double p1);
void    asset_dfps12T(ast_typ D,Int4 n,Int4 i,Int4 CB, b_type *B,double p1);
void	asset_error(const char *s);

/******************************* PUBLIC ******************************/
ast_typ	AssetSearch(ast_typ D);
ast_typ CreateAsset(char *DBS_NAME, a_type A);
void	NilAsset(ast_typ D);
ph_type NilAssetRtnPHeap(ss_type *data, ast_typ D);
void	PutAssetHeap(FILE *fptr,ast_typ D);
void	PutAssetIDs(FILE *fptr,ast_typ D);
void    PutInfoAsset(FILE *fptr, ast_typ D);
BooLean	SetAsset(char *command, ast_typ D);

/**************************** macro operations **********************/
#define	AssetPHeap(D)	((D)->pheap)
#define	AssetSeqSet(D)	((D)->P)
#define	AssetAlpha(D)	((D)->A)
#define MINDEPTH        2       /* minimum depth to store pattern */
#define MINSEGMENT      3

#define	MAXSEGMENT	100	/* this must be set <= 127 */
#define	USAGE_ASSET	"\nUSAGE: asset file [options]\n\
   options:\n\
     [-c<int>]    - set c_min equal to <int> (default = 4)\n\
                    c_min = minimum required # of matching segments\n\
     [-D<int>]	  - set d_max equal to <int> (default = 6)\n\
                    d_max = maximum # of pattern positions\n\
     [-d<int>]	  - set d_min equal to <int> (default = 4)\n\
                    d_min = minimum # of pattern positions\n\
     [-e<int>]    - maximum E-value needed to report pattern\n\
     [-f<int>]	  - create a scan file for alignments where the motif\n\
                    is present in at least <int>% of the sequences\n\
     [-h<int>]    - heap size (maximum number of patterns reported)\n\
     [-k<int>]	  - set k_max equal to <int> (default = 15)\n\
                    k_max = maximum pattern length\n\
     [-l]	  - DON'T eliminate low complexity sequences\n\
     [-n<int>]    - minimum required # of matching sequences for patterns\n\
     [-O<int>]	  - minimum segment block score (default = 100)\n\
     [-o<int>]	  - number of patterns shown (default = 20)\n\
     [-q]         - only output motifs to the scan file if they occur\n\
                    in the first sequence\n\
     [-s<int>]	  - s_min = minimum # of 1-residue positions in patterns\n\
                               suggested settings (default = 4):\n\
                       sequence relationships:     s_min:     speed:\n\
                          relatively close            5      fastest\n\
                         moderate to distant          4        :\n\
                            very distant              3      slowest\n\
     [-x<real>]   - sample E-value if log10(multiplication E-value) > <real>\n\n\
\n   It is best to eliminate closely related sequences prior to\n\
   analysis using the \"purge\" program with a cutoff score of\n\
   300 or less.\n"

#endif

