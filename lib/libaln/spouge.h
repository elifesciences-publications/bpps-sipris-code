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

/* spouge.h - code to call spouge procedures. */
#if !defined (SPOUGE)
#define SPOUGE
#include "stdinc.h"
#include "afnio.h"
#include "probability.h"
#include "smatrix.h"
#include "histogram.h"
#include "residues.h"
#include "sequence.h"
#include "smooth.h"

#if 0
#include "sls_pssm.h"

using namespace Sls
#endif

typedef struct {
        size_t  num_aa;			/* #(amino acids) [0...dimAminos-1] */
        double  *freq;			/* corresponding array of amino freq's */
        size_t  num_mtfs;		/* #(motifs) */
        size_t  *lengthMotif;		/* lengthMotif [0...dimMotif-1] */
        Int4    ***scoreAmino;		/* position-wise amino acid scores */
        size_t  maxLength;		/* max (protein length) */
        size_t	*dimscoreGap;		/* gap score dimensions: NULL no gaps */
        Int4    **scoreGap;		/* gap scores : NULL for no gapping */
        Int4    *threshold;		/* threshold [1...maxLength] */
	char	mode;			/* method of gapping */
	void	*neuwald;		/* New spouge routine */
	size_t	**sentinel;		/* sentinel [0...dimMotif-1][2] */
	Int4	sum_extends;		/* sum extensions **/
	BooLean	initialize;		// have values been initialized?
	Int4	*maxgap;		// maximum gap size observed for input
        Int4	*lengthMotif0;		/* lengthMotif [0...dimMotif-1] */
} spouge_type;
typedef spouge_type *js_type;

/********************************* private ********************************/
#ifdef __cplusplus
extern "C" {
#endif

#if 1
extern  void	*Vec_NewLocalDisjointNeuwald(size_t, size_t, double *, size_t,
			size_t *,size_t **, Int4 ***, size_t *, Int4 **);
extern  void	*Vec_NewGlobalNeuwald(size_t, size_t, double *, size_t,
			size_t *,size_t **, Int4 ***, size_t *, Int4 **);
extern  void	*Vec_NewLocalOverlapNeuwald(size_t, size_t, double *, size_t,
			size_t *,size_t **, Int4 ***, size_t *, Int4 **);
extern  void	*Vec_NewLocalCoreDisjointNeuwald(size_t , size_t ,
			double *, size_t , size_t *, size_t **,
			Int4 ***, size_t *, Int4 **);
extern  char    Vec_TestSentinel(    /* tests for valid sentinels */
			size_t , size_t *, size_t **);
extern  Int4	Vec_NeuwaldThreshold(void *,size_t,double);
extern  double	Vec_NeuwaldTail(void *,size_t,Int4);
extern  void	*Vec_FreeNeuwald(void *);
#endif

#ifdef __cplusplus
}
#endif

/********************************* private ********************************/
void    SmoothGapScoresSpouge(Int4 nl, Int4 nr, Int4 m, js_type S);
void	InitPsiGap (js_type S);
void    MedianSpouge(double *median, js_type S);
/********************************* PUBLIC ********************************/
js_type	MkSpouge(Int4 num_aa, double *freq, Int4 num_mtfs, smx_typ *sM,
	Int4 maxLength, Int4 **scoreGap, char mode);
void	PutHistSpouge(FILE *fp, Int4 t, double inc, js_type S);
void    NilSpouge(js_type S);
Int4    MaxGapScoreSpouge(js_type S);

Int4	SpougeThreshold(js_type S,double target,Int4 length);
double  SpougePvalue(js_type S,UInt8 dbs_leng, Int4 length, Int4 score);
/********************************* MACROS ********************************/
#define EndCoreSpouge(m,S)		((S)->sentinel[m-1][1]+1)
#define StartCoreSpouge(m,S)		((S)->sentinel[m-1][0]+1)
#define GapScoreSpouge(S)		((S)->scoreGap)
#define GapScoreSpouge1(r,S)		((S)->scoreGap[(r)])

#endif

