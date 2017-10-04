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

/************** fmodel.h - fragmented model abstract data type.***********/
#if !defined(FMODEL)
#define FMODEL
#include "stdinc.h"
#include "alphabet.h"
#include "random.h"
#include "probability.h"
#include "smatrix.h"
#include "blosum62.h"
#include "histogram.h"
#include "dheap.h"

/**************************** model ADT *****************************
        Product multinomial model for an ungapped element allowing 
	fragmentation.

r:  counts[r] freq[r]           
X( 0): 0      0.000             (X is a dummy residue)
A( 1): 806    0.064             freq[r] = counts[r]/tot_cnts
R( 2): 756    0.060             seqs:   H..GD..I.K
N( 3): 629    0.050                     H..AD..L.K
D( 4): 671    0.054                     H..RD..L.K
C( 5): 241    0.019                     H..RD..V.K
Q( 6): 577    0.046                     H..FD..I.T
E( 7): 701    0.056                     H..SD..I.S
G( 8): 691    0.055                     Y..RD..I.K
H( 9): 336    0.027                     C..RD..I.C
I(10): 680    0.054                     H..RD..L.A
L(11): 1376   0.110                     P..PE..L.K
K(12): 615    0.049			----+----+
M(13): 235    0.019			1   5   10
F(14): 471    0.038
P(15): 627    0.050                     totsites = 10
S(16): 1071   0.086                     length = 10
T(17): 704    0.056                     npseudo = 2
W(18): 158    0.013                     Ps[r] = npseudo * freq[r]
Y(19): 434    0.035
V(20): 730    0.058
tot_cnts: 12509

observed[pos][r]:

        pos     A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V 
     |   0    (nonsite observed frequencies)
     |   1    (null)
     :
     |  x-1   (null)
start+-> x      0  0  0  0  1  0  0  0  7  0  0  0  0  0  1  0  0  0  1  0
     |  x+1   (null)
     |	x+2   (null)
     |  x+3     1  5  0  0  0  0  0  1  0  0  0  0  0  1  1  1  0  0  0  0
     |  x+4     0  0  0  9  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0
     |  x+5   (null)
     |  x+6   (null)
     |  x+7     0  0  0  0  0  0  0  0  0  5  4  0  0  0  0  0  0  0  0  1
     |  x+8   (null)
 end +->x+9     1  0  0  0  1  0  0  0  0  0  0  6  0  0  0  1  1  0  0  0
     |  x+10  (null)
     :
     | maxlen (null)

********************************************************************/

typedef struct {
	a_type	A;		/* alphabet */
	Int4	length;		/* length of motif model */
	Int4	ncols;		/* number columns in motif model */
	Int4	maxlen;		/* maximum possible length of model */
	Int4	start;		/* pointer to first column */
	Int4	end;		/* pointer to last column */
	Int4	**observed;	/* observed(j,b)=# b's at pos j of site*/
	UInt4	*observedNS;	// observed counts for nonsites.
	UInt4	*counts;	// observed counts for nonsites.
	double	*targetNS;	/* nonsite target frequencies */
	double	**likelihood;	/* likelihood ratio for being in model */
	double	npseudo;	/* npseudo = # pseudo priors */
	double	*Ps;		/* Ps(nres) = # pseudo residues by type */
        double	*freq;          /* freq(b)= total frequency of b's */ 
        double	*sumlngamma;    /* sumlngamma(i)= column i sum of lngamma b's */ 
	double	*tmp_val;	/* temp value(maxlen) = temp double array */
        double	**blsm62p;      /* probabilities from blosum62 matrix */ 
	Int4	totsites;	/* number of sites in model */
	BooLean	update,recalc;  /* update the normalized freq */
} fmodel_type;
typedef fmodel_type *fm_type;

#define FMODEL_UNDEF	-9999
/********************************* private ********************************/
Int4	add_column_fmodel(Int4 *observed, Int4 pos, fm_type M);
void	center_model(Int4 pos, fm_type M);
void	fmodel_error(const char *s);
void    update_fmodel(fm_type M);

double  score_fmodel(register unsigned char r, register Int4 s,
        register double **likelihood);
double  likelihood_fmodel(register unsigned char *seq, register Int4 s, 
        register double **likelihood);
/********************************* PUBLIC ********************************/
#define	LikelihoodFModel(s,p,M)	( ((M)->update) ? update_fmodel(M), \
	likelihood_fmodel((s+p),((M)->end-(M)->start),((M)->likelihood+(M)->start)):\
	likelihood_fmodel((s+p),((M)->end-(M)->start),((M)->likelihood+(M)->start)))
#define	ScoreFModel(r,s,M)	( ((M)->update) ? update_fmodel(M), \
	score_fmodel((r),(s)-1,((M)->likelihood+(M)->start)):\
	score_fmodel((r),(s)-1,((M)->likelihood+(M)->start)))

Int4	AddColumnFModel(Int4 *observed, Int4 pos, fm_type M);
void    Add2FModel(unsigned char *seq, Int4 site, fm_type M);
Int4    WorstOnFModel(fm_type M);
Int4	ChoiceLemonFModel(fm_type M);
Int4	ChoiceOrangeFModel(fm_type M);
fm_type CopyFModel(fm_type M);
Int4	LemonFModel(fm_type M);
fm_type MkFModel(BooLean *null, Int4 length, Int4 maxlen, double npseudo,
        UInt4 *observedNS,UInt4 *counts,double *freq,a_type A);
Int4	MvColumnFModel(Int4 *observed, Int4 lemon, Int4 pos, fm_type M);
fm_type	NilFModel(fm_type M);
double  NormLikelihoodFModel(fm_type M);
BooLean NullSiteFModel(Int4 s,fm_type M);
Int4	NullSitesFModel(BooLean *null, fm_type M);
Int4    OrangeFModel(fm_type M);
double  ProbFModel(register unsigned char *seq, register Int4 pos, register double p,
        register fm_type M);
void    PutFModel(FILE *fptr, fm_type M);
void    PutFModelShort(FILE *fptr, fm_type M);
double  LnRatioFModel(Int4 *observed, Int4 d, fm_type M);
double  RatioFModel(Int4 *observed, Int4 d, fm_type M);
BooLean EnlargeFModel(Int4 maxlen, fm_type M);
BooLean	RmColumnFModel(Int4 pos, fm_type M);
Int4	RmColumnFModel2(Int4 pos, fm_type M);
void	RmFModel(unsigned char *seq, Int4 site, fm_type M);
void    SetPseudoFModel(double npseudo, fm_type M);
void    ShiftFModel(Int4 *observed, BooLean left, fm_type M);
Int4    ObservedFModel(Int4 **array, fm_type M);
Int4    *SeeColumnFModel(Int4 c, fm_type M);
float   *InfoFModel(fm_type M);

smx_typ GetSmatrixFModel(double pernats, fm_type M);
smx_typ SampleSmatrixFModel(double pernats, Int4 wt, fm_type M);
smx_typ SampleDirichletSmatrixFModel(double pernats, Int4 wt, fm_type M);
/********************************* MACROS ********************************/
#define IsColumnFModel(j,M)	((M)->observed[(j+(M)->start - 1)]!=NULL)
#define LenFModel(M)		((M)->length)
#define ContigFModel(M)		((M)->ncols==(M)->length)
#define PseudoFModel(M)		((M)->Ps)
#define nPseudoFModel(M)	((M)->npseudo)
#define nColsFModel(M)		((M)->ncols)
#define TotSitesFModel(M)	((M)->totsites)
#define MaxLenFModel(M)		((M)->maxlen)

#endif

