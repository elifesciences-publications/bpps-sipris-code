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

/************** wmodel.h - weighted model abstract data type.***********/
#if !defined(WMODEL)
#define WMODEL
#include "stdinc.h"
#include "alphabet.h"
#include "haussler.h"
#include "smatrix.h"
#include "blosum62.h"
#include "blosum45.h"
#include "random.h"

/**************************** model ADT *****************************
        Product multinomial model for an ungapped element allowing 
	fragmentation.

r:  counts[r] freq[r]           
X( 0): 0      0.000             (X is a dummy residue)
A( 1): 806    0.064             
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
S(16): 1071   0.086                     length = 5
T(17): 704    0.056                     npseudo = 2
W(18): 158    0.013                     N0[r] = npseudo * freq[r]
Y(19): 434    0.035
V(20): 730    0.058
tot_cnts: 12509

        site_freq[pos][r] = N0[r] + counts[r]
        pos   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V 
     |  1    0  0  0  0  1  0  0  0  7  0  0  0  0  0  1  0  0  0  1  0
     |	2 (null)
     |	3 (null)
     |  4   1  5  0  0  0  0  0  1  0  0  0  0  0  1  1  1  0  0  0  0
     |  5   0  0  0  9  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0
     |  6 (null)
     |  7 (null)
     |  8   0  0  0  0  0  0  0  0  0  5  4  0  0  0  0  0  0  0  0  1
     |  9 (null)
     |  10   1  0  0  0  1  0  0  0  0  0  0  6  0  0  0  1  1  0  0  0

        site_freq0[r] (nonsite model) == # r residues not in model + N0[r].

        site_freqN[0][r] = site_freq[0][r]/sum(site_freq[0][r=A..V]-N0[r=A..V])

 	factor = totsites 

        site_freqN[pos][r] = site_freq[pos][r]/(factor*site_freqN[0][r])

********************************************************************/

typedef struct {
	a_type	A;		/* alphabet */
	dmp_typ	D;		/* for 'd' method */
	smx_typ	smx;		/* smatrix for p-values */
	double	maxscore;	/* maximum possible score */
	Int4	length;		/* length of motif model */
	double	**site_freq;	/* site_freq(j,b)=# b's at pos j of site*/
	double	totsites;	/* number of sites in model */
	double	**likelihood;	/* likelihood ratio for being in model */
	double	pseudo;		/* pseudo = base pseudo count priors */
	double	npseudo;	/* npseudo = # pseudo priors */
	double	*N0;		/* N0(nres) = # pseudo residues by type */
        double	*freq;          /* freq(b)= total frequency of b's */ 
	double	*temp;		/* temp(nres) = temp double array */
	double	*tmp_val;	/* temp value(maxlen) = temp double array */
	char	method;		/* method to calc Staden pval (default: g)*/
	char	*null;		/* null site or other...?? */
	BooLean	update;		/* update the normalized freq */
} wmodel_type;
typedef wmodel_type *wm_type;

/********************************* PRIVATE ********************************/
#define WMODEL_MAX_SCORE 200.0
void    update_wmodel_freqN(wm_type M);
void    get_smx_wmodel(register wm_type M);
void	wmodel_error(char *s);
double  InfoColWModel(Int4 col, wm_type M);
/********************************* PUBLIC ********************************/
wm_type WModel(Int4 length,double npseudo,double *freq,a_type A);
wm_type MkWModel(char *null, Int4 length,double npseudo,double *freq,a_type A);
wm_type MergeWModels(Int4 N, Int4 *p, Int4 leng, double *freq, wm_type *M);
double  **RealScoresWModel(wm_type M);
double  ExpectedScoreWModel(BooLean *use, wm_type M);
void    InitWModel(wm_type M);
void    Add2WModel(unsigned char *seq, Int4 site, double w, wm_type M);
wm_type NilWModel(wm_type M);
Int4	MaxSegWModel(unsigned char *maxseq, wm_type M);
Int4	ScoreWModel(register unsigned char *seq, register Int4 pos, register wm_type M);
Int4	SubScoreWModel(unsigned char *seq, Int4 pos, Int4 start, Int4 end, wm_type M);
double  PvalWModel(register unsigned char *seq, register Int4 pos, register wm_type M);
smx_typ GetSMatrixWModel(wm_type M);
double	PutWModel(FILE *fptr, wm_type M);
BooLean NullSiteWModel(Int4 s,wm_type M);
Int4	CellScoreWModel(Int4 r, Int4 pos, wm_type M);
double  *ObservedWModel(Int4 j, wm_type M);
Int4    GetSegWModel(unsigned char *seg, wm_type M);
unsigned char   GetBackGroundWModel(wm_type M);
/********************************* MACROS ********************************/
#define LenWModel(M)		((M)->length)
#define SetMethodWModel(x,M)	((M)->method = (char) (x))
#define NumResWModel(j,b,M)	((M)->site_freq[(j)][(b)])
#define IsColumnWModel(j,M)	((M)->null[(j)] != '.')
#define TotSitesWModel(M)	((M)->totsites)

#endif
