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

/****************** sma.h - for storing msa's ***************/
#if !defined(_SMA)
#define _SMA
#include "stdinc.h"
#include "afnio.h"
#include "dheap.h"
#include "alphabet.h"
#include "residues.h"
#include "probability.h"
#include "dsets.h"
#include "dheap.h"
#include "random.h"
#include "histogram.h"
#include "seqset.h"
#include "sequence.h"
#include "block.h"
/*************************** ADT multaln *********************************/

typedef struct {
	Int4		*glength;   // gapped length of aligned block 
	char		**gnull;    // gapped null markers
	char		***gseq;    // gapped sequence

	double		*freq;
	Int4		nseqs;
	Int4		ntyps;
	Int4		totcol;
	BooLean		*ignore;
	Int4		*length,*cols;
	a_type		A;
	char		msaid[30],access[10],de[100];
	char		**null;
	float		*fieldmap,*single_map,*map_minus,lpr;
	float		**prob;
	char		**info,**id;
	Int4		**start,**end;
	unsigned char	***seq;  /** seq[n][t][s] **/
	BooLean		*use;   // use sequence n after purging?
} comultaln_type;
typedef comultaln_type *sma_typ;

/********************************* PRIVATE ********************************/
Int4    max_int_sma_array(Int4 *I, Int4 *L);
BooLean	*conserved_sma(Int4 *observed, double cutoff, sma_typ MA);
Int4    cluster_sma(Int4 cutoff, sma_typ MA);
BooLean *null_col_sma(Int4 t, double cbp_cut, sma_typ MA);
BooLean *purge_sma(Int4 cutoff, sma_typ MA);
BooLean *GoodSMA(double cutoff, Int4 *N, sma_typ MA);
sma_typ RmSMA(BooLean *keep, Int4 nseqs, sma_typ OMA);

/********************************* PUBLIC ********************************/
BooLean *ConservedSMA(Int4 t, Int4 s, sma_typ MA);
void    NetChargeSMA(FILE *fptr, sma_typ MA);
void	PutSelexSMA(FILE *fptr, e_type  *ListE, sma_typ MA);
void    PutClustalwSMA(FILE *fptr, Int4 blk, sma_typ MA);
void    RePutSeqsSMA(FILE *fptr, e_type *ListE, Int4 left, Int4 right, sma_typ MA);
BooLean	*FixSMA(double cutoff, sma_typ *NMA, sma_typ OMA);
sma_typ	ReadSMA(char *msafile);
sma_typ	ReadSMA(FILE *fp);
sma_typ	*MultiReadSMA(char *msafile, Int4 *n);
void	SMA2RTF(FILE *fptr, double cbp_cut, double infoLO, double infoHI, sma_typ MA);
void	SMA2RTF(FILE *fptr, double cbp_cut, double infoLO, double infoHI, 
	ss_type key_seq, sma_typ MA);
Int4    HspPurgeSMA(Int4 cutoff, sma_typ MA);
Int4    *GetGapsSMA(Int4 maxgap, float minprob, Int4 t1, Int4 t2, sma_typ MA);
BooLean	DiffSMA(FILE *fptr, sma_typ MA1, sma_typ MA2);
void    ReOrderSMA(FILE *fptr, char *file, sma_typ MA);
void	PutSMA(FILE *fptr, sma_typ MA);
void    PutTruncateBlkSMA(FILE *fptr, Int4 blk, Int4 remove, sma_typ MA);
BooLean	*PutPurgeSMA(FILE *fptr, Int4 cutoff, sma_typ MA);
void    PutGapsSMA(FILE *fp, sma_typ MA);
sma_typ TrimSMA(sma_typ MA);
Int4    CountsSMA(Int4 *counts, sma_typ MA);
Int4	ClusterSMA(sma_typ MA);
Int4	PurgeSMA(Int4 cutoff, sma_typ MA);
double  *WeightsSMA(sma_typ MA);
sma_typ	NilSMA(sma_typ MA);
void    PutDiffSMA(FILE *fptr, sma_typ MA1, sma_typ MA2);
Int4    GapLengthSMA(Int4 m, Int4 n, sma_typ MA);
Int4    *gnullInsrtLenSMA(Int4 blk, sma_typ MA);

float   **InfoSMA(Int4 purge_cutoff, sma_typ MA);
float   **ExcessInfoSMA(char *string, Int4 purge_cutoff, sma_typ MA);
/********************************* MACROS ********************************/
#define ignoreSMA(t,MA)		((MA)->ignore[(t)])
#define lengthSMA(t,MA)		((MA)->length[(t)])
#define glengthSMA(t,MA)	((MA)->glength[(t)])
#define fieldmapSMA(t,MA)	((MA)->fieldmap[(t)])
#define nullSMA(t,MA)		((MA)->null[(t)])
#define gnullSMA(t,MA)		((MA)->gnull[(t)])
#define seqSMA(t,n,MA)		((MA)->seq[(n)][(t)])
#define SeqSMA(n,MA)		((MA)->seq[(n)])
#define gseqSMA(t,n,MA)		((MA)->gseq[(n)][(t)])
#define seqidSMA(n,MA)		((MA)->id[(n)])
#define startSMA(t,n,MA)	((MA)->start[(n)][(t)])
#define endSMA(t,n,MA)		((MA)->end[(n)][(t)])
#define ntypSMA(MA)		((MA)->ntyps)
#define residueSMA(t,n,s,MA)	((MA)->seq[(n)][(t)][(s)])
#define AlphaSMA(MA)		((MA)->A)
#define nseqSMA(MA)		((MA)->nseqs)
#define freqSMA(r,MA)		((MA)->freq[(r)])
#define FreqSMA(MA)		((MA)->freq)
#define probSMA(s,t,MA)		((MA)->prob[(s)][(t)])
#define IsProbSMA(MA)		((MA)->prob != 0)
#define infoSMA(s,MA)		((MA)->info[(s)])
#define mapSMA(t,MA)		((MA)->map[(t)])
#define TotalColsSMA(MA)	((MA)->totcol)
#define DescriptionSMA(MA)	((MA)->de)

#endif

