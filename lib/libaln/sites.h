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

/****************** sites.h - sites abstract data type.***************/
#if !defined(SITES)
#define SITES
#include "stdinc.h"
#include "olist.h"
#include "afnio.h"
#include "seqset.h"
#include "dheap.h"
#include "alphabet.h"
#include "residues.h"
#include "probability.h"
#include "random.h"
#include "gss_typ.h"
#include "set_typ.h"
/*************************** ADT sites *********************************
Defines regions of sequences to different types of elements.  

ntyp = 6; 	t = ( A  B  C  D  E  F)
len_elem = 	      5  5  3  3  5  6
maxinc = min(len_elem) - 1 = 2

type[1]   =    A  F  B  E  D
pos[1]    =    3 18 30 40 50 
 
type[2]   =    E  D  C  B  F  A  D  D
pos[2]    =    6 13 21 24 33 42 52 55

type[3]   =    C  E  A  F  D  B
pos[3]    =    4 10 18 23 32 38

nseq = 3	( 1   2   3)
len_seq =	 57  59  51
nsites[A] =       1   1   1
nsites[B] =       1   1   1
nsites[C] =       0   1   1
nsites[D] =       1   3   1
nsites[E] =       1   1   1
nsites[F] =       1   1   1

type[1]:     ..Aaaaa..........Ffffff......Bbbbb.....Eeeee.....Ddd.....*
type[2]:     .....Eeeee..Ddd.....CccBbbbb....Ffffff...Aaaaa.....DddDdd..*
type[3]:     ...Ccc...Eeeee...AaaaaFfffff...Ddd...Bbbbb.........*
 	     ....:....|....:....|....:....|....:....|....:....|....:....|
                 5   10   15   20   25    30  35   40   45   50   55   60

Note: '.' = VACANT; lowercase  = BLOCKED(t); '*' = ENDTYPE_SITE.

 ********************************************************************/

typedef struct {
	Int4	ntyp;		/* ntyp = number of types of elements. */
	Int4	*len_elem;	/* len_elem(t) = the length of element t */
	Int4	**nsites;	/* nsites(t,n) = # of type t sites in seq n */
	unsigned short	***site_pos;	/* site_pos[t][n][s]: site positions */
	Int4	maxinc;		/* maxinc = length of shortest element - 1 */
	Int4	*totsites;	/* totsites(t) = total # of type t sites */
	Int4	*tmp;		/* temp. buffer for site positions */
	ol_type *pos;		/* ordered list of site positions in seq n */
	char	**type;		/* type(n,i) = type of ith element in seq n */

	gss_typ	gss;		/* gapped sequence set */
	Int4	nseq;		/* number of sequences */
	Int4	*len_seq;	/* length of sequences */
} sites_type;
typedef sites_type *st_type;

typedef struct {
	gss_typ	gss;			/* gapped sequence set */
	Int4	ntyp;			/* number of types of elements. */
	Int4	*len_elem;		/* the length of element t */
	Int4	**nsites;		/* number of type t sites in seq n */
	unsigned short	***site_pos;	/* positions of sites */
} sites_info_type;
typedef sites_info_type *sti_typ;

/********************************* PRIVATE ********************************/
#define VACANT			0
#define BLOCKED(t)		-(t)
#define MAX_NO_TYPESITES	125
#define	ENDTYPESITE		-120
#define	MAX_LENG_SITES		120

void	print_sites(unsigned short *L);
Int4    bubble_sort_sites(unsigned short *L);
void    mk_sites(Int4 ntyps, Int4 *len_elem, st_type S);
Int4    *get_site_freq(st_type S, Int4 t, Int4 d, BooLean no_blocked);
void	sites_error(const char *s);

/********************************* PUBLIC ********************************/
/*-------------------------------- archive -----------------------------*/
sti_typ ArchiveSites(st_type S);
st_type ExtractSites(sti_typ X);
sti_typ CopyArchiveSites(sti_typ S);
void    NilArchiveSites(sti_typ X);
/*-------------------------------- Sites -------------------------------*/
st_type MkSites(Int4 ntyps, Int4 *len_elem, ss_type data, Int4 open, Int4 extend,
	double pernats, Int4 leftflank, Int4 rightflank); 
st_type MakeSites(Int4 ntyps, Int4 *len_elem, gss_typ &gss0);
st_type	CopySites(st_type S);
st_type CreateNewSites(Int4 ntyps, Int4 *len_elem, st_type S);
st_type StartSites(Int4 ntyps, Int4 *sitelen, gss_typ& gss);
void    InitSites(st_type S);
void	NilSites(st_type S);

Int4    GapBetweenSites(Int4 n, Int4 t, st_type S);
Int4    *GetSiteFreq(st_type S, Int4 t, Int4 d);
Int4    *GetSiteFreq(st_type S, BooLean *skip, Int4 t, Int4 d);
Int4    *GetSiteFreq(st_type S, set_typ Set, Int4 t, Int4 d);
Int4    *AlwaysGetSiteFreq(st_type S,Int4 t,BooLean right,Int4 *pos,char *blocked);
Int4    GetEdgeBlocksSite(st_type S,Int4 t,BooLean right,Int4 *pos,char *blocked);
Int4	PosSites(Int4 n, Int4 *pos, st_type S);
Int4	PosTSites(Int4 t, Int4 n, Int4 *pos, st_type S);

BooLean ColinearSites(st_type S);
BooLean OccupiedSite(register Int4 t, register Int4 n, register Int4 site, 
        register st_type S);

// BASIC OPERATIONS ON SITES:
Int4    MkRoomForSite(Int4 blk, Int4 sq, Int4 site, st_type S);
BooLean IsBlockedSite(Int4 blk, Int4 sq, Int4 site, st_type S);
void    AddSite(Int4 t, Int4 n, Int4 site, st_type S);
void	GrowSites(Int4 t, st_type S);
void    ShiftSites(st_type S, Int4 t, BooLean left);
void    ShiftSitesM(st_type S, Int4 t, Int4 d);
void	ShrinkSites(Int4 t, st_type S);
BooLean ShuffleCLSites(st_type S);
BooLean AddRandomCLSites(Int4 n, st_type S);
void    VacateSites(Int4 n, st_type S);
void	VacateSite(Int4 t, Int4 n, Int4 site, st_type S);
void	ReplaceSeqSites(Int4 n, gsq_typ *gsq, st_type S);
gsq_typ *SwapSeqSites(Int4 n, gsq_typ *gsq, st_type S);
st_type FuseElementsSites(Int4 x, Int4 maxlen, st_type S);
st_type FuseElementsSites2(Int4 x, Int4 maxlen, st_type S);
st_type SplitElementSites(Int4 x, Int4 minlen, st_type S);
st_type SplitElementSites(Int4 x, Int4 length1, Int4 minlen, st_type S);
st_type AddElementSites(Int4 lenx, Int4 x, Int4 min_free,
        Int4 min_squeeze, st_type S);
void	InsertGapSites(Int4 n, Int4 pos, unsigned short gap, st_type S);
Int4    *GetPosSites(Int4 sq,st_type S);

// OUTPUT FOR SITES:
void    PutPhylipSites(FILE *fptr,st_type S,BooLean **off);
void    PutScanSites(FILE *fptr, Int4 t, st_type S, BooLean *off);
void    PutScanSitesProb(FILE *fptr, Int4 t, st_type S, double **prob, BooLean *off,
	double cutoff);
void    PutSites(FILE *fptr,Int4 t,st_type S,double **site_prob, BooLean *off);
void	PutTypeSites(FILE *fptr, st_type S);

/********************************* MACROS ********************************/
#define AlphabetSites(S)	((S)->gss.Alphabet())
#define SitesGSS(S)		(&((S)->gss))
#define SitesSeqSet(S)		((S)->gss.FakeSqSet())
#define SitesTrueSqSet(S)	((S)->gss.TrueSqSet())
#define CountsSites(S)		((S)->gss.Counts())
#define NSeqsSites(S)		((S)->nseq)
#define SitePos(t,n,k,S)	(((S)->nsites[(t)][(n)]) < 2?((S)->site_pos[(t)][(n)][(k)]):\
				  bubble_sort_sites((S)->site_pos[(t)][(n)]),\
					(S)->site_pos[(t)][(n)][(k)])
#define EndSitePos(t,n,k,S)	(((S)->site_pos[(t)][(n)][(k)])+((S)->len_elem[(t)])-1)
#define TypeSite(n,s,S)		((S)->type[(n)][(s)])
#define StartSite(n,s,S)	((S)->type[(n)][(s)] > 0)
#define nSites(t,n,S)		((S)->nsites[(t)][(n)])
#define SiteLen(t,S)		((S)->len_elem[(t)])
#define SiteLengths(S)		((S)->len_elem)
#define ArchivedSiteLen(t,S)	((S)->len_elem[(t)])
#define SeqLenSites(n,S)	((S)->len_seq[(n)])
#define nTypeSites(S)		((S)->ntyp)
#define OpenPos(n,s,S)		(!(S)->type[(n)][(s)])
#define BlockedSite(t,S)      	-(t)
#define MaxTypeSites(S)      	MAX_NO_TYPESITES

#endif

