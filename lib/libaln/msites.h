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
#if !defined(_MSITES_)
#define _MSITES_
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
} msites_type;
typedef msites_type *mst_type;

/********************************* PRIVATE ********************************/
#define VACANT			0
#define BLOCKED(t)		-(t)
#define MAX_NO_TYPESITES	125
#define	ENDTYPESITE		-120

void	print_msites(unsigned short *L);
Int4    bubble_sort_msites(unsigned short *L);
void    mk_msites(Int4 ntyps, Int4 *len_elem, mst_type S);
void	msites_error(char *s);

/********************************* PUBLIC ********************************/
/*-------------------------------- MSites -------------------------------*/
mst_type MkMSites(Int4 ntyps, Int4 *len_elem, ss_type data);
mst_type MakeMSites(Int4 ntyps, Int4 *len_elem, gss_typ &gss0);
void	NilMSites(mst_type S);

Int4    GapBetweenMSites(Int4 n, Int4 t, mst_type S);
Int4	PosMSites(Int4 n, Int4 *pos, mst_type S);
Int4	PosTMSites(Int4 t, Int4 n, Int4 *pos, mst_type S);

BooLean OccupiedMSite(register Int4 t, register Int4 n, register Int4 site, 
        register mst_type S);

// BASIC OPERATIONS ON SITES:
void    AddMSite(Int4 t, Int4 n, Int4 site, mst_type S);
void    VacateMSites(Int4 n, mst_type S);
void	VacateMSite(Int4 t, Int4 n, Int4 site, mst_type S);

// OUTPUT FOR SITES:
void    PutScanMSites(FILE *fptr, Int4 t, mst_type S, BooLean *off);
void    PutScanMSitesProb(FILE *fptr, Int4 t, mst_type S, double **prob, BooLean *off,
	double cutoff);
void    PutMSites(FILE *fptr,Int4 t,mst_type S,double **site_prob, BooLean *off);
void	PutTypeMSites(FILE *fptr, mst_type S);

/********************************* MACROS ********************************/
#define MSitesGSS(S)		(&((S)->gss))
#define MSitesSeqSet(S)		((S)->gss.FakeSqSet())
#define MSitesTrueSqSet(S)	((S)->gss.TrueSqSet())
#define CountsMSites(S)		((S)->gss.Counts())
#define NSeqsMSites(S)		((S)->nseq)
#define MSitePos(t,n,k,S)	(((S)->nsites[(t)][(n)]) < 2?((S)->site_pos[(t)][(n)][(k)]):\
				  bubble_sort_msites((S)->site_pos[(t)][(n)]),\
					(S)->site_pos[(t)][(n)][(k)])
#define EndMSitePos(t,n,k,S)	(((S)->site_pos[(t)][(n)][(k)])+((S)->len_elem[(t)])-1)
#define TypeMSite(n,s,S)		((S)->type[(n)][(s)])
#define StartMSite(n,s,S)	((S)->type[(n)][(s)] > 0)
#define nMSites(t,n,S)		((S)->nsites[(t)][(n)])
#define MSiteLen(t,S)		((S)->len_elem[(t)])
#define SeqLenMSites(n,S)	((S)->len_seq[(n)])
#define nTypeMSites(S)		((S)->ntyp)
#define OpenPos(n,s,S)		(!(S)->type[(n)][(s)])
#define BlockedMSite(t,S)      	-(t)
#define MaxTypeMSites(S)      	MAX_NO_TYPESITES

#endif

