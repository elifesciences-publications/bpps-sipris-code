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

/* scaninfo.h - scan information data type. */
#if !defined (SCANINFO)
#define SCANINFO
#include "stdinc.h"
#include "sequence.h"
#include "residues.h"

/************************** Scan Info Type ***************************

/*********************************************************************/

/*************************** Scan Info type **************************/
typedef struct scaninfo_link {
        e_type  		E;	/** database sequence **/
        float   		S;	/** sequence -log10(p-value) **/
        float   		*s;	/** segment log10(p-values) **/
        Int4			*site;	/** site of match **/
        unsigned char		rpts;	/** number of rpts **/
        unsigned char		nblk;	/** number of blocks **/
	Int4			id;	/** sequence id number **/
	char			method;	/** scanning method ('g','G',etc.) **/
	struct	scaninfo_link	*next;  /** link to next repeat **/
} scaninfo_type;

typedef scaninfo_type *sni_typ;

/*********************************************************************/

/******************************* private *******************************/
void    nil_scan_info(sni_typ I);

/******************************* PUBLIC *******************************/
sni_typ MakeScanInfo(e_type E, float score,Int4 nblk, Int4 rpts,
	float *blkscore, Int4 *site, Int4 seqid, char method);
sni_typ ReMakeScanInfo(sni_typ I, float score,Int4 nblk, Int4 rpts,
	float *blkscore, Int4 *site);
sni_typ MakeRptScanInfo(Int4 r, Int4 Nt_flank, Int4 Ct_flank, sni_typ I);
sni_typ AppendScanInfo(sni_typ I, sni_typ List);
sni_typ RmLastScanInfo(sni_typ List);
Int4    LengScanInfo(sni_typ List);
void    PutScanInfo(FILE *fptr,sni_typ I,Int4 *len,BooLean msaformat,a_type A);
void    NilScanInfo(sni_typ I);
e_type  NilScanInfoRtnE(sni_typ I);

/*********************************************************************/
#define MethodScanInfo(I)	((I)->method)
#define SeqIDScanInfo(I)	((I)->id)
#define RptsScanInfo(I)		((I)->rpts)
#define nBlksScanInfo(I)	((I)->nblk)
#define SeqScanInfo(I)		((I)->E)
#define ScoreScanInfo(I)	((I)->S)
#define SiteScanInfo(r,I)       ((I)->site[(r)])
#define scoreScanInfo(r,I)	((I)->s[(r)])

#endif

