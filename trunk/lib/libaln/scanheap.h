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

/* scanheap.h - scan hit min/max heap. */
#if !defined (SCANHEAP)
#define SCANHEAP
#include "afnio.h"
#include "alphabet.h"
#include "scaninfo.h"
#include "histogram.h"
#include "mheap.h"
#include "set_typ.h"

/************************** Scan Heap Type ***************************

/*********************************************************************/
typedef struct {
	sni_typ	*info;		/** information on hits **/
	mh_type	mH;		/** heap for best hits **/
	Int4    hpsz;		/** heap size **/
	Int4	nblk;		/** number of blocks **/
	h_type  HG;		/** histogram for aggregate score **/
	h_type  *hg;		/** histograms for blocks **/
	Int4	*len;		/** lengths of blocks **/
	/********** processed info **********************/
	double	Ecut;		/** segment Evalue cutoff **/
	double	neglog10Ecut;	/** -log10 Evalue cutoff **/
	sni_typ *sort;		/** sorted hit information **/
	BooLean	*okay;		/** is sequence match okay? **/
	Int4	*not_okay;	// number of blocks in sequence i not okay.
	Int4	max_not_ok;	// maximum number of bad blocks allowed.
	Int4	number;		/** length of array **/
	BooLean	*use;		/** which blocks are to be used **/
	double	use_cutoff;	/** mean cutoff to use block **/
} scanheap_type;

typedef scanheap_type *snh_typ;
/*********************************************************************/

/******************************* private *******************************/
void    process_scanheap(snh_typ H);
void    put_info_scan_heap(FILE *fptr, snh_typ H, BooLean msaformat, a_type A);

/******************************* PUBLIC *******************************/
snh_typ	MakeScanHeap(Int4 hpsz,Int4 nblk,Int4 *len, double Ecut);
snh_typ	MakeScanHeap(Int4 hpsz,Int4 nblk,Int4 *len, double Ecut,Int4 MaxNotOK);
void	NilScanHeap(snh_typ H);
sni_typ *RtnNilScanHeap(snh_typ H);
Int4    InsertScanHeap(sni_typ I, snh_typ H);
sni_typ	DelMinScanHeap(snh_typ H);
sni_typ	DelMaxScanHeap(snh_typ H);
Int4    PutSeqScanHeap(FILE *fptr, snh_typ H, a_type A);
Int4    PutTopSeqScanHeap(FILE *fptr, Int4 number, snh_typ H, a_type A);
Int4    PutFailSeqScanHeap(FILE *fptr, snh_typ H, a_type A);
void	PutInfoScanHeap(FILE *fptr, snh_typ H, a_type A);
void    PutFormatMSAScanHeap(FILE *fptr, snh_typ H, a_type A);
void    PutMSAScanHeap(FILE *fptr, snh_typ H, a_type A);
void	PutGapsScanHeap(FILE *fptr, snh_typ H);
Int4    AddSetScanHeap(snh_typ H, set_typ S, BooLean AddAll);
void	PutFAAlnScanHeap(FILE *fptr, snh_typ H, a_type A);
void	PutRptsScanHeap(FILE *fptr, Int4 left, Int4 right, snh_typ H, a_type A);
void    PutTopRptsScanHeap(FILE *fptr, Int4 left, Int4 right, Int4 number,
	snh_typ H, a_type A);
Int4    PutMinRptsScanHeap(FILE *fptr,Int4 left, Int4 right, Int4 min_rpt,
	snh_typ H,a_type A);
void    PutRptsScanHeapFull(FILE *fptr, Int4 left, Int4 right, Int4 number,
	UInt4 min_rpt, snh_typ H, a_type A);
void	PutSelexAlnScanHeap(FILE *fptr, snh_typ H, a_type A);
void	PutRasmolScanHeap(FILE *fptr, snh_typ H);
void	PutMaskSeqScanHeap(FILE *fptr, snh_typ H, BooLean neg, a_type A);
void    PutSeqBlkScanHeap(char *name, snh_typ H, Int4 Nflank, Int4 Cflank, 
	a_type A);
Int4    NumFailSeqScanHeap(snh_typ H);

Int4    AdjustBlksScanHeap(snh_typ H);
sni_typ *InfoScanHeap(Int4 *n, snh_typ H);
sni_typ	*DomainsScanHeap(Int4 *n, Int4 Nt_flank, Int4 Ct_flank, snh_typ H);
double  *WeightsScanHeap(snh_typ H, a_type A);
/*********************************************************************/
#define nScanHeap(H)		(((H)->number < 0)? ItemsInMheap((H)->mH):\
					(H)->number)
#define SizeScanHeap(H)		((H)->hpsz)
#define nblksScanHeap(H)	((H)->nblk)
#define BlkLenScanHeap(r,H)	((H)->len[(r)])
#define LengthsScanHeap(H)	((H)->len)

#define HistScanHeap(H)		((H)->HG)
#define histScanHeap(r,H)	((H)->hg[(r)])

#endif

