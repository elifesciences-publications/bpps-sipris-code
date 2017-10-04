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

#if !defined(_GSEQ_)
#define _GSEQ_
#include "sequence.h"
#include "alphabet.h"
#include "idp_typ.h"

/*************************** Gapped Sequence Class *****************************
  realE = original (real 
  indels[0] == N-terminal region deleted from gsq.
  indels[fake_len] == C-terminal region deleted from gsq.
  f2r[0] == 0 (always).
  f2r[1] == real sequence position at start of fake seq.
  f2r[fake_len] == real sequence position at end of fake seq.
  f2r[fake_len] + indels[fake_len] == end of real sequence.
  inserts[0] + (f2r[fake_len] -f2r[1] +1) + inserts[fake_len] = real seq length.

  f2r[f] for '-' residue == position in real seq. just before '---'.
         for lowercase residue == position in real seq. just before inserted region
				  of length inserts[f].
Example:
 (implied model:     Mmmmmmmmmmmmiimmmmmmmmmmmmmmmd   dmmmmmmmmmmmmmd )
 Real:      plwlbvrpfIEVIGKENICGApgIVASNHRSHLDPPVL-dee-GGILKHMRAIPLR-rainsg*
 Fake:               IEVIGKENICGA  IVASNHRSHLDPPVLx   xGGILKHMRAIPLRx      
 inserts:  9         000000000002  0000000000000003   000000000000006
 del_bit:  0         000000000000  0000000000000001   100000000000001
 f2r[f]:   0        10...                                        
 f:        0         1...                                         

inserts:
unsigned short Ins_Mask = 0x7FFF; // 0111 1111 1111 1111
unsigned short Del_Mask = 0x8000; // 1000 0000 0000 0000
                                     |
                                     delete boolean bit

insertion:  (gseq_ins_mask & indels[i])= following insert length (range 0 to 32767). 
deletion:   if(gseq_del_mask & indels[i])...

 *******************************************************************************/
class gsq_typ {
public: 
		gsq_typ( ){ realE=0; init( ); }
		gsq_typ(e_type E) { this->init( ); this->initialize(E); }
		gsq_typ(const gsq_typ& gsq){ realE=0; init(); copy(gsq); }
		   // constructor for 'gsq_typ gsq2=gsq;' or 'gsq_typ gsq2(gsq);
		gsq_typ(Int4 lf, Int4 rf, char *op, Int4 trace_len, Int4 st, e_type E,Int4 *pos)
		  { init(); initialize(lf,rf,op, trace_len, st, E,pos); }
		gsq_typ(char *op, e_type E,Int4 *pos)
		  { init(); realE=0; initialize(op,E,pos); }
		~gsq_typ( ){ Free(); }

	gsq_typ& operator=(const gsq_typ&);	// assignment operator.
	//======================= gsq_init.cc =============================
	Int4    read(FILE *,Int4 ,Int4 *,Int4 *,a_type);
	void	initialize(const gsq_typ &gsq){ copy(gsq); }
	void	initialize(e_type);
	void	initialize(Int4,Int4,char *,Int4,Int4,e_type);
	void	initialize(Int4,Int4,char *,Int4,Int4,e_type,Int4*);
	e_type	initialize(Int4,Int4,char *,Int4,Int4); // reinitialize values
	e_type	initialize(Int4,Int4,char *,Int4,Int4,Int4 *); // reinitialize values

	void    initialize(Int4 nBlks,Int4 *len,Int4 *pos,e_type E);
	void    initialize(char *operation, e_type E, Int4 *pos);

	void	initializeX(Int4,Int4,char *,Int4,Int4,e_type,Int4*);
	e_type	initializeX(Int4,Int4,char *,Int4,Int4,Int4 *); // reinitialize values
private:
	Int4    read(FILE *fp,Int4 *pos,a_type A);
	void	init( );
public:
	//======================= gsq_init.cc =============================

	//======================= gsq_put.cc =============================
	void 	PutFake(FILE *,a_type);  // print domain in fasta format
	void	Put(FILE *, a_type);
	void	PutX(FILE *, a_type);
	void    Put(FILE *, Int4 , a_type );
	void	Put_cma_format(FILE *fp,Int4 J0,Int4 nblk,Int4 *start, Int4 *blk_len,a_type A);

	//======================= gsq_operate.cc =============================
	char    *MvColsOperation(char *MvOper, char *operation);
	gsq_typ *MoveColumns(char *MvOper, Int4 nBlks,Int4 *len, Int4 *sites);
	gsq_typ *IronOut(Int4 nblk,Int4 *start, Int4 *blk_len,Int4 &X);
	gsq_typ *IronOut(Int4 nblk,Int4 *start, Int4 *blk_len);
	gsq_typ *IronOut2(Int4 nblk,Int4 *start, Int4 *blk_len,Int4 &X);
	gsq_typ *ConvertColsToInsert(Int4 theBlk, Int4 start, Int4 end, Int4 nBlks, 
			Int4 *len, Int4 *sites);

	gsq_typ *InsertColumns(Int4 theBlk, Int4 start, Int4 Length, Int4 nBlks,
                        Int4 *len, Int4 *sites);
	char    *AddColumnToOperation(Int4 Blk, Int4 start, Int4 Length, char *operation);

	gsq_typ *FuseBlks(Int4 theBlk, Int4 nBlks, Int4 *len, Int4 *sites);
	gsq_typ	*ToOneBlk(Int4 nblk,Int4 *start, Int4 *blk_len,Int4 &pos);
#if 1	// DEBUG...
	gsq_typ	*ToOneBlk2(Int4 nblk,Int4 *start, Int4 *blk_len,Int4 &pos);
	void	ConvertToOneBlock2(char *Op);
#endif

	gsq_typ	*RmBlock(Int4 Blk,Int4 nBlks, Int4 *len, Int4 *sites);

	gsq_typ	*ExtendBlk(Int4 Blk,Int4 AddLeft,Int4 AddRight,Int4 nBlks, Int4 *len, Int4 *sites);
	gsq_typ	*TrimBlk(Int4 Blk,Int4 start,Int4 end,Int4 nBlks, Int4 *len, Int4 *sites);

	gsq_typ	*SplitBlock(Int4 Blk, Int4 *pos,Int4 nBlks, Int4 *Len, Int4 *sites);
	gsq_typ	*SplitSingleBlock(Int4 *pos,Int4 site, Int4 Len);

#if 1	// New operations for mapgaps analyses.
	gsq_typ *ExtendAlnToEnds(Int4 nblk,Int4 *start, Int4 *blk_len,Int4 &pos);
#endif
	gsq_typ *RmOverHangs(Int4 nblk,Int4 *sites, Int4 *blk_len,Int4 start, Int4 end);
	char    *RmFlankingInOperation(char *operation);

private:	// operation array manipulations
	char    *TrimBlkInOperation(Int4 Blk,Int4 start,Int4 end,char *operation);
	char    *ExtendBlkInOperation(Int4 Blk,Int4 AddLeft,Int4 AddRight,char *operation);
	Int4	SplitBlkInOperation(Int4 Blk, Int4 *pos,char *operation);
	void	SplitBlkInOperation(Int4 *pos,char *operation);
	char    *FuseBlksInOperation(Int4 Blk, char *operation);
	void	ConvertToOneBlock(char *Op);
	char    *RmBlkInOperation(Int4 Blk,char *operation);
	char    *AddInsertToOperation(Int4 Blk, Int4 start, Int4 end, char *operation);
	Int4    RmInsertByDeleteOperation(char *operation);

	char    *FromSitesOperation(Int4 nBlks,Int4 *len,Int4 *pos,e_type E);
	//======================= gsq_operate.cc =============================
public:
	char    *Operation(Int4 nblk,Int4 *start, Int4 *blk_len);
	char    *Operation2(Int4 nblk,Int4 *start, Int4 *blk_len);
	char    *Operation3(FILE *fp, a_type A);
	char    *Operation4(FILE *fp, a_type A);

	//======================= gsq_typ.cc =============================

	Int4    FakeToReal(Int4 i);
	Int4    RealToFake(Int4 i); // CAUTION: this in not efficient!
	Int4    Region(char *,Int4,Int4,a_type);
	BooLean Gapped( ) const { return ((num_insert+num_del) > 0); }
	BooLean Identical(const gsq_typ&);  // Are sequences identical?
	Int4	NumIns( ){ return num_insert; }
	Int4	NumDel( ){ return num_del; }
	Int4	NumOpen( ){ return num_open; }
	Int4	NumExtend( ){ return (num_insert + num_del - num_open); }
	e_type	FakeSeq( ){ return fakeE; }
	e_type	TrueSeq( ){ return realE; }
	Int4    OverHangN( );
	Int4    OverHangC( );
	Int4    IronOut(a_type);
	// void	AddRFF(rff_typ *rff){ assert(realE); AddRFF2Seq(realE,rff); }

// Insertion & deletion routines....
	Int4    CheckForInsertsAtEnds( );
	BooLean	IsDeleted(UInt4);
	Int4	Insertion(UInt4);
	BooLean NewAlign(Int4,Int4,char *,Int4,Int4);
	Int4	InDels(UInt4,UInt4,Int4 *,Int4 *);
	void	InsertGap(Int4,unsigned short, e_type);
	e_type	InsertGap(Int4,unsigned short);
	void    MkGapless(e_type E);
	Int4	GapPenalty(UInt4,idp_typ*);
	Int4	GapPenalty(UInt4,Int4, idp_typ*);
	void	FindIndels(Int4 nblks, Int4 *start, Int4 *len, UInt4 &Nmm, 
			UInt4 &Nmi, UInt4 &Nmd, UInt4 &Nii,
        		UInt4 &Nim, UInt4 &Ndd, UInt4 &Ndm,
			UInt4 &Ndi,UInt4 &Nid,UInt4 &Nsd,UInt4 &Nsm);
	void    AddTransitions(Int4 nblks, Int4 *start, Int4 *len,
        		UInt4 *Nmm, UInt4 *Nmi, UInt4 *Nmd, UInt4 *Nii,
        		UInt4 *Nim, UInt4 *Ndd, UInt4 *Ndm,
        		UInt4 *Ndi,UInt4 *Nid,UInt4 *Nsd, UInt4 *Nsm, char SqWt);
	char    *Operation(FILE *fp, a_type A);
	// Int4    InsGapCost(UInt4,unsigned short,Int4 &,Int4 &);
	// Int4    *Recombine(gsq_typ&,gsq_typ&,Int4 *,Int4 *,Int4 *);
private:
#if 1	// 32-bits...
	typedef unsigned int indel_typ;
	typedef unsigned int f2r_typ;
							
	static const indel_typ gseq_del_mask = 0x80000000; // 1000 0000 0000 0000 0000 0000 0000 0000
	static const indel_typ gseq_ins_mask = 0x7FFFFFFF; // 0111 1111 1111 1111 1111 1111 1111 1111
#else
	typedef unsigned short indel_typ;
	typedef unsigned short f2r_typ;
	static const indel_typ gseq_del_mask = 0x8000; // 1000 0000 0000 0000
	static const indel_typ gseq_ins_mask = 0x7FFF; // 0111 1111 1111 1111
#endif
	void    	copy(const gsq_typ& gsq);
	BooLean		is_deleted(Int4 i){ return ((indels[i] & gseq_del_mask)? TRUE: FALSE); }
	indel_typ	deletion(Int4 i){ return (indels[i] & gseq_del_mask); }
	indel_typ	insert(Int4 i){ return (indels[i] & gseq_ins_mask); }
	void		Free();		// free memory...
	e_type		realE;		// original ungapped full sequence
	e_type		fakeE;		// gapped subsequence
	indel_typ	*indels; 	// insert length (note: USHRT_MAX='-').
	f2r_typ		*f2r;		// position in realE of fakeE residue; 0 for '-'.
	UInt4	num_insert,num_del,num_open;
	// unsigned short	delo,delx,inso,insx; // Indel information.
};

#endif

