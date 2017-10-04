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

#if !defined(_GDATA_)
#define _GDATA_
#include "sequence.h"
#include "seqset.h"
#include "gsq_typ.h"
#include "idp_typ.h"

/*********************** Gapped Sequence Set Class ***************************
   truedata = ungapped original sequence set.
   fakedata = gapped sequence set derived from truedata.
   if fakedata == NULL then uses truedata.
 *****************************************************************************/

class gss_typ {
public: 
		gss_typ( ){ init( ); truedata=NULL; }
		gss_typ(ss_type seqset,Int4,Int4,double,Int4,Int4);
		gss_typ(const gss_typ&);	// copy constructor
		~gss_typ(){ Free(); }
	gss_typ& operator=(const gss_typ&);
	void    initialize(const ss_type,Int4,Int4,double,Int4,Int4);
	ss_type	FromArray(Int4, gsq_typ **,Int4,Int4,double,Int4,Int4,char *,a_type);
	Int4	*Counts( ){ assert(truedata!=NULL); return CntsSeqSet(truedata); }
	ss_type	FakeSqSet(){if(fakedata) return fakedata; else return MkFakeSqSet();}
	ss_type	TrueSqSet(){ if(truedata == NULL) return NULL; else return truedata; }
	e_type	FakeSeq(Int4);
	e_type	TrueSeq(Int4 s){ assert(truedata!=NULL);
		if(s < 1 || s > number) return NULL; else return SeqSetE(s,truedata); }
	Int4	TrueSite(Int4, Int4);
	Int4	TrueSubSite(Int4 sq, Int4 p)
			{ return (TrueSite(sq,p)-OffSetSeq(TrueSeq(sq))); }
	Int4	TruePos(Int4, Int4);
	Int4	TrueFullLength(Int4);
	Int4    RealToFake(Int4,Int4);
	Int4	NumDel(Int4 J){ assert(J > 0 && J <= number); 
			if(Gs[J]== NULL) return 0; else return Gs[J]->NumDel(); }
	Int4	NumIns(Int4 J){ assert(J > 0 && J <= number); 
			if(Gs[J]== NULL) return 0; else return Gs[J]->NumIns(); }
	Int4    Region(Int4,char *,Int4,Int4);
	void	Replace(Int4, gsq_typ *);
	gsq_typ	*Swap(Int4, gsq_typ *);
	void	RmFake(Int4);
	Int4	Insertion(Int4, UInt4);
	BooLean	Gapped( ){ return ((num_open+num_extend) > 0); }
	BooLean	IsDeleted(Int4,UInt4);
	Int4	InDels(Int4,UInt4,UInt4,Int4 *,Int4 *);
	a_type	Alphabet() { return SeqSetA(truedata); }
	Int4	MaxTrueSqLen( ){ return MaxSeqSeqSet(truedata); }
	Int4	MaxFakeSqLen( ){ return MaxSeqSeqSet(FakeSqSet()); }
	Int4	MinTrueSqLen( ){ return MinSeqSeqSet(truedata); }
	Int4	MinFakeSqLen( ){ return MinSeqSeqSet(FakeSqSet()); }
	Int4	AveTrueSqLen( ){ return AveSeqSeqSet(truedata); }
	void	PutFA(char *);	// put gapped sequences in fasta format.
	void	PutFA(FILE *);
	void	Put(FILE *);
	void	Put(FILE *fp, Int4 line_len, Int4 s){ assert(truedata!=NULL);
			if(Gs[s]) Gs[s]->Put(fp,line_len, SeqSetA(truedata)); }
	void	Put(FILE *,Int4);
	char    *Operation(Int4 s){ return Operation(0,s); }
	char    *Operation(FILE *fp,Int4 s);
	BooLean Identical(Int4,const gsq_typ&);  // Are sequences identical?
	BooLean NewAlign(char *,Int4,Int4,Int4); // New alignment? 
	Int4	NumSeq( ){ assert(truedata!=NULL); return NSeqsSeqSet(truedata); }
	double	IndelPenalty();
	double	IndelPenalty(Int4 s);
	double	IndelPenalty(UInt4,UInt4,UInt4);
	void	SetIndelPenalty(Int4, Int4);
	void	SetPerNats(double pn){ pernats=pn; }
	Int4	NumOpen(Int4);
	Int4	NumExtend(Int4);
	Int4    GapPenalty(Int4,Int4,idp_typ*);
	Int4    GapPenalty(Int4,Int4,Int4,idp_typ*);
	Int4	GapOpen(){ return gapopen; }
	Int4	GapExtend(){ return gapextend; }
	// Int4    GapCost(Int4,UInt4,unsigned short);
	void	InsertGap(Int4, UInt4,unsigned short);
	double	PerNats(){ return pernats; }
	Int4	LeftFlank(){ return leftflank; }
	Int4	RightFlank(){ return rightflank; }
	void	SetLeftFlank(Int4 x){ assert(x >= 0); leftflank=x; }
	void	SetRightFlank(Int4 x){ assert(x >= 0); rightflank=x; }
	Int4    OverHangN(Int4);
	Int4    OverHangC(Int4);
	gsq_typ *ForcedGetGSQ(Int4 s);
	gsq_typ	*GetGSQ(Int4 s){ 
		   assert(truedata != NULL); assert(s >= 1 && s <= number);
        	   if(Gs[s] == NULL) return NULL; else return Gs[s];
		}
private:
	void    Free();
	ss_type	MkFakeSqSet( );
	void	init( );
	void	copy(const gss_typ&);
	Int4	number,num_open,num_extend;
	ss_type	fakedata,truedata; // gapped and full sequence set...
	Int4	gapopen,gapextend;
	Int4	leftflank,rightflank; // flanking regions for domain sampling.
	double	pernats;
	gsq_typ	**Gs;
};

#endif


