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

#include "gss_typ.h"

void	gss_typ::Free()
{
     if(fakedata != NULL) free(NilSeqSetRtnSeqs(fakedata));
     if(Gs != NULL){	// FakeSeqs are deleted.
	for(Int4 s=1; s<=number; s++) if(Gs[s] != NULL) delete []Gs[s];
	delete []Gs;
     }
     init( );
}

Int4    gss_typ::NumOpen(Int4 s)
{
	assert(truedata != NULL && s >= 1 && s <= number);
	if(Gs[s] == NULL){ return 0; } else return Gs[s]->NumOpen();
}

Int4    gss_typ::Insertion(Int4 s, UInt4 i)
{
	assert(truedata != NULL && s >= 1 && s <= number);
	if(Gs[s] == NULL){ return 0; } else return Gs[s]->Insertion(i);
}

Int4    gss_typ::NumExtend(Int4 s)
{
	assert(truedata != NULL && s >= 1 && s <= number);
	if(Gs[s] == NULL){ return 0; } else return Gs[s]->NumExtend();
}

BooLean	gss_typ::NewAlign(char *operation,Int4 trace_length,Int4 start,Int4 s)
{
	assert(truedata != NULL && s >= 1 && s <= number);
	if(Gs[s] == NULL) return TRUE; 
	else return Gs[s]->NewAlign(leftflank,rightflank,operation,trace_length,start);
}

BooLean	gss_typ::Identical(Int4 s,const gsq_typ& gsq)
{
	assert(truedata != NULL && s >= 1 && s <= number);
	if(Gs[s] == NULL){
	  if(gsq.Gapped()) return FALSE; else return TRUE;
	} else return Gs[s]->Identical(gsq);
}

void	gss_typ::PutFA(char *filename)
// Overload to allow specification of filename
{
	assert(truedata != NULL); 
	FILE *fp = open_file(filename,"","w"); PutFA(fp); fclose(fp); 
}

void	gss_typ::PutFA(FILE *fp)
{
	Int4	s;
	assert(truedata != NULL);
	if(fakedata != NULL) PutSeqSetEs(fp,fakedata);
	else for(s=1; s<=number; s++) PutSeq(fp,FakeSeq(s),SeqSetA(truedata));
}

char	*gss_typ::Operation(FILE *fp, Int4 s)
{
	assert(truedata != NULL && s >= 1 && s <= number);
	if(Gs[s] == NULL){
	   print_error("Operation input error");
	   return 0;
	} else { return Gs[s]->Operation(fp,SeqSetA(truedata)); }
}

void    gss_typ::Put(FILE *fp,Int4 s)
{
	assert(truedata != NULL && s >= 1 && s <= number);
	if(Gs[s] == NULL) PutSeq(fp,FakeSeq(s),SeqSetA(truedata));
	else { 
		Gs[s]->Put(fp,SeqSetA(truedata)); 
		Gs[s]->Put(fp,60,SeqSetA(truedata));
	}
}

void    gss_typ::Put(FILE *fp)
{
	assert(truedata != NULL);
	for(Int4 s=1; s<=number; s++){
		if(Gs[s] == NULL) PutSeq(fp,FakeSeq(s),SeqSetA(truedata));
		else Gs[s]->Put(fp,60,SeqSetA(truedata));
	}
}

e_type  gss_typ::FakeSeq(Int4 s)
{
	assert(truedata != NULL && s >= 1 && s <= number);
	if(Gs[s] == NULL) return SeqSetE(s,truedata);
	else return Gs[s]->FakeSeq( );
}

Int4    gss_typ::TrueFullLength(Int4 s)
{
	e_type	E;
	Int4	len;
	assert(truedata != NULL); assert(s >= 1 && s <= number);
	E = SeqSetE(s,truedata);
	if(Gs[s] == NULL){
		len = LenSeq(E);
	} else {
		Int4 start = Gs[s]->FakeToReal(1);
		Int4 end = LenSeq(Gs[s]->FakeSeq( ));
		end = Gs[s]->FakeToReal(end);
		len = end - start + 1;
	}
	return (len + OffSetSeq(E) + CtermExtendSeq(E));
}

Int4    gss_typ::TruePos(Int4 s, Int4 site)
// get the position in true sequence corresponding to site in fake seq.
{
	assert(truedata != NULL); assert(s >= 1 && s <= number);
	if(Gs[s] == NULL) return site; else return Gs[s]->FakeToReal(site);
}

void    gss_typ::InsertGap(Int4 s, UInt4 p,unsigned short gap)
// insert a gap at p in sequence s.
// WARNING: NEED TO ALLOW SEG'ED SEQUENCES TO GET GAP AS WELL!!!
{
	Int4	o,i,d,x;
	e_type	oldE,oldE2;
	assert(truedata != NULL); assert(s >= 1 && s <= number);
	if(Gs[s] == NULL){  // gsq == trueSq.
	   gsq_typ *gsq=new gsq_typ[1];
	   gsq->InsertGap(p,gap,SeqSetE(s,truedata));
	   Replace(s,gsq);	// this routine updates fakedata and gapspenalty.
	} else {
// TEMPORARY FIX TILL I GET InsGapCost() WORKING RIGHT...  
       	   num_open-=Gs[s]->NumOpen(); num_extend-=Gs[s]->NumExtend();
// o=Gs[s]->InsGapCost(p,gap,i,d);
// num_open+=o; num_extend+=i+d;
// if(s==20) Gs[s]->Put(stderr,SeqSetA(truedata));
	   oldE=Gs[s]->InsertGap(p,gap);
// if(s==20) fprintf(stderr,"p = %d; gap = %d\n",p,gap);
// if(s==20) Gs[s]->Put(stderr,SeqSetA(truedata));
       	   num_open+=Gs[s]->NumOpen(); num_extend+=Gs[s]->NumExtend();
	   if(fakedata != NULL){
	     oldE2 = ReplaceSeqSet(s,Gs[s]->FakeSeq( ),fakedata);
	     SetIndelPenaltySeqSet(IndelPenalty( ), fakedata);
	     assert(oldE==oldE2);
	   }
	   NilSeq(oldE);
	}
}

double  gss_typ::IndelPenalty( )
{ return (double) (num_open*gapopen + num_extend*gapextend)/pernats; }

double  gss_typ::IndelPenalty(Int4 s)
{
  if(Gs[s] == NULL) return 0.0;
  return (double) (Gs[s]->NumOpen()*gapopen + Gs[s]->NumExtend()*gapextend)/pernats;	
}

double  gss_typ::IndelPenalty(UInt4 del, UInt4 inso, UInt4 insx)
{ 
	double log_odds=(double)(inso*gapopen+insx*gapextend);  // insertion penalty.
	log_odds += (double)(del*(gapopen+gapextend))/2.0; // deletion penalty.
	log_odds /= pernats;	// express in nats.
	log_odds *= 0.43429448; // make base of log == 10;
	return log_odds;
}

Int4	gss_typ::GapPenalty(Int4 s, Int4 start, idp_typ *idp)
{ return GapPenalty(s, start, 1,idp); }

Int4	gss_typ::GapPenalty(Int4 s, Int4 start, Int4 blk, idp_typ *idp)
{
	assert(truedata != NULL); assert(s >= 1 && s <= number);
	if(Gs[s] == NULL) return 0; else return Gs[s]->GapPenalty(start,blk,idp);
}

Int4	gss_typ::InDels(Int4 s, UInt4 st,UInt4 e, Int4 *inso, Int4 *insx)
// returns the number of deletions and sets inso == # insertions
// and insx == total residued inserted.
{
	assert(truedata != NULL); assert(s >= 1 && s <= number);
	if(Gs[s] == NULL) { *inso=*insx=0; return 0; }
	else return Gs[s]->InDels(st,e,inso,insx);
}

BooLean	gss_typ::IsDeleted(Int4 s,UInt4 i)
// Is there a deletion at position i in sequence s?
{
	assert(truedata != NULL); assert(s >= 1 && s <= number);
	if(Gs[s] == NULL) return FALSE; else return Gs[s]->IsDeleted(i);
}

#if 0
Int4	gss_typ::GapCost(Int4 s,UInt4 p,unsigned short gap)
// Get the gap opening and extension penalty to insert a gap at p in s.
{
	Int4	o,x;
	assert(truedata != NULL); assert(s >= 1 && s <= number);
	if(Gs[s] == NULL){  // gsq == trueSq.
		// if(p == 0 || p == LenSeq(SeqSetE(s,truedata)) // ????
		return (gapopen + gap*gapextend); 
	} else { Gs[s]->InsGapCost(p,gap,o,x); return (o*gapopen+x*gapextend); }
}
#endif

Int4    gss_typ::RealToFake(Int4 sq, Int4 site)
// get the actual residue number of site in (sub)sequence sq.
{
	assert(truedata != NULL); assert(sq >= 1 && sq <= number);
	if(Gs[sq] == NULL) return (site - OffSetSeq(SeqSetE(sq,truedata)));
	else return (Gs[sq]->RealToFake(site - OffSetSeq(SeqSetE(sq,truedata))));
}

Int4    gss_typ::TrueSite(Int4 s, Int4 site)
// get the actual residue number of site in (sub)sequence s.
{
	assert(truedata != NULL); assert(s >= 1 && s <= number);
	if(Gs[s] == NULL) return (site + OffSetSeq(SeqSetE(s,truedata)));
	else return (Gs[s]->FakeToReal(site) + OffSetSeq(SeqSetE(s,truedata)));
}

Int4    gss_typ::OverHangC(Int4 s)
{
	assert(truedata != NULL); assert(s >= 1 && s <= number);
	if(Gs[s] == NULL) return 0;
	else return Gs[s]->OverHangC( );
	// else return (CtermExtendSeq(Gs[s]->FakeSeq( )) - CtermExtendSeq(SeqSetE(s,truedata)));
}

Int4    gss_typ::OverHangN(Int4 s)
{
	assert(truedata != NULL); assert(s >= 1 && s <= number);
	if(Gs[s] == NULL) return 0;
	else return Gs[s]->OverHangN( );
	// else return (OffSetSeq(Gs[s]->FakeSeq( ))-OffSetSeq(SeqSetE(s,truedata)));
}

Int4    gss_typ::Region(Int4 s,char *rtnaln, Int4 site, Int4 len_align)
{
	assert(truedata != NULL); assert(s >= 1 && s <= number);
	a_type	A=SeqSetA(truedata);
	if(Gs[s] == NULL){
	   Int4	i,j;
	   unsigned char *oseq=SeqPtr(SeqSetE(s,truedata));
	   for(i=site,j=0; j < len_align; j++,i++) rtnaln[j]=AlphaChar(oseq[i],A);
	   rtnaln[j]=0; return len_align;
	} else return Gs[s]->Region(rtnaln,site,len_align,A);
}

void	gss_typ::init( )
{
	fakedata=NULL; Gs=NULL; 
	number=gapopen=gapextend=0; pernats=1000.0;
	num_open=num_extend=0;
}

gss_typ::gss_typ(ss_type seqset,Int4 open,Int4 extend,double Pernats,Int4 left,Int4 right)
{ init(); truedata=NULL; initialize(seqset,open,extend,Pernats,left,right); }

void	gss_typ::initialize(const ss_type seqset,Int4 open,Int4 extend,double Pernats,
	Int4 left,Int4 right)
{
	assert(seqset!=NULL);
	if(truedata!=NULL) Free( ); 
	truedata=seqset; number=NSeqsSeqSet(truedata); 
	gapopen=open;  gapextend=extend; pernats=Pernats;
	leftflank=left; rightflank=right;
	Gs = new gsq_typ*[number+1];	
	for(Int4 s=0; s <= number; s++) Gs[s]=NULL;
}

ss_type	gss_typ::FromArray(Int4 N, gsq_typ **gsq,Int4 open,Int4 extend,
	double Pernats,Int4 left,Int4 right,char *name,a_type A)
// De Novo creation of gss -- need to run lngamma[N][b]?
{
   e_type *Elist;
   Free( );
   number=N; Gs = gsq; leftflank=left; rightflank=right;
   gapopen=open; gapextend=extend; pernats=Pernats;
   num_open=num_extend=0;
   NEW(Elist,number+3,e_type);
   for(Int4 s=1; s<=number; s++) {
     if(Gs[s] != NULL){ 
       	num_open+=Gs[s]->NumOpen(); num_extend+=Gs[s]->NumExtend();
	Elist[s] = Gs[s]->TrueSeq(); 
     } else print_error("gss_typ::FromArray( ): input error");
   }
   truedata=Array2SeqSet(Elist,number,name,A); // Elist "absorbed"
   return truedata;
}

gss_typ::gss_typ(const gss_typ& gss) { init( ); truedata=NULL; copy(gss); } 
// called for 'gss_typ gss2=gss;' or 'gss_typ gss2(gss);'.

gss_typ& gss_typ::operator=(const gss_typ& gss)
// this is called for gss_typ gsscopy; gsscopy=gss;  
{ if(this != &gss) { if(truedata!=NULL) Free( ); copy(gss); } return *this; }

void	gss_typ::copy(const gss_typ& gss)
// WARNING: private function that assumes 'this' has been initialized.
// This needs to be done before coping gss to 'this'.
{
   initialize(gss.truedata,gss.gapopen,gss.gapextend,gss.pernats,
		gss.leftflank,gss.rightflank); 
   num_open=num_extend=0;
   for(Int4 s=1; s<=number; s++) {
     if(gss.Gs[s] != NULL){ 
	Gs[s]=new gsq_typ[1]; *Gs[s]=*(gss.Gs[s]); 
       	num_open+=Gs[s]->NumOpen(); num_extend+=Gs[s]->NumExtend();
     }
   }
   assert(num_open==gss.num_open && num_extend==gss.num_extend);
   if(gss.fakedata!=NULL) MkFakeSqSet( );
}

ss_type	gss_typ::MkFakeSqSet( )
{
   e_type	*Elist;
   Int4		s;

   assert(truedata != NULL && fakedata==NULL);
   for(s=1; s<=number; s++) if(Gs[s] != NULL) break;
   if(s > number) return truedata;  // if no gapped sequences....
   NEW(Elist,number+3,e_type);
   for(s=1; s<=number; s++) {
	if(Gs[s] == NULL){ Elist[s] = SeqSetE(s,truedata); }
	else { Elist[s] = Gs[s]->FakeSeq(); }
   }
   fakedata=Array2SeqSet(Elist,number,NameSeqSet(truedata),
			SeqSetA(truedata)); // Elist absorbed 
   SetIndelPenaltySeqSet(IndelPenalty( ), fakedata);
   return fakedata;
}

void    gss_typ::SetIndelPenalty(Int4 open, Int4 extend)
{
	gapopen=open; gapextend=extend;
	if(fakedata!= NULL) SetIndelPenaltySeqSet(IndelPenalty( ), fakedata);
}

void	gss_typ::Replace(Int4 s, gsq_typ *gsq)
// replace sth sequence with gsq; return old e_type sequence.
{
   e_type	E,oldE,oldE2;
   gsq_typ	*gsqold=NULL;

   assert(truedata != NULL); assert(s >= 1 && s <= number);
// fprintf(stderr,"indel penalty=%g (%d open; %d extend)\n",IndelPenalty(),num_open,num_extend);
   if(Gs[s] == NULL){
	Gs[s] = gsq;
	if(fakedata != NULL){ oldE=SeqSetE(s,truedata); E=gsq->FakeSeq(); }
   } else {
	gsqold = Gs[s]; Gs[s] = gsq;
        num_open-=gsqold->NumOpen( ); num_extend-=gsqold->NumExtend( ); 
	if(fakedata != NULL){ oldE=gsqold->FakeSeq(); E=gsq->FakeSeq(); }
   }
   num_open+=gsq->NumOpen(); num_extend+=gsq->NumExtend();
   if(fakedata != NULL){
        oldE2 = ReplaceSeqSet(s,E,fakedata);
        assert(oldE == oldE2);
// fprintf(stderr,"indel penalty=%g (%d open; %d extend)\n",
//	IndelPenalty(),num_open,num_extend);
	SetIndelPenaltySeqSet(IndelPenalty( ), fakedata);
   }
   if(gsqold!= NULL) delete []gsqold; // destroys oldE.
}

gsq_typ *gss_typ::Swap(Int4 s, gsq_typ *gsq)
// replace sth sequence with gsq; return old qsq.
{
   e_type	E,oldE,oldE2;
   gsq_typ	*gsqold=0;

   assert(truedata != NULL); assert(s >= 1 && s <= number);
// fprintf(stderr,"indel penalty=%g (%d open; %d extend)\n",IndelPenalty(),num_open,num_extend);
   if(Gs[s] == NULL){
	Gs[s] = gsq;
	if(fakedata != NULL){ oldE=SeqSetE(s,truedata); E=gsq->FakeSeq(); }
   } else {
	gsqold = Gs[s]; Gs[s] = gsq;
        num_open-=gsqold->NumOpen( ); num_extend-=gsqold->NumExtend( ); 
	if(fakedata != NULL){ oldE=gsqold->FakeSeq(); E=gsq->FakeSeq(); }
   }
   num_open+=gsq->NumOpen(); num_extend+=gsq->NumExtend();
   if(fakedata != NULL){
        oldE2 = ReplaceSeqSet(s,E,fakedata);
        assert(oldE == oldE2);
// fprintf(stderr,"indel penalty=%g (%d open; %d extend)\n",
//	IndelPenalty(),num_open,num_extend);
	SetIndelPenaltySeqSet(IndelPenalty( ), fakedata);
   } return gsqold; 
}

void	gss_typ::RmFake(Int4 s)
// remove sth sequence and set to null;
{

   assert(truedata != NULL); assert(s >= 1 && s <= number);
   if(Gs[s] == NULL){ // then already removed; return;
	return;
   } else {
      e_type	E,oldE,oldE2;
      gsq_typ	*gsqold=NULL;
      gsqold = Gs[s]; Gs[s] = NULL;
      num_open-=gsqold->NumOpen( ); num_extend-=gsqold->NumExtend( ); 
      if(fakedata != NULL){ 
	oldE=gsqold->FakeSeq(); E=TrueSeq(s);
        oldE2 = ReplaceSeqSet(s,E,fakedata);
        assert(oldE == oldE2);
	SetIndelPenaltySeqSet(IndelPenalty( ),fakedata);
      }
      if(gsqold!= NULL) delete []gsqold; // destroys oldE.
   }
}

gsq_typ *gss_typ::ForcedGetGSQ(Int4 s)
{
        assert(truedata != NULL); assert(s >= 1 && s <= number);
        if(Gs[s] == NULL){  // gsq == trueSq.
           gsq_typ *gsq=new gsq_typ[1];
           gsq->MkGapless(SeqSetE(s,truedata));
           Replace(s,gsq);
        } return Gs[s];
}

