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

#include "gsq_typ.h"

// WARNING: NEED TO DEAL WITH MASKED REGIONS OF SEQUENCES!!!
// XSeqPtr( ) function call???

void    gsq_typ::MkGapless(e_type E)
// Make a gapless gsq_typ...
{
        Int4            g,f,r,N;

        if(fakeE != NULL) Free();
        N=LenSeq(E); realE=E;
        fakeE=CopySeq(E);
        num_del=0; num_insert=0; num_open=0;
        f2r = new f2r_typ[N+3];
        indels = new indel_typ[N+3];
        for(r=0; r <= LenSeq(E); r++){ indels[r]=0; f2r[r]=r; }
}

void	gsq_typ::InsertGap(Int4 p, unsigned short gap, e_type E) 
// Create a sequence with an insertion ('-') of length 'gap'.
// WARNING: CHECK TO MAKE SURE DUMMY RESIDUES AT N-TERMINAL ARE OK.
// C-terminal extension is always zero -- so don't need to set.
{
	Int4		g,f,r,N;

	assert(p <= LenSeq(E));
	if(fakeE != NULL) Free();
	N=LenSeq(E)+gap; realE=E; 
	fakeE=CopySeqPlus(E,p,gap);  // inserts 'gap' dummy residues at p
	num_del=gap; num_insert=0; num_open=1;
	f2r = new f2r_typ[N+3]; 
	indels = new indel_typ[N+3];
	for(f=r=0; r <= p; ){ indels[f]=0; f2r[f]=r; r++; f++;} r--;
	for(g=0; g< gap; g++){ f2r[f]=r; indels[f]=gseq_del_mask; f++; } 
	while(r < LenSeq(E)){ r++; indels[f]=0; f2r[f]=r; f++;} f--;
	assert(f==N);
}

Int4	gsq_typ::GapPenalty(UInt4 start, idp_typ *idp)
{ return GapPenalty(start, 1, idp); }

Int4	gsq_typ::GapPenalty(UInt4 start, Int4 blk, idp_typ *idp)
{
	assert(fakeE != NULL); assert(start > 0);
	Int4	penalty=0,ins,end;
	Int4	i,j,s,k,len=idp->Length(blk);
	char	state=' ';

	end=idp->EndBlk(blk);
	for(j=idp->StartBlk(blk),i=start; j <= end; j++,i++){
	   if(indels[i]){
		if(deletion(i)){
		   if(state == 'd') { penalty += idp->DelExtend(j); }
		   else {
			penalty += idp->DelOpen(j)+idp->DelExtend(j); state = 'd';
		   }
		} ins = insert(i);
		if(ins){
		  penalty += idp->InsOpen(j);
		  penalty += ins*idp->InsExtend(j);
#if 1	// NEED TO REPLACE THIS WITH A MORE EFFICIENT ROUTINE.
		  Int4 *tmp = idp->SeqPenalties(SeqPtr(realE),LenSeq(realE));
		  if(tmp){
		    for(k=1,s=f2r[i]; k <= ins; k++,s++); penalty += tmp[s]; 
		    free(tmp);
		  }
#endif
		  state = 'i';
		}
	   } else state = ' ';
	} return penalty;
}

Int4	gsq_typ::InDels(UInt4 s,UInt4 e,Int4 *inso,Int4 *insx) 
// returns the number of deletions and sets inso == # insertions
// and insx == total residued inserted.
{
	assert(fakeE != NULL); assert(e <= LenSeq(fakeE)); assert(s > 0);
	register Int4 i,io,ix,d;
	for(io=ix=d=0,i=s; i < e; i++){
	   if(indels[i]){
	      if(deletion(i)) d++;
	      if(insert(i)){ io++; ix += insert(i); }
	   }
	} if(deletion(i)) d++;  // see if last residue deleted.
	*inso = io; *insx = ix;
	return d;
}

BooLean	gsq_typ::IsDeleted(UInt4 p) 
// Is there a deletion at position p in the sequence?
{
	assert(fakeE != NULL); assert(p <= LenSeq(fakeE));
	return is_deleted(p);
}

Int4	gsq_typ::Insertion(UInt4 p) 
// returns the length of the insertion at position p
{
	assert(fakeE != NULL); assert(p <= LenSeq(fakeE));
	return insert(p);
}

Int4	gsq_typ::CheckForInsertsAtEnds( ) 
// 
{
	Int4	i=0;
	for(UInt4 p = 0; p <= LenSeq(fakeE); p++){
		i += Insertion(p);
	}
	return (num_insert - i);
}

e_type	gsq_typ::InsertGap(Int4 p,unsigned short gap) 
// Create a sequence with an inserted gaps ('-') of length 'gap'.
// If there was an insertion event at this position then 
// 'uninsert' this residue instead.
{
	Int4		r,f,n,N,ni,nd,no;
        e_type          oldE,nE;         // new gapped subsequence
        f2r_typ		g,*n2r; 
	unsigned char	*nsq,*rsq;
	indel_typ	i,*ins;

	assert(fakeE != NULL); assert(p <= LenSeq(fakeE));
	// 1. Create new gsq structural components...
	N=LenSeq(fakeE)+gap; ni=num_insert; nd=num_del; no=num_open;
	n2r = new f2r_typ[N+3]; ins = new indel_typ[N+3];
	nE=CopySeqPlus(fakeE,p,gap);  // puts 'gap' # of '-' residues at p
	nsq=SeqPtr(nE); rsq=SeqPtr(realE);

	if(p==0){	// CASE: N-terminal gap (no insert penalty).
	  assert(!deletion(0));
	  // Retract residues from realE 1st...
	  i=insert(0);
	  for(g=0,r=f2r[1]-1,n=gap; i>0 && g<gap; i--,g++){
		   nsq[n]=rsq[r]; n2r[n]=r; ins[n]=0; r--; n--;
	  }
	  // AdjustOffsetSeq(g,nE);   // OLD Adjust offset for insert.
	  // Let's hope this fixes the problem without introducing more...
	  AdjustOffsetSeq(-g,nE);   // Adjust offset for insert.
	  if(g < gap) no++;
	  for(nd+=gap-g; g<gap; g++){ n2r[n]=0; ins[n]=gseq_del_mask; n--; }
	  ins[0]=i; n2r[0]=0; f=0; n=gap;
	} else { // Retract looped-out regions from the right.
	  for(f=n=-1; f < (Int4) p; ){ n++; f++; ins[n]=indels[f]; n2r[n]=f2r[f];}
	  g=0; r=f2r[p];
	  if(!(i=deletion(p))){	// Not just extending existing gap.
	     for(ins[p]=0; i>0 && g < gap; i--,g++){ // Retract from realE 1st.
		 n++; r++; nsq[n]=rsq[r]; n2r[n]=r; ins[n]=0;
	     } ins[n]=i;
	     // Code was core dumping here (below) due to negative extension assert.
	     // Fixed by setting negative extensions to zero...
	     AdjustCtermExtendSeq(-g,nE); // Subtract length of real added back.
	     if(p < LenSeq(fakeE)){		// Gap inside sequence.
	        if(g==0){			// No insertion at p and not...
		  if(!deletion(f+1)) no++; // next to existing gap.
	        } else if(g < gap){		// Swallowed up insertion and...
		  if(deletion(f+1)) no--; // next to existing gap.
		} else {	// (g==gap) 	// Avoided gap & 
		  if(i==0) no--;		// eliminated insertion open penalty.
		} ni-=g;
	     } else if(g<gap){ no++; }		// C-terminal gap (no insert penalty).
	  }					// Next: Gap ('-') residues inserted.
	  for(nd+=gap-g; g < gap; g++){ n++; n2r[n]=r; ins[n]=gseq_del_mask; }
	}
	while(f < LenSeq(fakeE)){ n++; f++; ins[n]=indels[f]; n2r[n]=f2r[f]; }
	assert(n==N && f==LenSeq(fakeE));
	oldE=fakeE; fakeE=NULL; Free( ); // Move new structures over to gsq.
	fakeE=nE; indels=ins; f2r=n2r; num_insert=ni; num_del=nd; num_open=no;
	return oldE;
}

BooLean	gsq_typ::Identical(const gsq_typ& gsq) 
// Checks to see whether these sequences are identically gapped.
{
	Int4	i,N;
	unsigned char	*seq1,*seq2;

	if(realE != gsq.realE) return FALSE;
	N=LenSeq(fakeE);
	if(N != LenSeq(gsq.fakeE)) return FALSE;
	if(OffSetSeq(fakeE) != OffSetSeq(gsq.fakeE)) return FALSE;
	if(CtermExtendSeq(fakeE) != CtermExtendSeq(gsq.fakeE)) return FALSE;
	seq1=SeqPtr(fakeE); seq2=SeqPtr(gsq.fakeE); 
	for(i=0; i <= N; i++){
		if(seq1[i] != seq2[i]) return FALSE;
		if(indels[i] != gsq.indels[i]) return FALSE;
	}
	return TRUE;
}

gsq_typ& gsq_typ::operator=(const gsq_typ& gsq)
// called for gsq_typ gsq2; gsq2=gsq;
{ if (this != &gsq) { Free(); copy(gsq); } return *this; }

void	gsq_typ::Free()
{
	if(fakeE != NULL) NilSeq(fakeE);
	if(indels != NULL) delete []indels;
	if(f2r != NULL) delete []f2r;
	init( );
}

void	gsq_typ::copy(const gsq_typ& gsq) 
// private function that assumes 'this' has been initialized and
// copies gsq to 'this'.
{
	Int4 N = LenSeq(gsq.fakeE);
	realE=gsq.realE; fakeE=CopySeq(gsq.fakeE);
	f2r = new f2r_typ[N+3];
	indels = new indel_typ[N+3];
	for(Int4 i=N; i>=0; i--){
		indels[i]=gsq.indels[i]; f2r[i]=gsq.f2r[i];
	}
	num_insert=gsq.num_insert; num_del=gsq.num_del; num_open=gsq.num_open;
}

#if 0
Int4	*gsq_typ::Recombine(gsq_typ& gsq1, gsq_typ& gsq2, Int4 *pos1, Int4 *pos2, 
	Int4 *config)
/*****************************************************************************
 Set *this == to the recombinant of gsq1 & gsq2 at crossover points xoverpt1
  and xoverpt2.  Set pos == 

          1      2            3        4
        -===-----===---------===------===---    aln1
        ....     ....      ......     ......
            \.../    \..../      \.../          (crossover points)
        -----===------===----===--===-------    aln2
             -1       -2     -3    -4

  Int4 *config = [ n, 1,-1, 2,-2, 3,-4, 4] n = 7 = number of recombinant blocks.
 *****************************************************************************/
// Create 
{
	Free();
	assert(gsq1.realE==gsq2.realE);

	Int4	b,s,s1,s2,s0,site,j,nblks=config[0],*pos,offsetX,extendX;
	char	last,state,id[MAX_SEQ_DEFLINE_LENG + 100];
	char	descript[MAX_SEQ_DEFLINE_LENG + 101];	// current sequence == 1 or 2.
	unsigned char *fsq,*osq;
        indel_typ *ins;
        f2r_typ	  *of2r;

	realE=gsq1.realE;
	fsq = new unsigned char [LenSeq(gsq1.realE) + gsq1.num_del + gsq2.num_del + 3];
	NEW(pos,nblks+2,Int4);
	if(config[1] < 0){
		offsetX = OffSetSeq(gsq2.fakeE);
	} else {
		offsetX = OffSetSeq(gsq1.fakeE);
	}
	for(s=s1=s2=0,b=1; b <=nblks; b++){
	   if(config[b] < 0) state=2; else state=1; 
	   switch(state){
	    case 1: 
		osq=SeqPtr(gsq1.fakeE); s0=pos1[b-1];
		of2r=gsq1.f2r; ins=gsq1.indels;
		site=pos1[b];
	     break;
	    case 2: 
		osq=SeqPtr(gsq2.fakeE); s0=pos2[-b-1];
		of2r=gsq2.f2r; ins=gsq2.indels;
		site=pos2[-b];
	     break;
	    default: assert(state==1 || state ==2); break;
	   } // Now fill in sequence.
	   assert(FakeToReal(xopt1)==FakeToReal(xopt2));
	   for( ; s0 <= site
	   last=state;
        }
	if(last==2) extendX += CtermExtendSeq(gsq2.fakeE);
	else extendX += CtermExtendSeq(gsq2.fakeE);

	StrSeqID(id, DEFAULT_INFO_LINE_LEN_SEQ, realE); 
	StrSeqDescript(descript,DEFAULT_INFO_LINE_LEN_SEQ, realE);
	fakeE=MakeSeq(id,descript,offsetX,extendX,s,fsq);
	return pos;
}
#endif

Int4    gsq_typ::OverHangC( )
{ return (CtermExtendSeq(fakeE) - CtermExtendSeq(realE)); }	// this may be wrong!!!
// { return (CtermExtendSeq(realE) - CtermExtendSeq(fakeE)); }  // this instead?

Int4    gsq_typ::OverHangN( )
{ return (OffSetSeq(fakeE) - OffSetSeq(realE)); }

char	*gsq_typ::Operation(FILE *fp, a_type A)
// Return the series of operations needed to construct this seq using
// the initialize operation.
// WARNING: this is currently designed for single block models only!!!
// WARNING: THIS NEEDS MORE WORK!!! 
{
	unsigned char *rseq,*fseq;
	Int4	i,j,k,length,TrueLen;
	UInt4 ptr;
	char	c1,c2,operation,*Operation;

	assert(fakeE != NULL);
	TrueLen=LenSeq(realE);
	length=LenSeq(fakeE); rseq=SeqPtr(realE); fseq=SeqPtr(fakeE);
	NEW(Operation,LenSeq(fakeE) + LenSeq(realE) +3,char);
	if(fp) fprintf(fp,"     r f   ins f2r oper\n");
	Operation[0]='E'; ptr=1;
#if 0
	// indels at start...
	k = insert(0);
	for(ptr=1,j=f2r[0]; k > 0; k--){
                j++; c2 = tolower(AlphaChar(rseq[j],A));
		Operation[ptr]='i'; ptr++;
		fprintf(fp,"%3d: %c .    :   - %c\n",0,c2,operation);
	} 
#else	// add N-terminal extension
	// for(i = OverHangN( ); i > 0;  i--) Operation[ptr]='m'; ptr++; }
#endif
        for(i=1; i<=length; i++){
	    if(!deletion(i)){ // no deletion here
                c2 = AlphaChar(rseq[f2r[i]],A);
                c1 = AlphaChar(fseq[i],A);
		if(f2r[i] > 0 && f2r[i] <= TrueLen) operation='m';
		else operation='d';
		if(i==1)operation=toupper(operation);
		Operation[ptr]=operation; ptr++;
		if(fp) fprintf(fp,"%3d: %c %c  %3d %3d %c\n",i,c2,c1,
				insert(i),f2r[i],operation);
	    } else {				// deletion ='-'.
                c1 = AlphaChar(fseq[i],A);
		if(f2r[i] > 0 && f2r[i] <= TrueLen) operation='d';
		else operation='d';
		if(i==1)operation=toupper(operation);
		Operation[ptr]=operation; ptr++;
		if(fp) fprintf(fp,"%3d: . %c  %3d %3d %c\n",i,c1,
			insert(i), f2r[i],operation);
	    }
	    k = insert(i);
	    for(j=f2r[i]; k > 0; k--){
                j++; c2 = tolower(AlphaChar(rseq[j],A));
		operation='I';	// always assume insertions are within block for now...
		Operation[ptr]=operation; ptr++;
		if(fp) fprintf(fp,"%3d: %c .    :   - %c\n",i,c2,operation);
	    } 
        }
#if 1	// BUG FIX for insertions on ends...Do a better job latter...
	while(Operation[ptr-1] =='I'){ ptr--; } // back up over 'I' operations
#else
	// for(i = OverHangC( ); i > 0;  i--) Operation[ptr]='m'; ptr++; }
#endif
	Operation[ptr]='E'; ptr++; Operation[ptr]=0;
	if(fp){
	  fprintf(fp,"num_insert=%d; num_del=%d; num_open=%d.\n",
		num_insert,num_del,num_open);
	  fprintf(fp,"#### offset real = %d; offset fake = %d\n",
                OffSetSeq(realE),OffSetSeq(fakeE));
	  PutSeq(fp,fakeE,A);
	}
	return Operation;
}

Int4	gsq_typ::Region(char *rtnaln, Int4 site, Int4 len_align, a_type A)
// site is site in fake seq; want corresponding sites in original sequence.
{
	Int4		i,j,k,p;
	char		r;
	indel_typ	ins;
	unsigned char	*rseq=SeqPtr(realE);

	assert(fakeE != NULL); assert(site > 0 && site <= LenSeq(fakeE));
	for(j=k=0,i=site; j < len_align; ){
	   ins=insert(i); p = f2r[i];
	   if(deletion(i)) { rtnaln[k] = '-'; j++; k++; i++; }
	   else { rtnaln[k] = AlphaChar(rseq[p],A); j++; k++; i++; }
	   if(ins > 0 && j < len_align){ // insertion must be within block.
		while(ins > 0){
		   p++; r = AlphaChar(rseq[p],A);
		   rtnaln[k] = tolower(r); k++; ins--;
		}
	   }
	} i--;
	rtnaln[k] = 0; 
#if 0	// debug...
	if(1 || site==5){ std::cerr << "DEBUG\n"; 
		fprintf(stderr,"site = %d; len_align = %d; %s\n",site,len_align,rtnaln); 
	}
#endif
#if 1
	if(i > LenSeq(fakeE)){ 
	  Put(stderr,A); PutFake(stderr,A);
	  fprintf(stderr,"site = %d; len_align = %d\n",site,len_align);
	  fprintf(stderr,"%s: i = %d\n",rtnaln,i); 
	  fprintf(stderr,"len_align = %d (%s)\n",k,rtnaln); fflush(stderr);
	  assert(i <= LenSeq(fakeE));
	}
#endif
	return k;
}

Int4	gsq_typ::RealToFake(Int4 i)
// CAUTION: this in not efficient!
// get the position in fake seq corresponding to position i in realE.
// return 0 if nothing corresponds...
{ for(Int4 j=1; j<=LenSeq(fakeE); j++){ if(f2r[j] == i) return j; } return 0; }

Int4	gsq_typ::FakeToReal(Int4 i)
// Int4 gsq_typ::TrueSite(Int4 i)
// return the location in source sequence (realE) of residue i in fakeE 
//        (2)     10           15 (19)   25   30
// fakeE: ----+----+--........--+--xx+----+----+-
//         |(2)       (absent)     |(26)	      f2r[2]=2; f2r[19]=25;
// realE: ----+----+----+----+----+..----+----+--
//        1   5   10   15   20   25 26  30   35
{ assert(i >= 0 && i <= LenSeq(fakeE)); return (f2r[i]); }

BooLean	gsq_typ::NewAlign(Int4 leftflank,Int4 rightflank,char *operation,
	Int4 trace_length,Int4 start0)
/****************************************************************************
  See whether operation array corresponds to a new alignment for qsq.
  trace = "EDdmmmmmmmmmmmiiiMmmmmmmmmmmmmmmmmmmmmiiMmmmmmmmmmmE"
 ****************************************************************************/
{
	Int4		o,j,k,J0,flankL,flankR,start;
	char		state;

	if(fakeE == NULL) return TRUE;
	assert(realE != NULL);
	if(operation[1]=='i' || operation[0]!='E' || operation[trace_length-2]=='i'
		|| operation[trace_length-1]!='E'){
		fprintf(stderr,"operation = %s\n trace_length=%d(%d)\nstart=%d\n",
			operation,trace_length,strlen(operation),start0);
		print_error("gsq_typ::IsMatch( ) input error");
	}
	flankL=MINIMUM(Int4,leftflank,start0-1); 
	start=start0-flankL-1;	// start relative to real subsequence.
	for(J0=0,j=start,k=1; k<= flankL; k++){
		J0++; j++; 
		if(indels[J0]!=0 || f2r[J0]!=j) return TRUE;
	}
// std::cerr << "\nE"; 
        for(o=1,state='E'; operation[o] != 'E'; o++){
// std::cerr << operation[o]; 
            switch(operation[o]){
               case 'M': case 'm': J0++; j++;
		if(IsDeleted(J0) || f2r[J0] != j) return TRUE; break;
               case 'D': case 'd': // deletion in sequence relative to profile.
                  J0++; 
		  if(!deletion(J0) || f2r[J0] != j) return TRUE;
		break;
               case 'i': J0++; j++; 
		if(indels[J0] || f2r[J0] != j) return TRUE;
		break; // insert between profile blocks;
               case 'I': // Insert within a profile block; delete from seq.
		 {
		   assert(state != 'I'); 
		   indel_typ	i;
		   for(i=0; operation[o]=='I'; ){ j++; i++; o++; }
		   o--; assert(i < gseq_ins_mask);
		   if(insert(J0) != i) return TRUE;
		 } break;
               default: fprintf(stderr,"operations = %s\n",operation);
                 print_error("gsq_typ::IsMatch( ): input error"); break;
            }  state=operation[o];
        } // 3. Add rightflank region.
// std::cerr << "E\n"; 
        flankR = MINIMUM(Int4,rightflank,LenSeq(realE)-j);
        for(Int4 i=1; i < flankR; i++){
            J0++; j++;
	    if(indels[J0] || f2r[J0]!=j){
// fprintf(stderr,"indels(%d)=%d; f2r=%d; j=%d; flankR=%d; i=%d\n",J0,indels[J0],f2r[J0],j,flankR,i);
		return TRUE;
	    }
        }
	if(flankR > 0){ J0++; j++;
	  if(deletion(J0) || f2r[J0]!=j){
// fprintf(stderr,"indels(%d)=%d; f2r=%d; j=%d; flankR=%d; i=%d\n",J0,indels[J0],f2r[J0],j,flankR,i);
		return TRUE;
	  }
	}
	Int4 x=(LenSeq(realE)-j);
	assert(x < gseq_ins_mask);
        if(insert(J0) != (indel_typ) x) return TRUE;
	if(CtermExtendSeq(fakeE) != CtermExtendSeq(realE) + x) return TRUE;
	return FALSE;
}

void	gsq_typ::FindIndels(Int4 nblks, Int4 *start, Int4 *len, 
	UInt4 &Nmm, UInt4 &Nmi, UInt4 &Nmd, UInt4 &Nii,
	UInt4 &Nim, UInt4 &Ndd, UInt4 &Ndm,
	UInt4 &Ndi,UInt4 &Nid,UInt4 &Nsd,
                UInt4 &Nsm)
// returns the length of the insertion at position p
{
   assert(fakeE != NULL); 
   Int4	n,end;
   char	state;
   for(Int4 bk = 1 ; bk <= nblks; bk++){
	end=start[bk] + len[bk] -1;
	assert(end <= LenSeq(fakeE));
   	state='S';
	for(Int4 i=start[bk]; i <= end; i++){
           if(!deletion(i)){	// next state == 'M';
		switch(state){
		  case 'S': Nsm++; break;
		  case 'M': Nmm++; break;
		  case 'I': Nim++; break;
		  case 'D': Ndm++; break;
		  default: 
		     fprintf(stderr,"state = %c; sq=%s\n",state,SeqKey(realE));
		     print_error("gsq_typ::FindIndels 'M' input error"); 
		  break;
		} state='M';
	   } else {		// next state == 'D'
		  switch(state){
		   case 'S': Nsd++; break;
		   case 'M': Nmd++; break;
		   case 'D': Ndd++; break;
		   case 'I': Nid++; break;	// I->D transitions sometimes show up.
		   default: 
		     fprintf(stderr,"state = %c; sq=%s\n",state,SeqKey(realE));
		     print_error("gsq_typ::FindIndels 'D' input error"); 
		   break;
		  } state='D';
	   }
	   if(i != end && (n=insert(i))){	// is there a following insertion?
	   // Don't worry about insertions at ends...
		switch(state){
		   // e.g., "TaFfsv..." = previous insert followed by a match and more inserts.
		   // case 'I': Nim++;  Nmi++; Nii+=n-1; break;
			// case 'I' should not occur...
		   // case 'S': break; // insert or match...
		   case 'M': Nmi++; Nii+=n-1; break;
		   case 'D': Ndi++; Nii+=n-1; break;
		   default: 
		     fprintf(stderr,"state = %c; n=%d; bk = %d; i = %d; sq=%s\n",
			state,n,bk,i,SeqKey(realE));
		     for(Int4 j=start[bk]; j <= end; j++){
			 fprintf(stderr,"indels[%d]=%d\n",j,indels[j]);
		     }
		     print_error("gsq_typ::FindIndels 'I' input error"); 
		   break;
		} state='I';
	   }
	}
   } return;
}

char	*gsq_typ::Operation(Int4 nblk,Int4 *start, Int4 *blk_len)
// return operations need to create the gsq.
// 'i' == flank or gap; 'M' = start motif; 'D' = delete at start; 'I' = insert within motif.
{
	Int4	i,j,k,b,site,s,true_len=LenSeq(realE),fake_len=LenSeq(fakeE);
	UInt4	ptr;
	char    c,*Operation,state='F';
	NEW(Operation,LenSeq(fakeE) + LenSeq(realE) +3,char);
	Operation[0]='E'; ptr=1;
	// assert(blk_len[nblk+1]==0);

	for(s=1; s <= this->OverHangN(); s++){ Operation[ptr]='i'; ptr++; } 
	for(state='G',site=LenSeq(fakeE),b=0,i=1; i<=fake_len; i++,site++){
	   if(b < nblk && i == start[b+1]){ state='M'; b++; site=1; } else state='m';
	   if(site > blk_len[b]) state = 'G'; 
	   if(!deletion(i)){  // no deletion at this site.
                // if(f2r[i] > 0 && f2r[i] <= true_len) c='m'; else c='d';
                assert(f2r[i] > 0 && f2r[i] <= true_len);
		c='m';
		if(state == 'M') c=toupper(c);
		else if(state == 'G') c='i';
		Operation[ptr]=c; ptr++; 
	   } else {
                // assert(f2r[i] > 0 && f2r[i] <= true_len);
                // if(f2r[i] > 0 && f2r[i] <= true_len) c='d'; else c='d';
		c='d';
		if(state == 'M') c=toupper(c);
		// Operation[ptr]=c; ptr++;
		if(state !='G'){ Operation[ptr]=c; ptr++; }     // ignore deletions between blocks.
	   } 
	   // if(site == blk_len[b] && c=='d') state = 'G';
	   if(site == blk_len[b]) state = 'G';
	   k = insert(i);
           for(j=f2r[i]; k > 0; k--){ 	// insertion here.
		j++; c='I';
		if(state == 'G') c = tolower(c); 
		Operation[ptr]=c; ptr++;
           } 
        } Operation[ptr]='E'; ptr++;
	return Operation;
}

char    *gsq_typ::FromSitesOperation(Int4 nBlks,Int4 *len,Int4 *pos,e_type E)
{
        Int4    i,j,blk,s;
        char    *op; NEW(op,LenSeq(E)+10,char); op[0]='E';
        for(s=0,i=1,blk=0,j=1; i <= LenSeq(E); s++,i++){
           if(j <= nBlks && i == pos[j]){ blk++; s=1; op[i]='M'; j++; }
	   else if(blk > 0 && s <= len[blk]){ op[i]='m'; }
	   else op[i]='i';
        } j = pos[nBlks] + len[nBlks] -1; // j == end of last block.
	assert(j <= LenSeq(E));  // j not beyond the end of Seq E.
	op[i]='E'; i++; op[i]=0;
        return op;
}

