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

#include "cma_gmb.h"

#if 0
void	PutAlnSeqsCMSA(FILE *fp,cma_typ cma)
{
	Int4	sq,NN=NumSeqsCMSA(cma);
	set_typ	Set=MakeSet(NN+5);  ClearSet(Set);
	for(sq=1; sq <= NN; sq++) if(SeqIsInCMSA(sq,cma)) AddSet(sq,Set);
	PutInSetCMSA(fp,Set,cma);  NilSet(Set);
}
#endif

Int4	MoveAlignmentCMSA(set_typ Used, cma_typ From, cma_typ To)
// Move the sequences within the 'Used' set from the first to the second alignment.
{
	Int4	n,i,j,blk,sq,*sites;
	gsq_typ	*gsq,*gsq0;
	assert(nBlksCMSA(From) == nBlksCMSA(To));
	for(blk=1; blk <= nBlksCMSA(To); blk++) assert(LengthCMSA(blk,From) == LengthCMSA(blk,To));
	assert(TrueDataCMSA(From) == TrueDataCMSA(To));
	assert(NumSeqsCMSA(From) == NumSeqsCMSA(To));
	for(n=0,sq=1; sq <= NumSeqsCMSA(From); sq++){
		if(!MemberSet(sq,Used)) continue;
		gsq=gsqCMSA(sq,From); sites=GetPosSitesCMSA(sq,From);
		for(blk=1; blk <= nBlksCMSA(From); blk++) RmSiteCMSA(blk,sq,sites[blk],From);
		gsq0=SwapGsqCMSA(sq,gsq,To); SwapGsqCMSA(sq,gsq0,From);
		for(blk=1; blk <= nBlksCMSA(To); blk++) AddSiteCMSA(blk,sq,sites[blk],To);
		free(sites); n++;
	} return n;
}

cma_typ	RemoveTheseColumnsCMSA(char *Operation, cma_typ CMA)
// Remove the positions within CMA with a 'd' in corresponding Operation array.
{
	assert(nBlksCMSA(CMA) == 1);
	Int4	i,j,n,len=strlen(Operation)-2;
	Int4	start,end,o,Len=LengthCMSA(1,CMA);
	assert(Len==len);
	char	c;

	cma_typ xcma,rcma=CMA;
	for(i=Len; i > 0;  i--){
		start=end=0;
		c=tolower(Operation[i]);
		if(c == 'd'){
		   for(start=end=j=i; j > 0; ){ 
			j--; c=tolower(Operation[j]);
			if(c == 'd') start=j; else break; 
		   } i=start;
		   // fprintf(stderr,"start=%d; end=%d\n",start,end);
		   if(end == Len){
                      xcma=TrimBlkCMSA(rcma,1,0,end-start+1,1); // Blk=1; RmLeft=0; RmRight < 0; limit=1.
		      if(xcma==0){ if(rcma != CMA) NilCMSA(rcma); return 0; } // too many columns removed!
                   } else if(start==1){
                      xcma=TrimBlkCMSA(rcma,1,end,0,1); // Blk=1; RmLeft < 0; RmRight=0; limit=1.
		      if(xcma==0){ if(rcma != CMA) NilCMSA(rcma); return 0; } // too many columns removed!
		   } else {
	    	      xcma=ConvertColsToInsertsCMSA(rcma,1,start,end);
		  } if(rcma != CMA) NilCMSA(rcma); rcma=xcma;
	        }
	}
	if(rcma == CMA) return 0; else return rcma;
}

#if 0	// Moved to cmsa_operations.cc
char    *AddInsertToOperationArray2(Int4 start, Int4 end, char *operation)
// ========== Add an insertion to an operational array. ===========
{

        char state,*new_operation=0;
        Int4 o,no,column;

        Int4 trace_length=strlen(operation);
        NEW(new_operation,trace_length+5,char);
        new_operation[0]='E'; no=1;
	for(no=o=1; isalpha(operation[o]); o++){
		char c = operation[o];
		if(c == 'M' || c == 'D'){ break; }
		else { new_operation[no] = c; no++; }
	}
	if(!isalpha(operation[o])){
	    fprintf(stderr,"operation[%d]='%c'=%d\n",o,operation[o],operation[o]);
	    fprintf(stderr,"operation=%s\n",operation);
	    print_error("AddInsertToOperationArray2() input error");
	}
        // Int4 InsSize=end-start+1;
        // NEW(new_operation,trace_length+InsSize+5,char);
        for(column=1,state='E'; operation[o] != 'E'; o++){
          switch(operation[o]){
            case 'M':
            case 'm':
                if(column >= start && column <=end){
                           new_operation[no]='I';
                } else new_operation[no]=operation[o];
                no++; column++; break;
            case 'D':
            case 'd': // deletion in sequence relative to profile.
                if(column >= start && column <=end){
                           // do nothing in new_operation;
                } else { new_operation[no]=operation[o]; no++; }
                column++; break;
            case 'i': // insert is between profile blocks;
                new_operation[no]=operation[o]; no++;
                break;
            case 'I': // Insert ('-') within a profile block; delete from seq.
                new_operation[no]=operation[o]; no++;
                break;
            default:
            // fprintf(stderr,"operations = %s\n",operation);
            print_error("operation( ): input error"); break;
          }  state=operation[o];
        }
        new_operation[no]='E'; no++; new_operation[no]=0;
        // printf("new operation: %s\n",new_operation);
        return new_operation;
}

Int4	IronOutOperation(char *operation)
/*************************************************************************************
 Just switch ddddIIII to IIIIdddd to fix program...so can work with cma_gblastpgp...

    ddddIIIIm
 d: 01233333
 i: 00001234

  Simply move insertions to other side of deletions; don't do this 
 *************************************************************************************/
{
	Int4 length=strlen(operation);
	// char *new_operation;
	char state,last=' ';
	Int4 i,j,strt_i=0,strt_d=0,end;
	Int4 num_i=0,num_d=0,num_fix=0;

	assert(operation[0] == 'E');
	// NEW(new_operation,length+3,char);
	if(operation[1] == 'D'){ for(i=2; i < length; i++) if(operation[i] != 'd') break; }
	else i=1;
	last=operation[i-1];
	for( ; i < length; i++){
		state = operation[i];
		if(state != last){	// then 
		   if(state=='E'){
		     // if "..dE" or "..mE" then do nothing
		     // "..IE"
		     if(num_i > 0) print_error("IronOutOperation() input error 0");
		   } else if(state=='I' && last == 'd'){	// then start collect iron out data.
			// ...dI...
			num_i=1; strt_i=i;
		   } else if(last == 'I' && num_d > 0){	// then iron out...
			num_fix++;
#if 1
			// ..dI..Im or ..dI..Id
			end = strt_d + num_i;
			for(j=strt_d; j <  end; j++) operation[j]='I';
			end += num_d;
			for( ; j <  end; j++) operation[j]='d';
		        num_i=num_d=0; strt_i=0; strt_d=0;
			if(state=='d'){
				strt_d=i; num_d=1;
			}
#else		// change inserts into matches and remove deletions...
#endif
#if 0	// dddDIII
		   } else if(state == 'D' && operation[i+1] == 'I'){
			   operation[i]='M'; i++;
			   operation[i]='d'; last='M'; state='d'; continue;
#endif
		   } else {	// state 
		     num_i=num_d=0; strt_i=strt_d=0;
		     switch(state){
			case 'I': strt_i=i; num_i=1; break;
			case 'm': case 'M': break;
			case 'D': {
			    fprintf(stderr,"i=%d\n",i);
			    fprintf(stderr,"operation = %s\n",operation);
			    print_error("IronOutOperation() input error 1"); 
			  }break;
			case 'd': strt_d=i; num_d=1; break;
		     }
		   }
		} else {
		   switch(state){
			case 'I': num_i++; break;
			case 'm': case 'M': break;
			case 'd': num_d++; break;
			default: print_error("IronOutOperation() input error 2"); break;
		   }
		} last=state;
	} return num_fix;
}
#endif

void	PutShrinkFlanksCMSA(FILE *fp,Int4 left,Int4 right,cma_typ cma)
{
	gss_typ	*gss=gssCMSA(cma);
	Int4    start,end,End,i,j,J,sq,m,s,site,blk,len,N=gss->NumSeq();
	e_type	E,subE;
	a_type	AB=AlphabetCMSA(cma);

        fprintf(fp,"[%d_(%d)=%s(%d){go=%d,gx=%d,pn=%.1f,lf=%d,rf=%d}",
		cma->Level,nBlksCMSA(cma),NameCMSA(cma),N,
		gss->GapOpen(),gss->GapExtend(),gss->PerNats(),left,right);
	if(cma->alpha > 0.0 && cma->alpha <= 1.0){
            fprintf(fp,";BPPS=(%c,%d,%d:%.4f):\n",cma->set_mode,cma->A0,cma->B0,cma->alpha);
	} else { fprintf(fp,":\n"); }

	// 2a. Print column indicators:
	for(m=1; m <= nBlksCMSA(cma); m++){
	   fprintf(fp,"(%d)",LengthCMSA(m,cma));
	   fm_type fm=ModelCMSA(m,cma);
	   for(s=1; s <= LengthCMSA(m,cma); s++){
		if(NullSiteFModel(s,fm)) fprintf(fp,".");
		else fprintf(fp,"*");
	   } fprintf(fp,"\n");
	} fprintf(fp,"\n"); fprintf(fp,"\n");
	for(sq=1; sq <= N; sq++){
		gsq_typ *gsq=gsqCMSA(sq,cma);
		E=TrueSeqCMSA(sq,cma);
		blk=1; site=1; 
		start=TruePosCMSA(sq,blk,site,cma);
		blk=nBlksCMSA(cma); site=LengthCMSA(blk,cma); 
		end=TruePosCMSA(sq,blk,site,cma);

		start=MAXIMUM(Int4,start-left,1);
		end=MINIMUM(Int4,end+right,LenSeq(E));
		subE=MkSubSeq(start,end,E);
		
		Int4 *Sites=GetPosSitesCMSA(sq,cma);
		char *operation=gsq->Operation(nBlksCMSA(cma),Sites,LengthsCMSA(cma));  //
                fprintf(stderr,"%d: operation=%s\n",sq,operation);
		char *op,c; NEW(op,strlen(operation)+3,char);
		op[0]='E';
		for(i=1,j=start,J=1; J <= LenSeq(subE); j++,i++){
			op[i]=operation[j];
			c=operation[j];
			if(strchr("MmIi",c) != NULL) J++;
		}
		while((c=operation[j]) == 'd' || c == 'D'){ op[i]=c; i++; j++; }
		op[i]='E';
                fprintf(stderr,"%d: op=%s\n",sq,op);

		for(blk=1; blk <= nBlksCMSA(cma); blk++){
			fprintf(stderr,"%d: Site[%d]=%d ... ",sq,blk,Sites[blk]);
			Sites[blk]=Sites[blk]-start+1;
			fprintf(stderr," ... %d\n",Sites[blk]);
		}

		gsq_typ *gsq0 = new gsq_typ[1];
                gsq0->initialize(op,subE,Sites);    // copies positions to sites array.

		gsq0->Put_cma_format(fp,sq,nBlksCMSA(cma),Sites,LengthsCMSA(cma),AB);
                free(operation); free(op); free(Sites);
	} fprintf(fp,"_%d].\n",cma->Level); 
}

void    PutIndelsCMSA(FILE *fp, cma_typ cma)
{
	gss_typ *gss=gssCMSA(cma);
        for(Int4 bk = 1 ; bk <= nBlksCMSA(cma); bk++){
           h_type HG=Histogram("number of hits for each repeat",0,LengthCMSA(bk,cma),2);
           Int4 totins=0,totdel=0,totilen=0;
           Int4 nins,inslen,ins,del,pos[3];
           fprintf(stdout,"Block %3d:  ins(num)  del\n",bk);
           for(Int4 col=0; col < LengthCMSA(bk,cma); col++){
              inslen=del=nins=0;
              for(Int4 sq=1; sq <= NumSeqsCMSA(cma); sq++){
                  PosSiteCMSA(bk, sq, pos, cma);
                  ins = InsertionCMSA(sq, pos[1]+col,cma);
                  if(ins > 0){ nins++; inslen+=ins; }
                  if(IsDeletedCMSA(sq,pos[1]+col, cma)){ del++; }
              } fprintf(stdout,"%3d         %3d(%3d)  %3d\n",col+1,inslen,nins,del);
              if((col+1) < LengthCMSA(bk,cma)) IncdMHist(col+1,inslen,HG);
              totins+=nins; totdel+=del; totilen+=inslen;
           } fprintf(stdout,"\n");
           fprintf(stdout,"Total:      %3d(%3d)  %3d\n",totilen,totins,totdel);
           PutHist(stdout,60,HG); NilHist(HG);
        }
}

set_typ	IndelBreakPntsCMSA(Int4 Blk, Int4 MinHalfBlk,cma_typ cma, Int4 cut)
#if 0	//*********************************************************************
	returns the set of breakpoints at which to split the single block cma into pieces.
	Algorithm:
	Starting from the position with the most insertions add break points for blocks, 
	but disallowing blocks of length <= 2 x MinHalfBlk.
#endif
{
	assert(Blk > 0 && Blk <= nBlksCMSA(cma));
	Int4	nins,inslen,ins,del,pos[3],sq,N=NumSeqsCMSA(cma);
	Int4	i,j,max,col,end,*Nino,*Nins,Len=LengthCMSA(Blk,cma);
	Int4	MinBlkLen=2*MinHalfBlk;

	// 1. Store the # of insertions and inserted residues at each position on a heap.
	dh_type dH=dheap(Len+3,4);
	NEW(Nino,Len+5,Int4); NEW(Nins,Len+5,Int4); 
	for(col=0; col < Len; col++){
	    inslen=del=nins=0;
	    for(sq=1; sq <= N; sq++){
	        PosSiteCMSA(Blk, sq, pos, cma);
	        ins = InsertionCMSA(sq, pos[1]+col,cma);
	        if(ins > 0){ nins++; inslen+=ins; }
	        // if(IsDeletedCMSA(sq,pos[1]+col, cma)){ del++; }
	    } 
	    if(col < (Len - 1)) Nino[col+1]=nins; Nins[col+1]=inslen; 
	}

	// 2. Starting from the position with the most insertions; 
	h_type HG=0; // HG=Histogram("number of residues inserted",0,Len,1);
	set_typ	used=MakeSet(Len +4); ClearSet(used);
	for(i=1; i <= MinBlkLen; i++) AddSet(i,used);
	end=(Len-MinBlkLen+1);
	for(i=end; i <= Len; i++) AddSet(i,used);
	for(max=0,col=MinBlkLen; col < end; col++){	// ignore extensions on the the end!!
	    if(Nins[col] > 0){
		insrtHeap(col,(keytyp) -Nins[col],dH);
		if(HG) IncdMHist(col,Nins[col],HG);
		if(max < Nins[col]) max = Nins[col];
	    }
	} 
	if(HG){ PutHist(stderr,60,HG); NilHist(HG); }

	// 3. Starting from the position with the most insertions; 
	if(cut==0) cut=N/100;	// 1% of sequences.
	//Int4 cut=N/50;	// 2% of sequences.
	// Int4 cut=N/10;	// 10% of sequences.
	// Int4 cut=N/5;	// 20% of sequences.
	// Int4 cut=N/4;	// 25% of sequences.
	if(cut < 1) cut = 1;
	set_typ set=MakeSet(Len+5); ClearSet(set);
	while((i=delminHeap(dH)) != NULL){
	   if(!MemberSet(i,used) && Nins[i] >= cut){
		AddSet(i,set); end=i+MinBlkLen;
		for(j=(i-MinBlkLen+1); j < end; j++) AddSet(j,used);
	   }
	} free(Nins); free(Nino); Nildheap(dH); NilSet(used);
	return set;
}

cma_typ	SplitUpBlkCMSA(Int4 Blk, set_typ Set,cma_typ cma)
// Convert block 'Blk' into multiple blocks based on sites in Set.
{
	Int4	i,j,s,b,blk,*len,nBlk=0,*pos,Len;
	assert(Blk > 0 && Blk <= nBlksCMSA(cma));
   	a_type AB=AlphabetCMSA(cma);
	gss_typ *gss=gssCMSA(cma);
	NEW(len,nBlksCMSA(cma) + CardSet(Set) +4, Int4);
	for(blk=1; blk < Blk; blk++){ len[blk]=LengthCMSA(blk,cma); }
	Len=LengthCMSA(Blk,cma);
	NEW(pos,nBlksCMSA(cma) + CardSet(Set) +4, Int4);
	for(blk=Blk-1,s=1,j=1,i=1; i <= Len; s++,i++){
		if(MemberSet(i,Set)){ blk++; pos[blk]=j; len[blk]=s; j+=s; s=0; }
	} blk++; pos[blk]=j; len[blk]=s-1; nBlk=blk; 
	for(b=Blk+1,blk++; b <= nBlksCMSA(cma); b++,blk++){ len[blk]=LengthCMSA(b,cma); }
	nBlk=blk-1;
#if 0	// DEBUG...
	PutConfigCMSA(stderr,cma);
	for(blk=1; blk <= nBlk; blk++){
		if(blk == Blk) fprintf(stderr,"*");
		fprintf(stderr,"blk%d: pos=%d; len=%d\n", blk,pos[blk],len[blk]);
	} fflush(stderr); PutSet(stderr,Set);
#endif
   	cma_typ rcma=EmptyCMSA(nBlk,len,TrueDataCMSA(cma),gss->GapOpen(),
				gss->GapExtend(),PerNatsCMSA(cma),0,0);
	assert(rcma != 0);
	Int4	sq,*sites; NEW(sites,nBlk+5,Int4);
	for(sq=1; sq <= NumSeqsCMSA(cma); sq++){
	        // ========== X. Get operational array for sequence. ===========
		gsq_typ *gsq=gsqCMSA(sq,cma);
		Int4 *Sites=GetPosSitesCMSA(sq,cma);
#if 0	// debug
		gsq->Put_cma_format(stderr,sq,nBlksCMSA(cma),Sites,LengthsCMSA(cma),AB);
#endif
#if 1
		for(blk=1; blk <= nBlksCMSA(cma); blk++) sites[blk]=Sites[blk];
		gsq_typ *gsq0=gsq->SplitBlock(Blk,pos,nBlksCMSA(cma),LengthsCMSA(cma),sites);
#if 0	// debug
		gsq0->Put_cma_format(stderr,sq,nBlksCMSA(rcma),sites,LengthsCMSA(rcma),AB);
#endif
		free(Sites);
#else	// delete this after testing above...
		char *operation=gsq->Operation(nBlksCMSA(cma),Sites,LengthsCMSA(cma));	// 
		fprintf(stderr,"%d: operation=%s\n",sq,operation);
#if 0	// debug
		gsq->Put_cma_format(stderr,sq,nBlksCMSA(cma),Sites,LengthsCMSA(cma),AB);
#endif
		free(Sites);
		gsq->SplitBlkInOperation(Blk,pos,operation);
		fprintf(stderr,"  new operation=%s\n",operation);
		gsq_typ *gsq0 = new gsq_typ[1];
        	gsq0->initialize(operation,TrueSeqCMSA(sq,rcma),sites);    // copies positions to sites array.
		free(operation);
#endif
		ReplaceCMSA(sq,gsq0,rcma); // replace sequence gsq in CMSA & fmodel with gsq0.
		for(blk=1; blk <= nBlksCMSA(rcma); blk++) AddSiteCMSA(blk,sq,sites[blk],rcma);
#if 0	// debug
		free(sites); sites=GetPosSitesCMSA(sq,rcma);
		gsq0->Put_cma_format(stderr,sq,nBlksCMSA(rcma),sites,LengthsCMSA(rcma),AB);
		PutSeq(stderr,FakeSeqCMSA(sq,rcma),AB);
		PutSeq(stderr,TrueSeqCMSA(sq,rcma),AB);
#endif
	} free(len); free(sites); 
        // sprintf(str,"%s.insrt",argv[1]); WriteMtfCMSA(str, cma, NULL);
	return rcma;
}

cma_typ	ExtendBlkCMSA(cma_typ cma, Int4 Blk, Int4 AddLeft, Int4 AddRight)
// trim the ends of block 'Blk' but not below limit.
{
   	a_type AB=AlphabetCMSA(cma);
	ExtendFakeToRealCMSA(cma);      // need the extensions to be in fake seq.
	Int4	blk,*len,numfix;
	assert(Blk > 0 && Blk <= nBlksCMSA(cma));
	assert(AddLeft >= 0 && AddRight >= 0);
	NEW(len,nBlksCMSA(cma) +3, Int4);
	for(blk=1; blk <= nBlksCMSA(cma); blk++) len[blk]=LengthCMSA(blk,cma);
	len[Blk] = len[Blk] + AddLeft + AddRight;
	gss_typ *gss=gssCMSA(cma);
   	cma_typ rcma=EmptyCMSA(nBlksCMSA(cma),len,TrueDataCMSA(cma),gss->GapOpen(),
		gss->GapExtend(),PerNatsCMSA(cma),0,0);
	for(Int4 sq=1; sq <= NumSeqsCMSA(cma); sq++){
#if 1
        	if(!SeqIsInCMSA(sq,cma)) continue;
#endif
	        // ========== X. Get operational array for sequence. ===========
		char *operation=0; // gss->Operation(sq);
		gsq_typ *gsq=gsqCMSA(sq,cma);
		if(gsq){
		  Int4	*sites=GetPosSitesCMSA(sq,cma);
		  gsq_typ *gsq0=gsq->ExtendBlk(Blk,AddLeft,AddRight,nBlksCMSA(rcma),LengthsCMSA(rcma),sites);
		  ReplaceCMSA(sq,gsq0,rcma); // replace sequence s in CMSA & fmodel.
		  for(blk=1; blk <= nBlksCMSA(rcma); blk++) AddSiteCMSA(blk,sq,sites[blk],rcma);
#if 0
		  free(sites); sites=GetPosSitesCMSA(sq,rcma);
		  gsq0->Put_cma_format(stderr,sq,nBlksCMSA(rcma),sites,LengthsCMSA(rcma),AB);
#endif
		  free(sites); 
		} else assert(gsq != 0);
	} free(len);
        // sprintf(str,"%s.insrt",argv[1]); WriteMtfCMSA(str, cma, NULL);
	return rcma;
}

cma_typ	SplitUpBlocksCMSA(set_typ Set,cma_typ cma)
// Convert a single block into split up blocks based on indels.
{
	Int4	i,j,s,blk,*len,Blk=0,*pos,Len=LengthCMSA(1,cma);
	assert(nBlksCMSA(cma) == 1);	// for now...
	// assert(Blk > 0 && Blk < nBlksCMSA(cma));
   	a_type AB=AlphabetCMSA(cma);
	gss_typ *gss=gssCMSA(cma);
	NEW(len,CardSet(Set) +4, Int4);
	NEW(pos,CardSet(Set) +4, Int4);
	for(blk=0,s=1,j=1,i=1; i <= Len; s++,i++){
		if(MemberSet(i,Set)){ blk++; pos[blk]=j; len[blk]=s; j+=s; s=0; }
	} blk++; pos[blk]=j; len[blk]=s-1; Blk=blk; 
	for(blk=1; blk <= Blk; blk++){ fprintf(stderr,"blk%d: pos=%d; len=%d\n",
		blk,pos[blk],len[blk]); } fflush(stderr);
	PutSet(stderr,Set);
   	cma_typ rcma=EmptyCMSA(Blk,len,TrueDataCMSA(cma),gss->GapOpen(),
				gss->GapExtend(),PerNatsCMSA(cma),0,0);
	Int4	*sites; NEW(sites,Blk+5,Int4);
	for(Int4 sq=1; sq <= NumSeqsCMSA(cma); sq++){
	        // ========== X. Get operational array for sequence. ===========
		gsq_typ *gsq=gsqCMSA(sq,cma);
		for(blk=1; blk <= Blk; blk++){
			sites[blk]=TruePosCMSA(sq,pos[blk],cma);
			fprintf(stderr,"%d: blk%d = %d\n",sq,blk,sites[blk]);
		} 
#if 0
		Int4 *Sites=GetPosSitesCMSA(sq,cma);
		char *operation=gsq->Operation(1,Sites,LengthsCMSA(cma));	// 
	fprintf(stderr,"%d: operation=%s\n",sq,operation);
#if 1	// debug
		gsq->Put_cma_format(stderr,sq,nBlksCMSA(cma),Sites,LengthsCMSA(cma),AB);
#endif
		free(Sites);
		gsq->SplitBlkInOperation(pos,operation);
	fprintf(stderr,"  new operation=%s\n",operation);
	gsq_typ *gsq0 = new gsq_typ[1];
	for(blk=1; blk <= Blk; blk++){ fprintf(stderr,"%d: blk%d = %d\n",sq,blk,sites[blk]); } fflush(stderr);
        gsq0->initialize(operation,TrueSeqCMSA(sq,rcma),sites);    // copies positions to sites array.
	for(blk=1; blk <= Blk; blk++){ fprintf(stderr,"%d: blk%d = %d\n",sq,blk,sites[blk]); } fflush(stderr);
		ReplaceCMSA(sq,gsq0,rcma); // replace sequence s in CMSA & fmodel.
		for(blk=1; blk <= nBlksCMSA(rcma); blk++) AddSiteCMSA(blk,sq,sites[blk],rcma);
	free(operation);
#else
		// gsq_typ *SplitSingleBlock(Int4 *pos,Int4 site, Int4 Len);
        	// void    SplitBlkInOperation(Int4 *pos,char *operation);
		// gsq_typ *gsq0=gsq->SplitBlock(Blk,pos,sites);
		// gsq_typ *gsq0=gsq->SplitBlock(Blk,pos,nBlksCMSA(cma),LengthsCMSA(cma),sites);
		gsq_typ *gsq0=gsq->SplitBlock(1,pos,nBlksCMSA(cma),LengthsCMSA(cma),sites);
		sites[Blk+1]=0; 	
		ReplaceCMSA(sq,gsq0,rcma); // replace sequence s in CMSA & fmodel.
		for(blk=1; blk <= nBlksCMSA(rcma); blk++) AddSiteCMSA(blk,sq,sites[blk],rcma);
#endif
#if 1	// debug
		free(sites); sites=GetPosSitesCMSA(sq,rcma);
		gsq0->Put_cma_format(stderr,sq,nBlksCMSA(rcma),sites,LengthsCMSA(rcma),AB);
		PutSeq(stderr,FakeSeqCMSA(sq,rcma),AB);
		PutSeq(stderr,TrueSeqCMSA(sq,rcma),AB);
#endif
	} free(len); free(sites); 
        // sprintf(str,"%s.insrt",argv[1]); WriteMtfCMSA(str, cma, NULL);
	return rcma;
}

cma_typ	TrimBlkCMSA(cma_typ cma, Int4 Blk, Int4 RmLeft, Int4 RmRight, Int4 limit)
// trim the ends of block 'Blk' but not below limit.
{
   	a_type AB=AlphabetCMSA(cma);
	Int4	blk,*len,numfix;
	assert(Blk > 0 && Blk <= nBlksCMSA(cma));
#if 0
	assert(LengthCMSA(Blk,cma) >= (RmLeft + RmRight + limit));
#else
	// Int4 x=LengthCMSA(Blk,cma) - (RmLeft + RmRight);
	if(LengthCMSA(Blk,cma) < (RmLeft + RmRight + limit)){
#if 1	// check for this in calling environment.
		return 0;
#else
		assert(!((LengthCMSA(Blk,cma) < (RmLeft + RmRight + limit)))); 
		print_error("TrimBlkCMSA(): too many columns removed from alignment");
#endif
	}
#endif
	assert(RmLeft >= 0 && RmRight >= 0);
	gss_typ *gss=gssCMSA(cma);
	NEW(len,nBlksCMSA(cma) +3, Int4);
	for(blk=1; blk <= nBlksCMSA(cma); blk++) len[blk]=LengthCMSA(blk,cma);
	len[Blk] = len[Blk]-(RmLeft + RmRight);
   	cma_typ rcma=EmptyCMSA(nBlksCMSA(cma),len,TrueDataCMSA(cma),gss->GapOpen(),
		gss->GapExtend(),PerNatsCMSA(cma),0,0);
	for(Int4 sq=1; sq <= NumSeqsCMSA(cma); sq++){
		if(!SeqIsInCMSA(sq,cma)){ continue; }
	        // ========== X. Get operational array for sequence. ===========
		char *operation=0; // gss->Operation(sq);
		gsq_typ *gsq=gsqCMSA(sq,cma);
		Int4	*sites=GetPosSitesCMSA(sq,cma);
		Int4	start=RmLeft+1,end=len[Blk]+start-1;
		gsq_typ *gsq0=gsq->TrimBlk(Blk,start,end,nBlksCMSA(cma),LengthsCMSA(cma),sites);
#if 0
		fprintf(stderr,"TrimBlkCMSA(%d,%d,%d):\n",RmLeft,RmRight,limit);
		fprintf(stderr,"%d. operation=%s\n",sq,gsq->Operation(0,AB));
		fprintf(stderr,"%d. operation=%s\n",sq,gsq0->Operation(0,AB));
#endif
		ReplaceCMSA(sq,gsq0,rcma); // replace sequence s in CMSA & fmodel.
		for(blk=1; blk <= nBlksCMSA(rcma); blk++) AddSiteCMSA(blk,sq,sites[blk],rcma);
#if 0
		free(sites); sites=GetPosSitesCMSA(sq,rcma);
		gsq0->Put_cma_format(stderr,sq,nBlksCMSA(rcma),sites,LengthsCMSA(rcma),AB);
#endif
		free(sites); 
	} free(len);
        // sprintf(str,"%s.insrt",argv[1]); WriteMtfCMSA(str, cma, NULL);
	return rcma;
}

cma_typ	CreateRandomCMA(Int4 nblk, Int4 *len,ss_type data,double per_nats)
{
        cma_typ xcma,rcma=EmptyCMSA(nblk,len,data,19,2,per_nats,5000,5000);
        Int4	i,b,ln,sq,s,n = MaxSeqSeqSet(data);
        dh_type dH=dheap(n+2,3);
        for(sq=1; sq <= NSeqsSeqSet(data); sq++) {
                    ln=SqLenSeqSet(sq,data);
                    for(b=1; b <= nblk; b++){
                        ln -= len[b]; insrtHeap(b,(keytyp)Random(),dH);
                    }
                    for(s=b; ln > 0 ; s++,ln--) insrtHeap(s,(keytyp)Random(),dH);
                    for(b=s=1; (i=delminHeap(dH)) != NULL; ){
                        // WARNING: need to generalize beyond one colinear!
                        if(i<=nblk){
                                AddSiteCMSA(b,sq,s,rcma); s+=len[b]; b++;
                        } else s++;
                    }
        } Nildheap(dH);
	xcma=OneBlockCMSA(rcma); ExtendFakeToRealCMSA(xcma); NilCMSA(rcma);
	return xcma;
}

cma_typ	RmBlockCMSA(cma_typ cma, Int4 Blk)
{
   	a_type AB=AlphabetCMSA(cma);
	Int4	b,blk,*len,numfix;
	assert(Blk > 0 && Blk <= nBlksCMSA(cma));
	gss_typ *gss=gssCMSA(cma);
	NEW(len,nBlksCMSA(cma) +3, Int4);
	for(b=0,blk=1; blk <= nBlksCMSA(cma); blk++){
	   if(blk==Blk) continue; b++; len[b]=LengthCMSA(blk,cma);
	}
   	cma_typ rcma=EmptyCMSA(nBlksCMSA(cma)-1,len,TrueDataCMSA(cma),gss->GapOpen(),
		gss->GapExtend(),PerNatsCMSA(cma),0,0);
	for(Int4 sq=1; sq <= NumSeqsCMSA(cma); sq++){
	        // ========== X. Get operational array for sequence. ===========
		char *operation=0; // gss->Operation(sq);
		gsq_typ *gsq=gsqCMSA(sq,cma);
		Int4	*sites=GetPosSitesCMSA(sq,cma);
		Int4	start,end=len[Blk]; start=end+1;
		fprintf(stderr,"Blk %d: start...end=%d..%d\n",Blk,start,end);
		gsq_typ *gsq0=gsq->RmBlock(Blk,nBlksCMSA(cma),LengthsCMSA(cma),sites);
		ReplaceCMSA(sq,gsq0,rcma); // replace sequence s in CMSA & fmodel.
		for(blk=1; blk <= nBlksCMSA(rcma); blk++) AddSiteCMSA(blk,sq,sites[blk],rcma);
#if 1
		free(sites); sites=GetPosSitesCMSA(sq,rcma);
		gsq0->Put_cma_format(stderr,sq,nBlksCMSA(rcma),sites,LengthsCMSA(rcma),AB);
#endif
		free(sites); 
	} free(len);
        // sprintf(str,"%s.insrt",argv[1]); WriteMtfCMSA(str, cma, NULL);
	return rcma;
}

cma_typ	TransferColumnCMSA(Int4 Blk,BooLean right,Int4 numCols,cma_typ cma)
{
	Int4	limit=3,AddLeft,AddRight,RmLeft,RmRight;
	cma_typ	rcma=0,tmp_cma=0;
	if(right){	// Move the last column in Blk to the right.
		assert(Blk > 0 && Blk < nBlksCMSA(cma));
		assert(numCols <= (LengthCMSA(Blk,cma) - limit));
		RmLeft=0; RmRight=numCols;
		tmp_cma=TrimBlkCMSA(cma,Blk,RmLeft,RmRight,limit);
		if(tmp_cma==0) print_error("TransferColumnCMSA(): too many aligned columns");
		AddLeft=numCols; AddRight=0;
		rcma=ExtendBlkCMSA(tmp_cma, Blk+1, AddLeft, AddRight); 
	} else {	// move the last column in Blk to the left.
		assert(Blk > 1 && Blk <= nBlksCMSA(cma));
		assert(numCols <= (LengthCMSA(Blk,cma) - limit));
		RmLeft=numCols; RmRight=0;
		tmp_cma=TrimBlkCMSA(cma,Blk,RmLeft,RmRight,limit);
		if(tmp_cma==0) print_error("TransferColumnCMSA(): too many columns removed from alignment");
		AddLeft=0; AddRight=numCols;
		rcma=ExtendBlkCMSA(tmp_cma, Blk-1, AddLeft, AddRight); 
	} NilCMSA(tmp_cma);
	return rcma;
}

cma_typ	ShiftBlkCMSA(Int4 Blk,Int4 shift,cma_typ cma)
{
	Int4	limit=3,AddLeft,AddRight,RmLeft,RmRight;
	cma_typ	rcma=0,tmp_cma=0;
	assert(shift != 0);
	if(shift < 0){ AddLeft=1; AddRight=0; RmLeft=0; RmRight=1; }
	else { AddLeft=0; AddRight=1; RmLeft=1; RmRight=0; }
	tmp_cma=ExtendBlkCMSA(cma, Blk, AddLeft, AddRight); 
	rcma=TrimBlkCMSA(tmp_cma,Blk,RmLeft,RmRight,limit);
	if(rcma==0) print_error("ShiftBlkCMSA(): too many columns removed from alignment");
	NilCMSA(tmp_cma);
	return rcma;
}

cma_typ	TrimDownCMSA(cma_typ cma, char dms_mode, Int4 limit)
{
	Int4	blk,i,j,RmLeft,RmRight;
	double **SubLLR,d,D,dd,dmax,*deleted,RE,BS,BSps;
	cma_typ tmp_cma=cma,rcma=0;

print_error("TrimDownCMSA() currently not working");

	Int4    aa_per_io=1000,aa_per_do=1000,exp_ie=3,exp_de=1;
	// Int4    aa_per_io=200,aa_per_do=200,exp_ie=2,exp_de=1;
	double	pn=PerNatsCMSA(cma);
	fprintf(stderr,"PerNatsCMSA = %g\n",pn);
	pn=1000.0;
	ssx_typ *ssx = new ssx_typ(aa_per_io,aa_per_do,exp_ie, exp_de, pn,cma,dms_mode);
#if 0	// new method
	// ssx->UsePsiScoring();
	D = ssx->GapMap();
	dd = ssx->BildLLR('M');
	NEWP(SubLLR,nBlksCMSA(tmp_cma) +3, double);
	double TotalWtSq=ssx->TotalWtSeq();
	for(dmax=-9999999,blk=1; blk <= nBlksCMSA(tmp_cma); blk++){
	   	NEW(SubLLR[blk],LengthCMSA(blk,tmp_cma) +3, double);
#if 0
	        fprintf(stderr,"Rm blk %d: %.1f nats; oldmap = %.1f; bild LLR=%.1f.\n",
			blk,D - ssx->GapMap(0.0,blk),FieldRelMapCMSA(cma,blk),dd);
#endif
	        for(i=1; i <= LengthCMSA(blk,tmp_cma); i++){
		    d=ssx->GapMap(0.0,'V',blk,i); SubLLR[blk][i]=D-d;
		    RE=ssx->RelEntropy(blk,i);
		    BS=ssx->BildScore(blk,i);
		    BSps=BS/TotalWtSq;
		    // d=0; D=SubLLR[blk][i]=ssx->BildScore(blk,i);
		    // SubLLR[blk][i]=RE;
		    // SubLLR[blk][i]=BS;
		    SubLLR[blk][i]=BSps;
		    // if((D-d) > dmax) dmax=D-d;
		    if(SubLLR[blk][i] > dmax) dmax=SubLLR[blk][i];
#if 1
	            fprintf(stderr,"Col %d.%d = %.1f nats; bild=%.1f (%.2f nats/wtseq); %.1f%c deletions; RE=%.3f.\n",
			blk,i,D-d,BS,BSps,100.0*FractDeletionsCMSA(blk,i,tmp_cma),'%',RE);
#else
	            fprintf(stderr,"Column %d.%d:  bild = %.1f nats; %.1f%c deletions .\n",
			blk,i,ssx->BildScore(blk,i),100.0*FractDeletionsCMSA(blk,i,tmp_cma),'%');
#endif
		} fprintf(stderr,"\n");
	}
#endif
	double *AveDel,**Deleted = ssx->FractionDeleted();
	NEW(AveDel,nBlksCMSA(tmp_cma) +3, double);
	i=(Int4) ceil(dmax+2);
	d=(Int4) dmax/20.0;
	// h_type HG=Histogram("Column contributions to LLR (in nats)",-i,i,d);
	for(blk=1; blk <= nBlksCMSA(tmp_cma); blk++){
	        for(D=0,i=1; i <= LengthCMSA(blk,tmp_cma); i++){
		    // IncdHist(SubLLR[blk][i],HG);
#if 0
		    fprintf(stderr,"%d.%d: %.2f nats, %.3f%c\n",
				blk,i,SubLLR[blk][i],Deleted[blk][i],'%');
#endif
		    D+=Deleted[blk][i];
		} AveDel[blk]=D/(double)LengthCMSA(blk,tmp_cma);
		// fprintf(stderr,"  ave deleted = %.3f%c.\n\n",AveDel[blk],'%');
	} // PutHist(stderr,50,HG); NilHist(HG); // exit(1);
	delete ssx;
	float info_cut=0.5;
	double *Info,InfoCut=0.12;
	for(blk=1; blk <= nBlksCMSA(tmp_cma); blk++){
	     fprintf(stderr,"Motif %3d:\n",blk);
	     fm_type fm = ModelCMSA(blk,tmp_cma);
             Int4 lenM = LenFModel(fm);
	     Info=SubLLR[blk];
	     deleted=Deleted[blk];
             float *info = InfoFModel(fm);
	     for(RmLeft=RmRight=0,i=1; i <= lenM-limit; i++){
                if(Info[i] <= InfoCut) RmLeft=i; else break;	  
                // if(Info[i] < InfoCut || deleted[i] >= 0.50) RmLeft=i; else break;	  
                // if(info[i] < info_cut || deleted[i] >= 0.50) RmLeft=i; else break;	  
             } // fprintf(stderr,"\n");
             for(j=1,i=lenM; i > (limit+RmLeft); i--,j++){
                if(Info[i] <= InfoCut) RmRight=j; else break;	  
                // if(Info[i] < InfoCut || deleted[i] >= 0.50) RmRight=j; else break;	  
                // if(info[i] < info_cut || deleted[i] >= 0.50) RmRight=j; else break;	  
	     }
	     fprintf(stderr,"RmLeft=%d..%d..RmRight=%d\n",RmLeft,lenM,RmRight);
#if 0
	     if(RmLeft == 0 && RmRight == 0 && AveDel[blk] < 0.334){ continue; }
	     if(LengthCMSA(blk,tmp_cma) < (RmLeft + RmRight + limit) || AveDel[blk] >= 0.334){
	        rcma=RmBlockCMSA(tmp_cma,blk);
		free(SubLLR[blk]); free(Deleted[blk]); 
		for(Int4 x=blk; x < nBlksCMSA(tmp_cma); x++){
			SubLLR[x] = SubLLR[x+1]; Deleted[x] = Deleted[x+1];
			AveDel[x]=AveDel[x+1];
	        }
	     } else {
		rcma=TrimBlkCMSA(tmp_cma,blk,RmLeft,RmRight,limit); 
		if(rcma==0) print_error("TrimBlkCMSA(): too many columns removed from alignment");
	     }
#else
	     if(RmLeft == 0 && RmRight == 0){ continue; }
	     if(LengthCMSA(blk,tmp_cma) < (RmLeft + RmRight + limit)){
	        rcma=RmBlockCMSA(tmp_cma,blk);
		free(SubLLR[blk]); free(Deleted[blk]); 
		for(Int4 x=blk; x < nBlksCMSA(tmp_cma); x++){
			SubLLR[x] = SubLLR[x+1]; Deleted[x] = Deleted[x+1];
			AveDel[x]=AveDel[x+1];
	        }
	     } else {
		rcma=TrimBlkCMSA(tmp_cma,blk,RmLeft,RmRight,limit); 
		if(rcma==0) print_error("TrimBlkCMSA(): too many columns removed from alignment");
	     }
#endif
	     PutConfigCMSA(stderr,rcma);
	     if(tmp_cma != cma) NilCMSA(tmp_cma); tmp_cma=rcma;
	}
#if 1	// new method
	free(AveDel);
	for(blk =1; blk <= nBlksCMSA(tmp_cma); blk++) { free(SubLLR[blk]); free(Deleted[blk]); }
#endif
	return rcma;
}

cma_typ	OneBlockCMSA(cma_typ cma)
{
   a_type AB=AlphabetCMSA(cma);
   gss_typ *gss=gssCMSA(cma);
   Int4	len[5]; len[1]=TotalLenCMSA(cma);
   cma_typ rcma=EmptyCMSA(1,len,TrueDataCMSA(cma),gss->GapOpen(),gss->GapExtend(),
			PerNatsCMSA(cma),0,0);
   Int4 i,pos[9];  pos[2]=0;
   gss_typ *rgss=gssCMSA(rcma);
   for(i=1; i<= NumSeqsCMSA(cma); i++){
	gsq_typ *gsq=gsqCMSA(i,cma);
// BooLean	found=FALSE;
	if(gsq==0){
		gsq=gsqForceCMSA(i,cma); assert(gsq != 0); 
// found=TRUE; gsq->Put(stderr,AB); // exit(1);
	} // else if(gsq->NumOpen() == 0){ gsq->Put(stderr,AB); exit(1); }

	Int4	X,*sites=GetPosSitesCMSA(i,cma);
// gsq->Put_cma_format(stderr,i,nBlksCMSA(cma),sites,LengthsCMSA(cma),AB);
        gsq_typ *gsq0 = gsq->ToOneBlk(nBlksCMSA(cma),sites,LengthsCMSA(cma),X); 
        // gsq_typ *gsq0 = gsq->ToOneBlk2(nBlksCMSA(cma),sites,LengthsCMSA(cma),X); 
        ReplaceCMSA(i,gsq0,rcma); // replace sequence s in CMSA & fmodel.
        AddSiteCMSA(1,i,X,rcma);
#if 0
if(found) {
   gsq0=gsqCMSA(i,rcma); pos[1]=X; gsq0->Put_cma_format(stderr,i,1,pos,LengthsCMSA(rcma),AB);
   WriteCMSA("junk2.cma",rcma);
   gsq0->Put(stderr,AB); exit(1);
}
#endif
	free(sites);
   } // exit(1); 
   return rcma;
}

#if 0	// moved to cma_typ (cmsa_operations.cc)
cma_typ RmWrinklesCMSA(cma_typ cma)
// remove insertions next to deletions.
{
   assert(nBlksCMSA(cma)==1);
   a_type AB=AlphabetCMSA(cma);
   gss_typ *gss=gssCMSA(cma);
   Int4 len[5]; len[1]=LengthCMSA(1,cma);
   cma_typ rcma=EmptyCMSA(1,len,TrueDataCMSA(cma),gss->GapOpen(),gss->GapExtend(),
                        PerNatsCMSA(cma),0,0);
   Int4 i,pos[9];  pos[2]=0;
   gss_typ *rgss=gssCMSA(rcma);
   for(i=1; i<= NumSeqsCMSA(cma); i++){
        gsq_typ *gsq=gsqCMSA(i,cma);
        // gsq->Put(stderr,AB);
        // gsq->Put_cma_format(stderr,i,nBlksCMSA(cma),sites,LengthsCMSA(cma),AB);
        Int4    *sites=GetPosSitesCMSA(i,cma);
        Int4    X;
        gsq_typ *gsq0 = gsq->IronOut(nBlksCMSA(cma),sites,LengthsCMSA(cma),X);
        // gsq_typ *gsq0 = gsq->ToOneBlk(nBlksCMSA(cma),sites,LengthsCMSA(cma),X);
        // gsq0->Put_cma_format(stderr,i,1,pos,LengthsCMSA(rcma),AB);
        ReplaceCMSA(i,gsq0,rcma); // replace sequence s in CMSA & fmodel.
        AddSiteCMSA(1,i,X,rcma);
        free(sites);
   } return rcma;
}
#endif

cma_typ	FuseBlocksCMSA(Int4 Blk, cma_typ cma)
// Convert a region within a block into an insert (lower case in cma file).
{
	Int4	blk,*len;
	assert(Blk > 0 && Blk < nBlksCMSA(cma));
   	a_type AB=AlphabetCMSA(cma);
	gss_typ *gss=gssCMSA(cma);
	NEW(len,nBlksCMSA(cma) +3, Int4);
	for(blk=1; blk < Blk; blk++) len[blk]=LengthCMSA(blk,cma);
	len[Blk]=LengthCMSA(Blk,cma) + LengthCMSA(Blk+1,cma);
	for(blk=Blk+2; blk <= nBlksCMSA(cma); blk++) len[blk-1]=LengthCMSA(blk,cma);
   	cma_typ rcma=EmptyCMSA(nBlksCMSA(cma)-1,len,TrueDataCMSA(cma),gss->GapOpen(),
				gss->GapExtend(),PerNatsCMSA(cma),0,0);
	for(Int4 sq=1; sq <= NumSeqsCMSA(cma); sq++){
	        // ========== X. Get operational array for sequence. ===========
		gsq_typ *gsq=gsqCMSA(sq,cma);
		Int4	*sites=GetPosSitesCMSA(sq,cma);
		gsq_typ *gsq0=gsq->FuseBlks(Blk,nBlksCMSA(cma),LengthsCMSA(cma),sites);
		sites[nBlksCMSA(cma)]=0; 	// FuseBlks() sets sites == new sites.
		ReplaceCMSA(sq,gsq0,rcma); // replace sequence s in CMSA & fmodel.
		for(blk=1; blk <= nBlksCMSA(rcma); blk++) AddSiteCMSA(blk,sq,sites[blk],rcma);
#if 0	// debug
		free(sites); sites=GetPosSitesCMSA(sq,rcma);
		gsq0->Put_cma_format(stderr,sq,nBlksCMSA(rcma),sites,LengthsCMSA(rcma),AB);
#endif
		free(sites); 
	} free(len);
        // sprintf(str,"%s.insrt",argv[1]); WriteMtfCMSA(str, cma, NULL);
	return rcma;
}

void    IronOutCMSA(cma_typ cma)
// remove 'wrinkles' in a multiblock cmsa
{
        Int4    sq,blk,X;
        // a_type AB=AlphabetCMSA(cma);
        for(sq=1; sq <= NumSeqsCMSA(cma); sq++){
                // ========== X. Get operational array for sequence. ===========
                gsq_typ *gsq=gsqCMSA(sq,cma);
                Int4    *sites=GetPosSitesCMSA(sq,cma);
                gsq_typ *gsq0=gsq->IronOut(nBlksCMSA(cma),sites,LengthsCMSA(cma));
                ReplaceCMSA(sq,gsq0,cma); // replace sequence s in CMSA & fmodel.
                for(blk=1; blk <= nBlksCMSA(cma); blk++) AddSiteCMSA(blk,sq,sites[blk],cma);
#if 0   // debug
                sites=GetPosSitesCMSA(sq,cma);
                gsq0->Put_cma_format(stderr,sq,nBlksCMSA(cma),sites,LengthsCMSA(cma),AB);
#endif
                free(sites);
        }
}



//************************* hieraln.cc ******************************

#if 1	// copied to cma_typ; still needed here too.
static Int4 RtnMisMatch(char sq1[8], char sq2[8])
// Inner loop; most time intensive...
{
	register Int4   n=0;
	if(sq1[0] != sq2[0]) n++;
	if(sq1[1] != sq2[1]) n++;
	if(sq1[2] != sq2[2]) n++;
	if(sq1[3] != sq2[3]) n++;
	if(sq1[4] != sq2[4]) n++;
	if(sq1[5] != sq2[5]) n++;
	if(sq1[6] != sq2[6]) n++;
	if(sq1[7] != sq2[7]) n++;
	return n;
}
#endif

#if 0	// moved to cma_typ 

static Int4 *RtnWorstColumnsCMSA(set_typ Set, cma_typ cma)
// The maximum = ln 20 = 2.9957 (or ln 21 = 3.044).
// find worst positions...
{
        Int4    i,j,t,k,sq,r,n;
        double  p,sum,entropy,cnt[50],tcnt;
        a_type AB = AlphabetCMSA(cma);

	dh_type dH=dheap(TotalLenCMSA(cma)+3,4);
	Int4	*list; NEW(list,TotalLenCMSA(cma) +3,Int4);

     for(sum=0.0,n=0,t=1; t <= nBlksCMSA(cma); t++){
        for(k=1; k <= LengthCMSA(t,cma); k++){
          n++;
          for(tcnt=0.0,r=0; r <= nAlpha(AB); r++) cnt[r]=0.0;
          for(sq=1; sq <= NumSeqsCMSA(cma); sq++){
	      if(Set && !MemberSet(sq,Set)) continue;
              r=ResidueCMSA(t,sq,k,cma);
              if(r){ cnt[r]+=1.0; tcnt +=1.0; }
          }
          for(entropy=0.0,r=0; r <= nAlpha(AB); r++){
                p = cnt[r]/tcnt; 
                if(p > 0.0) entropy += -p * log(p);
          } sum += entropy;
	  insrtHeap(n,(keytyp)-entropy,dH);
        }
     }
     for(i=1; !emptyHeap(dH); i++){
	// entropy=-minkeyHeap(dH);
	n=delminHeap(dH);  list[i]=n; 
	// fprintf(stderr,"%d.%d %.3f\n",i,n,entropy);
     } Nildheap(dH); return list;
}

static UInt8 **RtnBlkAln(Int4 &NumBlks,Int4 &MaxMisMatch, Int4 percent_ident,set_typ InSet,cma_typ cma)
{
	Int4	b,i,j,k,N=NumSeqsCMSA(cma),score,nblk=nBlksCMSA(cma);
	Int4	s,sb,bl,si,sj,i8,pos[3],seti,setj;
	unsigned char	r;
        a_type  AB=AlphabetCMSA(cma);
	Int4	start=1;
	FILE	*efp=0; // efp=stderr;

	if(percent_ident < 0){ 	// then keep first sequence in cma output file...
		start=2; percent_ident = -percent_ident;
	} assert(percent_ident > 0 && percent_ident <= 100);

	// [0,1,2,3, 4,5,6,7] = 8 residues at a time.
	assert(sizeof(UInt4) != sizeof(UInt8)); // make sure this is a 64 bit machine...
	assert(sizeof(UInt8)==8); // make sure this is a 64 bit machine...
	char	*bsq; 
	Int4	col,TotLen=TotalLenCMSA(cma);
	double	d=(double) (TotLen*((double)(100-percent_ident)/100.0));
	double	D=(double) TotLen/8.0;
	MaxMisMatch= (Int4) floor(d); NumBlks=(Int4) ceil(D); 

	// 0. Store each fake sequence in 8 byte arrays.
	UInt8	**SQ; NEWP(SQ,N+3,UInt8);
#if 0	// Original code.
	for(i = start; i <= N; i++){
  	  if(!MemberSet(i,InSet)) continue;
	  NEW(SQ[i],(TotLen/7)+2,UInt8);
#if 0	  // purify uninitialized read error due to short seq?
	  unsigned char *isq = SeqPtrCMSA(i,cma);
	  for(bl=sb=0,b=1; b <= nblk; b++){
	    PosSiteCMSA(b,i,pos,cma); si=pos[1];
	    for(s=1; s <= LengthCMSA(b,cma); s++,si++,sj++){
		if(sb%8==0){ bsq=(char *) &SQ[i][bl]; bl++; sb=0; } 
		bsq[sb] = AlphaChar(isq[si],AB); sb++;
	    }
	  } // fprintf(stderr,"%d: %s\n",i,SQ[i]);
#else
	  for(bl=sb=0,b=1; b <= nblk; b++){
	    for(s=1; s <= LengthCMSA(b,cma); s++){
		if(sb%8==0){ bsq=(char *) &SQ[i][bl]; bl++; sb=0; } 
		r=ResidueCMSA(b,i,s,cma);
		bsq[sb] = AlphaChar(r,AB); sb++;
	    }
	  } // fprintf(stderr,"%d: %s\n",i,SQ[i]);
#endif
	}
#else	//====== NEW faster method.. ===========
	// sort the coluumns from least conserved first to hit mismatches sooner!
	assert(nblk == 1);
	Int4 *list=RtnWorstColumnsCMSA(InSet,cma);
	for(i=start; i <= N; i++){
  	  if(!MemberSet(i,InSet)) continue;
	  NEW(SQ[i],(TotLen/7)+2,UInt8);
	  for(bl=sb=0,j=1; j <= LengthCMSA(1,cma); j++){
		s=list[j];
		if(sb%8==0){ bsq=(char *) &SQ[i][bl]; bl++; sb=0; } 
		r=ResidueCMSA(1,i,s,cma); bsq[sb] = AlphaChar(r,AB); sb++;
	  } // fprintf(stderr,"%d: %s\n",i,SQ[i]);
	} free(list);
#endif
	if(efp) fprintf(efp,"NumBlks = %d; TotLen= %d; MaxMisMatch=%d\n",NumBlks,TotLen,MaxMisMatch);
	return SQ;
}

ds_type GetRepSetsStatic(Int4 start,Int4 N,Int4 NumBlks,Int4 MaxMisMatch,
				set_typ InSet, Int4 *edges,UInt8 **SQ)
// Don't make this static as the optimization option eliminates this and makes it slower!!!
{
	ds_type sets = DSets(N);
	register Int4	i,j,seti,setj,score,bl,sb;
	for(i = start; i < N; i++){
  	  if(!MemberSet(i,InSet)) continue;
	  seti=0;
	  // if(fp && i % 1000 == 0) fprintf(fp,"\r%.1f",100.0*((double)i/(double)N));
	  for(j=i+1; j <= N; j++){
  	     if(!MemberSet(j,InSet)) continue;
	     if(seti == 0) seti=findDSets(i,sets);
	     setj=findDSets(j,sets);
	     if(seti != setj){
	      UInt8 *SQ_i=SQ[i],*SQ_j=SQ[j];
	      for(score=0,bl=0,sb=0; bl <= NumBlks; bl++){
		if(SQ_i[bl] != SQ_j[bl]){
		   score += RtnMisMatch((char *)&SQ_i[bl],(char *)&SQ_j[bl]);
		   if(score > MaxMisMatch) break;
		}
	      }
	      if(score <= MaxMisMatch) { 
		   edges[i]++; edges[j]++;
		   seti=linkDSets(seti,setj,sets);
	      }
	     }
	  }
	} return sets;
}

set_typ	RtnFastRepSetCMSA(FILE *fp, Int4 percent_ident,set_typ InSet,cma_typ cma)
// return a representative set of sequences from cma.
{
	ds_type sets;
	Int4	b,i,j,k,N=NumSeqsCMSA(cma),score,nblk=nBlksCMSA(cma);
	Int4	s,sb,bl,si,sj,i8,pos[3],seti,setj;
        a_type  AB=AlphabetCMSA(cma);
	unsigned char	r;
	Int4	start=1;

	Int4	NumBlks=0,MaxMisMatch=0; 
	UInt8	**SQ=RtnBlkAln(NumBlks,MaxMisMatch,percent_ident,InSet,cma);
	if(fp) fprintf(fp,"NumBlks=%d; MaxMisMatch=%d\n",NumBlks,MaxMisMatch);
	if(percent_ident < 0) start=2; 
	assert(percent_ident > 0 && percent_ident <= 100);

	//============ 1. Cluster sequences into sets at the percent identity cutoff. ============
	
	Int4	*edges; NEW(edges,N+3,Int4);	// number of edges out of each node.
	h_type	HG=0,HG2=0,HG3=0;
	if(fp) HG=Histogram("number of edges",0,25,1.0);

#if 1	// new routine.
	sets=GetRepSetsStatic(start,N,NumBlks,MaxMisMatch,InSet,edges,SQ);
#else
	sets = DSets(N);
	for(i = start; i < N; i++){
  	  if(!MemberSet(i,InSet)) continue;
	  seti=0;
	  // if(fp && i % 1000 == 0) fprintf(fp,"\r%.1f",100.0*((double)i/(double)N));
	  for(j=i+1; j <= N; j++){
  	     if(!MemberSet(j,InSet)) continue;
	     if(seti == 0) seti=findDSets(i,sets);
	     setj=findDSets(j,sets);
	     if(seti != setj){
	      UInt8 *SQ_i=SQ[i],*SQ_j=SQ[j];
	      for(score=0,bl=0,sb=0; bl <= NumBlks; bl++){
		if(SQ_i[bl] != SQ_j[bl]){
		   score += RtnMisMatch((char *)&SQ_i[bl],(char *)&SQ_j[bl]);
		   if(score > MaxMisMatch) break;
		}
	      }
	      if(score <= MaxMisMatch) { 
		   edges[i]++; edges[j]++;
		   seti=linkDSets(seti,setj,sets);
	      }
	     }
	  }
	}
#endif
	if(fp){
	  for(i=1; i <= N; i++){
	    IncdHist((double)edges[i],HG);
	    // if(edges[i] > 70) PutSeq(stderr,TrueSeqCMSA(i,cma),AB);
	  } PutHist(fp,50,HG); NilHist(HG);
	}
	for(i = start; i <= N; i++) if(SQ[i]) free(SQ[i]);  free(SQ);

	// 2. Within each set pick the sequence with the highest score.
	double	bestprob,jprob,var;
	Int4	best;
	set_typ Set=0,SubSet=MakeSet(SetN(InSet)); ClearSet(SubSet);
	if(fp){
	  HG=Histogram("number of edge difference",0,25,1.0);
	  // Int4 inc=(Int4) ceil((double) CardSet(InSet)/100.0); 
	  HG2=Histogram("seqs set sizes",0,30,1);
	  HG3=Histogram("number of seqs in each set",0,100,1);
	  Set=MakeSet(N+3); // ClearSet(Set);
	}
        for(s=0,i=1; i <= N; i++){
  	  if(!MemberSet(i,InSet)) continue;
	  seti=findDSets(i ,sets);
	  if(i != seti) continue;	// skip non-canonical sequences
	  best=i;
	  // bestprob=GetTotalProbCMSA(i,cma);
	  bestprob=(double)edges[i];
	  if(fp){ ClearSet(Set); AddSet(i,Set); }
	  for(j=1; j <= N; j++){
  	        if(!MemberSet(j,InSet)) continue;
		if(j!=i){
		  setj=findDSets(j,sets);
		  if(seti == setj){
		        if(fp) AddSet(j,Set);
			// jprob=GetTotalProbCMSA(j,cma);
			jprob=(double) edges[j];
			if(jprob > bestprob){ best = j; bestprob=jprob; }
		  }
		}
	  } AddSet(best,SubSet); 
	  if(fp){
	     Int4 card=CardSet(Set);
	     if(card > 2){ s++; IncdMHist((double)s,card, HG3); }
	     IncdHist(card,HG2);
	     for(j=1; j <= N; j++){
		if(!MemberSet(j,Set)) continue;
		for(Int4 jj=1; jj <= N; jj++){
		  if(!MemberSet(jj,Set) || j == jj) continue;
		  IncdHist(abs(edges[j]-edges[jj]),HG);
  		}
	     }
	  }
	}
	if(fp){
	   PutHist(fp,50,HG); NilHist(HG); PutHist(fp,50,HG2); NilHist(HG2);
	   PutHist(fp,50,HG3); NilHist(HG3); NilSet(Set);
	} NilDSets(sets); free(edges);
	return SubSet;
}
#endif

set_typ	RtnRepSubSetCMSA(FILE *fp_err, Int4 percent_ident,set_typ Set,cma_typ cma)
// return a representative set of sequences from cma.
{
	ds_type sets;
	Int4	b,i,j,k,N=NumSeqsCMSA(cma);
	Int4	score,nblk=nBlksCMSA(cma);
	Int4	s,si,sj,pos[3],seti,setj;
	ss_type	data=DataCMSA(cma);
        a_type  A=AlphabetCMSA(cma);
	unsigned char	*isq,*jsq;
	Int4	start=1;
	assert(SetN(Set) >= N);

	if(percent_ident < 0){ 	// then keep first sequence in cma output file...
		start=2; percent_ident = -percent_ident;
	}
	assert(percent_ident > 0 && percent_ident <= 100);
	// 1. Cluster sequences into sets at the percent identity cutoff.
	Int4 total = TotalLenCMSA(cma);
	sets = DSets(N);
	for(i = start; i < N; i++){
	  if(!MemberSet(i,Set)) continue;
	  isq = SeqPtrCMSA(i,cma);
	  seti=findDSets(i,sets);
	  if(fp_err && i % 1000 == 0) fprintf(stderr,"\r%.1f",100.0*((double)i/(double)N));
	  for(j=i+1; j <= N; j++){
	     if(!MemberSet(j,Set)) continue;
	     setj=findDSets(j,sets);
	     if(seti != setj){
	      jsq = SeqPtrCMSA(j,cma);
	      for(score=0,b=1; b <= nblk; b++){
		PosSiteCMSA(b,i,pos,cma); si=pos[1];
		PosSiteCMSA(b,j,pos,cma); sj=pos[1];
		for(s=1; s <= LengthCMSA(b,cma); s++,si++,sj++){
			if(isq[si] == jsq[sj]) score++;
		}
	      } // score = (Int4) floor(((double)score*100.0/(double)total) +0.5);
	      score = (Int4) floor(((double)score*100.0/(double)total));
	      if(score >= percent_ident) seti=linkDSets(seti,setj,sets);
	     }
	  }
	}
	// 2. Within each set pick the sequence with the highest profile score.
	double	bestprob,jprob;
	Int4	best;
	set_typ	SubSet=MakeSet(SetN(Set)); ClearSet(SubSet);
        for(s=0,i=1; i <= N; i++){
	  if(!MemberSet(i,Set)) continue;
	  seti=findDSets(i,sets);
	  if(i==seti){	// is this the canonical sequence?
	     best=i;
	     bestprob=GetTotalProbCMSA(i,cma);
	     for(j=1; j <= N; j++){
	        if(!MemberSet(j,Set)) continue;
		if(j!=i){
		  setj=findDSets(j,sets);
		  if(seti == setj){
			jprob=GetTotalProbCMSA(j,cma);
			if(jprob > bestprob){ best=j; bestprob=jprob; }
		  }
		}
	     } AddSet(best,SubSet); s++;
	  }
	} NilDSets(sets);
	return SubSet;
}

cma_typ	InSetMkCMSA(set_typ set, cma_typ cma)
{
	a_type  AB=AlphabetCMSA(cma);
	FILE *fp=tmpfile(); PutInSetCMSA(fp,set,cma); rewind(fp);
	cma_typ rcma = ReadCMSA(fp,AB); fclose(fp);
	return rcma;
}

void	ExtendAlnToEndsCMSA(cma_typ cma)
{
}

cma_typ	LengthenBlksCMSA(cma_typ cma)
{
	Int4	blk,R,AddLeft=0,AddRight=0,min_half_blk=6,site;
	Int4    aa_per_io=500,aa_per_do=500,exp_ie=3,exp_de=1;
	char	dms_mode='S',C;
	ssx_typ *ssx=0,*rssx=0;
	double	DD,D,d,dd,pn=1000.0;
	cma_typ xcma=0,tmpcma=0,rcma=0;
	assert(nBlksCMSA(cma) == 1);

	//=========== 1. 
	set_typ bset=IndelBreakPntsCMSA(1,min_half_blk,cma);
	PutConfigCMSA(stderr,cma);
	if(CardSet(bset) > 0) xcma=SplitUpBlkCMSA(1,bset,cma); NilSet(bset);
	if(xcma == 0) return 0;
	PutConfigCMSA(stderr,xcma);
	dh_type dH=dheap(nBlksCMSA(xcma) *3,4);
	for(blk=1; blk <= nBlksCMSA(xcma); blk++){
		R=blk*2; insrtHeap(R,(keytyp) Random(),dH);
		R=blk*2+1; insrtHeap(R,(keytyp) Random(),dH);
	}
	char	method='B';	method=' ';	// virtual map...
	while(!emptyHeap(dH)){
	   R=delminHeap(dH); AddLeft=0; AddRight=0;
	   if(R%2 == 0){ blk = R/2; AddLeft=1; C='L'; site = 1; }
	   else { blk = (R-1)/2; AddRight=1; C='R'; site=LengthCMSA(blk,xcma) + 1; }
	   if(blk == 1 && AddLeft > 0) continue;
	   // if(blk == nBlksCMSA(xcma) && AddRight > 0) continue;
	   tmpcma=ExtendBlkCMSA(xcma,blk, AddLeft, AddRight);
	   if(tmpcma){
		fprintf(stderr,"Try adding column to block %d%c\n",blk,C);
   		ssx = new ssx_typ(aa_per_io,aa_per_do,exp_ie, exp_de, pn,xcma,dms_mode); 
		DD = ssx->GapMap(0,method); d = DD/(double) TotalLenCMSA(xcma);
		fprintf(stderr,"Original BILD LLR=%.3f (%.3f nats per column).\n",DD,d);
		delete ssx;
  		rssx = new ssx_typ(aa_per_io,aa_per_do,exp_ie, exp_de, pn,tmpcma,dms_mode); 
  		D = rssx->GapMap(0,method); d = D/(double) TotalLenCMSA(tmpcma);
double dd_cut=0.7;
		dd=rssx->FractionDeleted(blk,site);
		fprintf(stderr," fraction deleted in this position = %.2f\n",dd);
		fprintf(stderr," Modified BILD LLR=%.3f (%.3f nats per column).\n",D,d);
		delete rssx; // dd=0;
		if(D > DD && dd < dd_cut){
			fprintf(stderr,"Column added to block %d%c\n",blk,C);
			insrtHeap(R,(keytyp) Random(),dH);
			if(xcma!= cma) NilCMSA(xcma); xcma=tmpcma;
		} else { NilCMSA(tmpcma); tmpcma=0; }
	   }
	}
	if(xcma != cma){ rcma=OneBlockCMSA(xcma); NilCMSA(xcma); } else rcma=0;
	return rcma;
} 

cma_typ	gmb_typ::Optimize()
{
	cma_typ rcma=0;
#if 0
	double	map,bestmap,temp=300.0;
	do { 
		rcma=Sample("junk",'S',50,temp);
	} while(map > best_map);
#endif
	return rcma;
}

Int4    ModeLengthSeqSet(ss_type data)
// return the mode of the sequence length distribution.
{
        Int4    sq,i,M,N=NSeqsSeqSet(data);
        dh_type dH=dheap(N+2,4);
        for(sq=1; sq <= N; sq++) insrtHeap(sq,-(keytyp)SqLenSeqSet(sq,data),dH);
        M=MAXIMUM(Int4,1,N/2);
        for(i=1; sq=delminHeap(dH); i++){ if(i >= M) break; }
        Nildheap(dH); return SqLenSeqSet(sq,data);
}

#if 0	// moved to cmsa.h
Int4	NumInsertsCMSA(Int4 blk, Int4 site, cma_typ cma, Int4 &ins_res)
// return number of insertions following column 'site'.
{
	Int4	sq,n,nins=0; 
	for(ins_res=0,sq=1; sq <= NumSeqsCMSA(cma); sq++){
#if 1
		if(!SeqIsInCMSA(sq,cma)) continue;
#else
		assert(SeqIsInCMSA(sq,cma));
#endif
		n = InsertionCMSA(blk,sq,site,cma);
		if(n > 0){ nins++; ins_res += n; }
	} return nins;
}
#endif

Int4	*CommonColsCMSA(cma_typ cmaA, cma_typ cmaB, double *&DD,double *&Del,double frq_cut)
// return the set of columns in cmaA that also occur in cmaB.
{
	Int4	i,j,k,sq,N=NumSeqsCMSA(cmaA),lenA=LengthCMSA(1,cmaA),lenB=LengthCMSA(1,cmaB);
	//===================== 1. Check for input errors ==========================
	assert(nBlksCMSA(cmaA) == 1); assert(nBlksCMSA(cmaB) == 1);
	if(NumSeqsCMSA(cmaA) !=  NumSeqsCMSA(cmaB)){
		fprintf(stderr,"NumSeqs: %d != %d\n",NumSeqsCMSA(cmaA),NumSeqsCMSA(cmaB));
		print_error("FindCommonColumnsCMSA(): input alignments are inconsistent");
	} 
	for(sq=1; sq <= N; sq++){
	    e_type sqA=TrueSeqCMSA(sq,cmaA),sqB=TrueSeqCMSA(sq,cmaB);
	    if(!FastIdentSeqs(sqA,sqB)){
		a_type	AB=AlphabetCMSA(cmaA);
		PutSeq(stderr,sqA,AB); PutSeq(stderr,sqB,AB);
		AlnSeqSW(stderr,11,1,sqA,sqB,AB);
		fprintf(stderr,"Sequences in row %d differ.\n",sq);
		print_error("FindCommonColumnsCMSA(): input sequences are inconsistent");
	    }
	} assert(frq_cut >= 0.20 && frq_cut <= 0.50);

	//===================== 2. Allocate memory for analysis. ==========================
	Int4	n,x,y,z,M=MaxTrueSeqCMSA(cmaA),max_j;
	double	d,dd,max_d;
	set_typ SetB=MakeSet(lenB+3),*setB; NEW(setB,lenB+4,set_typ); 
	for(j=1; j <=lenB; j++){ setB[j]=MakeSet(N+4); }
	Int4	*ColI2J; NEW(ColI2J,lenA+4,Int4); 	// mapping of columns in cmaB to cmaA.
	if(lenA > lenB){ NEW(DD,lenA+4,double); NEW(Del,lenA+4,double); }
	else { NEW(DD,lenB+4,double); NEW(Del,lenB+4,double); }
	set_typ SetJ=MakeSet(lenB+3); ClearSet(SetJ);

	//================== 3. Compare seq alignments for A & B. ========================
	for(i=1; i <=lenA; i++){
	    ClearSet(SetB);
	    for(j=1; j <=lenB; j++){ ClearSet(setB[j]); }
	    for(sq=1; sq <= N; sq++){
		if(IsDeletedCMSA(1,sq,i,cmaA)) continue; 	// columns don't correspond...
	        x=TruePosCMSA(sq,1,i,cmaA);
	    	for(j=1; j <= lenB; j++){
		  if(IsDeletedCMSA(1,sq,j,cmaB)) continue;
	          y=TruePosCMSA(sq,1,j,cmaB);	// a match position....
		  if(y == x){ AddSet(sq,setB[j]); } // found a common alignment.
		  if(y >= x) break;	// already past position x in the sequence...
		}
	    }
	    for(max_j=0,max_d=dd=0.0,j=1; j <= lenB; j++){
		x=CardSet(setB[j]);
		d=(double)x/(double)N; dd += d;
		if(d >= 0.50){ //  then map j --> i;
			AddSet(j,SetJ); AddSet(j,SetB); ColI2J[i]=j; if(DD) DD[i]=d; break; 
		} else if(d > max_d){ max_d=d; max_j=j; }
		if(dd > 0.66) break;
	    }
#if 1
	    if(ColI2J[i]==0 && !MemberSet(max_j,SetB) && !MemberSet(max_j,SetJ) && max_d > frq_cut){
		AddSet(max_j,SetJ); AddSet(max_j,SetB); ColI2J[i]=max_j; DD[i]=max_d; 
	    }
#endif
	}
	
	//================== 4. Free sets. ========================
	NilSet(SetB); NilSet(SetJ);
	for(j=1; j <= lenB; j++){ NilSet(setB[j]);} free(setB);

	//================== 5. Compute # deletions at each position. ========================
        for(sq=1; sq <= N; sq++){
	    for(i=1; i <= lenA; i++){
		j=ColI2J[i]; 
		if(j==0) continue; assert(j > 0 && j <= lenB);
                if(!SeqIsInCMSA(sq,cmaA) || !SeqIsInCMSA(sq,cmaB)) continue;
                // if(IsDeletedCMSA(1,sq,i,cmaA) || IsDeletedCMSA(1,sq,j,cmaB)) Del[i]++;
                if(IsDeletedCMSA(1,sq,i,cmaA) && IsDeletedCMSA(1,sq,j,cmaB)) Del[i]++;
                // if(IsDeletedCMSA(1,sq,i,cmaA)) Del[i]++;
            }
	} return ColI2J;
}

cma_typ	SortByCMA(cma_typ cma, cma_typ scma)
// Sort the sequences in the scma file to correspond to order of cma file.
// Print out an error message if Template contains 
{
	ss_type data = TrueDataCMSA(cma);  // == True sequences...
        Int4    *list,i,j,I,J,Score,Number,N,n;
	a_type	AB=AlphabetCMSA(cma);
	BooLean	*Found;
	e_type	sqJ,sqI;
	char	strI[100],strJ[100];
#if 0
	Number=N=NumSeqsCMSA(cma);
	assert(N == NumSeqsCMSA(scma));
	NEW(list, N+3,Int4); NEW(Found, Number+4,BooLean);
	for(n=0,J=1; J <= N; J++){
	   sqJ = TrueSeqCMSA(J,cma);	// get template sequence
	   for(I = 1; I <= Number; I++){
		if(Found[I]) continue;
	   	sqI = TrueSeqCMSA(I,scma);	// get template sequence
		if(FastIdentSeqs(sqJ,sqI)){
		   assert(list[J] == 0);
		   list[J]=I; n++; Found[I]=TRUE; break;
		}
	   }
	   if(I > Number){	// corresponding sequence not found...
		PutSeq(stderr,sqJ,AB);
		fprintf(stderr,"I = %d; Number=%d\n",I,Number);
		print_error("SortByCMA() error 2: missing sequence");
	   }
	} 
	for(J = 1; J <= N; J++){
		if(!Found[J]){
		    fprintf(stderr,"missing seq. %d\n",J);
		    PutSeq(stderr,TrueSeqCMSA(J,cma),AB);
		    // print_error("SortByCMA() error 3: missing sequence");
		}
	}
#else
	Number=NumSeqsCMSA(cma); N=NumSeqsCMSA(scma);
	assert(N <= Number);
	ss_type dataI=DataCMSA(scma),dataJ=DataCMSA(cma);
	NEW(list, N+3,Int4); NEW(Found, Number+4,BooLean);
	for(n=0,J=1; J <= Number; J++){
	   sqJ = TrueSeqCMSA(J,cma);	// get template sequence
	   StrSeqID(strJ,40,SeqSetE(J,dataJ));
	   for(I = 1; I <= N; I++){
		if(Found[I]) continue;
	        StrSeqID(strI,40,SeqSetE(I,dataI));
		if(strncmp(strI,strJ,40) != 0) continue;
	   	sqI = TrueSeqCMSA(I,scma);	// get template sequence
		if(FastIdentSeqs(sqJ,sqI)){ n++; list[n]=I; Found[I]=TRUE; break; }
	   }
	} 
	for(I=1; I <= N; I++){
		if(!Found[I]){
		    fprintf(stderr,"missing seq. %d\n",I);
		    PutSeq(stderr,TrueSeqCMSA(I,scma),AB);
		    print_error("SortByCMA() error 3: missing sequence");
		}
	}
#endif
	FILE *fp=tmpfile();
	PutSelectOneCMSA(fp,0,list,scma); free(list);  free(Found);
	rewind(fp); cma_typ rcma=ReadCMSA(fp,AB); fclose(fp);
	return rcma;
}

#define TestETvsBPPS 1

void	ColumnQualityCMSA(ssx_typ *ssx)
{
	cma_typ	cma=ssx->RtnCMA();
	Int4	blk,i,j,RmLeft,RmRight;
	double **SubLLR,d,D,dd,dmax,*deleted,RE,BS,BSps,TotalBILD=0.0;
	double	pn=PerNatsCMSA(cma); fprintf(stderr,"PerNatsCMSA = %g\n",pn);
	pn=1000.0;
	// ssx->UsePsiScoring();
	D = ssx->GapMap();
	dd = ssx->BildLLR('M');
	NEWP(SubLLR,nBlksCMSA(cma) +3, double);
	h_type HG=Histogram("Column bild scores (in nats)",0,2000,20);
	double TotalWtSq=ssx->TotalWtSeq();
#if TestETvsBPPS	// Testing ET versus BPPS.
	dh_type dH=dheap(TotalLenCMSA(cma)+5002,4);
#endif
	for(dmax=-9999999,blk=1; blk <= nBlksCMSA(cma); blk++){
	   	NEW(SubLLR[blk],LengthCMSA(blk,cma) +3, double);
#if 0
	        fprintf(stderr,"Rm blk %d: %.1f nats; oldmap = %.1f; bild LLR=%.1f.\n",
			blk,D - ssx->GapMap(0.0,'V',blk),FieldRelMapCMSA(cma,blk),dd);
#endif
	        for(i=1; i <= LengthCMSA(blk,cma); i++){
		    // d=ssx->GapMap(0.0,'V',blk,i); SubLLR[blk][i]=D-d;
		    RE=ssx->RelEntropy(blk,i);
		    BS=ssx->BildScore(blk,i); TotalBILD += BS;
IncdHist(BS,HG);
		    BSps=BS/TotalWtSq;
		    Int4 ins_res,n_ins=NumInsertsCMSA(blk, i,cma,ins_res);
		    Int4 n_del=NumDeletionsCMSA(blk, i,cma);
		    // d=0; D=SubLLR[blk][i]=ssx->BildScore(blk,i);
		    // SubLLR[blk][i]=RE;
		    // SubLLR[blk][i]=BS;
		    SubLLR[blk][i]=BSps;
		    // if((D-d) > dmax) dmax=D-d;
		    // if(SubLLR[blk][i] > dmax) dmax=SubLLR[blk][i];
#if TestETvsBPPS	// Testing ET versus BPPS.
		    Int4 sq=1,x=TruePosCMSA(sq,blk,i,cma);
		    x = x + OffSetSeq(SeqSetE(sq,TrueDataCMSA(cma)));
#if 1
		    double Ent=ssx->Entropy(blk,i);
		    if(x > 0) insrtHeap(x,(keytyp)-Ent,dH);
#elif 0
		    if(x > 0) insrtHeap(x,(keytyp)-BS,dH);
#else
		    if(x > 0) insrtHeap(x,(keytyp)-RE,dH);
#endif
#endif
#if 1
	            fprintf(stderr,"%d.%d: bild=%.1f nats (%.2f npws); RE=%.3f; ins=%d (%d) del=%d (%.1f%c).\n",
			blk,i,BS,BSps,RE,n_ins,ins_res,n_del,100.0*FractDeletionsCMSA(blk,i,cma),'%');
		    if(ins_res > 0) ssx->PutPenalties(stderr,blk,i);
		    else if(n_del > 10) ssx->PutPenalties(stderr,blk,i);
		    else if(i == LengthCMSA(blk,cma)) ssx->PutPenalties(stderr,blk,i);
#else
	            fprintf(stderr,"Column %d.%d:  bild = %.1f nats; %.1f%c deletions .\n",
			blk,i,ssx->BildScore(blk,i),100.0*FractDeletionsCMSA(blk,i,cma),'%');
#endif
		} fprintf(stderr,"\n");
	} fprintf(stderr,"Total BILD score = %.2f (%.2f nats/column)\n",
		TotalBILD,TotalBILD/(double)TotalLenCMSA(cma)); 
#if TestETvsBPPS	// Testing ET versus BPPS.
	fprintf(stdout,"B=");
	for(j=1; !emptyHeap(dH) && j <= 50; j++){
		double score=-(double)minkeyHeap(dH);
	        assert((i=delminHeap(dH)) != 0);
		// fprintf(stderr,"%d.%d %.2f\n",j,i,score);
		fprintf(stdout,"%d,",i);
        } Nildheap(dH); fprintf(stdout,"\n\n");
	fprintf(stderr,"To see more, remove exit from cma_gmb.cc\n");
	exit(1);
#endif
	PutHist(stderr,50,HG); NilHist(HG); // exit(1);
	double *AveDel,**Deleted = ssx->FractionDeleted();
	NEW(AveDel,nBlksCMSA(cma) +3, double);
	// i=(Int4) ceil(dmax+2); d=dmax/20.0;
	// HG=Histogram("Column contributions to LLR (in nats)",-i,i,d);
	for(blk=1; blk <= nBlksCMSA(cma); blk++){
	        for(D=0,i=1; i <= LengthCMSA(blk,cma); i++){
		    // IncdHist(SubLLR[blk][i],HG);
#if 0
		    fprintf(stderr,"%d.%d: %.2f nats, %.3f%c\n",
				blk,i,SubLLR[blk][i],Deleted[blk][i],'%');
#endif
		    D+=Deleted[blk][i];
		} AveDel[blk]=D/(double)LengthCMSA(blk,cma);
		// fprintf(stderr,"  ave deleted = %.3f%c.\n\n",AveDel[blk],'%');
	} // PutHist(stderr,50,HG); NilHist(HG); // exit(1);
	float info_cut=0.5;
	double *Info,InfoCut=0.12;
	free(AveDel);
	for(blk =1; blk <= nBlksCMSA(cma); blk++) { free(SubLLR[blk]); free(Deleted[blk]); }
	free(SubLLR); free(Deleted); 
}

double	**BILD_ScoresCMSA(ssx_typ *ssx)
{
	cma_typ	cma=ssx->RtnCMA();
	Int4	blk,i,j;
	double **RtnBS,d,dmax,BS,BSps,TotalBILD=0.0;
	double TotalWtSq=ssx->TotalWtSeq();
	FILE	*efp=0; // efp=stderr;
	h_type HG=0;

	if(efp) HG=Histogram("Column bild scores (in nats)",0,2000,20);
	NEWP(RtnBS,nBlksCMSA(cma) +3, double);
	for(dmax=-9999999,blk=1; blk <= nBlksCMSA(cma); blk++){
	   	NEW(RtnBS[blk],LengthCMSA(blk,cma) +3, double);
	        for(i=1; i <= LengthCMSA(blk,cma); i++){
#if 0	// relative entropy
		    BS=ssx->RelEntropy(blk,i);
#else	// BILD scores.
		    BS=ssx->BildScore(blk,i); 
#endif
		    RtnBS[blk][i]=BS; TotalBILD += BS;
		    if(efp) IncdHist(BS,HG);
		    BSps=BS/TotalWtSq;
		    RtnBS[blk][i]=BSps;
		    if(RtnBS[blk][i] > dmax) dmax=RtnBS[blk][i];
	            if(efp) fprintf(stderr,"%d.%d: BILD=%.1f (%.2f npws).\n", blk,i,BS,BSps);
		} if(efp) fprintf(stderr,"\n");
	}
	if(efp){
	    fprintf(stderr,"Total BILD score = %.2f (%.2f nats/column)\n",
		TotalBILD,TotalBILD/(double)TotalLenCMSA(cma)); 
	    PutHist(stderr,50,HG); NilHist(HG); 
	} // assert(dmax >= 0.0);
	if(efp){
          i=(Int4) ceil(fabs(dmax)+2); d=dmax/20.0;
	  if(d <= 0.0) d=0.1; // assert(d >= 0.0);
	  HG=Histogram("Column BILD scores (in nats)",-i,i,d);
	  for(blk=1; blk <= nBlksCMSA(cma); blk++){
	        for(i=1; i <= LengthCMSA(blk,cma); i++){ IncdHist(RtnBS[blk][i],HG); } 
	  } PutHist(stderr,50,HG); NilHist(HG);
	} return RtnBS;
}

cma_typ	RmWorstColsCMSA(Int4 mincol,ssx_typ *ssx)
// Remove the lowest BILD scoring colums but retain mincol.
{
	cma_typ	xcma,rcma,cma=ssx->RtnCMA(); assert(nBlksCMSA(cma) == 1); 
	Int4	x,i,j,k,n,blk=1,N=LengthCMSA(blk,cma);
	double **BS=BILD_ScoresCMSA(ssx);
	a_type	AB=AlphabetCMSA(cma);

	set_typ	Set=MakeSet(N+3); ClearSet(Set);
	dh_type dH=dheap(N+2,4);
	h_type HG=Histogram("Column contributions to Contextual BILD scores(in nats)",-10,10,0.10);
	for(i=1; i <= N; i++){
		double	d,bs;
#if 0	// 
	   // switch(mode){
		n=NumberDeletionsCMSA(i, cma);	// remove columns with most deletions (n) and, if equal,
		insrtHeap(i,(keytyp)((BS[1][i]*0.000001)-(double)n),dH);  // with lower BILD scores.
#elif 1
		n=NumberDeletionsCMSA(i, cma);	// favor columns with fewer deletions (n) and...
		d=ssx->ContextBildScore(1,i); 
		insrtHeap(i,(keytyp)((d*0.000001)-(double)n),dH);  // delete lower BILD scores.
#elif 1
		n=NumberDeletionsCMSA(i, cma);	// favor columns with fewer deletions (n) and...
		insrtHeap(i,(keytyp)(BS[1][i]+((double)n/100000.0)),dH);  // focus on BILD scores.
#elif 1	// moved below to ssx_typ
		d=ssx->ContextBildScore(1,i); IncdHist(d,HG);
		insrtHeap(i,(keytyp)d,dH);  // focus on BILD scores.
#elif 1	// Contextual BILD scores...
		double frct=0.750;
		UInt4 WtCnts[3][30],sqwt;
		Int4	ins_res=0,nins=NumInsertsCMSA(1,i,cma,ins_res);
		unsigned char	r,del=AlphaCode('x',AB);
		enum context { left=1, right=2 };
		for(j=0; j <= nAlpha(AB); j++) WtCnts[left][j]=WtCnts[right][j]=0;
		Int4	sq,pos,nins=0,ins,ins_len=0,Pos[5],ndel=0;
		for(sq=1; sq <= NumSeqsCMSA(cma); sq++){
			e_type tE=TrueSeqCMSA(sq,cma); sqwt=ssx->RtnSeqWt(sq);

			pos=TruePosCMSA(sq,i,cma);
			unsigned char   R,RR,*isq=SeqPtrCMSA(sq,cma);	// FakeSeq ptr.
			PosSiteCMSA(1,sq,Pos,cma); j = Pos[1]+i-1;
			R=isq[j];  
                  	if(IsDeletedCMSA(sq,j,cma)){ RR=0; ndel++; } 
			else { RR=ResSeq(pos,tE); }
			
		// fprintf(stderr,"%d.%d: %c%d == %c%d?\n",i,sq,AlphaChar(R,AB),j,AlphaChar(RR,AB),pos);
			assert(R == RR);	// confirm that these correspond.
			// pos=FakeToRealCMA(sq,i,cma);
			if(pos == 1){
				WtCnts[left][del] += sqwt;
				// r=ResSeq(pos+1,tE); WtCnts[left][r] += sqwt;
				r=ResSeq(pos+1,tE); WtCnts[right][r] += sqwt;
			} else if(pos == LenSeq(tE)){
				r=ResSeq(pos-1,tE); WtCnts[left][r] += sqwt;
				// r=ResSeq(pos-1,tE); WtCnts[right][r] += sqwt;
				WtCnts[right][del] += sqwt;
			} else {
				r=ResSeq(pos-1,tE); WtCnts[left][r] += sqwt;
				r=ResSeq(pos+1,tE); WtCnts[right][r] += sqwt;
			}
		}
		double	TotalWtSq=ssx->TotalWtSeq();
		double	bld_lt=ssx->BildScore(WtCnts[left])/TotalWtSq;
		double	bld_rt=ssx->BildScore(WtCnts[right])/TotalWtSq;
		double  bs=ssx->BildScore(1,i)/TotalWtSq; 
		// double	d=BS[1][i] + frct*bld_rt + frct*bld_lt;
		double	d=bs + frct*bld_rt + frct*bld_lt;
// fprintf(stderr,"%d: BILD = %.3f; bld_lt = %.3f; bld_rt = %.3f; total = %.3f\n",i,BS[1][i],bld_lt,bld_rt,d);
fprintf(stderr,"%d: BILD = %.3f; bld_lt = %.3f; bld_rt = %.3f; total = %.3f; deletions = %d\n"
			,i,bs,bld_lt,bld_rt,d,ndel);
	if(nins > 0) fprintf(stderr,"    %d insertions\n",nins); 
		IncdHist(d,HG);
		insrtHeap(i,(keytyp)d,dH);  // focus on BILD scores.
#else
		insrtHeap(i,(keytyp)BS[1][i],dH); 
#endif
	} PutHist(stderr,50,HG); NilHist(HG); 
	for(n=N; n > mincol; ){
		assert((i=delminHeap(dH)) != 0); n--; AddSet(i,Set);
	} Nildheap(dH);
	for(blk=1; blk <= nBlksCMSA(cma); blk++) free(BS[blk]); free(BS);	
	char	*operation; NEW(operation,N+5,char); operation[0]='E'; operation[N+1]='E'; 
	for(j=1; j <= N; j++){
	    if(MemberSet(j,Set)) operation[j]='d'; else operation[j]='m';
	} rcma=RemoveTheseColumnsCMSA(operation,cma); free(operation); NilSet(Set);
	return rcma;
}

e_type	MKCsqSubCMSA(set_typ Set,cma_typ cma)
//  Return a consensus sequence for entire sequence;
{
        ss_type data=DataCMSA(cma);
        st_type S=SitesCMSA(cma);
        Int4    best,s,r,t,j,*cnts,d,lenM;
        char    *rtnseq;
        a_type  AB = AlphabetCMSA(cma);

   NEW(rtnseq,TotalLenCMSA(cma)+3,char);
   for(s=1,t=1; t <= nBlksCMSA(cma); t++){
    lenM=LengthCMSA(t,cma);
    // 2. generate consensus sequence.
    Int4 N=NSeqsSeqSet(data);
    for(d=0; d < lenM; d++){
        cnts = GetSiteFreq(S,Set,t,d);
        for(r=best=0,j=1; j <= nAlpha(AB); j++){
           if(cnts[j] > best) { best = cnts[j]; r = j; }
        } free(cnts);
        //fprintf(stderr,"best = %d; r = %c\n", best, AlphaChar(r,AB));
        rtnseq[s] = AlphaChar(r,AB); s++;
    }
   } e_type E=StringToSeq(rtnseq+1, "consensus seq", 1, AB);
   free(rtnseq); return E;
}

cma_typ	MkSubCMSA(set_typ Set, BooLean	AddCsq, cma_typ cma)
{
	a_type  AB=AlphabetCMSA(cma);
	Int4	i,j,N,sq,blk,*sites;
	e_type *EList; NEW(EList,CardSet(Set) + 5,e_type);
	j=0;
	if(AddCsq){ j++; EList[j]=MKCsqSubCMSA(Set,cma); }
	for(sq=1; sq <= NumSeqsCMSA(cma); sq++){
		if(!MemberSet(sq,Set)) continue;
		j++; EList[j]=CopySeq(TrueSeqCMSA(sq,cma));
	}
	ss_type data=Array2SeqSet(EList,j,NameCMSA(cma),AB); // data absorbs EList!!! 
	gss_typ *gss=gssCMSA(cma);
   	cma_typ rcma=EmptyCMSA(nBlksCMSA(cma),LengthsCMSA(cma),data,gss->GapOpen(),
		gss->GapExtend(),PerNatsCMSA(cma),0,0);
	gsq_typ	*gsq,*gsq2;
	j=0;
	if(AddCsq){
		gsq2 = new gsq_typ[1];
		Int4 *pos; NEW(pos,nBlksCMSA(rcma) + 5, Int4); pos[1]=1;
		for(i=1; i < nBlksCMSA(rcma); i++) pos[i+1]=pos[i] + LengthCMSA(i,rcma);
		gsq2->initialize(nBlksCMSA(rcma),LengthsCMSA(rcma),pos,TrueSeqCMSA(1,rcma));
		j++; ReplaceCMSA(j,gsq2,rcma); // replace sequence s in CMSA & fmodel.
		for(blk=1; blk <= nBlksCMSA(rcma); blk++) AddSiteCMSA(blk,j,pos[blk],rcma); free(pos);
	}
	for(sq=1; sq <= NumSeqsCMSA(cma); sq++){
		if(!MemberSet(sq,Set)) continue;
		gsq=gsqCMSA(sq,cma);
		gsq2 = new gsq_typ[1]; gsq2->initialize(*gsq);
		// gsq2= new gsq_typ(*gsq);	// core dumps due to array of [1] convention in gss_typ.
		// gsq->Put(stderr,AB); gsq2->Put(stderr,AB);
		j++; ReplaceCMSA(j,gsq2,rcma); // replace sequence s in CMSA & fmodel.
		sites=GetPosSitesCMSA(sq,cma);
		for(blk=1; blk <= nBlksCMSA(rcma); blk++) AddSiteCMSA(blk,j,sites[blk],rcma);
		free(sites);
	} return rcma;
}

#if 0
Int4	PutNewCMSA(FILE *fp,set_typ Set, BooLean put_csq, cma_typ cma)
// put in aligned sequence format (so can be used by other routines
{
	Int4	d,i,j,lenM,LenStr,AlignLen,J,J0,m,s,t,end;
	// Int4	netlen;
	Int4	insertions,deletions,matches,interblks,endseqs;
	st_type	S=SitesCMSA(cma);
	a_type	A=SeqSetA(DataCMSA(cma));
	e_type	fakeE,E;
	BooLean	**null;
	gss_typ& gss=*gssCMSA(cma);
	Int4	site,number=gss.NumSeq();
	Int4	NetNum;

	Int4	Num;
	if(Set == 0) Num = number;
	else for(Num=0,J=1; J <= number; J++) if(MemberSet(J,Set)) Num++;
	// if(cma->FullSeq) assert(J0 == PutFullSeqCMSA(fp,skip, cma));
	NetNum=Num;
	if(put_csq) NetNum++;
        fprintf(fp,"[%d_(%d)=%s(%d){go=%d,gx=%d,pn=%.1f,lf=%d,rf=%d}",
		cma->Level,nBlksCMSA(cma),NameCMSA(cma),NetNum,
		gss.GapOpen(),gss.GapExtend(),gss.PerNats(),
        	gss.LeftFlank(),gss.RightFlank());
	if(cma->alpha > 0.0 && cma->alpha <= 1.0){
            fprintf(fp,";BPPS=(%c,%d,%d:%.4f):\n",cma->set_mode,cma->A0,cma->B0,cma->alpha);
	} else { fprintf(fp,":\n"); }

	// 2a. Print column indicators:
	// null = NullCMSA(cma); null[m][s] == ?
	for(m=1; m <= nBlksCMSA(cma); m++){
	   fprintf(fp,"(%d)",LengthCMSA(m,cma));
	   fm_type fm=ModelCMSA(m,cma);
	   for(s=1; s <= LengthCMSA(m,cma); s++){
		if(NullSiteFModel(s,fm)) fprintf(fp,".");
		else fprintf(fp,"*");
	   } fprintf(fp,"\n");
	} fprintf(fp,"\n");
	fprintf(fp,"\n");
	//******  Print consensus sequence here:
	if(put_csq){
	  fprintf(fp,"$0=%d(%d):\n>%s consensus seq\n",
			TotalLenCMSA(cma),TotalLenCMSA(cma),NameCMSA(cma));
	  fprintf(fp,"{()");
	  for(m=1; m <= nBlksCMSA(cma); m++){
	     char    *csq=ConsensusSeqCMSA(m,cma);
	     fprintf(fp,"%s()",csq); free(csq);
	  }
	  fprintf(fp,"}*\n\n");
	}
	// create_gsq_cmsa(cma);	// create gsq objects for cma if not present...
	for(J0=J=1; J <= number; J++){ 
	   if(Set && !MemberSet(J,Set)) continue;
	   gsq_typ *gsq=gsqCMSA(J,cma);
	   assert(gsq != 0);	// this should not happen...
	   Int4    *sites=GetPosSitesCMSA(J,cma);
	   gsq->Put_cma_format(fp,J0,nBlksCMSA(cma),sites,LengthsCMSA(cma),A);
	   free(sites);
	   J0++;  // always number consecutively 
	} fprintf(fp,"_%d].\n",cma->Level); 
	return NetNum;
}
#endif

double  ResidueDiversityCMSA3(set_typ Set, cma_typ cma)
// what is the average number of residue types seen at each position?
{
     register sst_typ   sst;
     register Int4      k,i,sum,s;
     register unsigned char *seq;

     assert(nBlksCMSA(cma) == 1);
     for(sum=0,k=LengthCMSA(1,cma); k > 0; k--){
          for(sst=0,i=NumSeqsCMSA(cma); i > 0; i--){
              if(MemberSet(i,Set)){
		seq=SeqPtrCMSA(i,cma);
                // seq=SeqSeqSet(i,DataCMSA(cma));
		s = SitePos(1,i,1,SitesCMSA(cma))+k-1;
                sst=UnionSset(sst,SsetLet(seq[s]));
              }
          } for(i=nAlpha(AlphabetCMSA(cma)); i > 0; i--) if(MemSset(i,sst)) sum++;
     } return (double)sum/(double) TotalLenCMSA(cma);
}

double  AveEntropyCMSA(set_typ	Set, cma_typ cma)
// The maximum = ln 20 = 2.9957 (or ln 21 = 3.044).
{
        Int4    t,k,sq,r,n;
        double  p,sum,entropy,cnt[50],tcnt;
        a_type AB = AlphabetCMSA(cma);

     for(sum=0.0,n=0,t=1; t <= nBlksCMSA(cma); t++){
        for(k=1; k <= LengthCMSA(t,cma); k++){
          n++;
          for(tcnt=0.0,r=0; r <= nAlpha(AB); r++) cnt[r]=0.0;
          for(sq=1; sq <= NumSeqsCMSA(cma); sq++){
	      if(Set && !MemberSet(sq,Set)) continue;
              r=ResidueCMSA(t,sq,k,cma);
              if(r){ cnt[r]+=1.0; tcnt +=1.0; }
          }
          for(entropy=0.0,r=0; r <= nAlpha(AB); r++){
                p = cnt[r]/tcnt; 
                if(p > 0.0) entropy += -p * log(p);
          } sum += entropy;
        }
     } return (sum/(double)n);
}

set_typ	AutoPurgeCMSA(Int4 s, set_typ SubGrp,char *name,cma_typ cma)
// NodeCMA[1]=AutoPurgeCMSA(set_typ SubGrp,Hpt->SetName(1),TrueMainCMA)
{
	cma_typ	rcma;
	Int4	percent_ident,N;
	set_typ RepSet=0;
        a_type AB = AlphabetCMSA(cma);

	// double	d=ResidueDiversityCMSA3(SubGrp,cma);
	double	d=AveEntropyCMSA(SubGrp,cma),dd,D;
	d = exp(-d);	// expected residue frequency.
	// dd=ComputSeqWtsCMSA(cma,SubGrp); D=100*dd/(double) CardSet(SubGrp);
	// value == 0 (most conserved) to 1.0 (unconserved).
	N = CardSet(SubGrp);
	if(N > 100000){
	    if(d < 0.50) percent_ident=70; else if(d < 0.75) percent_ident=80; else percent_ident=85;
	} else if(N > 20000){	// 20,000..100,000
	    if(d < 0.50) percent_ident=80; else if(d < 0.75) percent_ident=85; else percent_ident=90;
	} else if(N > 4000){	// 4,000..20,000
	    if(d < 0.50) percent_ident=85; else if(d < 0.75) percent_ident=90; else percent_ident=95;
	} else if(N > 500){	// 500..4,000
	    if(d < 0.50) percent_ident=90; else if(d < 0.75) percent_ident=95; else percent_ident=98;
	} else if(N > 200){	// 200..500
	    if(d < 0.50) percent_ident=0; else if(d < 0.75) percent_ident=90; else percent_ident=95;
	} else { percent_ident=0; }
	if(percent_ident == 0) RepSet=CopySet(SubGrp); 
	else RepSet=RtnFastRepSetCMSA(0,percent_ident,SubGrp,cma);
#if 0
	fprintf(stderr,"%d: %s = %d; ave res freq = %.3f (-U%d); weighted #seqs =%.1f (%.1f%c).\n",s,name,
                CardSet(SubGrp),d,percent_ident,dd,D,'%');
#else
	fprintf(stderr,"%d: %s = %d seqs; ave res freq = %.3f (-U%d); retained seqs =%d.\n",
		s,name,N,d,percent_ident,CardSet(RepSet));
#endif
        // rcma=MkSubCMSA(RepSet,FALSE,cma); NilSet(RepSet);
	// return rcma;
	return RepSet;
}

cma_typ	ExtendColumnsCMSA(cma_typ cma, char dms_mode, double SqWtAdj, double minFrq2Del)
// this routine has problems!!
{
	assert(minFrq2Del > 0.0 && minFrq2Del < 1.0);
	double minFrq2Add=1.0 - minFrq2Del;
	cma_typ xcma=cma,tmp_cma=0;
	do {	// add columns.
	   if(tmp_cma != 0){ if(xcma != cma) NilCMSA(xcma); xcma=tmp_cma; tmp_cma=0; }
	   tmp_cma=ModifyColumnsCMSA(xcma, dms_mode, SqWtAdj,minFrq2Add,'A');
	} while(tmp_cma != 0);
#if 0
	tmp_cma=0;
	do {	// delete columns...
	   if(tmp_cma != 0){ if(xcma != cma) NilCMSA(xcma); xcma=tmp_cma; tmp_cma=0; }
	   tmp_cma=ModifyColumnsCMSA(xcma, dms_mode, SqWtAdj,minFrq2Del,'D');
	} while(tmp_cma != 0);
#endif
	if(xcma == cma) return 0;
	else return xcma;
}

cma_typ	ContractColumnsCMSA(cma_typ cma, char dms_mode, double SqWtAdj, double minFrq2Del)
{
	assert(minFrq2Del > 0.0 && minFrq2Del < 1.0);
	cma_typ xcma=cma,tmp_cma=0;
	do {	// delete columns...
	   if(tmp_cma != 0){ if(xcma != cma) NilCMSA(xcma); xcma=tmp_cma; tmp_cma=0; }
	   tmp_cma=ModifyColumnsCMSA(xcma, dms_mode, SqWtAdj,minFrq2Del,'D');
	} while(tmp_cma != 0);
	if(xcma == cma) return 0;
	else return xcma;
}

cma_typ	ModifyColumnsCMSA(cma_typ cma, char dms_mode, double SqWtAdj,double frq_cutoff,char mode)
// If mode=='A' (add) then convert conserved insert regions to aligned columns 
// If mode=='D' (delete) then convert columns with lots of deletions into insertions.
// if deletions are on the ends, then trim back the alignment?
{
print_error("ModifyColumnsCMSA() routine has problems! Exiting...");
	Int4	blk,i,j,RmLeft,RmRight,length=0,start=0,theBlk=0,ins_res,n_ins,n_del,begin;
	Int4    aa_per_io=1000,aa_per_do=1000,exp_ie=3,exp_de=1;
	double	DD,*deleted,max_frq=0.0,pn=PerNatsCMSA(cma); pn=1000.0;
	ssx_typ *ssx=new ssx_typ(aa_per_io,aa_per_do,exp_ie,exp_de,pn,cma,dms_mode); ssx->SetPriorWt(1);
	double TotalWtSq=ssx->TotalWtSeq();

	if(SqWtAdj > 0) ssx->SetSqWtAdjust(SqWtAdj);
	for(blk=1; blk <= nBlksCMSA(cma); blk++){
		if(mode=='A') begin=1; else begin=1;	// Disallow Deletions on the ends.
	        for(i=begin; i < LengthCMSA(blk,cma); i++){
		    n_del=NumDeletionsCMSA(blk, i,cma);
		    n_ins=NumInsertsCMSA(blk, i,cma,ins_res);
		    if(mode == 'A'){	// extend columns.
		      DD=(double)n_ins/(double)NumSeqsCMSA(cma);
		      if(DD >= frq_cutoff && DD > max_frq){
			max_frq=DD; start=i; theBlk=blk;
		    	DD= (double)ins_res/(double)n_ins; length=(Int4) ceil(DD);
		      }
		    } else {		// delete columns.
		      DD=(double)n_del/(double)NumSeqsCMSA(cma);
		      if(DD >= frq_cutoff && DD > max_frq){
			max_frq=DD; start=i; theBlk=blk; length=1;
		      }
		    }
#if 1
	            fprintf(stderr,"%d.%d%c RE=%.3f; ins=%d (%d) del=%d (%.1f%c).\n",
		    	blk,i,mode,ssx->RelEntropy(blk,i),n_ins,ins_res,n_del,
			100.0*FractDeletionsCMSA(blk,i,cma),'%');
		    if(ins_res > 10) ssx->PutPenalties(stderr,blk,i);
#endif
		} fprintf(stderr,"\n");
	}
	cma_typ rcma=0;
	if(max_frq > frq_cutoff){
	   if(mode=='A') rcma=InsertColumnsCMSA(cma,theBlk,start,length);
	   else if(start > 1) rcma=ConvertColsToInsertsCMSA(cma, theBlk, start,start);
	   else {
		rcma=TrimBlkCMSA(cma,theBlk,1,0,3);
		if(rcma==0) print_error("TrimBlkCMSA(): too many columns removed from alignment");
	   }
	} delete ssx;
	return rcma;
}

#if 0	// old code...
cma_typ	ModifyColumnsCMSA_OLD(cma_typ cma, char dms_mode, double SqWtAdj,double frq_cutoff,char mode)
{
	Int4	blk,i,j,RmLeft,RmRight,length=0,start=0,theBlk=0,ins_res,n_ins,n_del,begin;
	Int4    aa_per_io=1000,aa_per_do=1000,exp_ie=3,exp_de=1;
	double **SubLLR,d,D,FD,DD,dd,dmax,*deleted,RE,BS,BSps,max_frq=0.0;
	double	pn=PerNatsCMSA(cma); fprintf(stderr,"PerNatsCMSA = %g\n",pn); pn=1000.0;
	ssx_typ *ssx=new ssx_typ(aa_per_io,aa_per_do,exp_ie,exp_de,pn,cma,dms_mode); ssx->SetPriorWt(1);
	double TotalWtSq=ssx->TotalWtSeq();

	if(SqWtAdj > 0) ssx->SetSqWtAdjust(SqWtAdj);
	D = ssx->GapMap(0.0,'V'); dd = ssx->BildLLR('M');
	NEWP(SubLLR,nBlksCMSA(cma) +3, double);
	for(dmax=-9999999,blk=1; blk <= nBlksCMSA(cma); blk++){
	   	NEW(SubLLR[blk],LengthCMSA(blk,cma) +3, double);
	        fprintf(stderr,"Rm blk %d: %.1f nats; oldmap = %.1f; bild LLR=%.1f.\n",
			blk,D - ssx->GapMap(0.0,'V',blk),FieldRelMapCMSA(cma,blk),dd);
		if(mode=='A') begin=1; else begin=2;	// Disallow Deletions on the ends.
	        for(i=begin; i < LengthCMSA(blk,cma); i++){
		    d=ssx->GapMap(0.0,'V',blk,i); SubLLR[blk][i]=D-d;
		    RE=ssx->RelEntropy(blk,i); BS=ssx->BildScore(blk,i); BSps=BS/TotalWtSq;
		    n_del=NumDeletionsCMSA(blk, i,cma);
		    n_ins=NumInsertsCMSA(blk, i,cma,ins_res);
		    FD=FractDeletionsCMSA(blk,i,cma);

		    if(mode == 'A'){	// extend columns.
		      DD=(double)n_ins/(double)NumSeqsCMSA(cma);
		      if(DD >= frq_cutoff && DD > max_frq){
			max_frq=DD; start=i; theBlk=blk;
		    	DD= (double)ins_res/(double)n_ins; length=(Int4) ceil(DD);
		      }
		    } else {		// delete columns.
		      DD=(double)n_del/(double)NumSeqsCMSA(cma);
		      if(DD >= frq_cutoff && DD > max_frq){
			max_frq=DD; start=i; theBlk=blk; length=1;
		      }
		    }
		    SubLLR[blk][i]=BSps;
		    if(SubLLR[blk][i] > dmax) dmax=SubLLR[blk][i];
#if 1
	            fprintf(stderr,"%d.%d = %.1f nats; bild=%.1f (%.2f npws); RE=%.3f; ins=%d (%d) del=%d (%.1f%c).\n",
			blk,i,D-d,BS,BSps,RE,n_ins,ins_res,n_del,100.0*FD,'%');
		    if(ins_res > 10) ssx->PutPenalties(stderr,blk,i);
#endif
		} fprintf(stderr,"\n");
	}
#if 1
	cma_typ rcma=0;
	for(blk =1; blk <= nBlksCMSA(cma); blk++) free(SubLLR[blk]); 
	if(max_frq > frq_cutoff){
	   if(mode=='A') rcma=InsertColumnsCMSA(cma,theBlk,start,length);
	   else rcma=ConvertColsToInsertsCMSA(cma, theBlk, start,start);
	} delete ssx;
	free(SubLLR); 
	return rcma;
#else
	double *AveDel,**Deleted = ssx->FractionDeleted();
	NEW(AveDel,nBlksCMSA(cma) +3, double);
	i=(Int4) ceil(dmax+2);
	d=(Int4) dmax/20.0;
	h_type HG=Histogram("Column contributions to LLR (in nats)",-i,i,d);
	for(blk=1; blk <= nBlksCMSA(cma); blk++){
	        for(D=0,i=1; i <= LengthCMSA(blk,cma); i++){
		    IncdHist(SubLLR[blk][i],HG);
#if 0
		    fprintf(stderr,"%d.%d: %.2f nats, %.3f%c\n",
				blk,i,SubLLR[blk][i],Deleted[blk][i],'%');
#endif
		    D+=Deleted[blk][i];
		} AveDel[blk]=D/(double)LengthCMSA(blk,cma);
		// fprintf(stderr,"  ave deleted = %.3f%c.\n\n",AveDel[blk],'%');
	} PutHist(stderr,50,HG); NilHist(HG); // exit(1);
	delete ssx;
	free(AveDel);
	for(blk =1; blk <= nBlksCMSA(cma); blk++) { free(SubLLR[blk]); free(Deleted[blk]); }
	free(SubLLR); free(Deleted);
#endif
}
#endif

cma_typ AddThisCsqCMSA(e_type CsqE, cma_typ cma)
{
	a_type  AB=AlphabetCMSA(cma);
	FILE *fp=tmpfile(); PutWithCsqCMSA(fp,CsqE,cma); rewind(fp);
	cma_typ	rcma=ReadCMSA(fp,AB); fclose(fp);
	return rcma;
}

Int4	PutWithCsqCMSA(FILE *fp,e_type CsqE, cma_typ cma)
// put in aligned sequence format (so can be used by other routines
// == common language for alignments.
// == Psiblast too where master sequence defines a single block template.
{
	Int4	d,i,j,lenM,LenStr,J,J0,m,r,s,t;
	char	c;

	//=========== 1. Get input information: ==========
	gss_typ& gss=*gssCMSA(cma);
	Int4	site,number=gss.NumSeq(),NetNum=number;
	if(CsqE){ assert(LenSeq(CsqE) == TotalLenCMSA(cma)); NetNum++; }

	//=========== 2. Print priliminary information: ==========
	if(cma->FullSeq) assert(1 == PutFullSeqCMSA(fp,0,cma));
	fprintf(fp,"[%d_(%d)=%s(%d){go=%d,gx=%d,pn=%.1f,lf=%d,rf=%d}",
		cma->Level,nBlksCMSA(cma),NameCMSA(cma),NetNum,
		gss.GapOpen(),gss.GapExtend(),gss.PerNats(),
        	gss.LeftFlank(),gss.RightFlank());
	if(cma->alpha > 0.0 && cma->alpha <= 1.0){
            fprintf(fp,";BPPS=(%c,%d,%d:%.4f):\n",cma->set_mode,cma->A0,cma->B0,cma->alpha);
	} else { fprintf(fp,":\n"); }

	//=========== 3. Print column indicators: ==========
	for(m=1; m <= nBlksCMSA(cma); m++){
	   fprintf(fp,"(%d)",LengthCMSA(m,cma));
	   fm_type fm=ModelCMSA(m,cma);
	   for(s=1; s <= LengthCMSA(m,cma); s++){
		if(NullSiteFModel(s,fm)) fprintf(fp,"."); else fprintf(fp,"*");
	   } fprintf(fp,"\n");
	} fprintf(fp,"\n");
	fprintf(fp,"\n");

	//=========== 3. Print consensus sequence: ==========
   	a_type AB=AlphabetCMSA(cma);
	if(CsqE){
	  fprintf(fp,"$0=%d(%d):\n>%s consensus seq\n",
			TotalLenCMSA(cma),TotalLenCMSA(cma),NameCMSA(cma));
	  fprintf(fp,"{()");
	  for(j=0,m=1; m <= nBlksCMSA(cma); m++){
	     for(s=1; s <= LengthCMSA(m,cma); s++){
		j++; r=ResSeq(j,CsqE); fprintf(fp,"%c",AlphaChar(r,AB));
	     } fprintf(fp,"()");
	  } fprintf(fp,"}*\n\n");
	}

	//=========== 4. Print aligned sequences: ==========
	create_gsq_cmsa(cma);	// create gsq objects for cma if not present...
	for(J0=J=1; J <= number; J++){ 
	   gsq_typ *gsq=gsqCMSA(J,cma);
	   assert(gsq != 0);	// this should not happen...
	   Int4    *sites=GetPosSitesCMSA(J,cma);
	   gsq->Put_cma_format(fp,J0,nBlksCMSA(cma),sites,LengthsCMSA(cma),AB);
	   free(sites);
	   J0++;  // always number consecutively 
	} fprintf(fp,"_%d].\n",cma->Level); 
	return NetNum;
}

char	*FixTmplOperation(char *operation, e_type CsqE, Int4 &Start, Int4 &TrimNt, 
		Int4 &TrimCt)
// 1. Add insertions on ends...
// 2. Remove Deletions at start...(for now).
// 3. Trim N- and C-terminal ends if necessary.
// Use this for hieraln program as temporary fix for issues...
{
	char	c,*Op=0,*op=0; 
	Int4	i,j,o,d,x,y,z,overNt=0,overCt=0,delNt=0,delCt=0;
	BooLean	flag=TRUE;
	TrimNt=0; TrimCt=0;
	BooLean	debug=FALSE;

	// 1. Count overhangs and deletions on either end...
	overNt=Start-1; 
	assert(overNt >= 0);
	// for(x=1,z=1; x < Start; x++){ z++; Op[x]='i'; overNt++; }
	for(z=Start-1, y=1 ; (c=operation[y]) != 'E'; y++){
	    if(c != 'D' && c != 'd'){ z++; flag=FALSE; delCt=0;}
	    else { if(flag) delNt++; else delCt++; }
	} while(z < LenSeq(CsqE)){ z++; overCt++; }
	if(debug){
	   fprintf(stderr,"Start=%d; Len=%d; z=%d; overNt=%d; delNt=%d; overCt=%d; delCt=%d\n",
		Start,LenSeq(CsqE),z,overNt,delNt,overCt,delCt);
	}

	// 2. Remove Deletions at start...(for now).
	NEW(Op,strlen(operation)+overNt+overCt+8,char); Op[0]='E';
	// 2a. fix N-terminal end issues...
	d=delNt; o=overNt; 
	if(o == 0){		
	   strcpy(Op,operation); 
	   if(d > 0){ Op[1]='M'; Op[delNt+1]='d'; }	// e.g., "EDddmmm..." --> "EMdddmm..."
	} else if(o > 0 && d == 0 && operation[2]=='d'){	// e.g., "(iiiii)EMdddmmm.."
		Op[1]='M';
#if 1	// fix bug with "...mmmddddmE(ii)" --> "...mmmddddIImE" problem
		for(j=2; operation[j]=='d' && o > 0; ){ Op[j]='m'; o--; j++; }
		for(i=j; o > 0; o--,i++) Op[i]='I'; 
		strcat(Op,operation+j);
#endif	// back up to fill in deletions...
	} else if(o > d){	// e.g., "EiiiiiDdddmmm.."
	   if(d == 0){		// e.g., "EiiiiiMmmm.." --> "EMIIIIImmm..."
		for(Op[1]='M',i=2; o > 0; o--,i++) Op[i]='I'; 
		strcat(Op,operation+2);
	   } else { 		// e.g., "EiiiDdmmm.." --> "EMImmm..."
	     for(Op[1]='M',d--,o--, i=2 ; o > d; o--,i++) Op[i]='I';
	     for( ; d > 0; d--,o--,i++) Op[i]='m';
	     strcat(Op,operation+(delNt+1));
	   }
	} else if(o <= d){	// e.g., "EiiDddmmm..." overNt > 0 && delNt > 0
	     for(Op[1]='M', o--, x=2 ; o > 0; d--,o--,x++) Op[x]='m';
	     strcat(Op,operation+(overNt+1));
	} else print_error("FixTmplOperation(): this should not happen");
	// fprintf(stderr,"Op=%s\n",Op);
	Start=1; // In all cases resets start to 1.

	// 2b. fix C-terminal end issues...
	op=Op; NEW(Op,strlen(operation)+overNt+overCt+8,char); 
	// fprintf(stderr,"Op=%s\n",op);
	strcpy(Op,op);  y = strlen(Op)-2; // y = last operation before 'E'.
	free(op);
	d=delCt; o=overCt; x = y-d;		 // x = last 'm'. 
	assert(Op[x]=='m'); 
	if(o == 0){		
	   if(d > 0){ Op[y]='m'; Op[x]='d'; }	// last m to end: e.g., "...mmmdddE" --> "...mmdddmE"
	} else if(o > d){	
	   if(d == 0){		// e.g., "...mmmEiiii" --> "...mmIIIImE"
#if 1	// fix bug with "...mmmddddmE(ii)" --> "...mmmddddIImE" problem
		for(j=x-1; Op[j]=='d' && o > 0; ){ Op[j]='m'; o--; j--; }
#endif	// back up to fill in deletions...
	        for(i=x ; o > 0; o--,i++) Op[i]='I'; 
		Op[i]='m'; i++;
	   } else { 			// e.g., "...mmmdddEiiii" --> "...mmmImmmE"
	     for(i=x+1 ; o > d; o--,i++) Op[i]='I';
	     for( ; d > 0; d--,o--,i++) Op[i]='m';
	   } Op[i]='E'; i++; Op[i]=0;
	} else if(o <= d){	// e.g., "...mdddEii" --> "...mdmmE" (o > 0 && d > 0).
	     // fprintf(stderr,"Op[y]='%c'; Op[y+1]='%c'\n",Op[y],Op[y+1]);
	     for(i=y+1,j=1; j <= o; j++){ i--; Op[i]='m'; }	// backup... now d == o
	     i=y+1;
	     // for(i=y+1; o > 0; o--,i++) Op[i]='m'; 
	     Op[i]='E'; i++; Op[i]=0;
	} else print_error("FixTmplOperation(): this should not happen");
	return Op;
}

Int4	NumAlnSeqsCMSA(cma_typ cma)
{
	Int4    n,N=NumSeqsCMSA(cma),sq;
	for(n=0,sq=1; sq <= N; sq++){ if(SeqIsInCMSA(sq,cma)) n++; }
	return n;
}

void	PutSubGroupCsqScores(FILE *fp, set_typ Set, cma_typ cma)
{
       double	inc;
       a_type	AB=AlphabetCMSA(cma);
       Int4	j,sq,score;
       e_type	csqE=GetSeqAsCsqCMSA(Set,cma);
       for(score=0,j=1; j <= LenSeq(csqE); j++){
                 unsigned char r=ResSeq(j,csqE);
                 if(r) score += valAlphaR(r,r,AB);
       }
       inc=ceil((double)score/40.0);
       h_type  HG=Histogram("in-set sequence scores",-200,score,inc);
       h_type  xHG=Histogram("out-of-set sequence scores",-200,score,inc);
       for(sq=1; sq <= NumSeqsCMSA(cma); sq++){
	  if(!SeqIsInCMSA(sq,cma)) continue; 
          score=PseudoAlnScoreSqToCMSA(csqE,sq,cma);
	  if(MemberSet(sq,Set)){ IncdHist((double)score,HG); }
	  else { IncdHist((double)score,xHG); }
       } NilSeq(csqE);
       PutHist(fp,50,HG); NilHist(HG); PutHist(fp,50,xHG); NilHist(xHG);
}

void	PutSubGroupFootPrints(FILE *fp, set_typ Set, cma_typ cma)
// 
{
        h_type  HG=Histogram("in-set domain footprints",0,1000,10);
        h_type  xHG=Histogram("out-of-set domain footprints",0,1000,10);
	Int4    x,y,Blk=1,N=NumSeqsCMSA(cma),sq;
	for(sq=1; sq <= N; sq++){
	   if(!SeqIsInCMSA(sq,cma)) continue; 
	   x=TruePosCMSA(sq,Blk,1,cma);
	   y=TruePosCMSA(sq,Blk,LengthCMSA(1,cma),cma);
	   if(MemberSet(sq,Set)) IncdHist((double)(y-x+1),HG);
	   else IncdHist((double)(y-x+1),xHG);
	} PutHist(fp,50,HG); NilHist(HG); PutHist(fp,50,xHG); NilHist(xHG);
}

set_typ	*RtnTargetSizeSetsCMSA(Int4 &numSets, Int4 &percent_ident, cma_typ CMA, double MaxFraction)
// cluster sequences into sets with no more than MaxFraction in any one set.
{
	Int4    i,j,k,N=NumSeqsCMSA(CMA);
	double	dd;
        set_typ *Ordered,*CloseSq=0,InSet=MakeSet(N+4); FillSet(InSet);
	assert(MaxFraction <= 1.0);
        do {
            CloseSq=RtnFastClustersCMSA(numSets,percent_ident,InSet,CMA);
            if(percent_ident > 98) break;
            for(j=1; j <= numSets; j++){
                dd=(double) CardSet(CloseSq[j])/(double)N;
                if(dd > MaxFraction){
                  // fprintf(stderr,"%d%c %d seqs: ",percent_ident,'%',N);
                  for(k=1; k <= numSets; k++){
                        // fprintf(stderr,"%d ",CardSet(CloseSq[k]));
                        NilSet(CloseSq[k]);
                  } // fprintf(stderr,"\n");
                  free(CloseSq); CloseSq=0; percent_ident++; break;
                }
            }
        } while(CloseSq==0);
	dh_type dH=dheap(numSets+2,4);
	for(i=1; i <= numSets; i++) insrtHeap(i,(keytyp)-CardSet(CloseSq[i]),dH);
	NEW(Ordered,numSets+4,set_typ);
	for(j=1; !emptyHeap(dH); j++){
		assert((i=delminHeap(dH)) != 0); Ordered[j]=CloseSq[i];
                // fprintf(stderr,"%d ",CardSet(CloseSq[i]));
        } // fprintf(stderr,"\n"); 
	Nildheap(dH); free(CloseSq); NilSet(InSet);
	return Ordered;
}

set_typ	*RtnFastClustersCMSA(Int4 &numSets, Int4 percent_ident,set_typ InSet,cma_typ cma)
// return cliques of sequences within cma with percent_ident.
{
	ds_type sets;
	Int4	b,i,j,k,N=NumSeqsCMSA(cma),score,nblk=nBlksCMSA(cma);
	Int4	s,sb,bl,si,sj,i8,pos[3],seti,setj;
        a_type  AB=AlphabetCMSA(cma);
	unsigned char	r;
	Int4	start=1;

	if(percent_ident < 0){ 	// then keep first sequence in cma output file...
		start=2; percent_ident = -percent_ident;
	} assert(percent_ident > 0 && percent_ident <= 100);

	// [0,1,2,3, 4,5,6,7] = 8 residues at a time.
	assert(sizeof(UInt4) != sizeof(UInt8)); // make sure this is a 64 bit machine...
	assert(sizeof(UInt8)==8); // make sure this is a 64 bit machine...
	char	*bsq; 
	Int4	TotLen=TotalLenCMSA(cma);
	double	d=(double) (TotLen*((double)(100-percent_ident)/100.0));
	Int4	MaxMisMatch= (Int4) floor(d); 
	double	D=(double) TotLen/8.0;
	Int4	NumBlks=(Int4) ceil(D); 

	//================ 0. Store each fake sequence in 8 byte arrays. ===============
	UInt8	**SQ; NEWP(SQ,N+3,UInt8);
	for(i = start; i <= N; i++){
  	  if(!MemberSet(i,InSet)) continue;
	  NEW(SQ[i],(TotLen/7)+2,UInt8);
#if 0
	  unsigned char	*isq = SeqPtrCMSA(i,cma);
	  for(bl=sb=0,b=1; b <= nblk; b++){
	    PosSiteCMSA(b,i,pos,cma); si=pos[1];
	    for(s=1; s <= LengthCMSA(b,cma); s++,si++,sj++){
		if(sb%8==0){ bsq=(char *) &SQ[i][bl]; bl++; sb=0; } 
		bsq[sb] = AlphaChar(isq[si],AB); sb++;
	    }
	  }
#else
	  for(bl=sb=0,b=1; b <= nblk; b++){
	    for(s=1; s <= LengthCMSA(b,cma); s++){
		if(sb%8==0){ bsq=(char *) &SQ[i][bl]; bl++; sb=0; } 
		r=ResidueCMSA(b,i,s,cma);
		bsq[sb] = AlphaChar(r,AB); sb++;
	    }
	  } // fprintf(stderr,"%d: %s\n",i,SQ[i]);
#endif
	}
	// fprintf(stderr,"NumBlks = %d; TotLen= %d; MaxMisMatch=%d\n",NumBlks,TotLen,MaxMisMatch);

	//============ 1. Cluster sequences into sets at the percent identity cutoff. =============
	Int4 total = TotalLenCMSA(cma);
	
	Int4	*edges; NEW(edges,N+3,Int4);	// number of edges out of each node.
	sets = DSets(N);
	for(i = start; i < N; i++){
  	  if(!MemberSet(i,InSet)) continue;
	  seti=findDSets(i,sets);
	  for(j=i+1; j <= N; j++){
  	     if(!MemberSet(j,InSet)) continue;
	     setj=findDSets(j,sets);
	     if(seti != setj){
	      UInt8 *SQ_i=SQ[i],*SQ_j=SQ[j];
	      for(score=0,bl=0,sb=0; bl <= NumBlks; bl++){
		if(SQ_i[bl] != SQ_j[bl]){
		   score += RtnMisMatch((char *)&SQ_i[bl],(char *)&SQ_j[bl]);
		   if(score > MaxMisMatch) break;
		}
	      }
	      if(score <= MaxMisMatch)
		{ edges[i]++; edges[j]++; seti=linkDSets(seti,setj,sets); }
	     }
	  }
	} 
	for(i = start; i <= N; i++) if(SQ[i]) free(SQ[i]);  free(SQ);

	//============== 2. Create an array of sequence sets. ===========
	double	bestprob,var;
	Int4	best;
	set_typ *Set=0,SubSet=MakeSet(SetN(InSet)+1); ClearSet(SubSet);
	NEW(Set,N+3,set_typ);
        for(s=0,i=1; i <= N; i++){
  	  if(!MemberSet(i,InSet)) continue;
	  seti=findDSets(i ,sets);
	  if(i != seti) continue;	// skip non-canonical sequences
	  best=i;
	  // bestprob=GetTotalProbCMSA(i,cma);
	  bestprob=(double)edges[i];
	  ClearSet(SubSet); AddSet(i,SubSet); 
	  for(j=1; j <= N; j++){
  	        if(!MemberSet(j,InSet)) continue;
		if(j!=i){
		  setj=findDSets(j,sets);
		  if(seti == setj) AddSet(j,SubSet);
		}
	  } s++;
	  Set[s]=CopySet(SubSet);
	} NilSet(SubSet); NilDSets(sets); free(edges);
	numSets=s;
	return Set;
}

void	PutCsqScoresCMSA(FILE *fp,set_typ set,cma_typ cma)
{
       double RelStdDev,inc;
       a_type  AB=AlphabetCMSA(cma);
       Int4 j,sq,score;
       e_type  cE=MkConsensusCMSA(cma);
       for(score=0,j=1; j <= LenSeq(cE); j++){
                 unsigned char r=ResSeq(j,cE);
                 if(r) score += valAlphaR(r,r,AB);
       }
       inc=ceil((double)score/25.0);
       h_type  HG=Histogram("sequence scores",-200,score,inc);
       for(sq=1; sq <= NumSeqsCMSA(cma); sq++){
	  if(!MemberSet(sq,set)) continue;
          score=PseudoAlnScoreSqToCMSA(cE,sq,cma);
          IncdHist((double)score,HG);
       }
       RelStdDev=100.0*sqrt(VarianceHist(HG))/MeanHist(HG);
       // if(RelStdDev > 0.0 && RelStdDev < 10.0)
       fprintf(fp,"%s %.2f\n",NameCMSA(cma),RelStdDev);
       NilSeq(cE); PutHist(fp,50,HG); NilHist(HG);
       fprintf(fp,"RSD: %.2f\n",RelStdDev);
}

Int4    NumberDeletionsCMSA(Int4 Column, cma_typ cma)
// return the number of sequences with deletions at column position.
{
    if(nBlksCMSA(cma) != 1) print_error("NumberDeletionsCMSA() requires a single block");
    Int4        sq,hits,sq_hits,s,pos[4],len=LengthCMSA(1,cma),N=NumSeqsCMSA(cma);
    gss_typ *gss=gssCMSA(cma);
    a_type  A=AlphabetCMSA(cma);
    if(Column < 1 || Column > len) print_error("NumberDeletionsCMSA() column out of range!");

    for(hits=0,sq=1; sq<=N; sq++){
        if(IsDeletedCMSA(1,sq,Column,cma)){ hits++; }
    }
    return hits;
}

Int4	RmGappyColumnsCMSA(double cutoff, cma_typ &cma)
// remove columns with fraction 'cutoff' or more of deletions (e.g., 0.50 deletions or more)
// return new length of cma; call using 'iron' option
{
    Int4	Len,s,N=NumSeqsCMSA(cma),i,j,total=0,cycle;
    BooLean	*Delete;

#if 1	// fixes problem with inserting columns!
    ExtendFakeToRealCMSA(cma);
#endif
    NEW(Delete, LengthCMSA(1,cma)+3, BooLean);
    // for(s=2 ; s < LengthCMSA(1,cma); s++)	// Don't remove columns on ends!!
    for(s=1 ; s <= LengthCMSA(1,cma); s++)	// Don't remove columns on ends!!
    {
	Int4 n=NumberDeletionsCMSA(s, cma);
	double d=(double)n/(double)N;
	if(d >= cutoff){
		// fprintf(stderr,"%d: %d deletions (%.2lf)\n",s,n,d);
		total++; Delete[s]=TRUE;
	}
    } Len=LengthCMSA(1,cma);
    for(cycle=0, s = Len; s > 0; s--){
	if(Delete[s]){
	   cycle++;
	   Int4 end=s; 
	   while(Delete[s] && s > 0) s--;   s++;
	   // fprintf(stderr,"end=%d; s=%d\n",end,s);
	   assert(end >= s); 
	   // routines from cma_gmb.cc
	   cma_typ rcma=0;
	   if(end==Len){ rcma=TrimBlkCMSA(cma,1,0,(end-s+1), 1); }
	   else if(s==1){ rcma=TrimBlkCMSA(cma,1,(end-s+1),0, 1); }
	   else { rcma=ConvertColsToInsertsCMSA(cma,1,s,end); }
	   if(rcma){ NilCMSA(cma); cma=rcma; }
#if 0
	   fprintf(stderr,"\n********************** cycle %d **************************\n",cycle);
	   char str[200]; sprintf(str,"%s.cycle%d",NameCMSA(cma),cycle);
           // FILE *fp = open_file(str,".cma","w"); PutCMSA(fp,cma); fclose(fp);
#endif
	}
    } free(Delete);
    return total;
}

void	PutIndelBrkPtsCMSA(FILE *rfp,cma_typ cma)
{
	assert(nBlksCMSA(cma) == 1);
	Int4 i,j,N=NumSeqsCMSA(cma);
        set_typ Set=IndelBreakPntsCMSA(1,1,cma,N/10);        // > 10 % seqs.
        fprintf(rfp,"  %d blocks [",CardSet(Set));
        for(i=1,j=0; i <= LengthCMSA(1,cma);i++){
                j++;
                if(MemberSet(i,Set)){
                  Int4 ins_res,n_ins=NumInsertsCMSA(1,i,cma,ins_res);
                  fprintf(rfp," %d (%.1f)",j,(double)ins_res/(double)N); j=0;
                }
        } fprintf(rfp," %d]\n",j); NilSet(Set); 
}

cma_typ	*GetCoreCMSA(cma_typ cmaI, cma_typ cmaJ,double &ave_frq,double frq_cut,FILE *fp)
{
	int	i,j,k,n,NumADV;
	char	str[200];
	a_type	AB = AlphabetCMSA(cmaI);
	assert(nBlksCMSA(cmaI) == 1); assert(nBlksCMSA(cmaJ) == 1); 
	Int4	lenI=LengthCMSA(1,cmaI),lenJ=LengthCMSA(1,cmaJ);
	Int4	nI=NumSeqsCMSA(cmaI),nJ=NumSeqsCMSA(cmaJ); 
	if(nI != nJ) print_error("GetCoreCMSA() different numbers of sequences in the two cma_typ.");
	double	*dd=0,sum; // NEW(dd,lenI+4,double); this is allocated within CommonColsCMSA().
	double	*nDelI;
	Int4	*ColI2J=CommonColsCMSA(cmaI,cmaJ,dd,nDelI,frq_cut);
	set_typ SetI=MakeSet(lenI+4); ClearSet(SetI);
	set_typ SetJ=MakeSet(lenJ+4); ClearSet(SetJ);
	e_type	csq1=MkConsensusCMSA(cmaI),csq2=MkConsensusCMSA(cmaJ);
	for(sum=0.0,n=0,i=1; i <= lenI; i++){
		j=ColI2J[i]; 
		if(j > 0){
#if 1	// need to resolve these issues!!!
		   if(!MemberSet(i,SetI) && !MemberSet(j,SetJ)) 
			{ AddSet(i,SetI); AddSet(j,SetJ); n++; }
		   else fprintf(stderr,"WARNING: ambiguous columns (%d vs %d)\n",i,j);
#else
			{ AddSet(i,SetI); AddSet(j,SetJ); n++; }
#endif
#if 0
			unsigned char rI=ResSeq(i,csq1),rJ=ResSeq(j,csq2); 
			fprintf(stderr,"%d: %c%d --> %c%d (%.3f)\n",
				n,AlphaChar(rI,AB),i,AlphaChar(rJ,AB),ColI2J[i],dd[i]); 
#endif
			if(fp) fprintf(fp,"%d\t%.4f\t%.4f\n",
					n,dd[i],(double)nDelI[i]/(double)nI); 
			sum+=dd[i];
		}
	} free(dd); free(ColI2J); NilSeq(csq1); NilSeq(csq2); free(nDelI);
	if(n > 0) ave_frq=sum/(double)n; else ave_frq=0.0;
	cma_typ	*rcma=0;
	i=CardSet(SetI); j=CardSet(SetJ); 
	if(i != j){ fprintf(stderr,"# columns: align A = %d; align B = %d\n",i,j); }
	if(i > 0){
	  char    *operA; NEW(operA,lenI+5,char); operA[0]='E'; operA[lenI+1]='E';
	  char    *operB; NEW(operB,lenJ+5,char); operB[0]='E'; operB[lenJ+1]='E';
          for(i=1; i <= lenI; i++){ if(!MemberSet(i,SetI)) operA[i]='d'; else operA[i]='m'; }
          for(i=1; i <= lenJ; i++){ if(!MemberSet(i,SetJ)) operB[i]='d'; else operB[i]='m'; }
#if 0
	  fprintf(stderr,"operA=%s\n",operA); fprintf(stderr,"operB=%s\n",operB);
#endif
	  cma_typ cmaA=RemoveTheseColumnsCMSA(operA,cmaI); 
	  if(cmaA == 0) cmaA=CopyCMSA(cmaI);
	  cma_typ cmaB=RemoveTheseColumnsCMSA(operB,cmaJ); 
	  if(cmaB == 0) cmaB=CopyCMSA(cmaJ);
	  rcma; NEW(rcma,4,cma_typ); rcma[1]=cmaA; rcma[2]=cmaB;
	  Int4    lenA=LengthCMSA(1,cmaA),lenB=LengthCMSA(1,cmaB);
	  assert(lenA == lenB);
	  free(operA); free(operB);
	} NilSet(SetI); NilSet(SetJ);
	return rcma;
}

cma_typ	MvColumnsCMSA(ssx_typ *ssx,set_typ &SetSq)
// Move columns from one side of insertions to the other if this increases the BILD score.
{
	cma_typ	xcma,rcma=0,cma=ssx->RtnCMA(); assert(nBlksCMSA(cma) == 1); 
	Int4	x,i,j,k,n,blk=1,lenX=LengthCMSA(blk,cma),N=NumSeqsCMSA(cma);
	double **BS=BILD_ScoresCMSA(ssx);
	a_type	AB=AlphabetCMSA(cma);
	FILE	*efp=0; // efp=stderr;

	set_typ	SetR=MakeSet(lenX+3); ClearSet(SetR);
	set_typ	SetA=MakeSet(lenX+3); ClearSet(SetA);
	dh_type dH=dheap(lenX+2,4);
	// h_type HG=Histogram("Column contributions to Contextual BILD scores(in nats)",-10,10,0.10);
	Int4	*nINS,*resINS,ins_res,nins,*nDEL,ndel;
	NEW(nINS,lenX+5,Int4); NEW(resINS,lenX+5,Int4); NEW(nDEL,lenX+5,Int4); 
	for(i=1; i <= lenX; i++){
		nINS[i]=NumInsertsCMSA(1,i,cma,ins_res); resINS[i]=ins_res;
	}
	set_typ	*SetS; NEW(SetS,lenX+5,set_typ); 
	enum location { left=1, center=2, right=3 };
	double **BILD; NEWP(BILD,5,double);
	NEW(BILD[center],lenX+5,double); 
	NEW(BILD[left],lenX+5,double); NEW(BILD[right],lenX+5,double);
	assert(lenX > 2);
	for(i=1; i <= lenX; i++){
		// if(!(nINS[i] > 0 || nINS[i-1] > 0)) continue;
		if(nINS[i] > 0) SetS[i]=MakeSet(N+5); 
		double	d,frct=0.750;
		UInt4 WtCnts[5][30],sqwt;
		unsigned char	r,del=AlphaCode('X',AB);
		for(j=0; j <= nAlpha(AB); j++) WtCnts[left][j]=WtCnts[center][j]=WtCnts[right][j]=0;
		Int4	sq,pos,Pos[5]; 	
		if(SetS[i]) ClearSet(SetS[i]);
		double	TotalWtSq=ssx->TotalWtSeq();
		for(ndel=0,sq=1; sq <= N; sq++){
		   unsigned char   R,RR,*isq=SeqPtrCMSA(sq,cma);	// FakeSeq ptr.
		   e_type tE=TrueSeqCMSA(sq,cma); sqwt=ssx->RtnSeqWt(sq);

		   if(SetS[i]){ if(InsertionCMSA(1,sq,i,cma) > 0) AddSet(sq,SetS[i]); }
		   pos=TruePosCMSA(sq,i,cma);

		   // confirm that true and fake residues correspond.
		   PosSiteCMSA(1,sq,Pos,cma); j = Pos[1]+i-1; R=isq[j];  
                   // if(IsDeletedCMSA(sq,j,cma)){ RR=0; ndel++; } else { RR=ResSeq(pos,tE); } 
                   if(IsDeletedCMSA(1,sq,i,cma)){ RR=0; ndel++; } else { RR=ResSeq(pos,tE); } 
		   assert(R == RR);
	// gsq_typ *gsq=gsqCMSA(sq,cma); gsq->Put(stderr,AB); exit(1);
// sqwt=1;
		   WtCnts[center][R] += sqwt;
		   // fprintf(stderr,"%d.%d: %c%d == %c%d?\n",i,sq,AlphaChar(R,AB),j,AlphaChar(RR,AB),pos);
		   // pos=FakeToRealCMA(sq,i,cma);
		   if(pos == 1){			// start of sequence. "A--","A-C","--F"
                        if(i < lenX && IsDeletedCMSA(1,sq,i+1,cma)){ r=del; } 
			// else if(IsDeletedCMSA(1,sq,i,cma)){ r=del; }	// "?--" 2 deletions= 
			else { r=ResSeq(pos+1,tE); } 
			WtCnts[right][r] += sqwt; WtCnts[left][del] += sqwt;
		   } else if(pos == LenSeq(tE)){	// end of sequence.
                        if(IsDeletedCMSA(1,sq,i,cma)){			// "?-?"
                            if(i > 1 && IsDeletedCMSA(1,sq,i-1,cma)){ r=del; }	// "--?"
			    else r=ResSeq(pos,tE);  				// "X-?"; at pos is X;
                        } else if(i > 1 && IsDeletedCMSA(1,sq,i-1,cma)) r=del;	// "-X?" 
			else { r=ResSeq(pos-1,tE); } 			// "XY?" ; at pos-1 is X.
			WtCnts[left][r] += sqwt; WtCnts[right][del] += sqwt;
		   } else {
                        if(IsDeletedCMSA(1,sq,i,cma)){			// "?-?"
                            if(i > 1 && IsDeletedCMSA(1,sq,i-1,cma)){ r=del; }	// "--?"
			    else r=ResSeq(pos,tE);  				// "X-?"; at pos is X;
                        } else if(i > 1 && IsDeletedCMSA(1,sq,i-1,cma)) r=del;	// "-X?" 
			else { r=ResSeq(pos-1,tE); } 			// "XY?" ; at pos-1 is X.
			WtCnts[left][r] += sqwt;

                        if(i < lenX && IsDeletedCMSA(1,sq,i+1,cma)){ r=del; }
			// else if(IsDeletedCMSA(1,sq,i,cma)){ r=del; }	// "--" two deletions...
			else { r=ResSeq(pos+1,tE); } 
			WtCnts[right][r] += sqwt;
		   }
		}
// TotalWtSq=(double)N/1000;;
		if(1 || nINS[i] > 0) BILD[right][i]=ssx->BildScore(WtCnts[right])/TotalWtSq;
		if(1 || nINS[i-1] > 0) BILD[left][i]=ssx->BildScore(WtCnts[left])/TotalWtSq;
		// BILD[center][i]=ssx->BildScore(1,i)/TotalWtSq; 
		BILD[center][i]=ssx->BildScore(WtCnts[center])/TotalWtSq; 
#if 0
		fprintf(stderr,"%3d: left ",i); ssx->PutWtCnts(stderr,WtCnts[left]);
		// fprintf(stderr,"%3d:      ",i); ssx->PutWtCnts(stderr,1,i);
		fprintf(stderr,"%3d:      ",i); ssx->PutWtCnts(stderr,WtCnts[center]);
		fprintf(stderr,"%3d: right",i); ssx->PutWtCnts(stderr,WtCnts[right]);
#endif
		nDEL[i]=ndel;
	}
	double	last_bsL=0,bs,bsC,bsL,bsR,bld_lt,bld_rt;
	char	*operation=0; NEW(operation,lenX+5,char); operation[0]='E';
	BooLean	found=FALSE;
	for(i=1; i <= lenX; i++) operation[i]='.'; operation[i]='E';
	SetSq=MakeSet(N+5); ClearSet(SetSq);
	for(last_bsL=0,i=1; i <= lenX; i++){
	    // if(!(nINS[i] > 0 || nINS[i-1] > 0)) continue; 
	    // else an insertion before and/or after column i.
	    double	bild_cut=0.002,dd,DD,TotalWtSq=ssx->TotalWtSeq();
	    char	move='.';
	    bs=0.0; bsL=BILD[center][i-1]; bsR=BILD[center][i+1]; bsC=BILD[center][i]; 
	    bld_lt=BILD[left][i]; 		// to the left of & up against column i.
	    bld_rt=BILD[right][i-1];		// ro the right of & up agaisnt column i-1.
	    if(nINS[i-1] > 0){			// Insertion before column i: 
		dd=bld_lt - bsL; 
		if(dd > bild_cut){				// ...[i-1]xx(l)[i]...
		   move='R'; bs=bld_lt;			// move column i-1 right, to the left of i.
		   if(efp) fprintf(efp,"Move column %d (%.3f) right (to left of & up against column %d) (%.3f).\n",
				i-1,bsL,i,bs);
		} 
		dd=bld_rt - bsC; DD=bld_rt-bs; 
		if(dd > bild_cut && DD > bild_cut){		// ...[i-1](r)xx[i]...
		   move='L'; bs=bld_rt;			// move column i left, to right of i-1.
		   if(efp) fprintf(efp,"Move column %d (%.3f) left (to right of & up against column %d) (%.3f).\n",
			i,bsC,i-1,bs);
		}
		if(bs > 0){
		   if(move=='R'){
		      assert(i <= lenX);	
		      // core dumping here with i == lenX + 1 and next assertion failing; dbx says op[i-1]='.'!!
		      if(operation[i-1] != '.'){   // then need to decide to move either 'R' or 'L'.
		        assert(operation[i-1]=='L');
			if(last_bsL < bs){
			  operation[i-1]='R'; found=TRUE; 
			  assert(SetS[i-1]); UnionSet(SetSq,SetS[i-1]);
			} // FILE *ofp=open_file("pre_core_dump",".cma","w"); PutCMSA(ofp,cma); fclose(ofp);
		      } else {
		        operation[i-1]='R'; found=TRUE; assert(SetS[i-1]); UnionSet(SetSq,SetS[i-1]);
		      }
		   } else if(move=='L'){
			assert(i <= lenX); operation[i]='L'; found=TRUE; 
			assert(SetS[i-1]); UnionSet(SetSq,SetS[i-1]);
			last_bsL=bs;
		   }
		} else last_bsL=0.0;
	    } else last_bsL=0.0;
	    if(efp){ 
		 fprintf(efp,"%d: BILD = %.3f; bld_lt=%.3f; bld_rt=%.3f; deletions=%d (%.3f)\n",
			i,bsC,BILD[left][i],BILD[right][i],nDEL[i],(double) nDEL[i]/(double)N);
	         if(nINS[i] > 0) fprintf(efp,"   %d: %d insertions\n",i,nINS[i]); 
	    }
	    // IncdHist(d,HG);
	    // insrtHeap(i,(keytyp)d,dH);  // focus on BILD scores.
	} // PutHist(stderr,50,HG); NilHist(HG); 
	// fprintf(stderr,"operation = %s\n",operation);
	if(efp) fprintf(efp,"# seqs realigned = %d (%d)\n",CardSet(SetSq),N); // NilSet(SetSq);
	if(found) rcma=MoveColumnsCMSA(cma,operation); else rcma=0;
	Nildheap(dH); free(operation); 
	for(i=1; i <= lenX; i++){ if(SetS[i]) NilSet(SetS[i]); } free(SetS);
	free(BILD[center]); free(BILD[left]); free(BILD[right]); free(BS[1]); free(BS); 
	free(BILD); NilSet(SetA); NilSet(SetR); free(nINS); free(resINS); free(nDEL);
	return rcma;
}

char    *gsq_typ::MvColsOperation(char *MvOper, char *operation)
// ========== Move columns within an operational array. ===========
// NOTE: only need to worry about insertions within and at end of domains ('I' and 'i'). 
{
	Int4	mv_len=strlen(MvOper);
        Int4	o,i,j,no,column,blk,x,y,trace_length=strlen(operation);
        char	token,Token,*new_operation;
		NEW(new_operation,trace_length + 5,char); new_operation[0]='E';
        for(no=o=1,j=column=blk=0; operation[o] != 'E'; o++){
          token=operation[o]; assert(j < mv_len); 
          switch(token){
            case 'M': blk++; column=0;
            case 'm': new_operation[no]=token; no++; column++; j++; break;
            case 'D': blk++; column=0;
            case 'd': new_operation[no]=token; no++; column++; j++; break;
            case 'i': // insert is between profile blocks;
#if 0
		if(MvOper[j] == 'R'){	// Move firt column to the right: iii?IIX --> iii?XII
		} else if(MvOper[j+1] == 'L'){	// move next column to the left: mIII? --> m?III
#endif
		{ new_operation[no]=token; no++; }
                break;
            case 'I': // Insert ('-') within an aligned block.
		if(MvOper[j] == 'R'){	// Move last column to the right:      ?IIIX --> III?X
			no--; token=new_operation[no];
			if(column == 1) Token='i'; else Token='I';
			for( ; operation[o]=='I'; o++){
			     new_operation[no]=Token; no++;  // added same number...but shifted to left.
			} o--; new_operation[no]=token; no++;  	// last token; no now in sync with o.
		} else if(MvOper[j+1] == 'L'){	// move next column to the left: mIII? --> m?III
			if(MvOper[j+2]=='E') Token='i'; else Token='I';
			for(x=no ; operation[o]=='I'; o++){
				new_operation[no]=Token; no++;
			} new_operation[no] = Token; no++; 
			new_operation[x] = operation[o]; column++; j++;
        // fprintf(stderr,"new operation: %s\n",new_operation);
		} else { new_operation[no]=token; no++; } break;
            default:
            // fprintf(stderr,"operations = %s\n",operation);
            print_error("operation( ): input error"); break;
          }
        } new_operation[no]='E'; no++; new_operation[no]=0;
#if 0
	PutShortSeqID(stderr,realE); fprintf(stderr," operation: %s\n",operation);
        fprintf(stderr,"new operation: %s\n",new_operation);
#endif
	assert(no == trace_length);
        return new_operation;
}

gsq_typ *gsq_typ::MoveColumns(char *MvOper, Int4 nBlks, Int4 *len, Int4 *sites)
{
        Int4    numfix;
        char *new_operation,*operation=this->Operation(nBlks,sites,len);
        new_operation=MvColsOperation(MvOper,operation); free(operation);
        // do { numfix=this->RmInsertByDeleteOperation(new_operation); } while(numfix > 0);
        // fprintf(stderr,"new operation: %s\n",new_operation);
        gsq_typ *gsq0 = new gsq_typ[1];
        gsq0->initialize(new_operation,realE,sites);    // copies positions to sites array.
        free(new_operation);
        return gsq0;
}

cma_typ	MoveColumnsCMSA(cma_typ cma, char *operation)
// Insert new columns after start_ins (lower to upper case in cma file).
//      operation = E.R.......................R............L.............................E
// --> move columns   ^ to right              ^to right    ^ to left.
//     implies that there is an insertion before and after to-left and to-right columns, respectively.
{
	a_type  AB=AlphabetCMSA(cma); assert(nBlksCMSA(cma) == 1);
	Int4	blk=1,*len,numfix;
	gss_typ *gss=gssCMSA(cma); NEW(len,4, Int4); len[blk]=LengthCMSA(blk,cma);
   	cma_typ rcma=EmptyCMSA(nBlksCMSA(cma),len,TrueDataCMSA(cma),gss->GapOpen(),
			gss->GapExtend(),PerNatsCMSA(cma),0,0);
	for(Int4 sq=1; sq <= NumSeqsCMSA(cma); sq++){
	        // ========== X. Get operational array for sequence. ===========
		gsq_typ *gsq=gsqCMSA(sq,cma);
        	if(!SeqIsInCMSA(sq,cma)){    // subcma operations; check for vacant site ...afn: 11-26-2014.
			if(gsq==0){
			   gsq_typ *gsq0 = new gsq_typ[1];
                	   gsq0->initialize(TrueSeqCMSA(sq,cma));    // copies positions to sites array.
			   ReplaceCMSA(sq,gsq0,rcma);
			} continue;
		}
		Int4	*sites=GetPosSitesCMSA(sq,cma);
		gsq_typ	*gsq0=gsq->MoveColumns(operation,nBlksCMSA(cma),LengthsCMSA(cma),sites);
		ReplaceCMSA(sq,gsq0,rcma); // replace sequence s in CMSA & fmodel.
		AddSiteCMSA(blk,sq,sites[blk],rcma); free(sites); 
	} free(len);
	return rcma;
}

set_typ **ParticipatingColPairsCMSA(cma_typ cmaA, cma_typ cmaB, double frq_cut)
// return the set of columns in cmaA that also occur in cmaB.
{
	//=============== 1. Check alignment consistency =========================
	Int4	i,j,k,sq,N=NumSeqsCMSA(cmaA),lenA=LengthCMSA(1,cmaA),lenB=LengthCMSA(1,cmaB);
	assert(nBlksCMSA(cmaA) == 1); assert(nBlksCMSA(cmaB) == 1);
	if(NumSeqsCMSA(cmaA) !=  NumSeqsCMSA(cmaB)){
		fprintf(stderr,"NumSeqs: %d != %d\n",NumSeqsCMSA(cmaA),NumSeqsCMSA(cmaB));
		print_error("ParticipatingColPairsCMSA(): input alignments are inconsistent");
	} 
        if(lenA != lenB){
		print_error("ParticipatingColPairsCMSA(): input alignment lengths are inconsistent");
	}
	for(sq=1; sq <= N; sq++){
	    e_type sqA=TrueSeqCMSA(sq,cmaA),sqB=TrueSeqCMSA(sq,cmaB);
	    if(!FastIdentSeqs(sqA,sqB)){
		a_type	AB=AlphabetCMSA(cmaA);
		AlnSeqSW(stderr,11,1,sqA,sqB,AB);
		fprintf(stderr,"Sequences in row %d differ.\n",sq);
		print_error("ParticipatingColPairsCMSA(): input sequences are inconsistent");
	    }
	}

	//=============== 2. Check alignment consistency =========================
	Int4	n,x,xA,xB,yA,yB,y,z,M=MaxTrueSeqCMSA(cmaA),max_j;
	double	d,dd,max_d;
	set_typ	**RtnSet=0; NEWP(RtnSet, lenA+5,set_typ);
	for(i=1; i <= lenA; i++){
	   NEW(RtnSet[i], lenA+5,set_typ);
	   for(j=1; j <= lenA; j++){
	   	RtnSet[i][j]=MakeSet(N+5); ClearSet(RtnSet[i][j]);
	   }
	}
	for(sq=1; sq <= N; sq++){
	    for(i=1; i <=lenA; i++){
                if(IsDeletedCMSA(1,sq,i,cmaA) || IsDeletedCMSA(1,sq,i,cmaB)) continue;
	        xA=TruePosCMSA(sq,1,i,cmaA); xB=TruePosCMSA(sq,1,i,cmaB);
		// if(xA != xB) continue;
	    	for(j=i+1; j <= lenB; j++){
                   if(IsDeletedCMSA(1,sq,j,cmaA) || IsDeletedCMSA(1,sq,j,cmaB)) continue;
	           yA=TruePosCMSA(sq,1,j,cmaA); yB=TruePosCMSA(sq,1,j,cmaB);
		   // if(yA != yB) continue;
if(abs(xA - yA) < 8) continue;
if(abs(xB - yB) < 8) continue;
		   AddSet(sq,RtnSet[i][j]); AddSet(sq,RtnSet[j][i]);
		}
	    }
	}
	for(i=1; i <= lenA; i++){
	   for(j=i+1; j <= lenA; j++){
		x=CardSet(RtnSet[i][j]); y=CardSet(RtnSet[j][i]);
		d=(double)x/(double)N;
		if(x != y || d < frq_cut){
			// fprintf(stderr,"RtnSet[%d][%d]=%d (%3f)\n",i,j,x,d);
			// print_error("ParticipatingColPairsCMSA(): frequency of matches too low.");
		}
	   }
	} return RtnSet;
}

set_typ *IdenticallyAlignedCMSA(cma_typ cmaI, cma_typ cmaJ)
// return the set of columns in cmaI that also occur in cmaJ.
{
	//=============== 1. Check alignment consistency =========================
	if(NumSeqsCMSA(cmaI) !=  NumSeqsCMSA(cmaJ)){
		fprintf(stderr,"NumSeqs: %d != %d\n",NumSeqsCMSA(cmaI),NumSeqsCMSA(cmaJ));
		print_error("IdenticallyAlignedCMSA(): input alignments are inconsistent");
	} 
	assert(nBlksCMSA(cmaI) == 1); assert(nBlksCMSA(cmaJ) == 1);
	Int4	xI,xJ,x,i,j,k,sq,N=NumSeqsCMSA(cmaI),lenI=LengthCMSA(1,cmaI),lenJ=LengthCMSA(1,cmaJ);
	for(sq=1; sq <= N; sq++){
	    e_type sqI=TrueSeqCMSA(sq,cmaI),sqJ=TrueSeqCMSA(sq,cmaJ);
	    if(!FastIdentSeqs(sqI,sqJ)){
		a_type	AB=AlphabetCMSA(cmaI);
		AlnSeqSW(stderr,11,1,sqI,sqJ,AB);
		fprintf(stderr,"Sequences in row %d differ.\n",sq);
		print_error("IdenticallyAlignedCMSA(): input sequences are inconsistent");
	    }
	}

	//=============== 2. Find common columns =========================
        double  *DD,*Del,frq_cut=0.50,ave_frq=0,sum;
        Int4    *ColI2J=CommonColsCMSA(cmaI,cmaJ,DD,Del,frq_cut);

	//=============== 3. find sequences identical at each shared column =========================
	set_typ	*RtnSet=0; NEW(RtnSet,lenI+5,set_typ);
        for(i=1; i <= lenI; i++){
           if((j=ColI2J[i]) > 0){
	     RtnSet[i]=MakeSet(N+5); ClearSet(RtnSet[i]); 
	     for(sq=1; sq <= N; sq++){
                   if(IsDeletedCMSA(1,sq,i,cmaI) || IsDeletedCMSA(1,sq,j,cmaJ)) continue;
	           xI=TruePosCMSA(sq,1,i,cmaI); xJ=TruePosCMSA(sq,1,j,cmaJ); 
		   if(xI == xJ) AddSet(sq,RtnSet[i]); 
	     }
	   }
        } free(DD); free(Del); free(ColI2J); return RtnSet;
}

double	ComputSeqWtsCMSA(cma_typ CMA,set_typ Set)
// return an integer weight for each sequence in the cmsa.
{
        Int4    i,j,k,r,sq,b,n,totlen=TotalLenCMSA(CMA),N=NumSeqsCMSA(CMA);
        double  w,max,N2,d,*wt;
        unsigned char    *seq;
	UInt4	**nres,*ntyp;
	a_type	AB=AlphabetCMSA(CMA);

	NEW(wt,N+2,double); 
        NEWP(nres,totlen+2,UInt4); NEW(ntyp,totlen+2,UInt4);
        for(i=0; i<=totlen; i++) { NEW(nres[i],nAlpha(AB)+2,UInt4); }
        /*** 1. determine the number of residues at each position ***/
        for(sq=1; sq<=N; sq++){
	   if(Set && !MemberSet(sq,Set)) continue;
           for(j=1,b=1; b <= nBlksCMSA(CMA); b++){
		seq=GetAlnResInSiteCMSA(b,sq,CMA);
		if(seq == 0) continue;
                for(i=1; i<=LengthCMSA(b,CMA); i++){
                    r=seq[i]; 
		    if(r==0) continue;
		    if(nres[j][r] == 0) { ntyp[j]++; } 
		    nres[j][r]++; j++;
                }
           }
        }

        /*** 2. determine the sequence weights. ***/
        for(max=0.,sq=1; sq<=N; sq++){
	   if(Set && !MemberSet(sq,Set)) continue;
           for(w=0.0,j=1,b=1; b <= nBlksCMSA(CMA); b++){
		seq=GetAlnResInSiteCMSA(b,sq,CMA);
		if(seq == 0) continue;
                for(i=1; i<=LengthCMSA(b,CMA); i++){
                    r=seq[i]; 
		    if(r==0) continue;
                    w+= 1.0/ (double) (ntyp[j]*nres[j][r]); j++; 
                }
           } w /= (double) totlen; 
	   wt[sq]=w; max = MAXIMUM(double,w,max); 
        }
	h_type	HG=0;
#if 0	// For debugging code...
	HG = Histogram("number residue types", 0, 30,1.0);
        for(j=1; j<=totlen; j++){ IncdHist(ntyp[j],HG); }
	PutHist(stderr,60,HG);  NilHist(HG);

        /*** 3. normalize the weights. ***/
	HG = Histogram("sequence weights", 0, 100,2); 
	// HG = Histogram("sequence weights", 0, 100,1.0); 
#endif
        for(N2=0.0,sq=1; sq<=N; sq++){
	   if(Set && !MemberSet(sq,Set)) continue;
	   wt[sq] /= max; 
	   // if(SqWtAdjst < 1.0) wt[sq] = pow(wt[sq],SqWtAdjst);
	   N2 += wt[sq];
	   if(HG) IncdHist(100*wt[sq],HG);
        }
	if(HG){ PutHist(stderr,60,HG);  NilHist(HG); }
        for(i=0; i<=totlen; i++) { free(nres[i]); }
        free(nres); free(ntyp); free(wt); 
	return N2;
}

