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

char	*gsq_typ::RmFlankingInOperation(char *operation)
// Rmove sequence regions beyond the ends of the alignment (as for cma2fa routine).
{
	char *new_operation=0;
	Int4 blk,o,no,trace_length=strlen(operation);
        NEW(new_operation,trace_length+5,char);
	new_operation[0]='E'; no=1;
        for(o=1; operation[o] == 'i'; o++) ;
        for( ; operation[o] != 'E'; o++){ new_operation[no]=operation[o]; no++; }
        for(no-- ; new_operation[no] == 'i'; no--) ;
	no++; new_operation[no]='E'; no++; new_operation[no]=0;
        return new_operation;
}

gsq_typ *gsq_typ::RmOverHangs(Int4 nblk,Int4 *sites, Int4 *len,Int4 start,Int4 end)
// Rmove sequence regions beyond the ends of the alignment (as for cma2fa routine).
// this is needed for mapgaps.
{
	// 1. find the ends of the sequence.
	char *operation=this->Operation(nblk,sites,len);
	// fprintf(stderr,"operation=%s\n",operation);
	char	*nop=this->RmFlankingInOperation(operation);
	// fprintf(stderr," new operation=%s\n",nop);
	e_type  subE=MkSubSeq(start,end,realE);
        gsq_typ *gsq = new gsq_typ[1];
	gsq->initialize(nop,subE, sites);
	free(operation); free(nop);
	return gsq;
}

gsq_typ *gsq_typ::RmBlock(Int4 Blk,Int4 nBlks, Int4 *len, Int4 *sites)
{
        char *operation=this->Operation(nBlks,sites,len);
        fprintf(stderr,"operation=%s\n",operation);
        char *new_operation=this->RmBlkInOperation(Blk,operation);
        fprintf(stderr," new operation=%s\n",new_operation);
        Int4 numfix,numiters=0;
        do { numfix=this->RmInsertByDeleteOperation(new_operation); numiters++; } while(numfix > 0);
        if(numiters > 1) fprintf(stderr," new new operation=%s\n",new_operation);
        gsq_typ *gsq0 = new gsq_typ[1];
#if 0   // this is not needed...new sites copied to site in gsq->initialize();
        for(b=0,blk=1; blk <= nBlksCMSA(cma); blk++){
                if(blk == Blk) continue; b++; sites[b]=sites[blk];
        } sites[b+1]=0;
#endif
        gsq0->initialize(new_operation,realE,sites);    // copies positions to sites.
        free(new_operation);
        return gsq0;
}

Int4    gsq_typ::RmInsertByDeleteOperation(char *operation)
/*************************************************************************************
 Merge ..IIIdd.. and ..ddIII... to fix program...so can work with cma_gblastpgp...
	Actually need ...IIIddm... or ...mddIII...
 *************************************************************************************/
{
        Int4 length=strlen(operation);
        char state,last=' ';
        Int4 i,j,k,strt_j,strt_k,net,x;
        Int4 num_i=0,num_d=0,num_fix=0;

        assert(operation[0] == 'E');
        if(operation[1] == 'D'){ for(i=2; i < length; i++) if(operation[i] != 'd') break; }
        else i=1;
        last=operation[i-1];
        for( ; i < length; i++){
                state = operation[i];
                if(state=='E') break;
                if(state != last){      // then
                   BooLean found=FALSE;
                   if(state=='I' && last == 'd'){               // found ...dI...
                        for(num_i=0,j=i; operation[j] == 'I'; j++) num_i++;
                        for(num_d=0,j=i-1; operation[j] == 'd'; j--) num_d++;
#if 1	// if there are no matches before the deletion, then do nothing...
			if(operation[j] != 'M' && operation[j] != 'm'){
			   for(x=j; operation[x]=='I' || operation[x]=='i'; x--) ; // ..mIIddd.. is okay.
			   if(operation[x] != 'M' && operation[x] != 'm') continue;
			}
#endif
                   } else if(state=='i' && last == 'd'){                // found ...dI...
                        for(num_i=0,j=i; operation[j] == 'i'; j++) num_i++;
                        for(num_d=0,j=i-1; operation[j] == 'd'; j--) num_d++;
#if 1	// if there are no matches before the deletion, then do nothing...
			if(operation[j] != 'M' && operation[j] != 'm'){
			   for(x=j; operation[x]=='I' || operation[x]=='i'; x--) ; // ..mIIddd.. is okay.
			   if(operation[x] != 'M' && operation[x] != 'm') continue;
			}
#endif
                   } else if(state == 'd' && last == 'I'){      // found ...Id...
                        for(num_d=0,j=i; operation[j] == 'd'; j++) num_d++;
#if 1	// if there are no matches after the deletion, then do nothing...
			if(operation[j] != 'M' && operation[j] != 'm'){
			   for(x=j; operation[x]=='I' || operation[x]=='i'; x++) ; // ..dddIIm.. is okay.
			   if(operation[x] != 'M' && operation[x] != 'm') continue;
			}
#endif
                        for(num_i=0,j=i-1; operation[j] == 'I'; j--) num_i++;
                   } else { last=state; continue; }
                   if(num_i <= num_d){          // ..RGvs---KA.. or ..GKstl---LR..
                                strt_j = i - num_i; strt_k = i + num_i; net = num_i;
                   } else if(num_i > num_d){    // ..NVafp--LK..
                                strt_j = i - num_d; strt_k = i + num_d; net = num_d;
                   } // change 'I' to 'm'.
                   for(j=strt_j,k=0; k < net; j++,k++){ operation[j]='m'; }
                   for(k=strt_k; operation[k] != 'E'; k++,j++){
                                operation[j]=operation[k];
                   } operation[j]='E'; operation[j+1]=0;
                   i=strt_j;
                   state=operation[i]; num_fix++;
                   // fprintf(stderr,"state=%c; i = %d\n",state,i);
                   // if(state== 'I') print_error("RmInsertByDeleteOperation( ) input error");
                } last=state;
        }
        return num_fix;
}

gsq_typ *gsq_typ::IronOut(Int4 nblk,Int4 *start, Int4 *blk_len)
// for multi-block alignments.
{
        Int4 n;
        char *Op=this->Operation(nblk,start,blk_len);
        // fprintf(stderr,"%s\n",Op);
        // gsq->Put_cma_format(stderr,0,nBlksCMSA(cma),sites,LengthsCMSA(cma),AB);
        // fprintf(stderr,"Operation = %s\n",Op);
        do { n=this->RmInsertByDeleteOperation(Op); } while(n > 0);
        // fprintf(stderr,"Operation = %s\n",Op);
        gsq_typ *gsq0 = new gsq_typ[1];
        gsq0->initialize(Op,this->TrueSeq(),start);
        free(Op);
        return gsq0;
}

gsq_typ *gsq_typ::IronOut(Int4 nblk,Int4 *start, Int4 *blk_len,Int4 &X)
{
        Int4 n,i,pos[9];  pos[2]=0;
        char *Op=this->Operation(nblk,start,blk_len);
        // fprintf(stderr,"%s\n",Op);
        // gsq->Put_cma_format(stderr,i,nBlksCMSA(cma),sites,LengthsCMSA(cma),AB);
        // fprintf(stderr,"Operation = %s\n",Op);
	do { n=this->RmInsertByDeleteOperation(Op); } while(n > 0);
        // fprintf(stderr,"Operation = %s\n",Op);
        gsq_typ *gsq0 = new gsq_typ[1];
        gsq0->initialize(Op,this->TrueSeq(),pos);
        X = pos[1]; free(Op);
        return gsq0;
}

void	gsq_typ::ConvertToOneBlock(char *Op)
{
	Int4	end,m=0,i,j,len=strlen(Op);
	char	state;
	for(end=len-2; Op[end] == 'i'; end--) ;
	for(state='N',i=1; i <= end; i++){
	    if(state == 'N' && (Op[i] == 'M' || Op[i] == 'D')){
		state='A'; continue;
	    } else if(state == 'N'){ continue;
	    } else switch(Op[i]){
		case 'i': Op[i]='I'; break;
		case 'D': case 'M': Op[i]=tolower(Op[i]); break;
		case 'I': case 'm': case 'd': break; // do nothing.
		default: 
		  fprintf(stderr,"%s\n",Op);
		  fprintf(stderr,"error at position %d = %c\n",i,Op[i]);
		  print_error("ConvertToOneBlock() input error"); break;
	    }
	}
}

gsq_typ	*gsq_typ::ToOneBlk(Int4 nblk,Int4 *start, Int4 *blk_len,Int4 &X)
{
	Int4 i,pos[9];  pos[2]=0;
	char *Op=this->Operation(nblk,start,blk_len);
	// fprintf(stderr,"input: %s\n",Op);
	// gsq->Put_cma_format(stderr,i,nBlksCMSA(cma),sites,LengthsCMSA(cma),AB);
	this->ConvertToOneBlock(Op);
	// fprintf(stderr,"output: %s\n",Op);
        gsq_typ *gsq0 = new gsq_typ[1]; 
        gsq0->initialize(Op,this->TrueSeq(),pos);
	X = pos[1]; free(Op);
	return gsq0;
}

char	*gsq_typ::FuseBlksInOperation(Int4 Blk, char *operation)
// Fuse Blk with Blk + 1 converting gap residues to insertions...
{
	assert(Blk > 0); 
	char state,*new_operation=0;
	Int4 o,no,column,blk;
	Int4 trace_length=strlen(operation);
	NEW(new_operation,trace_length+5,char);
	new_operation[0]='E'; no=1;
	for(o=1,column=0,state='E',blk=0; operation[o] != 'E'; o++){
          switch(operation[o]){
            case 'M': 
		if(blk==Blk) new_operation[no]='m';  else new_operation[no]='M';
		blk++; no++; break;
            case 'm': new_operation[no]=operation[o]; no++; break;
            case 'D': 
		if(blk==Blk) new_operation[no]='d';  else new_operation[no]='D'; 
		blk++; no++; break;
            case 'd': new_operation[no]=operation[o]; no++; break;
            case 'i': // insert is between profile blocks;
		if(blk==Blk) new_operation[no]='I'; else new_operation[no]='i';  no++; break;
            case 'I': // Insert ('-') within a profile block; delete from seq.
		new_operation[no]=operation[o]; no++; break;
            default:
            // fprintf(stderr,"operations = %s\n",operation);
            print_error("operation( ): input error"); break;
          }  state=operation[o];
       	} assert(blk > Blk);
	new_operation[no]='E'; no++; new_operation[no]=0;
	return new_operation;
}

gsq_typ	*gsq_typ::FuseBlks(Int4 theBlk, Int4 nBlks, Int4 *len, Int4 *sites)
// gsq0=gsq->FuseBlks(nBlksCMSA(cma),sites,LengthsCMSA(cma));
{
	Int4	n;
	char *operation=this->Operation(nBlks,sites,len);
	// fprintf(stderr,"operation: %s\n",operation);
	// this->Put_cma_format(stderr,sq,nBlksCMSA(cma),sites,LengthsCMSA(cma),AB);
	char *new_operation=this->FuseBlksInOperation(theBlk,operation);
	// this->RmInsertByDeleteOperation(new_operation);
	do { n=this->RmInsertByDeleteOperation(new_operation); } while(n > 0);
	// fprintf(stderr,"new operation: %s\n",new_operation);
	gsq_typ *gsq0 = new gsq_typ[1];
	gsq0->initialize(new_operation,realE,sites);	// copies positions to sites array.
	free(new_operation); 
	return gsq0;
}

char	*gsq_typ::AddInsertToOperation(Int4 Blk, Int4 start, Int4 end, char *operation)
// ========== Add an insertion to an operational array. ===========
{
	assert(Blk > 0); 
	assert(start > 1 && start <= end);
	char *new_operation=0;
	Int4 o,no,column,blk;
	Int4 trace_length=strlen(operation);
	NEW(new_operation,trace_length+5,char);
	new_operation[0]='E'; no=1;
	for(o=1,column=0,blk=0; operation[o] != 'E'; o++){
          switch(operation[o]){
            case 'M': if(blk==Blk) assert(column-1 > end); blk++; column=1; 
            case 'm': 
		if(blk==Blk && column >= start && column <=end){
			   new_operation[no]='I'; 
		} else new_operation[no]=operation[o];
		no++; column++; break;
            case 'D': if(blk==Blk) assert(column-1 > end); blk++; column=1; 
            case 'd': // deletion in sequence relative to profile.
		if(blk==Blk && column >= start && column <=end){
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
          }  
       	} if(blk==Blk) assert(column-1 > end); 
	new_operation[no]='E'; no++; new_operation[no]=0;
	// printf("new operation: %s\n",new_operation);
	return new_operation;
}

gsq_typ	*gsq_typ::ConvertColsToInsert(Int4 theBlk, Int4 start, Int4 end, Int4 nBlks,
                        Int4 *len, Int4 *sites)
{
	Int4	numfix;
	char *operation=this->Operation(nBlks,sites,len);
	char *new_operation=AddInsertToOperation(theBlk,start,end,operation);
	do { numfix=this->RmInsertByDeleteOperation(new_operation); } while(numfix > 0);
	// this->RmInsertByDeleteOperation(new_operation);
	// printf("new operation(%d): %s\n",sq,new_operation);
	gsq_typ *gsq0 = new gsq_typ[1];
	gsq0->initialize(new_operation,realE,sites);	// copies positions to sites array.
	free(new_operation); free(operation);
	return gsq0;
}

#if 0
char	*gsq_typ::Operation2(Int4 nblk,Int4 *start, Int4 *blk_len)
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
		Operation[ptr]=c; ptr++;
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
#endif

gsq_typ *gsq_typ::IronOut2(Int4 nblk,Int4 *start, Int4 *blk_len,Int4 &X)
{
        Int4 n,i,pos[9];  pos[2]=0;
        char *Op=this->Operation(nblk,start,blk_len);
        // fprintf(stderr,"%s\n",Op);
        // gsq->Put_cma_format(stderr,i,nBlksCMSA(cma),sites,LengthsCMSA(cma),AB);
        // fprintf(stderr,"Operation = %s\n",Op);
        do { n=this->RmInsertByDeleteOperation(Op); } while(n > 0);
        // fprintf(stderr,"Operation = %s\n",Op);
        gsq_typ *gsq0 = new gsq_typ[1];
        gsq0->initialize(Op,this->TrueSeq(),start);
        X = pos[1]; free(Op);
        return gsq0;
}

void	gsq_typ::SplitBlkInOperation(Int4 *pos, char *operation)
{
        Int4	o,column,blk;
	assert(operation[0]=='E');
        for(o=1,column=0,blk=0; operation[o] != 'E'; o++){
          switch(operation[o]){
            case 'M': case 'D': 
		blk++; 
	    case 'm': case 'd':
		column++; 
		if(column == pos[blk]){ operation[o]=toupper(operation[o]); blk++; }
		break;
            case 'i': // insert is between profile blocks;
		break; // do nothing
            case 'I': // Insert ('-') within a profile block; delete from seq.
		if(column+1 == pos[blk]) operation[o]='i';
		break;
            default:
            // fprintf(stderr,"operations = %s\n",operation);
            print_error("SplitBlkInOperation( ): input error"); break;
          }  
        }  assert(pos[blk+1]==0);
	return;
}

char	*gsq_typ::TrimBlkInOperation(Int4 Blk,Int4 start,Int4 end,char *operation)
// ========== trim block in an operational array. ===========
{
        assert(Blk > 0);
        assert(start >= 1 && start <= end);
        char *new_operation=0;
        Int4 o,no,column,blk;
        Int4 trace_length=strlen(operation);
        NEW(new_operation,trace_length+5,char);
        new_operation[0]='E'; no=1;
        for(o=1,column=0,blk=0; operation[o] != 'E'; o++){
          switch(operation[o]){
            case 'M': blk++; if(blk==Blk) column=1;
            case 'm':
                if(blk==Blk && (column < start || column > end)){
                           new_operation[no]='i';
		} else if(blk==Blk && column==start){
                           new_operation[no]='M';
                } else new_operation[no]=operation[o];
                no++; column++; break;
            case 'D': blk++; if(blk==Blk) column=1; 
            case 'd': // deletion in sequence relative to profile.
                if(blk==Blk && (column < start || column > end)){
                           // do nothing in new_operation;
		} else if(blk==Blk && column==start){
                           new_operation[no]='D'; no++;
                } else { new_operation[no]=operation[o]; no++; }
                column++; break;
            case 'i': // insert is between profile blocks;
		new_operation[no]='i'; no++;
                break;
            case 'I': // Insert ('-') within a profile block; delete from seq.
                if(blk==Blk && (column <= start || column > end)){
                   new_operation[no]='i';	// next column == start
		} else new_operation[no]='I'; no++;
                break;
            default:
            // fprintf(stderr,"operations = %s\n",operation);
            print_error("operation( ): input error"); break;
          }
        } // if(blk==Blk) assert(column-1 > end);
        new_operation[no]='E'; no++; new_operation[no]=0;
        // printf("new operation: %s\n",new_operation);
        return new_operation;
}

gsq_typ *gsq_typ::TrimBlk(Int4 Blk,Int4 start,Int4 end,Int4 nblk, Int4 *len, Int4 *sites)
// sites input locations of blocks; sites are reinitialized to new sites.
{
	char debug=0;
	char *operation=this->Operation(nblk,sites,len);
	if(debug){
	  fprintf(stderr,"operation=%s\n",operation);
	  fprintf(stderr,"Blk %d: start...end=%d..%d\n",Blk,start,end);
	}
	char *new_operation=this->TrimBlkInOperation(Blk,start,end,operation); free(operation);
	if(debug) fprintf(stderr," new operation=%s\n",new_operation);
	Int4 numiters=0,numfix;
	do { numfix=this->RmInsertByDeleteOperation(new_operation); numiters++; } while(numfix > 0);
	if(debug && numiters > 1) fprintf(stderr," new new operation=%s\n",new_operation);
	gsq_typ *gsq0 = new gsq_typ[1];
	gsq0->initialize(new_operation,realE,sites);	// copies positions to sites.
	free(new_operation); 
	return gsq0;
}

char	*gsq_typ::RmBlkInOperation(Int4 Blk,char *operation)
// ========== trim block in an operational array. ===========
{
        assert(Blk > 0);
        char *new_operation=0;
        Int4 o,no,blk,trace_length=strlen(operation);
        NEW(new_operation,trace_length+5,char);
        new_operation[0]='E'; no=1;
        for(o=1,blk=0; operation[o] != 'E'; o++){
          switch(operation[o]){
            case 'M': blk++; 
            case 'm':
                if(blk==Blk){ new_operation[no]='i'; }
		else new_operation[no]=operation[o];
                no++; break;
            case 'D': blk++; 
            case 'd': // deletion in sequence relative to profile.
                if(blk != Blk){ new_operation[no]=operation[o]; no++; } break;
            case 'i': // insert is between profile blocks;
		new_operation[no]='i'; no++; break;
            case 'I': // Insert ('-') within a profile block; delete from seq.
                if(blk==Blk) new_operation[no]='i';	
		else new_operation[no]='I'; no++; break;
            default:
            // fprintf(stderr,"operations = %s\n",operation);
            print_error("operation( ): input error"); break;
          }
        } new_operation[no]='E'; no++; new_operation[no]=0;
        // printf("new operation: %s\n",new_operation);
        return new_operation;
}

Int4	gsq_typ::SplitBlkInOperation(Int4 Blk, Int4 *pos, char *operation)
{
        Int4    o,column,blk=0,x=Blk;
        assert(operation[0]=='E');
        for(o=1,column=-1,blk=0; operation[o] != 'E'; o++){
          switch(operation[o]){
            case 'M': case 'D':
                blk++; if(blk == Blk) column=0;
            case 'm': case 'd':
                if(column >= 0) column++;
                if(column == pos[x]){ operation[o]=toupper(operation[o]); x++; blk++; }
                break;
            case 'i': // insert is between profile blocks;
                break; // do nothing
            case 'I': // Insert ('-') within a profile block; delete from seq.
                if(column+1 == pos[x]) operation[o]='i';
                break;
            default:
            // fprintf(stderr,"operations = %s\n",operation);
            print_error("SplitBlkInOperation( ): input error"); break;
          }
        }
	if(pos[blk+1] != 0){
		fprintf(stderr,"blk=%d; pos[blk+1]=%d; column=%d\n",blk,pos[blk+1],column);
		 assert(pos[blk+1]==0);
	}
        return blk;
}

gsq_typ *gsq_typ::SplitBlock(Int4 Blk, Int4 *pos,Int4 nBlks, Int4 *Len, Int4 *sites)
{
	char *operation=this->Operation(nBlks,sites,Len);	// 
	// fprintf(stderr,"operation=%s\n",operation);
	this->SplitBlkInOperation(Blk,pos,operation);
	// fprintf(stderr,"  new operation=%s\n",operation);
	gsq_typ *gsq0 = new gsq_typ[1];
        gsq0->initialize(operation,realE,sites);    // copies positions to sites array.
	return gsq0;
}

static void debug_extend_blk(Int4 blk,Int4 no, Int4 o, char *op, char *nop)
{
	if(o > 5){
           fprintf(stderr," Blk %d op: '%c%c%c%c%c'\n",blk,op[o-4],op[o-3],op[o-2],op[o-1],op[o]);
           fprintf(stderr,"    new op: '%s'\n",nop+(no-5));
	}
}

char	*gsq_typ::ExtendBlkInOperation(Int4 Blk,Int4 AddLeft,Int4 AddRight,char *operation)
// ========== trim block in an operational array. ===========
// 'I' = Insert ('-') within a block; deleted from fake seq.
// 'd' = Deletion within a block.
// 'i' = insert between profile blocks.
// 'M' Add to C-terminal end: Nins== 3; "mmiiiM" --> "mmmiiM" 
//        Nins==0; "mmM" --> "mmdM" or "mdM" --> "mddM
// 'M' Add to N-terminal end: Nins== 3; "mmiiiM" --> "mmiiMm".
//        Nins==0; "mmM" --> "mmDm" or "mdM" --> "mdDm"
// 'D' Add to C-terminal end: Nins == 3; "mmiiiD" --> "mmmiiD"
//        Nins==0; "mmD" --> "mmdD" or "mdD" --> "mddD"
// 'D' Add to N-terminal end: Nins == 3; "mmiiiD" --> "mmmiiMd".
//        Nins==0; "mmD" --> "mmDd" or "mdD" --> "mdDd"
{
        assert(Blk > 0);
        assert(AddLeft >= 0 && AddRight >= 0);
        char	token,*new_operation=0;
        Int4	i,j,o,no,blk,Nins,trace_length=strlen(operation);
        NEW(new_operation,trace_length + AddLeft + AddRight +5,char);
	new_operation[0]='E';
	BooLean NeedToAdd,debug=TRUE; debug=FALSE;

        for(no=1,o=1,Nins=0,blk=0; (token=operation[o]) != 'E'; o++){
	  NeedToAdd=TRUE; 
          switch(token){
            case 'M': case 'D': 
		if(blk==Blk && AddRight){  // Add to C-terminal end... 
		   if(Nins > 0){ new_operation[no-Nins]='m'; } 
		   else { new_operation[no]='d'; no++; }
		   new_operation[no]=token; no++; NeedToAdd=FALSE;
		   if(debug) debug_extend_blk(blk,no,o,operation,new_operation);
		} blk++; 
		if(blk==Blk && AddLeft){   // Add to N-terminal end...
		   if(Nins > 0){ new_operation[no-1]='M'; }
		   else { new_operation[no]='D'; no++; }
		   new_operation[no]=tolower(token); no++;  NeedToAdd=FALSE;
		   if(debug) debug_extend_blk(blk,no,o,operation,new_operation);
		} if(NeedToAdd){ new_operation[no]=token; no++; }
		Nins=0; break;
            case 'i': new_operation[no]='i'; no++; Nins++; break;
            case 'I': case 'd': case 'm': new_operation[no]=operation[o]; no++; break;
            default: print_error("gsq_typ::ExtendBlkInOperation( ): input error"); break;
          }
	} 
	if(blk == Blk && AddRight){  // Add to C-terminal end... 
	   if(Nins > 0){ new_operation[no-Nins]='m'; }
	   else { new_operation[no]='d'; no++; } 
	   if(debug) debug_extend_blk(blk,no,o,operation,new_operation);
	} new_operation[no]='E'; no++; new_operation[no]=0;
        return new_operation;
}

gsq_typ *gsq_typ::ExtendBlk(Int4 Blk,Int4 AddLeft,Int4 AddRight,Int4 nblk, Int4 *len, Int4 *sites)
// sites input locations of blocks; sites are reinitialized to new sites.
{
	char *new_operation=0,*operation=this->Operation(nblk,sites,len);
	char debug=1; debug=0;
	if(debug){
	  fprintf(stderr,"operation=%s\n",operation);
	  fprintf(stderr,"Blk %d: +%d...+%d.\n",Blk,AddLeft,AddRight);
	}
	do {
		new_operation=this->ExtendBlkInOperation(Blk,AddLeft,AddRight,operation);
		free(operation); operation=0;
		AddLeft = MAXIMUM(Int4,AddLeft-1,0); AddRight = MAXIMUM(Int4,AddRight-1,0);
		if(AddLeft > 0 || AddRight > 0) { operation = new_operation; }
	} while(AddLeft > 0 || AddRight > 0);
	if(debug){ fprintf(stderr," new operation=%s\n",new_operation); }
	Int4 numiters=0,numfix;
	do { numfix=this->RmInsertByDeleteOperation(new_operation); numiters++; } while(numfix > 0);
	if(debug && numiters > 1) fprintf(stderr," new new operation=%s\n",new_operation);
	gsq_typ *gsq0 = new gsq_typ[1];
	gsq0->initialize(new_operation,realE,sites);	// copies positions to sites.
	free(new_operation);
	return gsq0;
}

char	*gsq_typ::Operation3(FILE *fp, a_type A)
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
			insert(i),f2r[i],operation);
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

char    *gsq_typ::AddColumnToOperation(Int4 Blk, Int4 start, Int4 Length, char *operation)
// ========== Add a column to an operational array. ===========
{
        assert(Blk > 0);
        assert(start > 0);
        Int4 o,no,column,blk,x,add,trace_length=strlen(operation);
        char token,last='E',*new_operation; NEW(new_operation,trace_length + Length + 5,char);
        new_operation[0]='E';
        for(no=1,o=1,add=0,column=0,blk=0; operation[o] != 'E'; o++){
	  token=operation[o];
          switch(token){
            case 'M': blk++; column=1;
            case 'm':
                if(blk==Blk && column == start){
		   add=Length; new_operation[no]=token; no++;
		} else if(blk==Blk && add > 0){
		  while(add > 0){ new_operation[no]='d'; no++; add--; }
		  new_operation[no]=token; no++;
		} else { new_operation[no]=token; no++; }
                column++; break;
            case 'D': blk++; column=1;
            case 'd': // deletion in sequence relative to profile.
                if(blk==Blk && column == start){
		   add=Length; new_operation[no]=token; no++;
		} else if(blk==Blk && add > 0){
		  while(add > 0){ new_operation[no]='d'; no++; add--; }
		  new_operation[no]=token; no++;
		} else { new_operation[no]=token; no++; }
                column++; break;
            case 'i': // insert is between profile blocks;
		assert(add <= 0);
                new_operation[no]=token; no++;
                break;
            case 'I': // Insert ('-') within a profile block; delete from seq.
		if(add > 0){ new_operation[no]='m'; no++; add--; }
                else { new_operation[no]=token; no++; }
                break;
            default:
            // fprintf(stderr,"operations = %s\n",operation);
            print_error("operation( ): input error"); break;
          } last=token;
        }
#if 0	// this need not hold; delete it...
	if(blk==Blk && column-1 <= Length){
            	fprintf(stderr,"operation = %s\n",operation);
            	fprintf(stderr,"new operation = %s\n",new_operation);
            	fprintf(stderr,"column=%d; Length=%d; Blk=%d; start=%d\n",
			column,Length,Blk,start);
		assert(column-1 > Length);
	}
#else
	assert(no <= (trace_length + Length)); 
#endif
        new_operation[no]='E'; no++; new_operation[no]=0;
        // printf("new operation: %s\n",new_operation);
        return new_operation;
}

char *ReverseOperation(char *operation)
{
        Int4    i,j,len;
        char    *rev_op;
        len=strlen(operation)-1;
        NEW(rev_op,len + 5,char);
        rev_op[0]='E'; rev_op[len]='E';
        for(i=1, j=len-1; i < len ; i++,j--){
                if(operation[j] == 'I') rev_op[i] = operation[j];
                else rev_op[i] = tolower(operation[j]);
        }
        for(i=1; i <= len; i++){
                if (rev_op[i]=='m'){ rev_op[i]='M'; break; }
                else if (rev_op[i]=='d'){ rev_op[i]='D'; break; }
        }
        if(i > len){
            assert(i <= len);
            print_error("ReverseOperation(): error in operation array");
        } return rev_op;
}

gsq_typ *gsq_typ::InsertColumns(Int4 theBlk, Int4 Start, Int4 Length, Int4 nBlks,
                        Int4 *len, Int4 *sites)
{
        Int4    numfix,start;
        char *new_operation,*operation=this->Operation(nBlks,sites,len);
        if(Length < 0){ // then insert column(s) before Start.
          assert(nBlks == 1); Length = -Length;
          // PutShortSeqID(stderr,realE);
          //fprintf(stderr," %d: operation=%s\n",Start,operation);
          char *ro=ReverseOperation(operation); free(operation);
          // fprintf(stderr,"rev_opera=%s\n",ro);
          start=len[1]-Start + 1;
          // if(this->IsDeleted(Start)) fprintf(stderr,"%d = '-'\n",Start);
          operation=AddColumnToOperation(theBlk,start,Length,ro); free(ro);
          // fprintf(stderr,"Add column %d opera=%s\n",start,operation);
          new_operation=ReverseOperation(operation); free(operation);
          // fprintf(stderr,"new_opera=%s\n\n",new_operation);
        } else {
          new_operation=AddColumnToOperation(theBlk,Start,Length,operation); free(operation);
        }
        do { numfix=this->RmInsertByDeleteOperation(new_operation); } while(numfix > 0);
        // fprintf(stderr,"new operation: %s\n",new_operation);
        gsq_typ *gsq0 = new gsq_typ[1];
        gsq0->initialize(new_operation,realE,sites);    // copies positions to sites array.
        free(new_operation);
        return gsq0;
}

