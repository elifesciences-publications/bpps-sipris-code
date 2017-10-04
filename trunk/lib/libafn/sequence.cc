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

#include "sequence.h"

const unsigned char default_overlap = 20;

static unsigned char *alloc_seq(Int4 length)
{
	unsigned char	*seq;
	Int4	i;

	MEW(seq,length+2*default_overlap+4, unsigned char);
	seq += default_overlap;
	for(i=-default_overlap; i<=0; i++) seq[i]=0;
	for(i=length + default_overlap; i > length; i--) seq[i]=0;
	return seq;
}

static void dealloc_seq(unsigned char *seq)
{ seq -= default_overlap; free(seq); }

void	LabelSeq(e_type E)
{
	if(E->kingdom == 0){
		assert(E->phylum == 0);
		E->phylum=AllocString("Unknown"); E->kingdom='u';
	} else E->kingdom=tolower(E->kingdom);
}

void	UnLabelSeq(e_type E)
{
	if(E->kingdom == 0){
		assert(E->phylum == 0);
		E->phylum=AllocString("Unknown"); E->kingdom='U';
	} else E->kingdom=toupper(E->kingdom);
}

long	IsFastaFormatSeq(FILE *fp,long &MaxLen)
// Check to see whether this file contains fasta sequences in valid format.
// (Instead of a fasta format alignment.)
// if valid fasta file the return number of fasta sequences;
// if appears to be fasta alignment file then return negative of number of sequences
// else return 0;
// set MaxLen to the longest sequence found.
{
	long	len=0,N=0,min_len=LONG_MAX,max_len=0;
	long	gaps=0;
	char	c;

	// scan over space characters to first defline:
	do { c=fgetc(fp); } while(isspace(c)); if(c==EOF) return 0;
	if(c!='>') return 0; 
  do {	// c=='>' at this point...

	// scan over defline:
	do { c=fgetc(fp); } while(c!='\n' && c!=EOF);
	if(c==EOF) return 0;

	// scan over sequences:
	N++; len=0;
	while(c=fgetc(fp)){
	   if(isalpha(c)){ len++; }
	   else if(c == '-'){ gaps++; len++; }
	   else if(isspace(c)) ;   // do nothing.
	   else if(c=='>' || c==EOF){	// done with this sequence..
	     if(len ==0) return 0; // defective file.
	     else {
	       if(len > max_len)max_len=len;
	       if(len < min_len)min_len=len;
	     }
	     break;	// do another loop;
	   } else return 0;	// error occurred 
	}
   } while(c != EOF);
   if(gaps){
	if(max_len==min_len){ MaxLen=max_len; return -N; } // fasta format alignment.
	else return 0;			// invalid input.
   } else {
	MaxLen=max_len; return N;
   }
}

void	AllocXSeq(e_type E)
{
	if(E->S != E->X) print_error("AllocXSeq( ): already allocated!");
	E->X = alloc_seq(E->n); E->xed=TRUE; 
}

e_type	AllocSeq(UInt4 I, Int4 length)
{
	e_type	E;

	NEW(E,1,sequence_type);
	E->S = alloc_seq(length);
	E->overlap=default_overlap; E->offset=E->extend=0;
	E->kingdom=0; E->phylum=0; E->genetic_code=0;
	E->info=NULL;  E->n = length;
	E->next = 0;
        E->I = I; E->X = E->S; E->xed = FALSE; E->rff=0;
        return E;
}

void	UnXSeq(e_type E)
/*** remove 'X'ed out regions ****/
{
	if(E==NULL) seq_error("UnXSeq( ) can't UnX null sequence.");
	if(E->xed){ dealloc_seq(E->X); E->X=E->S; E->xed = FALSE; } 
}

void	SeverNextSeqList(e_type E) { E->next = 0; }

Int4	ConcatSeqList(e_type headE, e_type tailE)
// attach TailE to headE and return the length of headE...
{
	e_type endE,next,last;
	Int4 i;
	assert(headE != 0 && tailE != 0);
	for(i=1,last=headE,next=NextSeq(headE); next;  next=NextSeq(next)){
		last=next; i++;
	} assert(last->next == 0);
	last->next = tailE;
	return i;
}

void	NilSeq(e_type E)
{
	if(E!=NULL){
		if(E->xed) dealloc_seq(E->X);
		dealloc_seq(E->S);
		if(E->info!=NULL)free(E->info); 
		if(E->rff) delete E->rff; 
		if(E->phylum) free(E->phylum);
		free(E);
	}
}

BooLean StringInSeqInfo(char *phrase, e_type E)
// return true if info line for sequence contains phrase
{ if(strstr(E->info,phrase)) return TRUE; else return FALSE; }

e_type	MkSeq(char *defline, Int4 len, unsigned char *seq)
{
	e_type	E;
	Int4 	i;

	E = AllocSeq(1, len);
	E->info = get_info_seq(E, defline);
	for(i=1;i<=(Int4)E->n; i++) E->S[i]= seq[i];
	return E;
}

e_type  MakeSeq(char *id, char *descript, Int4 offset, Int4 extend, Int4 len, 
	unsigned char *seq)
{ return MakeSeq(id,descript, offset, extend, len, seq,0,0); }

void	TaxAssignSeq(char *phylum, char kingdom, e_type E)
{
	if(E->phylum) free(E->phylum);
	E->phylum=AllocString(phylum);
	E->kingdom=kingdom;
}

e_type  MakeSeq(char *id, char *descript, Int4 offset, Int4 extend, Int4 len, 
	unsigned char *seq,char *phylum, char kingdom)
{
	Int4	lenstr=strlen(id) + strlen(descript) + 20;
	char	*defline;
	e_type  E;

	defline = new char[lenstr+100];
	if(phylum) sprintf(defline,"%s {|%d(%d)|<%s(%c)>}%s",
				id,offset,extend,phylum,kingdom,descript);
	else sprintf(defline,"%s {|%d(%d)|}%s",id,offset,extend,descript);
	E = MkSeq(defline, len, seq);
	delete []defline;
	return E;
}	

Int4    ReadSeqStr(FILE *fptr, register char *defline, register unsigned char *seq,
        Int4 max, a_type A)
// This is only used by transitive blast routine; can probably be eliminated...
{
        register Int4	length,i;
	register char	c='x',last_c='\n';

	do {
	   if((c=fgetc(fptr)) == '>' && last_c=='\n') break;
	   if(c == EOF) return 0;
	   else if(!isspace(c)) print_error("ReadSeqStr( ): fasta input error"); 
	   else last_c=c;
	} while(TRUE); 
        for(i=0; (c=fgetc(fptr))!=EOF; i++){
                if(c=='\n') break; 
                else if(i < DEFAULT_INFO_LINE_LEN_SEQ && isprint(c)) defline[i] = c;
                else { while((c=fgetc(fptr))!='\n' && c != EOF) ; break;}
	} defline[i] = '\0';

	while((c=fgetc(fptr)) == '+'){ // scan over additional information.
	   ScanOverRFF(fptr); // will just skip over this information...
	} if(ungetc(c,fptr) == EOF) print_error("ReadSeqStr( )-> fasta input error");

        for(last_c=0,length=1; (c=fgetc(fptr))!=EOF; ){
                if(c == '>') {
			if(last_c != '\n') print_error("ReadSeqStr( )-> fasta input error");
			ungetc(c,fptr); return (Int4) (length-1); 
		} else if(isalpha(c)) seq[length++] = AlphaCode(c,A);
		// else if(c == '-') seq[length++] = AlphaCode(c,A);
		else if(c == '-') seq[length++] = AlphaCode('X',A);
                if(length >= max) {
                        fprintf(stderr,">%s\n",defline);
                        print_error("input sequence too long");
                } last_c=c;
        } 
	if(length <= 1) print_error("ReadSeqStr( )-> fasta input error");
	return (Int4) (length-1);
}

e_type	ReadSeqStrII(FILE *fptr, UInt4 I, register char *defline, 
	register unsigned char *seq, Int4 max, a_type A)
{
        register Int4	length,i;
	register char	c='x',last_c='\n';
	e_type		E;
	rff_typ	*rff=0;

	do {
	   if((c=fgetc(fptr)) == '>' && last_c=='\n') break;
	   if(c == EOF) return 0;
	   else if(!isspace(c)) {
		do {
		   if(c == EOF) break;
		   else fprintf(stderr,"%c",c);
		} while((c=fgetc(fptr)) != '\n');
		fprintf(stderr,"\n");
		print_error("ReadSeqStrII( ): fasta input error1"); 
	   }
	   else last_c=c;
	} while(TRUE); 
        for(i=0; (c=fgetc(fptr))!=EOF; i++){
                if(c=='\n') break; 
                else if(i < DEFAULT_INFO_LINE_LEN_SEQ && isprint(c)) defline[i] = c;
                else { while((c=fgetc(fptr))!='\n' && c != EOF) ; break;}
	} defline[i] = '\0';

	while((c=fgetc(fptr)) == '+'){ // scan over additional information.
	   ungetc(c,fptr); // put back '+';
	   rff = new rff_typ(fptr);
	} if(ungetc(c,fptr) == EOF) print_error("ReadSeqStrII( )-> fasta input error2");

        for(last_c=0,length=1; (c=fgetc(fptr))!=EOF; ){
                if(c == '>') {
			if(last_c != '\n') print_error("ReadSeqStrII( )-> fasta input error3");
			ungetc(c,fptr); break;
		} else if(isalpha(c)) seq[length++] = AlphaCode(c,A);
		// else if(c == '-') seq[length++] = AlphaCode(c,A);
		else if(c == '-') seq[length++] = AlphaCode('X',A);
                if(length >= max) {
                        fprintf(stderr,">%s\n",defline);
                        print_error("input sequence too long");
                } last_c=c;
        } 
	if(length <= 1) print_error("ReadSeqStrII( )-> fasta input error4");
	// fprintf(stderr,"length=%d\n",length-1);
	E = StringToSeq2(seq, length-1, defline, I);
	if(rff){
		E->rff = rff;
		if(rff->LengthS() != E->n) print_error("ReadSeqStrII( )-> rff input error5"); 
	} return E;
}

e_type	ReadSeq(FILE *fptr, Int4 I, Int4 size, a_type A)
{
	char	id[MAX_SEQ_DEFLINE_LENG+4];
	unsigned char	*seq;
	Int4	length,i;
	char	c='x',last_c='\n';
	e_type	E;
	rff_typ	*rff=0;

	do {
	   if((c=fgetc(fptr)) == '>' && last_c=='\n') break;
	   if(c == EOF) return 0;
	   else if(!isspace(c)) print_error("ReadSeq( ): fasta input error 1"); 
	   else last_c=c;
	} while(TRUE); 
        for(i=0; (c=fgetc(fptr))!=EOF; i++){ 
		if(c=='\n') { break; }
		else if(i< DEFAULT_INFO_LINE_LEN_SEQ && isprint(c)) id[i] = c;
                else { while((c=fgetc(fptr))!='\n' && c != EOF) ; break;}
	} id[i] = '\0';
	E = AllocSeq(I,size); seq = E->S;

	while((c=fgetc(fptr)) == '+'){ // scan additional information.
	   ungetc(c,fptr); // put back '+';
	   rff = new rff_typ(fptr);
	} if(ungetc(c,fptr) == EOF) print_error("ReadSeq( )-> fasta input error 2");
// fprintf(stderr,"%d: c = %c\n",I,c);

        for(last_c=0,length=1; (c=fgetc(fptr))!=EOF; ){ 
// fprintf(stderr,"length=%d: c = %c\n",length,c);
                if(c == '>') {
			if(last_c != '\n') print_error("ReadSeq( ): fasta input error 3");
			else ungetc(c,fptr); break;
		} else if(isalpha(c)) { seq[length++] = AlphaCode(c,A); }
		// else if(c == '-') seq[length++] = AlphaCode(c,A);
		else if(c == '-') seq[length++] = AlphaCode('X',A);
		else last_c=c;
        } length--;
	if(length != size) {
		fprintf(stderr,"I=%d; size = %d; length = %d\n",I,size,length);
		fprintf(stderr,"possible database mung!\n");
		assert(length == size);
		seq_error("size and length inconsistency");
	}
	E->info = get_info_seq(E, id);
	if(rff){
		E->rff = rff;
		if(rff->LengthS() != E->n) print_error("ReadSeq( )-> rff input error"); 
	} return E;
}

BooLean	NonNullIdentSeqs(register e_type E1, register e_type E2)
/** If E1 has the same sequence as E2 return TRUE; else return FALSE **/
{
	register Int4	s,r1,r2;

	if(E1->n != E2->n) return FALSE;
	for(s=1; s<= (Int4) E1->n; s++){
		r1 = E1->X[s]; r2 = E2->X[s]; 
		if(r1 != 0 && r2 != 0 && r1 != r2) return FALSE;
	}
	return TRUE;
}

e_type	StringToSeq2(unsigned char *seq, Int4 length, char *info, UInt4 I)
/*** convert a character string to a sequence **/
{
	char	id[MAX_SEQ_DEFLINE_LENG+4];
	unsigned char	*seq0;
	Int4	i,c;
	e_type	E;

        for(i=0; (c=info[i])!=0; i++){ 
		if(c=='\n') break; 
		else if(i > DEFAULT_INFO_LINE_LEN_SEQ || !isprint(c)) break;
		else id[i] = c;
	} 
	if(i > 0 && id[i-1] != ' '){ id[i]=' '; i++; } 
	id[i] = '\0';
	E = AllocSeq(I,length); seq0 = E->S;
        for(i=1; i <= length; i++) seq0[i] = seq[i];
	E->info = get_info_seq(E,id);
	return E; 
}

e_type	StringToSeq(char *buffer, char *info, UInt4 I, a_type A)
/*** convert a character string to a sequence **/
{
	char	id[MAX_SEQ_DEFLINE_LENG+4];
	unsigned char	*seq;
	Int4	length,i,c;
	e_type	E;

        for(i=0; (c=info[i])!=0; i++){ 
		if(c=='\n') break; 
		else if(i > DEFAULT_INFO_LINE_LEN_SEQ || !isprint(c)) break;
		else id[i] = c;
	} 
	if(i > 0 && id[i-1] != ' '){ id[i]=' '; i++; } 
	id[i] = '\0';
        for(i=0; (c=buffer[i])!=0; i++) if(!isalpha(c)) break; 
	E = AllocSeq(I,i); seq = E->S; 
        for(i=0; (c=buffer[i])!=0; ){ 
	    if(isalpha(c)) { i++; seq[i] = AlphaCode(c,A); }
	    else { break; }
        } length = i;
	E->info = get_info_seq(E, id);
	return E; 
}

Int4	SeqToString(register char *buffer, e_type E, register a_type A)
/*** convert a sequence to a character string **/
{
	register Int4	i;
	register unsigned char *s;
	
	s = SeqPtr(E); 
	for(s++, i=0; i< LenSeq(E); i++) buffer[i] = AlphaChar(s[i],A);
	buffer[i]=0;
	return i;
}

void    PutDeleteSeq(FILE *fptr,Int4 start_del, Int4 end_del, e_type E,a_type A)
{
        Int4    i,j;

        if(LenSeq(E) <= (end_del - start_del + 1)) return;
        fprintf(fptr,">del(%d..%d)_",start_del,end_del); PutSeqInfo(fptr,E);
        for(i=j=1; (Int4)j <= (Int4)LenSeq(E); j++) {
          if(j < start_del || j > end_del){
           fprintf(fptr,"%c",AlphaChar(ResSeq(j,E),A));
           if(i%50 == 0) fprintf(fptr,"\n"); i++;
          }
        } fprintf(fptr,"\n\n");
}

BooLean	PutSuperSeq(FILE *fp, Int4 left, Int4 right, e_type sE, e_type oE,
	a_type A)
{
	Int4	s,start,end,len_s;
	unsigned char *ss,*os;
	BooLean	flag;

	if(sE->n > oE->n) return FALSE;
	ss=SeqPtr(sE) ; end = LenSeq(oE) - LenSeq(sE);
	len_s=LenSeq(sE);
	for(start=0; start <=end; start++){
	   os=SeqPtr(oE) + start; flag=TRUE;
	   for(s=1; s<= (Int4) len_s; s++){
		if(os[s] != ss[s]) { flag=FALSE; break; }
	   }
	   if(flag){
		PutSubSeq(fp,start-left,len_s+start+right,oE,A);
		return TRUE;
	   } 
	}
	return FALSE;
}

#if 0
BooLean	OverlappingSeq(Int4 *offset, Int4 MinOverlap, e_type E1, e_type E2)
{ return OverlappingSeq(offset,MinOverlap,0,E1,E2); }

BooLean	OverlappingSeq(Int4 *offset, Int4 MinOverlap, Int4 MaxMisMatch, e_type E1, e_type E2)
/*************************************************************
return TRUE if E1 & E2 are overlapping fragments of the same sequence.

Start:					end1=21; end2=24
       	E1  ----+----+----+----+-     	start1=end1 -MinOverlap +1 = 21-5 +1=17;
                            |||||	start2 = 1;
	E2                  ----+----+----+----+----           
		:	:	:		(*offset= start1 - start2 = 17 -1 = 16)
		:	:	:		( add 16 to second seq)
	E1  ----+----+----+----+-	
	    |||||||||||||||||||||		start1=start2=1 offset = 1 - 1 = 0
	E2  ----+----+----+----+----
		:	:	:	
		:	:	:	
End:					
	E1                     ----+----+----+----+-
                               |||||	
	E2  ----+----+----+----+----    start2=end2 - MinOverlap +1 = 24 -5 + 1 = 20       
					offset = 1 - 20 = -19 (add 19 to first seq)

 *************************************************************/
{
	Int4	start1,start2,end1,end2;
	Int4	s1,s2,e1,e2,stop1,stop2;
	unsigned char	*sq1,*sq2;
	Int4    mismatches;

	end1=LenSeq(E1); end2=LenSeq(E2);
	sq1=SeqPtr(E1); sq2=SeqPtr(E2);
	stop1=end1 - MinOverlap +1;
	stop2=end2 - MinOverlap +1;
	for(start1=stop1 ; start1 >= 1; start1--){
	   for(start2=1; start2 <= stop2; start2++){
		for(s1=start1, s2=start2; s1 <= end1 && s2 <= end2; s1++,s2++){
			mismatches=0;
	       		if(sq1[s1] != sq2[s2]){ 
			   mismatches++;
			   if(sq1[s1] != 0 && sq2[s2] != 0 && mismatches > MaxMisMatch) break; 
			}	// allow for missing residues in pdb files...
		}
		if(s1 > end1 || s2 > end2){ *offset= start1 - start2;  return TRUE; }
	   }
	}
	return FALSE;
}
#endif

BooLean	OverlappingSeq(Int4 &offset, e_type E1, e_type E2, Int4 MinOverlap, Int4 MaxMisMatch)
/*************************************************************
return TRUE if E1 & E2 are overlapping fragments of the same sequence.
This sets offset > 0 if query starts before pdbseq else it sets dcaOS <= 0.
It also allows for 'X' residues in pdb files...

Start:					end1=21; end2=24; MinOverlap=5;
       	E1  ----+----+----+----+-     	start1=end1-MinOverlap+1 = 21-5+1=17;
                            |||||	start2 = 1;
	E2                  ----+----+----+----+----           
		:	:	:		(offset= start1 - start2 = 17 -1 = 16)
		:	:	:		( add 16 to second seq)
	E1  ----+----+----+----+-	
	    |||||||||||||||||||||		start1=start2=1 offset = 1 - 1 = 0
	E2  ----+----+----+----+----
		:	:	:	
		:	:	:	
End:					
	E1                     ----+----+----+----+-
                               |||||	
	E2  ----+----+----+----+----    start2=end2 - MinOverlap +1 = 24 -5 + 1 = 20       
					offset = 1 - 20 = -19 (add 19 to first seq)
 (move to sequence.cc eventually)
 *************************************************************/
{
	Int4	start1,start2,end1,end2;
	Int4	s1,s2,e1,e2,stop1,stop2;
	unsigned char	*sq1,*sq2;
	Int4    mismatches=0;

	end1=LenSeq(E1); end2=LenSeq(E2); sq1=SeqPtr(E1); sq2=SeqPtr(E2);
	stop1=end1 - MinOverlap +1; stop2=end2 - MinOverlap +1;
	// for(start1=stop1,start2=1; start1 >= 1 && start2 <= stop2; start1--,start2++)
        for(start1=stop1 ; start1 >= 1; start1--)
	   {
           for(start2=1; start2 <= stop2; start2++){

		mismatches=0;
		for(s1=start1, s2=start2; s1 <= end1 && s2 <= end2; s1++,s2++){
	       		if(sq1[s1]==0 || sq2[s2] == 0) continue;	// one or two == 'X'.
	       		if(sq1[s1] != sq2[s2]){ 
			   mismatches++;
			   if(sq1[s1] != 0 && sq2[s2] != 0 && mismatches > MaxMisMatch) break; 
			}	// allow for missing residues in pdb files...
		}
		if(s1 > end1 || s2 > end2){ offset= start1-start2; return TRUE; }
#if 0
		else {
		   fprintf(stderr,"start1=%d; end1=%d; offset=%d; start2=%d; end2=%d\n",
			start1,end1,start1-start2,start2,end2);
		}
#endif
	   }
	} return FALSE;
}

BooLean	IsSameSeqID(e_type E1, e_type E2)
{
	char	*info1=E1->info,*info2=E2->info;
	Int4	i;
	for(i=0; info1[i] && info2[i]; i++){
		if(info1[i] != info2[i]) return FALSE;
		if(info1[i] == ' ' && info2[i] == ' ') return TRUE;
	} if(info1[i]==0 && info2[i]==0) return TRUE;
	else return FALSE;
}

char    IsSameSeq(e_type E1, e_type E2,Int4 *Start,Int4 MinOverlap,BooLean IgnoreX)
// find out whether or not E1 and E2 are the same sequence.
// return 1 if seq E1 lacks an N-terminal extension.
// return 2 if seq E2 lacks an N-terminal extension.
{
        Int4    st,sp,sb,end,lenSb,lenSp,test,NumX,Overlap;
        char    offset_seq=0;
        unsigned char   *sup,*sub;      // superseq and subseq
        BooLean IsSub=FALSE;

   if(LenSeq(E1) < MinOverlap || LenSeq(E2) < MinOverlap) return 0;
   for(test=1; test <= 2; test++){
      if(test==1){	// 1. start from sequence E1 and inch along seq. E2.
   	sup=SeqPtr(E2); sub=SeqPtr(E1); 
	lenSp=LenSeq(E2); lenSb=LenSeq(E1); offset_seq=1;
	end=LenSeq(E2)-LenSeq(E1)+1;
      } else { 		// 2. if failed: start from sequence E2 and inch along seq. E1.
	sup=SeqPtr(E1); sub=SeqPtr(E2);
	lenSp=LenSeq(E1); lenSb=LenSeq(E2); offset_seq=2;
	end=LenSeq(E1)-LenSeq(E2)+1;
      }
      if(IgnoreX){ // ignore 'X' residues (in pdb files)
        for(st=1; st<= end; st++){
	  if((lenSp-st+1) < MinOverlap) break;
          for(NumX=0,IsSub=TRUE,sp=st,sb=1; sb <= lenSb && sp <= lenSp; sb++,sp++){
            if(sup[sp]==0 || sub[sb]==0){ NumX++; continue; } 
	    if(sup[sp] != sub[sb]){ IsSub=FALSE; break; }
          }
          if(IsSub){
		Overlap=lenSb-st+1-NumX;
		if(Overlap < MinOverlap) break;	// can only get shorter from here...
		*Start=st-1; return offset_seq; 
	  }
        }
      } else {
        for(st=1; st<= end; st++){
	  if((lenSp-st+1) < MinOverlap) break;
          for(IsSub=TRUE,sp=st,sb=1; sb <= lenSb && sp <= lenSp; sb++,sp++){
            if(sup[sp] != sub[sb]){ IsSub=FALSE; break; }
          }
          if(IsSub){
		Overlap=lenSb-st+1;
		if(Overlap < MinOverlap) break;	// can only get shorter from here...
		*Start=st-1; return offset_seq; 
	  }
        }
      }
   } return 0;
}

static Int4 is_same_seq_fast(register unsigned char *sup, register unsigned char *sup_end,
	register unsigned char *sub, register unsigned char *sub_end)
// returns the overlap between sequences...
{
	register Int4 NumX=0;

#if 0
	for( ; sup != sup_end && sub != sub_end; sup++,sub++){
            if(*sup==0 || *sub==0){ NumX++; continue; } 
	    if(*sup != *sub){ return -1; } 
        } return NumX;
#else
	while(sup != sup_end && sub != sub_end){
	    if(*sup != *sub){
		if(*sup!=0 && *sub!=0){ return -1; } else NumX++; 
	    } sup++,sub++; 
	} return NumX;
#endif
}

char    IsSameSeqFast(e_type E1, e_type E2,Int4 *Start,Int4 MinOverlap)
{ Int4 NumX; return IsSameSeqFast(E1, E2,Start,&NumX, MinOverlap); }

char    IsSameSeqFast(e_type E1, e_type E2,Int4 *Start,Int4 *RtnNumX,Int4 MinOverlap)
// find out whether or not E1 and E2 are the same sequence.
// return 1 if seq E1 lacks an N-terminal extension.
// return 2 if seq E2 lacks an N-terminal extension.
// Sets Start to the position in one sequence corresponding to the start of the other.
{
        Int4    st,end,lenSb,lenSp,NumX;
        unsigned char   *sup,*sub;      // superseq and subseq

	*RtnNumX=0;
	if(LenSeq(E1) < MinOverlap || LenSeq(E2) < MinOverlap) return 0;
	
	end=LenSeq(E2)-MinOverlap+1;
	//          st=end                            st=1
	//    E2 ----+---------+   <-- to...from <--   +---------+-------  E2
	//           |--MinOL--|  	               |--MinOL--|
	//    E1     +---------+-----                  +---------+------   E1
	if(end > 0){
   	  sup=SeqPtr(E2); lenSp=LenSeq(E2);	sub=SeqPtr(E1); lenSb=LenSeq(E1); 
          for(st=1; st<= end; st++){
	    // if((lenSp-st+1) < MinOverlap) break;	// guarranteed to be long enough.
	    NumX=is_same_seq_fast(sup+st,sup+lenSp+1,sub+1,sub+lenSb+1);
	    if(NumX >=0){			// sequences match!
	      // if((lenSb-st+1-NumX) < MinOverlap) break;	// can only get shorter from here...
	      // if((end-st < NumX)) return 0;  // short perfect match; assume no perfect match below.
	      if((end-st < NumX)) continue; // could be due to a string of X residues on end.
	      *Start=st-1; *RtnNumX=NumX; return 1; 
	    }
          }
	}
	//  E2  +---------+-------                         +---------+-------  E2 
	//      |--MinOL--|                                |--MinOL--|
	//  E1  +---------+------  --> from..to -->  ------+---------+         E1
	//     st=1                                        st=end
	end=LenSeq(E1)-MinOverlap+1;
	if(end > 0){
	  sup=SeqPtr(E1); lenSp=LenSeq(E1);	sub=SeqPtr(E2); lenSb=LenSeq(E2); 
          for(st=1; st<= end; st++){
	    // if((lenSp-st+1) < MinOverlap) break;
	    NumX=is_same_seq_fast(sup+st,sup+lenSp+1,sub+1,sub+lenSb+1);
	    if(NumX >=0){			// sequences match!
	      // if((lenSb-st+1-NumX) < MinOverlap) break;	// can only get shorter from here...
	      //if((end-st < NumX)) return 0;	// perfect match but too short...
	      if((end-st < NumX)) continue;	// could be due to a string of X residues on end.
	      *Start=st-1; *RtnNumX=NumX; return 2; 
	    }
	  }
	} return 0;
}

char    IsSameSeq0(e_type E1, e_type E2,Int4 *Start,BooLean IgnoreX)
// find out whether or not E1 and E2 are the same sequence.
// return 1 if seq. E1 lacks an N-terminal extension.
// return 2 if seq. E2 lacks an N-terminal extension.
{
        Int4    st,sp,sb,end,nSb,test;
        char    offset_seq=0;
        unsigned char   *sup,*sub;      // superseq and subseq
        BooLean IsSub=FALSE;

   // print_error("WARNING: this procedure can return a sequence of length 1");
   for(test=1; test <= 2; test++){
      if(test==1){	// 1. start from sequence E1 and inch along seq. E2.
   	sup=SeqPtr(E2); sub=SeqPtr(E1); nSb=LenSeq(E1); offset_seq=1;
	end=LenSeq(E2)-LenSeq(E1)+1;
      } else { 		// 2. if failed: start from sequence E2 and inch along seq. E1.
	sup=SeqPtr(E1); sub=SeqPtr(E2); nSb=LenSeq(E2); offset_seq=2;
	end=LenSeq(E1)-LenSeq(E2)+1;
      }
      if(IgnoreX){ // ignore 'X' residues (in pdb files)
        for(st=1; st<= end; st++){
          for(IsSub=TRUE,sp=st,sb=1; sb <= nSb; sb++,sp++){
            if(sup[sp] && sub[sb] && sup[sp] != sub[sb]){ IsSub=FALSE; break; }
          }
          if(IsSub){ *Start=st-1; return offset_seq; }
        }
      } else {
        for(st=1; st<= end; st++){
          for(IsSub=TRUE,sp=st,sb=1; sb <= nSb; sb++,sp++){
            if(sup[sp] != sub[sb]){ IsSub=FALSE; break; }
          }
          if(IsSub){ *Start=st-1; return offset_seq; }
        }
      }
   } return 0;
}

char	IsSubSeq(e_type E1, e_type E2)
{ Int4	st; char rtn = IsSubSeq(E1,E2,&st,FALSE); return rtn; }

char	IsSubSeq(e_type E1, e_type E2,Int4 *Start,BooLean IgnoreX)
// If E1 is a subseq of E2 return 1; else if E2 < E1 return 2; 
// if E1 == E2 return 3; else return 0.
{
	Int4	st,sp,sb,end,nSb;
	char	subsq=0;
	unsigned char	*sup,*sub;	// superseq and subseq
	BooLean	IsSub=FALSE;

	if(LenSeq(E1) <= LenSeq(E2)){
	  sup=SeqPtr(E2); sub=SeqPtr(E1); nSb=LenSeq(E1); 
	  if(LenSeq(E1)==LenSeq(E2)) subsq=3; else subsq=1; 
	  end=LenSeq(E2)-LenSeq(E1)+1;
	} else { 
	  sup=SeqPtr(E1); sub=SeqPtr(E2); nSb=LenSeq(E2); subsq=2; 
	  end=LenSeq(E1)-LenSeq(E2)+1;
	}
   if(IgnoreX){	// ignore 'X' residues (in pdb files)
	for(st=1; st<= end; st++){
	  for(IsSub=TRUE,sp=st,sb=1; sb <= nSb; sb++,sp++){
	    if(sup[sp] && sub[sb] && sup[sp] != sub[sb]){ IsSub=FALSE; break; }
	  } 
	  if(IsSub){ *Start=st-1; return subsq; }
	}
   } else {
	for(st=1; st<= end; st++){
	  for(IsSub=TRUE,sp=st,sb=1; sb <= nSb; sb++,sp++){
	    if(sup[sp] != sub[sb]){ IsSub=FALSE; break; }
	  } 
	  if(IsSub){ *Start=st-1; return subsq; }
	}
   }
	return 0;
}

BooLean FastIdentSeqs(register e_type E1, register e_type E2)
/** If E1 has the same sequence as E2 return TRUE; else return FALSE **/
{
        register Int4    s;
        if(E1->n != E2->n) return FALSE;
        for(s=1; s<= E1->n; s++) if(E1->S[s] != E2->S[s]) return FALSE;
        return TRUE;
}

BooLean	IdentSeqs(e_type E1, e_type E2)
/** If E1 has the same sequence as E2 return TRUE; else return FALSE **/
{
	Int4	s;

	if(E1->n != E2->n) return FALSE;
	for(s=1; s<= (Int4) E1->n; s++){
		if(E1->X[s] != E2->X[s]) return FALSE;
	}
	return TRUE;
}

e_type	MergeSeqs(e_type *E)
/** merge array of sequences into one Int4 sequence **/
{
	e_type mE;
	Int4  length,s,i,n;
	UInt4 I;

	if(E[0]==NULL) return NULL;
	for(I=length=n=0; E[n] != NULL; n++){
		length += E[n]->n;
		I = MAXIMUM(UInt4 ,I,E[n]->I);
	} I++;
	if(n==1) return CopySeq(E[0]);
	mE = EmptySeq(I,length);
	NEW(mE->info,DEFAULT_INFO_LINE_LEN_SEQ,char);
	for(i=0; i<DEFAULT_INFO_LINE_LEN_SEQ; i++) mE->info[i] = E[0]->info[i];
	for(s=1,n=0; E[n] != NULL; n++){
	   for(i=1; i<= (Int4) E[n]->n; i++,s++){
		mE->S[s] = E[n]->S[i];
	   }
	}
	return mE;
}

e_type	MkEmptySeq(Int4 I, char *info, Int4 length)
{
	e_type	E=EmptySeq(I, length);
	NEW(E->info,DEFAULT_INFO_LINE_LEN_SEQ,char);
	for(Int4 i=0; i<DEFAULT_INFO_LINE_LEN_SEQ && info[i]!=0; i++) E->info[i] = info[i];
	return E;
}

e_type	EmptySeq(Int4 I, Int4 length)
/*  create and return an empty sequence of length */
{
	Int4		i;
	e_type	E = AllocSeq(I,length);

	for(i=0; i<=length; i++) E->S[i]=0;
        return E;
}

Int4    GetLongFastaInfo(char *DBS_NAME, Int4 max, Int4 **pcounts, Int4 **psize, a_type A)
/*****************************************************************
 Get information on the number of sequences, lengths (=size) and residue
    compositions for a fasta file
 *****************************************************************/
{
        Int4    i,length,*counts;
        Int4    *size;
        char    c=' ',last_c='\n';
        FILE    *fptr;

        if((fptr = fopen(DBS_NAME,"r")) == NULL) {
                fprintf(stderr,"Could not open file \"%s\"\n",DBS_NAME);
                print_error("File does not exist.");
        }
        if(setvbuf(fptr,NULL,_IOFBF, BUFSIZ*200))
                print_error("GetFastaInfo() buffer error");
        if(max > 0) NEW(size,max+2,Int4);
        NEW(counts,nAlpha(A)+3,Int4);
#if 0	// OLD
        while((c=fgetc(fptr))!=EOF){
		if(c=='>') break; 
	}
#endif
        while(c != EOF){
           if((c=fgetc(fptr)) == '>' && last_c=='\n') break;
           else if(!isspace(c)) print_error("File not in fasta format!");
           else last_c=c;
        }
        if(c==EOF) print_error("File not in fasta format!");

        for(i=1,length=0;c!=EOF;length=0,i++) {
                if(c=='>'){
		   	while((c=fgetc(fptr))!=EOF){ if(c=='\n') break; }
			while((c=fgetc(fptr)) == '+'){ // scan additional information.
			   ScanOverRFF(fptr); // will just skip over rff file...
			} // if(ungetc(c,fptr) == EOF) print_error("GetFastaInfo( )-> input error");
		}
                while(c!='>') {
                   if(isalpha(c)) { length++; counts[(AlphaCode(c,A))]++; }
		   // else if(c == '-'){ length++; counts[(AlphaCode(c,A))]++; }
		   else if(c == '-'){ length++; counts[(AlphaCode('X',A))]++; }
                   else if(!isspace(c)) {
                        fprintf(stderr,"sequence %d: %c undefined ",i,c);
                        fprintf(stderr," (ignored)\n");
                   } last_c=c;
                   if((c=fgetc(fptr))==EOF) break;
                }
		if(last_c != '\n' && c != EOF){
			 print_error("fasta input error: last_c != rtn && c != EOF");
                } else if(length==0){
			 print_error("fasta input error: length==0");
		}
                if(max > 0){
                  if(i > max){
                        fprintf(stderr,"%d sequences exceeds limit of %d\n",i,max);
                        print_error("reset maximum for GetFastaInfo( )");
                  } size[i] = length;
                }
        } i--;
        fclose(fptr);
        *pcounts = counts;
        if(max > 0) *psize = size;
        return i;
}

Int4	GetFastaInfo(char *DBS_NAME, Int4 max, Int4 **pcounts, 
	unsigned short **psize, a_type A)
/*****************************************************************
 Get information on the number of sequences, lengths (=size) and residue 
    compositions for a fasta file 
 *****************************************************************/
{
	Int4	i,j,length,*counts;
	unsigned short	*size;
	char	c=' ',last_c='\n';
        FILE    *fptr;

        if((fptr = fopen(DBS_NAME,"r")) == NULL) {
                fprintf(stderr,"Could not open file \"%s\"\n",DBS_NAME);
                print_error("File does not exist.");
        }
	if(setvbuf(fptr,NULL,_IOFBF, BUFSIZ*200))
		print_error("GetFastaInfo() buffer error");
	if(max > 0) NEW(size,max+2,unsigned short);
	NEW(counts,nAlpha(A)+3,Int4);
#if 0	// OLD
        while((c=fgetc(fptr))!=EOF){ if(c=='>') break; }
#endif
        while(c != EOF){
           if((c=fgetc(fptr)) == '>' && last_c=='\n') break;
           else if(!isspace(c)) print_error("File not in fasta format!");
           else last_c=c;
        }
	if(c==EOF) print_error("File not in fasta format!");

#if 0
	fprintf(stderr,"Determining database composition...");
#endif
        for(i=1,length=0;c!=EOF;length=0,i++) {
                if(c=='>'){
			while((c=fgetc(fptr))!=EOF){ if(c=='\n') break; }
			while((c=fgetc(fptr)) == '+'){ // scan additional information.
			   ScanOverRFF(fptr);
			} // if(ungetc(c,fptr) == EOF) print_error("GetFastaInfo( )-> input error");
		}
                while(c!='>') {
                   if(isalpha(c)) { length++; counts[(AlphaCode(c,A))]++; }
		   // else if(c == '-'){ length++; counts[(AlphaCode(c,A))]++; }
		   else if(c == '-'){ length++; counts[(AlphaCode('X',A))]++; }
                   else if(!isspace(c)) {
                        fprintf(stderr,"sequence %d: %c undefined ",i,c);
                        fprintf(stderr," (ignored)\n");
#if 0
                        fprintf(stderr," (replaced with null character: %c)\n",
				AlphaChar(UndefAlpha(A),A));
			length++; counts[UndefAlpha(A)]++;
#endif
                   } last_c=c;
                   if((c=fgetc(fptr))==EOF) break;
                }
		if(last_c != '\n' && c != EOF) {
			fprintf(stderr,"last_c=%c (sequence #%d)\n",last_c,i); 
		 	fprintf(stderr,"%c",c); 
			for(j=0; ((c=fgetc(fptr))!=EOF); j++){
				fprintf(stderr,"%c",c); 
				if(j > 200) break;
			} 
		 	fprintf(stderr,"\n"); 
			assert(last_c == '\n');
			print_error("GetFastaInfo(): fasta input error: last_c != rtn");
		} else if(length==0){
			print_error("GetFastaInfo(): fasta input error: length = 0");
		}
		if(max > 0){
                  if(i > max){
			fprintf(stderr,"%d sequences exceeds limit of %d\n",i,max);
			print_error("reset maximum for GetFastaInfo( )");
		  } size[i] = length;
		}
        } i--; 
        fclose(fptr);
#if 0
	fprintf(stderr," (%d sequences)\n",i); 
#endif
	*pcounts = counts; 
	if(max > 0) *psize = size;
	return i;
}

BooLean	AddCountsSubSeq(UInt4 *counts,Int4 start,Int4 end,e_type E)
// add residue counts within sequence region to counts array.
{
	Int4 j,r;

	assert(start > 0 && end <= LenSeq(E) && start <= end);
	for(j=start; j <= end; j++){ r=E->S[j]; counts[r]++; }
	return TRUE;
}

double	*FreqResSeq(e_type E, double *freq, a_type A)
{
	Int4 j,r;

	for(r = 0; r <= nAlpha(A); r++) freq[r] = 0.0;
	for(j=1; j <= (Int4)E->n; j++) { r=E->S[j]; freq[r]+=1.0; }
	for(r = 0; r <= nAlpha(A); r++) freq[r] /= (double) E->n;
	return freq;
}

void	NumResSeq(register e_type E,register UInt4 *num, a_type A)
/** return the number of residues in sequence. **/
{
	register Int4 	j;
	for(j=0; j <= nAlpha(A); j++) num[j]=0;
	for(j=1; j <= (Int4)E->n; j++) num[E->S[j]]++; 
}

Int4	*NumResSeq(e_type E,a_type A)
/** return residue frequencies **/
{
	Int4 	j,r,*num;

	NEW(num,nAlpha(A)+3,Int4);
	for(j=1; j <= (Int4)E->n; j++) { r=E->S[j]; num[r]++; }
	return num;
}

Int4	NonXResSeq(e_type E,a_type A)
/** return residue frequencies **/
{ Int4 N=0,j; for(j=1; j <= (Int4)E->n; j++) { if(E->S[j] > 0) N++; } return N; }

char	*copy_seq_info(e_type E)
{
	char	*info;
	assert(E->info != NULL);
	Int4 n=strlen(E->info); NEW(info,n+20,char);
	for(Int4 i=0;i<=n; i++) info[i]= E->info[i];
	return info;
}

void    AdjustCtermExtendSeq(Int4 diff, e_type E)
{ 
#if 0	// OLD: was core dumping...(on 4/21/03)
	assert((Int4)E->extend + diff >= 0); E->extend += diff; 
#else	// Not sure that this fixes problem, but give it a try...need to check!
	if((Int4)E->extend + diff >= 0){
		E->extend += diff; 
	} else {	// no need for an extension...
		E->extend = 0; // set to zero if negative extension...
	}
#endif
}

void	AdjustOffsetSeq(Int4 diff, e_type E)
{ assert((Int4)E->offset + diff >= 0); E->offset += diff; }

e_type	CopySeqPlus(e_type E, UInt4 p, UInt4 N)
// copy sequence E with N 'X' residues added at position p.
{
	e_type	E2;
	Int4	i,j,n;

	assert(E!=NULL);
	E2 = AllocSeq(E->I, E->n+N);
	E2->info=copy_seq_info(E);
	E2->offset = E->offset; E2->extend = E->extend;
	E2->kingdom=E->kingdom; 
	if(E->phylum) E2->phylum = AllocString(E->phylum);
	E2->genetic_code = E->genetic_code;
	E2->X=E2->S;
	for(i=0; i < p; ){ i++; E2->S[i]= E->S[i]; }
	for(j=i,n=0;n<N; n++){ j++; E2->S[j]= 0; }
	while(i<=(Int4)E->n){ i++; j++; E2->S[j]= E->S[i]; } 
	if(E->xed){ 
		AllocXSeq(E2);
		for(i=0; i < p; ){ i++; E2->X[i]= E->X[i]; }
		for(j=i,n=0;n<N; n++){ j++; E2->X[j]= 0; }
		while(i<=(Int4)E->n){i++; j++; E2->X[j]= E->X[i]; } 
	}
	return E2;
}

e_type	CopySeq(register e_type E)
{
	register e_type	E2;
	register Int4	i,len;

	assert(E!=NULL);
	E2 = AllocSeq(E->I, E->n);
	if(E->info != NULL){
		len = strlen(E->info);
		NEW(E2->info,len + 50,char);
		for(i=0;i<=len; i++) E2->info[i]= E->info[i];
	}
	E2->offset = E->offset; E2->extend = E->extend;
	E2->kingdom=E->kingdom; E2->genetic_code=E->genetic_code; 
	if(E->phylum) E2->phylum = AllocString(E->phylum);
	for(i=E->n;i>=1; i--) E2->S[i]= E->S[i];
	if(E->xed){
		AllocXSeq(E2); 
		for(i=1;i<=(Int4)E->n; i++) E2->X[i]= E->X[i];
	} else E2->X=E2->S;
	if(E->rff){ E2->rff = new rff_typ(E->rff); }
	return E2;
}

void	AddRFF2Seq(register e_type E, rff_typ *rff)
{ assert(rff); if(E->rff) delete E->rff; E->rff = rff; }

void	PutSeqRFF(FILE *fp,e_type E) { if(E->rff != 0) E->rff->Put(fp); }

e_type  CopySeq(register e_type E, char *rffstr)
{
	e_type E2= CopySeq(E);
	if(E2->rff) delete E2->rff;
	E2->rff= new rff_typ(rffstr);
	return E2;
}

/************************* Randomization Routines ***********************/
e_type	RandomSeq(Int4 length, Int4 I, double *freq, a_type A)
/* return a randomly generated new sequence of length with residue 
   frequency " freq" */
{
	Int4 	i,c;
	double	r;
	e_type	E;

	// h_type HG=Histogram("random numbers",0,100,4);;
	E = EmptySeq(I, length);
	for(i=1;i<=(Int4)E->n; i++) {
		r = (double) Random()/(double) RANDOM_MAX; /* 0 <= r <= 1 */
		// IncdHist(r, HG);
		for(E->S[i]=0,c=1; c <= nAlpha(A); c++){
			r=r-freq[c];
			if(r <= 0.0){ E->S[i]=c; break; }
		}
	}
	// PutHist(stdout,60,HG); NilHist(HG);
	return E;
}

char	RandomResSeq(register e_type E)
/* return a randomly obtained residue in sequence E */
{
	register Int4	i;
	register double	r;

	do {
		r = (double) Random()/(double) RANDOM_MAX; /* 0 <= r <= 1 */
		i = (Int4) (r*E->n + 1); 
	} while(i < 1 && i > (Int4) E->n);
	return E->S[i];
}

unsigned char *ShuffleSeqArray(Int4 length, unsigned char *seq)
// randomly permute the sequence seq from position 0 to length-1;
// this uses a heap & does not move 0's. 
{
	unsigned char	*S;
	Int4 	i;
	dh_type	H;

	NEW(S,length+3,unsigned char);
	H = dheap(length+1,4);
	for(i=1; i <= length; i++) insrtHeap(i,((keytyp)Random()),H);
	for(i=1;i <= length; i++) S[i] = seq[delminHeap(H)];
	Nildheap(H); 
	return S;
}

e_type	ReverseSeq(e_type E) // reverse the sequence E
{
	unsigned char	r;
	Int4 	j,i,os,ex;

	if(E->xed) UnXSeq(E);
	for(j=E->n,i=1; i < j; i++,j--){
		r = E->S[j]; E->S[j]=E->S[i]; E->S[i]=r;
	}
	// if(E->info!=NULL) free(E->info); E->info = NULL; 
	os=E->offset; ex=E->extend;
	E->offset=ex; E->extend=os;
	if(E->rff){ delete E->rff; E->rff=0; }
	return E;
}

e_type	ShuffleSubSeq(Int4 start, Int4 end, e_type E)
// randomly permute the sequence E from positions start to end;
// this uses a heap & does not move 'X's. 
{
	unsigned char	*S;
	Int4 	i;
	dh_type	H;

	if(start >= end) return E;
	if(E->xed) UnXSeq(E);
	NEW(S,E->n+3,unsigned char);
	H = dheap(E->n+1,4);
	start = MAXIMUM(Int4,start,1);
	end = MINIMUM(Int4,end,E->n);
	for(i=start; i <= end; i++){
		S[i]=E->S[i];
		if(S[i] != 0){ insrtHeap(i,((keytyp)Random()),H); }
	}
	for(i=start;i <= end; i++){
		if(E->S[i] != 0){ E->S[i] = S[delminHeap(H)]; }
	}
	Nildheap(H); free(S); 
	if(E->info!=NULL)free(E->info); E->info = NULL; 
	E->offset=E->extend=0;
	if(E->rff){ delete E->rff; E->rff=0; }
	return E;
}

e_type	ShuffleSeq(e_type E){ return ShuffleSubSeq(1,E->n,E); }

Int4	MaxSegSeq(e_type E)
/** return the length of the longest high complexity segment in E **/
{
	Int4 	j,len,max;

	for(max=len=0,j=1; (Int4)j <= (Int4)E->n; j++) {
	   if(E->X[j] > 0) len++;
	   else { if(len > max) max=len; len=0; }
	}
	if(len > max) max=len;
	return max;
}

BooLean	EST_Seq(e_type E)
{
	char *result=strstr(E->info,"_EST");
	if(result==0) return FALSE;
	else {
		char *space=strchr(E->info,' ');
		if(space > result) return TRUE; else return FALSE;
	}
}

BooLean	PdbSeq(e_type E)
{
	char *result=strstr(E->info,"|pdb|");
	if(result==0) return FALSE;
	else {
		char *space=strchr(E->info,' ');
		if(space > result) return TRUE; else return FALSE;
	}
}

BooLean	SwissSeq(e_type E)
{	
#if 0
	if(strstr(E->info,"|sp|") == NULL) return FALSE; else return TRUE; 
#else
	char *result=strstr(E->info,"|sp|");
	if(result==0) return FALSE;
	else {
		char *space=strchr(E->info,' ');
		if(space > result) return TRUE; else return FALSE;
	}
#endif
}

Int4    NumXSeq(e_type E)
{
	Int4 n,i;
	for(n=0,i=1; i <= LenSeq(E); i++) if(ResSeq(i,E)==0) n++; 
	return n;
}

/************************* Print Routines ***********************/
void	PutSeqInfo2(FILE *fptr,e_type E)
{
	if(E->info !=NULL){
		fprintf(fptr,"%s\n", E->info); 
	} else {
		fprintf(fptr,"random%d \n",(Int4)E->I);
	}
}

void    ChangeInfoSeq(char *new_info, e_type E)
{
	if(E->info !=NULL) free(E->info);
	assert(isprint(new_info[0]));
	E->info = AllocString(new_info);
}

void    StrSeqID(char *str, Int4 n, e_type E)
{
	Int4	i;
	char	c,str2[40];

	if(E->info !=NULL){
	     for(i=0; (c=E->info[i])!=0; i++){
		if(i==n || isspace(c)){ str[i]=0; return; }
		else { str[i]=c; }
	     } str[i]=0; return;
	} else { sprintf(str2,"random%d ",(Int4)E->I); strncat(str,str2,n); }
}

void    StrSeqDescript(char *str,Int4 n, e_type E)
// get the sequence description.
{
	Int4 	i,j;
	char	c;

	if(E->info !=NULL){
	     // i=0; while(!isspace(E->info[i])) i++;	// modified for short descriptions...
	     i=0; while(!isspace(E->info[i])){ 
		if(!isprint(E->info[i])){ str[0]=0; return; }
		i++; 
	     }
	     if(E->info[i]==0){ str[0]=0; return; }  // no description given...
	     for(j=0,i++; (c=E->info[i])!=0; i++,j++){
		if(j==n || c=='\n'){ str[j]=0; return; }
	    	else { str[j]=c; }
	     } str[j]=0; return;
	} else sprintf(str," ");
}

void	StrSeqInfo(char *str,e_type E)
{	
	Int4	i;
	BooLean	flag; 
	char	c,str2[MAX_SEQ_DEFLINE_LENG + 100];

	if(E->info !=NULL){
	   if(E->offset || E->extend || E->phylum){
	     str[0]=0;
	     for(flag=TRUE,i=0; (c=E->info[i])!=0; i++){
	       if(flag && isspace(c)) {
		 if(E->offset || E->extend){
		  if(E->phylum){
		    if(E->genetic_code){
			 sprintf(str2," {|%u(%u)|<%s(%c%d)>}", E->offset,E->extend,
						E->phylum,E->kingdom,E->genetic_code);
		    } else {
			 sprintf(str2," {|%u(%u)|<%s(%c)>}",
					E->offset,E->extend,E->phylum,E->kingdom);
		    }
		  } else sprintf(str2," {|%u(%u)|}", E->offset,E->extend);
	         } else if(E->phylum){
		   if(E->genetic_code){
			sprintf(str2," {<%s(%c%d)>}",E->phylum,E->kingdom,E->genetic_code);
		   } else {
			sprintf(str2," {<%s(%c)>}",E->phylum,E->kingdom);
		   }
		 }
		 strcat(str,str2);
		 flag = FALSE;
	       } else { sprintf(str2,"%c",c); strcat(str,str2); }
	     }
	   } else { sprintf(str,"%s", E->info); }
	} else {
	   if(E->offset || E->extend){
		sprintf(str,"random%d {|%d(%d)|}",
		   (Int4)E->I,E->offset,E->extend);
	   } else sprintf(str,"random%d ",(Int4)E->I);
	}
}

void	PutSeqInfo(FILE *fptr,e_type E) { PutSubSeqInfo(fptr,0,E->n,E); }

void	PutSubSeqInfo(FILE *fptr,Int4 start, Int4 end, e_type E)
{	
	Int4	i,j;
	BooLean	flag; 
	char	c,special_info[200];
	UInt4	o,e;

	start = MAXIMUM(Int4,start,1);
	end = MINIMUM(Int4,end,E->n);
	o = (start-1 + E->offset);
	e = (E->n-end + E->extend);
	if(E->info !=NULL){
	   if(o || e || E->phylum){
	     for(flag=TRUE,j=i=0; ((c=E->info[i])!=0) && j < DEFAULT_INFO_LINE_LEN_SEQ; j++,i++){
	       if(flag && isspace(c)) { //****** stop here and output info ******//
#if 0	// OLD
	         if(o || e){
		  if(E->phylum) {
		    if(E->genetic_code){
			fprintf(fptr," {|%d(%d)|<%s(%c%d)>}",o,e,
				E->phylum,E->kingdom,E->genetic_code);
		    } else {
			fprintf(fptr," {|%d(%d)|<%s(%c)>}",
				o,e,E->phylum,E->kingdom);
		    }
		  } else fprintf(fptr," {|%d(%d)|}", o,e);
	   	 } else if(E->phylum){
		    if(E->genetic_code){
			fprintf(fptr," {<%s(%c%d)>}",E->phylum,E->kingdom,E->genetic_code);
		    } else {
			fprintf(fptr," {<%s(%c)>}",E->phylum,E->kingdom);
		    }
		 } flag = FALSE;
#else 	// NEW
	         if(o || e){
		  if(E->phylum) {
		    if(E->genetic_code){
			sprintf(special_info," {|%u(%u)|<%s(%c%d)>}",o,e,
				E->phylum,E->kingdom,E->genetic_code);
		    } else {
			sprintf(special_info," {|%u(%u)|<%s(%c)>}",
				o,e,E->phylum,E->kingdom);
		    }
		  } else sprintf(special_info," {|%u(%u)|}",o,e);
	   	 } else if(E->phylum){
		    if(E->genetic_code){
			sprintf(special_info," {<%s(%c%d)>}",
				E->phylum,E->kingdom,E->genetic_code);
		    } else {
			sprintf(special_info," {<%s(%c)>}",E->phylum,E->kingdom);
		    }
		 } flag = FALSE;
		 j += strlen(special_info);
		 fprintf(fptr," %s",special_info); flag = FALSE;
#endif
	       } else fprintf(fptr,"%c",c);
	     } fprintf(fptr,"\n");
	   } else { fprintf(fptr,"%s\n", E->info); }
	} else {
	   if(o || e){
		fprintf(fptr,"random%d {|%d(%d)|}\n",(Int4)E->I,o,e);
	   } else {
		fprintf(fptr,"random%d \n",(Int4)E->I);
	   }
	}
}

void	ShortenSeqInfo(e_type E)
// eliminate path information from seq info (if present). 
{ 
	Int4	k,i,start;
	char	c,*s;
	char	str[MAX_SEQ_DEFLINE_LENG+10];

	if(E->info !=NULL) {
	   str[MAX_SEQ_DEFLINE_LENG]=0;
	   for(start=i=0; E->info[i] && !isspace(E->info[i]); i++){
		if(E->info[i] == '/') start = i+1;
		if(isspace(E->info[start])) print_error("SeqInfo error in ShortenSeqInfo()");
	   } sprintf(str,"%s",E->info + start);
		   if(str[MAX_SEQ_DEFLINE_LENG] != 0) print_error("ShortenSeqInfo() error");
	   free(E->info); E->info=AllocString(str);
       	}
}

void	PutShortSeqID(FILE *fptr,e_type E)
{ 
	Int4	k,i;
	char	c,*s;
	char	str[MAX_SEQ_DEFLINE_LENG+10];

	if(E->info !=NULL) {
	   for(i=0; E->info[i] && !isspace(E->info[i]); i++){
		str[i]=E->info[i];
	   } str[i]=0;
	   if(i > 0 && str[i-1]=='|') str[i-1]=0;
	   if(strstr(str,"|") == NULL) fprintf(fptr,"%s",str);
	   else if((s=strstr(str,"pdb|")) != NULL){
	     for(s+=4; (!isspace(*s) && *s!=0); s++) fprintf(fptr,"%c",*s);
	   } else {
	     c = str[0];
	     for(k=0,i=0; c && !isspace(c); k++,c=str[k]){
		if(c == '|') i = k+1;
	     } 
	     if(isspace(str[i])){ // '|'s at end of line
		for(i--; i >= 0 && str[i] == '|'; ){ k=i; i--; }
		//  k = i-1; i-=2; // OLD...
		while(i >= 0 && str[i] != '|'){ i--; }
		i++;
	     }
	     for( ; i < k; i++){ fprintf(fptr,"%c", str[i]); }
	   }
       	} else fprintf(fptr,"random%d ",(Int4)E->I);
}

void	PutSeqID(FILE *fptr,e_type E)
{ 
	Int4	k;

	if(E->info !=NULL) {
	   for(k=0; !isspace(E->info[k]); k++){
		if(E->info[k]==0) break;
		fprintf(fptr,"%c", E->info[k]);
	   }
       	} else fprintf(fptr,"random%d ",(Int4)E->I);
}

void	PutSeqID2(FILE *fptr,e_type E, Int4 len)
{
	Int4	k;

	if(E->info !=NULL) {
	   for(k=0; TRUE; k++){
		if(E->info[k]==0 || isspace(E->info[k])) { 
			for( ;k < len; k++) fprintf(fptr," "); 
			break;
		} else if(k >= len) break;
		fprintf(fptr,"%c", E->info[k]);
	   }
       	} else fprintf(fptr,"random%-4d    ",(Int4)E->I);
       	fprintf(fptr," ");
}

void    PutXSeq(FILE *fptr,e_type E,a_type A)
{
	Int4 	j;

	fprintf(fptr,">"); PutSeqInfo(fptr,E);
	if(E->rff){ E->rff->Put(fptr); }
	for(j=1; (Int4)j <= (Int4)E->n; j++) {
	   if(E->X[j]==0) fprintf(fptr,"x");
	   else fprintf(fptr,"%c",AlphaChar(E->X[j],A));
	   if(j%10 == 0) fprintf(fptr," ");
	   if(j%50 == 0) fprintf(fptr,"\n");
	}
	fprintf(fptr,"\n\n");
}

void    PutSeqMtfMask(FILE *fptr,Int4 M, Int4 *len, Int4 *p0, Int4 *p1,
	e_type E1,a_type A)
// Print out masked sequence for threading through a PDB structure 
{
	Int4 	i,j,k,m;

	fprintf(fptr,">"); PutSeqInfo(fptr,E1);
	for(i=m=1; m <= M; m++){
	  while(i < p0[m]){ fprintf(fptr,"x"); i++; }
	  for(j=p1[m], k=1; k <= len[m]; k++,j++) {
	     fprintf(fptr,"%c",AlphaChar(E1->S[j],A));
	  }
	  i+=len[m];
	  fprintf(fptr,"\n");
	}
	while(i < p0[m]){ fprintf(fptr,"x"); i++; }
	fprintf(fptr,"\n\n");
}

void    PutSeq(FILE *fptr,e_type E,a_type A)
{ PutSeq(fptr,0,E,A); }

void    PutSeq(FILE *fptr,Int4 subID, e_type E,a_type A)
{ PutSeq(fptr,subID, E,A,TRUE); }

void    PutSeq(FILE *fptr,Int4 subID, e_type E,a_type A,BooLean Rtns)
{
	Int4 	i,j;

	// added phylogenetic info...
	char	c;
	Int4 o = (E->offset);
	Int4 e = (E->extend);
	BooLean	flag;
	fprintf(fptr,">"); 
	if(o != 0 || e != 0){
	   for(flag=TRUE,i=0; (c=E->info[i])!=0; i++){
	    if(flag && isspace(c)) {
		if(subID > 0){ fprintf(fptr,"%d",subID); }
		if(E->phylum) {
		    if(E->genetic_code){
			fprintf(fptr," {|%d(%d)|<%s(%c%d)>}", o,e,
					E->phylum,E->kingdom,E->genetic_code);
		    } else {
			fprintf(fptr," {|%d(%d)|<%s(%c)>}", o,e,E->phylum,E->kingdom);
		    }
		} else fprintf(fptr," {|%d(%d)|}", o,e);
		flag = FALSE;
	    } else fprintf(fptr,"%c",c);
	   }
	   if(flag){
		fprintf(fptr," {|%d(%d)|", o,e);
		if(E->phylum) fprintf(fptr,"<%s(%c)>}",E->phylum,E->kingdom);
		else fprintf(fptr,"}");
	   } fprintf(fptr,"\n");
	} else if(E->phylum) {
	  if(E->info == NULL) {
       		fprintf(fptr,"random%d \n",(Int4)E->I);
	  } else {
	   for(flag=TRUE,i=0; (c=E->info[i])!=0; i++){
	    if(flag && isspace(c)) {
		if(subID > 0){ fprintf(fptr,"%d",subID); }
		if(E->phylum) {
		    if(E->genetic_code){
			fprintf(fptr," {<%s(%c%d)>}",
					E->phylum,E->kingdom,E->genetic_code);
		    } else {
			fprintf(fptr," {<%s(%c)>}",E->phylum,E->kingdom);
		    }
		} else fprintf(fptr," ");
		flag = FALSE;
	    } else fprintf(fptr,"%c",c);
	   } if(flag && E->phylum) fprintf(fptr," {<%s(%c)>}",E->phylum,E->kingdom);
	   fprintf(fptr,"\n");
	 }
	} else {
#if 1	// allow unique ids for est_other_txaa for indexfa
	  if(E->info == NULL) { fprintf(fptr,"random%d \n",(Int4)E->I); }
	  else {
	   for(flag=TRUE,i=0; (c=E->info[i])!=0; i++){
	    if(flag && isspace(c)) {
		if(subID > 0){ fprintf(fptr,"%d",subID); }
		fprintf(fptr," ");
		flag = FALSE;
	    } else fprintf(fptr,"%c",c);
	   } fprintf(fptr,"\n");
	 }
#else
		PutSeqInfo(fptr,E);
#endif
	}
	if(E->rff){ E->rff->Put(fptr); }
	for(j=1; (Int4)j <= (Int4)E->n; j++) {
	   fprintf(fptr,"%c",AlphaChar(E->S[j],A));
	   if(Rtns){
	   	if(j%10 == 0) fprintf(fptr," ");
	   	if(j%50 == 0) fprintf(fptr,"\n");
	   }
	}
	fprintf(fptr,"\n\n");
}

void    MaskSeq(Int4 start, Int4 end, e_type E)
/** set the region from start to end to "null" residues (0) **/
{
	Int4	i;
	unsigned char	*s;
	
	if(start < 1 || start > end || end > (Int4) LenSeq(E)){
		seq_error("MaskSeq( ) - out of range.");
	} 
	if(!E->xed) AddXArraySeq(E);
	s=E->X;
	for(i=start; i <= end; i++) s[i] = 0;
}

void	AddXArraySeq(register e_type E)
/** add an x-able array to E (xed=null) **/
{
	register Int4	i;
	
	if(E->X == E->S){
	   AllocXSeq(E);
	   for(E->X[0]=0,i=1; i<=(Int4)E->n; i++) E->X[i]=E->S[i];
	}
}

char    *get_info_seq(e_type E, char *info)
/*** retrieves pertinent information from fasta defline and converts defline
    into useable format; also retrieves subsequence information.  ***/
{
	char	*id,*p,*p0;
	Int4	i,n,code;
	char	kingdom,phylum[100];

	if(info == NULL) return NULL;
// fprintf(stderr,"DEBUG1\n");
	n = MINIMUM(UInt4, strlen(info), DEFAULT_INFO_LINE_LEN_SEQ);
	NEW(id,n+10,char);
	E->kingdom=0; E->phylum=0; E->genetic_code=0;
	if((p=strstr(info," {|")) == NULL){	// no offset info.
	  E->offset=0; E->extend=0;
	  if((p=strstr(info," {<")) == NULL){
	     for(i=0, p=info; isprint(*p) && i < n; i++,p++) { id[i] = *p; }
	     id[i] = 0; return id;
	  } else {	// Taxonomic information probable...
	    if(sscanf(p," {<%[^(](%c%d)>}",phylum,&kingdom,&code) != 3) {
	      p=strstr(info," {<"); 
	      if(sscanf(p," {<%[^(](%c)>}",phylum,&kingdom) != 2) {
		// No tax info either, though weird defline...
	        for(i=0, p=info; isprint(*p) && i < n; i++,p++) { id[i] = *p; }
	        id[i] = 0; return id;
	      } else { E->kingdom=kingdom; E->phylum=AllocString(phylum); }
	    } else {	// Taxonomic information...
		E->genetic_code=(unsigned char)code;
		E->kingdom=kingdom; E->phylum=AllocString(phylum);
	    }
	  }
	} else {	// Offset information probable.
	  if(sscanf(p," {|%u(%u)|<%[^(](%c%d)>}",
			&E->offset,&E->extend,phylum,&kingdom,&code) != 5) {
	    p=strstr(info," {|"); 
	    if(sscanf(p," {|%u(%u)|<%[^(](%c)>}",&E->offset,&E->extend,phylum,&kingdom) != 4) {
	      p=strstr(info," {|");
	      if(sscanf(p," {|%u(%u)|}",&E->offset,&E->extend) != 2) {
	        E->offset=0; E->extend=0;
	        p=strstr(info," {|");
	        if(sscanf(p," {|%u|}",&E->offset) != 1) { 
		  E->offset = 0;
	   	  for(i=0, p=info; isprint(*p) && i < n; i++,p++) { id[i] = *p; }
	   	  id[i] = 0; return id;
	        }
	      }
	    } else { E->kingdom=kingdom; E->phylum=AllocString(phylum); }
	  } else {
		E->genetic_code=(unsigned char)code; 
		E->kingdom=kingdom; E->phylum=AllocString(phylum); 
	  }
	}
	for(i=0,p=info; *p != ' '; p++,i++){ id[i] = *p; } id[i] = ' '; i++;
	p0=p;
	if((p=strstr(p0,"|}")) == NULL){ p=strstr(p0,">}"); }
	for(p+=2; isprint(*p) && i < n; i++,p++) { id[i] = *p; }
	if(i > 0 && id[i-1] != ' '){ id[i]=' '; i++; }
	id[i] = 0; i++;
	return id;
}

e_type	MkSubSeq(Int4 start, Int4 end, e_type E)
/*** make a subsequence of E starting at s and ending at e ***/
{
	Int4	len,i,k;
	char	str[MAX_SEQ_DEFLINE_LENG+10],str2[100];
	UInt4	o,e;
	e_type	sE;

	if(start < 1 || end > E->n){
		fprintf(stderr,"MkSubSeq(): start=%d; end = %d\n",start,end);
		fprintf(stderr,"SubSeq out of range\n");
		assert(!(start < 1 || end > E->n));
	}
	o = (start-1 + E->offset);
	e = (E->n-end + E->extend);
	if(E->info == NULL) {
       		sprintf(str,"random%d {|%u(%u)|}\n",
			(Int4)E->I,o,e);
	} else {
	  for(i=0; isgraph(E->info[i]); i++){ 	// while printing char except space.
		str[i] = E->info[i];
		if(i > DEFAULT_INFO_LINE_LEN_SEQ) break;
	  }
	  str[i] = 0; i++;
	  sprintf(str2," {|%u(%u)|}", o,e);
	  strcat(str,str2);
	  k = strlen(str);
	  if(E->info[i] != 0) strncat(str,E->info + i, DEFAULT_INFO_LINE_LEN_SEQ- k);
	} 
	len = end - start + 1;
	sE=MkSeq(str, len, (E->S + start - 1));
	sE->kingdom = E->kingdom; 
	if(E->phylum) sE->phylum = AllocString(E->phylum); else sE->phylum = 0;
	if(E->rff){ sE->rff = new rff_typ(E->rff,start,end); }
	return sE;
}

void    PutSubSeq(FILE *fptr, Int4 start, Int4 end, e_type E, a_type A)
// put the subseq in fasta format...
{
	Int4	len,i;
	char	c;
	BooLean	flag;
	UInt4	o,e;

	start = MAXIMUM(Int4,start,1);
	end = MINIMUM(Int4,end,E->n);
	o = (start-1 + E->offset);
	e = (E->n-end + E->extend);
	fprintf(fptr,">"); 
	if(E->info == NULL) {
	   if(o == 0 && e == 0) fprintf(fptr,"random%d \n",(Int4)E->I);
	   else fprintf(fptr,"random%d {|%d(%d)|}\n",(Int4)E->I,o,e);
	} else {
	  if(o == 0 && e == 0) {
	   if(E->phylum){
	     for(flag=TRUE,i=0; (c=E->info[i])!=0; i++){
	      if(flag && isspace(c)) {
		if(E->genetic_code){
		  fprintf(fptr," {<%s(%c%d)>}",E->phylum,E->kingdom,E->genetic_code);
		} else fprintf(fptr," {<%s(%c)>}",E->phylum,E->kingdom);
		flag = FALSE;
	      } else fprintf(fptr,"%c",c);
	     } fprintf(fptr,"\n");
	   } else fprintf(fptr,"%s\n",E->info);
	  } else {
	   for(flag=TRUE,i=0; (c=E->info[i])!=0; i++){
	    if(flag && isspace(c)) {
		if(E->phylum) {
		    if(E->genetic_code){
			fprintf(fptr," {|%d(%d)|<%s(%c%d)>}", o,e,
					E->phylum,E->kingdom,E->genetic_code);
		    } else {
			fprintf(fptr," {|%d(%d)|<%s(%c)>}", o,e,E->phylum,E->kingdom);
		    }
		} else fprintf(fptr," {|%d(%d)|}", o,e);
		flag = FALSE;
	    } else fprintf(fptr,"%c",c);
	   }
	   fprintf(fptr,"\n");
	  }
	}
	if(E->rff){ E->rff->Put(fptr,start,end); }
	for(len=0,i=start; i <= end; i++){
	   if(i > 0 && i <= (Int4) E->n) {
	   	fprintf(fptr,"%c", AlphaChar(E->S[i],A)); len++;
	   	if(len%70==0) fprintf(fptr,"\n");
	   }
	}
	fprintf(fptr,"\n\n");
}

void    PutSeqRegion(FILE *fptr,Int4 start, Int4 length, e_type E, a_type A)
{ PutSeqRegion2(fptr,start, length, E, 10, A); }

void    PutSeqRegionFormatMSA(FILE *fptr,Int4 start, Int4 length, e_type E, 
	double p, a_type A)
{
	Int4	i,offset,r,end;
	char	c;

	offset = OffSetSeq(E);
	end = start + length -1;
        fprintf(fptr,"     ");
        for(i=start; i <= end; i++){
            if(i < 1 || i > (Int4) LenSeq(E)) fprintf(fptr," ");
            else{
               r = ResSeq(i,E);
               c = AlphaChar(r,A);
               fprintf(fptr,"%c", c);
            }
        }
        fprintf(fptr," %d-%d ",start+offset,end+offset);
        PutShortSeqID(fptr,E);
        fprintf(fptr," %.2f",p);
        fprintf(fptr,"\n");
}

void    PutSeqRegion2(FILE *fptr,Int4 start, Int4 length, e_type E, 
	Int4 flank, a_type A)
{
	Int4	e,end,i;

	i = MAXIMUM(Int4,1,(start+E->offset));
	fprintf(fptr,"%4d  ",i);
	e = start + length - 1;
	if(flank > 0){
	  end = e + flank;
	  for(i=start-flank; i <= end; i++){
		if(i < 1 || i > (Int4) E->n) fprintf(fptr," ");
		else if(i == e) {
			fprintf(fptr,"%c ", AlphaChar(E->S[i],A));
		} else if(i == start) {
			fprintf(fptr," %c", AlphaChar(E->S[i],A));
		} else {
			fprintf(fptr,"%c", AlphaChar(E->S[i],A));
		}
	  }
	} else {
	  for(i=start; i <= e; i++){
		if(i < 1 || i > (Int4) E->n) fprintf(fptr," ");
		else fprintf(fptr,"%c", AlphaChar(E->S[i],A));
	  }
	}
	e = MINIMUM(Int4,e,LenSeq(E));
	fprintf(fptr," %4d",e+E->offset);
}

Int4	get_diagonal_ends_seq(unsigned char *seq1, unsigned char *seq2, 
	char **R, Int4 n, Int4 *begin, Int4 *end)
{
	register Int4 min=INT4_MAX,sum=0,score=-9999;
	Int4	min_n=0;

	if(n <= 0) { *begin = n+1; *end = n; return 0; }
	score = min = sum = R[seq1[n]][seq2[n]];
	*begin = *end = min_n = n;
	n--;
	while(n > 0){
        	if(min > (sum += R[seq1[n]][seq2[n]])){
			min = sum; min_n = n-1;
		} 
		if(score < (sum-min)){
			score = (sum-min);
			*end = min_n; *begin = n;
		} n--;
        };
	/*** printf("\tEND\n\n"); /****/
	return score;
}

Int4	FindMaxWordSeq(Int4 q_start, Int4 s_start, Int4 N, e_type qE, 
	e_type sE, Int4 *word_score, a_type A)
// Find the first residue of the maximum word in ungapped alignment 
// of len N between sequence qE with sE starting at q_start & s_start.
// returns the offset from q_start & s_start for starting residue.
// word length = 3;
{
	Int4	os,s0,s1,s2,offset,max;
	register Int4 	i,s,score;
	unsigned char   *qsq=SeqPtr(qE)+q_start,*ssq=SeqPtr(sE)+s_start;

	assert(N >= 3);
	assert((q_start+N-1) <= LenSeq(qE));
	assert((s_start+N-1) <= LenSeq(sE));
	
	s0 = valAlphaR(qsq[0],ssq[0],A);
	s1 = valAlphaR(qsq[1],ssq[1],A);
	s2 = valAlphaR(qsq[2],ssq[2],A);
	max = score = s0+s1+s2; offset = 0;
	for(os =1, i = 3; i < N; i++, os++){
	   s = valAlphaR(qsq[i],ssq[i],A);
	   score = score + s - s0;
	   if(score > max){ max=score; offset = os; }
	   s0 = s1; s1 = s2; s2 = s;
	} *word_score = max;
	return offset;
}

Int4	CenterSeqHSP(Int4 offset, e_type E1, e_type E2, a_type A)
/* return the center of the high scoring pair with offset **/
{
	Int4	v,n,n1,n2,e,start,end,center;
	unsigned char	*seq1,*seq2;

	n1 = LenSeq(E1); n2 = LenSeq(E2);
	v = offset;
	if(v < (1-n2) || v >= n1) seq_error("offset out of range");
	seq1 = XSeqPtr(E1); seq2 = XSeqPtr(E2);
	if(v > 0){		/** case 1: offset **/
	   n = MINIMUM(Int4,n1-v,n2);
	   get_diagonal_ends_seq(seq1+v,seq2,AlphaR(A),n,&start,&e);
	   start = v+start; end = v+e-1;
	   
	} else {
	   n = MINIMUM(Int4,n2+v,n1);
	   get_diagonal_ends_seq(seq2-v,seq1,AlphaR(A),n,&start,&e);
	   end = e - 1;
	}
	n = (end - start + 1)/2;
	center = start + n;
	return center;
}

Int4	PutDiagonalSeq(FILE *fptr, Int4 offset, e_type E1, e_type E2, a_type A)
/*******************************************************************
 print the MSP for diagonal at offset between seq1 and seq2.

 *******************************************************************/
{
	Int4	v,i,j,n,n1,n2,score=0,b,e,start,end,w=40;
	unsigned char	*seq1,*seq2;

/** ALPHABETA=A; /*****/
	n1 = LenSeq(E1); n2 = LenSeq(E2);
	v = offset;
	if(v < (1-n2) || v >= n1) seq_error("offset out of range");
	seq1 = XSeqPtr(E1); seq2 = XSeqPtr(E2);
	if(v > 0){
	   n = MINIMUM(Int4,n1-v,n2);
	   score=get_diagonal_ends_seq(seq1+v,seq2,AlphaR(A),n,&start,&e);
	   fprintf(fptr,"\n\n");
	   for(b=start; b <= e; b+=w){
		end = MINIMUM(Int4,e,b+w-1);
		fprintf(fptr,"%4d ",v+b+OffSetSeq(E1));
		for(i=v+b,j=b; j <= end; i++,j++)
			fprintf(fptr,"%c",AlphaChar(seq1[i],A));
	  	fprintf(fptr," %4d ",v+end+OffSetSeq(E1));
		PutSeqID(fptr,E1);
	  	fprintf(fptr,"\n     ");
		for(i=v+b,j=b; j <= end; i++,j++)
			if(seq1[i]==seq2[j]) fprintf(fptr,":");
			else if(valAlphaR(seq1[i],seq2[j],A) >= 0)
							fprintf(fptr,".");
			else fprintf(fptr," ");
		fprintf(fptr,"\n%4d ",b+OffSetSeq(E2));
		for(j=b; j <= end; j++)
			fprintf(fptr,"%c", AlphaChar(seq2[j],A));
		fprintf(fptr," %4d ",end+OffSetSeq(E2));
		PutSeqID(fptr,E2);
		fprintf(fptr,"\n\n");
	   }
	   fprintf(fptr,"   score = %d\n\n",score);
	} else {
	   n = MINIMUM(Int4,n2+v,n1);
	   score=get_diagonal_ends_seq(seq2-v,seq1,AlphaR(A),n,&start,&e);
	   fprintf(fptr,"\n\n");
	   for(b=start; b <= e; b+=w){
		end = MINIMUM(Int4,e,b+w-1);
		fprintf(fptr,"%4d ",b+OffSetSeq(E1));
		for(j=b; j <= end; j++)
			fprintf(fptr,"%c",AlphaChar(seq1[j],A));
		fprintf(fptr," %4d ",end+OffSetSeq(E1));
		PutSeqID(fptr,E1);
	  	fprintf(fptr,"\n     ");
		for(i=b-v,j=b; j <= end; i++,j++)
			if(seq2[i]==seq1[j]) fprintf(fptr,":");
			else if(valAlphaR(seq2[i],seq1[j],A) >= 0)
							fprintf(fptr,".");
			else fprintf(fptr," ");
		fprintf(fptr,"\n%4d ",b-v+OffSetSeq(E2));
		for(i=b-v,j=b; j <= end; i++,j++)
			fprintf(fptr,"%c", AlphaChar(seq2[i],A));
	  	fprintf(fptr," %4d ",end-v+OffSetSeq(E2));
		PutSeqID(fptr,E2);
		fprintf(fptr,"\n\n");
	   }
	   fprintf(fptr,"   score = %d\n\n",score);
	}
	return score;
}

e_type  ReadGuideSeqFA(char *infile, Int4 I, BooLean **ignore, a_type A)
{
        char    c;
        Int4    length,s;
	FILE	*fptr;
	BooLean *ig;
	e_type	E;

	E = ReadSeqFA(infile, I, A);
        fptr = open_file(infile,"","r");
        length=LenSeq(E);
	NEW(ig, length+2, BooLean);
        while((c=fgetc(fptr))!=EOF){ if(c=='>') break; }
        if(c=='>') while((c=fgetc(fptr))!=EOF){ if(c=='\n') break; }
        for(s=0; c!='>'; ) {
             if(isalpha(c)) {
                s++;  
		if(islower(c)) ig[s] = TRUE;
		else ig[s] = FALSE;
             } else if(!isspace(c)) seq_error("ReadGuideSeqFA( ) fatal error");
             if((c=fgetc(fptr))==EOF) break;
        }
        fclose(fptr);
	*ignore = ig;
	return E;
}

sst_typ *SST_FromSeq(e_type keyE)
{
        Int4    j,r,Length=LenSeq(keyE);
        sst_typ *xsst; NEW(xsst,Length +4, sst_typ);
        for(j=1; j <= Length; j++){ r=ResSeq(j,keyE); xsst[j] = SsetLet(r); }
        return xsst;
}

e_type	*ReadSeqFileFA(char *infile, a_type A, Int4 max_in_seq)
{
	Int4	number,*counts,J;
	unsigned short  *nsize;
	e_type	*E;
	FILE	*fp;
	
	number = GetFastaInfo(infile, max_in_seq, &counts, &nsize, A);
	NEW(E,number+3,e_type);
	fp = open_file(infile,"","r");
	for(J=1; J <= number; J++) E[J] =ReadSeq(fp,J,nsize[J],A);
	free(counts); free(nsize); fclose(fp);
	return E;
}

e_type  ReadSeqFA(char *infile, Int4 I, a_type A)
/** read a single sequence in from a fasta file ***/
{
        char    c=' ',last_c='\n';
        Int4    length;
	FILE	*fptr;
	e_type	E;

        fptr = open_file(infile,"","r");
        length=0;
        // while((c=fgetc(fptr))!=EOF){ if(c=='>') break; } // old...
        while(c != EOF){
           if((c=fgetc(fptr)) == '>' && last_c=='\n') break;
           else if(!isspace(c)) print_error("ReadSeqFA( ) fasta input error");
           else last_c=c;
        }
        if(c=='>') while((c=fgetc(fptr))!=EOF){ if(c=='\n') break; }

	while((c=fgetc(fptr)) == '+'){ // scan additional information.
	   ScanOverRFF(fptr);
	} // if(ungetc(c,fptr) == EOF) print_error("ReadSeqFA( )-> fasta input error");

        while(c!='>') {
             if(isalpha(c)) { length++; }
	     else if(c == '-') length++;
	     else if(!isspace(c)) {
                  fprintf(stderr,"seq %d: illegal character -> %c",I,c);
                  while((c=fgetc(fptr)) != EOF) {
                          fprintf(stderr,"%c",c);
                          if(c == '\n') break;
                  }
                  fprintf(stderr,"\n");
                  seq_error("fatal error.");
             } last_c=c;
             if((c=fgetc(fptr))==EOF) break;
        }
	if(last_c != '\n' && c != EOF) print_error("fasta input error");
	else if(length==0) print_error("fasta input error");
        fclose(fptr);
        fptr = open_file(infile,"","r");
	E = ReadSeq(fptr, I, length, A);
        fclose(fptr);
	return E;
}

char    IsSameSeqFastX(e_type E1, e_type E2,Int4 *Start,Int4 MinOverlap)
{ Int4 NumX; return IsSameSeqFastX(E1, E2,Start,&NumX, MinOverlap); }

char    IsSameSeqFastX(e_type E1, e_type E2,Int4 *Start,Int4 *RtnNumX,Int4 MinOverlap)
// find out whether or not E1 and E2 are the same sequence.
// return 1 if seq E1 lacks an N-terminal extension.
// return 2 if seq E2 lacks an N-terminal extension.
// Sets Start to the position in one sequence corresponding to the start of the other.
{
        Int4    st,end,lenSb,lenSp,NumX;
        unsigned char   *sup,*sub;      // superseq and subseq

	// Ignore X's on either ends...
	Int4 start1=0,end1=0;
        for(start1=1; ResSeq(start1,E1)==0; start1++) ;
        for(end1=LenSeq(E1); ResSeq(end1,E1)==0; end1--) ;

	Int4 start2=0,end2=0;
        for(start2=1; ResSeq(start2,E2)==0; start2++) ;
        for(end2=LenSeq(E2); ResSeq(end2,E2)==0; end2--) ;

	*RtnNumX=0;
	if(end1 < MinOverlap || end2 < MinOverlap) return 0;
	
	end=end2-MinOverlap+1;
	//          st=end                            st=1
	//    E2 ----+---------+   <-- to...from <--   +---------+-------  E2
	//           |--MinOL--|  	               |--MinOL--|
	//    E1     +---------+-----                  +---------+------   E1
	if(end > 0){
   	  sup=SeqPtr(E2); lenSp=LenSeq(E2);	sub=SeqPtr(E1); lenSb=LenSeq(E1); 
          for(st=start2; st<= end; st++){
	    // if((lenSp-st+1) < MinOverlap) break;	// guarranteed to be long enough.
	    NumX=is_same_seq_fast(sup+st,sup+lenSp+1,sub+1,sub+lenSb+1);
	    if(NumX >=0){			// sequences match!
	      // if((lenSb-st+1-NumX) < MinOverlap) break;	// can only get shorter from here...
	      // if((end-st < NumX)) return 0;  // short perfect match; assume no perfect match below.
	      if((end-st < NumX)) continue; // could be due to a string of X residues on end.
	      *Start=st-1; *RtnNumX=NumX; return 1; 
	    }
          }
	}
	//  E2  +---------+-------                         +---------+-------  E2 
	//      |--MinOL--|                                |--MinOL--|
	//  E1  +---------+------  --> from..to -->  ------+---------+         E1
	//     st=1                                        st=end
	end=end1-MinOverlap+1;
	if(end > 0){
	  sup=SeqPtr(E1); lenSp=LenSeq(E1);	sub=SeqPtr(E2); lenSb=LenSeq(E2); 
          for(st=start1; st<= end; st++){
	    // if((lenSp-st+1) < MinOverlap) break;
	    NumX=is_same_seq_fast(sup+st,sup+lenSp+1,sub+1,sub+lenSb+1);
	    if(NumX >=0){			// sequences match!
	      // if((lenSb-st+1-NumX) < MinOverlap) break;	// can only get shorter from here...
	      //if((end-st < NumX)) return 0;	// perfect match but too short...
	      if((end-st < NumX)) continue;	// could be due to a string of X residues on end.
	      *Start=st-1; *RtnNumX=NumX; return 2; 
	    }
	  }
	} return 0;
}

void	seq_error(const char *s){fprintf(stderr,"Seq: %s\n",s);exit(1);}

