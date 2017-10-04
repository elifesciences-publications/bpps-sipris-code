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

#include "seqset.h"

ss_type	MkXSeqSet(char *filename,a_type A) 
/** return NULL if can't xnu **/
{ return xnu_seqset(seqset(filename,A)); }

ss_type	PSegSeqSet(ss_type data){ return xnu_seqset(data); }

ss_type	SeqSet1(char *filename, e_type E,a_type A)
/* create a seqset of one sequence E */
{
	ss_type	P;
	e_type	E2 = CopySeq(E);
	Int4	s,r;

	NEW(P,1,seqset_type);
	P->indel_penalty =0.0;
	P->name = AllocString(filename); P->A = A;
	P->nent = 1; 
	NEW(P->entity,2,e_type); 
	NEW(P->counts,nAlpha(A)+1,Int4);
	P->entity[1]=E2;
	for(s=1; s<= (Int4) LenSeq(E2); s++){
		r = ResSeq(s,E2);
		P->counts[r]++;
	}
	P->max_leng = P->min_leng = P->total = LenSeq(E2);
	P->tfreq = NULL;
	calcseqsetfreq(P);
	P->xnu = FALSE;
	return (P);
}

void	RenameSeqSet(char *name, ss_type P)
{ char *tmp=AllocString(name); free(P->name); P->name=tmp; }

Int4	AveSeqSeqSet(ss_type data)
{ return (Int4) ceil((double)data->total/(double)data->nent); }

ss_type	Array2SeqSet(e_type *EList, Int4 N, char *filename,a_type A) 
/* create a seqset from the input list with alphabet A. */
/* WARNING: EList is absorbed into the new object. */
{
	ss_type	P;
	Int4	i,length,s,r;
	e_type	E;

	NEW(P,1,seqset_type);
	P->max_num_seqs = MAX_NUMBER_SEQS;
	P->indel_penalty =0.0;
	NEW(P->name,strlen(filename)+2,char);
	strcpy(P->name,filename); P->A = A;
assert(N > 0);
	if(N > P->max_num_seqs) {
		   fprintf(stderr,"greater than %d sequences in input file\n",
			P->max_num_seqs);
		   seqset_error("too many sequences; reset MAX_NUMBER_SEQS");
	} else if(N <= 0) seqset_error("no sequences in input List."); 
	P->nent = N;
	NEW(P->entity,P->nent+1,e_type); 
	NEW(P->counts,nAlpha(A)+1,Int4); P->total = 0;
	P->max_leng = 0; P->min_leng = INT4_MAX; 
	for(i=1; i<=P->nent; i++){
	   P->entity[i]=E=EList[i];
	   length = LenSeq(E);
	   P->max_leng = MAXIMUM(Int4,P->max_leng,length);
	   P->min_leng = MINIMUM(Int4,P->min_leng,length);
	   for(s=1;s<=(Int4)LenSeq(E);s++){r=ResSeq(s,E);P->counts[r]++;}
	   P->total += LenSeq(E);
	   if(P->total > MAX_NUMBER_RESIDUES){
                seqset_error("MAX_NUMBER_RESIDUES overflow error. (exiting...)");
	   }
	}
	P->tfreq = NULL; calcseqsetfreq(P); P->xnu = FALSE;
	free(EList);
	return (P);
}

void    ReNumberSeqSet(ss_type data)
// Renumber the sequences in data from 1 to N sequentially.
{
       Int4 N = NSeqsSeqSet(data);
       for(Int4 i=1; i<= N; i++) EqSeqI(i,SeqSetE(i,data)); 
}

BooLean	IsInSeqSet(e_type E, ss_type data)
// Determines whether a sequence is in data based on its actual sequence.
{
       Int4 N = NSeqsSeqSet(data);
       for(Int4 i=1; i<= N; i++) if(IdentSeqs(E,SeqSetE(i,data))) return TRUE;
       return FALSE;
}

ss_type	AddSeqSet(e_type E, ss_type P)
/** add the nth sequence from P **/
{
	e_type	E2,*entity;
	Int4	i,s,r,k,n;

	n = P->nent+1;
	NEW(entity,n+3,e_type); 
	entity[n]=CopySeq(E); EqSeqI(n,entity[n]);
	for(i=1; i <= P->nent; i++) entity[i]=P->entity[i];
	free(P->entity); P->entity = entity; P->nent++;
	for(s=1; s<= (Int4) LenSeq(entity[n]); s++){
		r = ResSeq(s,entity[n]); P->counts[r]++;
	}
	P->total = P->max_leng = 0; 
	P->min_leng = INT4_MAX;
        for(s=1;s<=P->nent; s++){
		E2 = P->entity[s];
		k = LenSeq(E2);
		if(P->max_leng < k) P->max_leng = k;
		if(P->min_leng > k) P->min_leng = k;
		P->total += LenSeq(E2);
	   if(P->total > MAX_NUMBER_RESIDUES){
                seqset_error("MAX_NUMBER_RESIDUES overflow error. (exiting...)");
	   }
	} calcseqsetfreq(P);
	return P;
}

ss_type	RemoveSeqSet(Int4 n, ss_type P)
/** remove the nth sequence from P **/
{
	e_type	E2;
	Int4	s,r,k;

	E2 = P->entity[n];
	for(s=1; s<= (Int4) LenSeq(E2); s++){
		r = ResSeq(s,E2); P->counts[r]--;
	}
	NilSeq(E2);
        for(s=n;s < P->nent;s++){
	   P->entity[s] = P->entity[s+1];
	}
	P->nent--; 
	P->total = P->max_leng = 0; 
	P->min_leng = INT4_MAX;
        for(s=1;s<=P->nent; s++){
		E2 = P->entity[s];
		k = LenSeq(E2);
		if(P->max_leng < k) P->max_leng = k;
		if(P->min_leng > k) P->min_leng = k;
		P->total += LenSeq(E2);
	}
	calcseqsetfreq(P);
	return P;
}

void    SetIndelPenaltySeqSet(double penalty, ss_type P)
{ P->indel_penalty=penalty; }

BooLean	IsInSeqSet(ss_type data, char *seqid)
// Determines whether a sequence is in data based on its seqid.
{
       char str[108];
       Int4 N = NSeqsSeqSet(data);
       for(Int4 i=1; i<= N; i++) {
             StrSeqID(str,100,SeqSetE(i,data));
             if(strncmp(seqid,str,100) == 0){  return TRUE; }
       } return FALSE;
}

BooLean	*AreInSeqSet(ss_type data, ss_type keydata)
// Determines which sequences in data are in keydata based on seqid.
{
       char str[108],str2[108];
       Int4 N = NSeqsSeqSet(data);
       Int4 N2 = NSeqsSeqSet(keydata);
       BooLean	*arein;
       NEW(arein,N+3,BooLean); 
       for(Int4 i=1; i<= N; i++) {
	     arein[i] = FALSE;
             StrSeqID(str,100,SeqSetE(i,data));
             for(Int4 j=1; j<= N2; j++){
                 StrSeqID(str2,100,SeqSetE(j,keydata));
                 if(strncmp(str2,str,100) == 0){  arein[i] = TRUE; break; }
             }
       } return arein;
}

ss_type	RmSeqSet(e_type E, ss_type P)
/** remove all sequences in P that are identical to E **/
/** WARNING: THIS APPEARS NOT TO BE CALLED BY ANYTHING! **/
{
	e_type	E2;
	Int4	s,r,n,k;

        for(n=1;n<=P->nent;){
		E2 = P->entity[n];
		if(IdentSeqs(E2, E)){
		   for(s=1; s<= (Int4) LenSeq(E2); s++){
			r = ResSeq(s,E2); P->counts[r]--;
		   }
		   NilSeq(E2);
		   P->entity[n] = P->entity[P->nent];
		   P->nent--; 
		} else n++;
	}
	P->total = P->max_leng = 0; 
	P->min_leng = INT4_MAX;
        for(n=1;n<=P->nent; n++){
		E2 = P->entity[n];
		k = LenSeq(E2);
		if(P->max_leng < k) P->max_leng = k;
		if(P->min_leng > k) P->min_leng = k;
		P->total += LenSeq(E2);
	}
	calcseqsetfreq(P);
	return P;
}

ss_type	SeqSet(char *filename,a_type A) { return seqset(filename,A); }

ss_type	MkSeqSet(char *filename,a_type A) { return seqset(filename,A); }

ss_type	MergeSeqSets(ss_type data1, ss_type data2)
// returns a ss_type that concatenates data2 after data1.
{
	FILE    *fp=tmpfile();
	PutSeqSetEs(fp,data1); PutSeqSetEs(fp,data2); rewind(fp);
        ss_type data=SeqSet_fptr(fp,SeqSetA(data1)); fclose(fp);
	return data;
}

ss_type	SeqSet_fptr(FILE *fptr,a_type A) { return fptr_seqset(fptr,A); }

ss_type	fptr_seqset(FILE *fptr,a_type A) 
/* create a seqset from the input file with using alphabet A. */
{
	ss_type	P;
	Int4	i,s,r,*nsize;
	e_type	E;

	NEW(P,1,seqset_type);
	P->max_num_seqs = MAX_NUMBER_SEQS;
	NEW(nsize,P->max_num_seqs+2,Int4);
	P->indel_penalty =0.0;
	NEW(P->name,25,char);
	strcpy(P->name,"temp_file"); P->A = A;
	P->nent = count_seqset_entities(fptr,P,nsize); 
	NEW(P->entity,P->nent+1,e_type); 
	NEW(P->counts,nAlpha(A)+1,Int4); P->total = 0;
	rewind(fptr);
	for(i=1; i<=P->nent; i++){
	   E = ReadSeq(fptr,i,nsize[i],A);
	   P->entity[i]=E;
	   for(s=1;s<=(Int4)LenSeq(E);s++){r=ResSeq(s,E);P->counts[r]++;}
	   P->total += LenSeq(E);
	   if(P->total > MAX_NUMBER_RESIDUES){
                seqset_error("MAX_NUMBER_RESIDUES overflow error. (exiting...)");
	   }
	}
	P->tfreq = NULL; calcseqsetfreq(P); P->xnu = FALSE;
	free(nsize);
	return (P);
}

ss_type MakeSeqSet(char *name,Int4 max_len,a_type A)
{
	ss_type	P;
	FILE	*fptr;
	Int4	i,s,r,max_defline=MAX_SEQ_DEFLINE_LENG;
	e_type	E;
	char	*defline;
	UInt4 I;
	unsigned char *seq;
	sqlnk_typ head,tmp,last;

	NEW(seq,max_len+2,unsigned char); NEW(defline,max_defline+4,char);
	NEW(P,1,seqset_type);
	NEW(P->name,strlen(name)+2,char); strcpy(P->name,name); 
	P->A = A; P->indel_penalty =0.0;
	fptr = OpenSeqSetFile(P);
	for(head=0,P->nent=0,I=1; TRUE; I++){
#if 1
		E = ReadSeqStrII(fptr,I,defline,seq,max_len,A);
		if(E == 0) break;	// i.e., end of file...
		NEW(tmp,1,sqlnk_type);
		tmp->E = E;
#endif
#if 0
		Int4 len = ReadSeqStr(fptr,defline,seq,max_len,A);
		if(len==0) break;  	// i.e., end of file...
		NEW(tmp,1,sqlnk_type);
		tmp->E = StringToSeq2(seq, len, defline, I);
#endif
		if(head == 0){ head = last = tmp; }
		else { last->next = tmp; last = tmp; }
		P->nent++;
	} fclose(fptr); free(defline); free(seq);
	assert(P->nent > 0);
	NEW(P->entity,P->nent+1,e_type); 
	NEW(P->counts,nAlpha(A)+1,Int4); P->total = 0;
	P->max_leng = 0; P->min_leng = INT4_MAX;
	for(last=head,i=1; i<=P->nent; i++){
	   E = last->E;
	   P->entity[i]=E;
	   if(LenSeq(E) > P->max_leng) P->max_leng=LenSeq(E);
	   if(LenSeq(E) < P->min_leng) P->min_leng=LenSeq(E);
	   for(s=1;s<=(Int4)LenSeq(E);s++){r=ResSeq(s,E);P->counts[r]++;}
	   P->total += LenSeq(E);
	   if(P->total > MAX_NUMBER_RESIDUES){
                seqset_error("MAX_NUMBER_RESIDUES overflow error. (exiting...)");
	   }
	   last = last->next; free(head); head = last;
	} P->tfreq = NULL; calcseqsetfreq(P); P->xnu = FALSE;
	return (P);
}

ss_type	seqset(char *filename,a_type A) 
/* create a seqset from the input file with segment length k and */
/* alphabet A. */
{
	ss_type	P;
	FILE	*fptr;
	Int4	i,s,r,*nsize;
	e_type	E;

	NEW(P,1,seqset_type);
	P->max_num_seqs = MAX_NUMBER_SEQS;
	NEW(nsize,P->max_num_seqs+2,Int4);
	P->indel_penalty =0.0;
	NEW(P->name,strlen(filename)+2,char);
	strcpy(P->name,filename); P->A = A;
	fptr = OpenSeqSetFile(P);
	P->nent = count_seqset_entities(fptr,P,nsize); 
        fclose(fptr);
	NEW(P->entity,P->nent+1,e_type); 
	NEW(P->counts,nAlpha(A)+1,Int4); P->total = 0;
	fptr = OpenSeqSetFile(P);
	for(i=1; i<=P->nent; i++){
	   E = ReadSeq(fptr,i,nsize[i],A);
	   P->entity[i]=E;
	   for(s=1;s<=(Int4)LenSeq(E);s++){r=ResSeq(s,E);P->counts[r]++;}
	   P->total += LenSeq(E);
	   if(P->total > MAX_NUMBER_RESIDUES){
                seqset_error("MAX_NUMBER_RESIDUES overflow error. (exiting...)");
	   }
	}
	fclose(fptr); P->tfreq = NULL; calcseqsetfreq(P); P->xnu = FALSE;
	free(nsize);
	return (P);
}

FILE    *OpenSeqSetFile(ss_type P)
{
        FILE    *fptr;
        if((fptr = fopen(P->name,"r")) == NULL) {
                fprintf(stderr,"Could not open file \"%s\"\n",P->name);
                seqset_error("File does not exist!\n");
        }
        return fptr;
}

e_type	*NilSeqSetRtnSeqs(ss_type P)
// sequence structures are borrowed; don't 
{
	e_type	*entity=P->entity;

	free(P->name); free(P->counts);
	if(P->tfreq != NULL) free(P->tfreq); 
	free(P);
	return entity;
}

e_type	ReplaceSeqSet(Int4 n, e_type E, ss_type P)
// replaces sequece n in P with E and returns oldE.
// returns NULL if unsuccessful. 
{
	e_type	E2;
	Int4	s,r,k;

	if(n < 1 || n > P->nent) return NULL;
	// 1. Remove old sequence.
	E2 = P->entity[n];
	for(s=1; s<= (Int4) LenSeq(E2); s++){
		r = ResSeq(s,E2); P->counts[r]--;
	}
	P->total -= LenSeq(E2); 
	// 2. Add new sequence.
	if(P->xnu) ProcessSeqPSeg(17,2.2,2.5,100,E,P->A); 
	P->total += LenSeq(E); 
	if(P->total > MAX_NUMBER_RESIDUES){
                seqset_error("MAX_NUMBER_RESIDUES overflow error. (exiting...)");
	}
	P->entity[n] = E;
	for(s=1; s<= (Int4) LenSeq(E); s++){ 
		r = ResSeq(s,E); P->counts[r]++; 
	}
	P->max_leng = 0; P->min_leng = INT4_MAX;
        for(s=1;s<=P->nent; s++){
		k = LenSeq(P->entity[s]);
		if(P->max_leng < k) P->max_leng = k;
		if(P->min_leng > k) P->min_leng = k;
	}
	calcseqsetfreq(P);
	return E2;
}

void	ReverseSeqSet(ss_type P)
// reverse the sequences in P. 
{ for(Int4 i=1; i<=P->nent; i++) ReverseSeq(P->entity[i]); }

ss_type	CopySeqSet(ss_type P)
{
	ss_type	P2;
	Int4	i,r;
	a_type	A=P->A;

	NEW(P2,1,seqset_type);
	P2->indel_penalty = P->indel_penalty;
	P2->name = AllocString(P->name); 
	// P2->A = CopyAlphabet(A);
	P2->A = P->A;
	P2->nent = P->nent; P2->xnu = P->xnu;
	P2->min_leng = P->min_leng;
	P2->max_leng = P->max_leng;
	NEW(P2->entity,P2->nent +2, e_type);
	NEW(P2->counts,nAlpha(A)+1,Int4);
	NEW(P2->tfreq,nAlpha(A)+1,double);
	for(r=0; r<=nAlpha(A); r++){
		P2->counts[r]=P->counts[r];
		P2->tfreq[r]=P->tfreq[r];
	}
	P2->total=P->total;
	for(i=1; i<=P->nent; i++){
	   if(P->entity[i] != NULL){
		P2->entity[i]=CopySeq(P->entity[i]); 
	   }
	}
	return P2;
}

ss_type	NilSeqSet(ss_type P)
{
	Int4 i;
	for(i=1; i<=P->nent; i++){
		if(P->entity[i] != NULL) NilSeq(P->entity[i]); 
	}
	free(P->name); free(P->entity); free(P->counts);
	if(P->tfreq != NULL) free(P->tfreq); 
	free(P);
	return (ss_type) NULL;
}

ss_type	xnu_seqset(ss_type P)
/*** remove low complexity regions from SeqSet using seg ***/
{
        Int4	i;
	e_type	E;
	a_type	A=P->A;

	if(P->xnu) return P;
	else P->xnu = TRUE;
        for(i=1;i<=NSeqsSeqSet(P);i++){
	   E = P->entity[i];
           if(LenSeq(E) != 0) {
		ProcessSeqPSeg(17, 2.2,2.5,100,E,A); 
		// I might want to fix this better??? 
	   }

	}
	calcseqsetfreq(P);
	return P;
}

double  LogL0SeqSet(ss_type P)
{
        double  *freq,L0,n,r;
        Int4     b;

        freq = tFreqSeqSet(P);
        for(L0=0.0, b=1; b<= nAlpha(P->A); b++){
                    if(CountsSeqSet(b,P) > 0){
                        n = (double) CountsSeqSet(b,P);
                        r = freq[b];
                        L0 += n * log(r);
                    }
        }
        return (1.4427*L0);
}

Int4	*LengthsSeqSet(ss_type P)
/* returns an array containing the sequence lengths */
{
	Int4	*len_seq,n;

	NEW(len_seq,NSeqsSeqSet(P) +1,Int4);
	for(n=1; n<= NSeqsSeqSet(P); n++)len_seq[n]=SqLenSeqSet(n,P);
	return len_seq;
}

/******************** Counting and Numbering Operations *******************/
Int4     count_seqset_entities(FILE *fptr,ss_type P, Int4 nsize[])
{
        Int4 i=0,j,length,c; 

	P->max_leng = 0; P->min_leng = INT4_MAX; 
	while((c=fgetc(fptr))!=EOF){ if(c=='>') break; }
        for(i=1,length=0;c!=EOF;length=0,i++) { 
		if(c=='>'){
		  while((c=fgetc(fptr))!=EOF){ if(c=='\n') break; }
        	  while((c=fgetc(fptr)) == '+'){ // scan over additional information.
		    ScanOverRFF(fptr);
		    // while((c=fgetc(fptr))!=EOF){ if(c=='\n') break; }
        	  } // okay to get first character in sequence here...
		}
		while(c!='>') {
		   if(isalpha(c)){ length++; }
		   else if(c == '-'){
			fprintf(stderr,"WARNING: '-' character! (Convert to 'X'?)\n");
			length++; 
		   } else if(!isspace(c)) {
			fprintf(stderr,"seq %d: illegal character -> %c",i,c);
			for(j=0; (c=fgetc(fptr)) != EOF; j++) {
				if(c == '\n') break;
				fprintf(stderr,"%c",c);
				if(j > 10 || isspace(c)) break;
			}  
			fprintf(stderr,"\n");
			seqset_error("input file error - fatal.");
		   } 
           	   if((c=fgetc(fptr))==EOF) break; 
	     	}
		if(i > P->max_num_seqs) {
		   fprintf(stderr,"greater than %d sequences in input file\n",
			P->max_num_seqs);
		   seqset_error("too many sequences; reset MAX_NUMBER_SEQS");
		}
	   	P->max_leng = MAXIMUM(Int4,P->max_leng,length);
		P->min_leng = MINIMUM(Int4,P->min_leng,length);
		nsize[i] = length;
	}
	i--;
	if(i <= 0) seqset_error("no sequences in input file."); 
        return i;
}

void	PutLengthsSeqSet(FILE *fptr,ss_type P)
{
        h_type	H=Histogram("sequence lengths", 0, P->max_leng+1,10.0);
        for(Int4 i=1; i<=NSeqsSeqSet(P); i++) IncdHist(LenSeq(SeqSetE(i,P)),H);
        PutHist(fptr,60,H); NilHist(H); 
}

ss_type	PutSeqSet(FILE *fptr,ss_type P)
{
	fprintf(fptr,"\n  input file:\n");
	fprintf(fptr,"\tname: \"%s\"\n\ttotal sequences: %d (Seg = %d)",
			P->name,P->nent,P->xnu);
	fprintf(fptr,"\n\tsequence lengths: %d-%d residues\n",
		       P->min_leng,P->max_leng);
	fprintf(fptr,"\tindel_penalty: %.2f\n",P->indel_penalty);
	fprintf(fptr,"\n\ttotal residues: %d\n",P->total);
	return P;
}

/*********************** Put SeqSet Entities Operations **********************/
ss_type	PutSeqSetEs(FILE *fptr,ss_type P)
/* print all sequence entities in seqset using fasta format */
{
	Int4     i;
	for(i=1;i<=P->nent; i++) {
 	   if(SeqI(P->entity[i]) !=0) PutSeqSetE(fptr, i,P);
	}
	fprintf(fptr,"\n\n");
	return P;
}

void	PutTrimmedSeqSet(FILE *fptr, Int4 trim, ss_type P)
// put the sequences trimming trim residues from the N-terminus
// if trim < 0 then trim from C-terminus
{
   e_type	E;
   Int4		i;
   if(trim > 0){
      for(i=1;i<=P->nent; i++) {
	E = P->entity[i];
	PutSubSeq(fptr, 1+trim, LenSeq(E), E, P->A);
      }
   } else if(trim < 0){
      for(i=1;i<=P->nent; i++) {
	E = P->entity[i];
	PutSubSeq(fptr, 1, LenSeq(E)+trim, E, P->A);
      }
   } else PutSeqSetEs(fptr,P);
}

ss_type	PutSeqSetE(FILE *fptr, Int4 i, ss_type P) /* print the ith entity */
{
	e_type	E;

	if(i <= P->nent && i > 0) { E = P->entity[i]; PutSeq(fptr,E,P->A); }
	return P;
}

ss_type	PutSeqSetPIDs(FILE *fptr, ss_type P)
/* print entity ids for selected entities */
{
	e_type	E;
	Int4     i;
	for(i=1;i<=P->nent;i++) {
	   E = P->entity[i];
 	   if(SeqI(E) != 0){
		fprintf(fptr,"#%-3d ",(Int4)SeqI(E));
		PutSeqInfo(fptr,E);
	   }
	}
	fprintf(fptr,"\n\n");
	return P;
}

/********************** Frequency Operations ********************/

ss_type	calcseqsetfreq(ss_type P)
/* calculate the residue frequencies for seqset */
{
    Int4		s;
    a_type	A=P->A;

    if(P->tfreq==NULL) NEW(P->tfreq,nAlpha(A)+1,double);
    for(s=0;s<=(Int4) nAlpha(A);s++) {
	P->tfreq[s] = (double) P->counts[s]/(double) P->total;
    }
    return P;
}

double	SeqSetEntropy(ss_type P)
{
	Int4	i;
	double	*freq,H;

	freq = tFreqSeqSet(P);
	for(H=0.0,i = 1; i <= nAlpha(P->A); i++){
		if(freq[i] > 0.0) H += freq[i] * log(freq[i]);
	} return (-1.442695041*H);
}

ss_type	PutSeqSettFreqs(FILE *fptr,ss_type P)
{
	Int4 i; double T=0.0;

	fprintf(fptr,"RES    %-6s %-s\n","NUM","FREQ");
	for(i=0;i<=nAlpha(P->A);T+=P->tfreq[i],i++)
	    fprintf(fptr,"%c(%2d): %-6d %-2.3f\n",
			AlphaChar(i,P->A),i,P->counts[i],P->tfreq[i]);
	fprintf(fptr,"TOTAL: %-6d %-2.3f\n\n", P->total, T);
	return P;
}

/************************* Randomization Routines ***********************/

ss_type	ShuffleSeqSet2(ss_type P)
{ 
	Int4 r,s,i,n,item;
	dh_type H;
	e_type	E;
	char	*S;

	for(n=0,i=1;i<=P->nent;i++) { E=P->entity[i]; n += LenSeq(E); }
	H = dheap(n+2,4);
	NEW(S,n+2,char);
	for(item=i=1;i<=P->nent;i++) {
		E=P->entity[i];
		for(s=1; s<= (Int4) LenSeq(E); s++){
			r = ResSeq(s,E);
			insrtHeap(item,((keytyp)Random()),H);
			S[item++]=r;
		}
	}
	for(i=1;i<=P->nent;i++) {
		E=P->entity[i];
		for(s=1; s<= (Int4) LenSeq(E); s++){
			item=delminHeap(H);
			if(item==0) seqset_error("shuffleSeqSet2 error");
			r=S[item];
			EqSeq(s,r,E);
		}
		EqSeqI(i,P->entity[i]);
	}
	Nildheap(H); free(S);
	return P;
}

ss_type	ShuffleSeqSet(ss_type P)
{ 
	Int4 i;
	for(i=1;i<=P->nent;i++) {
		ShuffleSeq(P->entity[i]);
		EqSeqI(i,P->entity[i]);
	}
	return P;
}

void	seqset_error(const char *s)
{fprintf(stderr,"Seq_Set: %s\n",s); exit(1);}

