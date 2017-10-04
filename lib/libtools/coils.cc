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

/* coils.c - coiled coils prediction program */
#include "coils.h"

csn_typ MkCoils(double pseudo, char *snfile, Int4 *counts, a_type A, 
	char method)
{
	csn_typ F;

	NEW(F,1,coils_type);
	NEW(F->M,MAX_NUM_MODELS+1,wm_type);
	if(snfile!=NULL) F->snfile = AllocString(snfile);
	else F->snfile = NULL;
	F->shuffle=FALSE;
	F->dfp = NULL;
	F->recomb = NULL;
	F->A = A; 
	F->method = method; 
	if(F->snfile != NULL) ReadCoils(F,pseudo,F->snfile,counts);
	return F;
}

void	NilCoils(csn_typ F)
{
	Int4	i;

	if(F->snfile != NULL){
	  for(i=1;i<=F->N; i++) { NilWModel(F->M[i]); }
	  free(F->M); free(F->freq); 
	  free(F->snfile);
	}
	if(F->dfp != NULL) fclose(F->dfp);
	if(F->recomb != NULL) fclose(F->recomb);
	free(F);
}

Int4	MotifInfoCoils(char *snfile, Int4 *counts, a_type A)
/** Get counts for conserved regions and return the number of blocks **/
{
	FILE	*fptr;
	Int4	N,m,r,i,length,pos;
	char	c = '\n',*null;

   NEW(null,MAX_BLOCK_LENGTH+2,char);
   for(i=0; i<=nAlpha(A); i++) counts[i]=0;
   if((fptr = fopen(snfile,"r")) == NULL) {
        fprintf(stderr,"Could not open file \"%s\"\n",snfile);
        print_error("File does not exist!");
   }
   for(c=' ',length = 0,N=0,m=1;c!=EOF && c!='@'; m++){
        for(; isspace(c); c=fgetc(fptr))	/** GO TO FIRST CHARACTER **/;
	if(c == '*'){			/** fragmentation row **/
		length = 0;
		do {
                   if(c=='*') {
			length++; null[length] = c;
		   } else if(c=='.') {
			length++; null[length] = c;
		   } else if(!isspace(c)) {
                        fprintf(stderr,"illegal character -> %c",c);
                        fprintf(stderr,"\n");
			print_error("fatal error.");
                   }
		   if(length >= MAX_BLOCK_LENGTH) 
                   		print_error("reset MAX_BLOCK_LENGTH");
                } while((c=fgetc(fptr))!='\n' && c!=EOF && c!='@');
	} else print_error("model file input error1");
	if(length == 0) print_error("input error");
        for(; isspace(c); c=fgetc(fptr))	/** GO TO FIRST CHARACTER **/;
	if(isalpha(c)){				/** = segment row **/
            for(i=1,pos=0;c!=EOF && c!='@' && c!='*';i++,pos=0) {
		do{
                   if(c=='*') { break; }
		   if(isalpha(c)) {
			if(islower(c)) c=toupper(c);
		        if(++pos >= MAX_BLOCK_LENGTH) 
                   		print_error("reset MAX_BLOCK_LENGTH");
			r = AlphaCode(c,A); counts[r]++;
		   } else if(!isspace(c)) {
                        fprintf(stderr,"seq %d: illegal character -> %c\n",i,c);
			print_error("fatal error.");
                   }
                } while((c=fgetc(fptr))!='\n' && c!=EOF && c!='@');
                if(i >= MAXSCN_BLOCK_SIZE)
                   print_error("aligned segments > MAXSCN_BLOCK_SIZE");
		else if(pos == 0)	i--;
		else if(pos != length) 
                   print_error("input segments have inconsistent lengths.");
	    }
	} else print_error("model file input error2");
        i--; 
	if(i < 1){
	    if(N==0) print_error("zero segments in snfile.");
	    else m--;	/** no segments in this pass **/
	} else if(i > 0) N++;
   }
   free(null); fclose(fptr);
   return N;
}

void	ReadCoils(csn_typ F, double pseudo,char *snfile, Int4 *counts)
{
	FILE	*fptr;
	Int4	m,i,length,number,pos,total;
	unsigned char	**seq;
	char	c = '\n',*null;
	BooLean	use_null=TRUE, verbose=FALSE;
	a_type	A = F->A;

   if(isupper(F->method)){ use_null=FALSE; F->method=tolower(F->method); }
   NEW(null,MAX_BLOCK_LENGTH+2,char);
   NEWP(seq,MAXSCN_BLOCK_SIZE+2,unsigned char);
   for(i=0; i<=MAXSCN_BLOCK_SIZE; i++) NEW(seq[i],MAX_BLOCK_LENGTH+2,unsigned char);
   NEW(F->freq,nAlpha(A)+2,double);
   for(total=i=0; i<=nAlpha(A); i++) total += counts[i];
   for(i=0; i<=nAlpha(A); i++) {
	F->freq[i]= (double)counts[i]/(double)total;
	if(verbose) fprintf(stderr,"freq[%d] = %d/%d = %g\n",
		i,counts[i],total,F->freq[i]);
   }
   F->total = (double) total;
   if((fptr = fopen(snfile,"r")) == NULL) {
        fprintf(stderr,"Could not open file \"%s\"\n",snfile);
        print_error("File does not exist!");
   }
   for(c=' ',length = 0,F->N=0,m=1;c!=EOF && c!='@'; m++){
        if(verbose) fprintf(stderr,"\n"); /****/
        for(; isspace(c); c=fgetc(fptr))	/** GO TO FIRST CHARACTER **/;
	if(c == '*'){			/** fragmentation row **/
		length = 0;
		do {
                   if(c=='*') {
			length++; null[length] = c;
			if(verbose) fprintf(stderr,"*");
		   } else if(c=='.') {
			length++; null[length] = c;
			if(verbose) fprintf(stderr,".");
		   } else if(!isspace(c)) {
                        fprintf(stderr,"illegal character -> %c",c);
                        fprintf(stderr,"\n");
			print_error("fatal error.");
                   }
		   if(length >= MAX_BLOCK_LENGTH) 
                   		print_error("reset MAX_BLOCK_LENGTH");
                } while((c=fgetc(fptr))!='\n' && c!=EOF && c!='@');
	} else print_error("model file input error3");
	if(verbose) fprintf(stderr," length = %d\n",length);
	if(length == 0) print_error("input error");
        for(; isspace(c); c=fgetc(fptr))	/** GO TO FIRST CHARACTER **/;
	if(isalpha(c)){				/** = segment row **/
            for(i=1,pos=0;c!='*' && c!=EOF && c!='@' && c!='*';i++,pos=0) {
		do{
                   if(c=='*') { break; }
		   if(isalpha(c)) {
			if(islower(c)) c=toupper(c);
                        if(verbose) fprintf(stderr,"%c",c);/***/
		        if(++pos >= MAX_BLOCK_LENGTH) 
                   		print_error("reset MAX_BLOCK_LENGTH");
			seq[i][pos] = AlphaCode(c,A);
		   } else if(!isspace(c)) {
                        fprintf(stderr,"seq %d: illegal character -> %c",i,c);
                        fprintf(stderr,"\n");
			print_error("fatal error.");
                   }
                } while((c=fgetc(fptr))!='\n' && c!=EOF && c!='@');
		if(verbose) {
		    if(c == '\n' && pos > 0) fprintf(stderr," #%d\n",i);
		}
                if(i >= MAXSCN_BLOCK_SIZE)
                   print_error("aligned segments > MAXSCN_BLOCK_SIZE");
		else if(pos == 0)	i--;
		else if(pos != length) 
                   print_error("input segments have inconsistent lengths.");
	    }
	} else print_error("model file input error4");
        i--; number = i; 
	if(i < 1){
	    if(F->N==0) print_error("zero segments in snfile.");
	    else m--;	/** no segments in this pass **/
	} else if(i > 0){
		F->N++;
		if(F->N > MAX_NUM_MODELS) print_error("increase MAX_NUM_MODELS");
		if(use_null) F->M[m]=MkWModel(null,length,pseudo,F->freq,A);
		else F->M[m]=MkWModel(NULL,length,pseudo,F->freq,A);
		SetMethodWModel(F->method,F->M[m]);
		for(i=1; i<=number; i++){
			Add2WModel(seq[i],1,1.0,F->M[m]);
/*****
for(j=1; j<=length; j++){ printf("%c", AlphaChar(seq[i][j],A));
if(i%50 == 49) printf("\n"); } printf("\n\n");
/*****/
		}
		if(verbose){
		   fprintf(stderr,"\n\n");
		   PutWModel(stderr,F->M[m]); 
		   fprintf(stderr,"\n\n"); 
		}
	}
   }
   if(verbose) fprintf(stderr,"\n\t %d motif models\n",F->N);
   for(i=0; i<=MAXSCN_BLOCK_SIZE; i++) free(seq[i]);
   free(seq); free(null); 
   fclose(fptr);
}

/******************************** KUZIO *************************************/

static FloatHi PPCC (FloatHi score)
{
  FloatHi   meang = 0.77;
  FloatHi   stdg  = 0.20;
  FloatHi   meanc = 1.63;
  FloatHi   stdc  = 0.24;
  FloatHi   constg, constc;
  FloatHi   evalg, evalc;
  FloatHi   Gg, Gc;
  FloatHi   prob;

  constg = 1 / (stdg * 2 * 3.14159);
  constc = 1 / (stdc * 2 * 3.14159);
  evalg = (((score - meang) / stdg) * ((score - meang) / stdg)) / 2;
  evalc = (((score - meanc) / stdc) * ((score - meanc) / stdc)) / 2;
  evalg = exp (evalg);
  evalc = exp (evalc);
  Gg = constg / evalg;
  Gc = constc / evalc;
  prob = Gc / ((30 * Gg) + Gc);
  return prob;
}

static FloatHi nroot (FloatHi fv, Int4 wv)
{
  FloatHi  fi = 1.0;

  fi = fi / (FloatHi) wv;
  fv = pow (fv, fi);
  return fv;
}

static void PCC (SeqPortPtr spp, Int4 start, Int4 end, Int4 window,
                 FloatHi scr[24][42], Char res[24], FloatHiPtr fptr)
{
  Int4        i, n, col;
  Int4        stop;
  Char        aa;
  Boolean     flagMatch;
  FloatHi     temp;
  FloatHiPtr  fheadptr;

  fheadptr = fptr;
  for (i = start; i < end; i++)
    *fptr++ = 0.0;
  fptr = fheadptr;
  stop = end - window;
  for (i = start; i < stop; i++) {
    SeqPortSeek (spp, i, SEEK_SET);
    temp = 1.0;
    for (col = 0; col < window; col++) {
      aa = SeqPortGetResidue (spp);
      flagMatch = FALSE;
      n = 0;
      while (res[n] != '\0') {
        if (aa == res[n]) { flagMatch = TRUE; break; }
        n++;
      }
      if (flagMatch) temp *= scr[n][col];
    }
    temp = nroot (temp, window);
    temp = PPCC (temp);
    if (temp < 0.0) temp = 0.0;
    temp *= 100.0;
    fheadptr = fptr;
    for (n = 0; n < window; n++) {
      if (temp > *fptr) *fptr = temp;
      fptr++;
    }
    fptr = fheadptr;
    fptr++;
  }
  return;
}

static Int4 ReadPCC (CharPtr res, FloatHi scr[24][42], Int4 window)
{
  FILE    *fin;
  Int2    i, n, val, count;
  Int2    numcol;
  FloatHi score;
  Char    buff[256];
  CharPtr bptr;
  static Char c[] = {'A','B','C','D','E','F','G','H','I','K','L','M',
                     'N','P','Q','R','S','T','V','W','X','Y','Z','*'};

  numcol = 7;

  if ((fin = FileOpen("/users/neuwald/COILS/KSpcc.dat", "r")) == NULL) {
    if ((fin = FileOpen("KSpcc.dat", "r")) == NULL) { return 0; }
  }

  val = 0;
  while ((FileGets (buff, sizeof (buff), fin)) != NULL) {
    if (buff[0] == ';') continue;
    res[val] = c[val];
    bptr = buff;
    for (n = 0; n < numcol; n++) {
      sscanf (bptr, "%lf", &score); scr[val][n] = score; bptr += 7;
    }
    val++;
    if (val == 24) break;
  }
  count = val;

  for (n = numcol; n < window; n++) {
    if (n == 42) break;
    val = n % numcol;
    for (i = 0; i < 24; i++) { scr[i][n] = scr[i][val]; }
  }

  FileClose (fin);
  return (Int4) count;
}

/******************************** KUZIO *************************************/

Uint1Ptr GetSequenceWithDenseSeg2(DenseSegPtr dsp, Boolean query, Int4Ptr start, 
	Int4Ptr length)
/*** taken from blast code... GetSequenceWithDenseSeg( )  ***/
{
        BioseqPtr bsp;
        Int4 index, offset;
        SeqIdPtr id;
        SeqPortPtr spp;
        Uint1Ptr buffer;

        if (dsp == NULL) return NULL;
        if (query == TRUE) { offset = 0; id = dsp->ids; }
        else { offset = 1; id = dsp->ids->next; }

        *start = dsp->starts[offset];
        *length = 0;
        for (index=0; index<dsp->numseg; index++) {
                if (dsp->starts[offset+2*index] != -1)
                        *length += dsp->lens[index];
        }
        bsp = BioseqLockById(id);
        spp = SeqPortNew(bsp, *start, (*start)+(*length)-1, 
			Seq_strand_unknown, Seq_code_ncbistdaa);
        buffer = (Uint1Ptr) MemNew((*length)*sizeof(Uint1));
        for (index=0; index<*length; index++)
                buffer[index] = SeqPortGetResidue(spp);
        spp = SeqPortFree(spp);
        BioseqUnlock(bsp);
        return buffer;
}

SeqEntryPtr	MySeqToSEP(unsigned char *seq, Int4 len_seq, char *info, 
		Int4 id, a_type A)
/**** Jinghui's method for creating an NCBI bioseq *****/
{
	FILE			*fp;
	SeqEntryPtr		sep;
	Int4			i;

	fp = tmpfile();
	if(info != NULL) fprintf(fp,">lcl|%d|%s\n", 10000 + id,info);
	else fprintf(fp,">lcl|%d \n", 10000 + id);
/***
for(i=1; i<=len_seq; i++) fprintf(stdout,"%c",AlphaChar(seq[i],A));
/** debug **/
	for(i=1; i<=len_seq; i++) fprintf(fp,"%c",AlphaChar(seq[i],A));
	fprintf(fp,"\n\n");
	rewind(fp);
	if((sep = FastaToSeqEntry(fp, FALSE)) == NULL) 
				print_error("error in MySeqToSEP( )");
	fclose(fp);
	return sep;
}

