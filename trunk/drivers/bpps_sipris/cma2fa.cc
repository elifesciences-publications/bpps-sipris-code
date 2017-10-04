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

#include "cmsa.h"

#define	USAGE_START_NEW	"USAGE: cma2fa fafile [options]\n\
   options:\n\
     -A            - Use PutFastaAlnCMSA() routine\n\
     -P=<int>      - Permute cma file at column position <int>\n\
     -P=<int>      - Permute cma file at split point <int>\n\
     -F            - Use PutFastaCMSA() routine (default)\n\
     -fasta        - Use PutFastaCMSA() routine (with a consensus added)\n\
\n\n"

#define	USAGE_START	"\n\
===========================================================\n\
  cma2fa Copyright 2017 The University of Maryland\n\
\n\
 For information contact:\n\
   Andrew F. Neuwald\n\
   Institute for Genome Sciences\n\
   Department of Biochemistry & Molecular Biology\n\
   University of Maryland School of Medicine, Baltimore\n\
   E-mail: aneuwald@som.umaryland.edu\n\
===========================================================\n\
USAGE: cma2fa fafile [options]\n\
     Note: this program eliminates extensions on ends (by necessity)\n\
\n\n"

void	PermuteCMSA(FILE *fp, Int4 SplitSite, cma_typ cma)
// output a permuted version of input cma file
// permution occurs directly AFTER the SplitSite 
{
	a_type	A=AlphabetCMSA(cma);

        // 1. check input parameters 
	if(nBlksCMSA(cma) > 1) print_error(USAGE_START);
	if(SplitSite < 1 || SplitSite >= LengthCMSA(1,cma)) print_error(USAGE_START);

	// 2. add consensus to cma and make a copy.
	cma_typ	cma1=AddConsensusCMSA(cma);
	cma_typ cma2=CopyCMSA(cma1);
	Int4 N = NumSeqsCMSA(cma1);

	// 3. create an N- and C-terminal alignments
	// BooLean TrimCMSA(Int4 blk, unsigned short RmLeft, unsigned short RmRight, cma_typ cma)
	BooLean okay=TrimCMSA(1, 0, LengthCMSA(1,cma) - SplitSite,cma1);
	if(!okay) print_error("PermuteCMSA() error 1");
	okay=TrimCMSA(1,SplitSite,0,cma2);
	if(!okay) print_error("PermuteCMSA() error 2");

	// 4. print out the two halves in fasta format into tmpfiles
	FILE *fp1,*fp2;
	fp1=tmpfile();
	Int4 len1 = PutFastaCMSA(fp1,cma1); 
	Int4 n1 = len1+3;
	char *buffer1; NEW(buffer1,len1 +10,char);
	fprintf(stderr,"alignment length 1 = %d\n",len1);

	fp2=tmpfile();
	Int4 len2 = PutFastaCMSA(fp2,cma2); 
	char *buffer2; NEW(buffer2,len2 +10,char);
	Int4 n2 = len2+3;
	fprintf(stderr,"alignment length 2 = %d\n",len2);

	// 5. Read in the two halves in fasta format and print out permuted seqs.
	rewind(fp1); rewind(fp2);
	e_type E;
	for(Int4 n=1; n <= N; ){
		char c1=fgetc(fp1);
		char c2=fgetc(fp2);
		if(c1 == '>' && c2 == '>'){
			while((c1=fgetc(fp1) != '\n')); 
			while((c2=fgetc(fp2) != '\n')); 
		} else if(c1 == '\n' && c2 == '\n'){
		} else if(isprint(c1) && isprint(c2)){
		   ungetc(c1,fp1); ungetc(c2,fp2);
		   fgets(buffer1,len1+1,fp1); fgets(buffer2,len2+1,fp2);
			E= TrueSeqCMSA(n,cma1);
			printf(">"); PutSeqID(fp,E); 
			if(PhylumSeq(E)){
                  		fprintf(fp," {<%s(%c)>} (permuted)\n%s%s\n\n",
					PhylumSeq(E),kingdomSeq(E),buffer2,buffer1);
                	} else printf(" (permuted)\n%s%s\n\n",buffer2,buffer1);
			n++;
		} else if(c1 == EOF && c2 == EOF) break;
		else print_error("Permute input error 3");
	}
}

int	main(Int4 argc,char *argv[])
{ 
	Int4	arg,permute_pt;
	Int4    time1;
	char	str[300],mode='f';
	// char	str[300],mode='F';
	cma_typ	cma;
	a_type	A;
	UInt4   seed=7061950;
	Int4	aln_len;

	time1=time(NULL); 
	if(argc < 2) print_error(USAGE_START);
	for(arg = 2; arg < argc; arg++){
	   if(argv[arg][0] != '-') print_error(USAGE_START);
	   switch(argv[arg][1]) {
             case 'P': mode = 'P'; 
		if(sscanf(argv[arg],"-P=%d",&permute_pt) != 1) print_error(USAGE_START);
		break;
             case 'A': mode = 'A'; break;
             case 'F': mode = 'F'; break;
             case 'f': 
		if(strcmp("-fasta",argv[arg]) != 0) print_error(USAGE_START);
		mode = 'f'; 
		break;
	     default: print_error(USAGE_START);
	   }
	}
	if(seed == 7061950) seed = (UInt4) time(NULL);
	sRandom(seed);
	A = MkAlpha(AMINO_ACIDS,PROT_BLOSUM62);
	sprintf(str,"%s.cma",argv[1]);
	cma=ReadCMSA2(str,A);
	if(!cma) print_error("cma file read error");
	switch(mode) {
	  case 'A': PutFastaAlnCMSA(stdout,cma); break;
	  case 'F': PutFastaCMSA(stdout,cma); break;
	  case 'f': 
		{
		cma_typ cma0 = AddConsensusCMSA(cma);
		aln_len = PutFastaCMSA(stdout,cma0); 
		fprintf(stderr,"alignment length = %d\n",aln_len);
		}
	  break;
	  case 'P': 
		PermuteCMSA(stdout, permute_pt,cma);
	  break;
	  default : print_error(USAGE_START); break;
	}
	if(cma) TotalNilCMSA(cma);
	// if(data != NULL) NilSeqSet(data);
	NilAlpha(A);
	fprintf(stderr, "\ttime: %d seconds (%0.2f minutes)\n",
                        time(NULL)-time1,(float)(time(NULL)-time1)/60.0);
}

