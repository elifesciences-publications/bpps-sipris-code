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
#include "seqset.h"
#include "residues.h"
#include "blosum62.h"

#define USAGE_START "Usage: driver_seqset input_seqs [options]\n\
       -x          dummy\n\n"

int     main(int argc, char *argv[])
{
        int	arg;
        Int4    i;
        a_type  A;
        ss_type data = NULL;
        e_type  Seq;
        char	str[MAX_SEQ_DEFLINE_LENG+100],*buffer=0;
	FILE	*fptr=stderr;

        if(argc < 2) print_error(USAGE_START);
        for(arg = 2; arg < argc; arg++){
          if(argv[arg][0] == '-'){
           switch(argv[arg][1]) {
             case 'x': break;
             default : print_error(USAGE_START);
           }
	  }
        }
        A = MkAlpha(AMINO_ACIDS,PROT_BLOSUM62); /****/
        data = MakeSeqSet(argv[1],200000, A);
	PutSeqSet(stderr,data);
	PutLengthsSeqSet(stderr,data);
	PutSeqSetEs(stderr,data);
	PutSeqSettFreqs(stderr,data);


	fprintf(fptr,"<TaxSeqs(%d;%d)={",NSeqsSeqSet(data),NSeqsSeqSet(data));
	for(i=1; i < NSeqsSeqSet(data); i++){
                    Seq = SeqSetE(i,data);
                    if(KingdomSeq(Seq)) fprintf(fptr,"%s(%c1)=%d;",PhylumSeq(Seq),
                        KingdomSeq(Seq),LenSeq(Seq));
                    else fprintf(fptr,"NotKnown(X1)=%d;",LenSeq(Seq));
	}
	Seq = SeqSetE(i,data);
	if(KingdomSeq(Seq)) fprintf(fptr,"%s(%c1)=%d};\n\n",
                        PhylumSeq(Seq),KingdomSeq(Seq),LenSeq(Seq));
	else fprintf(fptr,"unknown(X)=%d};\n\n",LenSeq(Seq));
	NEW(buffer,MaxSeqSeqSet(data)+50,char);
	for(i=1; i <= NSeqsSeqSet(data); i++){
                Seq = SeqSetE(i,data);
                StrSeqID(str,MAX_SEQ_DEFLINE_LENG,Seq);
                fprintf(fptr,">%s {<%s(%c)>}",str,PhylumSeq(Seq),KingdomSeq(Seq));
                StrSeqDescript(str,MAX_SEQ_DEFLINE_LENG,Seq);
                fprintf(fptr,"%s\n",str);
                SeqToString(buffer, Seq, A);
                for(Int4 j=0; j < LenSeq(Seq); j++){
                        if(j%60 == 59) fprintf(fptr,"\n");
                        fprintf(fptr,"%c",buffer[j]);
                } fprintf(fptr,"\n");

	} fprintf(fptr,">.\n"); 
        if(buffer) free(buffer);
        NilSeqSet(data); NilAlpha(A);
}


