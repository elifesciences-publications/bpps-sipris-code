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

#include "psimatrix.h"

BooLean	TestPsiMatrix(  )
{
	/*** PNTR EditScriptPtr	esp; /***/
	EditScriptPtr	*esp; /****/
	SeqEntryPtr sep;
	BioseqPtr bsp;
	FILE *fp;
	SeqLocPtr slp;
	SeqPortPtr spp;
	Int2	i=1,len;
	posSearchItems  *posSearch;
	char	seq1[] = "HPVYNMVMSDVTISCTSLEKEKREEVHKYVQMMGGRVYRDLNVSVT";
	char	seq2[] = "EDFKCPIFLGCIICVTGLCGLDRKEVQQLTVKHGGQYMGQLKMNEC";
	
	fp = tmpfile();
	fprintf(fp,">lcl|10001 test seq\n%s\n\n",seq1);
	fprintf(fp,">lcl|10002 test seq\n%s\n\n",seq2);
	rewind(fp);
	while((sep = FastaToSeqEntry(fp, FALSE)) != NULL) {
	     bsp = (BioseqPtr) sep->data.ptrvalue;
	     slp = SeqLocIntNew(0, bsp->length-1, 0, bsp->id);
	     len = SeqLocLen(slp);
	     printf("length seq %d = %d\n",i++,SeqLocLen(slp));
	     spp = SeqPortNew(bsp, 0, len-1, 0, 0);
	     SeqPortTell(spp);
	     spp = SeqPortFree(spp);
	}
#if 0
	posAllocateMemory(posSearch, 20, 100, 40);
#endif
	return TRUE;
}

/************************** Some functions ******************************
NLM_EXTERN SeqLocPtr LIBCALL FindSpliceSites(SeqEntryPtr sep, Boolean findOnProtein);NLM_EXTERN SeqFeatPtr LIBCALL FindCodingRegion(SeqEntryPtr sep);

/***********************************************************************/

/************************** From Jinghui ******************************
In the tools directory of the ncbi toolkit, there is file simutil.h.
The Webb structure that I talked to you about is like this:

#define INS ((Uint1)1)
#define DEL ((Uint1)2)

/* EditScript is the internal representation of an alignment in SIM3. *\
typedef struct edit_script {
        Uint1 op_type;  /* SUB, INS, or DEL *\
        Int4 num;       /* Number of operations *\
        struct edit_script PNTR next;
} EditScript, PNTR EditScriptPtr;

The function that I told you about is:
/* *       convert EditScript Structure into a Seq-align * *\
SeqAlignPtr ConvertEditScript PROTO((EditScriptPtr esp, SeqLocPtr slp1,
SeqLocPtr slp2, Int4 offset_1, Int4 offset_2, Boolean show_end_gap));

In this function, slp1, slp2, defines a location on the two sequences.
you can set offset_1, offset_2 to 0. The parameter show_end_gap is an option
to decide when the end of the alignment is gap, wether to include it in
Seq-align.

Hope this will be useful to you.

(PS. to make a Seq-loc, there is a function:
SeqLocIntNew, which takes the start, stop, strand and the id of the
sequence)

write your sequence into a temp file:

>lcl|10001
GATCCTAATTTTAAGCTTTTGCAATTCTTATTTTTTTTTCAAAAAAAGTAAAATAATAAAATGGAAAAGA
AGCCCGGAACCACACCCCCCCCCACCACCCAGGGCCCTCCCCAACACACACATATTTAAGAATAAGAGAA
AATGCCCTTGGGACATTCTCACACCTTTAGGAGAAGTAACTTTCACAATGCTTATATTTAAAATTCCTCC
ACCACTATCCTCCTTTTTTTTTTCTTCATTTCTTGACATAATTGGGCAGAGTAGGCAGAGGCATATTTCT


SeqEntryPtr sep;
BioseqPtr bsp;
FILE *fp;
SeqLocPtr slp;

while((sep = FastaToSeqEntry(fp, FALSE)) != NULL)
{
     bsp = sep->data.ptrvalue;
     slp = SeqLocIntNew(0, bsp->length-1, 0, bsp->id);
}

/**********************************************************************/

