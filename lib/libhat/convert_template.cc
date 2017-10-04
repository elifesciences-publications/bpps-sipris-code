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
#include "hat_typ.h"
#include "residues.h"

#define	USAGE_START	"USAGE: convert_tmplt template file [options]\n\
    Realign the sequences in file.cma based on the file template.cma\n\
      Note: the first sequence in file.cma must occur in the template aligment\n\n\
    Thus, the run_convert program realigns a (more extensive) subgroup template alignment \n\
      based on a (less extensive) supergroup template alignment.  This allows modeling of \n\
      subgroups sharing more extensive sequence similarity with one template (e.g., the AGC kinases)\n\
      that can then be incorporated into a supergroup template alignment for analysis \n\
      of a larger (super)group (e.g., all EPKs).  \n\
\n\n"

/**************************** Global Variables ******************************/
int	convert_template_main(Int4 argc,char *argv[])
{ 
	Int4	arg,i,j,s,cutoff=-999,blk=0,lenrm,mingap;
	a_type	A;
	UInt4   seed=7061950;
	cma_typ	cmsa,cmsa2,cma=0;
	cma_typ	tmplt=0;
	FILE	*fp;
	char	str[300],mode=' ';

	Int4 time1=time(NULL); 
	if(argc < 3) print_error(USAGE_START);
	for(arg = 3; arg < argc; arg++){
	   if(argv[arg][0] != '-') print_error(USAGE_START);
	   switch(argv[arg][1]) {
             case 'x': mode = 'x'; break;
	     default: print_error(USAGE_START);
	   }
	}
	if(seed == 7061950) seed = (UInt4) time(NULL);
	sRandom(seed);
	A = MkAlpha(AMINO_ACIDS,PROT_BLOSUM62);
	sprintf(str,"%s.cma",argv[1]);

	tmplt=ReadCMSA2(str,A);
	if(!tmplt) print_error("cma file read error");
	Int4 Number=NumSeqsCMSA(tmplt);

#if 1
        Int4 J,NumberCMAs;
        fp=open_file(argv[2],".cma","r");
        cma_typ *MULTI_CMA=MultiReadCMSA(fp,&NumberCMAs,A); fclose(fp);
    for(Int4 f=1; f<=NumberCMAs;f++) {
#if 1		// temporary fix for overhangs on ends of sequences...
             cma=MULTI_CMA[f];
#else
	     FILE *tfp=tmpfile();
	     Int4 aln_len = PutFastaCMSA(tfp,MULTI_CMA[f]);
	     rewind(tfp); cma=FastaToCMSA(tfp,A); fclose(tfp);
#endif
#else
	sprintf(str,"%s.cma",argv[2]); cma=ReadCMSA2(str,A);
	if(!cma) print_error("cma file read error");
	{
#endif
	cma_typ *IN_CMA=0;
	NEW(IN_CMA,Number+3,cma_typ);
	BooLean okay=FALSE;
	e_type csq=FakeSeqCMSA(1,cma);
	for(Int4 i=1; i < Number; i++){
		e_type trueE = TrueSeqCMSA(i+1,tmplt);
		if(IdentSeqs(trueE,csq)){
			okay=TRUE;
			IN_CMA[i]=cma; break;
		}
	}
	if(!okay){
	    fprintf(stderr,"In input file %d\n",f);
	    PutSeq(stderr,csq,A);
// continue;
	    print_error("FATAL: Template does not contain cma file concensus sequence");
	}

	cma_typ *OUT_CMA=ConversionViaTemplateCMSA3(tmplt,IN_CMA);
	for(Int4 II=1; II < Number; II++){
	   if(OUT_CMA[II]){
#if 0
		char str2[500];
            	sprintf(str2,"%s_New.cma",argv[2],II);
            	WriteCMSA(str2,OUT_CMA[II]); break;
#else
		PutCMSA(stdout,OUT_CMA[II]); break;
#endif
           }
	} free(OUT_CMA);
	if(cma) TotalNilCMSA(cma);
    }
#if 1
	free(MULTI_CMA);
#endif
	if(tmplt) TotalNilCMSA(tmplt);
	// if(data != NULL) NilSeqSet(data);
	NilAlpha(A);
	fprintf(stderr,
		"\ttime: %d seconds (%0.2f minutes)\n",
                        time(NULL)-time1,(float)(time(NULL)-time1)/60.0);
}

