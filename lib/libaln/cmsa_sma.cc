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

cma_typ	SMA2CMSA(char *fafile, sma_typ MA)
{ return SubSMA2CMSA(fafile, NULL, MA); }

cma_typ	SubSMA2CMSA(char *fafile, BooLean *remove, sma_typ MA)
{
	a_type	A;
	ss_type	data;

        A = MkAlpha(AMINO_ACIDS,PROT_BLOSUM62); data = SeqSet(fafile,A);
	return SubSMAtoCMSA(data, remove, MA);
}

cma_typ	SMAtoCMSA(ss_type data, sma_typ MA)
{ return SubSMAtoCMSA(data, NULL, MA); }

cma_typ	SubSMAtoCMSA(ss_type data, BooLean *remove, sma_typ MA)
{
        Int4    s,*len,t,n,ntyp,t2;
        st_type S;
        e_type  E;
        a_type  A;
        BooLean **null;
        char    *null0;
        cma_typ cmsa;

	A = SeqSetA(data);
	ntyp=ntypSMA(MA);
	if(remove != NULL){
	   for(t=1; t<= ntypSMA(MA); t++){
		if(remove[t]) { ntyp--; }
	   }
	}
        NEW(len,ntypSMA(MA)+2,Int4);
        NEWP(null,ntypSMA(MA)+2,BooLean);
        for(t2=0,t=1; t<= ntypSMA(MA); t++) {
	  if(remove == NULL || !remove[t]){
	    t2++;
            len[t2]=lengthSMA(t,MA);
            NEW(null[t2],len[t2]+2,BooLean);
            null0 = nullSMA(t,MA);
            for(n=1; n<=len[t2]; n++) {
                switch(null0[n-1]){
                 case '.': null[t2][n] = TRUE; break;
                 case '!': null[t2][n] = FALSE; break;
                 case '*': null[t2][n] = FALSE; break;
                 case '^': null[t2][n] = FALSE; break;
                 default: print_error("input error"); break;
                }
            }
	  }
        }
        S = MkSites(ntyp,len,data,10,5,5.0,9,9);
        for(n=1; n<=nseqSMA(MA); n++){
           E = SeqSetE(n,data);
           for(t2=0,t=1; t<= ntypSMA(MA); t++){
	     if(remove == NULL || !remove[t]){
		t2++;
                s = startSMA(t,n,MA);
                s -= OffSetSeq(E);
                AddSite(t2,n,s,S);
	     }
           }
        }
        cmsa=MakeCMSA(S, null);
        for(t=1; t<= ntyp; t++) free(null[t]);
        free(null); free(len);
        return cmsa;
}

