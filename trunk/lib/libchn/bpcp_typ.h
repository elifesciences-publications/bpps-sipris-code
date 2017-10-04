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

#if !defined (_BPCP_TYP_)
#define _BPCP_TYP_

#include "set_typ.h"
#include "cmsa.h"
#include "hpt_typ.h"

class bpcp_typ {         // Bayesian partitioning checkpoint type
public:
	bpcp_typ( ){ Init(); }
	bpcp_typ(Int4 sd, Int4 ns,set_typ *st, Int4 n, cma_typ *s, hpt_typ *h){
		Init(sd,ns,st,n,s,h); 
	}
	~bpcp_typ( ){ Free(); }
private:
	void Free(){ }
	Int4	seed;
	set_typ	*sets;
	Int4	NumSets;
	hpt_typ	*hpt;
	cma_typ	*sma;
	Int4	NumSMA;
	void	Init(Int4 sd=0, Int4 ns=0, set_typ *st=0, Int4 n=0,cma_typ *s=0, hpt_typ *h=0){
		seed=sd; sets=st; sma=s; hpt=h; NumSMA=n; NumSets=ns;
		if(n > 0){  // then check input for consistency.
		   Int4	*P; assert(hpt->IsTree(P)); free(P);
		   assert(NumSets==hpt->NumSets());
		   assert(NumSets==NumSMA+1);
		}
	}
public:
	void Put(char *outfile){
		FILE *fp=open_file(outfile,".seed","w"); fprintf(fp,"-seed=%d",seed); fclose(fp);
		fp=open_file(outfile,".sma","w");
    		for(Int4 i=1; i<=NumSMA; i++){ PutCMSA(fp,sma[i]); } fclose(fp);
		fp=open_file(outfile,".sets","w"); WriteSets(fp,NumSets,sets); fclose(fp);
		fp=open_file(outfile,".hpt","w"); hpt->Put(fp); fclose(fp);
	}
	void Write(FILE *fp){
          assert(fp != NULL);
	  assert(NumSMA != 0);
          WriteSets(fp,NumSets,sets);
          fprintf(fp,"seed=%d; N=%d.\n",seed,NumSMA);
          for(Int4 i=1; i <= NumSMA; i++){ PutCMSA(fp,sma[i]); }
          hpt->Put(fp,TRUE,FALSE,TRUE); 
	}
	cma_typ *RtnSMA(Int4 &N){ N=NumSMA; return sma; }
	hpt_typ *RtnHpt(){ return hpt; }
	Int4	RtnSeed(){ return seed; }
	set_typ *RtnSets(Int4 &n){ n=NumSets; return sets; }
	void Read(FILE *fp, a_type AB){	
          assert(fp != NULL);
	  assert(NumSMA == 0);
	  char str[202];
          fprintf(stderr,"Reading checkpoint file\n");
          sets=ReadSets(fp,NumSets);
          if(fgets(str,200,fp) == NULL) print_error("*.chk file input error 1.");
          if(sscanf(str,"seed=%d; N=%d.",&seed,&NumSMA) != 2) print_error("*.chk file input error 2.");
          NEW(sma, NumSMA+3, cma_typ);
          for(Int4 i=1; i<=NumSMA; i++){
               if((sma[i]=ReadCMSA(fp,AB))==NULL) print_error("*.chk file input error 3.");
          }
          hpt=new hpt_typ(fp);
          if(NumSets != hpt->NumSets()) print_error("FATAL: *.hpt and *.sets files are inconsistent");
	}
};

#endif

