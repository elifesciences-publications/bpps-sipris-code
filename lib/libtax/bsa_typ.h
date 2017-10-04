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

#if !defined (_BSA_TYP_)
#define _BSA_TYP_

#include "cmsa.h"
#include "residues.h"
#include "histogram.h"
#include "swaln.h"
#include "dheap.h"
#include "mheap.h"
#include "dsets.h"
#include "tax_typ.h"
#include "set_typ.h"

class bsa_typ {         // build seed alignments
public:
        bsa_typ( ){ assert(!"Illegal constructor"); }
        bsa_typ(cma_typ in_cma, Int4 argc, char *argv[]) { Init(in_cma,argc,argv); }
	Int4    FindSeedAlns(Int4 iter,char *infilename,Int4 num_random, cma_typ cma);
	Int4    FindSeedAlns(Int4 iter,char *infilename,cma_typ cma)
		{ return FindSeedAlns(iter,infilename,1+ (N/3), cma); }
	void	DecrementIter(){ Iter--; }
	set_typ	RtnCopyTriedSeeds(){ return CopySet(TriedSeeds); }
	cma_typ	RtnCsqCMSA(char *name){
		    FILE *fp=tmpfile(); PutCsqCMSA_File(fp,name,input_cma); rewind(fp);
		    cma_typ csq_cma=ReadCMSA(fp,AB); fclose(fp); return csq_cma;
		}
	BooLean	IsEmpty() { return EmptyMheap(mH); }
	cma_typ GetSeedAln( );
	cma_typ GetSeedAln(Int4 N, set_typ &RtnSet){
		   RtnSet=MakeSet(N+1); ClearSet(RtnSet);
		   cma_typ rtn=GetSeedAln( );
		   Int4 M=SetN(DisplaySet);	// DisplaySet initialized by GetSeedAln();
		   for(Int4 i=1; i < M; i++){
			if(MemberSet(i,DisplaySet)) AddSet(i,RtnSet);
		   }
		   return rtn;
		}
        ~bsa_typ( ){ Free( ); }
private:
	BooLean	RandomlyInsert;
        //****************** Multiple category routines: ******************
        void    Init(cma_typ in_cma, Int4 argc, char *argv[]);
	void	Free();

	void    PutCsqCMSA_File(FILE *fp, cma_typ cma){ PutCsqCMSA_File(fp,NameCMSA(cma),cma); }
	void    PutCsqCMSA_File(FILE *fp, char *name, cma_typ cma);
	double  fastest_sq2sq_percent_identCMSA(register Int4 , 
			register unsigned char *, register unsigned char *);
	void    PutChnFileCMSA(FILE *fp, cma_typ cmaQ, char *main_set_file);
	void    PutChnFileCMSA(FILE *fp, cma_typ cmaQ, cma_typ cmaM);

	mh_type mH;
	Int4	ReloadHeap( );
	unsigned char **PtrFakeSq;	// ptr to fake sequences to speed up search.

	ds_type	sets;
	e_type	*SeqList;
	Int4	LinkDisjointSets(Int4 setI, Int4 setJ)
		{
		   assert(setI != setJ);
		   Int4 setIJ= linkDSets(setI,setJ,sets);
                   assert(setIJ == setI || setIJ == setJ);
                   e_type SeqI = SeqList[setI], SeqJ = SeqList[setJ];
                   Int4 len_list= ConcatSeqList(SeqI,SeqJ);     // attaches SeqJ list to end of SeqI list.
		   SeqList[setIJ]=SeqI;
		   // fprintf(stderr,"list length = %d\n",len_list); // len_list is for SeqI only!!!
		   return setIJ;
		}
	Int4    GetDisplaySet(Int4 setI,ss_type data);

	Int4	current_rank;
	Int4	Iter;
	Int4    *StoreI,*StoreJ,*Rank;

	e_type	*CSQ;
	BooLean	UsableCsq(e_type csq)	// does this consensus correspond to a usable subgroup?
		{
                   Int4 x,rx,n;
                   for(x=1; x < NumDisplaySets; x++){
                          for(n=0, rx = 1; rx <= LenSeq(csq); rx++) if(ResSeq(rx,csq) == ResSeq(rx,CSQ[x])) n++; 
                          if(SeedDiversity < ((double)n/ (double)LenSeq(csq))) return FALSE; 
                   } return TRUE;
		}

	// parameter settings:
	void	InitDefaults( ),InitAsNull( ),InitFlags( );
	double	SeedDiversity,MinSeqIdent;
	Int4	MaxNumDisplaySets;
	Int4	PrintSize,MaxNumPrint,Factor;
	BooLean	PrintLots;

	Int4	NumDisplaySets; 		// Keep track so can check for seed diversity.
	set_typ TriedSeeds,DisplaySet,SeedSet;
	cma_typ	CreateDisplayCMA(Int4 NumPrinted, cma_typ cma);

	char	*main_set_file;
	Int4	N,time1;
        cma_typ input_cma;
        a_type  AB;
};

#endif

