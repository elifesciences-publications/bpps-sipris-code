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

#if !defined (_GSM_TYP_)
#define _GSM_TYP_
#include "gibbs.h"
#include "binomial.h"
#include "purge.h"
#include "goscan.h"
#include "msaheap.h"
#include "set_typ.h"
#include "oscan2msa.h"
#include "watchgibbs.h"
#include "cluster.h"
#include "blosum62.h"
#include "cma_gmb.h"
#include "ppg_typ.h"

class gsm_typ {
public: 
	gsm_typ(int argc, char *argv[], a_type);
	gsm_typ(char *USAGE,int argc, char *argv[], a_type);
	gsm_typ(char *USAGE,int argc, char *argv[], a_type, ss_type);
	BooLean	use_msa(){ return input_msa; } // start with an msa search.
	BooLean	Gapped(){ return use_gseq; } // use gapped sequences. 
	BooLean	doBreed(){ return use_breed; } // apply a genetic algorithm.
	BooLean	GibbsSearch(){ return (cma_purge > 0); } // optimize using cma_alignment.
	void	Gibbs(cma_typ cma){ gibbs(cma); optimize(); } // optimize using cma_alignment.
	void	RefineInputAln(cma_typ cma);	// assume MAPGAPS created...
	cma_typ	optbreed(cma_typ); 
	void	purge(); 	// purge database 

	BooLean	doImprove(){ return improve; } // start with an msa search
	cma_typ	Improve(cma_typ); 
	cma_typ	BestCMA( ){ return bestmsa; } 
	Int4	Maxrun(){ return maxrun; }
	Int4	AlignMode(){ return align_mode; }
	BooLean	Go(){ return go; }
	~gsm_typ();
	//********************* gsm_plus.cc *************************
	Int4	*StrtLens,StrtBlks;
	void	SetStartBlks(Int4 b, Int4 *L){
			assert(b > 0 && L != 0); StrtBlks=b; 
			if(StrtLens) free(StrtLens); NEW(StrtLens, b+3, Int4);
			for(b=1; b <= StrtBlks; b++) StrtLens[b]=L[b];
		}
	//********************* gsm_plus.cc *************************
	void	SetMaxIter(Int4 x){ assert(x > 0); MaxIter=x; }
	void	SetTestMode(){ TestMode=TRUE; }
	BooLean	TestMode;
	void	DoWrite(){ write=TRUE; }
	Int4	RtnMinPurgedSetSize(){ return minseq; }
private:
	BooLean	write;
	Int4	MaxIter;
	//********************* gsm_init.cc *************************
	void	Init(char *USAGE,int argc, char *argv[], a_type ab, ss_type D=0);
	void	make_binomials(ss_type);
	void	read_input_data();
	void	free_input_data();
	void    InitDefaults();
	void    InitAsNull();
	//********************* gsm_init.cc *************************

	//********************* gsm_srch.cc *************************
	static const Int4 MAX_IN_SEQS=20000000;
public: 
	void	ResetParameters(Int4 b,Int4 c, Int4 mnb,Int4 mxb){
		   assert(b > 0); assert(c > 2); assert(mnb > 0); assert(mnb <= mxb); 
			aveblk=b; avecol=c; minblk=mnb; maxblk=mxb; fix=TRUE;
		}
	void	search();  // scan search
	e_type	**search(Int4,Int4,Int4,BooLean,char *); 	// scan search
	void	gapped_srch();	// new 7/8/09 (AFN)
	BooLean	doGappedSrch(){ return improve; } // start with an msa search
private:
	BooLean	do_gapped_srch; 		// NEW Gapped search options 7/8/09 (AFN)
	Int4	Gap_o,Gap_x;
	// void	setup_next_srch(); not defined or used...
	//********************* gsm_srch.cc *************************

	//********************* gsm_smpl.cc *************************
	double  the_gibbs_sampler(gs_type G,FILE *rfp=0);
	// ^ modified version of what is in libaln/gibbs.cc
	BooLean	SmplBlksCols;
	//********************* gsm_smpl.cc *************************

	//********************* gsm_typ.cc *************************
public: 
#if 1	// This ignores block combinatorics when computing LLR (used as key for maH).
	double  LogLikeRatio(cma_typ cma){ return log_like_ratio(cma); }
#else	// This includes block combinatorics.
	double  LogLikeRatio(cma_typ cma){ return RelMapCMSA(cma); }
#endif
	double  log_like_ratio(cma_typ cma);
	Int4	align(); 	// create and optimize alignment
	cma_typ	*Align(Int4 b,Int4 c, Int4 mnb,Int4 mxb,Int4 len_mode); // for gismo++
	Int4	Align2(Int4 b,Int4 c, Int4 mnb,Int4 mxb,Int4 len_mode); // for gismo++
	void	gapped_align();	// create and optimize gapped alignment
	cma_typ	*RtnBestCMA();

	// double	RunPublicGibbs(cma_typ *M){ strcpy(options,"-t1 "); return core_gibbs(M,'N',-1.0,0.0); }
private:
	// customized gibbs sampler:
	// double	sim_anneal_gibbs(cma_typ *M, char mode, double temp) { return core_gibbs(M,mode,temp,0.0); }
	double  core_gibbs(cma_typ *M, char mode='N', double temp=-1.0, double minmap=0.0);
public:
	//=====================
	cma_typ	AddEndBlocks(cma_typ icma){
	   double lpr0,lpr,d,dd;
	   Int4	minleng=5,x,end,i;
	   cma_typ rcma=0,xcma=icma,bcma=0;
	   lpr0=RelMapCMSA(xcma);
	   for(i=1,x=nBlksCMSA(xcma); i <= 2; i++,x=0){
	     rcma=AddBlkCMSA(x,minleng,xcma);
	     if(rcma){
                if(x==0) fprintf(stderr,"N-terminal block added.\n");
		else fprintf(stderr,"C-terminal block added.\n");
		this->use_gseq=TRUE; this->gapopen=1000; this->gapextend=100;
           	this->pernats=1000;
		SetPenaltyCMSA(gapopen,gapextend,rcma);
           	cma_typ tcma=this->cmsa; this->cmsa=0;
	        // sprintf(options,"-t1 -g -l%d ",limit);
	        sprintf(options,"-t10 -g -l%d ",20);
	        lpr = this->core_gibbs(&rcma,'D',300,0.0);  this->cmsa=tcma;
		d=FieldRelMapCMSA(rcma,x+1);	// new block = 1 or = nBlk+1.
		if(d <= 0.0){
	            sprintf(options,"-t10 -g -l%d ",20);
	            lpr = this->core_gibbs(&rcma,'S',100,0.0);  this->cmsa=tcma;
		    d=FieldRelMapCMSA(rcma,x+1);
		}
	        // if(lpr <= lpr0 || d <= 0.0){ NilCMSA(rcma); rcma=0; }
	        if(lpr <= lpr0){ NilCMSA(rcma); rcma=0; }
		else { if(bcma) NilCMSA(bcma); xcma=bcma=rcma; rcma=0; lpr0=lpr; }
             }
	   } return bcma;
	}
	BooLean	AddEndBlks(cma_typ &icma){
	   cma_typ rcma=0,xcma=icma;
	   double d,dd;
           while(rcma=this->AddEndBlocks(xcma)){
             if(xcma != icma) NilCMSA(xcma); xcma=rcma; rcma=0;
           } 
           if(xcma != icma){
#if 0
		x=nBlksCMSA(xcma); d=FieldRelMapCMSA(xcma,x);
		if(d <= 0.0){
		   rcma=DeleteBlkCMSA(x,xcma); NilCMSA(xcma); xcma=rcma; rcma=0;
		}
#endif
		NilCMSA(icma); icma=xcma; return TRUE; 
	   } else return FALSE;
	}
private:
        cma_typ bestmsa,cmsa,cmsa_in;
	BooLean	gibbs(cma_typ ma2,char tweakmode='S');
	BooLean	use_gseq,dont_delete;
	Int4	gapopen,gapextend;
	double	indel_penalty,pernats;

	Int4	SampleConfig();
	BooLean	create_align(); 	
	Int4	create_cma(Int4 &num_success, Int4 max_futile=20);
	void	Optimize();
	BooLean	Record( );
	BooLean	Refine();
	BooLean	SplitUp(Int4);	// split up a raw alignment from MAPGAPS.

	Int4	MinSeqLen;
	double	trim_info;

	void	RefineMSAHeap();
	//********************* gsm_typ.cc *************************

	//********************* gsm_blks.cc *************************
public:
	BooLean	DeletePoor(char tweakmode='N');
private:
	BooLean	Add(Int4);
	BooLean	Delete(char tweakmode='N');
	BooLean	Split(Int4);
	BooLean	Fuse(Int4);
	double	Breed();
	void	create_population();
	void	Recombine();
	void	optimize();	// with block operations.

	BooLean	SFR(Int4,Int4); 
	BooLean SFL(Int4,Int4); 
	BooLean FSR(Int4);
	BooLean FSL(Int4);
	BooLean SSF(Int4, Int4); // NOT VERY USEFUL...
	BooLean SFF(Int4, Int4); // NOT VERY USEFUL...
	//********************* gsm_blks.cc *************************

	//********************* parameters **************************
	Int4	COL,BLK;
	char	*Argv[3],method,mode,*guide,aafreq;

	BooLean	use_breed,input_msa,weight,combine,repeats,force,fix,go,flag,noseg;
        BooLean report_gaps,useSA,UseRepset,UseLabeled,mask_nonglobular,improve;

	double	improve_cut;
	Int4	max_motif_width;
	Int4	*num,number,item,*counts,Run,maxbreed,align_mode;
        Int4    avecol,aveblk,minseq,cycle,maxlenModels;
        Int4    maxseq,inc,maxrun,cutoff,over,under;
        Int4    maxcycle,mhpsz,maxrpts,maxseqleng; // for gapped scan
        Int4    left_flank,right_flank,minblk,maxblk;
        Int4    cardB,oldcardB;	// number of sequence detected
        char    str[300],*name,options[300];
	Int4	limit;
        float   minmap;
        double  map,bestmap,pseudo,Ecut,ecut,temperature;
        double  *freq,repeatEval,breedVar;
        UInt4	seed;
	time_t	time1;
        UInt4 min_rpt;
        unsigned short  *nsize;
        wg_type WG;		// WatchGibbs object
        a_type  AB;
        ss_type data,in_data;
	BooLean	input_data;
        mah_typ maH;
        set_typ DBS_hits;	// keep track of database hits
        gd_type Guide;
	Int4	patience;	// how Int4 am I willing to wait for an improvement?
	Int4	cma_purge;	// realignment purge cutoff in percent identity.

	Int4	MinGoodBlks,MaxBadBlks;
	// gss_typ	gss2,gss;
	gss_typ	*gss;
};

#endif


