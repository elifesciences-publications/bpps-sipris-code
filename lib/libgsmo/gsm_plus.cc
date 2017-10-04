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

#include "gsm_typ.h"
#include "gmb_typ.h"
#include "cma_gmb.h"
#include "gpr_typ.h"

void	GetParameters(Int4 stage, double &bild_cut, double &bild_add, double &Temp,
		double &afrq, double &rfrq, Int4 &aa_per_io, Int4 &aa_per_do,
		Int4 &exp_ie, Int4 &exp_de, double &prior_wt,char &dms_mode);

cma_typ	*RtnBestCore(FILE *rfp, cma_typ	*rtn_cma,Int4 StopHere)
{
        cma_typ *one_cma,rcma;
	double	shared[20][20],CL,scr,bst_scr;
	Int4	i,j,k,bst_core,num_shared;
        cma_typ corecma[20][20],*cor_cma,*strt_cma=0;
        for(num_shared=0,i=1; rtn_cma[i]; i++) {     // convert to one block.
		if(rtn_cma[i]) num_shared++;
	} assert(num_shared < 20);
        NEW(one_cma,num_shared + 3, cma_typ);
        for(i=1; i <= num_shared; i++) {     // convert to one block.
	  assert(rtn_cma[i]);
	  if(rfp){
              fprintf(rfp,"%d: GISMO LLR = %0.2f; (%d blocks; %d columns).\n   ",
                i,UnGappedRelMapCMSA(rtn_cma[i]),nBlksCMSA(rtn_cma[i]),NumColumnsCMSA(rtn_cma[i]));
	      PutConfigCMSA(rfp,rtn_cma[i]); fflush(rfp);
	  } for(j=1; j <= num_shared; j++){ shared[i][j]=0; corecma[i][j]=0; }
#if 0	// Copying cma is not creating gss type properly!
	  if(nBlksCMSA(rtn_cma[i]) > 1) rcma=OneBlockCMSA(rtn_cma[i]); 
	  else rcma=CopyCMSA(rtn_cma[i]);
#else
	  rcma=OneBlockCMSA(rtn_cma[i]); 
#endif
	  ExtendFakeToRealCMSA(rcma);

#if 1	// Do additional sampling to introduce gaps.
	  double	bc,ba,T,afq,rfq,prior_wt,MaxFrctn=0.05;
	  Int4		aa_per_io,aa_per_do,exp_ie,exp_de,MinSticky=2; 
	  char		dms_mode;
	  GetParameters(0,bc,ba,T,afq,rfq,aa_per_io,aa_per_do,exp_ie,exp_de,prior_wt,dms_mode);
	  gmb_typ *gmb = new gmb_typ(aa_per_io,aa_per_do,exp_ie,exp_de,rcma,dms_mode,prior_wt);
	  fprintf(rfp," "); gmb->SampleStickyTogether(rfp, MinSticky, MaxFrctn,0,200,'R');
	  delete gmb;
#endif
	  one_cma[i]=rcma; 
        }

        //================= Use best common core alignment. ===============
        for(i=1; i < num_shared; i++) {
          assert(rtn_cma[i]); assert(one_cma[i]);
          for(j=i+1; j <= num_shared; j++) {
                assert(rtn_cma[j]); assert(one_cma[j]);
                double  frq_cut=0.50,ave_frq=0;
		cor_cma=GetCoreCMSA(one_cma[i],one_cma[j],ave_frq,frq_cut);
		if(cor_cma){
                   corecma[i][j]=cor_cma[1]; corecma[j][i]=cor_cma[2]; free(cor_cma);
		  // fprintf(stderr,"%d.%d: %g\n",i,j,ave_frq);
                  // CL=(double) LengthCMSA(1,corecma[i][j][1]) + ave_frq;
                  CL=(double) LengthCMSA(1,corecma[i][j])*ave_frq;
                  // CL=ave_frq;
		  // WARNING: Important to take shared score over all alignments!!!
                  shared[i][i] += CL; shared[j][j] += CL;
                  shared[i][j]=shared[j][i] = CL;
		} else { shared[i][j]=shared[j][i] = 0.0; }
          }
        } NEW(strt_cma, num_shared+3,cma_typ);
	dh_type dH=0;
        for(bst_core=0,bst_scr=0,i=1; i < num_shared; i++){
          for(j=i+1; j <=num_shared; j++) {
		if(corecma[i][j] == 0) continue;
                if((CL=shared[i][j]) > bst_scr){
		   bst_scr=CL;
                   if(shared[i][i] >= shared[j][j]){ bst_core=i; } else { bst_core=j; }
                } 
          }
        } fprintf(rfp,"\n best core = %d (%.3f score; %.2f aggregate).\n   :",
                                        bst_core,bst_scr,shared[bst_core][bst_core]);
        for(j=1; j <= num_shared; j++) fprintf(rfp," %8d",j); fprintf(rfp,"\n");
	cma_typ Rtn_cma[20];
	dH=dheap(num_shared+4,3);
        for(i=1; i <= num_shared; i++) {
                fprintf(rfp,"%3d:",i);
                for(bst_core=0,bst_scr=0,j=1; j <= num_shared; j++){
		   fprintf(rfp," %8.2f",shared[i][j]); 
		   if(i != j && shared[i][j] > bst_scr){ bst_core=j; bst_scr=shared[i][j]; }
		}
                fprintf(rfp," (%d cols; %d blks; %d best)\n",
                        NumColumnsCMSA(rtn_cma[i]),nBlksCMSA(rtn_cma[i]),bst_core);
#if 1
		CL=shared[i][i];	// this setting seems to work best; was used for FprH7 run...
#else	// This appears to degrade the alignment quality.
		CL=shared[i][i] + NumColumnsCMSA(rtn_cma[i]);   // favor more columns to avoid truncations?
#endif
		insrtHeap(i,(keytyp)-CL,dH);
		Rtn_cma[i]=rtn_cma[i]; rtn_cma[i]=0;
#if 0	// return core only...this doesn't seem to help.
                if(one_cma[i]) NilCMSA(one_cma[i]); one_cma[i]=0;
		if(bst_core > 0){
		   one_cma[i]=corecma[i][bst_core]; corecma[i][bst_core]=0; 
		   assert(one_cma[i]);
		}
#endif
        } fprintf(rfp,"\n"); fflush(rfp); // free(one_cma); // free(rtn_cma); 
	fprintf(rfp,"best seed alns:");
	for(j=1; (i=delminHeap(dH)) != 0; j++){
		if(j > StopHere){
                   if(Rtn_cma[i]) NilCMSA(Rtn_cma[i]); Rtn_cma[i]=0;
                   // if(rtn_cma[i]) NilCMSA(rtn_cma[i]); rtn_cma[i]=0;
                   if(one_cma[i]) NilCMSA(one_cma[i]); one_cma[i]=0;
		} else {
		   fprintf(rfp," %d(%.2f)",i,(shared[i][i]+(double) NumColumnsCMSA(Rtn_cma[i])));
		   rtn_cma[j]=Rtn_cma[i]; strt_cma[j]=one_cma[i]; 
		}
	} Nildheap(dH); fprintf(rfp,"\n\n"); free(one_cma);
        for(i=1; i <= num_shared; i++) {
           for(j=1; j <= num_shared; j++) if(corecma[i][j]) NilCMSA(corecma[i][j]); 
	}
// dsw=rssx->RtnWtNumSeqs();
	return strt_cma;
}


// UnLabelSeq(e_type E); LabelSeq(e_type E); LabeledSeq(E);

// GetParameters(stage,bild_cut,bild_add,Temp,afrq,rfrq,aa_per_io,aa_per_do,exp_ie,exp_de,prior_wt);

void	GetParameters(Int4 stage, double &bild_cut, double &bild_add, double &Temp,
		double &afrq, double &rfrq, Int4 &aa_per_io, Int4 &aa_per_do,
		Int4 &exp_ie, Int4 &exp_de, double &prior_wt,char &dms_mode) 
{
	switch(stage){	// "TfF" works  better than "OTf"...
#if 0	// FprJ4[a-h]f run...
	     case 0: bild_cut=-4.2; bild_add=-4.0; Temp=200; dms_mode='T'; 
	  	  aa_per_io=20; aa_per_do=150; exp_ie=1;exp_de=1;
		  // add columns with > 30 % matches; rm with > 40% deletions (but not used).
		  afrq=0.30; rfrq=0.40; prior_wt=0.90; break;

	     case 1: bild_cut=-3.7; bild_add=-3.5; Temp=100; dms_mode='T'; 
	  	  aa_per_io=20; aa_per_do=150; exp_ie=1;exp_de=1;
		  // add columns with > 30% matches; rm with > 40% deletions.
		  afrq=0.30; rfrq=0.40; prior_wt=1.00; break;

	     case 2: bild_cut=-3.2; bild_add=-3.0; Temp=0; dms_mode='f'; 
	  	  aa_per_io=20; aa_per_do=170; exp_ie=1;exp_de=1;
		  // add columns with > 50% matches; rm with > 52% deletions.
		  afrq=0.50; rfrq=0.52; prior_wt=1.00; break;

	     case 3: bild_cut=-1.7; bild_add=-1.5; Temp=0; dms_mode='F'; 
	  	  aa_per_io=20; aa_per_do=200; exp_ie=1;exp_de=1;
		  // add columns with > 60% matches; rm with > 42% deletions.
		  afrq=0.60; rfrq=0.42; prior_wt=1.00; break;

	     case 4: bild_cut=0.0; bild_add=0.0; Temp=0; dms_mode='F'; 
	  	  aa_per_io=20; aa_per_do=200; exp_ie=1;exp_de=1;
		  afrq=0.60; rfrq=0.45; prior_wt=1.00; break;

#elif 0		// FprK7[t-z]f
	     case 0: bild_cut=-4.7; bild_add=-4.5; Temp=200; dms_mode='T'; 
	  	  aa_per_io=20; aa_per_do=150; exp_ie=1;exp_de=1;
		  // add columns with > 30 % matches; rm with > 40% deletions (but not used).
		  afrq=0.25; rfrq=0.35; prior_wt=0.90; break;

	     case 1: bild_cut=-4.2; bild_add=-4.0; Temp=100; dms_mode='T'; 
	  	  aa_per_io=20; aa_per_do=150; exp_ie=1;exp_de=1;
		  // add columns with > 30% matches; rm with > 40% deletions.
		  afrq=0.30; rfrq=0.40; prior_wt=1.00; break;

	     case 2: bild_cut=-2.2; bild_add=-2.0; Temp=0; dms_mode='f'; 
	  	  aa_per_io=20; aa_per_do=170; exp_ie=1;exp_de=1;
		  afrq=0.50; rfrq=0.52; prior_wt=1.00; break;

	     case 3: bild_cut=-0.2; bild_add=0.0; Temp=0; dms_mode='F'; 
	  	  aa_per_io=20; aa_per_do=200; exp_ie=1;exp_de=1;
		  afrq=0.50; rfrq=0.52; prior_wt=1.00; break;

	     case 4: bild_cut=0.2; bild_add=0.2; Temp=0; dms_mode='F'; 
	  	  aa_per_io=20; aa_per_do=200; exp_ie=1;exp_de=1;
		  afrq=0.50; rfrq=0.52; prior_wt=1.00; break;

#elif 1		// FprJ9[t-z]f
	     case 0: bild_cut=-4.7; bild_add=-4.5; Temp=200; dms_mode='T'; 
	  	  aa_per_io=20; aa_per_do=120; exp_ie=1;exp_de=1;
		  // add columns with > 30 % matches; rm with > 40% deletions (but not used).
		  afrq=0.25; rfrq=0.35; prior_wt=0.90; break;

	     case 1: bild_cut=-4.2; bild_add=-4.0; Temp=100; dms_mode='T'; 
	  	  aa_per_io=20; aa_per_do=120; exp_ie=1;exp_de=1;
		  // add columns with > 30% matches; rm with > 40% deletions.
		  afrq=0.30; rfrq=0.40; prior_wt=1.00; break;

	     case 2: bild_cut=-2.7; bild_add=-2.5; Temp=0; dms_mode='f'; 
	  	  aa_per_io=20; aa_per_do=130; exp_ie=1;exp_de=1;
		  afrq=0.50; rfrq=0.52; prior_wt=1.00; break;

	     case 3: bild_cut=-1.2; bild_add=-1.0; Temp=0; dms_mode='F'; 
	  	  aa_per_io=20; aa_per_do=150; exp_ie=1;exp_de=1;
		  afrq=0.50; rfrq=0.52; prior_wt=1.00; break;

#if 0
	     case 4: bild_cut=0.2; bild_add=0.2; Temp=0; dms_mode='F'; 
	  	  aa_per_io=20; aa_per_do=150; exp_ie=1;exp_de=1;
		  afrq=0.50; rfrq=0.52; prior_wt=1.00; break;
#else
	     case 4: bild_cut=-0.2; bild_add=0.0; Temp=0; dms_mode='F'; 
	  	  aa_per_io=20; aa_per_do=150; exp_ie=1;exp_de=1;
		  afrq=0.50; rfrq=0.52; prior_wt=1.00; break;
#endif

#elif 1		// FprJ9[t-z]f
	     case 0: bild_cut=-4.7; bild_add=-4.5; Temp=200; dms_mode='T'; 
	  	  aa_per_io=20; aa_per_do=150; exp_ie=1;exp_de=1;
		  // add columns with > 30 % matches; rm with > 40% deletions (but not used).
		  afrq=0.25; rfrq=0.35; prior_wt=0.90; break;

	     case 1: bild_cut=-4.2; bild_add=-4.0; Temp=100; dms_mode='T'; 
	  	  aa_per_io=20; aa_per_do=150; exp_ie=1;exp_de=1;
		  // add columns with > 30% matches; rm with > 40% deletions.
		  afrq=0.30; rfrq=0.40; prior_wt=1.00; break;

	     case 2: bild_cut=-3.2; bild_add=-3.0; Temp=0; dms_mode='f'; 
	  	  aa_per_io=20; aa_per_do=170; exp_ie=1;exp_de=1;
		  afrq=0.50; rfrq=0.52; prior_wt=1.00; break;

	     case 3: bild_cut=-1.7; bild_add=-1.5; Temp=0; dms_mode='F'; 
	  	  aa_per_io=20; aa_per_do=200; exp_ie=1;exp_de=1;
		  afrq=0.60; rfrq=0.42; prior_wt=1.00; break;

	     case 4: bild_cut=0.0; bild_add=0.0; Temp=0; dms_mode='F'; 
	  	  aa_per_io=20; aa_per_do=200; exp_ie=1;exp_de=1;
		  afrq=0.60; rfrq=0.45; prior_wt=1.00; break;
#else
	     case 0: bild_cut=-4.7; bild_add=-4.5; Temp=200; dms_mode='T'; 
	  	  aa_per_io=20; aa_per_do=120; exp_ie=5; exp_de=1; // for FprJ4[t-z]f
	  	  // aa_per_io=20; aa_per_do=120; exp_ie=1; exp_de=1;    // for FprJ5[a-h]f
		  // add columns with > 30 % matches; rm with > 40% deletions (but not used).
		  // afrq=0.30; rfrq=0.40; prior_wt=0.60; break;
		  afrq=0.25; rfrq=0.35; prior_wt=0.60; break;

	     case 1: bild_cut=-4.2; bild_add=-4.0; Temp=100; dms_mode='T'; 
	  	  aa_per_io=20; aa_per_do=120; exp_ie=1;exp_de=1;
		  // add columns with > 30% matches; rm with > 40% deletions.
		  afrq=0.30; rfrq=0.40; prior_wt=0.90; break;

	     case 2: bild_cut=-3.2; bild_add=-3.0; Temp=0; dms_mode='f'; 
	  	  aa_per_io=20; aa_per_do=150; exp_ie=1;exp_de=1;
		  afrq=0.50; rfrq=0.50; prior_wt=1.00; break;

	     case 3: bild_cut=-1.7; bild_add=-1.5; Temp=0; dms_mode='F'; 
	  	  aa_per_io=20; aa_per_do=150; exp_ie=1;exp_de=1;
		  afrq=0.60; rfrq=0.42; prior_wt=1.00; break;

	     case 4: bild_cut=0.0; bild_add=0.0; Temp=0; dms_mode='F'; 
	  	  aa_per_io=20; aa_per_do=150; exp_ie=1;exp_de=1;
		  afrq=0.70; rfrq=0.35; prior_wt=1.00; break;
#endif

	     default: print_error("run_gambit() this should not happen"); break;
	}
}

cma_typ	run_gambit(FILE *rfp,char *name, Int4 iter, ssx_typ *issx, Int4 StageStart, Int4 Stages)
{
	Int4	time1=time(NULL);
	cma_typ	in_cma=issx->RtnCMA();
	char	dms_mode=' ',str[100];
	double	dd,d,D,prior_wt=0,pn;
	Int4    i,j,x,similarity=0,aa_per_io=0,aa_per_do=0,exp_ie=0,exp_de=0;
	issx->GetParameters(aa_per_io,aa_per_do,exp_ie,exp_de,prior_wt,dms_mode);
	Int4	FloorI2I=issx->RtnFloorI2I(), FloorD2D=issx->RtnFloorD2D();
	// sets transition probability priors...
	Int4	MinSticky=2,stage;
	char	StMd=' ';
	
	a_type	AB=AlphabetCMSA(in_cma);
	cma_typ	final_cma,rcma,ocma,cma,xcma,bst_cma=0;

	Int4 timeR,timeI=time(NULL); cma=CopyCMSA(in_cma);
	//============ 3c. Optimize the MSA using GAMBIT =============
	assert(exp_de > 0); sprintf(str,"%s_gmb_tmp",name); 
	// double *xfq = &sfrq[0]; // can change for different stages...

	//============ 3d. Add & remove columns & realign sequences using GAMBIT =============
	char	change=' ';
	BooLean	keep_going,did_rm=TRUE,did_add=FALSE,did_mv=FALSE;
	double bild_cut=-2.0,bild_add=-2.0,Temp=300,afrq=0.50,rfrq=0.50;
	// Using dmp='f' or 'F' takes about twice as long versus 'T'.

	if(Stages > 4) print_error("run_gambit() Stages input error");
	//   afrq=0.60 --> add columns with > 60% matches...
	//   rfrq=0.40 --> remove columns with more than 40% deletions.
        for(stage=StageStart; stage <= Stages; stage++){
	  timeR = time(NULL); 
	  GetParameters(stage,bild_cut,bild_add,Temp,afrq,rfrq,
				aa_per_io,aa_per_do,exp_ie,exp_de,prior_wt,dms_mode);
	  //****************************************************************************
	  fprintf(rfp,"   ----- stage %d ----- (dmp: '%c').\n",stage,dms_mode); 
	  fprintf(rfp,"   GAMBIT: expected io=%d, do=%d, ie=%d, de=%d; prior_wt=%.2f; floor I2I=%d; D2D=%d.\n",
			aa_per_io,aa_per_do,exp_ie,exp_de,prior_wt,FloorI2I,FloorD2D);
	  {
	      Int4	nAddRms,nAddRmMv=0;
	      BooLean	did_mv;
	      gpr_typ gpr(name,iter,stage);
	      gpr.SetParameters(aa_per_io,aa_per_do,exp_ie,exp_de,prior_wt,dms_mode);
	      gpr.SetAddRmCols(bild_cut,bild_add,Temp,afrq,rfrq); gpr.SetSticky(MinSticky,StMd);
	      gpr.DoNotWrite();
	      if(stage == 0 && Stages==0){
		gpr.DoAddColumns(rfp,cma); continue; // first add columns (for initial purging).
	      } 
	      if(stage == 1){	// purged set...
		did_mv=gpr.DoMvColumns(rfp,cma); 
		if(did_mv) nAddRmMv++;
		nAddRms=gpr.DoAddRmColumns(rfp,cma,1);        // need to first add columns (for 1 cycle).
		if(nAddRms==0){
		   gpr.DoStickySample(rfp,cma);
		   nAddRms=gpr.DoAddRmColumns(rfp,cma,1);
		   if(nAddRms > 0) nAddRmMv += nAddRms;
		} else nAddRmMv += nAddRms; 
		gpr.DoPurgeSampling(rfp,cma); // gpr.DoSingleSample(rfp,cma); 
	      } else if(stage > 1){	// do presampling...Skipping this helps for PH domain!!!
		Int4 maxRnds=3;
		if(stage == 3) gpr.DoStickySample(rfp,cma);
		did_mv=gpr.DoMvColumns(rfp,cma); 
		if(did_mv) nAddRmMv++;
		for(did_mv=TRUE,i=1; (did_mv && i <= 3); i++){
		    nAddRms=gpr.DoAddRmColumns(rfp,cma,maxRnds); 
		    if(nAddRms > 0){
			nAddRmMv += nAddRms; did_mv=gpr.DoMvColumns(rfp,cma); 
			if(did_mv) nAddRmMv++; 
		    } maxRnds=1;
		} if(stage <= 3){
		   if(gpr.DoPurgeSampling(rfp,cma) > 0 && stage == 4) gpr.DoStickySample(rfp,cma);
		}
		if(stage == 4){ gpr.DoSingleSample(rfp,cma); }
	      } else {	// stage == 0;
		  gpr.DoWriteX(cma); 
		  // gpr.DoStickySample(rfp,cma);  // AddColumns does sticky sampling already...
		  nAddRmMv=gpr.DoAddColumns(rfp,cma); 
		  gpr.DoPurgeSampling(rfp,cma); 
	      }
	      if(nAddRmMv == 0) gpr.DoStickySample(rfp,cma);
          } if(stage > 0) fprintf(rfp," time: %d secs (%0.2f mins)\n\n",
				time(NULL)-timeR,(float)(time(NULL)-timeR)/60.0); fflush(rfp);
        } return cma;
}

static char GISMO_USAGE[]= "\n\
 \n\
   GISMO version 2.0 (9-8-2015)\n\
   for accurate alignment of large numbers of diverse protein sequences\n\
\n\
   Copyright (c) 2016 The University of Maryland and Andrew F. Neuwald \n\
   http://gismo.igs.umaryland.edu \n\
\n\n Usage: gismo fastafile [options] \n\
\n\
  options:\n\
\t-fast         - Run the sampler in fast mode (h=H=stages=2)\n\
\t-h=<int>      - Number of initial block-based alignments = N (2 to 18: default 10)\n\
\t-H=<int>      - Number of candidate HMM alignments (N to 9: default 5)\n\
\t-L            - Mask low complexity regions\n\
\t-maxseq=<int> - Set the maximum number of input sequences (default: 20,000)\n\
\t-rtf          - Also output a MS-Word-viewable rich text file (*.rtf)\n\
\t-seed=<int>   - Random seed\n\
\t-stages=<int> - Specify the number of sampling stages (2-4: default 4)\n\
\t                    \n";

#if 0
\t-C           - sample on closely-related clusters\n\
\t-C<int>      - sample on Clusters with \% similarity cutoff = <int>\n\
\t-p=<int1>:<int2>:<int3>:<int4>  - priors for computing transition probabilities\n\
                  int1 = average residue spacing between insertions (within blocks) (default: 200)\n\
                  int2 = average length of insertion extensions (default: 5)\n\
                  int3 = average residue spacing between deletions (default: 200)\n\
                  int4 = average length of deletion extensions (default: 2)\n\
\t-n            - mask potential nonglobular regions using 'seg x 45 3.4 3.75'\n\
\t-wt=<real)   - Set the weight to place on prior indel transition probs. (default: 1*data)\n\
\t-W           - output a *.see file to plot progress of sampler\n\
\t-strategy=<int>  - choose a sampling strategy for testing (range: 2-4 stages)\n\
\t                    \n";
#endif

cma_typ	run_gismo_plus(int argc, char *argv[],a_type AB)
{
	// Int4    i,j,x,arg,similarity=0,aa_per_io=35,aa_per_do=170,exp_ie=1,exp_de=1;
	Int4    i,j,x,arg,similarity=0,aa_per_io=20,aa_per_do=150,exp_ie=1,exp_de=1;
	Int4	time1=time(NULL);
 	ss_type Data=0,data=0;
	char	dms_mode=' ',str[100];
	Int4	max_in_seq=20000;
	Int4	Stages=4;	// works best for PH domain...
	// Int4	Stages=3;	// works for PH domain...??
	BooLean	RunGambit=FALSE,noseg=TRUE,fast=FALSE;
	double	prior_wt=0;
	cma_typ	final_cma,rcma,ocma,cma,xcma,bst_cma=0;
	UInt4	seed=18364592;
	Int4	ave_blk=10,avecol=7,minblk=3,maxblk=40,iter,maxiter=3;
	Int4	mhpsz=10,NumCandidates=5;
	double	dd,d,D,last_map=0,map,bst_map=0; 
	BooLean	output_rtf=FALSE,debug=FALSE;	// debug=TRUE;
// mhpsz=18; NumCandidates=9;
	FILE	*efp=0; // efp=stderr;

	//================ 0. get input parameters and file ===================
	if(argc < 2) print_error(GISMO_USAGE);
	FILE	*ofp,*rfp,*fptr = open_file(argv[1],".cmd","w");
	char	**ArgV; NEWP(ArgV,argc+20,char);
	for(i = 0; i < argc; i++) { fprintf(fptr,"%s ",argv[i]); ArgV[i]=argv[i]; } 
        for(arg = 2; arg < argc; arg++){
	   if(argv[arg][0] != '-') print_error(GISMO_USAGE);
           switch(argv[arg][1]) {
	      case 'C':
                if(argv[arg][2] != 0) similarity=IntOption(argv[arg],'C',10,100,GISMO_USAGE);
                else similarity=-1;  argv[arg][1]=' ';
		break;
	      case 'f':
		  if(strcmp("-fast",argv[arg]) == 0){ fast=TRUE; argv[arg][1]=' '; }
		  else { print_error(GISMO_USAGE); } break;
	      case 'H':
		  if(sscanf(argv[arg],"-H=%d",&NumCandidates) == 1){
			if(NumCandidates < 2 || NumCandidates > 9) print_error(GISMO_USAGE); 
			argv[arg][1]=' ';
		  } else print_error(GISMO_USAGE); 
		  break;
	      case 'h':
		  if(sscanf(argv[arg],"-h=%d",&mhpsz) == 1){
			if(mhpsz < 2 || mhpsz > 18) print_error(GISMO_USAGE); 
			argv[arg][1]=' ';
		  } else print_error(GISMO_USAGE); 
		  break;
	      case 'L': noseg=FALSE; argv[arg][1]=' '; break;
	      case 'm':
		  if(sscanf(argv[arg],"-maxseq=%d",&max_in_seq) == 1){
			if(max_in_seq < 20 || max_in_seq >= INT4_MAX) print_error(GISMO_USAGE); 
			argv[arg][1]=' ';
		  } else print_error(GISMO_USAGE); 
		  break;
#if 0
	      case 'p':
		  if(sscanf(argv[arg],"-p=%d:%d:%d:%d",&aa_per_io,&exp_ie,&aa_per_do,&exp_de) == 4){
			argv[arg][1]=' ';
		  } break;
#endif
	      case 'r':
		  {
		    if(strcmp("-rtf",argv[arg]) == 0){ output_rtf=TRUE; argv[arg][1]=' '; }
		    else { print_error(GISMO_USAGE); }
		  } break;
	      case 's':
		  if(sscanf(argv[arg],"-strategy=%d",&Stages) == 1){
			if(Stages < 2 || Stages > 4) print_error(GISMO_USAGE); 
		  } else if(sscanf(argv[arg],"-seed=%u",&seed) == 1){
			argv[arg][1]=' ';
		  } else print_error(GISMO_USAGE); break;
	      case 'w':
		  if(sscanf(argv[arg],"-wt=%lf",&prior_wt) == 1){
			if(prior_wt < 0.01 || prior_wt > 1000) print_error(GISMO_USAGE); 
			argv[arg][1]=' ';
		  } // else print_error(GISMO_USAGE); 
		  break;
	      default: break;	// ignore the rest...
	   }
	}
	if(mhpsz < NumCandidates) mhpsz=NumCandidates;
	if(fast) { mhpsz=NumCandidates=2; Stages=2; }
        Int4    number,*counts,J;
        unsigned short  *nsize;

	// Read in the sequences as an array to allow full and purged alignments.
        number = GetFastaInfo(argv[1], max_in_seq, &counts, &nsize, AB);
        e_type  *Seqs; NEW(Seqs,number+3,e_type);
        e_type  *tmpSeqs; NEW(tmpSeqs,number+3,e_type);
        FILE    *ifp= open_file(argv[1],"","r");
        for(J=1; J <= number; J++) tmpSeqs[J]=Seqs[J]=ReadSeq(ifp,J,nsize[J],AB);
        free(counts); free(nsize); fclose(ifp);
	data=Array2SeqSet(tmpSeqs,number,argv[1],AB);	// tmpSeqs is freed by Array2SeqSet().
	if(!noseg) data=PSegSeqSet(data); 

   	if(seed == 18364592){ seed = (UInt4) time(NULL); fprintf(fptr,"-seed=%d",seed);}
	if(prior_wt==0) prior_wt=1;
	fprintf(fptr,"\n"); fclose(fptr); sRandom(seed);
	if(debug){ sprintf(str,"%s_gmb",argv[1]); rfp=open_file(str,".log","w"); }
	else rfp=stderr;
	// PutLengthsSeqSet(stderr,data);
        Int4    M=ModeLengthSeqSet(data);
	fprintf(stderr,"mode=%d; ave = %d\n",M,AveSeqSeqSet(data));
        ave_blk = (Int4) ceil((double)M/10.0);        // A block of 5 residue every 10 residues.
	if(ave_blk > maxblk) ave_blk=maxblk;
	int	Argc=argc; ArgV[1]=0;
	Int4	timeGsm=time(NULL);
	ssx_typ	*issx=0,*jssx=0; dms_mode='F';
	gmb_typ	*gmb=0;
	cma_typ	*rtn_cma=0;
	gsm_typ *gsm=0;
	Int4	timeGp=0;

	// Purge the Data input set...
	//============= 1. Run gsm->Align() for one full set... ==================
	sprintf(str,"%s_prg",argv[1]); ArgV[1]=AllocString(str); RenameSeqSet(argv[1],data);
	i=Argc=argc; ArgV[i]=AllocString("-DoNotSeed"); Argc++;
	sprintf(str,"-h2"); i++; ArgV[i]=AllocString(str); Argc++;

	gsm = new gsm_typ(GISMO_USAGE,Argc,ArgV,AB,data); gsm->SetMaxIter(1);
        M=ModeLengthSeqSet(data); maxblk=MAX_NO_TYPESITES;
	rtn_cma=gsm->Align(ave_blk,avecol,minblk,maxblk,M);
	if(rtn_cma == 0) print_error("gismo++: failed to find a significant alignment");
	fprintf(rfp,"   GISMO PURGE: %d secs (%0.2f mins)\n",
			time(NULL)-timeGsm,(float)(time(NULL)-timeGsm)/60.0); fflush(rfp);

	//================  2. Purge the input set... ======================
	timeGp=time(NULL);
	for(i=1; rtn_cma[i]; i++) {
	   //================= 2. Output block-based alignment information ==================
	   fprintf(rfp,"%d: GISMO LLR = %0.2f; (%d blocks; %d columns).\n   ",
		i,UnGappedRelMapCMSA(rtn_cma[i]),nBlksCMSA(rtn_cma[i]),NumColumnsCMSA(rtn_cma[i]));
	   PutConfigCMSA(rfp,rtn_cma[i]); 
	   if(i==1){
	     cma=OneBlockCMSA(rtn_cma[1]); ExtendFakeToRealCMSA(cma);
	     // Do additional sampling to introduce gaps.
	     double	bc,ba,T,afq,rfq,prior_wt,MaxFrctn=0.05;
	     Int4		MinSticky=2; 
	     GetParameters(0,bc,ba,T,afq,rfq,aa_per_io,aa_per_do,exp_ie,exp_de,prior_wt,dms_mode);
	     gmb = new gmb_typ(aa_per_io,aa_per_do,exp_ie,exp_de,cma,dms_mode,prior_wt);
	     fprintf(rfp," "); gmb->SampleStickyTogether(rfp, MinSticky, MaxFrctn,0,200,'R');
	     delete gmb;
	   }
	} fprintf(stderr,"\tadd block time: %d seconds (%0.2f minutes)\n",
                        time(NULL)-timeGp,(float)(time(NULL)-timeGp)/60.0);
	
	//================  3. free up gismo ... ======================
	for(i=1; rtn_cma[i]; i++) NilCMSA(rtn_cma[i]); free(rtn_cma);
Int4 MinSetSize=gsm->RtnMinPurgedSetSize();
	delete gsm; gsm=0; 	// calls NilCMSA(rtn_cma[i]);  Does not free data...
	free(ArgV[1]); ArgV[1]=0;
	for(i=argc; i < Argc; i++){ free(ArgV[i]); ArgV[i]=0; }

#if 0	// run gambit to improve purging...This doesn't seem to help much...
        issx = new ssx_typ(aa_per_io,aa_per_do,exp_ie,exp_de,1000,cma,dms_mode);
	rcma=run_gambit(rfp,argv[1],0,issx,0,0); delete issx; 
	assert(rcma != 0); NilCMSA(cma); cma=rcma; rcma=0;
#endif
	Int4 Nset,Nsq=NumSeqsCMSA(cma),percent_ident=60;
#if 0	// old method...
	set_typ InSet=MakeSet(NumSeqsCMSA(cma)+4); FillSet(InSet);
        set_typ Set=RtnFastRepSetCMSA(efp, percent_ident,InSet,cma); Nset=CardSet(Set);
	fprintf(rfp,"   file (\"%s\"): (%d/%d removed; %d remain).\n",
                                                NameCMSA(cma),Nsq-Nset,Nsq,Nset);
#else	// new method
	set_typ Set=0,InSet=MakeSet(NumSeqsCMSA(cma)+4); FillSet(InSet);
	for(Nset=0,percent_ident=60; percent_ident <= 90 && Nset < MinSetSize; percent_ident+=5){
	   if(Set != 0) NilSet(Set);
           Set=RtnFastRepSetCMSA(efp, percent_ident,InSet,cma); Nset=CardSet(Set);
	   fprintf(rfp,"   file (\"%s\")(%d%c): (%d/%d removed; %d remain).\n",
                                     NameCMSA(cma),percent_ident,'%',Nsq-Nset,Nsq,Nset);
	}
#endif
        e_type  *RepSeqs; NEW(RepSeqs,number+3,e_type);
	for(i=0,J=1; J <= number; J++){
		if(MemberSet(J,Set)){ i++; RepSeqs[i]=Seqs[J]; }
	} assert(i==Nset);
	Data=data;  data=Array2SeqSet(RepSeqs,Nset,argv[1],AB); // RepSeqs is freed by Array2SeqSet()
	NilCMSA(cma); 	// don't need this anymore...
	NilSet(InSet);

	timeGsm=time(NULL);
	if(ArgV[1]) free(ArgV[1]); 
	sprintf(str,"%s_gsm",argv[1]); ArgV[1]=AllocString(str); RenameSeqSet(str,data);
	i=Argc=argc; ArgV[i]=AllocString("-DoNotSeed"); Argc++;
	sprintf(str,"-h%d",mhpsz); i++; ArgV[i]=AllocString(str); Argc++;

	//================ 1. Run GISMO sampler (longest step) ===================
	gsm = new gsm_typ(GISMO_USAGE,Argc,ArgV,AB,data);	// MaxIter set to default (5?).
        M=ModeLengthSeqSet(data); maxblk=MAX_NO_TYPESITES;
	rtn_cma=gsm->Align(ave_blk,avecol,minblk,maxblk,M);
	if(rtn_cma == 0) print_error("gismo++: failed to find a significant alignment");
	if(gsm){ delete gsm; gsm=0; }
	if(ArgV[1]) free(ArgV[1]); 
	for(i=argc; i < Argc; i++) if(ArgV[i]) free(ArgV[i]); free(ArgV);
	fprintf(rfp,"   GISMO SEEDS: %d secs (%0.2f mins)\n",
			time(NULL)-timeGsm,(float)(time(NULL)-timeGsm)/60.0); fflush(rfp);
	// pre-sampling step...
	assert(mhpsz <= 18);
	Int4	num_shared=0,bst_core;
	timeGp=time(NULL);
#if 0
	for(i=1; rtn_cma[i]; i++) {
	   //================= 2. Output block-based alignment information ==================
	   fprintf(rfp,"%d: GISMO LLR = %0.2f; (%d blocks; %d columns).\n   ",
		i,UnGappedRelMapCMSA(rtn_cma[i]),nBlksCMSA(rtn_cma[i]),NumColumnsCMSA(rtn_cma[i]));
	   PutConfigCMSA(rfp,rtn_cma[i]); 
	   num_shared++; 
	}
#endif
	cma_typ	*start_cma=RtnBestCore(rfp,rtn_cma,NumCandidates); 	// NilCMSA(start_cma);
	for(num_shared=0,i=1; rtn_cma[i]; i++) num_shared++;
	// for(i=1; rtn_cma[i]; i++) NilCMSA(rtn_cma[i]); free(rtn_cma); rtn_cma=start_cma;
	fprintf(rfp,"\tblock-based time: %d seconds (%0.2f minutes)\n\n",
                        time(NULL)-timeGp,(float)(time(NULL)-timeGp)/60.0);

	//=========== Find shared core for the gismo block-based alignments. ============
	double	shared[20][20],CL,scr,bst_scr;
	cma_typ	*one_cma; NEW(one_cma,num_shared + 3, cma_typ);
	Int4	*self_cols;  NEW(self_cols,num_shared+3,Int4);
	double	*self_core;  NEW(self_core,num_shared+3,double);
	for(i=1; rtn_cma[i]; i++) {	// convert to one block.
	  fprintf(rfp,"%d: GISMO LLR = %0.2f; (%d blocks; %d columns).\n",
		i,UnGappedRelMapCMSA(rtn_cma[i]),nBlksCMSA(rtn_cma[i]),NumColumnsCMSA(rtn_cma[i]));
	  fflush(rfp); 
	  rcma=start_cma[i];
	  // one_cma[i]=rcma; continue;
	  aa_per_io=20; aa_per_do=100; prior_wt=0.50; // for stages 0-1.
          issx = new ssx_typ(aa_per_io,aa_per_do,exp_ie,exp_de,1000,rcma,dms_mode);
	  issx->SetPriorWt(prior_wt);
	  one_cma[i]=run_gambit(rfp,argv[1],i,issx,0,1); delete issx; 

	  //============= Compute the core that the original and optimized alignments share.
	  double  frq_cut=0.50,ave_frq=0;
          cma_typ *tcma=GetCoreCMSA(one_cma[i],rcma,ave_frq,frq_cut);
	  if(tcma != 0){
	    j=LengthCMSA(1,tcma[1]); dd=(double)j * ave_frq;
	    // fprintf(rfp,"  Shared core: %d columns; score = %.3f\n",j,dd); 
	    self_cols[i]=j; self_core[i]=dd;
	    NilCMSA(tcma[1]); NilCMSA(tcma[2]); free(tcma);
	  } else { self_cols[i]=0; self_core[i]=0.0; }
	  NilCMSA(rcma);
	  shared[i][i]=0.0;
	} free(start_cma);

	//================= Find the shared cores between alignments. ===============
	cma_typ	strt_cma=0;
	for(i=1; i < num_shared; i++) {
	  assert(rtn_cma[i]); assert(one_cma[i]);
	  for(j=i+1; j <=num_shared; j++) {
	     assert(rtn_cma[j]); assert(one_cma[j]);
	     Int4    lenI=LengthCMSA(1,one_cma[i]);
	     double  *DD,*Del,frq_cut=0.50,ave_frq=0,sum;
	     Int4    n,*ColI2J=CommonColsCMSA(one_cma[i],one_cma[j],DD,Del,frq_cut);
	     for(sum=0.0,n=0,x=1; x <= lenI; x++){
		if(ColI2J[x] > 0){ n++; sum+=DD[x]; }
	     } free(DD); free(Del); free(ColI2J);
	     if(n > 0) ave_frq=sum/(double)n; else ave_frq=0.0;
	     CL=(double) n*ave_frq;
	     shared[i][i] += CL; shared[j][j] += CL;
	     shared[i][j]=shared[j][i]=CL;
	  }
	}
	//================= Pick the best starting (core) alignment. ===================
	//================= Pick the best starting (core) alignment. ===================
	//================= Pick the best starting (core) alignment. ===================
	// pick core alignment with best overall score with other alignments.
	for(bst_core=0,bst_scr=0,i=1; i <= num_shared; i++){
	   // CL=shared[i][i];
	   // CL=shared[i][i] + self_core[i];
	   CL=shared[i][i] + NumColumnsCMSA(one_cma[i]);
	   if(CL >= bst_scr){ bst_scr=CL; bst_core=i; }
	}
#if 1	// find all columns shared by at least one other alignment.
	for(bst_scr=0,i=bst_core,j=1; j <= num_shared; j++){
	   if(i != j && (CL=shared[i][j]) >= bst_scr){ bst_scr=CL; }
	}
     {
	if(strt_cma) NilCMSA(strt_cma); 
	i=bst_core; bst_cma=one_cma[i]; assert(one_cma[i]); 
	Int4    lenI=LengthCMSA(1,one_cma[i]);
	set_typ SetI=MakeSet(lenI+4); ClearSet(SetI);
	for(j=1; j <= num_shared; j++){
	     if(i==j) continue;
	     double  *DD,*Del,frq_cut=0.50,ave_frq,sum;
	     Int4    n,*ColI2J=CommonColsCMSA(one_cma[i],one_cma[j],DD,Del,frq_cut);
	     for(sum=0.0,n=0,x=1; x <= lenI; x++){
		if(ColI2J[x] > 0){ n++; sum+=DD[x]; AddSet(x,SetI); }
	     } free(DD); free(Del); free(ColI2J);
	     if(n > 0) ave_frq=sum/(double)n; else ave_frq=0.0;
	}
	char    *operA; NEW(operA,lenI+5,char); operA[0]='E'; operA[lenI+1]='E';
	for(x=1; x <= lenI; x++){ if(!MemberSet(x,SetI)) operA[x]='d'; else operA[x]='m'; }
	strt_cma=RemoveTheseColumnsCMSA(operA,bst_cma); free(operA);
	if(strt_cma == 0){ strt_cma=CopyCMSA(one_cma[bst_core]); }
	NilSet(SetI);
     }
	fprintf(rfp,"\n   trimmed alignment down to %d shared columns\n",NumColumnsCMSA(strt_cma));
#else	// use the entire starting alignment instead of the core...
	strt_cma=CopyCMSA(one_cma[bst_core]); // one_cma[bst_core]=0;
#endif
	//=============== print out starting core information... ===============
	fprintf(rfp,"\n best core = %d (%.3f score; %.2f aggregate; %d cols).\n   :",
		bst_core,bst_scr,shared[bst_core][bst_core],NumColumnsCMSA(strt_cma));
	for(j=1; j <= num_shared; j++) fprintf(rfp," %8d",j); fprintf(rfp,"\n");
	for(i=1; i <= num_shared; i++) {
		fprintf(rfp,"%3d:",i); assert(one_cma[i]);
	        for(j=1; j <= num_shared; j++){ fprintf(rfp," %8.2f",shared[i][j]); }
		fprintf(rfp," (%d/%d cols; %.3f (s)core; %d blks --> %d cols)\n",
			self_cols[i],NumColumnsCMSA(rtn_cma[i]),self_core[i],
			nBlksCMSA(rtn_cma[i]),NumColumnsCMSA(one_cma[i]));
		NilCMSA(one_cma[i]); NilCMSA(rtn_cma[i]);
	} fprintf(rfp,"\n"); fflush(rfp); free(rtn_cma); free(one_cma);
	free(self_cols); free(self_core);
	//========================================================================
	//========================================================================
// dsw=rssx->RtnWtNumSeqs();

	//========================================================================
	//================== The gambit part of the program ======================
	//========================================================================
	//================= Make a new full size cma from strt_cma. ==============
	aa_per_io=20; aa_per_do=150; prior_wt=1.0;	// stage 1.5.
	double T,bc,ba,afq,rfq;	// these are not used here...
GetParameters(1,bc,ba,T,afq,rfq,aa_per_io,aa_per_do,exp_ie,exp_de,prior_wt,dms_mode);
	gmb = new gmb_typ(aa_per_io,aa_per_do,exp_ie,exp_de,strt_cma,dms_mode,prior_wt);
        // if(AdjSqWt > 0.0) gmb->SetSqWtAdjust(AdjSqWt);
	rcma=gmb->CreateFullCMSA(Data,Set);
        /// ofp = open_file(argv[1],"_finalB.cma","w"); PutCMSA(ofp,rcma);  fclose(ofp);
        delete gmb;
        NilCMSA(strt_cma);    // does not delete data.
	e_type *xSeqs=NilSeqSetRtnSeqs(data); free(xSeqs); // free array only!!! Seqs in use by Data...
	data=Data; strt_cma=rcma;  // Data is in here...

#if 1	// Do additional sampling only on sequences not in Set.
	GetParameters(1,bc,ba,T,afq,rfq,aa_per_io,aa_per_do,exp_ie,exp_de,prior_wt,dms_mode);
	gmb = new gmb_typ(aa_per_io,aa_per_do,exp_ie,exp_de,strt_cma,dms_mode,prior_wt);
	gmb->SetExcludeSet(Set);
	for(i=1; i <= 2; i++){
		Int4	MinSticky=2; 
		double	MaxFrctn=0.05;
		gmb->SampleStickyTogether(rfp, MinSticky, MaxFrctn,0,100,'R');
		// double  lpr=gmb->RandomTogether(MinSticky, MaxFrctn, temp);
	} gmb->UnSetExcludeSet( ); NilSet(Set);
        delete gmb;
	GetParameters(4,bc,ba,T,afq,rfq,aa_per_io,aa_per_do,exp_ie,exp_de,prior_wt,dms_mode);
#endif

	cma_typ	*saved_cma; NEW(saved_cma,num_shared+5,cma_typ);
#if 1
	//================= 3. Run GAMBIT on MSAs returned by GISMO ==================
        issx = new ssx_typ(aa_per_io,aa_per_do,exp_ie,exp_de,1000,strt_cma,dms_mode);
	issx->SetPriorWt(prior_wt);
	saved_cma[1]=run_gambit(rfp,argv[1],bst_core,issx,2,Stages); // run stages 2-4.
	sprintf(str,"%s_gmb%d",argv[1],bst_core); RenameCMSA(str,saved_cma[1]);
	// ofp=open_file(str,".cma","w"); PutCMSA(ofp,saved_cma[1]); fclose(ofp);
	delete issx; NilCMSA(strt_cma);
#elif 1
        issx = new ssx_typ(aa_per_io,aa_per_do,exp_ie,exp_de,1000,strt_cma,dms_mode);
	saved_cma[1]=run_gambit(rfp,argv[1],bst_core,issx,2,2); // run stage 2 only.
	sprintf(str,"%s_gmb%d",argv[1],bst_core); RenameCMSA(str,saved_cma[1]);
	// ofp=open_file(str,".cma","w"); PutCMSA(ofp,saved_cma[1]); fclose(ofp);
	delete issx; NilCMSA(strt_cma); strt_cma=saved_cma[1]; saved_cma[1]=0;

	//================= 3. Run stage 3 GAMBIT ==================
	aa_per_io=30; aa_per_do=170; // for stage 3;
        issx = new ssx_typ(aa_per_io,aa_per_do,exp_ie,exp_de,1000,strt_cma,dms_mode);
	saved_cma[1]=run_gambit(rfp,argv[1],bst_core,issx,3,3); // starts & ends at stage 3.
	delete issx; NilCMSA(strt_cma); 
	strt_cma=saved_cma[1];  saved_cma[1]=0;

	//================= 3. Run GAMBIT on final alignment set ==================
	aa_per_io=30; aa_per_do=200; // for stage 4;
        issx = new ssx_typ(aa_per_io,aa_per_do,exp_ie,exp_de,1000,strt_cma,dms_mode);
	// saved_cma[1]=run_gambit(rfp,argv[1],bst_core,issx,4,Stages); // run stages 3-4.
	saved_cma[1]=run_gambit(rfp,argv[1],bst_core,issx,4,4); // run stage 4.
	sprintf(str,"%s_gmb%d",argv[1],bst_core); RenameCMSA(str,saved_cma[1]);
	ofp=open_file(str,".cma","w"); PutCMSA(ofp,saved_cma[1]); fclose(ofp);
	delete issx; NilCMSA(strt_cma);
#else	
	//================= 3. Run GAMBIT on MSAs returned by GISMO ==================
        issx = new ssx_typ(aa_per_io,aa_per_do,exp_ie,exp_de,1000,strt_cma,dms_mode);
	saved_cma[1]=run_gambit(rfp,argv[1],bst_core,issx,2,2); // run stage 2 only.
	sprintf(str,"%s_gmb%d",argv[1],bst_core); RenameCMSA(str,saved_cma[1]);
	// ofp=open_file(str,".cma","w"); PutCMSA(ofp,saved_cma[1]); fclose(ofp);
	delete issx; NilCMSA(strt_cma); strt_cma=saved_cma[1]; saved_cma[1]=0;

	//=============== Purge the new set for stage 2... ======================
	cma=strt_cma;
	Int4	timeP=time(NULL); percent_ident=66;
	gmb = new gmb_typ(aa_per_io,aa_per_do,exp_ie,exp_de,cma,dms_mode,prior_wt);
	strt_cma=gmb->CreatePurgedCMSA(percent_ident,Seqs,Set,data);
	delete gmb; NilCMSA(cma);
	fprintf(rfp,"\t2nd purge time: %d seconds (%0.2f minutes)\n",
                        time(NULL)-timeP,(float)(time(NULL)-timeP)/60.0);

	//================= 3. Run GAMBIT on new purged alignment ==================
        issx = new ssx_typ(aa_per_io,aa_per_do,exp_ie,exp_de,1000,strt_cma,dms_mode);
	saved_cma[1]=run_gambit(rfp,argv[1],bst_core,issx,3,3); // starts & ends at stage 3.
	delete issx; NilCMSA(strt_cma); 
	strt_cma=saved_cma[1];  saved_cma[1]=0;

	//================= Make a new full size cma from final_cma.
	gmb = new gmb_typ(aa_per_io,aa_per_do,exp_ie,exp_de,strt_cma,dms_mode,prior_wt);
	rcma=gmb->CreateFullCMSA(Data,Set); delete gmb; NilCMSA(strt_cma);    // does not delete data.
	xSeqs=NilSeqSetRtnSeqs(data); free(xSeqs); // free array only!!! Sequences in use by Data...
	data=Data; strt_cma=rcma;  // Data is in here...
	NilSet(Set);
	
	//================= 3. Run GAMBIT on final full alignment set ==================
        issx = new ssx_typ(aa_per_io,aa_per_do,exp_ie,exp_de,1000,strt_cma,dms_mode);
	// saved_cma[1]=run_gambit(rfp,argv[1],bst_core,issx,4,Stages); // run stages 3-4.
	saved_cma[1]=run_gambit(rfp,argv[1],bst_core,issx,4,4); // run stage 4.
	sprintf(str,"%s_gmb%d",argv[1],bst_core); RenameCMSA(str,saved_cma[1]);
	ofp=open_file(str,".cma","w"); PutCMSA(ofp,saved_cma[1]); fclose(ofp);
	delete issx; NilCMSA(strt_cma);
#endif
	//************** compare the files to find the best alignment **********************
	//************** compare the files to find the best alignment **********************
	//************** compare the files to find the best alignment **********************
	{
	   ssx_typ *rssx,**ssx; NEWP(ssx,mhpsz+5,ssx_typ);
	   double  bst_lpr=-9999999999.9,pn=1000; 
	   Int4	   MinCol=INT4_MAX,len;
	   // dms_mode='M';	// 110 components
	   // dms_mode='L';	// 134 components  doesn't seem to help...
	   //-------------- Find minimum core length after trimming columns. ---------------
	   for(iter=1; saved_cma[iter]; iter++){
#if 0
		cma=saved_cma[iter]; assert(pn == PerNatsCMSA(cma)); i=NumColumnsCMSA(cma);
	   	Int4 nColsRm=RmGappyColumnsCMSA(0.20,cma); saved_cma[iter]=cma; j=NumColumnsCMSA(cma);
	   	fprintf(rfp,"  %d: %d/%d columns removed (%d cols remain)\n",iter,nColsRm,i,j);
		if((len=LengthCMSA(1,cma)) < MinCol) MinCol=len;
#else	// fix problem with overly aggressive 
		cma=saved_cma[iter]; assert(pn == PerNatsCMSA(cma)); i=NumColumnsCMSA(cma);
		for(double gap_cut=0.20; gap_cut <= 0.50; ){
		  cma_typ zcma=CopyCMSA(cma); i=NumColumnsCMSA(zcma);
	   	  Int4 nColsRm=RmGappyColumnsCMSA(gap_cut,zcma); j=NumColumnsCMSA(zcma);
	   	  fprintf(rfp,"  %d: %d/%d columns removed (%d cols remain)\n",iter,nColsRm,i,j);
		  gap_cut+=0.10;
		  if(j < 35 && gap_cut <= 0.50) NilCMSA(zcma);
		  else {
		     NilCMSA(cma); cma=zcma; saved_cma[iter]=zcma; 
		     if((len=LengthCMSA(1,cma)) < MinCol) MinCol=len; break;
		  } 
		} 
#endif
	   } //--------------- 
	   double  nWtSq,MaxWtSq=0; 
	   //========== New heuristics for finding the optimum structural alignment. =========
	   cma_typ *core_cma; NEW(core_cma,iter+5,cma_typ);
#if 0	// Retain only 90% of the highest-scoring core columns.
	   MinCol = (Int4) ceil(0.90*(double)MinCol -0.49);	
	   MinCol=MAXIMUM(Int4,MinCol,10);  
#endif
	   //----------------- shrink msa down to the common core. ------------------
	   for(iter=1; saved_cma[iter]; iter++){
		xcma=saved_cma[iter]; assert(pn == PerNatsCMSA(xcma));
                ssx_typ *tssx = new ssx_typ(aa_per_io,aa_per_do,exp_ie,exp_de,pn,xcma,dms_mode);
                rcma=RmWorstColsCMSA(MinCol,tssx);	// create alignments with same # columns
		if(rcma) core_cma[iter]=rcma; else core_cma[iter]=CopyCMSA(xcma);
		delete tssx;
	   }
	   //--------------- Get the sequence weights; save the maximum weight. ----------------
	   for(iter=1; saved_cma[iter]; iter++){
		cma=core_cma[iter]; assert(pn == PerNatsCMSA(cma));
                ssx[iter] = new ssx_typ(aa_per_io,aa_per_do,exp_ie,exp_de,pn,cma,dms_mode);
		if((nWtSq=ssx[iter]->RtnWtNumSeqs()) > MaxWtSq) MaxWtSq=nWtSq;
	   } bst_cma=0; 
	   //--------------- compare the adjusted BILD scores for the core alignments. ----------------
	   for(iter=1; core_cma[iter]; iter++){
		rcma=core_cma[iter]; // ExtendFakeToRealCMSA(rcma);
		rssx = ssx[iter];
                double ddd,dd,dsw,gp,DD;
#if 1
		ddd=rssx->GapMap();
                // dd=rssx->AdjstdBildLLR(MaxWtSq);	// Adjust for different # weighted seqs.
                dd=rssx->AdjstDirichletLLR(MaxWtSq);	// Adjust for different # weighted seqs.
                dsw=rssx->RtnWtNumSeqs(); gp=rssx->RtnIndelPenalty();
                fprintf(rfp,"LLR=%.2f; AdjstLLR(%d cols)=%.2f; SqWt=%.1f; penalty=%.2f; NetLLR=%.2f; dmp='%c'.\n",
                        ddd,MinCol,dd,dsw,gp,dd+gp,dms_mode); fflush(rfp);
		// dd += gp;
#else
		dd=rssx->GapMap();
#endif
                delete rssx; ssx[iter]=0;   // delete rssx first and then rcma to avoid a core dump!
		if(dd > bst_lpr){
		    if(bst_cma) NilCMSA(bst_cma); 
		    bst_cma=saved_cma[iter]; bst_lpr=dd;
		} else NilCMSA(saved_cma[iter]); saved_cma[iter]=0;
#if 0
	   	sprintf(str,"%s_gmbY%d",argv[1],iter); RenameCMSA(str,rcma);
	   	ofp=open_file(str,".cma","w"); PutCMSA(ofp,rcma); fclose(ofp);
#elif 0
	   	sprintf(str,"%s_gmb%d",argv[1],iter); RenameCMSA(str,rcma);
	   	ofp=open_file(str,".cma","w"); PutCMSA(ofp,rcma); fclose(ofp);
#endif
		NilCMSA(rcma); core_cma[iter]=0;
	   } free(saved_cma); free(core_cma); free(ssx);
	}
	// This bst_cma file is getting corrupted at some point after 
	sprintf(str,"%s_gmb",argv[1]); 
	// ofp=open_file(str,".cma","w"); PutCMSA(ofp,bst_cma); fclose(ofp);
	//************** end: compare the files to find the best alignment **********************
	{	// sample the best at the end.
	   sprintf(str,"%s_bst",argv[1]); 
	   char *name=AllocString(str);
#if 0	// Shutting off for now; doesn't improve much, but might be able to improve...
	   cma=SampleSubAlignsCMSA(str,bst_cma,dms_mode);
	   sprintf(str,"%s_final.cma",name); WriteCMSA(str,cma);
	   gmb=new gmb_typ(aa_per_io,aa_per_do,exp_ie,exp_de,cma,dms_mode,prior_wt);
           //============ sample sequences in final alignment one-at-a-time ================
           gmb->Sample(str,'S',similarity,150.0,0.0,50,rfp); delete gmb;
           sprintf(str,"%s_finalA.cma",name); WriteCMSA(str,cma);
	   Int4 nColsRm=RmGappyColumnsCMSA(0.20,cma);
	   fprintf(stderr,"\n %d columns removed\n",nColsRm);
           sprintf(str,"%s_finalB.cma",name); WriteCMSA(str,cma); free(name);
	   final_cma=cma; 
#else
#if 0	// did this earlier...
	   Int4 nColsRm=RmGappyColumnsCMSA(0.20,bst_cma);
	   fprintf(stderr,"\n %d columns removed\n",nColsRm);
           sprintf(str,"%s_finalB.cma",name); WriteCMSA(str,bst_cma); free(name);
#endif
	   final_cma=bst_cma; 
#endif
	   free(name);
	} if(rfp != stderr) fclose(rfp); free(Seqs);
#if 0	//================= Make a new full size cma from final_cma.
	gmb = new gmb_typ(aa_per_io,aa_per_do,exp_ie,exp_de,final_cma,dms_mode,prior_wt);
        // if(AdjSqWt > 0.0) gmb->SetSqWtAdjust(AdjSqWt);
        // std::cerr << "DEBUG X" << std::endl;

        rcma=gmb->CreateAlign(Data);	// samples at Temperature = 0.
        ofp = open_file(argv[1],"_finalB.cma","w"); PutCMSA(ofp,rcma);  fclose(ofp);
        delete gmb;
        TotalNilCMSA(final_cma);  // deletes data as well.
	final_cma=rcma;		  // Data is in here...
#endif
	fprintf(stderr,"\ttime: %d seconds (%0.2f minutes)\n",
                        time(NULL)-time1,(float)(time(NULL)-time1)/60.0);
#if 1
	if(!debug){
	  cma_typ cma0 = AddConsensusCMSA(final_cma);
	  ofp = open_file(argv[1],".fa","w"); PutFastaCMSA(ofp,cma0); fclose(ofp);
	} else { sprintf(str,"%s_bst_finalA.cma",argv[1]); WriteCMSA(str,bst_cma); } 
	sprintf(str,"%s.cma",argv[1]); WriteCMSA(str,bst_cma);
#endif
	if(output_rtf) { return final_cma; } // final_cma and AB and data freed in calling environment...
	else { TotalNilCMSA(final_cma); return 0; }
}

