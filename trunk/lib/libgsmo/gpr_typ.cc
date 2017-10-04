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
#include "gpr_typ.h"

Int4	gpr_typ::DoPurgeSampling(FILE *rfp, cma_typ &cma)
{
	double	d,dd,T,pfrq,WtSq,lpr0,lpr,bst_lpr;
	cma_typ	bcma=0,rcma=0;
	Int4	rtn=0;

	timeR=time(NULL);  
	if(gmb) delete gmb; gmb=0; 
        gmb=new gmb_typ(aa_per_io,aa_per_do,exp_ie,exp_de,cma,dms_mode,prior_wt);
	gmb->SetPriorWt(prior_wt); gmb->SetMaxIter(3);

	if(stage == 0){
#if 1
	  char	Strategy[20]=" RX   ";
                                    //   R    X            
	  double	sfrq[20]={ 0.0,0.15,0.20,0.0,0.0,0.0,0.0};
	  lpr0=gmb->RtnMap();
	  rcma=gmb->SampleStickyTogether(rfp,MinSticky,Strategy,sfrq,Temp,StMd);
	  lpr=gmb->RtnMap(); assert(rcma == cma); 	// Don't call P in Strategy...
	  if(lpr > lpr0) { rtn++; routine='s'; if(write) this->DoWrite(cma); }
#else
	  this->DoStickySample(rfp,cma);
#endif
	} if(1) {
	  Int4	dt,i,stop=6,percent_id=30,timeI; // bcma=CopyCMSA(cma); 
#if 0
	  if(stage == 1){ T=200; pfrq=0.95; stop=2; }
	  else if(stage == 2){ T=100;  pfrq=0.95; stop=3; }
	  else if(stage == 3){ T=50;  pfrq=0.90; stop=3; }
	  else { T=0; pfrq=0.85; stop=2; }
	  for(d=1.0,i=1; (d > 0 && i <= stop); pfrq-= 0.05,T=0,i++){
	    timeI=time(NULL); // lpr0=gmb->RtnMap();
	    if(!gmb->SamplePurged(percent_id,pfrq,T,lpr0,lpr,stage)) continue; 
	    cma=gmb->RtnCMSA();	// WARNING: cma has changed...
	    fprintf(rfp,"   SamplePurged (%.2f)(%d cols): LLR = %.2f (%.2f);",
			pfrq,LengthCMSA(1,cma),lpr,lpr-lpr0);
	    d=lpr-lpr0; dt=time(NULL)-timeI; WtSq=gmb->RtnWtNumSeqs(); 
            if(dt < 60) fprintf(rfp," %.1f wt.sq. (%d secs)(%.1f K)\n", WtSq,dt,T);
            else fprintf(rfp," %.2f WtSq (%0.2f mins)(%.1f K)\n", WtSq,(float)(dt)/60.0,T);
	    if(lpr > lpr0){
		rtn++; routine='p'; if(write) this->DoWrite(cma); // NilCMSA(bcma); bcma=CopyCMSA(cma); 
	    } else break; 
	  }
#elif 1
	  if(stage == 0){ T=150; stop=1; }
	  else if(stage == 1){ T=75; stop=2; }
	  else if(stage == 2){ T=0; stop=1; }
	  else if(stage == 3){ T=0; stop=1; }
	  else { T=0; stop=1; }
	  Int4 minsize=50;
	  for(d=1.0,i=1; (d > 0 && i <= stop); T=0,i++){
	    timeI=time(NULL); // lpr0=gmb->RtnMap();
	    if(!gmb->SampleByLayers(minsize,T,lpr0,lpr,stage)) continue;
	    cma=gmb->RtnCMSA();	// WARNING: cma has changed...
	    fprintf(rfp,"   SamplePurged (%d)(%d cols): LLR = %.2f (%.2f);",
			minsize,LengthCMSA(1,cma),lpr,lpr-lpr0);
	    d=lpr-lpr0; dt=time(NULL)-timeI; WtSq=gmb->RtnWtNumSeqs(); 
            if(dt < 60) fprintf(rfp," %.1f wt.sq. (%d secs)(%.1f K)\n", WtSq,dt,T);
            else fprintf(rfp," %.2f WtSq (%0.2f mins)(%.1f K)\n", WtSq,(float)(dt)/60.0,T);
	    if(lpr > lpr0){
		rtn++; routine='p'; if(write) this->DoWrite(cma); // NilCMSA(bcma); bcma=CopyCMSA(cma); 
	    } else break; 
	  }
#else
	  if(stage == 0){ T=200; } else if(stage == 1){ T=150; } else if(stage == 2){ T=75; }
	  else if(stage == 3){ T=0; } else { T=0; }
	  timeI=time(NULL); // lpr0=gmb->RtnMap();
	  if(gmb->SampleByLayers(50,T,lpr0,lpr,stage)){
	    cma=gmb->RtnCMSA();	// WARNING: cma has changed...
	    fprintf(rfp,"   SamplePurged (%.2f)(%d cols): LLR = %.2f (%.2f);",
			pfrq,LengthCMSA(1,cma),lpr,lpr-lpr0);
	    d=lpr-lpr0; dt=time(NULL)-timeI; WtSq=gmb->RtnWtNumSeqs(); 
            if(dt < 60) fprintf(rfp," %.1f wt.sq. (%d secs)(%.1f K)\n", WtSq,dt,T);
            else fprintf(rfp," %.2f WtSq (%0.2f mins)(%.1f K)\n", WtSq,(float)(dt)/60.0,T);
	    if(lpr > lpr0){
		rtn++; routine='p'; if(write) this->DoWrite(cma); // NilCMSA(bcma); bcma=CopyCMSA(cma); 
	    } 
	  }
#endif
	} delete gmb; gmb=0; // NilCMSA(cma); cma=bcma; // revert back to best cma (bcma).
	fprintf(rfp,"   time (DoPurgeSampling): %d secs",time(NULL)-timeR); 
	fprintf(rfp," (%0.2f mins)\n",(float)(time(NULL)-timeR)/60.0); fflush(rfp);
	return rtn;
}

cma_typ	gpr_typ::DoStickySample(FILE *rfp, cma_typ &cma)
// Do initial broad sampling to get unstuck...
{
	// this->DoAddRmColumns(rfp,cma,2);	// do for 2 cycles.
	char	Strategy[20]=" RTNPXDVSW    ";  // 'D' ruins the alignment?
	// double	s1frq[20]={ 0.0,0.15,0.10,0.10,0.49,0.15,0.15,0.05,0.15,0.08,0.0,0,0,0.0};
                             //   R    T    N    P    X    D    V    S    W      
	double	s0frq[20]={ 0.0,0.15,0.10,0.10,0.45,0.15,0.08,0.05,0.10,0.15,0.0,0,0,0.0};
	double	s1frq[20]={ 0.0,0.15,0.10,0.10,0.30,0.15,0.08,0.05,0.08,0.15,0.0,0,0,0.0};
	double	s2frq[20]={ 0.0,0.10,0.10,0.10,0.20,0.10,0.08,0.05,0.07,0.15,0.0,0,0,0.0};
	double	s3frq[20]={ 0.0,0.05,0.10,0.10,0.15,0.10,0.08,0.05,0.06,0.15,0.0,0,0,0.0};
	double	s4frq[20]={ 0.0,0.01,0.10,0.10,0.10,0.10,0.08,0.05,0.05,0.05,0.0,0,0,0.0};
	double	d,dd,lpr,best,*sfrq;
	if(stage == 0) sfrq=s0frq; 
	else if(stage==1) sfrq=s1frq; 
	else if(stage==2) sfrq=s2frq;
	else if(stage==3) sfrq=s3frq;
	else sfrq=s4frq;
	timeR=time(NULL);  Cycle=1;
	if(gmb) delete gmb; gmb=0; 
        gmb=new gmb_typ(aa_per_io,aa_per_do,exp_ie,exp_de,cma,dms_mode,prior_wt);
	gmb->SetPriorWt(prior_wt); 
	gmb->SetMaxIter(3);
	cma_typ rcma=gmb->SampleStickyTogether(rfp,MinSticky,Strategy,sfrq,Temp,StMd);
	delete gmb; gmb=0; cma=rcma;
	routine='s'; if(write) this->DoWrite(rcma);
	fprintf(rfp,"   StickySample: %d secs",time(NULL)-timeR); 
	fprintf(rfp," (%0.2f mins)\n",(float)(time(NULL)-timeR)/60.0); fflush(rfp);
	return rcma;
}

cma_typ	gpr_typ::DoSingleSample(FILE *rfp, cma_typ &cma)
// Do initial broad sampling to get unstuck...
{
	// this->DoAddRmColumns(rfp,cma,2);	// do for 2 cycles.
	// char	Strategy[20]=" ROROR    ";  
	char	Strategy[20]=" RO    ";  
                            //   R    O    R    O    R
	double	sfrq[20]={ 0.0,0.01,0.01,0.01,0.01,0.01,0.0,0,0,0.0};
	double	d,dd,lpr,best;
	timeR=time(NULL);  Cycle=1;
	if(gmb) delete gmb; gmb=0; 
        gmb=new gmb_typ(aa_per_io,aa_per_do,exp_ie,exp_de,cma,dms_mode,prior_wt);
	gmb->SetPriorWt(prior_wt); 
	gmb->SetMaxIter(1);
	cma_typ rcma=gmb->SampleStickyTogether(rfp,MinSticky,Strategy,sfrq,Temp,StMd);
	delete gmb; gmb=0; cma=rcma;
	routine='s'; if(write) this->DoWrite(rcma);
	fprintf(rfp,"   SingleSample: %d secs",time(NULL)-timeR); 
	fprintf(rfp," (%0.2f mins)\n",(float)(time(NULL)-timeR)/60.0); fflush(rfp);
	return rcma;
}

BooLean	gpr_typ::DoMvColumns(FILE *rfp, cma_typ &cma, Int4 printCycle)
{
	BooLean	did_mv=FALSE;
	timeR=time(NULL); Cycle=0; sprintf(str,"%s_gmb_tmp",name);
	cma_typ	rcma=0,bcma=0;
	char	Strategy[20]=" RTNXVS  ";  // 'D' ruins the alignment!!!
                                           //  R    T    N    X    V    S  
	double	d,dd,lpr,best,sfrq[20]={ 0.0,0.12,0.10,0.10,0.10,0.10,0.10,0.0,0.0};
        bcma=cma; cma=0; SaveBestCMSA(bcma); 	// saves the best column config as well.
	gmb_typ	*gmb0=new gmb_typ(aa_per_io,aa_per_do,exp_ie,exp_de,bcma,dms_mode,prior_wt);
        gmb0->SetPriorWt(prior_wt);

	BooLean	improved=FALSE;
	assert(nBlksCMSA(bcma) == 1);
	set_typ SetSq=0; best=gmb0->RtnMap(); 
        do {
          rcma=MvColumnsCMSA(gmb0->RtnSSX(),SetSq);
          if(rcma) {
		Cycle++; assert(rcma != cma); 
	        if(gmb) delete gmb; gmb=0; 
		if(cma) NilCMSA(cma); cma=rcma; rcma=0; 
                gmb=new gmb_typ(aa_per_io,aa_per_do,exp_ie,exp_de,cma,dms_mode,prior_wt,gmb0->RtnSSX());
		gmb->SetMaxIter(3);
                       // sprintf(str,"%s_gmb%d",argv[1],i);  RenameCMSA(str,cma); 
		double dd=(double)CardSet(SetSq)/(double)NumSeqsCMSA(cma);
		d=gmb->RtnMap();	// perhaps d > best?
		if(d > best) SaveBestCMSA(cma);	// need this to avoid infinite loop.
	 	if(dd < 0.10){
		  dh_type dH=dheap(NumSeqsCMSA(cma)+5,4);
		  for(Int4 sq=1; sq <= NumSeqsCMSA(cma); sq++){
			if(MemberSet(sq,SetSq)) insrtHeap(sq,(keytyp)Random(),dH);
		  } gmb->SampleTogether(SetSq,Temp,lpr,dH); Nildheap(dH);  rcma=cma;
		  if(d < lpr) SaveBestCMSA(cma);  // sampling improved even more?
		  // if 2 seqs sampled back in same way then d==best; lpr < d == best;
		} else { rcma=gmb->SampleStickyTogether(rfp,MinSticky,Strategy,sfrq,Temp,StMd); cma=rcma; }
		gmb->RestoreBest();	// this restores the best; could be original configuration.
		lpr=gmb->RtnMap(); NilSet(SetSq);
		// WARNING: sequence weights could change with new gmb_typ...
		if(lpr > best){
		   did_mv=TRUE; best=lpr; improved=TRUE; 
		   delete gmb0; gmb0=0; NilCMSA(bcma); bcma=cma; cma=0; SaveBestCMSA(bcma); 
		   gmb0=new gmb_typ(aa_per_io,aa_per_do,exp_ie,exp_de,bcma,dms_mode,prior_wt);
		   gmb0->SetPriorWt(prior_wt); rcma=bcma;
		} else { if(gmb) delete gmb; gmb=0; NilCMSA(cma); cma=rcma=0; break; } // leave do loop...
          } else NilSet(SetSq);
	} while(rcma); 	// avoid infinite loops...need to see which column moved...
	Int4 cycle=Cycle; Cycle=printCycle;
	char	index=' ';
	if(did_mv) index='o';
	if(improved){ routine='m'; if(write) this->DoWrite(bcma); }
	fprintf(rfp," %c MvColumns: %d cycles %d secs",index,cycle,time(NULL)-timeR); 
	fprintf(rfp," (%0.2f mins)\n",(float)(time(NULL)-timeR)/60.0); fflush(rfp);
	if(gmb0) delete gmb0;
	if(gmb) delete gmb; gmb=0; if(cma) NilCMSA(cma); cma=bcma;
	return did_mv;
}

void	gpr_typ::DoWrite(cma_typ cma)
// print out intermediate files...
{
	// sprintf(str,"%s%dstg_gmb%d",name,iter,stage); RenameCMSA(str,cma);
	sprintf(str,"%s%d%c%d_gmb%d",name,iter,routine,Cycle,stage); RenameCMSA(str,cma);
	FILE *ofp=open_file(str,".cma","w"); PutCMSA(ofp,cma); fclose(ofp);
	sprintf(str,"%s_gmb_tmp",name); 
}

#define	UseStdSSX 0
Int4	gpr_typ::DoAddRmColumns(FILE *rfp, cma_typ &cma,Int4 MaxCycles)
{
     char	*ao,change=' ';
     double 	pn,d,dd,lpr,best;
     Cycle=0; sprintf(str,"%s_gmb_tmp",name);
     char	Strategy[20]=" RTNPXRVCDSW   ";	// stage == 1 (purged set).
                            //   R    T    N    P    X    R    V    C    D    S    W               
     double	sfrq[20]={ 0.0,0.20,0.10,0.10,0.30,0.15,0.20,0.20,0.20,0.08,0.15,0.20,0.0,0.0,0.0,0.0};

     char	StrategyA[20]=" RXTSDNRV  ";	// stage == 2(full).
                             //   R    X    T    S    D    N    R    V
     double	sfrqA[20]={ 0.0,0.15,0.10,0.10,0.10,0.08,0.10,0.10,0.15,0.0,0.0};

     char	StrategyB[20]=" RXVSWT  ";  	// stage == 3..4(full).
                             //   R    X    V    S    W    T      
     double	sfrqB[20]={ 0.0,0.10,0.15,0.15,0.10,0.15,0.05,0.0,0.0};
     ssx_typ	*std_ssx=0,*ssx=0;

     Int4	nAddRms=0,old_ncols;
     BooLean	keep_going,did_rm,did_add;
     cma_typ	std_cma=0,rcma=0,bcma=cma; cma=0; SaveBestCMSA(bcma);
#if UseStdSSX
     std_cma=CopyCMSA(bcma); pn=PerNatsCMSA(std_cma);
     std_ssx= new ssx_typ(aa_per_io,aa_per_do,exp_ie,exp_de,pn,std_cma,dms_mode);
     std_ssx->SetPriorWt(prior_wt); 
#endif
     gmb_typ	*gmb0=new gmb_typ(aa_per_io,aa_per_do,exp_ie,exp_de,bcma,dms_mode,prior_wt,std_ssx);

     for(keep_going=TRUE,Cycle=1; keep_going && Cycle <= MaxCycles; Cycle++){
       keep_going=did_add=FALSE; timeR=time(NULL); old_ncols = LengthCMSA(1,bcma);
	Int4 ntimes = 0;
       do {	//-------------- 3d.i. Add columns... -------------------
	// best=gmb0->RtnMap(); // weights can change...
	rcma=gmb0->AddColumns(bild_add,TRUE,0,afrq); 	// IMPORTANT: allow up to 50% gaps!
	ao=gmb0->RtnAddOp(); if(ao){ free(ao); } // sets AddOp=0 for gmb0;
       	if(rcma){
	   ntimes++;
	   if(gmb) delete gmb; gmb=0; 	 	// delete prior to NilCMSA!!
	   if(rcma != cma){ if(cma) NilCMSA(cma); cma=rcma; } rcma=0; 
#if 1	// avoid very long runs... 12/23/2015...
	   if(ntimes > 3) break;
#endif
#if UseStdSSX
           gmb=new gmb_typ(aa_per_io,aa_per_do,exp_ie,exp_de,cma,dms_mode,prior_wt,std_ssx);
#else
           gmb=new gmb_typ(aa_per_io,aa_per_do,exp_ie,exp_de,cma,dms_mode,prior_wt,gmb0->RtnSSX());
#endif
	   RenameCMSA(str,cma); gmb->SetPriorWt(prior_wt); gmb->SetMaxIter(2); SaveBestCMSA(cma);
	   best=gmb->RtnMap();
     	   if(stage < 2){
	   	if(Cycle <= 1){ 
		   // rcma=gmb->SampleStickyTogether(rfp,MinSticky,0.10,0,Temp,'R'); cma=rcma;
		} rcma=gmb->SampleStickyTogether(rfp,MinSticky,Strategy,sfrq,Temp,StMd);
     	   } else if(stage < 3){ rcma=gmb->SampleStickyTogether(rfp,MinSticky,StrategyA,sfrqA,Temp,StMd);
	   } else { rcma=gmb->SampleStickyTogether(rfp,MinSticky,StrategyB,sfrqB,Temp,StMd); }
#if 0	// Adjust for sequence weights...
	   gmb->RestoreBest(); cma=rcma; lpr = gmb->RtnMap(); fflush(rfp); 
	   ssx=gmb0->RtnSSX();
	   Int4	nWtSq,MaxWtSq=ssx->RtnWtNumSeqs();
	   ssx=gmb->RtnSSX();
	   if((nWtSq=ssx->RtnWtNumSeqs()) > MaxWtSq) MaxWtSq=nWtSq;
	   dd=ssx->AdjstDirichletLLR(MaxWtSq);    // Adjust for different # weighted seqs.
	   ssx=gmb0->RtnSSX();
	   d=ssx->AdjstDirichletLLR(MaxWtSq);    // Adjust for different # weighted seqs.
	   if(dd > d){
	   }
#else
	   // gmb->RestoreBest(); assert(cma==rcma); lpr = gmb->RtnMap(); fflush(rfp); 
	   gmb->RestoreBest(); cma=rcma; lpr=gmb->RtnMap(); fflush(rfp); 
#endif
	   if(lpr > best){
       		fprintf(rfp,"    columns: %d --> %d; LLR: %.2f --> %.2f\n",
				LengthCMSA(1,bcma),LengthCMSA(1,cma),best,lpr);
		keep_going=TRUE; did_add=TRUE; best=lpr; nAddRms++;
		delete gmb0; gmb0=0; NilCMSA(bcma); bcma=cma; cma=0;
#if UseStdSSX
	        gmb0=new gmb_typ(aa_per_io,aa_per_do,exp_ie,exp_de,bcma,dms_mode,prior_wt,std_ssx);
#else
	        gmb0=new gmb_typ(aa_per_io,aa_per_do,exp_ie,exp_de,bcma,dms_mode,prior_wt,gmb->RtnSSX());
#endif
	   } else { if(gmb) delete gmb; gmb=0; NilCMSA(cma); cma=rcma=0; } // leave do loop...
	}
       } while(rcma);
       if(did_add) change='+'; else { change=' '; }
       fprintf(rfp," %c AddColumns(%d->%d): bild add %.1f frq %.2f. %d secs",
		change,old_ncols,LengthCMSA(1,bcma),bild_add,afrq,time(NULL)-timeR); 
       fprintf(rfp," (%0.2f mins)\n",(float)(time(NULL)-timeR)/60.0); fflush(rfp);
       routine='a'; if(did_add && write) this->DoWrite(bcma); fflush(rfp); 

       // if(Cycle >= MaxCycles) break;

       did_rm=FALSE; timeR=time(NULL); old_ncols = LengthCMSA(1,bcma);
       do {	//-------------- 3d.ii. Remove columns... -------------------
	// best=gmb0->RtnMap();
        rcma=gmb0->RmColumns(bild_cut,TRUE,rfrq); 	// remove columns with > 50% deletions.
	ao=gmb0->RtnAddOp(); if(ao){ free(ao); } // sets AddOp=0 for gmb0;
        if(rcma){
	   if(gmb) delete gmb; gmb=0; 
	   if(rcma != cma){ if(cma) NilCMSA(cma); cma=rcma; } rcma=0; 
	   SaveBestCMSA(cma);
#if UseStdSSX
	   gmb=new gmb_typ(aa_per_io,aa_per_do,exp_ie,exp_de,cma,dms_mode,prior_wt,std_ssx); 
#else
	   gmb=new gmb_typ(aa_per_io,aa_per_do,exp_ie,exp_de,cma,dms_mode,prior_wt,gmb0->RtnSSX()); 
#endif
	   RenameCMSA(str,cma); gmb->SetMaxIter(2); 
	   best=gmb->RtnMap();	// LPR will go down if remove columns!!!
     	   if(stage < 2){ rcma=gmb->SampleStickyTogether(rfp,MinSticky,Strategy,sfrq,Temp,StMd);
     	   } else if(stage < 3){ rcma=gmb->SampleStickyTogether(rfp,MinSticky,StrategyA,sfrqA,Temp,StMd);
	   } else { rcma=gmb->SampleStickyTogether(rfp,MinSticky,StrategyB,sfrqB,Temp,StMd); }
	   // gmb->RestoreBest(); assert(cma==rcma); lpr = gmb->RtnMap(); fflush(rfp); 
	   gmb->RestoreBest(); cma=rcma; lpr = gmb->RtnMap(); fflush(rfp); 
	   if(lpr > best){
       		fprintf(rfp,"    columns: %d --> %d; LLR: %.2f --> %.2f\n",
				LengthCMSA(1,bcma),LengthCMSA(1,cma),best,lpr);
		keep_going=TRUE; did_rm=TRUE; best=lpr; nAddRms++; 
		delete gmb0; gmb0=0; NilCMSA(bcma); bcma=cma; cma=0;
#if UseStdSSX
	        gmb0=new gmb_typ(aa_per_io,aa_per_do,exp_ie,exp_de,bcma,dms_mode,prior_wt,std_ssx);
#else
	        gmb0=new gmb_typ(aa_per_io,aa_per_do,exp_ie,exp_de,bcma,dms_mode,prior_wt,gmb->RtnSSX());
#endif
		gmb0->SetPriorWt(prior_wt); 
		char *ao=gmb0->RtnAddOp(); if(ao){ free(ao); }
	   } else { if(gmb) delete gmb; gmb=0; NilCMSA(cma); cma=0; rcma=0; } // leave do loop...
	}
       } while(rcma);
       if(did_rm) change='-'; else { change=' '; }
       fprintf(rfp," %c RmColumns(%d->%d): bild rm %.1f frq %.2f. %d secs",
		change,old_ncols,LengthCMSA(1,bcma),bild_cut,rfrq,time(NULL)-timeR); 
       fprintf(rfp," (%0.2f mins)\n",(float)(time(NULL)-timeR)/60.0); fflush(rfp);
       if(did_rm){ routine='r'; if(did_rm && write) this->DoWrite(bcma); }
     }
     if(gmb0) delete gmb0; 
     if(std_ssx) delete std_ssx;
     if(std_cma) NilCMSA(std_cma);  
     if(gmb) delete gmb; gmb=0; if(cma) NilCMSA(cma); cma=bcma; 
     // return rcma;
     return nAddRms;
}

#define DBG_Pts 1	// find out why adding so many columns...
Int4	gpr_typ::DoAddColumns(FILE *rfp, cma_typ &cma)
{
     char	*ao,change=' ';
     double 	d,dd,lpr,best;
     Cycle=0; sprintf(str,"%s_gmb_tmp",name);
     char	Strategy[20]=" RXSTWN   ";
                            //   R    X    S    T    W    N                  
     double	sfrq[20]={ 0.0,0.10,0.10,0.20,0.10,0.08,0.10,0.0,0.0};

     Int4	nAddRms=0,old_ncols,Cycle=0;
     BooLean	did_add;
     cma_typ	rcma=0,bcma=cma; cma=0; SaveBestCMSA(bcma);
     gmb_typ	*gmb0=new gmb_typ(aa_per_io,aa_per_do,exp_ie,exp_de,bcma,dms_mode,prior_wt);

     did_add=FALSE; timeR=time(NULL); old_ncols = LengthCMSA(1,bcma);
     do {	//-------------- 3d.i. Add columns... -------------------
	Cycle++;
	if(Cycle > 5){ if(gmb) delete gmb; gmb=0; if(cma) NilCMSA(cma); cma=rcma=0; } // quit...
	else {
#if DBG_Pts
	  best=gmb0->RtnMap();
#else
#endif
	  rcma=gmb0->AddColumns(bild_add,TRUE,0,afrq); 	// IMPORTANT: allow up to 50% gaps!
	  ao=gmb0->RtnAddOp(); if(ao){ free(ao); } // sets AddOp=0 for gmb0;
       	  if(rcma){
	   if(gmb) delete gmb; gmb=0; nAddRms++;  	// delete prior to NilCMSA!!
	   if(rcma != cma){ if(cma) NilCMSA(cma); cma=rcma; } rcma=0; 
           gmb=new gmb_typ(aa_per_io,aa_per_do,exp_ie,exp_de,cma,dms_mode,prior_wt,gmb0->RtnSSX());
#if DBG_Pts
#else
	   best=gmb->RtnMap();
#endif
	   RenameCMSA(str,cma); gmb->SetMaxIter(2); SaveBestCMSA(cma);
     	   rcma=gmb->SampleStickyTogether(rfp,MinSticky,Strategy,sfrq,Temp,StMd);
	   gmb->RestoreBest(); cma=rcma; lpr=gmb->RtnMap(); fflush(rfp); 
	   if(lpr > best){
		did_add=TRUE; best=lpr;
		delete gmb0; gmb0=0; NilCMSA(bcma); bcma=cma; cma=0;
	        gmb0=new gmb_typ(aa_per_io,aa_per_do,exp_ie,exp_de,bcma,dms_mode,prior_wt,gmb->RtnSSX());
		gmb0->SetPriorWt(prior_wt); 
	   } else { if(gmb) delete gmb; gmb=0; NilCMSA(cma); cma=rcma=0; } // leave do loop...
	  }
	}
     } while(rcma);
     if(did_add) change='+'; else { change=' '; }
     fprintf(rfp," %c AddColumns(%d->%d): bild add %.1f frq %.2f. %d secs",
		change,old_ncols,LengthCMSA(1,bcma),bild_add,afrq,time(NULL)-timeR); 
     fprintf(rfp," (%0.2f mins)\n",(float)(time(NULL)-timeR)/60.0); fflush(rfp);
     routine='A'; if(did_add && write) this->DoWrite(bcma); fflush(rfp); 
     if(gmb0) delete gmb0;  
     if(gmb) delete gmb; gmb=0; if(cma) NilCMSA(cma); cma=bcma; 
     return nAddRms;
}

