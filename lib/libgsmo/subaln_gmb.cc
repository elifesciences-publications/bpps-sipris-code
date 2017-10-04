#include "gmb_typ.h"

void	PutSmplSubAlnsCMSA(char *name, cma_typ  cma, cma_typ rcma, cma_typ CMA,
			cma_typ finalCMA, set_typ *Set, Int4 XX, char *AddOp, char dms_mode)
{
	Int4    i,s,aa_per_io=1000,aa_per_do=1000,exp_ie=1,exp_de=1;
	char	str[200];
        a_type  AB=AlphabetCMSA(cma);
        FILE    *fp;
	//--------------------- output subalign to debug -------------------------
	sprintf(str,"%s_%d",name,XX); 
	FILE *wfp=open_file(str,".cma","w");
	ReNameCMSA("AlnA",cma); PutInSetCMSA(wfp,Set[XX],cma); 
	// fclose(wfp); sprintf(str,"%s_%do",name,XX); wfp=open_file(str,".cma","w");
	ReNameCMSA("AlnB",CMA); PutInSetCMSA(wfp,Set[XX],CMA); fclose(wfp);

	e_type CsqE=GetSeqAsCsqCMSA(Set[XX],CMA);
	e_type csqE=GetSeqAsCsqCMSA(Set[XX],cma);
	gmb_typ *xgmb=new gmb_typ(aa_per_io,aa_per_do,exp_ie,exp_de,CMA,dms_mode); 
	FILE *afp = open_file(str,".aln","w");
	// PutSubGroupFootPrints(afp,Set[XX],CMA);
	PutSubGroupCsqScores(afp,Set[XX],CMA);
	// PutSubGroupCsqScores(afp,Set[XX],cma);
	Int4 Start;
	char *op=xgmb->GlobalAlignSeq(csqE,Start);
        fprintf(afp,"\n%s\n",op);
        xgmb->PutGappedAln(afp,csqE,op,Start); free(op); delete xgmb;
	AlnSeqSW(afp,11,1,CsqE, csqE,AB);
	if(AddOp){
	      xgmb=new gmb_typ(aa_per_io,aa_per_do,exp_ie,exp_de,cma,dms_mode); 
	      fprintf(afp,"%s",AddOp); fflush(afp);
              xgmb->PutGappedAln(afp,CsqE,AddOp,1); delete xgmb;
	} NilSeq(csqE); NilSeq(CsqE); fclose(afp);
#if 1
	if(AddOp){
	      if(rcma){
        	wfp=open_file(str,"_gld.cma","w");
          	PutInSetCMSA(wfp,Set[XX],rcma); fclose(wfp);
        	wfp=open_file(str,"_full.cma","w");
          	PutInSetCMSA(wfp,Set[XX],cma); fclose(wfp);
	      }
              FILE *wfp=open_file(str,"_gold.cma","w");
              PutInSetCMSA(wfp,Set[XX],finalCMA); fclose(wfp);
	}
#endif
}

#if 0
// can turn this back on if I decide to use this...
cma_typ	SampleSubAlignsCMSA(char *name, cma_typ  cma, char dms_mode)
{
	//===================== 0. Set parameters ========================
	Int4    i,j,s,similarity=0;
        Int4    aa_per_io=1000,aa_per_do=1000,exp_ie=1,exp_de=1;
        char    method='S',mode='S',str[200];
        // method = alignment method: 'E' = gapped blast; 'S' = \% similarity (default S).
        a_type  AB=AlphabetCMSA(cma);
        FILE    *fp;
        double  temp=-1.0,prior_wt=0,AdjSqWt=0.0;
	Int4	percent_id=40;
	double	max_fract=0.15;
	FILE *rfp=open_file(name,".log","w");
	cma_typ	xcma,rcma=0;
        double bild_cut=0.0,BildCut=0.0,dd,MaxFrctn=0.33;
	Int4	sq=1,N=NumSeqsCMSA(cma),MinSticky=2;
	Int4	numSets,*oldpos,Begin;

	//================== 1. Find related subsets for subalignment =====================
#if 1	// default...
	set_typ *CloseSq=RtnTargetSizeSetsCMSA(numSets,percent_id,cma,max_fract);
	Begin=1;
#else	// using omcBPPS defined sets...
        FILE *sfp = open_file(name,".sets","r");
	assert(sfp != 0);
	set_typ *CloseSq=ReadSets(sfp,numSets); fclose(sfp);
	Begin=2;
#endif
	//===================== 2. Subalign the larger subsets ========================
	Int4	MinSubSetSize=10;
	cma_typ finalCMA=CopyCMSA(cma);
	for(Int4 XX=Begin; XX <= numSets; XX++){
	  if(CardSet(CloseSq[XX]) < MinSubSetSize) break;
	  cma_typ CMA=CopyCMSA(cma);
          gmb_typ *gmb=new gmb_typ(aa_per_io,aa_per_do,exp_ie,exp_de,cma,dms_mode); 
	  if(prior_wt > 0) gmb->SetPriorWt(prior_wt); gmb->SetMaxIter(3);
	  //------------ 2a. Find additional sequences closely related to this set -----------------
	  Int4 *TestSet=gmb->RtnTestSet(CloseSq[XX]);

	  //----------------- 2b. Retain only sequences in the set -------------------
  	  ssx_typ *SSX=gmb->RtnSSX(); 
	  for(sq=1; sq<= N; sq++){
	    if(MemberSet(sq,CloseSq[XX])) continue;
	    oldpos=SSX->RemoveFromAlign(sq); assert(oldpos != 0); free(oldpos);
  	  } 
	  //------------------ 2b. Iteratively add columns & realign ----------------------
	  char *AddOp=0; sprintf(str,"%s_gmb",name); 
	  do {
	     RenameCMSA(str,cma); BildCut=-1.0;
	     // gmb->SampleStickyTogether(rfp,MinSticky,MaxFrctn,str);
	     // gmb->SampleStickyTogether(rfp,MinSticky,MaxFrctn,str,0.0,'W'); gmb->Sample(str,'S',similarity);

	     //.................. 2b.i. Add columns .....................
             rcma=gmb->AddColumns(BildCut,TRUE,CloseSq[XX]);	// widen alignment...
             // rcma=gmb->AddColumns(BildCut,FALSE,CloseSq[XX]);	// don't widen alignment...
             if(rcma){	// if columns were added...then realign...
		fprintf(rfp,"%d.add: %d -> %d columns; %d -> ",
			XX,LengthCMSA(1,cma),LengthCMSA(1,rcma),CardSet(CloseSq[XX]));
		//.................. Save added column positions .....................
		if(AddOp) free(AddOp); AddOp=gmb->RtnAddOp(); // keep track of added columns..
		delete gmb; NilCMSA(cma); cma=rcma;  // rcma=0;
		//.................. Add related sequences not in the set .................
                gmb=new gmb_typ(aa_per_io,aa_per_do,exp_ie,exp_de,cma,dms_mode); 
		if(prior_wt > 0) gmb->SetPriorWt(prior_wt); gmb->SetMaxIter(3); 
   		gmb->PutSampleOutSqIn(stderr,0.0,CloseSq[XX],TestSet);	// can add to CloseSq set.
		fprintf(rfp,"%d seqs.\n",CardSet(CloseSq[XX])); fflush(rfp);
		//.................. realign by sampling .....................
		gmb->SetAddOp(AddOp); 
		if(CardSet(CloseSq[XX]) >= 50){
	        	gmb->SampleStickyTogether(rfp,MinSticky,MaxFrctn,str); fflush(rfp);
		} gmb->Sample(str,'S',0,0.0);
#if 0	// WARNING: need to keep track of removed columns...
		//.................. remove poor columns? .....................
		bild_cut=0.0; xcma=gmb->RmColumns(bild_cut);
                if(xcma){
		   delete gmb; NilCMSA(cma); rcma=cma=xcma; 
                   gmb=new gmb_typ(aa_per_io,aa_per_do,exp_ie,exp_de,cma,dms_mode); 
		   if(prior_wt > 0) gmb->SetPriorWt(prior_wt); gmb->SetMaxIter(3); 
		}
#endif
             }
          } while(rcma);
	  // gmb->SetMaxIter(3); gmb->Sample(str,'S',similarity,0.0);
	  delete gmb; free(TestSet);

	  //-------------- 2c. Revert to original columns & move gsq to final alignment -----------------
	  if(AddOp){
	    rcma=RemoveTheseColumnsCMSA(AddOp,cma); 
	    MoveAlignmentCMSA(CloseSq[XX],rcma,finalCMA);	
	  } else { rcma=0; }

	  //------------------- 2d. Output subalign to debug -----------------------
#if 0
	  PutSmplSubAlnsCMSA(name,cma,rcma,CMA,finalCMA,CloseSq,XX, AddOp, dms_mode);
#endif
	  if(rcma) NilCMSA(rcma); rcma=0;
	  NilCMSA(cma); cma=CMA; // revert to saved CMA...
   	} fclose(rfp);
#if 0	// moved to calling environment...
	sprintf(str,"%s_final.cma",name); WriteCMSA(str,finalCMA); 
// exit(1);
	gmb_typ *gmb=new gmb_typ(aa_per_io,aa_per_do,exp_ie,exp_de,finalCMA,dms_mode); 
	//============ sample sequences in final alignment one-at-a-time ================
	gmb->Sample(str,'S',similarity,150.0); delete gmb;
	sprintf(str,"%s_finalA.cma",name); WriteCMSA(str,finalCMA); 
#endif
	return finalCMA;
}
#endif

