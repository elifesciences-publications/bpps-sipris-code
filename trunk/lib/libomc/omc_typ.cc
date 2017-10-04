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

#include "omc_typ.h"

void	omc_typ::PutCheckPoint(char *string,BooLean change)
// main routine to perform operations on mcs_typ mcs.
{
	Int4	*P=0,id,sq,i,j,pI,I,J,r,R,m,n,p,x,cID=0;
	Info	*info = new Info; info->Init(MaxNumNodesPlus); 
	// NEW(info->sma,MaxNumNodesPlus + 3,cma_typ); // created in info->Init();
	hpt_typ	*hpt=0;
	char	*name=AllocString(string);
	sst_typ **xsst=mcs->RtnCopyOfSSTs();
	char    **pttrn; NEWP(pttrn,MaxNumNodesPlus+3,char);

	set_typ *set=mcs->CopyOfBestTreeSets();
	// NEW(info->set,MaxNumNodesPlus + 3,set_typ);  // created by info->Init();
	hpt=Hpt->Copy();  info->hpt=hpt;  info->T=0; 
// fprintf(stderr,"taking checkpoint\n");
	for(j=1; j < hpt->NumSets(); j++){
		if(change) {
		   id=Hpt->ItoSetID(j); sprintf(str,"Set%d",id);
		} else { sprintf(str,"%s",Hpt->ElmntSetName(j)); }
		info->sma[j]=mcs->RtnCsqAsCMSA(j,str);  // won't redefine xsst[i].
		pttrn[j]=SST2ArgStrAlpha(xsst[j],mcs->RtnLengthMainCMSA(),AB);
		info->set[j]=set[j];
	} info->set[hpt->NumSets()] = set[j]; // info->set=set; set=0;
	for(i=1; i <= Hpt->NumBPPS(); i++){ if(xsst[i]) free(xsst[i]); } free(xsst); 
     	for(i=1; i <= hpt->NumBPPS(); i++){
// fprintf(stderr,"pttrn[%d]=%s\n",i,pttrn[i]);
            if(pttrn[i]){
		char *tmp[3]; tmp[0]=pttrn[i]; hpt->ReSetArgv(i,1,tmp); free(pttrn[i]); 
		if(change){
		   id=hpt->ItoSetID(i); sprintf(str,"Set%d",id); hpt->ReNameGroup(i,str);
		} 
	    }
     	} 
	// hpt->Put(stderr,TRUE,FALSE,TRUE);
// fprintf(stderr,"DEBUG: file name = %s\n",name);

        FILE *fp=0;
#if 1	// *.chk single output file.
        fp=open_file(name,".chk","w"); 
	{
	  bpcp_typ bpcp(RandomSeed,hpt->NumSets(),info->set,hpt->NumSets()-1,info->sma,hpt);
	  bpcp.Write(fp);
	} fclose(fp); fp=0; 
#else	// output old format files.
	fp=open_file(name,".hpt","w"); hpt->Put(fp,TRUE,FALSE,TRUE); fclose(fp); fp=0; 
        fp=open_file(name,".sets","w"); WriteSets(fp,hpt->NumSets(),info->set); fclose(fp); fp=0; 
        fp=open_file(name,".sma","w");
	for(i=1; i < info->hpt->NumSets(); i++){ PutCMSA(fp,info->sma[i]); } fclose(fp); fp=0; 
        fp=open_file(name,".seed","w"); fprintf(fp,"-seed=%d\n",RandomSeed); fclose(fp); fp=0;
#endif
	if(info->hpt){ 
	  if(info->sma){
		for(i=1; i < info->hpt->NumSets(); i++){
		   if(info->sma[i]) TotalNilCMSA(info->sma[i]); info->sma[i] = 0;
		} 
	  } if(info->set){
		for(i=1; i <= info->hpt->NumSets(); i++){
		   if(info->set[i]) NilSet(info->set[i]); info->set[i] = 0;
		}
	  } delete info->hpt; 
	} info->Free(); delete info; 
	if(set) free(set); free(pttrn); if(P) free(P); free(name);
}

Int4	omc_typ::DeletableNodes(mcs_typ *x_mcs,double minLLR, Int4 minCard)
// Delete nodes with LLRs < minLLR
{
	Int4	i,Ndel=0;
	x_mcs->RestoreBest();
	double  *Map=x_mcs->RtnSubLPR( );
	for(i=2; i <= Hpt->NumBPPS(); i++){	// 1. Eliminate empty leaf nodes...
	    if(Hpt->TypeOfSet(i) != '!') continue;
	    if(x_mcs->SetCard(i) < minCard){ Ndel++; }
	    else if(Map[i] < minLLR){ Ndel++; }
	} return Ndel; 
}

void    omc_typ::PutAsBestHpt(mcs_typ *x_mcs)
{
	assert(x_mcs != 0);
	FILE *fp=open_file(infile,"_bst.out","w"); x_mcs->PutHyperPartition(fp); fclose(fp);
}

Int4	omc_typ::ToggleOutFileName(BooLean NoSubID)
{
      const char letter[]=" ABCDEFGHIJKLMNOPQRSTUVWXYZ   ";
      if(outfilename) free(outfilename); Iteration++;
      Int4 id = delminHeap(FileIdHeap); assert(id > 0 && id <= 26);
      if(NoSubID) outfilename=AllocString(infile);
      else {
        sprintf(str,"%s_%c",infile,letter[id]); outfilename=AllocString(str);
      } return  id;
}

void    omc_typ::DeleteMCS(mcs_typ *x_mcs)
{
        Int4 id = x_mcs->FileID; assert(id >= 0 && id <=26);
        if(mcs == x_mcs){ mcs=0; Hpt=0; }
        if(id > 0){  // if id == 0 then deleting current best_mcs.
              assert(!memHeap(id,FileIdHeap)); insrtHeap(id,id,FileIdHeap);
        } delete x_mcs;
}

void    omc_typ::RestoreFinalBest( )
// Restore best and set best_mcs == 0.
{
       CalcLLR(mcs); // in case current is better than best_cms
       if(best_mcs != 0){
	  fprintf(stderr,"Restoring best_mcs and resetting to null (%.2f =? %.2f)\n",
		best_mcs->CalcTotalLPR( ),Best_LPR);
          DeleteMCS(mcs); mcs=best_mcs; Hpt=mcs->GetHpt( ); best_mcs=0; Best_LPR=0.0;
          mcs->StoreBest(); mcs->RestoreBest();
       } else print_error("Failed to find a significant hierarchy.");
}

double  omc_typ::RestoreVeryBest( )
// Restore best but save best_mcs.
{
       CalcLLR(mcs); // in case current is better than best_cms
       if(best_mcs != 0){
	  // sprintf(str,"%s_chk",infile); PutCheckPoint(str); // DEBUG...
	  fprintf(stderr,"Restoring best_mcs but retaining it (%.2f =? %.2f)\n",
					best_mcs->CalcTotalLPR( ),Best_LPR);
          DeleteMCS(mcs); mcs=best_mcs; Hpt=mcs->GetHpt( );
          mcs->RestoreBest(); mcs->StoreBest();
          mcs_typ *tmp_mcs = this->RtnCopy('C'); // make a copy...
          best_mcs=mcs; mcs=tmp_mcs; Hpt=mcs->GetHpt( );
	  assert(Best_LPR == best_mcs->CalcTotalLPR( ));
          return mcs->RtnBestLPR();
       } else mcs->RestoreBest();
}

double  omc_typ::CalcLLR(mcs_typ *&xmcs){
      xmcs->RestoreBest();
      double d=xmcs->CalcTotalLPR( );
      if(d > Best_LPR){
        Int4 n=xmcs->NumRawFailedNodes(MinimumSetSize,MinimumLLR);
	if(n == 0){
         fprintf(stderr,"Copying mcsID=%d to best_mcs (LLR=%.2f).\n",xmcs->FileID,d);
	 BooLean checkthis=TRUE;
	 if(best_mcs == 0) checkthis=FALSE; // samples negative columns at this point???
#if 0	// I believe this generates the temporary files *_A.out, *_B.out etc.
         xmcs->PutHyperPartition(); // report...
#endif
#if 0
         mcs_typ *tmp=0;
         if(mcs != xmcs){ tmp=mcs; mcs=xmcs; Hpt=mcs->GetHpt( ); }
	 // Hpt->PutAsTree(FILE *fp);
         if(best_mcs != 0) DeleteMCS(best_mcs); // need to delete first!!!
         best_mcs=Operate('B',0,0); // == RtnCopy('B'); 
	 assert(best_mcs != 0); Best_LPR=d; 
         best_mcs->PutHyperPartition( ); // report... <outfilename>_bst.out
	 best_mcs->PutCheckPoint( );
         if(tmp){ xmcs=mcs; mcs=tmp; Hpt=mcs->GetHpt( ); }
	 if(checkthis && !CheckBestCopy(xmcs)){
		best_mcs->PutHyperPartition(stderr);
		xmcs->PutHyperPartition(stderr);
		print_error("omc_typ::CalcLLR() error: best_mcs check failed!");
	 } fprintf(stderr,"Done creating mcs (FileID=%d) as best.\n",best_mcs->FileID);
#else
	 if(CheckSSTvsCSQ(xmcs)){
           if(best_mcs != 0) DeleteMCS(best_mcs); // need to delete first!!!
           if(mcs != xmcs){   //  xmcs is different from mcs; retain mcs as is.
	   	best_mcs=xmcs; Best_LPR=d; 
           	mcs_typ *tmp=mcs; mcs=best_mcs; Hpt=mcs->GetHpt( ); // Hpt->PutAsTree(FILE *fp);
	   	xmcs=Operate('B',0,0); // == RtnCopy('B'); 
		assert(xmcs != 0); mcs=tmp; // set xmcs equal to copy and mcs back to original.
	   } else {	// xmcs == mcs; set mcs = copy of itself!
	   	best_mcs=mcs; Best_LPR=d; mcs=Operate('B',0,0); // == RtnCopy('B'); 
		assert(mcs != 0); xmcs=mcs;
	   } Hpt=mcs->GetHpt( );  
	   if(checkthis && !CheckBestCopy(xmcs)){
		best_mcs->PutHyperPartition(stderr); xmcs->PutHyperPartition(stderr);
		print_error("omc_typ::CalcLLR() error: best_mcs check failed!");
	   } fprintf(stderr,"Done creating mcs (FileID=%d) as best.\n",best_mcs->FileID);
	   PutAsBestHpt(best_mcs);
           // best_mcs->PutHyperPartition( ); // report... <outfilename>_bst.out
#if 0	// I believe that this creates the *_A.hpt and *_A.cma intermediate files.
	   best_mcs->PutCheckPoint( );
#endif
	 } else print_error("omc_typ::CalcLLR() error: best_mcs pattern inconsistent with consensus!");
#endif
	} else {
	   // xmcs->PutHyperPartition(stderr);
	   fprintf(stderr,"LLR from %.2f to %.2f but with %d failed nodes (not saved)!\n", Best_LPR,d,n);
	}
      } return d;
}

mcs_typ	*omc_typ::CreateMCS(hpt_typ *hpt,set_typ *set, cma_typ *in_sma,
					double Temp,char mode, Int4 pID,Int4 cID)
// create and return a mcs_typ object for sampling...
// pID (if positive) is the ID for the set, over the subtree of which sampling should occur.
// mode: 'A' = AddLeaf; 'B' = StoreBest; 'C' = Copy; 'D' = DeleteNode; 'O' = OptimizeDisplay; 
// 'S' = SortHpt; 
{
     static	UInt4 calls=0; 
     Int4	id,x,i,j,*P;
     time_t	timeX=time(NULL);
     set_typ	st=0;
     mcs_typ	*rtn_mcs;

     if(mode != 'I' && mode != 'i') fprintf(stderr,"============== CreateMCS( ) '%c' operation ==============\n",mode);
     id=ToggleOutFileName(); 
     this->SetDefaultArguments( );
#if 1	// DEBUG.
     if(mode == 'B') mcs->PutHyperPartition(stderr);
#endif
     if(mode=='c' || mode=='b' || mode=='r'){
	 rtn_mcs = new mcs_typ(TrueMainCMA,MainCMA,MainHSW,hpt,in_sma,Argc,Argv);
     } else rtn_mcs = new mcs_typ(TrueMainCMA,MainCMA,MainHSW,hpt->NumSets(),set,hpt,in_sma,Argc,Argv);
     rtn_mcs->SetStringency(MinimumLLR,MinimumSetSize,stable);
// rtn_mcs->PutHyperPartition(stderr);
     rtn_mcs->UnLabelAllSeqs( );  // don't want this here...
     // rtn_mcs->DisOwnInSMA(); // let mcs_typ destructor free up in_sma.
     rtn_mcs->FileID=id; rtn_mcs->IsTreeHpt=TRUE; 
     rtn_mcs->DoEvolve();
     rtn_mcs->NoFailureMode=TRUE; // if node configuration subLPR <= 0.0, then keep nodes anyway.
                                   // Don't accept higher LPR due to node deletion along.
     rtn_mcs->SaveBest=TRUE;       // Start saving the best configuration immediately.
     switch (mode){
	case ' ': Sample(Temp,2,2,2,2,rtn_mcs,0.03); break; // Optimize over all nodes.
	case 'S': { Sample(Temp,2,2,2,2,rtn_mcs,0.03); } break; // created from a sorted Hpt.
	case 'a': {	
	    CheckID4NewMCS(rtn_mcs,pID); CheckID4NewMCS(rtn_mcs,cID); assert(hpt->IsTree(P)); free(P);
     	    st=MakeSet(hpt->NumSets()+1); ClearSet(st); 
	    i=hpt->SetIDtoI(pID); AddSet(i,st); i=hpt->SetIDtoI(cID); AddSet(i,st); 
	    i=hpt->NumSets(); AddSet(i,st); // sample over node + parent + random.
	    // if(mode != 'a') rtn_mcs->SampleColumns( ); 
#if 0	// DEBUG...
rtn_mcs->PutPttrnVsConsSeq(stderr,"omc->CreateMCS() debug 0");
rtn_mcs->CheckPttrnCsqMatch("omc->CreateMCS() debug 1");       // ZZZZZZZZZZZZ
fprintf(stderr,"******** CardSet(st) = %d *********\n",CardSet(st));
PutSet(stderr,st);
#endif
	    Sample(Temp,2,2,2,2,rtn_mcs,0.05,st); NilSet(st);
	    rtn_mcs->PutHyperPartition(stderr);
	  } break;
	case 'A': {	// sample over parent (pID), child (cID) (and, for 'A', reject) nodes only.
	    CheckID4NewMCS(rtn_mcs,pID); CheckID4NewMCS(rtn_mcs,cID); assert(hpt->IsTree(P)); free(P);
     	    st=MakeSet(hpt->NumSets()+1); ClearSet(st); 
	    i=hpt->SetIDtoI(pID); AddSet(i,st); i=hpt->SetIDtoI(cID); AddSet(i,st); 
	    // if(mode=='a'){ i=hpt->NumSets(); AddSet(i,st); } // sample over node + parent + random.
	    // Sample(Temp,2,3,2,2,rtn_mcs,0.05,st); NilSet(st);
	    // rtn_mcs->SampleColumns( ); Sample(Temp,2,2,2,2,rtn_mcs,0.05,st); NilSet(st);
	    rtn_mcs->SampleColumns( ); Sample(Temp,2,2,2,2,rtn_mcs,0.20,st); NilSet(st);
	  } break;
	case 'd': { rtn_mcs->SampleColumns( ); } break; 
	case 's':	// Sort with sampling over subtree only...
	case 'D': {		// sample over subtree of node pID.
	    i=hpt->SetIDtoI(pID);  CheckID4NewMCS(rtn_mcs,pID); assert(hpt->IsTree(P)); 
	    if(P[i] == 0){ st=0; }
     	    else { st=MakeSet(hpt->NumSets()+1); this->ReSetSubTree(st,i,hpt); }
	    // x=hpt->NumSets(); AddSet(x,st);
	    // Sample(Temp,2,2,2,2,rtn_mcs,0.03,st); if(st) NilSet(st); free(P);
	    Sample(Temp,2,2,2,2,rtn_mcs,1.5,st); if(st) NilSet(st); free(P);
	  } break; 
	case 'i': { 
	    CheckID4NewMCS(rtn_mcs,pID); CheckID4NewMCS(rtn_mcs,cID); assert(hpt->IsTree(P)); free(P);
     	    st=MakeSet(hpt->NumSets()+1); ClearSet(st); 
	    i=hpt->SetIDtoI(pID); AddSet(i,st); // i=hpt->SetIDtoI(cID); AddSet(i,st); 
     	    rtn_mcs->NoFailureMode=FALSE; // if node configuration subLPR <= 0.0, then keep nodes anyway.
	    rtn_mcs->SampleColumns( ); Sample(Temp,2,2,2,2,rtn_mcs,0.95,st); NilSet(st);
	    rtn_mcs->StoreBest(); rtn_mcs->NoFailureMode=TRUE;
	  } break; 
	case 'I': { 
     	   rtn_mcs->NoFailureMode=FALSE; // if node configuration subLPR <= 0.0, then keep nodes anyway.
	   rtn_mcs->StoreBest(); rtn_mcs->NoFailureMode=TRUE;
	  } break; 
	case 'r':  // rtn_mcs->SampleColumns( ); // rtn_mcs->Sample(2,2,2,2);
	case 'c':  // rtn_mcs->Sample(2,2,2,2);
	case 'b':  // rtn_mcs->Sample(2,2,2,2);
	case 'R': 
	case 'B': 	// copy for best.
	case 'u': 	// unscramble only!
	case 'O': 
	case 'C': {
     	   rtn_mcs->CalcTotalLPR(0,FALSE); // Need to Calculate LPR prior to storing TotalLPR as BestLPR;
	   // rtn_mcs->DidRestoreBest=FALSE;
     	   rtn_mcs->NoFailureMode=FALSE; // if node configuration subLPR <= 0.0, then keep nodes anyway.
	   rtn_mcs->StoreBest(); rtn_mcs->NoFailureMode=TRUE;
	   // rtn_mcs->CheckInputSets(); rtn_mcs->ChecksOut();  rtn_mcs->ConsistencyCheck(); // debug...
	   // rtn_mcs->PutHyperPartition(stderr);
	  } break; 
	default: break;	// do nothing...
     } rtn_mcs->RestoreBest(); 
     // rtn_mcs->PutHyperPartition(stderr);
     rtn_mcs->CalcTotalLPR( ); // Need to Calculate LPR prior to outputting PutHyperPartition();
     if(mode != 'I' && mode != 'i') fprintf(stderr,"   >>>>>>>> time = %d seconds <<<<<<<<\n",
			(Int4) difftime(time(NULL),timeX));
     return rtn_mcs;
}

BooLean	omc_typ::SampleMCS(mcs_typ *x_mcs,Int4 ID)
// Sample over the passed in mcs object and, if it is an improvement, use it.
{
	double	d,D;
	Int4	x,n,m;
	d = x_mcs->RtnBestLPR(); D = mcs->RtnBestLPR();
	fprintf(stderr,"============ Best Lpr=%.2f; Test Lpr=%.2f ============\n",D,d);
	if(d > D){
            mcs_typ *tmp_mcs=mcs; mcs=x_mcs; Hpt=mcs->GetHpt();
            if((x=this->DeleteNodes( )) > 0){ fprintf(stderr,"%d rejected nodes deleted.\n",x); }
            n=Hpt->SetIDtoI(ID);
fprintf(stderr,"Node %d ('Set%d') == 0 ?\n",n,ID);
	    // mcs->PutHyperPartition(stderr);
            if(n > 0){	// The Inserted node was not deleted; delete old mcs and use the new....
	      mcs->RestoreBest(); d=mcs->RtnBestLPR(); // mcs is new.
	      m=mcs->NumRawFailedNodes(MinimumSetSize,MinimumLLR);
	      if(d > D &&  m==0){
		// fprintf(stderr,"Sampling step (node %d) sucessful.\n",n);
		fprintf(stderr,"!!!!!!!!!!!!!!!!!! Sampling step (node %d) sucessful. !!!!!!!!!!!!!!!!!!!!\n",n);
		mcs->PutHyperPartition(stderr); 
#if 0	// I believe that this generates the temporary files... *_A.out *_B.out etc.
	     mcs->PutHyperPartition(); 
#endif
		DeleteMCS(tmp_mcs); 
		this->CalcLLR(mcs);
	        return TRUE; 
	      }
	    } // else the added node was deleted; revert back to original mcs (== tmp_mcs).
	    DeleteMCS(mcs); mcs=tmp_mcs; Hpt=mcs->GetHpt(); // mcs->PutHyperPartition(stderr); 
	    return FALSE; 
        } else { DeleteMCS(x_mcs); return FALSE; }
}

Int4	omc_typ::Put(BooLean SkipCheckPoint,FILE *fpmma,FILE *fphpt,FILE *ptrnfp)
{
	FILE *fp;

	fprintf(stderr,"############### Printing results ###############\n");
	mcs->RenameInfile(infile);	// use input file name without toggle letter.
        // mcs->PutHyperPartition( );
        if(SaveSets){		// save sets...
           set_typ *sets=mcs->CopyOfSeqSets();
           Int4 NumSets=mcs->RtnNumElmntSetCMA( );
	   // sprintf(str,"_%d.sets",Iteration);
           // fp=open_file(infile,str,"w");
	   if(verbose){ fp=open_file(infile,"_new.sets","w"); WriteSets(fp,NumSets,sets); fclose(fp); }
	   fp=0; 
	   for(Int4 i=1; i <= NumSets; i++) NilSet(sets[i]); free(sets);
	   if(verbose){ fp=open_file(infile,".cntrb","w"); mcs->PutMapContributions(fp); fclose(fp); }
        } // for further analysis...
	if(verbose){
	  fp=open_file(infile,"_ptn.lpr","w"); mcs->PutPttrnLLRs(fp); fclose(fp); fp=0;
	  // fp=open_file(infile,"_new.sma","w"); mcs->PutDisplayCMA(fp); fclose(fp); fp=0;
          // PutSet(stderr,DisplaySet[NewNode]); // exit(1);
	  fp=open_file(infile,".diagnose","w"); 
	  if(goodHG) PutHist(fp,60,goodHG);
	  if(badHG) PutHist(fp,60,badHG);
	  if(testHG) PutHist(fp,60,testHG);
 	  if(triesHG) PutHist(fp,60,triesHG);
	  fclose(fp); fp=0;
	}
	if(!SkipCheckPoint){ sprintf(str,"%s",infile); PutCheckPoint(str); }
	if(OutPutRTF){
	    mcs->SetFontSize(font_size); mcs->SetPageFormat(page_format);
	    if(NthSeqForDisplay > 0) mcs->SetNthSeqForDisplay(NthSeqForDisplay);
	    fprintf(stderr,"Printing rich text format contrast alignments...\n");
	}
#if 1
	char **names=0;
	if(verbose==FALSE) { mcs->VerboseOff(); }
	if(NameSeedAln){
	   std::cerr << NameSeedAln << std::endl;
	   names=GetNames( ); 
	   if(use_usr_sma) mcs->SetUserDisplaySet();
	   if(OutPutRTF) this->PrintRTF(FALSE);
	   mcs->Put(FALSE,names,FALSE,fpmma);  
	   hpt_typ *tmp_hpt=Hpt->Copy();
	   for(Int4 g=1; g <= Hpt->NumSets(); g++){
		if(names[g]){ tmp_hpt->ReNameSet(g,names[g]); free(names[g]); names[g]=0; }
	   }
	   if(fphpt){ tmp_hpt->Put(fphpt); }
	   else { fp=open_file(infile,"_new.hpt","w"); tmp_hpt->Put(fp,TRUE,FALSE,TRUE); fclose(fp); }
	   fp=0; delete tmp_hpt; 
	} else {
	   if(use_usr_sma) mcs->SetUserDisplaySet();
	   if(SecretCode || OutPutRTF) this->PrintRTF(FALSE);	// Always output if SecretCode > 0
	   if(SecretCode != 1){
	     if(SecretCode == 3) mcs->Put(FALSE,0,FALSE,fpmma,ptrnfp); 
	     else mcs->Put(TRUE,0,FALSE,fpmma);  // creates mcs_typ output files including <infile>_new.mma
	     if(fphpt){ Hpt->Put(fphpt); }
	     else { fp=open_file(infile,"_new.hpt","w"); Hpt->Put(fp,TRUE,FALSE,TRUE);  fclose(fp); fp=0; }
	   }
	   // fp=open_file(infile,"_tree.hpt","w"); Hpt->PutSorted(fp);  fclose(fp); fp=0;
	}
#else
	if(use_usr_sma) mcs->SetUserDisplaySet();
	if(OutPutRTF) this->PrintRTF(FALSE);
	mcs->Put(FALSE);  // creates mcs_typ output files including <infile>_new.mma
#endif
        fprintf(stderr,"done printing results\n");
}

char	**omc_typ::GetNames()
// Get names corresponding to a curated hierarchy; pick highest scoring consensus seq to name.
{
	char	**names=0;
	Int4	Number=0,n,g;
	Int4	sq,N=NumSeqsCMSA(TrueMainCMA);
	if(NameSeedAln == 0) return 0;
	FILE *fp=open_file(NameSeedAln,"","r");
	cma_typ *SeedCMA=MultiReadCMSA(fp,&Number,AB); fclose(fp);
	Int4 *BestScore=0,*BestG=0,*BestScoreG=0;
	NEW(BestScoreG,Number +5,Int4); NEW(BestScore,Number +5,Int4); NEW(BestG,Number +5,Int4);
	NEWP(names,Hpt->NumSets()+5,char);
	names[1]=AllocString(NameCMSA(SeedCMA[1]));	// name the root as root!
	for(g=2; g <= Hpt->NumSets(); g++){
	  for(sq=1; sq <= N; sq++){
	     if(!mcs->IsInSet(sq,g)) continue;
	     for(n=1; n <= Number; n++){
		Int4 Score=PseudoAlnScoreTwoCMSA(1,SeedCMA[n],sq,TrueMainCMA);
                if(Score > BestScore[n]){ BestScore[n]=Score; BestG[n]=g; }
	     }
          } 
	}
	for(n=2; n <= Number; n++){
	    if(BestG[n] != 0){
		g=BestG[n]; 
		if(BestScoreG[g] < BestScore[n]){
		   if(names[g]) free(names[g]);
		   names[g]=AllocString(NameCMSA(SeedCMA[n]));
		   std::cerr << n << "-" << g << "." << names[g] << std::endl;
		}
	    } TotalNilCMSA(SeedCMA[n]);
	} free(BestScore); free(BestG); free(BestScoreG);
	return names;
}

void    omc_typ::PrintRTF(BooLean updateCSQ,BooLean SaveChnFiles,Int4 KeyStart,Int4 KeyEnd,char *KeySetName)
{
        double **value=0,**tmp; 
	NEWP(value,Hpt->NumSets()+3, double); NEWP(tmp,Hpt->NumSets()+3, double);
        Int4    i,j,n,d,depth=0;
#if 1
        for(n=2; n < Hpt->NumSets(); n++){
            depth=Hpt->NodeDepth(n);
            assert(depth <= 10);
            if(Hpt->TypeOfSet(n) == '!'){        // this is a leaf node.
                 value[n]=this->FindCrossConserved(n,0);
                 for(d=2; d < depth; d++){
                    if(tmp[d]){
		       if(value[n]==0) NEW(value[n], mcs->RtnLengthMainCMSA()+3,double);
                       for(i=1; i <= mcs->RtnLengthMainCMSA(); i++){
                          if(tmp[d][i] > value[n][i]){
                               value[n][i]=tmp[d][i];
                               if(value[n][0] < tmp[d][i]) value[n][0] = tmp[d][i];
                          }
                       }
                    }
                 }
            } else { if(tmp[depth]) free(tmp[depth]); tmp[depth]=this->FindCrossConserved(n,0); }
        }
        for(n=2; n < Hpt->NumSets(); n++) if(tmp[n]) free(tmp[n]); free(tmp);
#endif
	if(use_usr_sma) mcs->SetUserDisplaySet();
        mcs->PutRTF(updateCSQ,SaveChnFiles,KeyStart,KeyEnd,KeySetName,value);
        free(value); // value[x] will be freed by mcs->PutRTF();
}


