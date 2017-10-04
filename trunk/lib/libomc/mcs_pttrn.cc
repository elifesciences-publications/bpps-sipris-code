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

#include "mcs_typ.h"
#include "blosum62.h"
#include <iostream>

#if 1
Int4	mcs_typ::LoadUpBestColumns(Int4 NumToKeep)
// NumToKeep==0 then set NumToKeep = current number...
{
	Int4	i,j,n;
	char tmp_str[3003],str0[100];
        tmp_str[0]=0;

	// set che_typ max_num_columns to Length.
	this->UpdateCSQ(stderr);
	Int4 Length=LenSeq(che[1]->KeyE());
	fprintf(stderr,"Loading up with %d best columns...\n",NumToKeep);

	// Remove all columns.
	
	// load up with the best.
	for(n=2; n<= Hpt->NumBPPS(); n++){
	   for(j=1; j <= Length; j++){ 
	      // che[n]->RemoveColumn(j); 
	      sst_typ bsst=che[n]->ForceBestColumn(j,sst[n][j]);
	      // fprintf(stderr," pattern added to group %d column %d \n",n,j);
#if 0
                  i++;
                  if(i > 1) strncat(tmp_str,",",200);
                  PutSST(str0,bsst,AB);  // copies pattern to str.
                  strncat(tmp_str,str0,22);
                  sprintf(str0,"%d",j);
                  strncat(tmp_str,str0,22);
                  assert(tmp_str[3000]==0);  // don't over-extend string
#endif
	   }
	}
	// remove all but the target number whether or not they are positive.
	fprintf(stderr,"removing all but the target number whether or not they are positive\n");
	for(i=2; i <= Hpt->NumBPPS(); i++) this->RmWorstColumn(i,NumToKeep);
	fprintf(stderr,"Done\n");
}
#endif

Int4	mcs_typ::LoadUpColumns( )
// Load up columns to reinitialize FD-table. (afn: 3-12-2013).
{
	Int4 mc,n,s,len;
	for(n=2; n<= Hpt->NumBPPS(); n++){
	    assert(NumDisplayCMA == Hpt->NumBPPS());
	    mc=che[n]->RtnMaxNumColumns();
	    if(che[n]->NumColumns( ) >= mc) continue;
            len=LengthCMSA(1,DisplayCMA[n]);
	    che[n]->SetMaxNumColumns(len);
	    e_type keyE=che[n]->KeyE( ); 
            for(s=1; s <=len; s++){
	       if(che[n]->IsPatternPos(s)) continue;
               unsigned char r = ResSeq(s,keyE);
	       sst_typ sstJ=SsetLet(r);
	       assert(che[n]->AddColumn(s,sstJ));
	    } che[n]->SetMaxNumColumns(mc);
	    mc=len/2 + 1;
fprintf(stderr," ************ %d. columns = %d ************\n",n,mc);
// exit(1);
	    while(che[n]->NumColumns( ) > mc){
		if(che[n]->RemoveWorstColumn( )==FALSE){
		    fprintf(stderr,"%d: %d columns; %d maxcols.\n",n,che[n]->NumColumns( ), mc);
		    print_error("mcs_typ::LoadUpColumns( ): fatal - che[n]->RemoveWorstColumn( ) failed!");
		    print_error("Failed to find a significant hierarchy given the input constraints and seed.");
		}
	    }
	} 
}

void	mcs_typ::UpdateCSQ(FILE *fp)
// Only relevant to do this prior to column sampling; otherwise does not change the model...
{
	Int4	n,j,ri,ris,Length;
	e_type  keyE=0;

	if(!Evolve) return;
	if(fp) fprintf(fp,"=========== Updating evolving consensus sequence ===========\n");
	UnLabelAllSeqs( );	// this only happens on first call to this function.
	for(n=1; n<= Hpt->NumBPPS(); n++){
	    // if(n==2) PutSeq(stderr,che[n]->KeyE( ),AB);
	    che[n]->UpdateConsensus(); // if(n==2) PutSeq(stderr,che[n]->KeyE( ),AB);
	    keyE=che[n]->KeyE( ); 
	    bpps_typ *pps=che[n]->BPPS();  
	    if(pps->UpdateRho(sst[n],keyE)) DidRestoreBest=FALSE;   // now different from the best.
	} this->CheckPttrnCsqMatch("Error in UpdateCSQ()");
}

BooLean	mcs_typ::IsConflict(Int4 k, Int4 better, Int4 worse)
// Is there a conflict between the patterns at position k in categories n and m?
// if(che[n]->PttrnPos(k)){ insrtHeap(n,(keytyp)-(LPR_Wt[n]*SubLPR[n][k]),dH); } 
{
	Int4	n=better,m=worse;
	assert(n != m);			// This should not happen.
	sst_typ qst_n=che[n]->BPPS()->RtnSST(k); assert(qst_n);
	sst_typ qst_m=che[m]->BPPS()->RtnSST(k); assert(qst_m);
        char relFG=RelateFGs[n][m],relBG=RelateBGs[n][m];
	sst_typ	ist=IntersectSset(qst_n,qst_m);

	if(ist == 0) return FALSE; // patterns are distinct --> okay (no conflict).
	if(relBG == '0') return FALSE; // backgrounds are distinct --> okay.
	if(qst_n == qst_m){	// identical patterns...
		if(relFG != '0') return TRUE; // same pattern disallowed if FGs intersect.
		// if(relBG != '0') return TRUE; // 
		else if(relBG != '=') return FALSE; // Allow BGs to intersect...
		// else if(relBG != '0') return TRUE; // 
		else return FALSE; 
	} else if(SubSset(qst_m,qst_n)){

	} else if(SubSset(qst_n,qst_m)){

	} else {	// patterns overlap 
		assert(ist != 0);
	} return FALSE;
}

BooLean	mcs_typ::IsConflict0(Int4 k, Int4 better, Int4 worse)
// Is there a conflict between the patterns at position k in categories n and m?
// if(che[n]->PttrnPos(k)){ insrtHeap(n,(keytyp)-(LPR_Wt[n]*SubLPR[n][k]),dH); } 
{
	Int4 n=better,m=worse;

	assert(n != m);			// This should not happen.
	sst_typ qst_n=che[n]->BPPS()->RtnSST(k); assert(qst_n);
	sst_typ qst_m=che[m]->BPPS()->RtnSST(k); assert(qst_m);
        char relFG=RelateFGs[n][m],relBG=RelateBGs[n][m];
	// fprintf(stderr,"entering IsConflict(%d vs %d) col %d: relFG=%c; relBG=%c\n",n,m,k,relFG,relBG);
	if(qst_n == qst_m){	// identical patterns...
		// fprintf(stderr,"qst_n == qst_m = "); PutSST(stderr,qst_n,AB); fprintf(stderr,".\n");
	        if(StrictIndepend){
		    if(relFG == '0') return FALSE;	// foreground partitions disjoint=OK.
		    else return TRUE;
		}
		// if(relFG != '0' && relBG != '0')
		if(relFG != '0' || relBG != '0')
		{
 		   // fprintf(stderr,"Incompatible patterns (%d vs %d) %d:\n",n,m,m);
		   return TRUE;
		} else return FALSE; 
	} else if(relFG == '0'){ return FALSE; // if FGs distinct then pattern positions irrelevant.
#if 1	// Fix problem with subset constraints...
	} else if(relFG == '<'){   // FG_n is a proper subset of FG_m
		if(SubSset(qst_m,qst_n)) return TRUE; // qst_m < qst_n --> a conflict...
	} else if(relFG == '>'){ // FG_n is a proper superset of FG_m
		if(SuperSset(qst_m,qst_n)) return TRUE; // qst_m > qst_n --> a conflict...
#endif
	} else {	// i.e., qst_n != qst_m && relFG != '0'.
	  // fprintf(stderr,"qst_n != qst_m && relFG != '0'.\n");
	  sst_typ ist=IntersectSset(qst_n,qst_m);
	  if(ist == 0) return FALSE;	// patterns distinct; no problem.
	  sst_typ ust=UnionSset(qst_n,qst_m);
	  switch (relFG){		// based on relationship to 
	   case '0': return FALSE;		// FGs don't overlap; do nothing; okay.
	   case '=': 			// identical FG sets.
	     if(StrictIndepend) return TRUE;
	     if(relBG != '0') return TRUE; else return FALSE;
	   case '+': 			// FG sets overlap but not identical or subsets.
	        if(StrictIndepend) return TRUE;
		if(relBG != '0'){	
		  if(ist != 0) return TRUE; else return FALSE;
		} else return FALSE;
	   case '<':			// FG_n is a proper subset of FG_m
	    if(qst_n == ist){		  // qst_n is a proper subset of qst_m.
		assert(qst_n != qst_m);
		return FALSE; 		 // do nothing... okay.
	    } else {			  // patterns overlap but !(qst_n < qst_m).
	       if(StrictIndepend) return TRUE;
	       if(relBG != '0') return TRUE;	// BGs overlap
	       else return FALSE;	// else do nothing as BGs are distinct...
	    } break;
	   case '>':				// FG_n is a proper superset of FG_m
	    if(qst_m == ist){		  	// qst_m is a proper superset of qst_n	
		assert(qst_n != qst_m);
		return FALSE; 			// do nothing... okay.
	    } else {
	       if(StrictIndepend) return TRUE;
	       if(relBG != '0') return TRUE;		// BGs overlap
	       else return FALSE;		// else do nothing as BGs are distinct...
	    } break;
	   default: print_error("mcs_typ::IsConflict(). This should not happen (1).\n");
	    break;
	  }
	} return FALSE;
}

sst_typ	*mcs_typ::PruneSST(Int4	n, Int4 k, sst_typ *isst)
// remove patterns that won't work...to speed up sampling.
{
	sst_typ *osst;
	Int4	i,j,num=0,r;
	UInt4	totFG,totBG,*size,s;
	double	fg,bg,d,*FG,*BG,*best;
	BooLean	okay;

	assert(n > 0 && n <= Hpt->NumBPPS());
        for(num=0,i=1; isst[i]; i++) num++;
	NEW(osst,num + 3, sst_typ);
#if 1
	NEW(FG,num + 3, double); NEW(BG,num + 3, double); NEW(size,num+3,UInt4);
	NEW(best,30,double); // best ratio by size of pattern.
	for(i=1; isst[i]; i++){
	    size[i]=CardSST(isst[i],AB);
	    for(totFG=totBG=0,r=1; r <= nAlpha(AB); r++){
		totFG += che[n]->GetResWtFG(k,r); totBG += che[n]->GetResWtBG(k,r);
		if(MemSset(r,isst[i])){ 
		   FG[i] += (double) che[n]->GetResWtFG(k,r);
		   BG[i] += (double) che[n]->GetResWtBG(k,r);
		}
	    } FG[i] = FG[i]/(double) totFG; BG[i] = BG[i]/(double) totFG;
	}
	for(s=1; s <= 20; s++) best[s]=-9999999999999999.9;
	for(i=1; isst[i]; i++){
	    d=FG[i]/BG[i]; 
	    if(best[size[i]] < d){ best[size[i]]=d;  }
	}
	for(j=0,i=1; isst[i]; i++){
	    d=10*(FG[i]/BG[i]);
	    if(d >= best[size[i]]){ j++; osst[j]=isst[i];  }
	}
	free(FG); free(BG); free(size); free(best); 
#else
	for(j=0,i=1; isst[i]; i++){
	    for(totFG=totBG=0,r=1; r <= nAlpha(AB); r++){
		totFG += che[n]->GetResWtFG(k,r); totBG += che[n]->GetResWtBG(k,r);
	    }
	    for(okay=TRUE,r=1; r <= nAlpha(AB); r++){
		if(!MemSset(r,isst[i])) continue;
		fg=(double) che[n]->GetResWtFG(k,r)/(double) totFG;
		bg=(double) che[n]->GetResWtBG(k,r)/(double) totBG;
		d= fg/bg;
		if(d < 0.5){ okay=FALSE; break; }
	    } if(okay){ j++; osst[j]=isst[i]; }
	}
#endif
	if(j == 0){ free(osst); return 0; }
	else return osst;
}

#define DebugMissingPatterns 0	// can only use with a specific NAT analysis.

double	mcs_typ::SampleColumns(BooLean UseNegCol,char mode)
#if 0	//*****************************************************************
// Operations: remove or add columns; move sequences up or down.
set options:
   A & B != 0 (Intersecting but 	-->	
   	A < B	  (subset)	-->	only subpatterns allowed.
	A !< B	  (not subset)	-->	no common pattern positions allowed.
   A & B = 0 (disjoint)		-->	All common pattern positions allowed...
 					(As long as transitive non-existent)
	(A | C) & (B | D) != 0 	-->	C and D intersect. 
	Implies that FG BG sampling is coupled between sets.
#endif	//*****************************************************************
{
    Int4	i,j,k,K,m,n,x;
    double	temp,D,d;
    static Int4 toggle=0;

// cerr << "Debug 1\n";
    //==============================================================================
    //============== 0. compute LPR and update consensus sequences. ================
    //==============================================================================
    D=CalcTotalLPR(efp); // computes SubLPR[n][k] for all; Do prior to UpdateCSQ()!! see below.
    CheckPttrnCsqMatch("debug 4");
    if(Evolve) UpdateCSQ(); 	// This may cause the Csq and sst patterns to be out of sync!!!
    				// Don't want CalcTotalLPR() to store as best if out of sync!!!
    // CheckPttrnCsqMatch("debug 5");  // debug... ZZZZZZZZZZZZZZZ

    //==============================================================================
    //==== 1. Sample in all columns that most positively contribute to the LPR. ====
    //==============================================================================
    Int4 Length=LenSeq(che[1]->KeyE());
    temp = temperature; 
    if(temp > 300) temp=300;
// for(n=1; n<= Hpt->NumBPPS(); n++) che[n]->PutSubLPRs(stderr);  
    for(j=1; j <=1; j++){  // do this twice to converge on stable alpha...
     if(efp) fprintf(efp,"loading up columns....\n");
     for(n=1; n<= Hpt->NumBPPS(); n++){
	if(IsFailedBPPS[n]) continue;
// che[n]->PutSubLPRs(stderr);
	bpps_typ *bpps=che[n]->BPPS();
	assert(LenSeq(che[n]->KeyE()) == Length);  // make sure that all same # columns!
	m=che[n]->RtnMinNumColumns();
	che[n]->SetMinNumColumns(0);		// set the minimum to zero...
	che[n]->SetMaxNumColumns(Length);	// redundant; also set above...
        for(k=1; k <= Length; k++){
	    che[n]->RemoveColumn(k);	// want to recalculate LPR here.
#if 0	// prune the pattern set to speed up sampling.  THIS DOESN'T SEEM TO HELP!!!
	    sst_typ *osst=PruneSST(n,k,sst[n][k]);
	    if(osst==0) continue;
	    BooLean okay=che[n]->SamplePattern(-99999999999999.9,k,osst,temp); free(osst);
	    if(okay)
#else
#if DebugMissingPatterns
	    // WARNING: THIS ASSUMES A TREE!!!
	    if(n==9 && (k==57 || k==18 || k==36)){
		che[n]->SpeakUp( );
	    }
#endif
	    if(che[n]->SamplePattern(-99999999999999.9,k,sst[n][k],temp))
#endif
	    {
#if DebugMissingPatterns
		if(n==9 && (k==57 || k==18 || k==36)){
		   che[n]->BeQuiet( );
		   sst_typ xsst= bpps->RtnSST(k);
		   fprintf(stderr,"%d. ",n);
		   PutSST(stderr,xsst,AB);
		    double  *xSubLpr=che[n]->SubMap();
		   fprintf(stderr,"%d (%.2f)\n",k,xSubLpr[k]);
		}
#endif
		double  *sub_lpr=che[n]->SubMap();
		// if(!UseNegCol && sub_lpr[k] < -10.0)	// changed for testing...
		// if(sub_lpr[k] < -10.0)
		if(sub_lpr[k] <= 0.0)
		{
		   if(efp) fprintf(efp,"Added negative col %d from analysis %d (%.2f;%.3f).\n",
			k,n,sub_lpr[k],bpps->Alpha( ));
		   // if(j == 1) che[n]->RemoveColumn(k);
		   // alpha changes when add or remove columns and this changes LPR at other sites!
		   // null_LPR depends only on priors --> priors effect zero cutoff.
		   // a negative column can increase the LPR by changing alpha...
		}
	    } else assert(!che[n]->PttrnPos(k)); // if failed then there should be no pattern.
	} che[n]->SetMinNumColumns(m);
// che[n]->PutSubLPRs(stderr); exit(1);
	// { che[n]->PutSubLPRs(stderr);  exit(1); }
	// che[n]->SetMaxNumColumns(MaxNumCol[n]);
// che[n]->PutSubLPRs(stderr);
     }
    }
    CheckPttrnCsqMatch("debug 5");	// ZZZZZZZZZZZZ

   //==============================================================================
   // 2. ****************** Create a pattern column max heap... ****************
   //==============================================================================
   Int4 *BestToWorstCol,*BestToWorstCat,TotalBestToWorst;
   {	//****************** begin of pch scope. **********************
     double key;
     pch_typ pch = pch_typ(Hpt->NumBPPS(),Length);
     // 3. Compute the subLPR for all categories and positions...
     CalcTotalLPR(efp,FALSE); // computes SubLPR[n][k] for all.
     if(efp) fprintf(efp,"sampling columns....\n");
     // 4. Insert all positive columns into pch.
     for(j=0,n=1; n<= Hpt->NumBPPS(); n++){
	if(IsFailedBPPS[n]) continue;
	for(k=1; k <= Length; k++){
	   if(che[n]->PttrnPos(k)){
#if 0	 // Eventually sample a key proportional to the pattern probabilities.
	// 1. need to normalize LPRs; sample one and put on the heap...
	//    ... then renormalize and resample, etc. until all have been sampled.
	      key = SubLPR[n][k];
#else	// currently just taking the actual prob. as the key.
	      key = SubLPR[n][k];
#endif
	      // if(UseNegCol || key > -10.0){ j++; pch.Insert(j,key,n,k); } 
	      // if(key > -10.0) { j++; pch.Insert(j,key,n,k); }
	      if(key > 0.0) { j++; pch.Insert(j,key,n,k); 
#if DebugMissingPatterns
		if(n==9 && (k==57 || k==18 || k==36)){ // for NAT analysis...
		   bpps_typ *bpps=che[n]->BPPS();
		   sst_typ xsst= bpps->RtnSST(k);
		   fprintf(stderr,"insert onto heap: %d. ",n);
		   PutSST(stderr,xsst,AB);
		   fprintf(stderr,"%d (%.2f)\n",k,key);
		}
#endif
	      } else {
		if(efp) fprintf(efp,"Removing negative col %d from analysis %d (%f).\n",k,n,key);
		// che[n]->RemoveColumn(k); // recalculates LPR!!! changes SubLPR[n][k] !!!
		che[n]->RmColButKeepLPR(k); // removes column without recalculating LPR.
	      }
	   } else assert(SubLPR[n][k] == 0.0);
	}
     } TotalBestToWorst=pch.NumItems(); assert(TotalBestToWorst == j);

     //==============================================================================
     //======== 3. Order all pattern-positions by contribution to total LPR; ========
     //==============================================================================
     // fprintf(stderr,"Ordering pattern-positions by contribution to total LPR.\n");
     NEW(BestToWorstCol,TotalBestToWorst+3,Int4);
     NEW(BestToWorstCat,TotalBestToWorst+3,Int4);
     for(j=0; !pch.Empty(); ){
	i=pch.DeleteMax(key,n,k);	// pch.DeleteMax(&key,&n,&k);
	j++; BestToWorstCol[j]=k; BestToWorstCat[j]=n;
	// fprintf(stderr,"%d. key = %.3f; k = %d; n = %d\n",j,key,k,n);
     } assert(TotalBestToWorst == j);
   } //********************** end of pch scope. **************************

    //==============================================================================
    //======== 4. "Add in" columns starting from the best to the worst.  ===========
    //==============================================================================
    // fprintf(stderr,"Adding columns starting from the best to the worst (%d).\n",TotalBestToWorst); 
#if 1  // alternative method that disallows identical patterns at the same position...
    //------- 4a. Remove column pattern-positions that are identical for distinct BPPS. 
    // keep track of identical positions and patterns
    ds_type djs=DSets(TotalBestToWorst+3);
    Int4 id,jd,djs_i,djs_j;
    set_typ SetN=MakeSet(Hpt->NumBPPS()+3);
    set_typ SetI=MakeSet(TotalBestToWorst+3); ClearSet(SetI);
    Int4 *NtoI;
    for(k=1; k <= Length; k++){

      //=== 4a. Find the subgroups with a pattern at position k. ===
      ClearSet(SetN); NEW(NtoI,Hpt->NumBPPS() +3,Int4);
      for(i=1; i <= TotalBestToWorst; i++){
	K=BestToWorstCol[i]; assert(K <= Length);
	if(K != k) continue;
	n=BestToWorstCat[i]; assert(n <= Hpt->NumBPPS());
	AddSet(n,SetN); NtoI[n]=i; 
	assert(!MemberSet(i,SetI)); AddSet(i,SetI);
#if DebugMissingPatterns
		if(n==9 && (k==57 || k==18 || k==36)){ // for NAT analysis...
		   bpps_typ *bpps=che[n]->BPPS();
		   sst_typ xsst= bpps->RtnSST(k);
		   fprintf(stderr,"removing best from heap: %d. ",n);
		   PutSST(stderr,xsst,AB);
		    double  *xSubLpr=che[n]->SubMap();
		   fprintf(stderr,"%d (%.2f)\n",k,xSubLpr[k]);
		}
#endif
      } // fprintf(stderr,"k = %d\n",k);

      //=== 4b. Put those subgroups with the same pattern into the same disjoint set. ===
      for(n=1; n < Hpt->NumBPPS(); n++){
	if(MemberSet(n,SetN)) i=NtoI[n]; else continue;
        // fprintf(stderr,"n = %d; i = %d; ",n,i);
	djs_i = findDSets(i,djs);
	sst_typ qst_n=che[n]->BPPS()->RtnSST(k); assert(qst_n);
	// PutSST(stderr,qst_n,AB); fprintf(stderr,"%d.\n",k); 
        for(m=n+1; m <= Hpt->NumBPPS(); m++){
	  if(MemberSet(m,SetN)) j=NtoI[m]; else continue;
	  djs_j = findDSets(j,djs);
          // fprintf(stderr,"m = %d; j = %d\n",m,j);
	  sst_typ qst_m=che[m]->BPPS()->RtnSST(k); assert(qst_m);
	  // PutSST(stderr,qst_m,AB); fprintf(stderr,"%d.\n",k); 
	  if(qst_n == qst_m){     // identical patterns at the same position
		if(djs_i != djs_j) djs_i = linkDSets(djs_i,djs_j,djs);
#if DebugMissingPatterns
		if((k==57 || k==18 || k==36) && (n==9 || m==9)){ // for NAT analysis...
		   bpps_typ *bpps=che[n]->BPPS();
		   fprintf(stderr,"same patterns: %d & %d. ",n,m);
		   PutSST(stderr,qst_n,AB); fprintf(stderr,"%d\n",k);
		}
#endif
	  }
        } 
      } free(NtoI);
    } NilSet(SetN);  ClearSet(SetI);

    //=== 4c. Convert disjoint set to a linked list matching subgroups. ===
    Int4 *djsID; NEW(djsID,TotalBestToWorst + 4, Int4);
    Int4 *djsNxt,last; NEW(djsNxt,TotalBestToWorst + 4, Int4);
//         0  1  2  3  4  5  6  7  8  9		== djsNxt[ ]
//        |0| 3| 5| 7| 6| 8| 0| 9| 0| 0|          (0 = end of list.)
//            1---->3---------->7---->9.
//               2------->5------->8. 
//                     4---->6.
    for(i=1; i <= TotalBestToWorst; i++){ djs_i=findDSets(i,djs); djsID[i]=djs_i; }
    for(i=1; i <= TotalBestToWorst; i++){
	id=djsID[i]; assert(id != 0); 
	if(djsNxt[i] != 0) continue;	// this one already filled in.
	for(last=i,j=i+1; j <= TotalBestToWorst; j++){
	   if(djsID[j]==id){ djsNxt[last]=j; last=j; }
	}
    } NilDSets(djs); free(djsID);
#if 0	// DEBUG...
    for(i=1; i <= TotalBestToWorst; i++){
	if(MemberSet(i,SetI)) continue; else AddSet(i,SetI);
	n=BestToWorstCat[i]; k=BestToWorstCol[i];  
	if(djsNxt[i]==0) continue;
	fprintf(stderr,"%d(%d:%d)",i,k,n);
	for(j=djsNxt[i]; j != 0; j=djsNxt[j]){
	   m=BestToWorstCat[j]; K=BestToWorstCol[j];  
	   assert(k==K);
	   fprintf(stderr,"->%d(%d:%d)",j,K,m); AddSet(j,SetI);
	} fprintf(stderr,".\n");
	// if(i > 200) break;
    }
#endif
    NilSet(SetI);
#endif

   //==============================================================================
   //===================== 5. Remove remaining columns. ==========================
   //==============================================================================
    Int4 *TotalAddedCol,*TotalAddedCat;
    NEW(TotalAddedCol,Hpt->NumBPPS()+3,Int4); NEW(TotalAddedCat,Length+3,Int4);
    for(j=1; j <= TotalBestToWorst; j++){
	n=BestToWorstCat[j]; 
	if(n == 0) continue; 			// this column has already been removed.
        //=== 5a. if the maximum reached for subgroup n, then remove the remaining columns. ===
	if(TotalAddedCol[n] == MaxNumCol[n]){	
	  for(i=j; i <= TotalBestToWorst; i++){	
	    if(BestToWorstCat[i] == n){			// is this category n?
		k=BestToWorstCol[i]; 
		che[n]->RmColButKeepLPR(k); // removes column without recalculating LPR.
		BestToWorstCat[i] = 0; BestToWorstCol[i]=0;
	    }
	  }
	} else {					// still more left to add.
        //=== 5b. Otherwise add another column. ===
	  assert(TotalAddedCol[n] < MaxNumCol[n]);
	  k=BestToWorstCol[j]; assert(k != 0);
	  if(TotalAddedCat[k] > 0){	// if already added a column here then check compatibility.
	     for(i=j-1; i > 0; i--){		// find the previous categories at this position.
		if(BestToWorstCol[i] == k){		// is this the same column?
		  m=BestToWorstCat[i];			// m is a previous (better) category
		  assert(m != 0 && m != n);		// this should not happen...
		  // fprintf(stderr,"cat %d vs %d at col %d.\n",n,m,k);
		  if(IsConflict(k,m,n)){		// is there a conflict between m & n at k?
#if DebugMissingPatterns
		    if((k==57 || k==18 || k==36) && (n==9)){ // for NAT analysis...
			bpps_typ *bpps=che[n]->BPPS();
			sst_typ xsst=bpps->RtnSST(k);
			fprintf(stderr,"removing pattern: %d. ",n);
			PutSST(stderr,xsst,AB); fprintf(stderr,"%d\n",k);
		    }
#endif
		        // fprintf(stderr,"conflict found!\n\n");
			che[n]->RmColButKeepLPR(k); // removes column without recalculating LPR.
			BestToWorstCat[j] = 0; BestToWorstCol[j]=0;
			break;	// no need to check further...one conflict is enough.
		  } // else fprintf(stderr,"conflict not found.\n\n");
		}					// else skip over these...
	     }
	     if(i == 0){ TotalAddedCat[k]++; TotalAddedCol[n]++; }  // else add column k for n.
	  } else { TotalAddedCat[k]++; TotalAddedCol[n]++; }	// add first pattern at k...
#if 0	// don't want to remove all because some identical patterns may not have a conflict!!!
	  if(BestToWorstCat[j] != 0){ // then need to remove equivalent pattern/positions. 
		if(efp) fprintf(efp,"%d(%d:%d)",j,k,n);
		for(x=djsNxt[j]; x != 0; x=djsNxt[x]){
		     m=BestToWorstCat[x]; K=BestToWorstCol[x];  
		     if(efp) fprintf(efp,"->%d(%d:%d)",x,K,m);
		     if(m == 0) continue;	// item was previously deleted.
		     assert(k==K);  
#if DebugMissingPatterns
		    if((k==57 || k==18 || k==36) && (m==9)){ // for NAT analysis...
			bpps_typ *bpps=che[m]->BPPS();
			sst_typ xsst=bpps->RtnSST(k);
			fprintf(stderr,"removing pattern: %d. ",m);
			PutSST(stderr,xsst,AB); fprintf(stderr,"%d\n",k);
		    }
#endif
		     che[m]->RmColButKeepLPR(k); // removes column without recalculating LPR.
		     BestToWorstCat[x] = 0; BestToWorstCol[x]=0;
		} if(efp) fprintf(efp,".\n");
	  }
#endif
	}
    } free(TotalAddedCol); free(TotalAddedCat); free(djsNxt);

   //==============================================================================
   //============================= 5. Free memory. ================================
   //==============================================================================
   free(BestToWorstCol); free(BestToWorstCat); 
   d=CalcTotalLPR(efp); // only recalculate LPR at very end.
   toggle++; if(toggle > PrintToggle) toggle=0;
   if(toggle==0) fprintf(stderr,"Sampling columns...LPR: %.2f to %.2f\n",D,d);
   // CheckPttrnCsqMatch("debug 6");	// ZZZZZZZZZZZZ
   return d; // only recalculate LPR at very end.
}

#if 0	// old code (10-30-2013) replaced by above...
double	mcs_typ::SampleColumns(BooLean UseNegCol,char mode)
#if 0	//*****************************************************************
// Operations: remove or add columns; move sequences up or down.
set options:
   A & B != 0 (Intersecting but 	-->	
   	A < B	  (subset)	-->	only subpatterns allowed.
	A !< B	  (not subset)	-->	no common pattern positions allowed.
   A & B = 0 (disjoint)		-->	All common pattern positions allowed...
 					(As long as transitive non-existent)
	(A | C) & (B | D) != 0 	-->	C and D intersect. 
	Implies that FG BG sampling is coupled between sets.
#endif	//*****************************************************************
{
    Int4	i,j,k,m,n;
    double	temp,D,d;
    static Int4 toggle=0;

    D=CalcTotalLPR(efp); // computes SubLPR[n][k] for all; Do prior to UpdateCSQ()!! see below.
    CheckPttrnCsqMatch("debug 4");
    if(Evolve) UpdateCSQ(); 	// This may cause the Csq and sst patterns to be out of sync!!!
    				// Don't want CalcTotalLPR() to store as best if out of sync!!!
    // CheckPttrnCsqMatch("debug 5");  // debug...
    //==== 1. Sample in all columns that most positively contribute to the LPR. ====
#if 0	// This is due to columns with patterns LACKING matches to foreground sequences.
    if(D < -999999.0){	// goes away if minimum # columns not imposed or if evolve csq/pattern.
	this->PutHyperPartition(stderr);
	for(n=1; n<= Hpt->NumBPPS(); n++){
		if(SubLPR[n][0] < -99999.0) che[n]->PutSubLPRs(stderr);
	}
    }
#endif
    Int4 Length=LenSeq(che[1]->KeyE());
    temp = temperature; 
    if(temp > 300) temp=300;
// for(n=1; n<= Hpt->NumBPPS(); n++) che[n]->PutSubLPRs(stderr);  
    // for(j=1; j <=2; j++){  // do this twice to converge on stable alpha...
    for(j=1; j <=1; j++){  // do this twice to converge on stable alpha...
     if(efp) fprintf(efp,"loading up columns....\n");
     for(n=1; n<= Hpt->NumBPPS(); n++){
	if(IsFailedBPPS[n]) continue;
// che[n]->PutSubLPRs(stderr);
	bpps_typ *bpps=che[n]->BPPS();
	assert(LenSeq(che[n]->KeyE()) == Length);  // make sure that all same # columns!
	m=che[n]->RtnMinNumColumns();
	che[n]->SetMinNumColumns(0);
	che[n]->SetMaxNumColumns(Length);	// redundant; also set above...
        for(k=1; k <= Length; k++){
	    che[n]->RemoveColumn(k);	// want to recalculate LPR here.
	    if(che[n]->SamplePattern(-99999999999999.9,k,sst[n][k],temp)){
		double  *sub_lpr=che[n]->SubMap();
		// if(!UseNegCol && sub_lpr[k] < -10.0)	// changed for testing...
		// if(sub_lpr[k] < -10.0)
		if(sub_lpr[k] <= 0.0)
		{
		   if(efp) fprintf(efp,"Added negative col %d from analysis %d (%.2f;%.3f).\n",
			k,n,sub_lpr[k],bpps->Alpha( ));
		   // if(j == 1) che[n]->RemoveColumn(k);
		   // alpha changes when add or remove columns and this changes LPR at other sites!
		   // null_LPR depends only on priors --> priors effect zero cutoff.
		   // a negative column can increase the LPR by changing alpha...
		}
	    } else assert(!che[n]->PttrnPos(k)); // if failed then there should be no pattern.
	} che[n]->SetMinNumColumns(m);
// che[n]->PutSubLPRs(stderr); exit(1);
	// { che[n]->PutSubLPRs(stderr);  exit(1); }
	// che[n]->SetMaxNumColumns(MaxNumCol[n]);
// che[n]->PutSubLPRs(stderr);
     }
    }
   // 2. ****************** Create a pattern column max heap... ****************
   Int4 *BestToWorstCol,*BestToWorstCat,TotalBestToWorst;
   {	//****************** begin of pch scope. **********************
     double key;
     pch_typ pch = pch_typ(Hpt->NumBPPS(),Length);
     // 3. Compute the subLPR for all categories and positions...
     CalcTotalLPR(efp,FALSE); // computes SubLPR[n][k] for all.
     if(efp) fprintf(efp,"sampling columns....\n");
     // 4. Insert all positive columns into pch.
     for(j=0,n=1; n<= Hpt->NumBPPS(); n++){
	if(IsFailedBPPS[n]) continue;
	for(k=1; k <= Length; k++){
	   if(che[n]->PttrnPos(k)){
	      key = SubLPR[n][k];
	      // if(UseNegCol || key > -10.0){ j++; pch.Insert(j,key,n,k); } 
	      // if(key > -10.0) { j++; pch.Insert(j,key,n,k); }
	      if(key > 0.0) { j++; pch.Insert(j,key,n,k); }
	      else {
		if(efp) fprintf(efp,"Removing negative col %d from analysis %d (%f).\n",k,n,key);
		// che[n]->RemoveColumn(k); // recalculates LPR!!! changes SubLPR[n][k] !!!
		che[n]->RmColButKeepLPR(k); // removes column without recalculating LPR.
	      }
	   } else assert(SubLPR[n][k] == 0.0);
	}
     } TotalBestToWorst=pch.NumItems(); assert(TotalBestToWorst == j);

     // 5. Order all pattern-positions by contribution to total LPR;
     // fprintf(stderr,"Ordering pattern-positions by contribution to total LPR.\n");
     NEW(BestToWorstCol,TotalBestToWorst+3,Int4);
     NEW(BestToWorstCat,TotalBestToWorst+3,Int4);
     for(j=0; !pch.Empty(); ){
	i=pch.DeleteMax(key,n,k);	// pch.DeleteMax(&key,&n,&k);
	j++; BestToWorstCol[j]=k; BestToWorstCat[j]=n;
	// fprintf(stderr,"%d. key = %.3f; k = %d; n = %d\n",j,key,k,n);
     } assert(TotalBestToWorst == j);
   } //********************** end of pch scope. **************************

    // 6. "Add in" columns starting from the best to the worst. 
    // fprintf(stderr,"Adding columns starting from the best to the worst (%d).\n",TotalBestToWorst); 
    Int4 *TotalAddedCol,*TotalAddedCat;
    NEW(TotalAddedCol,Hpt->NumBPPS()+3,Int4); NEW(TotalAddedCat,Length+3,Int4);
    for(j=1; j <= TotalBestToWorst; j++){
	n=BestToWorstCat[j]; 
	if(n == 0) continue; 			// this column has already been removed.
	if(TotalAddedCol[n] == MaxNumCol[n]){	// if the maximum have been added, then remove the rest.
	  for(i=j; i <= TotalBestToWorst; i++){	// remove all remaining columns in this category.
	    if(BestToWorstCat[i] == n){			// is this the same category?
		k=BestToWorstCol[i]; 
		che[n]->RmColButKeepLPR(k); // removes column without recalculating LPR.
		// che[n]->RemoveColumn(k);
		BestToWorstCat[i] = 0; BestToWorstCol[i]=0;
	    }
	  }
	} else {					// still more left to add.
	  assert(TotalAddedCol[n] < MaxNumCol[n]);
	  k=BestToWorstCol[j];
	  assert(k != 0);
	  if(TotalAddedCat[k] > 0){	// if already added a column here then check compatibility.
	     for(i=j-1; i > 0; i--){		// find the previous categories at this position.
		if(BestToWorstCol[i] == k){		// is this the same column?
		  m=BestToWorstCat[i];			// m is a previous (better) category
		  assert(m != 0 && m != n);		// this should not happen...
		  // fprintf(stderr,"cat %d vs %d at col %d.\n",n,m,k);
		  if(IsConflict(k,m,n)){		// is there a conflict between m & n at k?
		        // fprintf(stderr,"conflict found!\n\n");
			// che[n]->RemoveColumn(k);	// then remove the worst pattern...
			che[n]->RmColButKeepLPR(k); // removes column without recalculating LPR.
			BestToWorstCat[j] = 0; BestToWorstCol[j]=0;
			break;	// no need to check further...one conflict is enough.
		  } // else fprintf(stderr,"conflict not found.\n\n");
		}					// else skip over these...
	     }
	     if(i == 0){ TotalAddedCat[k]++; TotalAddedCol[n]++; }  // else add column k for n.
	  } else { TotalAddedCat[k]++; TotalAddedCol[n]++; }	// add first pattern at k...
	}
    }
   // 7. Free memory.
   free(BestToWorstCol); free(BestToWorstCat); free(TotalAddedCol); free(TotalAddedCat);
   d=CalcTotalLPR(efp); // only recalculate LPR at very end.
   toggle++; if(toggle > PrintToggle) toggle=0;
   if(toggle==0) fprintf(stderr,"Sampling columns...LPR: %.2f to %.2f\n",D,d);
   return d; // only recalculate LPR at very end.
}
#endif


