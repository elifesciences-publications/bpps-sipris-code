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

#include "gmb_typ.h"

BooLean	gmb_typ::SampleClustersOfClusters(set_typ ClstSet,Int4 **cluster, Int4 NumClst, Int4 **offset,
				double Temp,double &MAP)
// sample a gapped alignment for a cluster of sequences.
{
     Int4	i,j,n,start,strt,trace_length,score,sq,ii,*tmppos,size;
     BooLean	*Modified=0,found=FALSE;
     Int4	os,x,z,N,*clstr,***oldpos,***newpos,card=CardSet(ClstSet); 
     cma_typ	CMA=SSX->RtnCMA(); 
     double	*newmap,*oldmap; MAP=SSX->GapMap(Temp); 
     e_type	csqE=0;
     gsq_typ *gsq,***ogsq; 
     NEW(clstr, card+3,Int4); 
     NEWPP(oldpos, card+3, Int4); 
     NEWP(ogsq,card+2,gsq_typ *); NEW(oldmap,card+2,double); NEW(newmap,card+2,double);
     // PutSet(stderr,ClstSet);
     //======= 1. Remove seq clusters from the alignment.  ==========
     for(x=0,n=1; n <= NumClst; n++){
	   size = cluster[n][0];
	   assert(size > 0);
	   // if(size <= 0) continue;	// Then why was this included!!?? free memory read.
	   if(MemberSet(n,ClstSet)){
		x++; NEW(ogsq[x],size+2,gsq_typ *); clstr[x]=n;
		oldpos[x]=SSX->RemoveFromAlign(cluster[n]); 
		for(j=0,ii=1; cluster[n][ii] != 0; ii++) j++;
		assert(j == size);
	   }
     } N=x;
     //========== 2. Compute original LLRs ('oldmap') for each cluster. ============
     for(x=1; x <= N; x++){ 
	n=clstr[x]; 
	for(ii=1; (sq=cluster[n][ii]) != 0; ii++) SSX->AddToAlign(sq,oldpos[x][ii]); 
	oldmap[x]=SSX->GapMap(Temp);
	SSX->RmFromAlign(cluster[n]);	// does not return positions...
     }
     // fprintf(stderr,"N=%d; NumClst=%d\n",N,NumClst);

     //========== 3. Sample a gapped alignment for each master sequence csqE. ===========
     // SSX->InitNDL(Temp);	// initialize JHL parameters? Not used...?
     BooLean success=FALSE;
     gss_typ *gss=gssCMSA(CMA);
     NEW(Modified,N+3,BooLean); // Modified=FALSE;
     Int4 LenHMM=TotalLenCMSA(CMA);
     for(x=1; x <= N; x++){	// N = Number of clusters...
	n=clstr[x]; 
	csqE=CsqQuery[n]; 
 	if(csqE){	// implies multiple sequences in cluster.
	 // fprintf(stderr,"%d(%d): size=%d\n",x,n,cluster[n][0]);
	 this->GetMultipleTraces(csqE,Temp); // OPS->Put(stderr); // create ops_typ *OPS.
	 // ========== 3b. Align sequences based on the master. ===========
	 for(ii=1; (sq=cluster[n][ii]) != 0; ii++) {
	   e_type sbjE=gss->TrueSeq(sq); 	
	   assert(n == clstr[x]); os=offset[n][ii]; 
	   char *operation,*op,*bst_op=0;
	   Int4 bst_strt=0,bst_scr=0,tmp_scr=0; 
	   for(Int4 J=1; J <= OPS->N; J++){	// find the best consensus repeat among alternatives.
		operation=OPS->operation[J]; start=OPS->start[J]; strt=start-os;

		if((strt + LenHMM -1) <= 5) continue;	// domain maps prior to sbjE!
		else if((LenSeq(sbjE) - start) < 5) continue; // domain maps after sbjE!

		op=ConvertOperation(operation,sbjE,start,os);
		tmp_scr=this->AlignLenOperation(op);
		if(strt < 1) z=1; else z=strt;
		// tmp_scr=this->ScoreSeqAlnSMatrix(op,sbjE,z);
		if(strt >= LenSeq(sbjE)){
		  // PutDebugClocl(strt,os,sq,op,sbjE);
		  // fprintf(stderr,"3 ===========> %d: op=%s\n",sq,op); fflush(stderr);
		  free(op); op=0; continue;
		} else if(strt < 1){ 
		    // PutDebugClocl(strt,os,sq,op,sbjE);
	            // fprintf(stderr,"%d: strt=%d; os=%d; lenHMM=%d.\n",sq,start-os,os,LenHMM);
		    strt=1;	// starts < 1 already modeled in the operation array ("EDdddmm...");
		    if(bst_scr < tmp_scr){
			if(bst_op) free(bst_op);
			bst_strt=strt; bst_scr=tmp_scr; bst_op=op; op=0; 
		    } else { 
		    	// fprintf(stderr,"O ===========> %d: op=%s\n",sq,operation); 
		    	// fprintf(stderr,"4 ===========> %d: op=%s\n",sq,op); 
			free(op); op=0; continue; 
		    }
		} else if(bst_scr < tmp_scr){
		    if(bst_op) free(bst_op);
		    bst_strt=strt; bst_scr=tmp_scr; bst_op=op; op=0; 
		} else { free(op); op=0; }
	   } op=bst_op; strt=bst_strt;
	   if(op == 0){	// nothing works for this sequence...!!!
	      // fprintf(stderr,"%d: strt=%d; os=%d; score=%d\n",sq,strt,os,bst_scr);
	      SSX->AddToAlign(sq,oldpos[x][ii]); 
	   } else { 
	      assert(op[0]=='E');
	      NEW(tmppos,nBlksCMSA(CMA)+2,Int4);  // for new sites.
	      gsq = new gsq_typ[1]; // NEEDS TO BE AN ARRAY FOR gss_typ!!!
	      trace_length=strlen(op);
	      gsq->initialize(gss->LeftFlank(),gss->RightFlank(),op,trace_length,strt,sbjE,tmppos);
	      if(gss->Identical(sq,*gsq)){ 
		 BooLean the_same=TRUE;
                 for(Int4 b=1; b<= nBlksCMSA(CMA); b++)
                     { if(oldpos[x][ii][b] != tmppos[b]){ the_same=FALSE; break; } }
                 if(the_same){ delete []gsq; }
                 else { Modified[x]=TRUE; ogsq[x][ii]=SwapGsqCMSA(sq,gsq,CMA); }
	      } else { Modified[x]=TRUE; ogsq[x][ii]=SwapGsqCMSA(sq,gsq,CMA); }
	      SSX->AddToAlign(sq,tmppos); free(tmppos); free(op); 
	   }
	 } if(OPS) delete OPS; OPS=0;	// end of for(ii=1; (sq=cluster[n][ii]) != 0; ...) loop
	} else {	// Single sequence in the cluster.
	    assert(cluster[n][0] == 1); sq=cluster[n][1];
	    e_type  sbjE = gss->TrueSeq(sq);
	    char    *operation=SSX->GapAlnTrace(sbjE,Temp,start,score);
	    trace_length=strlen(operation);
	    gsq = new gsq_typ[1]; NEW(tmppos,nBlksCMSA(CMA)+3,Int4);
            gsq->initialize(gss->LeftFlank(),gss->RightFlank(),operation,trace_length,start,sbjE,tmppos);
	    if(gss->Identical(sq,*gsq)){
		 BooLean the_same=TRUE;
                 for(Int4 b=1; b<= nBlksCMSA(CMA); b++)
                     { if(oldpos[x][1][b] != tmppos[b]){ the_same=FALSE; break; } }
                 if(the_same){ delete []gsq; }
                 else { Modified[x]=TRUE; ogsq[x][1]=SwapGsqCMSA(sq,gsq,CMA); }
	    } else { Modified[x]=TRUE; ogsq[x][1]=SwapGsqCMSA(sq,gsq,CMA); }
	    SSX->AddToAlign(sq,tmppos); free(tmppos); free(operation); 
	}
     }	// end of for(x=1; x <= N; x++)...loop
     NEWPP(newpos, card+3, Int4); 
     //========== 4. Sample one of each of the two alignments. =================
     for(x=1; x <= N; x++) newpos[x]=SSX->RemoveFromAlign(cluster[clstr[x]]); 
     //============= compute new LLRs ('newmap') for each cluster. ============
     for(x=1; x <= N; x++){ 
	n=clstr[x]; 
	for(ii=1; (sq=cluster[n][ii]) != 0; ii++) { SSX->AddToAlign(sq,newpos[x][ii]); }
	newmap[x]=SSX->GapMap(Temp);
	SSX->RmFromAlign(cluster[n]);	// does not return positions...
     }
     //=============== Add them all back in again. ===============
     for(x=1; x <= N; x++){
	for(ii=1; (sq=cluster[clstr[x]][ii]) != 0; ii++){
	   SSX->AddToAlign(sq,newpos[x][ii]); free(newpos[x][ii]);
	} free(newpos[x]);
     } free(newpos);

     //========== 5. Sample either old or new alignment for each cluster. =============
     for(x=1; x <= N; x++){	// N = Number of clusters...
	n=clstr[x]; 
	//if(found){ this->PutClusterAln("junk_debug",cluster[n],CMA); } 
	if(Modified[x] && this->SampleTransition(newmap[x],oldmap[x],Temp)){
	   for(ii=1; cluster[n][ii] != 0; ii++) {
		free(oldpos[x][ii]); 
		if(ogsq[x][ii]) delete []ogsq[x][ii]; 
	   } success=TRUE;
	} else {	// reject the new cmsa file.
	   SSX->RmFromAlign(cluster[n]);
	   for(ii=1; (sq=cluster[n][ii]) != 0; ii++) {
		if(ogsq[x][ii]){ ReplaceCMSA(sq,ogsq[x][ii],CMA); }
		SSX->AddToAlign(sq,oldpos[x][ii]); free(oldpos[x][ii]); 
	   }
	} free(ogsq[x]); free(oldpos[x]); 
     } free(oldmap); free(oldpos); free(ogsq); free(clstr); free(Modified); free(newmap);
     MAP=SSX->GapMap(Temp);
     return success;
}

double	gmb_typ::ClustersOfClustersSampling(char mode,double Temp,Int4 percent_ident)
// sample with gaps using simulated annealing...
// 	char	method='E'; // gapped E-value
{
   	Int4	time1=time(NULL);
	cma_typ CMA=SSX->RtnCMA(); assert(CMA != 0); 
	Int4	n,i,j,k,r,s,sq,t,N=NumSeqsCMSA(CMA);
	double	dd,best,map;
	FILE	*efp=0; // efp=stderr;
	dh_type dH;

	SetPseudoToMapCMSA(CMA); best=SSX->GapMap(Temp); SaveBestCMSA(CMA); // Save initial alignment.

	//============== 1a. cluster the sequences into closely related sets. ============
	Int4    Nset,**cluster,**offset,numSets;
	cluster = this->GetClusters(efp,offset,Nset);
	if(Nset < this->MinNumClstrs){
	   fprintf(stderr,"too few (%d < %d) clusters",Nset,this->MinNumClstrs);
           for(n=1; n <= Nset; n++){ free(cluster[n]);  if(offset[n]) free(offset[n]); }
           free(cluster); free(offset);
           return GambitSingles(mode,Temp);
	}
#if 1
	set_typ *CloseSq=RtnTargetSizeSetsCMSA(numSets,percent_ident,CMA,0.25);
#else
	set_typ	*CloseSq=0,InSet=MakeSet(N+4); FillSet(InSet);
	do {
	    CloseSq=RtnFastClustersCMSA(numSets,percent_ident,InSet,CMA);
	    if(percent_ident > 98) break;
       	    for(j=1; j <= numSets; j++){
		dd=(double) CardSet(CloseSq[j])/(double)N;
		if(dd > 0.5){
		  fprintf(stderr,"%d%c %d seqs: ",percent_ident,'%',N);
       	    	  for(k=1; k <= numSets; k++){ 
			fprintf(stderr,"%d ",CardSet(CloseSq[k]));
			NilSet(CloseSq[k]);
		  } fprintf(stderr,"\n");
		  free(CloseSq); CloseSq=0; percent_ident++; break; 
		}
	    }
	} while(CloseSq==0); NilSet(InSet);
#endif

	//============== 1b. cluster the sequences into more distant sets. ============
	set_typ	Used=MakeSet(SetN(CloseSq[1])); ClearSet(Used);
	set_typ	*ClstSet=0; NEW(ClstSet,numSets +4,set_typ);
	Int4	size,*Size; NEW(Size,numSets +4,Int4);
       	for(j=1; j <= numSets; j++){
#if 0	// DEBUG...
	   fprintf(stderr,"========== Set %d: ========\n  ",j);
	   PutSet(stderr,CloseSq[j]);
#endif
	   ClstSet[j]=MakeSet(Nset + 4); ClearSet(ClstSet[j]);
       	   for(size=0,n=1; n <= Nset; n++){
	   	s=cluster[n][1];
		if(MemberSet(s,CloseSq[j])){
		    AddSet(n,ClstSet[j]);
		    size+=cluster[n][0];
		    assert(!MemberSet(n,Used)); AddSet(n,Used);
#if 0	// DEBUG...
	  	    fprintf(stderr,"cluster %d (%ld): \n  ",n,CsqQuery[n]);
       		    for(i=1; i <= cluster[n][0]; i++){ 
			sq=cluster[n][i]; fprintf(stderr," %d",sq);
			// assert(MemberSet(sq,CloseSq[j])); // some may be in different sets!!
		    } fprintf(stderr,"\n  ");
#endif
		} 
	   } // fprintf(stderr,"(%d seqs)\n",size);
	   NilSet(CloseSq[j]);
	   Size[j]=size;
	   // if(CardSet(ClstSet[j]) <= 0){ }
	} free(CloseSq);
	fprintf(stderr,"CardSet(Used)=%d; Nset=%d; numSets=%d\n",CardSet(Used),Nset,numSets);
	NilSet(Used); 
	// exit(1);

	//============== 2. Order from largest to smallest clusters. ============
	dH = dheap(numSets+5,3); fprintf(stderr,"."); 
       	for(j=1; j <= numSets; j++){
		if(Size[j]==0) continue;
		insrtHeap(j,(keytyp)-Size[j],dH); // insert by cluster size...
	}

	//============== 3. Sample larger to smaller clusters. ============
       	for(n=1; (j=delminHeap(dH)) != 0; n++){ 
	        // fprintf(stderr,"cluster %d(%d) (%d seqs)\n",n,j,Size[j]);
#if 0	// Test inverse progressive alignment strategy (crude hieraln?);
	  if(Size[j] > 2){
		sprintf(str,"Junk%d",j);
		FILE *dfp=open_file(str,".cma","w");
		set_typ SqSet=MakeSet(NumSeqsCMSA(CMA)+4); ClearSet(SqSet);
		Int4	Sq;
     		for(Int4 xx=1; xx <= Nset; xx++){
	   	   if(MemberSet(xx,ClstSet[j])){
	   	   	size = cluster[xx][0]; assert(size > 0);
			for(Int4 jj=1; jj <= size; jj++){
			   Sq=cluster[xx][jj]; assert(Sq != 0);
			   AddSet(Sq,SqSet);
			}
		   }
		} PutInSetCMSA(dfp,SqSet, CMA); fclose(dfp); NilSet(SqSet);
	  }
#endif
	  if(this->SampleClustersOfClusters(ClstSet[j],cluster,Nset,offset,Temp,map)){
            if(map > best){ 
#if 0
	        fprintf(stderr,"Clstr-Of-Clstr %d(%d) (%d seqs): ",n,j,Size[j]);
		fprintf(stderr,"map improves from %.1f to %.1f (%.1f K)\n",best,map,Temp);
#endif
         	best=map; SaveBestCMSA(SSX->RtnCMA());  // Save this configuration.
            } 
	  } 
	} free(Size);

        //======================= 4. Free up memory. =======================
       	for(j=1; j <= numSets; j++){ NilSet(ClstSet[j]); } free(ClstSet);
        for(n=1; n <= Nset; n++){
	  free(cluster[n]);
	  if(offset[n]) free(offset[n]);
	  if(0 && j>=200){ 	// DEBUG!!!
       	    while((n=delminHeap(dH)) != 0){ free(cluster[n]); if(offset[n]) free(offset[n]); } 
	    break;	// WARNING: purify will report memory leak with this...
	  }
	} free(cluster); free(offset); Nildheap(dH); 
	// InitMAPCMSA(SSX->RtnCMA()); // reset to optimum configuration.

        //======================= 5. compute and return LPR. =======================
	double map0=SSX->GapMap(Temp);
	SSX->RestoreBest(); // reset to optimum CMA configuration.
	map=SSX->GapMap(Temp);
	fprintf(stderr,"map1 = %.1f; restored map = %.1f; best = %.1f\n",map0,map,best);
#if 0
	fprintf(stderr,"\tpost-cluster time: %d seconds (%0.2f minutes)\n",
                        time(NULL)-time1,(float)(time(NULL)-time1)/60.0);
	fprintf(stderr,"final temperature = %2f K; LLR= %.1f\n", Temp,map);
#endif
//exit(1);
	return best;
}

