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

#define debug_indels 0

//************************ Single Sequence Sampling Routines ***********************
BooLean	gmb_typ::SampleSingle(Int4 sq,double Temp,double &MAP)
// sample insertions and deletions in sequence sq.
{
	Int4	i,start,trace_length,score,*newpos,*oldpos;
	cma_typ CMA=SSX->RtnCMA(); 
	MAP=SSX->GapMap(Temp);

#if debug_indels
	// if(sq == 15184 || sq==2099 || sq==134) this->PutCMA_SqIndels(stderr,sq);
	if(sq==8765) this->PutCMA_SqIndels(stderr,sq);
#elif 0
	if(sq == 15184)
	{ FILE *tfp=open_file("junk_15184",".cma","w"); PutCMSA(tfp,CMA); fclose(tfp); }
	if(sq==56)
	{ FILE *tfp=open_file("junk_56",".cma","w"); PutCMSA(tfp,CMA); fclose(tfp); }
#endif
	// if(sq == 15184) { this->PutCMA_SqAln(stderr,sq); }

	//======== 1. Remove sq from alignment, sample scoring matrix and realign. =====
	assert(SeqIsInCMSA(sq,CMA));
	oldpos=SSX->RemoveFromAlign(sq); 
	gss_typ *gss=gssCMSA(CMA);
	e_type	sbjE = gss->TrueSeq(sq); 	
	SSX->InitNDL(Temp);	// initialize JHL parameters?
	char	*operation=SSX->GapAlnTrace(sbjE,Temp,start,score);
	trace_length=strlen(operation);

	//========== 2. Create a gapped sequence using 'operation'. ===========
	gsq_typ	*ogsq,*gsq = new gsq_typ[1]; NEW(newpos,nBlksCMSA(CMA)+3,Int4);  
        gsq->initialize(gss->LeftFlank(),gss->RightFlank(),
				operation,trace_length,start,sbjE,newpos);
#if debug_indels	// DEBUG...
	// if(sq == 469)
	// if(sq == 15184 || sq==2099 || sq==134)
	// if(sq == 1235)
	// if(sq == 7 || sq==664|| sq==1110 || sq == 1566 || sq==5216 || sq==5952)
#elif 0
	if(sq == 56)
	{
		std::cerr << "========= sampled ==========" << std::endl;
		fprintf(stderr,"start=%d; score=%d\n",start,score); 
		fprintf(stderr,"operation=%s\n",operation);
		this->PutGappedAln(stderr,sbjE,operation,start);
		this->PutCMA_SqIndels(stderr, sq, newpos, gsq);
		char *ooper=gsq->Operation3(stderr,AB);
		fprintf(stderr,"ooper=%s\n",ooper);
		this->PutGappedAln(stderr,sbjE,ooper,start);
		free(ooper);
		// gsq->Put(stderr,AB);
		// std::cerr << "========= AlignSeq ==========" << std::endl;
		// this->AlignSeq(stderr,'S',sbjE,Temp);
		fprintf(stderr,"old gsq:\n");
		ogsq=gss->GetGSQ(sq); ogsq->Put(stderr,60,AB);
		fprintf(stderr,"new gsq:\n");
		gsq->Put(stderr,60,AB);
		// exit(1);
	}
#endif
	free(operation); 

	//========== 3. If qsq is unchanged then retain old alignment and return. =========
	if(gss->Identical(sq,*gsq)){
#if 0
	if(sq == 56){
		// gsq->Put(stderr,AB);
		fprintf(stderr,"sequence 56 new alignment is identical\n");
		fprintf(stderr,"oldpos=%d; newpos=%d\n",oldpos[1],newpos[1]);
		gsq_typ *xgsq=gss->GetGSQ(sq);
		xgsq->Put(stderr,60,AB);
		exit(1);
	}
#endif
	   BooLean the_same=TRUE;
	   for(Int4 b=1; b<= nBlksCMSA(CMA); b++)
	      { if(oldpos[b] != newpos[b]){ the_same=FALSE; break; } }
	   if(the_same){
		SSX->AddToAlign(sq,oldpos); free(newpos); free(oldpos); delete []gsq; 
	        return FALSE; 
	   }
	}

	//========== 4. Otherwise sample one of the two alignments. =========
	ogsq=SwapGsqCMSA(sq,gsq,CMA); SSX->AddToAlign(sq,newpos); free(newpos); 
	// ExtendFakeToRealCMSA(sq,CMA);
#if 0
	gsq_typ *gsq0=gsqCMSA(sq,CMA);
        newpos=GetPosSitesCMSA(sq,CMA);
        gsq0->Put_cma_format(stderr,sq,nBlksCMSA(CMA),newpos,LengthsCMSA(CMA),AB);
	free(newpos);
#endif
        double	newmap=SSX->GapMap(Temp); 
	// if(sq == 15184) { this->PutCMA_SqAln(stderr,sq); }
#if debug_indels	// DEBUG...
	// if(sq == 469)
	// if(sq == 15184 || sq==2099 || sq==134)
	if(sq == 8765)
#elif 0
	if(sq == 56)
	{
		this->PutCMA_SqIndels(stderr,sq);
		fprintf(stderr,"MAP=%.2f; newMAP=%.2f.\n",MAP,newmap);
		gsq_typ *xgsq=gss->GetGSQ(sq);
		xgsq->Put(stderr,60,AB);
		exit(1);
	}
#endif
	if(SampleTransition(newmap, MAP, Temp)){
		free(oldpos); delete []ogsq; MAP=newmap; return TRUE;
	} else {
		SSX->RmFromAlign(sq); assert(ogsq); ReplaceCMSA(sq,ogsq,CMA);
		SSX->AddToAlign(sq,oldpos); free(oldpos);
	// if(sq == 15184) this->PutCMA_SqAln(stderr,sq);
		return FALSE; 
	}
}

double	gmb_typ::GambitSingles(char mode,double temp)
// sample with gaps using simulated annealing...
{
	cma_typ CMA=SSX->RtnCMA(); assert(CMA != 0); 
	Int4	n,i,j,r,s,t,N=NumSeqsCMSA(CMA);
	double	best,map;
	FILE	*efp=0; // efp=stderr;

	SetPseudoToMapCMSA(CMA); best=SSX->GapMap(temp); 
	// fprintf(stderr,"best map = %.2f; Temp = %.1f\n",best,temp);
	SaveBestCMSA(CMA); // don't discard previous map for this one
	if(efp) fprintf(stderr,"best map = %.2f; Temp = %.1f\n",best,temp);
// #if debug_indels
	dh_type dH = dheap(N+4,3); // fprintf(stderr,"."); 
       	for(n = 1; n <= N; n++){
#if 1	// working on sampling over subalignment...
	   if(!SeqIsInCMSA(n,CMA)) continue;
#else
	   assert(SeqIsInCMSA(n,CMA));
#endif
	   insrtHeap(n,((keytyp)Random()),dH);
	}
	// chkeyHeap(15184,(keytyp) INT4_MIN, dH);
	// chkeyHeap(15485,(keytyp) INT4_MIN, dH);
       	for(j=1; (n=delminHeap(dH)) != 0; j++){ 
	  // fprintf(stderr,"--------------- %d: sq=%d ------------------\n",j,n);
	  if(SampleSingle(n,temp,map)){
            if(map > best){ 
		if(efp) fprintf(stderr,"%d(%d): map improves from %.1f to %.1f (%.1f K) gap_LLR=%.2f\n",
			j,n,best,map,temp,SSX->RtnIndelPenalty());
         	best=map; SaveBestCMSA(SSX->RtnCMA());
            }
	  }
// fprintf(stderr,"%d(%d): map changes from %.1f to %.1f (%.1f K)\n",j,n,best,map,temp); exit(1);
	  if(0 && j >= 5){ 	// DEBUG...
       	    while((n=delminHeap(dH)) != 0) ;  
	    break;	// WARNING: purify will report memory leak with this...
	  }
	} // fprintf(stderr,"final temperature = %.1f K\n %.1f\n", temp,best);
	Nildheap(dH); 
	// InitMAPCMSA(SSX->RtnCMA()); 
	// SSX->UpdateWts(); // don't update these as may reject correct, but uncertain alignments.
	double map0=SSX->GapMap(temp);
	SSX->RestoreBest(); // reset to optimum CMA configuration.
	map=SSX->GapMap(temp);
	if(efp) fprintf(stderr,"map1 = %.1f; restored map = %.1f; best = %.1f\n",map0,map,best);
	return best;
}

BooLean gmb_typ::SampleTransition(double nmap, double omap, double temp)
/** short routine to be used for sampling **/
{
        double  ratio,rand;

	if(temp <= 0.0){ if(nmap > omap)  return TRUE; else return FALSE; }
	temp=300.0/temp;   // 1/temp = 1.0 for "room temperature" and infinite for absolute zero.
	// fprintf(stderr,"Temperature = %.1f\n",temp);
	assert(temp >= 0.5);  // <= 600 K from calling...
	if(temp < 30.0){  if(nmap > omap)  return TRUE; else return FALSE; } // lowest sampling temp=30.
        // if(fabs(nmap-omap) > 500) print_error("error in SampleTransition( )");
        if((nmap-omap) >= 20)  return TRUE;
        if((omap-nmap) >= 20)  return FALSE;
        if(omap < nmap){
                ratio = exp(nmap-omap); // ratio of new to old == nLR/oLR == likelihood ratio.
                if(temp != 1.0) ratio = pow(ratio,temp);
#if 1
                ratio = ((double)Random()/(double)RANDOM_MAX) * (ratio+1);
#else
                rand = ((double)Random()/(double)RANDOM_MAX);
		fprintf(stderr,"oLLR=%.2f; nLLR=%.2f; Temp = %.2f; ratio=%.2f; rand=%.2f\n",
			omap,nmap,temp,ratio,rand);
		ratio = rand * (ratio+1);
#endif
                if(ratio >= 1.0) return TRUE; else return FALSE;
        } else {        /** omap >= nmap **/
		assert(temp > 0.10);
                ratio = exp(omap-nmap); /** ratio of old to new **/
                if(temp != 1.0) ratio = pow(ratio,temp);
#if 1
                ratio = ((double)Random()/(double)RANDOM_MAX) * (ratio+1);
#else
                ratio = ((double)Random()/(double)RANDOM_MAX) * (ratio+1);
		fprintf(stderr,"oLLR=%.2f; nLLR=%.2f; Temp = %.2f; ratio=%.2f\n",omap,nmap,temp,ratio);
#endif
                if(ratio > 1.) return FALSE; else return TRUE;
        }
}

//************************ Clustered sampling routines ***********************
BooLean	gmb_typ::SampleCluster(Int4 *cluster, Int4 *offset,e_type csqE,double Temp,double &MAP)
// sample a gapped alignment for a cluster of sequences.
{
	Int4	start,strt,trace_length,score,sq,ii,*newpos,**oldpos,size=cluster[0];
	BooLean	Modified=FALSE,found=FALSE;
	cma_typ CMA=SSX->RtnCMA(); 
	double	newmap; MAP=SSX->GapMap(Temp); 

	//======= 1. Remove seq cluster from the alignment and sample matrix.  ==========
	gsq_typ *gsq,**ogsq; NEW(ogsq,size+2,gsq_typ *);
	oldpos=SSX->RemoveFromAlign(cluster); 

	//========== 2. Sample a gapped alignment for master sequence csqE. ===========
#if 0	// DEBUG...
	SSX->InitNDL(Temp);	// initialize JHL parameters?
#endif
	//================ New 'Repeated-Domain' method. =================
	this->GetMultipleTraces(csqE,Temp);	// find all repeats....
	//  OPS->Put(stderr);
	// char	*operation=SSX->GapAlnTrace(csqE,Temp,start,score);

// fprintf(stderr,"%d: operation=%s\n",sq,operation);

	// ========== 3. Align sequences based on the master. ===========
	gss_typ	*gss=gssCMSA(CMA);  
	Int4 os;
	// Int4 LenHMM=TotalLenCMSA(CMA);
	for(ii=1; (sq=cluster[ii]) != 0; ii++) {
	   e_type sbjE=gss->TrueSeq(sq); 	
	   os=offset[ii];
	   char *operation,*op,*bst_op=0;
	   Int4 bst_strt=0,bst_scr=0,tmp_scr,x;
	   for(Int4 J=1; J <= OPS->N; J++){
	      operation=OPS->operation[J];
	      assert(operation[1]=='M' || operation[1]=='D');
	      start=OPS->start[J];
	      strt=start-os; 

              // if((strt + LenHMM -1) <= 5) continue;   // domain maps prior to sbjE!
              // else if((LenSeq(sbjE) - start) < 5) continue; // domain maps after sbjE!

	      op=ConvertOperation(operation,sbjE,start,os);
	      tmp_scr=this->AlignLenOperation(op);
	      if(strt < 1) x=1; else x=strt;
	      // tmp_scr=this->ScoreSeqAlnSMatrix(op,sbjE,x);
	      if(strt < 1){ 	// starts < 1 already modeled in the operation array ("EDdddmm...");
		if(bst_scr < tmp_scr){
		  if(bst_op) free(bst_op);
		  bst_strt=1; bst_scr=tmp_scr; bst_op=op; 
		} else free(op);  
	      } else if(strt >= LenSeq(sbjE)){
		  // fprintf(stderr,"%d: os = %d; strt=%d\n",sq,os,strt);
                  // fprintf(stderr,"op = %s\n",op); PutSeq(stderr,sbjE,AB);
		  free(op);
	      } else if(bst_scr < tmp_scr){
		  if(bst_op) free(bst_op);
		  bst_strt=strt; bst_scr=tmp_scr; bst_op=op; 
	      } else free(op); op=0;
	   } op=bst_op; strt=bst_strt;
#if 0	// DEBUG...
	   fprintf(stderr,"%d: os = %d; strt=%d\n",sq,os,strt);
	   fprintf(stderr,"op = %s\n",op); PutSeq(stderr,sbjE,AB);
#endif
	   if(op==0){
	      // fprintf(stderr,"%d: strt=%d\n",sq,strt);
	      SSX->AddToAlign(sq,oldpos[ii]); 
	   } else {
              assert(op[0]=='E');
	      NEW(newpos,nBlksCMSA(CMA)+2,Int4);  // for new sites.
	      gsq = new gsq_typ[1]; // NEEDS TO BE AN ARRAY FOR gss_typ!!!
	      trace_length=strlen(op);
	      gsq->initialize(gss->LeftFlank(),gss->RightFlank(),op,trace_length,strt,sbjE,newpos);
	      if(gss->Identical(sq,*gsq)){
	          BooLean the_same=TRUE;
	   	  for(Int4 b=1; b<= nBlksCMSA(CMA); b++)
	      	     { if(oldpos[ii][b] != newpos[b]){ the_same=FALSE; break; } }
	   	  if(the_same){ delete []gsq; }
	          else { Modified=TRUE; ogsq[ii]=SwapGsqCMSA(sq,gsq,CMA); }
	      } else { Modified=TRUE; ogsq[ii]=SwapGsqCMSA(sq,gsq,CMA); }
	      SSX->AddToAlign(sq,newpos); free(newpos); free(op); 
	   }
	} if(OPS) delete OPS; OPS=0;

	//========== 4. Sample one of the two alignments. =================
	if(Modified) newmap=SSX->GapMap(Temp);
	//if(found){ this->PutClusterAln("junk_debug",cluster,CMA); } 
	if(Modified && this->SampleTransition(newmap, MAP, Temp)){
	   for(ii=1; (sq=cluster[ii]) != 0; ii++) {
		if(ogsq[ii]) delete []ogsq[ii]; free(oldpos[ii]); 
	   } free(oldpos); MAP=newmap; free(ogsq); return TRUE;
	} else {	// reject the new cmsa file.
	   SSX->RmFromAlign(cluster);
	   for(ii=1; (sq=cluster[ii]) != 0; ii++) {
		if(ogsq[ii]){ ReplaceCMSA(sq,ogsq[ii],CMA); }
		SSX->AddToAlign(sq,oldpos[ii]); free(oldpos[ii]); 
	   } free(ogsq); free(oldpos); return FALSE; 
	}
}

double	gmb_typ::GambitClusters(char mode,Int4 similarity,double Temp)
// sample with gaps using simulated annealing...
// 	char	method='E'; // gapped E-value
{
	cma_typ CMA=SSX->RtnCMA(); assert(CMA != 0); 
	Int4	n,i,j,r,s,t,N=NumSeqsCMSA(CMA);
	double	best,map;
	FILE	*efp=0; // efp=stderr;
	dh_type dH;

	SetPseudoToMapCMSA(CMA); best=SSX->GapMap(Temp); SaveBestCMSA(CMA); // Save initial alignment.

	//============== 1. cluster the sequences into closely related sets. ============
	Int4    Nset,**cluster,**offset;
	cluster = this->GetClusters(efp,offset,Nset);
// fprintf(stderr,"Need to comment out exit(1) in gmb_typ::GambitClusters() in gmb_smpl.cc\n\n");
// exit(1);
	if(Nset < this->MinNumClstrs){
	  if(efp) fprintf(stderr,"too few (%d < %d) clusters",Nset,this->MinNumClstrs); 
       	  for(n=1; n <= Nset; n++){ free(cluster[n]);  if(offset[n]) free(offset[n]); }
	  free(cluster); free(offset); 
          return GambitSingles(mode,Temp);
	}

	//============== 2. Order from largest to smallest clusters. ============
   	Int4	time1=time(NULL);
	dH = dheap(Nset+5,3); if(efp) fprintf(stderr,"."); 
       	for(n=1; n <= Nset; n++){
		// this->PutClusterIDs(stderr,n,cluster[n]);
		if(cluster[n][0] == 1) insrtHeap(n,((keytyp)abs(Random())),dH);
		else insrtHeap(n,((keytyp)-cluster[n][0]),dH); // insert by group size...
	}

	//============== 3. Order from largest to smallest clusters. ============
       	for(j=1; (n=delminHeap(dH)) != 0; j++){ 
	  BooLean success;
	  if(cluster[n][0] > 1){	// use only for multiple sequences...
	     assert(offset[n]); assert(CsqQuery[n]);
	     success=this->SampleCluster(cluster[n],offset[n],CsqQuery[n],Temp,map);
	  } else success=SampleSingle(cluster[n][1],Temp,map);
	  // fprintf(stderr,"cluster %d(%d) (%d seqs).\n",j,n,cluster[n][0]);
	  if(success){
            if(map > best){ 
	        if(efp){
		   fprintf(stderr,"cluster %d(%d) (%d seqs): ",j,n,cluster[n][0]);
		   fprintf(stderr,"map improves from %.1f to %.1f (%.1f K)\n",best,map,Temp);
		} best=map; SaveBestCMSA(SSX->RtnCMA());  // Save this configuration.
            } 
	  } free(cluster[n]);
	  if(offset[n]) free(offset[n]);
	  if(0 && j>=200){ 	// DEBUG!!!
       	    while((n=delminHeap(dH)) != 0){ free(cluster[n]); if(offset[n]) free(offset[n]); } 
	    break;	// WARNING: purify will report memory leak with this...
	  }
	} free(cluster); free(offset);
	Nildheap(dH); 
	// InitMAPCMSA(SSX->RtnCMA()); // reset to optimum configuration.

	double map0=SSX->GapMap(Temp);
	SSX->RestoreBest(); // reset to optimum CMA configuration.
	map=SSX->GapMap(Temp);
	if(efp){
	  fprintf(stderr,"map1 = %.1f; restored map = %.1f; best = %.1f\n",map0,map,best);
	  fprintf(stderr,"\tpost-cluster time: %d seconds (%0.2f minutes)\n",
                        time(NULL)-time1,(float)(time(NULL)-time1)/60.0);
	  fprintf(stderr,"final temperature = %2f K; LLR= %.1f\n", Temp,map);
	} return best;
}

cma_typ	gmb_typ::Sample(char *name, char mode,Int4 similarity,double startTemp,FILE *rfp)
{ return this->Sample(name, mode,similarity,startTemp,0.0,30,rfp); }

cma_typ	gmb_typ::Sample(char *name, char mode,Int4 similarity,double startTemp,
		double endTemp, double incTemp,FILE *rfp)
{
   Int4		iter=0;
   double	Temperature=startTemp,map,bestmap=-999999999999999.9;
   cma_typ	CMA=0;
   Int4		T,timeS=time(NULL);
   FILE		*efp=0; // efp=stderr;
	
   do {
        if(similarity < 0) map=GambitClusters(mode,similarity,Temperature);
        else if(similarity > 0){
	  if(1 || iter==0) map=ClustersOfClustersSampling(mode,Temperature,similarity);
	  // if(iter % 5==0) map=ClustersOfClustersSampling(mode,Temperature,similarity);
	  // else if(iter % 5==4) map=GambitSingles(mode,Temperature);  // no sense in doing this really.
	  else map=GambitClusters(mode,similarity,Temperature);
	} else map=GambitSingles(mode,Temperature);
        if(bestmap < map){
		CMA=SSX->RtnCMA(); SaveBestCMSA(CMA);
		if(name){
                  // sprintf(str,"%s.tmp.cma",name); WriteCMSA(str,CMA);
		  if(0 && NumSeqsCMSA(CMA) < 10000){
                    sprintf(str,"%s.new",name); WriteMtfCMSA(str, CMA, NULL);
		  }
		} bestmap=map; if(efp) fprintf(stderr,"########## map %.1f ###########\n",map);
        } else if(Temperature <= endTemp) break;
	else Temperature=endTemp; // drop down to zero if failed to improve.
	if(startTemp < endTemp) break;
	Temperature = MAXIMUM(double,0.0,Temperature - incTemp); // 30 == 10 cycles to zero.
        iter++; if(efp) fprintf(stderr,"********** end iter %d ***********\n",iter);
	// break; // SSX->UpdateSeqWts(); // gets updated by SSX->RestoreBest();
  } while(iter < MaxIter);
  this->RestoreBest();
  if(efp) fprintf(stderr,"********** end iter %d ***********\n",iter);
  CMA=SSX->RtnCMA();  T=time(NULL);
  if(rfp){
	fprintf(rfp,"   Sample Seqs (%d cols): LLR = %.2f. %d secs (%0.2f mins)(%.1f..%.1f K)\n",
		NumColumnsCMSA(CMA),bestmap,T-timeS,(float)(T-timeS)/60.0,startTemp,endTemp); fflush(rfp);
  } return CMA;
}

