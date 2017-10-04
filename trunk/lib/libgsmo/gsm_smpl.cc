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

double	gsm_typ::the_gibbs_sampler(gs_type G, FILE *rfp)
{
	Int4	n,iter,i,j,r,s,t,N=NumSeqsCMSA(G->cmsa);
	Int4	nruns=G->nruns,nconverg=G->nconverge;
	Int4	run,number,num_fail;
	double	TotLike,best,lastbest,map,map2,T;
	BooLean	flag;
	dh_type dH;
	char	stage;
	/**** DEBUG *****/
	BooLean	debug=FALSE; // debug=TRUE;
	UInt4	time1;
	Int4	*ncols,c;
	/**** DEBUG *****/
	BooLean	*moved;
	Int4	*sq_moves,sqmv,tot_moves,cmvN,tot_try;
	double	*prob,prb,max_prb,min_prb,ave_prb,max_max_prb;
	double	oldMap,newMap;
	keytyp	key;
   ppg_typ *ppg= new ppg_typ(G);

   FILE *fp=0;
   if(debug) {
	NEW(ncols, nBlksCMSA(G->cmsa)+2, Int4);
	time1=time(NULL);
	fp = open_file(NameSeqSetCMSA(G->cmsa),".prn","a");
   }

   //========================== 1. Set up memory ============================
   // fprintf(stderr,"nconv = %d\n", nconverg); 
   assert(N < INT4_MAX && N > 0);
   map=lastbest=-DBL_MAX; 
   ppg->Update('c');
   dH = dheap(N,3);
   NEW(sq_moves,N+2,Int4); NEW(prob,N+2,double);
   mem_typ memS(N,'S'); // Sites: number of actions == number of sequences.
   mem_typ memG(N,'G'); // Gaps: number of actions == number of gapped sequences

   unsigned short x=0,X=0;
   if(N < (Int4) (USHRT_MAX-20)) X= (unsigned short) N; else X = (USHRT_MAX-20);
   X = MINIMUM(unsigned short,X,nBlksCMSA(G->cmsa)*1000);	// no need to make this so huge!
   x = MAXIMUM(unsigned short,X/5,5);
   if(x >= X) X = x + 1;
   mem_typ memB(nBlksCMSA(G->cmsa),'B',x,X,0.10); // BLOCKS: # actions == # blocks.
   // mem_typ memB(nBlksCMSA(G->cmsa),'B',5,35,0.10); // BLOCKS: # actions == # blocks.
   mem_typ memC(5,'C',3,15,0.33); // Columns: 5 types of operations.
   memC.Clear(4); // set cardinality of these to zero.
   moved = new BooLean [nBlksCMSA(G->cmsa)+2]; 

   //========================== 2. Perform sampling ============================
   for(run = 1; (nruns <= 0 || run <= nruns); run++){
	if(debug) {
		if(run<=0) fprintf(fp,"************* tweaking **************\n");  
		else fprintf(fp,"************* new run **************\n");  
	}
	// fprintf(stderr,"."); 
	best=-DBL_MAX; 
	if(nruns <= 0) { 
	    iter=0; G->stage=stage=1; 
            cmvN=NumColumnsCMSA(G->cmsa); 
	    SetPseudoToMapCMSA(G->cmsa);
	    if(nruns == 0){
		best=map=lastbest=TotLikeGibbs(G); TotLike=0.0; 
		SaveBestCMSA(G->cmsa); // don't discard previous map for this one
	    } else { map=0.01; TotLike=0.0; iter=-2; SaveBestCMSA(G->cmsa); }
	} else {
		map=TotLike=0.0; 
		ppg->Update('i'); 
		G->stage=stage=0; iter=-3; cmvN=nBlksCMSA(G->cmsa); 
	}
	ppg->Update('p'); flag=FALSE;
	max_max_prb=-999999999.0; 
        /********************* MAIN Gibbs sampling loop ************************/
	for(number=num_fail=0; iter <= nconverg; iter++) {
        /********************* 2a. SLIDE COLUMNS SAMPLING ROUTINES *********************/
	  if(SmplBlksCols && map > 0.0){  
	  // if(map > 0.0 || !memC.Empty(4))  // at start only.
if(debug) std::cerr << "sliding columns loop\n";
	    for(t=1; t < nBlksCMSA(G->cmsa); t++){
		newMap=oldMap=0.0; T = 1.0;	// T = temperature.  // T = 0.5;
		// if(iter > 5) T = G->temp; else T = 0.5;
		if(dfsSlideColLtCMSA(t,&oldMap,&newMap,G->cmsa,T)){
			memC.Success(4); flag=TRUE;
         		if(newMap > best){ number = 0; best = newMap; }
		} else {
		   memC.Failure(4); 
         	   if(oldMap > best){ number = 0; best = oldMap; }
         	   if((nruns >= 0 || stage > 1) && map < best){ map=best; SaveBestCMSA(G->cmsa); }
		   newMap=oldMap=0.0;
		   if(dfsSlideColRtCMSA(t,&oldMap,&newMap,G->cmsa,T)){
			  memC.Success(4); flag=TRUE;
         		  if(newMap > best){ number = 0; best=newMap; }
		   } else { memC.Failure(4); if(oldMap > best){ number = 0; best=oldMap; } }
		}
         	if((nruns >= 0 || stage > 1) && map < best){ map=best; SaveBestCMSA(G->cmsa); }
	    }
	  }
	  /********************* 2b. COLUMN SAMPLING. **********************/
	  if(SmplBlksCols && map > 0.0){
	    if(debug) std::cerr << "starting column sampling loop\n";
	    for(n = 1; n <= cmvN; n++){
   	      if(!G->test[1]){
	        do { t=random_integer(nBlksCMSA(G->cmsa))+1;
	        } while(memB.Chance(t) < SampleUniformProb());
		switch(random_integer(4)){	// return 0..3
		  case 0: 
		    if(ContigBlkCMSA(t,G->cmsa)){
			if(ShiftCMSA(G->cmsa,t)){ memC.Success(1); flag=TRUE; }
			else { memC.Failure(1); }
		    } break;
		  case 1:
		    if(SampAddColTempCMSA(t,G->temp,G->cmsa)){ memC.Success(2); flag=TRUE; }
		    else { memC.Failure(2); } break;
		  case 2: // do more often.
		    if(SampRmColTempCMSA(t,G->temp,G->cmsa)){ memC.Success(3); flag=TRUE; }
		    else { memC.Failure(3); } break;
		  default: 
		    if(SampleEndColCMSA(t,G->temp,G->cmsa)){
				memC.Success(5); EnlargeGibbs(G); flag=TRUE;
		    } else { memC.Failure(5); } break;
		}
              }
	    }
	  } else for(t = 1; t <= nBlksCMSA(G->cmsa); t++){
	      if(ContigBlkCMSA(t,G->cmsa) && memC.Chance(1) >= SampleUniformProb()){
		if(ShiftCMSA(G->cmsa,t)){ memC.Success(1); flag=TRUE; }
	      }
	  }
	  if(flag) { ppg->Update('p'); flag=FALSE; }
	  // for(t=1; t<=nBlksCMSA(G->cmsa); t++){ std::cerr << "BLOCKS:";  memB.Put(stderr,t); }
	  // for(i=1; i<=4; i++) if(memC.Chance(i) > 0.5){ std::cerr << "COLUMNS:";  memC.Put(stderr,i); }
	  // if(G->test[1]) for(i=1; i<=4; i++)if(!memC.Empty(i)){ std::cerr << "COLUMNS:";  memC.Put(stderr,i); }
if(debug) if(stage > 0 && !memC.Empty(5)){ std::cerr << "END COLUMNS:";  memC.Put(stderr,5); }

	  if(G->mod_temp && G->temp0 > 0){
		// if(G->temp0 > 75.) SetTemperatureGibbs(G->temp0-1.0, G);
		if(G->temp0 > 75.){ SetTemperatureGibbs(G->temp0-2.0, G);
		    // fprintf(stderr,"*************** temperature = %2f K **************\n", G->temp0);
		}
	  } else {
	    /********************* 2d. BLOCK-BASED SEQ. SAMPLING. ********************/
if(map > 0.0){  	//=================== sticky sampling... =======================
#if 0
	  Int4	card,sq,size=(Int4) ceil((double)N*0.10);  // ten percent of the sequences.
	  set_typ Set=MakeSet(N+5); ClearSet(Set);
	  for(n=1;(card=CardSet(Set)) < N || n <= 12; n++){
	       if(n >= 12){
       	    	 for(sq = 1; sq <= N; sq++){
		   if(!MemberSet(sq,Set)){ ppg->Propagate(sq,&prb,moved); AddSet(sq,Set); break; }
		 } n=12;
	       } else { ppg->PropagateClusters(Set,size,&prb); }

               if((TotLike=TotLikeGibbs(G)) > best){ 
		best=TotLike;
		if((nruns >= 0 || stage > 1) && map < best){
		  map=best; SaveBestCMSA(G->cmsa);
		  if(G->mod_temp && G->temp0 > 0)
		     { fprintf(stderr,"map = %.2f @ %.2f K)\n", map,G->temp0); }
		}
	       } else {
		if(map > 0 && TotLike < 0.8*map){  // then start over...
		     if(AlignedCMSA(G->cmsa))	// if cmsa->best != null; start over...
			{ InitMAPCMSA(G->cmsa); ppg->Update('p'); }
	         }
	       } if(CardSet(Set) == N) break;
	   } NilSet(Set);
#elif 0
	   Int4 percent_ident=50; 
	   double MaxFrctn=0.20;
	   ppg->PropagateRepSets(percent_ident,MaxFrctn,TotLike);
           if((TotLike=TotLikeGibbs(G)) > best){ 
	     best=TotLike;
	     if((nruns >= 0 || stage > 1) && map < best){
		  map=best; SaveBestCMSA(G->cmsa);
		  if(G->mod_temp && G->temp0 > 0){ fprintf(stderr,"map=%.2f (%.2f K)\n",map,G->temp0); }
	     }
	   } else {
		if(map > 0 && TotLike < 0.8*map){  // then start over...
		     if(AlignedCMSA(G->cmsa)){ InitMAPCMSA(G->cmsa); ppg->Update('p'); }
	        }
	   }
#else
	   double  dd,MaxFrctn=0.20;
	   Int4    percent_ident=50,NumSets,sq,N=NumSeqsCMSA(G->cmsa);
	   set_typ *set=RtnTargetSizeSetsCMSA(NumSets, percent_ident,G->cmsa,MaxFrctn);
	   dh_type dH=0; 
	   for(i=1; i <= NumSets; i++){
	      x=CardSet(set[i]);
              if(x==1){ Int4 *Sq=ListSet(set[i]); ppg->Propagate(Sq[0],&prb,moved);  free(Sq); }
              else {
                for(dH=dheap(N+5,4),sq=1; sq <= N; sq++){
                   if(MemberSet(sq,set[i]))	// align worst first...
                      { insrtHeap(sq,(keytyp)GetTotalProbCMSA(sq, G->cmsa),dH); }
                } ppg->MultiPropagate(set[i],&dd,dH); Nildheap(dH);
              }
              if((TotLike=TotLikeGibbs(G)) > best){ 
	        best=TotLike;
	        if((nruns >= 0 || stage > 1) && map < best){
		  map=best; SaveBestCMSA(G->cmsa);
		  if(G->mod_temp && G->temp0 > 0){ fprintf(stderr,"map=%.2f (%.2f K)\n",map,G->temp0); }
	        }
	      } else {
		if(map > 0 && TotLike < 0.8*map){  // then start over...
		     if(AlignedCMSA(G->cmsa)){ InitMAPCMSA(G->cmsa); ppg->Update('p'); }
	        }
	      } NilSet(set[i]);
	   } free(set);
#endif
} else { 	//=================== sample one-at-a-time ======================
	    max_prb=-999999999.0; min_prb=+999999999.0;
   	    tot_try=tot_moves=0;
       	    for(n = 1; n <= N; n++) insrtHeap(n,(keytyp) Random(),dH); 
       	    while((n=delminHeap(dH)) != NULL){ 
#if 0
	      if(memS.Chance(n) >= SampleUniformProb())
	      {
         	if((sqmv=ppg->Propagate(n,&prb,moved)) > 0)
		{
		   for(t=1; t <= nBlksCMSA(G->cmsa); t++){
			if(moved[t]) memB.Success(t); else memB.Failure(t);
// memB.Put(stderr,t); 
		   } memS.Success(n);
		} else memS.Failure(n);
// if(memS.Chance(n) <= 0.50) memS.Put(stderr,n); 
// memS.Put(stderr,n); 
		if(max_prb < prb){
         	  if((max_prb=prb) > max_max_prb){
		    max_max_prb = max_prb;
         	    if((TotLike=TotLikeGibbs(G)) > best){ 
         		num_fail=number=0; best=TotLike; 
         		if((nruns >= 0 || stage > 1) && map < best){ 
			   map=best; SaveBestCMSA(G->cmsa); 
			   if(G->mod_temp && G->temp0 > 0){
				fprintf(stderr,"map = %.2f @ %.2f K\n", map,G->temp0);
			   }
			}
	            } else {
			num_fail++;
			if(map > 0 && TotLike < 0.8*map){  // then start over...
			   if(AlignedCMSA(G->cmsa)){ InitMAPCMSA(G->cmsa); ppg->Update('p'); }
			}
	            }
         	  }
         	} min_prb=MINIMUM(double,min_prb,prb);
         	sq_moves[n] = sqmv; prob[n] = prb;
         	tot_moves+=sqmv; tot_try+=nBlksCMSA(G->cmsa);
	     }
#else
             if((sqmv=ppg->Propagate(n,&prb,moved)) > 0){
         	if((TotLike=TotLikeGibbs(G)) > best){ 
         	   num_fail=number=0; best=TotLike; 
         	   if((nruns >= 0 || stage > 1) && map < best){ 
		     map=best; SaveBestCMSA(G->cmsa); 
		     if(G->mod_temp && G->temp0 > 0)
			{ fprintf(stderr,"map = %.2f @ %.2f K\n", map,G->temp0); }
		   }
	        } else if(map > 0 && TotLike < 0.8*map){  // then start over...
			if(AlignedCMSA(G->cmsa)){ InitMAPCMSA(G->cmsa); ppg->Update('p'); }
	        }
	     }
#endif
          } ave_prb=((max_prb-min_prb)/4.0) + min_prb;
}
	}
	/********************* checkpoint for map ************************/
	if(G->mod_temp && G->temp0 > 0){
	   if(G->temp0 > 75.){
		SetTemperatureGibbs(G->temp0-2.0, G);
	        // fprintf(stderr,"*************** temperature = %2f K **************\n", G->temp0);
	   }
	}
	if(iter > 0){ // simulated annealing option 
	     TotLike = TotLikeGibbs(G);
#if 0
	     fprintf(stderr,
		"** stage %d: blocks = %d; number = %d; TotLike = %g \n",
				stage,nBlksCMSA(G->cmsa), number,TotLike);
	     fprintf(stderr,
		"  max = %.2f; min %.2f; ave_prob = %.2f\n", max_prb,min_prb,ave_prb);
	     fprintf(stderr, " %.2f seq_moves (%d/%d)\n",
			(double)tot_moves/(double)tot_try,tot_moves,tot_try);
#endif
	     if(TotLike > best){ 
	       if(debug) {
		ConfigCMSA(ncols, G->cmsa);
		fprintf(fp,"%4d  %5.2f ",iter,TotLike);
		for(c=0, t=1; t<=nBlksCMSA(G->cmsa); t++){	 
			fprintf(fp,"%6.1f ",FieldRelMapCMSA(G->cmsa,t));
			c+=ncols[t];
		}
		fprintf(fp,"%3d ",c);
		for(t=1; t<=nBlksCMSA(G->cmsa); t++){ fprintf(fp,"%3d ",ncols[t]); }
		fprintf(fp,"\n");
	       }
		num_fail=number=0; best=TotLike; 
		// if(stage > 0 && map < best){ map=best; SaveBestCMSA(G->cmsa); }
		if((nruns >= 0 || stage > 1) && map < best){
			map=best; SaveBestCMSA(G->cmsa); }
	     } else {
		num_fail++;
		if(TotLike < 0.8*map){  // then start over...
		 if(AlignedCMSA(G->cmsa)){
if(debug) std::cerr << "starting over with previous best alignment...\n\n"; 
		  InitMAPCMSA(G->cmsa); ppg->Update('p');
		 }
		}
	     }
	}
	if(debug) fprintf(stderr,
		"********** %d: cols=%d; blks=%d; num_fail=%d; LLR=%.1f *********\n",
		run,NumColumnsCMSA(G->cmsa),nBlksCMSA(G->cmsa),num_fail, map);
/********************* check for convergence ************************/
	if(iter <= 0) number=0;
	// if(++number > G->limit || num_fail > NumSeqsCMSA(G->cmsa)) {
	if(++number > G->limit) {
		if(stage == 3) {
			G->stage++; stage++; break; 
		} else if(stage == 2) { 
			G->stage++; stage++; number=G->limit-3; G->mod_temp=2; G->temp0=0;
		} else if(stage==1) { G->stage=stage=2; 
			number=G->limit-3; 
			// max_max_prb=-9999999999999.; map = 0.0;
			// number=0; // best=map=0.0;
		} else {  /** stage == 0 **/
		   G->stage=stage=1; number=0; 
		   cmvN=NumColumnsCMSA(G->cmsa);
		   SetPseudoToMapCMSA(G->cmsa);
		   best=TotLikeGibbs(G); 
		   if((nruns >= 0 || stage > 1) && map < best){
			map=best; SaveBestCMSA(G->cmsa); }
		   if(map < G->stage1_minmap){ map=0.0; break; }
		   if(nruns > 0 && G->test[1]) break;
                   if(best < 0){ /** THEN ABORT **/ break; }
		}
	}
#if 0
	  if(nruns > 0 && G->test[1]){
		if(memC.Empty(1) && memC.Empty(2) && memC.Empty(3)){
		   SetPseudoToMapCMSA(G->cmsa);
		   best=TotLikeGibbs(G); 
		   if((nruns >= 0 || stage > 1) && map < best){
			map=best; SaveBestCMSA(G->cmsa); }
		   break; 
		}
	  }
#endif
      } // end of 'for(number=num_fail=0; iter <= nconverg; iter++)'
	if(best==lastbest || nruns <= 0) break; 
	if(lastbest < best) lastbest = best;
    } // end of 'for(run = 1; (nruns <= 0 || run <= nruns); run++)'
    // if(G->mod_temp && G->temp0 > 0) fprintf(stderr,"final temperature = %2f K\n", G->temp0);
    Nildheap(dH); free(sq_moves); free(prob);
    if(debug) {
        if(map > 10.0) fprintf(stderr," %.1f\n", map);
	fclose(fp); 
	free(ncols);
	if(nruns > 0){
	    fp = open_file(NameSeqSetCMSA(G->cmsa),".log","a");
	    fprintf(fp,"%5.1f %3d %3ld %3d %3d\n",
		map,iter,time(NULL)-time1, c,nBlksCMSA(G->cmsa));
	    fclose(fp);
  	}
    }
    delete [] moved;
    delete ppg;
    return map;
}


