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

BooLean	gsm_typ::Add(Int4 minlength)
// TRY ADDING BLOCKS with minlength contiguous columns
{
	BooLean improved=FALSE;
	cma_typ	ma2;

	fprintf(stderr,"%d blocks; try adding blocks...\n",nBlksCMSA(cmsa));
	for(Int4 end=nBlksCMSA(cmsa), t=0; t <= end; t++){
	   if((ma2=AddBlkCMSA(t, minlength, cmsa))!= NULL){
	     fprintf(stderr,"block added between motifs %d & %d?\n",t,t+1);
	     if(gibbs(ma2)) improved=TRUE;
           } 
	}
	return improved;
}

BooLean	gsm_typ::DeletePoor(char tweakmode)
// Delete block b from the cmsa alignment...
{
	Int4	b,n=0,len;
	cma_typ	cma=cmsa,xcma;
	double	d;
	if(cma){
	   for(b=nBlksCMSA(cma); b > 0; b--){
		if(nBlksCMSA(cma) <= 1) break;
		len=LengthCMSA(b,cma); d=FieldRelMapCMSA(cma,b);
		// if(d <= 0.0 || (d < 10.0 && len < 4))
		if(d <= 0.0){
		   xcma=DeleteBlkCMSA(b,cma); 
		   if(cma!=cmsa) NilCMSA(cma); cma=xcma; xcma=0;  n++;
		}
	   }
	}
	if(n > 0){
	   // fprintf(stderr,"%d blocks deleted. Tweaking alignment...\n",n);
	   map=RelMapCMSA(cma); NilCMSA(cmsa); cmsa=cma; bestmap=0.0; 
	   SaveBestCMSA(cmsa); Record( ); gibbs(NULL,tweakmode);  // tweak the best alignment.
	   // fprintf(stderr,"Done Tweaking alignment\n");
	   return TRUE;
	} return FALSE;
}

BooLean	gsm_typ::Delete(char tweakmode)
/********************* TRY DELETING BLOCKS *************************/
{
	BooLean improved=FALSE;
	cma_typ	maD;

	fprintf(stderr,"delete blocks?\n");
        while(nBlksCMSA(cmsa) > 1 && (maD=DeleteBlkCMSA(cmsa))!=NULL){
	     // fprintf(stderr,"Block deleted. Tweaking alignment...\n");
  	     map = RelMapCMSA(maD); NilCMSA(cmsa); cmsa=maD; 
	     SaveBestCMSA(cmsa); Record( ); gibbs(NULL,tweakmode);  // tweak the best alignment.
	     // fprintf(stderr,"Done Tweaking alignment\n");
	     improved = TRUE;
        } return improved;
}


BooLean	gsm_typ::Split(Int4 minlength)
/********** TRY Splitting BLOCKS with >= minlength columns ***********/
{
	BooLean improved=FALSE;
	double	netlen,p;
	cma_typ	maS;

	fprintf(stderr,"Try splitting blocks...(%d total blks)\n",nBlksCMSA(cmsa));
	for(Int4 end=nBlksCMSA(cmsa), t=1; t <= end; t++){
	  if(LengthCMSA(t,cmsa) < minlength) continue;
	  netlen= (double) (LengthCMSA(t,cmsa) - 3);
	  p = (netlen*netlen + 1.0)/400.0;
	  // len=8: p = 0.065; len=12: p = 0.20; len=18: p = 0.56; len=22: p = 1.0.
	  if(netlen >= minlength && p >= SampleUniformProb( )){
	   if((maS=SplitBlkCMSA(t, minlength, cmsa))!= NULL){
	     fprintf(stderr,"block %d split in two.\n",t);
	     if(gibbs(maS,'n')) improved = TRUE;
	     // if(gibbs(maS,'S')) improved = TRUE;
           }
	  }
	} return improved;
}

BooLean	gsm_typ::Fuse(Int4 maxlengs)
// TRY Fusing BLOCKS with >= 16 columns.
// Sample fusions proportional to number of intervening residues in gaps.
{
	cma_typ	maF;
	double	avegap,p;
	BooLean improved=FALSE;

	fprintf(stderr,"Fuse blocks?\n");
	for(Int4 end=nBlksCMSA(cmsa), t=1; t <  end; t++){
	  avegap=(double)ResBetweenCMSA(t,cmsa)/(double)NumSeqsCMSA(cmsa);
	  p = 1.1/(avegap*avegap + 0.1); // always try fusing blocks if avegap <= 1.0
	  // avegap <=1.0 --> p = 1.0; avegap=2.0 --> p = 0.25; avegap=10 --> p = 0.01;
	  if(p >= SampleUniformProb( )){
	   if((maF=FuseBlksCMSA(t, maxlengs, cmsa))!= NULL){
	    fprintf(stderr,"blocks %d & %d fused.\n",t,t+1);
	    if(gibbs(maF)) improved = TRUE;
	   }
	  }
        }
	return improved;
}

double	gsm_typ::Breed( )
{
   double	Map,map1,map2,*mapn,*mapo;
   cma_typ	tmp_ma,ma,ma1,ma2,*New,*old;
   Int4		N,n,nn,no,i,j,o,c=1;
   BooLean	*recomb;
   Int4		last_time;

   last_time=time(NULL);
   if(use_gseq) strcpy(options,"-t1 -g -l1 ");
   else strcpy(options,"-t1 -l1 ");
   N = 2*SizeMSAHeap(maH); 
   NEW(New,N+2,cma_typ); NEW(old,N+2,cma_typ);
   NEW(mapn,N+2,double); NEW(mapo,N+2,double); NEW(recomb,N+2,BooLean);
   PurgeMSAHeap(maH); // remove all but one of those items with identical keys
   for(n=0,c=1; !EmptyMSAHeap(maH); c++){
	for(no=0,o=1; o <= n; o++){  // n==0 --> skipped the first time around 
	   if(recomb[o]){ NilCMSA(New[o]); } // discard parent alignments 
	   else { ++no; old[no] = New[o]; mapo[no]=mapn[o]; }
	}
        for(nn=0; (ma=DelMinMSAHeap(&Map, maH)) != NULL; )
		{ ++nn; New[nn]=ma; mapn[nn] = Map; recomb[nn]=FALSE; }
	for(n=nn,o=1; o <= no; o++) // also add previous recombinants to heap
		{ ++n; New[n]=old[o]; mapn[n]=mapo[o]; recomb[n]=FALSE; }
	fprintf(stderr,
	   "******** %d: %d new, %d old, %d total alignments ********\n",
		c,nn,no,n);
	if(no==0) nn--;
	for(i=1; i <= nn; i++){  /** try every combination of recombinants **/
	   ma1 = New[i]; map1=mapn[i];
	   for(j=i+1; j <= n; j++){
	     ma2 = New[j]; map2=mapn[j];
	     fprintf(stderr,"x");
	     if(DataCMSA(ma1) != DataCMSA(ma2)) ma=GRecombineCMSA(ma1,ma2); 
	     else ma=RecombineCMSA(ma1,ma2);
	     if(ma != NULL){
		 // if get a recombinant then try to refine it. 
		 recomb[i]=recomb[j]=TRUE;
		 fprintf(stderr,"RecombineCMSA( ) map = %g\n",map2=RelMapCMSA(ma));
		 if(map2 > 0.0){
		   Map=core_gibbs(&ma,'S',150);
		   if(Map > bestmap){
	        	fprintf(stderr,"!!!map improved from %.1f to %.1f\n", map2,Map);
			tmp_ma=cmsa; cmsa=ma; map=Map; 
			Record( ); last_time=time(0);
			cmsa=tmp_ma; 
		   } 
		   if(!KeyInMSAHeap(Map, maH)){
			if(InsertMSAHeap(ma, Map, maH) == NULL) NilCMSA(ma);
			else { fprintf(stderr,"\thybrid map=%f\n",Map); }
		   } else NilCMSA(ma);
		 } else NilCMSA(ma);
	     } 
   	   }
	}
	fprintf(stderr,"\n");
	if(ConvergedMSAHeap(maH) < 2.0) break; // I haven't tested this at all!!!!
	
	if((time(0)-last_time) > patience) break; // give up after 30 minutes...
	// ConvergedMSAHeap(maH);
   }
   for(o=1; o <= n; o++) {
	/** after breeding move the parent population to the heap **/
	if(InsertMSAHeap(New[o],mapn[o],maH)==NULL) NilCMSA(New[o]);
   }
   PurgeMSAHeap(maH);
   // get the Map for the best alignment.
   ma=DelMinMSAHeap(&Map,maH); 
   if(InsertMSAHeap(ma,Map,maH)==NULL) NilCMSA(ma);
   // get the map for the best alignment.
   free(New); free(old); free(mapo); free(mapn); free(recomb);
   return map;
}

void	gsm_typ::create_population()
// BUILD INITIAL POPULATION (MSAHEAP).
{
   Int4			futile_attempts=0;
   UInt4	time0,time0b;

   Int4	t,aln,NumAln;
   cma_typ	*list,lastcma=0;
   double	avelen=0.0,*lpr;

   maH=MkMSAHeap(mhpsz);
   item=NULL; time0b=time(NULL); 
   if(cmsa != NULL){ NilCMSA(cmsa); cmsa= NULL; }
   for(cycle=1; cycle <= maxcycle; cycle++){
	time0=time(NULL); 
	cmsa = NULL; gibbs(NULL);
	if(cmsa != NULL && map > 0) {  // THEN ADD ALIGNMENT TO HEAP.
	// if(create_align() && map > 0) // ADD ALIGNMENT TO HEAP.
	    // if(cycle==1) Refine( );  // Delete bad blocks and tweak alignment.
	    // if(cycle==1) Delete( );  // Delete and tweak alignment.
#if 0
	    map0=ReAlignBestCMSA(cmsa);
	    if(map0 > map){ map=map0; Record( ); }
	    else InitMAPCMSA(cmsa);
#endif
	    // Delete('N');  // Delete and tweak alignment.
	    if(!dont_delete) Delete('n');  // Delete and tweak alignment.
	    // Delete('n');  // Delete and tweak alignment.
#if 0
	    if(lastcma!=0){
		cmaGR=GRecombineCMSA(lastcma, cmsa);
		if(cmaGR != NULL){
		   if(use_gseq) sprintf(options,"-t1 -g -l%d ",limit);
		   else sprintf(options,"-t1 -l%d ",limit);
        	   map0=core_gibbs(&cmaGR);
		   if(map0 > map){
			map=map0; NilCMSA(cmsa); cmsa=cmaGR; Record( ); 
		   } else NilCMSA(cmaGR);
		}
	    }
#endif
	    if(InsertMSAHeap(cmsa,map,maH)==NULL){ NilCMSA(cmsa); cmsa=0; }
	    ConvergedMSAHeap(maH);
	    if(BestItemMSAheap(maH) != item){ // i.e., new optimum alignment.
		     item = BestItemMSAheap(maH);
		     // SeeMSAHeap(item,maH); 
	    }
	    lastcma=cmsa;
	    futile_attempts=0;
	} else {
		if(cmsa!=NULL) NilCMSA(cmsa); cmsa=NULL; 
		cycle--; futile_attempts++; 
		if(futile_attempts > 20) print_error("failed to find a motif");
		else if(futile_attempts % 10 == 0) { 
		   if(!fix){
			minblk--; maxblk--;
			minblk = MAXIMUM(Int4,minblk,1); 
			if(maxblk < minblk) maxblk = minblk;
		   }
		}
	}
        fprintf(stderr,"\ttime: %d seconds (%0.2f minutes)\n",
                time(NULL)-time0,(float)(time(NULL)-time0)/60.0);
     }
     fprintf(stderr,"\tTotal time: %d seconds (%0.2f minutes)\n",
                time(NULL)-time0b,(float)(time(NULL)-time0b)/60.0);
#if 1 // DO SOME POST PROCESSING HERE...
     NumAln=nMSAHeap(maH);
     NEW(list,NumAln+3,cma_typ);
     NEW(lpr,NumAln+3,double);
     for(aln=1; nMSAHeap(maH) > 0; aln++){ 
	assert((cmsa=DelMinMSAHeap(&map, maH)) != NULL);
fprintf(stderr,"******************** refining alignment %d *********************\n",aln);
	for(t=1; t<=nBlksCMSA(cmsa); t++) avelen+=LengthCMSA(t,cmsa);
	avelen/=(double)nBlksCMSA(cmsa);
#if 1
	if(use_gseq) sprintf(options,"-t1 -g -l%d ",limit);
	else sprintf(options,"-t1 -l%d ",limit);
        // map=core_gibbs(&cmsa,'S',300);
        map=core_gibbs(&cmsa);
	Record( );
#else
	if(use_gseq) sprintf(options,"-t1 -g -l%d ",limit);
	else sprintf(options,"-t1 -l%d ",limit);
	Delete();	// map?  did this above already...!!!
	Record( );
#endif
#if 1
        if(align_mode == 2){
	  if(Split(8)){
		Record( );
		fprintf(stderr,"Splitting block improved map\n"); 
	  }
	  if(Delete()){
		Record( );
		fprintf(stderr,"Deleting block improved map\n"); 
	  }
	  // if(avelen <= 15){ Delete(); Fuse(50); Split(8); }
	  // else { Delete(); Split(8); Fuse(50); }
        } else if(align_mode == 3){
	  if(SFR(8,40)){	// slide pieces to the right.
		Record( );
		fprintf(stderr,"Split Fuse to right improved map\n"); 
	  }
	  if(SFL(8,40)){	// slide pieces to the right.
		Record( );
		fprintf(stderr,"Split Fuse to right improved map\n"); 
	  }
        }
#endif
#if 0
	map0=ReAlignBestCMSA(cmsa);
	if(map0 > map){ map=map0; Record( );}
	else InitMAPCMSA(cmsa);
#endif
	list[aln]=cmsa; lpr[aln]=map;
     }
     for(aln=1; aln <= NumAln; aln++){
	assert(InsertMSAHeap(list[aln],lpr[aln],maH)!=NULL);
     }
     free(list); free(lpr);
#endif
}

void	gsm_typ::Recombine()
// GENETIC RECOMBINATION OF ALIGNMENT POPULATION.
{
   Int4		min;
   cma_typ	ma2;
   double	map2;

     if(maH==NULL) print_error("gsm_typ::Recombine() input error");
     if(use_gseq) sprintf(options,"-t1 -g "); else sprintf(options,"-t1 ");
     fprintf(stderr,"******** breed ********\n");
     map=Breed( );
     while(nMSAHeap(maH) > 3){ 
	if((ma2 = DelMaxMSAHeap(&map, maH)) != NULL) NilCMSA(ma2); 
     }
     assert(bestmsa != NULL); min = nBlksCMSA(bestmsa); 
     if(!fix){ minblk = MAXIMUM(Int4,min-under,2); maxblk = min+over;}
     fprintf(stderr,"--> End breed ");
     double runtime=difftime(time(NULL),time1);
     fprintf(stderr," time thus far: %0.1f seconds (%0.2f minutes)\n",runtime,runtime/60.0);
     cmsa = DelMinMSAHeap(&map2, maH); NilMSAHeap(maH); maH=NULL;
}

void gsm_typ::optimize() 
// FINAL OPTIMIZATION OF BEST ALIGNMENT.
{
   Int4		j;
   cma_typ	ma2,ma3;
   double	map2;

   if(align_mode == 2){
	map = bestmap;
	Split(8); while(Fuse(50)){ if(!Split(8)) break; }
	do {
	    if(SFL(8, 50)) FSR(60); else if(!FSR(60)) break; 
	    if(SFR(8, 50)) FSL(60); else if(!FSL(60)) break; 
	} while(TRUE);
   }
   if(use_gseq) strcpy(options,"-t1 -g "); else strcpy(options,"-t1 ");
   ma2 = CopyCMSA(cmsa);
   fprintf(stderr,"starting temperature = %.2f K\n", 300.);
   for(j=1; j<=10; j++){
     map2=core_gibbs(&ma2,'S',300);
     fprintf(stderr,"map = %.2f\n", map2);
     if((ma3 = RecombineCMSA(cmsa, ma2)) != NULL){
   	map2=core_gibbs(&ma3,'S',150);
        NilCMSA(ma2); ma2=ma3; ma3=NULL;
     } 
     if(map2 > bestmap){ 
	NilCMSA(cmsa); cmsa=ma2; bestmap=map2; PutAlnCMSA(name,cmsa,Guide); 
	ma2=NULL; break;
     } 
   }
   if(ma2 != NULL) NilCMSA(ma2);

#if 0
     map = core_gibbs(&cmsa,'S',150);
     if(map > bestmap){ bestmap=map; PutAlnCMSA(name,cmsa,Guide); }
#endif

}

//==================== TEST NEW ROUTINES =====================

BooLean	gsm_typ::SFR(Int4 minlength, Int4 maxlengs)
//   Split & Fuse operation (to the right): 2 -> 3 -> 2 blocks.
//
// --[____|\\\\]----[____]--   ====> --[____]----[/////|____]--  
//        :---->                                 
//
// decide whether to fuse first & split second or vice versa
{
   cma_typ	maS,maSF;
   Int4		t,numblks = nBlksCMSA(cmsa);
   BooLean	improved=FALSE;
#if 1
static Int4	ncalls=0;
#endif

   if(numblks > 1){
     fprintf(stderr,"%d blocks; try sliding halves of blocks right...\n",numblks);
     for(t=1; t < numblks; t++){  // requires at least two blocks 
	// split & then fuse (also try other way?).
	if((LengthCMSA(t,cmsa)+LengthCMSA(t+1,cmsa)) > (minlength+maxlengs))
		continue;
	if((maS=SplitBlkCMSA(t, minlength, cmsa)) != NULL){
	   if((maSF=FuseBlksCMSA(t+1, maxlengs, maS)) != NULL){
#if 0
		PutAlnCMSA("junk0", cmsa,NULL);
		PutAlnCMSA("junk1", maS,NULL);
		PutAlnCMSA("junk2", maSF,NULL);
#endif
		if(gibbs(maSF)){ 
			improved = TRUE;
#if 0
			PutAlnCMSA("junk3", cmsa,NULL);
			exit(1);
#endif
		}
	   }
	   if(gibbs(maS)) improved = TRUE;
	}
        numblks=nBlksCMSA(cmsa);
     }
   } // for completeness do another split operation at the end.
#if 0
   if((maS=SplitBlkCMSA(numblks, minlength, cmsa)) != NULL){
	   if(gibbs(maS)) improved = TRUE;
   }
#endif
   return improved;
}

BooLean	gsm_typ::SFL(Int4 minlength, Int4 maxlengs)
//   Split & Fuse operation (to the right): 2 -> 3 -> 2 blocks.
// TRY Splitting BLOCKS with >= minlength columns and fusing with next block
//
// --[____]----[\\\\|____]--   ====> --[____|/////]----[____]--  == !slide_right
//             <----:                                 
{
   cma_typ	maS,maSF;
   Int4		t,numblks = nBlksCMSA(cmsa);
   BooLean	improved=FALSE;

   if(numblks > 1){
     fprintf(stderr,"Try sliding halves of blocks left...\n");
     for(t=numblks; t > 1; t--){
	if((LengthCMSA(t,cmsa)+LengthCMSA(t+1,cmsa)) > (minlength+maxlengs))
		continue;
	if((maS=SplitBlkCMSA(t,minlength,cmsa)) != NULL){
	   if((maSF=FuseBlksCMSA(t-1, maxlengs, maS)) != NULL){
		if(gibbs(maSF)) improved = TRUE;
	   }
	   if(gibbs(maS)) improved = TRUE;
	} 
     }
   }	// for completeness do another split operation at the end.
#if 0
   if((maS=SplitBlkCMSA(1, minlength, cmsa)) != NULL){
	   if(gibbs(maS)) improved = TRUE;
   }
#endif
   return improved;
}

BooLean	gsm_typ::FSR(Int4 maxlength)
// FuseSplit operation:   2 -> 1 -> 2 blocks (net = 0)
// ***************************************************************************
//                                  F(t)           
//   msa   --[____]----[\\\\\\]--   ===>  --[____|\\\\\\]---  maF
//             t          t+1                    t        
//                                                 |  
//                                                 V  S(t)
//                                     
//                                        --[____|\\]--[\\\]-  maFS
//                                               t      t+1     
// ***************************************************************************
{
   cma_typ	maF,maFS;
   Int4		t,len1,len2,numblks = nBlksCMSA(cmsa);
   BooLean	improved=FALSE;

   if(numblks < 2) return FALSE;
   fprintf(stderr,"Try FSR operation ...\n");
   for(t=1; t < numblks; t++){
      len1 = LengthCMSA(t,cmsa); len2 = LengthCMSA(t+1,cmsa);
      if((len1+len2) <= maxlength){
	if((maF=FuseBlksCMSA(t,len1+len2,cmsa)) != NULL){
	   if((maFS= SplitBlkCMSA(t,len1,maF)) != NULL){
		if(gibbs(maFS)) improved = TRUE; 
	   } else print_error("gsm_typ::FSR( ): This should not happen");
	   if(gibbs(maF)) improved = TRUE;
	} else print_error("gsm_typ::FSR( ): This should not happen");
        numblks=nBlksCMSA(cmsa);
      }
   }
   return improved;
}

BooLean	gsm_typ::FSL(Int4 maxlength)
// Same as FSL only going to the left. 2 -> 3 -> 4 -> 3
{
   cma_typ	maF,maFS;
   Int4		t,len1,len2,numblks = nBlksCMSA(cmsa);
   BooLean	improved=FALSE;

   if(numblks < 2) return FALSE;
   fprintf(stderr,"Try FSL operation ...\n");
   for(t=numblks-1; t >= 1; t--){
      len1 = LengthCMSA(t,cmsa); len2 = LengthCMSA(t+1,cmsa);
      if((len1+len2) <= maxlength){
	if((maF=FuseBlksCMSA(t,len1+len2,cmsa)) != NULL){
	   if((maFS= SplitBlkCMSA(t,len1,maF)) != NULL){
		if(gibbs(maFS)) improved = TRUE; 
	   } else print_error("gsm_typ::FSL( ): This should not happen");
	   if(gibbs(maF)) improved = TRUE;
	} else print_error("gsm_typ::FSL( ): This should not happen");
      }
   }
   return improved;
}

BooLean	gsm_typ::SSF(Int4 t, Int4 minlength)
// SplitSplitFuse operation:   2 -> 3 -> 4 -> 3
// ***************************************************************************
//                                     S(t)           
//   msa   --[__|//]----[\\\|____]--   ===>  --[__|//]--[\\\]--[____]--  maS
//             t-1          t                    t-1     t      t+1
//                                                       |  
//                                                       V  S(t-1)
//                                     F(t)
//   maSSF --[__]--[//|\\]--[____]--   <===  --[__]-[//]-[\\\]-[____]--  maSS
//             t-1    t       t+1               t-1  t    t+1    t+2  
// ***************************************************************************
{
      cma_typ	maS,maSS,maSSF;
      Int4	len1,len2,numblks = nBlksCMSA(cmsa);
      BooLean	improved=FALSE;

      if(numblks < 2 || t < 2 || t > numblks) return FALSE;
      fprintf(stderr,"Try SSF operation on block %d...\n",t);
      len1 = LengthCMSA(t-1,cmsa); len2 = LengthCMSA(t,cmsa);
      if(len1 >= minlength && len2  >= minlength){
	if((maS = SplitBlkCMSA(t,minlength,cmsa)) != NULL){
	   if((maSS = SplitBlkCMSA(t-1,minlength,maS)) != NULL){
		maSSF=FuseBlksCMSA(t,len1+len2,maSS);
	   	if(maSSF==NULL) print_error("gsm_typ::SSF( ): This should not happen");
		if(gibbs(maSSF)) improved = TRUE; 
	   } else print_error("gsm_typ::SSF( ): This should not happen");
	   // gibbs(maSS); 
	   NilCMSA(maSS); 
	} else print_error("gsm_typ::SSF( ): This should not happen");
	// gibbs(maS);
	NilCMSA(maS); 
      }
      numblks=nBlksCMSA(cmsa);
      return improved;
}

BooLean	gsm_typ::SFF(Int4 t, Int4 minlength)
// SplitFuseFuse operation:   3 -> 4 -> 3 -> 2
// ***************************************************************************
//                                      S(t)           
//   msa   --[___]---[\\|//]---[___]--  ===>  --[___]--[\\]-[//]--[___]--  maS
//            t-1       t       t+1              t-1     t   t+1   t+2
//                                                       |  
//                                                       V  F(t-1)
//                                      F(t)
//  msaSFF --[___|\\]-------[//|___]--  <===  --[___|\\]----[//]--[___]--  maSF
//             t-1             t                   t-1       t     t+1     
// ***************************************************************************
{
   cma_typ	maS,maSF,maSFF;
   Int4		len1,len2,len3,numblks=nBlksCMSA(cmsa);
   BooLean	improved=FALSE;

   if(numblks < 3 || t < 2 || t >= numblks) return FALSE;
   fprintf(stderr,"Try SFF operation on block %d ...\n",t);
   len1 = LengthCMSA(t-1,cmsa); len2 = LengthCMSA(t,cmsa); len3 = LengthCMSA(t+1,cmsa);
   if(len2 >= minlength){
	if((maS = SplitBlkCMSA(t,minlength,cmsa)) != NULL){
	   if((maSF = FuseBlksCMSA(t-1,len1+len2,maS)) != NULL){
		maSFF=FuseBlksCMSA(t,len2+len3,maSF);
	   	if(maSFF==NULL) print_error("gsm_typ::SFF( ): This should not happen");
		if(gibbs(maSFF)){ // implies that maSFF == NULL;
			improved = TRUE;
		} 
	   } else print_error("gsm_typ::SFF( ): This should not happen");
	   // gibbs(maSF);
	   NilCMSA(maSF); 
	} else print_error("gsm_typ::SFF( ): This should not happen");
	// gibbs(maS);
	NilCMSA(maS); 
   }
   return improved;
}

cma_typ	gsm_typ::optbreed(cma_typ cma0)
{ 
  char		Str[300];
  double	Temperature=200;
  gss_typ	*gssX;

  if(cma0 == 0){ sprintf(Str,"%s.cma",name); cma0=ReadCMSA2(Str,AB); }
  SetPenaltyCMSA(gapopen,gapextend,cma0);
  sprintf(name,"%s%d",Argv[1],Run);  
  FILE *fptr=open_file(name,"","w");
  PutSeqSetEs(fptr,TrueDataCMSA(cma0));
  fclose(fptr);
  assert(cmsa==0); cmsa = cma0; map = RelMapCMSA(cmsa);
  // gssX=gssCMSA(cmsa); make_binomials(gssX->FakeSqSet());
  Record( ); cmsa=0; // sets bestmsa to cmsa;
  if(use_gseq) sprintf(options,"-t1 -g -l%d ",limit);
  else sprintf(options,"-t1 -l%d ",limit);
  assert((maH=MkMSAHeap(mhpsz)) != 0);
  if(InsertMSAHeap(cma0,bestmap,maH)==0) print_error("error in gsm_typ::optbreed");
  fprintf(stderr,"starting Temperature = %.2f K (%.2f degrees F)\n",
                Temperature,((9.*Temperature)/5.)-459);
  for(Int4 n=1; n <= maxcycle; n++){
	cmsa = CopyCMSA(bestmsa);
	map=core_gibbs(&cmsa,'S',temperature);
	if(map > bestmap){ Record( ); }
	if(InsertMSAHeap(cmsa,map,maH)==NULL){ 
   		gssX=gssCMSA(cmsa);     // gssX not owned by cma.
		NilCMSA(cmsa); cmsa=0;
   		gssX->~gss_typ(); 
	} else cmsa=0;
	if(ConvergedMSAHeap(maH) < breedVar) break;
  }
  map=Breed( ); 
  cmsa=DelMinMSAHeap(&map, maH); NilMSAHeap(maH);
  if(map > bestmap){ Record( ); }
  // cmsa = RecombineCMSA(cma0,cmsa);
  free_input_data(); // fixes memory leak: AFN 5/23/01
  return cmsa;	// need to free this up in calling environment.
}

#if 0
cma_typ	gsm_typ::optbreed(cma_typ msa0)
{ 
	double  optmap,hotmap,Map,oldmap,newmap,Temperature=300,tempmap,map2;
	Int4	r,n,i,*Num;
	cma_typ msa,oldmsa,tempmsa;
	char	Str[300],Str2[300],Str0[300];
	Int4	Aveblk,Avecol,t,numblks,num_aln=mhpsz;

  if(msa0 == NULL){ sprintf(Str,"%s.cma",name); msa0=ReadCMSA2(Str,AB); }
  SetPenaltyCMSA(gapopen,gapextend,msa0);
  sprintf(name,"%s%d",Argv[1],Run);  
  if(cmsa != msa0) cmsa = CopyCMSA(msa0);
  map = RelMapCMSA(cmsa);
  oldmsa = cmsa;
  optmap = RelMapCMSA(oldmsa);
  Aveblk = nBlksCMSA(oldmsa);
  Avecol = (NumColumnsCMSA(oldmsa)/Aveblk) + 1;

  st_type S=SitesCMSA(oldmsa);
  NEW(Num,nBlksCMSA(oldmsa)+2,Int4);
  for(i=1; i<= nBlksCMSA(oldmsa); i++) Num[i] = SiteLen(i,S);
  if(use_gseq) strcpy(options,"-t1 -g -T1 ");
  else strcpy(options,"-t1 -T1 ");
  free(Num);

  maH=MkMSAHeap(mhpsz); assert(maH != NULL);
  if(InsertMSAHeap(cmsa,map,maH)==NULL)
		print_error("input error in gsm_typ::optbreed");
  sprintf(Str0,"%s.best",Argv[1]);
  fprintf(stderr,"starting Temperature = %.2f K (%.2f degrees F)\n",
                Temperature,((9.*Temperature)/5.)-459);
  for(n=1; n <= num_aln; n++){
	msa = CopyCMSA(msa0);
	map = hotmap = 0.0;
	core_gibbs(&msa,'D',150);
	tempmap = map; tempmsa= cmsa; map = hotmap; cmsa= msa;
	switch(align_mode){
	  case 0:	// fast & dirty mode do simulated annealing only.
#if 0
std::cerr << "SAMPLE HOTGIBBS\n\n";
	     // map=hotmap=HotGibbs(options,&msa,Temperature);
	     map=hotmap=core_gibbs(&msa,'H',Temperature);
#endif
std::cerr << "SAMPLE SIM_ANEAL_GIBBS\n\n";
	     core_gibbs(&msa,'D',150);
	     map=hotmap=RelMapCMSA(msa);
	   break;
	  case 1:	// fast mode: only split & fuse once.
		Split(8); Fuse(50); break;
		// Fuse(50); break;
	  case 2: Split(8); Fuse(50); Split(8); break;
	  case 3: Split(8); while(Fuse(50)){ if(!Split(8)) break; } break;
	  case 4:		// 
		SFR(8, 50);	// slide pieces to the right.
		FSL(60);	// Nudge going back to the left.
		do {
		    if(SFL(8, 50)) FSR(60); else break; 
		    // else if(!FSR(60)) break; 
		    if(SFR(8, 50)) FSL(60); else break;
		    // else if(!FSL(60)) break;
		} while(TRUE);
	   break;
	  case 5:
		Split(8); while(Fuse(50)){ if(!Split(8)) break; }
		do {
		    if(SFL(8, 50)) FSR(60); else if(!FSR(60)) break; 
		    if(SFR(8, 50)) FSL(60); else if(!FSL(60)) break; 
		} while(TRUE);
	   break;
	  case 6: 	// = abhfbestB.msa
	    for(t=2,numblks=nBlksCMSA(cmsa); t <= numblks; t++){
		map = 0.0;  // forces SSF to accept altered alignment
	  	if(SSF(t,8)) { map = 0.0; SFF(t,8); }
		Split(8); while(Fuse(50)) if(!Split(8)) break;
		numblks=nBlksCMSA(cmsa); 
	    } break;
	  case 7: case 8: case 9:
	  default: print_error("gsm_typ::align() input error"); break;
        } 
	hotmap = map; msa = cmsa; map = tempmap; cmsa = tempmsa;
        Map = RelMapCMSA(msa);
	if(Map > optmap){ 
		optmap=Map; PutAlnCMSA(Str0,msa,NULL); 
	}
#if 1
	fprintf(stderr,"%d: hot map = %f; Map = %f\n",n,hotmap,Map);
#endif
	if(InsertMSAHeap(msa,Map,maH)==NULL){ NilCMSA(msa); }
	ConvergedMSAHeap(maH);
  }
  sprintf(Str,"%s.best",NameCMSA(msa0));
  Map=Breed( );

  msa = DelMinMSAHeap(&map2, maH); NilMSAHeap(maH);
  cmsa = RecombineCMSA(msa0, msa);
  if(cmsa != NULL) { Map = RelMapCMSA(cmsa); NilCMSA(msa); msa=NULL; }
  else { cmsa = msa; msa=NULL; }
#if 0	// PROBABLY NEED TO DEALLOCATE gss_typ!!!
   s_typ *gssX=gssCMSA(cmsa_in);     // not owned by cmsa_in.
   gssX->~gss_typ();
#endif
  return cmsa;
}
#endif

BooLean	gsm_typ::create_align()
// create an alignment and store it in cmsa.
{
	cmsa = NULL; 
	gibbs(NULL);
	if(cmsa == NULL) return FALSE;
#if 0
	else return TRUE;
#else
	BooLean	S,F,D;
	Int4	t;
	double	avelen=0.0;
	for(t=1; t<=nBlksCMSA(cmsa); t++) avelen+=LengthCMSA(t,cmsa);
	avelen/=(double)nBlksCMSA(cmsa);
	switch(align_mode){
	  case 0: break; // fast & dirty mode don't do anything
	  case 1: Split(8); break;
	  // case 1: Delete(); Add(7); Split(8); Fuse(50); break;
	  case 2: Add(5); Split(8); break;
	  // case 2: Split(8); Delete(); break;
	  case 3: Delete(); Split(8); break;
	  case 4: Split(8); Fuse(50); Split(8); break;
	  case 5: Split(8); while(Fuse(50) || Split(8)); break;
	  case 6: Split(8); while(Fuse(50) || Delete() || Split(8)); 
	   break;
	  case 7: 
	     if(avelen <= 9){
		do{ D=Delete(); F=Fuse(50); S=Split(8); } while(S || F || D);
	     } else if(avelen >= 16){
		do{ S=Split(8); D=Delete(); F=Fuse(50); } while(S || F || D);
	     } else do{ F=Fuse(50); S=Split(8); D=Delete(); } while(S || F || D);
	   break;
	  case 8: case 9:
	  default: print_error("gsm_typ::align() input error"); break;
        }
	// gibbs(NULL); // tweak cmsa
	return TRUE;
#endif
}

