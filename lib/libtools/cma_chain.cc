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

#include "cma_chain.h"

#define MAX_CMA_CHAIN_SETS	100

set_typ  *RunClustCMSA(FILE *fp,Int4 *num_subsets, Int4 **subset_increase,
	Int4 score_cut, cma_typ cmsa)
{
	ds_type sets;
	Int4	i,j,set1,score;
	time_t	time1=time(NULL); 

	assert(nBlksCMSA(cmsa) ==1);
	ss_type data = TrueDataCMSA(cmsa);
	ss_type fakedata = DataCMSA(cmsa);
	a_type	A=SeqSetA(data);
	Int4 N = NSeqsSeqSet(data);

	sets = DSets(N);
	Int4 **ScoreMTRX;
	NEWP(ScoreMTRX,N+3,Int4);
	for(i=1; i <= N; i++) NEW(ScoreMTRX[i],N+3,Int4);
	dh_type dH = dheap((N*N/2)+N+3,4);
	Int4 hpsz=0;

	for(i=1; i < N; i++){
	  for(j=i; j <= N; j++){
	    // very fast pseudo scores as clustering heuristic...
	    score = PseudoAlnScoreCMSA(i,j,cmsa);
	    // fprintf(stderr,"%d vs %d: %d\n",i,j,score);
	    ScoreMTRX[i][j]=ScoreMTRX[j][i]=score;
	    hpsz++;
	    if(score > 0) insrtHeap(hpsz,(keytyp)-score,dH);
	  } fprintf(stderr,"%d\n",i);
	}
	if(fp){
		for(i=1; i <= N; i++){
		   for(j=1; j <= N; j++){
			fprintf(fp," %5d",ScoreMTRX[i][j]);
		   } fprintf(fp,"\n");
		}
	}
	double	runtime=difftime(time(NULL),time1);
	fprintf(stderr,"\ttime: %0.1f seconds (%0.2f minutes)\n",runtime,runtime/60.0);
	score = (Int4) -minkeyHeap(dH);
	h_type HG=Histogram("cluster size as a function of cutoff score",0,score,1.0);
	h_type HG2=Histogram("change in cluster size as a function of cutoff score",0,score,5.0);
	fprintf(stderr,"create histogram\n");
	Int4 lastscore=0,lastsize=1;
#if 0	// Faster routine:
	// sort edges by score and link sets using sorted edges rather then 
	// iterating through the score matrix so many times...
	// need two arrays of unsigned shorts for storing i and j nodes.
#endif
	Int4 max_increase=0,cumulative=0,num_stored=0;
	set_typ *category=0;
	Int4	size_set;
#if 1
	set_typ	*cluster;
	NEW(cluster,N+2,set_typ);
	for(i=1; i <= N; i++) { cluster[i] = MakeSet(N+1); AddSet(i,cluster[i]); }
	set_typ LastQuerySet = MakeSet(N+1); AddSet(1,LastQuerySet);
#endif
	Int4 *set_increase;
	while(NSetsDSets(sets) > 1 && !emptyHeap(dH) &&
				(score = (Int4) -minkeyHeap(dH)) > 0){
	   delminHeap(dH);
	   while(score == lastscore){ 	// find next score cutoff
	    if(emptyHeap(dH)) break;
	    score = (Int4) -minkeyHeap(dH); delminHeap(dH); 
	   } // fprintf(stderr,"score = %d\n",score);
	   for(i=N; i > 1; i--){		// cluster at current score cutoff
		register Int4 Set1,J,Score,*mtrx;
		register ds_type Sets;
		register Int4 Set0,Set2;
		register set_typ TmpSet;

		Sets = sets; Set1=findDSets(i,sets);
		Score=score; mtrx=ScoreMTRX[i]; 
		for(J = i-1; J > 0; J--){
		  if(mtrx[J] == Score){	// higher scores are already merged.
		    Set2 = findDSets(J,Sets);
		    if(Set1!=Set2){
		      UnionSet(cluster[Set1], cluster[Set2]); // sets set1 to set1 U set2.
		      Set0=linkDSets(Set1,Set2,Sets); 
		      if(Set0 == Set2){ 
			TmpSet=cluster[Set2]; cluster[Set2]=cluster[Set1];  
			cluster[Set1]=TmpSet; Set1=Set2;
		      } else Set1=Set0;
		    }
		  }
		}
	   }
	   set1 = findDSets(1,sets);
	   size_set = CardSet(cluster[set1]);
	   Int4 increase_in_size = size_set-lastsize;
	   if(size_set > 1){
	    if(increase_in_size > 0){
	     if(num_stored == 0){	// create arrays for returning sets and set increases.
	       NEW(category,score+3,set_typ);
	       NEW(set_increase,score+3,Int4);
	       *subset_increase = set_increase;
	       num_stored=score;
	     }
	     fprintf(stderr,"set1=%d; score = %d; size_set = %d; increase = %d\n",
					set1,score,size_set,increase_in_size);
	     if(max_increase < increase_in_size) max_increase = increase_in_size;
	     set_increase[score] = increase_in_size;
	     category[score]=MakeSet(N+1);
	     CopySet(category[score],cluster[set1]); // return new query set.
	     IntersectNotSet(category[score],LastQuerySet); // == new sequences...
	     CopySet(LastQuerySet,cluster[set1]); // This should work.
	     // UnionSet(LastQuerySet,cluster[set1]); // This should work.
	     // CopySet(category[score],LastQuerySet); // return new query set.
	    }
	    for(Int4 tmpscore=lastscore-1; tmpscore > score; tmpscore--){
	   	IncdMHist(tmpscore,lastsize,HG);
	    } IncdMHist(score,size_set,HG);
	    IncdMHist(score,size_set-lastsize,HG2);
	   }
	   lastscore=score; lastsize=size_set;
	} PutHist(stdout,60,HG); NilHist(HG);
	PutHist(stdout,60,HG2); NilHist(HG2);
	for(i=1; i <= N; i++){
	   free(ScoreMTRX[i]); NilSet(cluster[i]);
	} free(cluster); free(ScoreMTRX); NilSet(LastQuerySet);
	NilDSets(sets); Nildheap(dH);
	*num_subsets = num_stored;

	// set zero set to cluster at input cutoff score.
	for(i=score_cut; i <= num_stored; i++){
	   if(category[i]){ // i.e., set exists.
		category[0] = category[i]; break;
	   }
	} if(num_stored > 0 && category[0] == 0) category[0] = category[i];

	runtime=difftime(time(NULL),time1);
	fprintf(stderr,"\ttime: %0.1f seconds (%0.2f minutes)\n",runtime,runtime/60.0);
	// return telescoping sets of related sequences.
	return category;
}

