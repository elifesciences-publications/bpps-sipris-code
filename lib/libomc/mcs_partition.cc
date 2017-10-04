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

void	mcs_typ::PartitionByInputSetCMSA( )
// Partition the alignment cma based on sequence scores against the seed cma files
{
	// 0. Create the Main alignment by adding in Random sequences.
	assert(passed_in_cma); assert(TrueMainCMA==passed_in_cma);
	Int4 i,J;
	for(i=1; i <= NumDisplayCMA; i++){ InitSet[i] = passed_in_sets[i]; }
        if(i == num_passed_in_sets) InitSet[i] = passed_in_sets[i];  
	else InitSet[i] = CopySet(RandomSet); // random set not passed in.
	ComputeSeedCMAs();
}

void	mcs_typ::RandomlyPartitionSets( )
{
	ComputeSeedCMAs();
	Int4	sq,i,M=SizeTrueMain,N=SizeMain,R=Hpt->NumSets();
	for(i=1; i <= Hpt->NumSets(); i++){ InitSet[i] = MakeSet(N+1); ClearSet(InitSet[i]); }
	InitSet[i] = MakeSet(N+1); ClearSet(InitSet[i]);  // One more For Random Sequence set
	InitSet[0]=InitSet[Hpt->NumSets()]; 	// Let 0 point to Random sequence set as well
	for(sq=1; sq <= M; sq++){
	  i=random_integer(Hpt->NumSets()-1) + 1;	// don't put these into random set.
	  AddSet(sq,InitSet[i]); 	// Set representation.
	} for(sq=M+1; sq <= N; sq++) AddSet(sq,InitSet[0]); 	// These go into the random Set.
	InitSet[0]=0;	// Don't need this pointer anymore.
}

void	mcs_typ::PartitionBySeedAlnCMSA(Int4 NumSeedAln, cma_typ *seed_cma, char *TypeOfSet)
// Create an initial partitioning of the aligned input sequences (cma) based on seed alignment consensus scores.
{
	// //****** 2a. If only one BPPS category, then call the following. ******
	// if(NumSeedAln==1 && NumBPPS == 1){ PartitionSingleSeedAlnCMSA(NumSeedAln, seed_cma,TypeOfSet); return; }
        Int4    i,j,n,s,I,J,N=NumSeqsCMSA(MainCMA);
	FILE	*ofp=0;
	//****** 2b. Add a consensus sequence to each of the input seed alignments. ******
	ComputeSeedCMAs();
	//****** 3a. Compute the maximum score between any two seed alignment consensus sequences. ******
	if(ofp){
	  fprintf(stdout,"   ");
	   for(j=1; j <= NumSeedAln; j++) fprintf(stdout," %5d",j); fprintf(stdout,"\n");
	}
	Int4 *MaxScoreNonID; NEW(MaxScoreNonID,NumSeedAln+3,Int4);
	for(i=1; i <= NumSeedAln; i++){
	  e_type Ei=FakeSeqCMSA(1,SeedCMA[i]);
	  if(ofp) fprintf(stdout,"%2d: ",i);
	  MaxScoreNonID[i]=INT4_MIN;
	  assert(MaxScoreNonID[i] < 0);
	  for(j=1; j <= NumSeedAln; j++){
	     e_type Ej=FakeSeqCMSA(1,SeedCMA[j]);
	     assert(LenSeq(Ei) == LenSeq(Ej));
	     Int4 score=0;
	     for(s=1; s <= LenSeq(Ej); s++){
		unsigned char ri=ResSeq(s,Ei);
		unsigned char rj=ResSeq(s,Ej);
		score += valAlphaR(ri,rj,AB);
	     } if(ofp) fprintf(stdout,"%5d ",score);
	     if(i != j && TypeOfSet[j] != '?' && score > MaxScoreNonID[i]) MaxScoreNonID[i]=score;
	  } if(ofp) fprintf(stdout," %s%c\n",NameCMSA(SeedCMA[i]),TypeOfSet[i]);
	}
	if(ofp){
	  fprintf(stdout,"---");
	  for(j=1; j <= NumSeedAln; j++) fprintf(stdout,"------",j); fprintf(stdout,"\n");
	  fprintf(stdout,"   ");
	  for(j=1; j <= NumSeedAln; j++) fprintf(stdout," %5d",MaxScoreNonID[j]); fprintf(stdout,"\n");
	}

	//****** 4. Find the Best Partition for each of the input sequences. ******
	Int2	*Partition; NEW(Partition, N+3, Int2);		
	h_type HG=Histogram("partitions", 0,NumSeedAln,1.0);
	h_type sHG=Histogram("best scores", -100,5000,25.0);
	// if(key_seq < 1 || key_seq > N) print_error("PartitionBySeedAlnCMSA( ) input error");
	Int4	M=NumSeqsCMSA(TrueMainCMA);
	for(J=1; J <= M; J++){
	   Int4 BestScore = INT4_MIN;
	   Int2	BestPartition=0;
	   for(n=1; n <= NumSeedAln; n++){
		Int4 Score=PseudoAlnScoreTwoCMSA(1,SeedCMA[n],J,TrueMainCMA);
		if((TypeOfSet[n] == '?')){ 	// Compare score for all Misc sets...
		    if(Score > BestScore){ BestScore=Score; BestPartition=n; }
		} else if(Score >= MaxScoreNonID[n] && Score > BestScore)
		    { BestScore=Score; BestPartition=n; }  // For others, only if > MaxScore...
	   }
	   if(BestPartition == 0){
		print_error("Hyperpartition lacks a miscellaneous set (indicator: '?')");
	   	// assert(BestPartition != 0);
	   }
	   Partition[J]=BestPartition;
	   IncdHist((double)BestPartition,HG); IncdHist((double)BestScore,sHG);
	}
        if(efp){ PutHist(stdout,50,HG); PutHist(stdout,50,sHG); }
	NilHist(sHG); NilHist(HG); 

	//****** 5. Create initialization sets for each partition. ******
	for(i=1; i <= NumSeedAln; i++){ InitSet[i] = MakeSet(N+1); ClearSet(InitSet[i]); }
	InitSet[i] = MakeSet(N+1); ClearSet(InitSet[i]);  // One more For Random Sequence set
	InitSet[0]=InitSet[i]; 	// Let 0 point to Random sequence set as well
	HG=Histogram("partitions", 0,NumSeedAln,1.0);
	for(n=1; n <= NumSeedAln; n++){
	   Int4 Cnt=0;
	   for(J=1; J <= M; J++){
		if(Partition[J] == n){
		   AddSet(J,InitSet[n]); 	// Set representation.
		   IncdHist((double)n,HG); Cnt++; 
		}
	   }
	   // assert(Cnt > 0);	// if none get put in will not work...
	   if(0 && Cnt <= 0){
		fprintf(stderr,"Partition initialization failed (seed alignment: '%s')\n",
			NameCMSA(seed_cma[n]));
		fprintf(stderr,"input alignment lacks sequences matching this seed alignment\n");
		print_error(" Try modifying the hyperpartition and/or seed alignment.");
	   }
#if 0
	   RenameCMSA(NameCMSA(seed_cma[n]),cma);
#endif
	} free(Partition);
	for(J=M+1; J <= N; J++) AddSet(J,InitSet[0]); 	// These go into the random Set.
        if(efp) PutHist(stdout,50,HG); NilHist(HG); free(MaxScoreNonID);
	InitSet[0]=0;	// Don't need this pointer anymore.
	return;
}

void	mcs_typ::ComputeSeedCMAs()
{
	Int4	i,j,s,score;  
	assert(NumDisplayCMA < MAX_NUM_ELMENTARY_SETS);
	for(Int4 i=1; i <= NumDisplayCMA; i++){ SeedCMA[i] = AddConsensusCMSA(DisplayCMA[i]); } 
	assert(SeedCMA[1]);		 // ComputeMinSeed2CsqScores()
	Int4	r1,rj,Length=LengthCMSA(1,SeedCMA[1]);
	char *TypeOfSet=Hpt->TypeOfSet();
	for(i=1; i <= NumDisplayCMA; i++){
		assert(SeedCMA[i]);
		assert(LengthCMSA(1,SeedCMA[i]) == Length);
         	e_type Ei=FakeSeqCMSA(1,SeedCMA[i]);	// assume the first is a consensus seq.
		assert(Ei);
		assert(i <= NumDisplayCMA);
         	MinSeed2CsqScore[i]=INT4_MAX;
		for(j=2; j <= NumSeqsCMSA(SeedCMA[i]); j++){
		   e_type Ej=FakeSeqCMSA(j,SeedCMA[i]);
		   for(score=0,s=1; s <= Length; s++){
			r1=ResidueCMSA(1,1,s,SeedCMA[i]);
			rj=ResidueCMSA(1,j,s,SeedCMA[i]);
			score += valAlphaR(r1,rj,AB);
		   }
		   if(TypeOfSet[i] != '?' && score < MinSeed2CsqScore[i]) MinSeed2CsqScore[i]=score;
         	}
        }
}

Int4	mcs_typ::Seq2SeedCsqScore(Int4 sq,Int4 set)
{
	Int4	s,r,R,score=0;  assert(SeedCMA[1]);
#if 0
	// char *TypeOfSet=Hpt->TypeOfSet();
        e_type Ei=FakeSeqCMSA(1,SeedCMA[set]);
	e_type Ej=FakeSeqCMSA(sq,MainCMA);
	if(LenSeq(Ej) != LenSeq(Ei)){
		s=SitePos(1,sq,1,SitesCMSA(MainCMA));
		Int4 r=ResidueCMSA(1,sq,1,MainCMA);
		fprintf(stderr,"%d; site=%d; r=%d; residue='%c'\n",sq,s,r,AlphaChar(r,AB));
		if(1 || LenSeq(Ei) != LengthCMSA(1,TrueMainCMA)) PutSeq(stderr,Ei,AB);
		if(1 || LenSeq(Ej) != LengthCMSA(1,TrueMainCMA)) PutSeq(stderr,Ej,AB);
		fprintf(stderr,"Input alignment formating error: aligned lengths inconsistent.\n");
		// assert(LenSeq(Ej) == LenSeq(Ei));
		print_error("Seq2SeedCsqScore() fatal error");
	}
	for(s=1; s <= LenSeq(Ej); s++) score += valAlphaR(ResSeq(s,Ei),ResSeq(s,Ej),AB);
#else
	assert(LengthCMSA(1,MainCMA) == LengthCMSA(1,SeedCMA[set]));
	for(s=1; s <= LengthCMSA(1,MainCMA); s++){
		r=ResidueCMSA(1,sq,s,MainCMA); R=ResidueCMSA(1,1,s,SeedCMA[set]);
		score += valAlphaR(r,R,AB);
#if 0	// DEBUG...
	   if(sq==2439){
		fprintf(stderr,"%d: site=%d; r=%d; residue='%c'\n",sq,s,r,AlphaChar(r,AB));
		fprintf(stderr,"seed csq: site=%d; R=%d; residue='%c'\n",s,R,AlphaChar(R,AB));
	   }
	}
	if(sq==2439){
        	e_type Ei=FakeSeqCMSA(1,SeedCMA[set]);
		e_type Ej=FakeSeqCMSA(sq,MainCMA);
		PutSeq(stderr,Ei,AB);
		PutSeq(stderr,Ej,AB);
		 exit(1);
	} return score;
#else
	} return score;
#endif
#endif
}

Int4    *mcs_typ::SortByScoreCMSA(FILE *fp, char mode, Int4 &first_best, double cut, cma_typ cma, Int4 set)
// return a list of sq numbers for
{
	Int4	*list,i,N,Item,*ItemToSq,NumSq=NumSeqsCMSA(cma);
	keytyp	first_best_key=0;
	h_type	HG=0;
	BooLean	found=FALSE;
	set_typ Set=InitSet[set];

	if(nBlksCMSA(cma) > 1){
	    // assert(nBlksCMSA(cma) > 1);
	    print_error("-scores option requires single blk cma file");
	}
	N=CardSet(Set);
	NEW(list,N+3,Int4); NEW(ItemToSq,N+4,Int4);
	mh_type mH=Mheap(N+2,3);
	for(i=1; i <= NumSq; i++){
	   if(MemberSet(i,Set)){	// is sequence in the seed set?
		Int4 score=i; 	// For Random set, order by sequence position = i.
		if(set != Hpt->NumSets()) score=Seq2SeedCsqScore(i,set);
		Item = InsertMheap((keytyp)score,mH);
		assert(Item > 0 && Item <= N+2); ItemToSq[Item]=i;
	   }
	}
        if(Verbose && ItemsInMheap(mH) > 0){
	   Int4 start= (Int4) floor(MinKeyMheap(mH));
	   Int4 end= (Int4) ceil(MaxKeyMheap(mH));
	   if(start==end) end++;
	   double inc=ceil((double)(end-start+1)); inc = inc/30.0;
	   char tmpstr[100];
	   sprintf(tmpstr,"ungapped csq scores for set %d (\"%s\")",set,Hpt->ElmntSetName(set));
	   HG=Histogram(tmpstr,start,end,inc);
        }
	for(i=0; ItemsInMheap(mH) > 0; ){
		Item=MinItemMheap(mH);
		keytyp key=MinKeyMheap(mH);
		if(HG) IncdHist((double)key,HG);
		assert(DelMinMheap(mH)==Item);
		i++; list[i]=ItemToSq[Item];
		if(!found && key >= cut){ found=TRUE; first_best=i; first_best_key=key; }
	}
	if(!found) first_best=N+1;
	if(HG){
	     if(fp) {
		fprintf(fp,"%s:\n num-to-label: %d (score=%.1f)\n",
			NameCMSA(cma),N-first_best+1,first_best_key);
		PutHist(fp,60,HG); 
	     } NilHist(HG); 
	} NilMheap(mH); free(ItemToSq);
	return list;
}

Int4	mcs_typ::GetElmntSets( )
// *.cma files; *.mma file; *.sma file.
{
	char	*name;
	Int4	i;

	cma_typ	*OldDisplayCMA=0;
	if(passed_in_sma){	// Make sure not to clobber passed_in_sma files...
		NEW(OldDisplayCMA,NumDisplayCMA + 3, cma_typ);
		for(i=0; i <= NumDisplayCMA; i++) OldDisplayCMA[i]=DisplayCMA[i]; 
		free(DisplayCMA);
	} else OldDisplayCMA=DisplayCMA;
	NEW(DisplayCMA,NumDisplayCMA + 3, cma_typ);
	// Hpt->NumSets()=Hpt->NumSets(); NumBPPS=Hpt->NumBPPS(); // moved further up...
	for(i=1; i <= Hpt->NumSets(); i++){
	   char mode='t'; // 'T' = gapped; 't' = ungapped scores.
	   Index1stBest[i] = 0;
#if 0	// DEBUG.
	 if(Hpt->NumSets() == 3){
	   FILE *tfp=open_file("junk",".cma","w");
	   PutCMSA(tfp,MainCMA); fclose(tfp);
	 }
#endif
	   WorstToBest[i]=SortByScoreCMSA(stderr,mode,Index1stBest[i],(double)MinSeed2CsqScore[i],MainCMA,i);
	   name=Hpt->ElmntSetName(i);
	   Int4 InSet=0;
	   if(Hpt->FixedElmntSet(i)){	// e.g., Random=20000.
		Index1stBest[i]=1; // don't sample sets with == definition...
	   	if(efp) fprintf(stderr,"%d.%s ***** WILL NOT BE SAMPLED *****\n", i,name);
	   } else if((InSet=CardSet(InitSet[i])) > 0){
#if 0
		Int4 tenthway = (Int4) ceil(0.90*InSet);
		assert(tenthway > 0);
	   	if(Index1stBest[i] > tenthway) Index1stBest[i] = tenthway; 
		// don't let it go beyond 90% of each input subset .
#endif
		// let empty sets exist?
	   } else if(Hpt->TypeOfSet(i) == '?'){ 	// a Misc set...
	   	// fprintf(stderr,"(\"%s?\"): Index1stBest[%d] = %d\n",name,i,Index1stBest[i]);
		Index1stBest[i]=1;	// doesn't matter as set is empty
	   } else {
		fprintf(stderr,"Too few sequences (%d) assigned to subgroup %d ('%s').\n",
			InSet,i,Hpt->ElmntSetName(i));
		// print_error("fatal error: Consider modifying the hyperpartition.");
	   }
	   if(efp) fprintf(stdout,"(%s): Index1stBest[%d] = %d\n",name,i,Index1stBest[i]);
	   if(Index1stBest[i] == 0) print_error("fatal error in mcs_typ::GetElmntSets()");
	   if(i > NumDisplayCMA) continue;	// temporary kluge; put non-display sets last.
	   BooLean found=FALSE;
	   for(Int4 jj=1; jj <= NumDisplayCMA; jj++){
		if(OldDisplayCMA[jj] && strcmp(name,NameCMSA(OldDisplayCMA[jj]))==0){
		  if(!found){
			DisplayCMA[i] = OldDisplayCMA[jj]; OldDisplayCMA[jj]=0; found=TRUE;
		  } else print_error("*.sma file contains identically-labeled alignments");
		}
	   } if(!found){
		 fprintf(stderr,"'%s': ",name);
		 // fprintf(stderr,"display set not found within *.sma input file!\n");
		 print_error("display set not found within *.sma input file!");
	   } else if(efp) fprintf(efp,"%d.%s == %d.%s\n",i,name,i,NameCMSA(DisplayCMA[i]));
	} 
	free(OldDisplayCMA);
// fprintf(stderr,"DEBUG OomcBPPS()\n");
	return Hpt->NumSets();
}

