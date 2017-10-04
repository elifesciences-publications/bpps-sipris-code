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

#include "lha_typ.h"
#include "tax_typ.h"

#define	USAGE_START2	"USAGE: bpps 3 <file_name> target_node [options]\n\
    Input files: \n\
       <file_name>.hpt: hyperpartion (tree) file indicating the relationships between nodes in the hierarchy\n\
       <file_name>.tpl: a set of template alignments showing how to align between levels in the hierarchy\n\
         There is one template for each node along the lineage except for the leaf node\n\
       <file_name>.cma: multiple input cma files, one for each node along the lineage\n\
         The first sequence in each cma file must occur in the corresponding template\n\
    Output files: \n\
       <file_name>_out.mma: lineage hierarchical alignment cma file.\n\
       <file_name>_out.sma: lineage hierarchical alignment seed alignment file.\n\
       <file_name>_out.hpt: lineage hierarchical alignment hyperpartition (same as input).\n\
\n\n"

Int4	PutRepSeqsCMSA(FILE *fp, cma_typ cma)
{
	Int4 Number,KeySeq=0,put_best=50;
	BooLean	DistinctPhyla=TRUE,*skip;
	double prob,best;
	FILE	*efp=0; // efp=stderr;
	dh_type dH=NULL;
	Int4    Score,I,J,best_J=0,N=NumSeqsCMSA(cma);
	if(nBlksCMSA(cma) != 1) print_error("this option requires only one blk");
	if(N < 2){ print_error("only one sequence in cma file"); }
	best=-9999999999.0;
	for(J=1; J <= N; J++){
              	prob = GetProbCMSA(1,J,cma);
		// prob = GetGappedProbCMSA(1,J,cma);
		if(prob > best) { best=prob; KeySeq=J; }
	}
	// if(KeySeq > N) print_error("-Best option input error (key seq > number seqs");
	if(KeySeq > N) fprintf(stderr,"-Best option (key seq > number seqs)\n");
	if(efp) fprintf(stderr,"KeySeq=%d\n",KeySeq);
	dH=dheap(N+2,4); NEW(skip,N+3,BooLean); 
	Int4 SelfScore=PseudoAlnScoreCMSA(KeySeq, KeySeq,cma);
	insrtHeap(KeySeq,-((keytyp)SelfScore + 0.01),dH); 
	if(KeySeq==J) fprintf(stderr,"KeySeq=%d; Self Score = %d\n",J,Score);
	for(J=1; J <= N; J++){
		skip[J]=TRUE;
		if(J==KeySeq) continue;
		e_type sE=TrueSeqCMSA(J,cma);
		if(DistinctPhyla){
		  if(KingdomSeq(sE) == 'U') continue;	// kingdom is unknown...
		  if(KingdomSeq(sE) == 'X') continue;	// kingdom is unknown...
		  if(KingdomSeq(sE) == 0) continue;	// kingdom is unknown...
		}
		Score=PseudoAlnScoreCMSA(KeySeq, J,cma);
		insrtHeap(J,-(keytyp)Score,dH); 
	}

	Int4	NumPhyla=0,NumHits;
	BooLean	*Printed=0;
	Int4    *phyla=GetPhylaSeqSet(efp, &NumPhyla, TrueDataCMSA(cma));
	NEW(Printed,NumPhyla+3,BooLean);
	Int4 *list; NEW(list, N+3,Int4); 

	for(NumHits=I=0; I < put_best && !emptyHeap(dH); ){
               	assert((J=delminHeap(dH)) != 0);
		if(!Printed[phyla[J]]) {
			skip[J]=FALSE; NumHits++;
			// fprintf(stderr,"Seq Num: %d\n",J);
			if(DistinctPhyla) Printed[phyla[J]]=TRUE;
			I++; list[I]=J;
		}
	} Nildheap(dH); 
	if(efp) fprintf(stderr,"I=%d; N=%d; KeySeq=%d; NumHits=%d.\n",I,N,KeySeq,NumHits);
	if(NumHits == 0){ J=KeySeq; skip[J]=FALSE; NumHits++; I++; list[I]=J; }
	for(J=1; J <= N; J++){ if(skip[J]){  I++; assert(I <= N); list[I]=J; } }
	assert(I == N);
	if(NumHits > 0) PutSelectOneCMSA(fp,skip,list,cma);
	free(list); free(skip); free(phyla); free(Printed);
	return NumHits;
}


int	run_hybrid(int argc,char *argv[],FILE *mfp,FILE *hfp,FILE *stfp,FILE *smfp)
/******************************************************************************
  Hybrid procedure: Extends out the alignment to include insert regions within a specified
  node by removing inserts and adding gaps to the background alignments.
  0. Takes Hieraln output as input (so eventually combine these into one procedure).
  1. Should make this so that starts with child of root and goes down to a leaf for
     a particular lineage; thus expand at each level.
  2. Modify Hpt to include lineage only nodes.
  3. Modify Hieraln to retain patterns (GISMO sampling over inserted columns only?).
  4. One lineage output for each leaf node? (Or limit to specified leaf nodes?)
 ******************************************************************************/
{ 
	Int4	arg,i,j,k,s,x,y,z,n,target;
	a_type	AB;
	UInt4   seed=7061950;
	cma_typ	cma=0;
	FILE	*fp,*efp=0;
	char	str[300],mode=' ';

	//================== 0. Get arguments. ===================
	Int4 time1=time(NULL); 
	if(argc < 3) print_error(USAGE_START2);
	target=atol(argv[2]);
	for(arg = 3; arg < argc; arg++){
	   // fprintf(stderr,"argv[%d] = '%s'\n",arg,argv[arg]);
	   if(argv[arg][0] != '-') print_error(USAGE_START2);
	   switch(argv[arg][1]) {
	     case 's': if(sscanf(argv[arg],"-seed=%d",&seed)!=1) print_error(USAGE_START2); 
		break;
             case 'x': mode = 'x'; break;
	     default: print_error(USAGE_START2);
	   }
	}
        sprintf(str,"%s_%d",argv[1],target);
	if(efp){ 
		fp = open_file(str,".cmd","w");
		for(i = 0; i < argc; i++) { fprintf(fp,"%s ",argv[i]); } 
	}
	if(seed == 7061950) {  // not provided by user
                seed=(UInt4) time(NULL)/2; // seed = (UInt4) time(NULL);
                if(efp) fprintf(fp,"-seed=%d",seed);
        } if(efp) { fprintf(fp,"\n"); fclose(fp); } sRandom(seed);

	//================== 1. Read input files. ===================
	AB = MkAlpha(AMINO_ACIDS,PROT_BLOSUM62);
	Int4 NumberTPLs=0,J,NumberCMAs;
        FILE *tfp=open_file(argv[1],".tpl","r");
        cma_typ *MULTI_TPL=MultiReadCMSA(tfp,&NumberTPLs,AB); fclose(tfp);
        fp=open_file(argv[1],".cma","r");
        cma_typ *MULTI_CMA=MultiReadCMSA(fp,&NumberCMAs,AB); fclose(fp);

	//================== 1. Find the target lineage hierarchy. ===================
	Int4	NumHpts=0;
	hpt_typ *TargetHpt=0,*hpt,*Hpt=MultiReadHpt(NumHpts,argv[1]);
	if(target < 2 || target >= Hpt->NumSets()){ 
	        fprintf(stderr,"Target (%d) is out of range (2..%d)\n",target,Hpt->NumSets());
		print_error(USAGE_START2);
	}
	char	*TargetName=Hpt->SetName(target);
	for(hpt=Hpt; hpt; hpt=hpt->RtnNext()){
		if(hpt->RtnFileName() == 0) print_error("hpt file input error");
		if(strcmp(hpt->RtnFileName(),TargetName) == 0){
		   TargetHpt=hpt; if(efp) hpt->Put(stderr); break;
		}
	} if(TargetHpt==0) print_error("target hierarchy not found!");
	if(NumberCMAs != (Hpt->NumSets() -1)){
	   fprintf(stderr,"NumberCMAs = %d + 1 != Hpt->NumSets() = %d\n",NumberCMAs,Hpt->NumSets());
	   print_error("NumberCMAs != (Hpt->NumSets() -1)");
	}

	//================== 1. Read Sequence Sets. ===================
        fp=open_file(argv[1],".sets","r");
	Int4	*NumSets=0;	NEW(NumSets, Hpt->NumSets()+3, Int4);
        set_typ	 **Sets=0,*sets=0;  NEWP(Sets,Hpt->NumSets()+3,set_typ);
	for(i=1; (sets=ReadSets(fp,n)) != 0; i++){	// skip root; start at 2.
	    Sets[i]=sets; NumSets[i]=n;
	    if(efp){
		fprintf(stderr,"%d.%s (%d): ",i,Hpt->SetName(i),SetN(Sets[i][1]),n);
	        for(j=1; j <= n; j++) fprintf(stderr," %d",CardSet(Sets[i][j]));
	        if(i==target) fprintf(stderr," <-- target");
	        fprintf(stderr,"\n");
	    }
	} fclose(fp);

	//================== 2. Check input files. ===================
	Int4	NumParent,*Parent=0,*S2I,*I2J,*J2S;
	if(!Hpt->IsTree(Parent)) print_error("FATAL: hpt file does not correspond to a tree");
	for(NumParent=0,i=1; i < Hpt->NumSets(); i++){ if(Hpt->IsParentNode(i)) NumParent++; }
	if(NumParent != NumberTPLs){
	   fprintf(stderr,"NumberParent=%d; NumberTPLs = %d\n",NumParent,NumberTPLs);
	   print_error("FATAL: tpl files != cma files - 1");
	}

	//======== 3. Find all Leaf nodes.
	set_typ Leaves=Hpt->MkLeafSet();

	//======== 4. Find the lineage for each node.
	Int4	child,**Lineage; NEWP(Lineage,Hpt->NumSets() +3, Int4);
	set_typ	*SubTree; NEW(SubTree,Hpt->NumSets() +3, set_typ);
	for(s=1; s < Hpt->NumSets(); s++){   
		NEW(Lineage[s],Hpt->NumSets() +3, Int4);
		SubTree[s]=Hpt->MkSubTreeSet(s);
		for(i=0,j=s; j != 0; j=Parent[j]){ i++; Lineage[s][i]=j; } Lineage[s][0]=i;
	}

	//======== 5. Get the target lineage node set.
	set_typ	LineageSet=Hpt->MkLineageSet(target);
	fprintf(stderr,"Target node (%d.%s) lineage set:",target,Hpt->ElmntSetName(target));
	PutSet(stderr,LineageSet);

	//======= 6. See which nodes correspond to the input alignments. =========
	//=======    and which template files correspond to each node. =========
	NEW(S2I,Hpt->NumSets() + 3, Int4); NEW(I2J,Hpt->NumSets() + 3, Int4);
	NEW(J2S,Hpt->NumSets() + 3, Int4);
S2I[1]=1; I2J[1]=1; J2S[1]=1;
	for(i=2; i <= NumberCMAs; i++){
	   char *name,*Name=NameCMSA(MULTI_CMA[i]);
	   for(s=2; s < Hpt->NumSets(); s++){	// Find mapping between cma and hpt.
		name=Hpt->ElmntSetName(s);
		if(strcmp(name,Name)==0){
		   S2I[s]=i;
		   for(z=s; !MemberSet(z,LineageSet); z=Parent[z]) ;	// Find a target lineage node.
		   assert(z > 0); Name=Hpt->ElmntSetName(z);
		   if(z == target && MemberSet(z,Leaves)) break;  
	   	   for(j=1; j <= NumberTPLs; j++){	// find the right template file.
	              name=NameCMSA(MULTI_TPL[j]);
		      if(strcmp(name,Name)==0){
			I2J[i]=j; if(J2S[j]) assert(J2S[j] == z); else J2S[j]=z; break; 
		      }
		   }
		   if(j > NumberTPLs){ print_error("FATAL: tpl inconsistent with hpt and cma"); } 
		   break;	// found corresponding template file at this point.
		}
	   }
	   if(s >= Hpt->NumSets()){
		fprintf(stderr,"%d. \"%s\": no match! (j=%d)\n",i,Name,I2J[i]);
		print_error("FATAL: cma file does not correspond to a hpt node");
	   } else if(I2J[i]){	// this will be zero for target node.
		if(efp) fprintf(stderr,"%d(%d) (\"%s\"): template #%d = %d.\"%s\".\n",
			i,s,NameCMSA(MULTI_CMA[i]),I2J[i],J2S[I2J[i]],NameCMSA(MULTI_TPL[I2J[i]]));
	   } else if(efp) fprintf(stderr,"%d (\"%s\"): target node (no template).\n",
			i,NameCMSA(MULTI_CMA[i]));
	}

	//======== 7. Find the node sets for each node in the target lineage.
	for(child=0,s=target; s != 0; child=s,s=Parent[s]){
		if(child){
		    for(i=target; i != 0; i=Parent[i]){
			IntersectNotSet(SubTree[s],SubTree[i]); // SubTree[s] is modified.
			if(i == child) break;
		    }
		}
	}

	//======== 7b. Print out for debugging.
     if(efp){
	fprintf(stderr,"\nNodes to be collapsed into the lineage node (1st on list):\n");
	for(s=target; s != 0; s=Parent[s]){ PutSet(stderr,SubTree[s]); }
	fprintf(stderr,"\nLineages for each (first) node:\n");
	for(s=1; s < Hpt->NumSets(); s++){   
		for(i=1; Lineage[s][i]; i++) fprintf(stderr,"%d ",Lineage[s][i]); 
		if(MemberSet(s,SubTree[target])) fprintf(stderr," <- %d",target);
		if(s==target) fprintf(stderr," (target)");
		fprintf(stderr," (depth = %d)\n",Lineage[s][0]);
	}
	fprintf(stderr,"\nKey nodes in target lineage:\n");
	for(z=target; z != 0; z=Parent[z]){
	  fprintf(stderr,"Node %d subTree: ",s); PutSet(stderr,SubTree[z]);
	}
     }
	//======= 8. Collapse nodes to match nearest node in the target lineage!
	//      Save these cma files for subsequence expansion into a hybrid alignment.
	cma_typ *MultiCMA; NEW(MultiCMA, NumberCMAs+3, cma_typ);
	MultiCMA[1]=CopyCMSA(MULTI_CMA[1]);	// assume root is first!
	for(z=target; z != 0; z=Parent[z]){
	  if(efp) fprintf(stderr,"************ node %d in target lineage: ************\n",z);
	  if(CardSet(SubTree[z]) == 1){	// this is the target set for a leaf node.
		assert(MemberSet(z,LineageSet));
		if(z==1) continue;
		assert(MultiCMA[z]==0);
		i=S2I[z]; MultiCMA[z]=CopyCMSA(MULTI_CMA[i]);
		continue;
	  } 
	  // for(x=Hpt->NumSets(); x > 0; x--)
	  for(x=2; x < Hpt->NumSets();  x++)
	  {
assert(MultiCMA);
	     if(MemberSet(x,SubTree[z])){	
		if(x==z){	// already at nearest lineage node...use it as is.
		  assert(MultiCMA[x]==0); 
		  if(MULTI_CMA[S2I[x]]) MultiCMA[x]=CopyCMSA(MULTI_CMA[S2I[x]]);
		  continue; 
		}
	        if(efp) fprintf(stderr,"\n=== node %d in node %d subtree: ===\n",x,z);
		s=x; i=S2I[x]; cma_typ	key_cma=MULTI_CMA[i];
		if(key_cma==0) continue;
		Int4 St;
		i=S2I[x]; j=I2J[S2I[s]]; 
		e_type tE=0,cE=TrueSeqCMSA(1,key_cma);
		cma_typ *out_cma,*in_cma; NEW(in_cma,Hpt->NumSets()+3,cma_typ);
		for(y=1; y < NumSeqsCMSA(MULTI_TPL[j]); y++){
		   	tE=TrueSeqCMSA(y+1,MULTI_TPL[j]);
			if(IsSubSeq(tE,cE,&St,FALSE)){ in_cma[y]=CopyCMSA(key_cma); break; }
		} assert(y < NumSeqsCMSA(MULTI_TPL[j]));
		out_cma=ConversionViaTemplateCMSA3(MULTI_TPL[j],in_cma); 
		assert(out_cma[y] != 0); 
		if(efp){
	          fprintf(stderr,"in: %d.%s; len=%d\n",s,NameCMSA(in_cma[y]),LengthCMSA(1,in_cma[y]));
	          fprintf(stderr,"out: %d.%s; len=%d (parent = %d)\n",s,
			NameCMSA(out_cma[y]),LengthCMSA(1,out_cma[y]),Parent[s]);
		}
		// WARNING: data set used for both in and out cma; RenameCMSA will free name!!!
		sprintf(str,"%s",NameCMSA(in_cma[y])); RenameCMSA(str,out_cma[y]);
		key_cma=out_cma[y]; NilCMSA(in_cma[y]); free(out_cma); free(in_cma);
		assert(MultiCMA[x]==0); MultiCMA[x]=key_cma;
		   // s=Parent[s];
		// while(!MemberSet(s,SubTree[target]));
		// PutCMSA(tfp,MultiCMA[x]);
	     } 
	  }
	}
#if 1	// Debug...
  if(mode != 'x'){
assert(MultiCMA);
        sprintf(str,"%s_%d",argv[1],target);
        tfp=open_file(str,"_shrink.cma","w"); 
	for(s=1; s < Hpt->NumSets();  s++) if(MultiCMA[s] != 0) PutCMSA(tfp,MultiCMA[s]);
	fclose(tfp);
  }
#endif

	//======= 9. Expand each SubTree set member's alignment to match the target node. 
	//   Nodes not in the target subtree only need to have deletions added.
	cma_typ *OutCMA; NEW(OutCMA, Hpt->NumSets()+3, cma_typ); 
	for(s=1; s < Hpt->NumSets(); s++){   
	   if(MemberSet(s,SubTree[target])){	// these are already of the right length.
	        if(efp){ fprintf(stderr,"s=%d element of ",s); PutSet(stderr,SubTree[target]); }
		OutCMA[s]=CopyCMSA(MultiCMA[s]); // PutCMSA(tfp,MultiCMA[s]); 
		SetLevelCMSA(1,OutCMA[s]);
	   } else for(z=target; z != 0; z=Parent[z]){
	     if(0 && z==target){	// these are already in final form...
		assert(!MemberSet(s,SubTree[z])); continue; // this should not be reached.
	     }
	     if(MemberSet(s,SubTree[z])){	
	        if(efp){
	           fprintf(stderr,"\n======== s=%d.%s element of (%d.%s) ",
				s,Hpt->ElmntSetName(s),z,Hpt->ElmntSetName(z));
	           PutSet(stderr,SubTree[z]);
		}
		if(s==z){
		  for(j=1; MULTI_TPL[j]; j++){ if(J2S[j]==z){ break; } }
		  if(efp) fprintf(stderr,"s=%d; j=%d\n",s,j);
		} else {
#if 1
		  i=S2I[s]; j=I2J[i];
		  if(efp) fprintf(stderr,"s=%d; i=%d; j=%d\n",s,i,j);
#else
		  char *name,*Name=NameCMSA(MultiCMA[z]);
		  for(j=1; MULTI_TPL[j]; j++){
			name=NameCMSA(MULTI_TPL[j]);
			if(strcmp(name,Name)==0); break;
		  } // i=S2I[z]; j=I2J[i]; if(s==1) j=1;  else if(i==1) j=1;
#endif
		}
	        if(efp) fprintf(stderr,"template = %d.%s.\n",j,NameCMSA(MULTI_TPL[j]));
		e_type csqE=TrueSeqCMSA(1,MultiCMA[target]);
		OutCMA[s]=GrowViaTemplateCMSA(MULTI_TPL[j],MultiCMA[s],target,csqE);
		SetLevelCMSA(2,OutCMA[s]);
		break;
	     }
	   }
	}

	//======== 5. Print out the omcBPPS input files for target node. ============
        sprintf(str,"%s_%d",argv[1],target);
        tfp=0; // tfp=open_file(str,".cma","w");
        FILE *sfp=0;
	if(smfp==0) sfp=open_file(str,".sma","w"); else sfp=smfp;
	set_typ	SetXX,*OutSets=0; NEW(OutSets,Hpt->NumSets() +3, set_typ);
	set_typ	*SetMerge; NEW(SetMerge,Hpt->NumSets() +3, set_typ);
	cma_typ *MergedCMA; NEW(MergedCMA,Hpt->NumSets() +3, cma_typ);
	Int4	i_mrg=0,xx,II=0,NumOutSets=0,size_set=NumSeqsCMSA(OutCMA[1]) + Hpt->NumberRandom() + 1;
	for(Int4 p=1; p < Hpt->NumSets(); p++){	// Find a target lineage node.
	   if(!MemberSet(p,LineageSet)) continue;
	   assert(OutCMA[p]); sets=Sets[p]; n=NumSets[p];
	   fprintf(stderr,"Number Seqs in OutCMA[p=%d]=%d; SetN(sets[1]) = %d; n=%d\n",
		p,NumSeqsCMSA(OutCMA[p]),SetN(sets[1]),n);
	   assert(NumSeqsCMSA(OutCMA[p]) <= SetN(sets[1]));
	   if(p==target){
		char *Name=AllocString(NameCMSA(OutCMA[p])); 
		for(j=p,i=1; i <= n; j++,i++){ 
		  ReNameCMSA(Hpt->SetName(j),OutCMA[p]);
		  if(CardSet(sets[i]) > 0){ 
			i_mrg++; SetMerge[i_mrg]=CopySet(sets[i]); MergedCMA[i_mrg]=OutCMA[p];
			if(tfp) PutInSetCMSA(tfp,sets[i],OutCMA[p]); 
		  }
	          SetXX=MakeSet(size_set); ClearSet(SetXX); 
		  for(xx=CardSet(sets[i]); xx > 0;  xx--){ II++; AddSet(II,SetXX); }
		  NumOutSets++; OutSets[NumOutSets]=SetXX; 
#if 0
		  if(j==p){ PutRepSeqsCMSA(sfp,OutCMA[p]); }
		  else {	// WARNING: what if set is empty!!??
	            e_type csqE=MkConsensusCMSA(sets[i],OutCMA[p]);
		    // PutSeq(stderr,csqE,AB);
	            cma_typ tmp_cma=RtnConSqAsCMSA(Hpt->SetName(j),0,csqE,AB); 
		    assert(LengthCMSA(1,OutCMA[p]) == LengthCMSA(1,tmp_cma));
	            PutCMSA(sfp,tmp_cma); TotalNilCMSA(tmp_cma); NilSeq(csqE);
		  }
#else
		  if(efp) fprintf(stderr,"  %d: %s (%s)\n",j,Hpt->SetName(j),NameCMSA(OutCMA[j]));
		  PutRepSeqsCMSA(sfp,OutCMA[j]);
#endif
		} ReNameCMSA(Name,OutCMA[p]); free(Name);
		  // PutCMSA(tfp,OutCMA[p]);	// then print the whole set.
	   } else {	// then print out everything except the subtree of c.
	        set_typ setP=CopySet(sets[1]); ClearSet(setP);
		for(j=p,i=1; i <= n; j++,i++){ 
		   if(MemberSet(j,SubTree[p])){	// for each node in parent subtree - target sublineage.
			UnionSet(setP,sets[i]);		// Add these to setP.
		   } 
		}
		if(CardSet(setP) > 0){
		  i_mrg++; SetMerge[i_mrg]=CopySet(setP); MergedCMA[i_mrg]=OutCMA[p];
		  if(tfp) PutInSetCMSA(tfp,setP,OutCMA[p]); 
		}
#if 1
	          SetXX=MakeSet(size_set); ClearSet(SetXX); 
		  for(xx=CardSet(setP); xx > 0;  xx--){ II++; AddSet(II,SetXX); }
		  NumOutSets++; OutSets[NumOutSets]=SetXX; 
#endif
		NilSet(setP);
		// WARNING: what if set is empty!!??
	        e_type csqE=MkConsensusCMSA(OutCMA[p]);
	        cma_typ tmp_cma=RtnConSqAsCMSA(Hpt->SetName(p),0,csqE,AB);
	        PutCMSA(sfp,tmp_cma); TotalNilCMSA(tmp_cma); NilSeq(csqE);
	   }
	} if(tfp) fclose(tfp); if(smfp==0) fclose(sfp);  else smfp=sfp;
#if 1
        if(mfp==0){ tfp=open_file(str,".mma","w"); } else tfp=mfp;
	PutMergedCMSA(tfp,(unsigned short)i_mrg,SetMerge,MergedCMA,0);
	for(i=1; i < Hpt->NumSets(); i++){
		if(SetMerge[i]) NilSet(SetMerge[i]);
	} free(SetMerge); free(MergedCMA);
	if(mfp==0) fclose(tfp); else mfp=tfp;
#else 
        if(mfp){ PutMergedCMSA(mfp,(unsigned short)i_mrg,SetMerge,MergedCMA,0); }
	else {
	   tfp=open_file(str,".mma","w"); 
	   PutMergedCMSA(tfp,(unsigned short)i_mrg,SetMerge,MergedCMA,0); fclose(tfp);
	}
	for(i=1; i < Hpt->NumSets(); i++){ if(SetMerge[i]) NilSet(SetMerge[i]); }
	free(SetMerge); free(MergedCMA);
#endif
	NumOutSets++; OutSets[NumOutSets]=MakeSet(size_set); ClearSet(OutSets[NumOutSets]);
	fprintf(stderr,"II = %d; NumSeqsCMSA=%d\n",II,NumSeqsCMSA(OutCMA[1]));
	TargetHpt->ChangeNumRandom(size_set-(II+1));
	for(xx=1; xx <= TargetHpt->NumberRandom(); xx++){ II++; AddSet(II,OutSets[NumOutSets]); }
	if(stfp==0) sfp=open_file(str,".sets","w"); else sfp=stfp;
	WriteSets(sfp,NumOutSets,OutSets); if(stfp == 0) fclose(sfp); else stfp=sfp;

        // tfp=open_file(str,".hpt","w"); TargetHpt->Put(tfp); fclose(tfp);
	if(hfp==0) tfp=open_file(str,".hpt","w"); else tfp=hfp;
	TargetHpt->Put(tfp,FALSE); if(hfp==0) fclose(tfp); else hfp=tfp;

#if 0
	//======== 5. Print out the hybrid alignments for each node. ============
        tfp=open_file(argv[1],".cma","w");
	for(i=2; i <= NumberCMAs; i++){		// i=2 (skip consensus sequence...
	     fprintf(stderr,"%d.%s + ",i-1,NameCMSA(MultiCMA[i-1])); 
	     fprintf(stderr,"%d.%s; tpl = ",i,NameCMSA(MultiCMA[i])); 
	     fprintf(stderr,"%d.%s\n",i,NameCMSA(MultiTPL[i])); 
	     if(i==target) PutHybridCMSA(tfp,tfp,MultiCMA[i],MultiCMA[1],MultiTPL[i],0,0); 
	     else PutHybridCMSA(tfp,0,MultiCMA[i],MultiCMA[1],MultiTPL[i],i-1,MultiCMA); 
	} fclose(tfp); tfp=open_file(argv[1],"_out.cma","r");

	Int4 NumberOut;
	cma_typ *MultiOut=MultiReadCMSA(tfp,&NumberOut,AB); fclose(tfp);
	sprintf(str,"_out.mma");
        fp=open_file(argv[1],str,"w"); PutMergedCMSA(fp,NumberOut,MultiOut); fclose(fp);
	for(j=1; j <= NumberOut; j++) TotalNilCMSA(MultiOut[j]); free(MultiOut);
#endif

#if 0	//=========== 6. Output seed sequence alignments. ================
	tfp=tmpfile();
	i=NumberCMAs; PutHybridCMSA(tfp,0,MultiTPL[i],MultiCMA[1],MultiTPL[i],0,0); 
	rewind(tfp); 
	cma_typ tmp_cma=ReadCMSA(tfp,AB); fclose(tfp);
        tfp=open_file(argv[1],"_out.sma","w");
	PutSinglesCMSA(tfp,tmp_cma); fclose(tfp); 
	TotalNilCMSA(tmp_cma);
#endif
	
	//=========== 7. Free up memory. ================
	free(Parent); free(S2I); free(I2J); free(J2S);

	for(i=1; Sets[i]; i++){ 
	   for(j=1; j <= NumSets[i]; j++) NilSet(Sets[i][j]);  free(Sets[i]);
	} free(Sets); free(NumSets);
	for(i=1; i <= Hpt->NumSets(); i++){
	   if(OutSets[i]) NilSet(OutSets[i]);
	   if(Lineage[i]) free(Lineage[i]); 
	   if(SubTree[i]) NilSet(SubTree[i]);
	} free(Lineage); free(SubTree); free(OutSets);
	for(i=1; i <= NumberCMAs; i++){
		if(LevelCMSA(OutCMA[i]) == 1) NilCMSA(OutCMA[i]); 
		else TotalNilCMSA(OutCMA[i]);
	} free(OutCMA); 
	for(i=1; i <= NumberCMAs; i++) NilCMSA(MultiCMA[i]); free(MultiCMA); 
	for(i=1; i <= NumberCMAs; i++) TotalNilCMSA(MULTI_CMA[i]); free(MULTI_CMA); 
	for(i=1; i <= NumberTPLs; i++) TotalNilCMSA(MULTI_TPL[i]); free(MULTI_TPL); 
	// for(i=1; i <= NumberSMAs; i++) TotalNilCMSA(MULTI_SMA[i]); free(MULTI_SMA);
	delete Hpt; NilSet(Leaves); NilAlpha(AB); NilSet(LineageSet);
	fprintf(stderr, "\ttime: %d seconds (%0.2f minutes)\n",
                        time(NULL)-time1,(float)(time(NULL)-time1)/60.0);
}


