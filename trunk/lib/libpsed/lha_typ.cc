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

e_type  MkCsqInSetCMSA(set_typ set,cma_typ cma)
{
}

set_typ	*MkSubSets(set_typ *set, Int4 NumSets, set_typ Merge)
/*****************************************************************
  1. Merge all of the sequences belonging to an internal node subtree 
     into a single set.
  2. Find the sequences within the merged set belonging to each node
     in the subtree.
  3. Use this in hieraln to keep track of subset within main alignments.

print out only the sequences 
  WriteSets(fp,Num,sets); ReadSets(fp,Num,sets);
 *****************************************************************/
{
	Int4	i,j,k,s,n,N,M;
	assert(NumSets >= 1);
	set_typ	OutSet=CopySet(set[1]); ClearSet(OutSet);
	set_typ	*RtnSet;

	//******* 1. Make sure input sets are disjoint. ********
	for(i=1; i < NumSets; i++){
	    for(j=i+1; j <= NumSets; j++) assert(CardInterSet(set[i],set[j])==0);
	}

	//******* 2. Create main merged set. ********
	for(i=0,s=1; s <= NumSets; s++){
	  if(MemberSet(s,Merge)){ i++; UnionSet(OutSet,set[s]); }
	}  NEW(RtnSet, i+3,set_typ);
	N=CardSet(OutSet); N++; // add another item to accomodate consensus seq. at start.

	//******* 3. Create merged subsets. ********
	for(M=0,s=1; s <= NumSets; s++){
	  if(MemberSet(s,Merge)){ M++; RtnSet[M]=MakeSet(N+3); ClearSet(RtnSet[M]); }
	}

	//******* 4. Add items to each subset based on main set. **********
	for(k=1,i=1; i < SetN(set[1]); i++){	// start with k=1 to accommodate consensus.
	   if(!MemberSet(i,OutSet)) continue;
	   for(j=0,s=1; s <= NumSets; s++){
		if(!MemberSet(s,Merge)) continue;
		j++;  assert(j <= M);
		if(MemberSet(i,set[s])){ k++; AddSet(k,RtnSet[j]); break; }
		// i,k not in other sets because they are disjoint!!
	   } assert(k <= N);
	} RtnSet[0]=OutSet; 
	return RtnSet;
}

cma_typ MkSubSetCMSA(set_typ *NodeSet, hpt_typ *Hpt, Int4 s, BooLean  AddCsq, cma_typ cma, set_typ **rtnSets)
/*****************************************************
  Keep track of sequences belonging to each Set.
  Merge in the background sets later...?
  Need lineage and 
 *****************************************************/
{
        Int4    i,j,N,sq,blk,*sites;
	//========  6a. Generate subgroup alignment with concensus. ========
        fprintf(stderr,"************** child = %d **************\n",s);
        set_typ subtree=Hpt->MkSubTreeSet(s); PutSet(stderr,subtree);
        fprintf(stderr,"Number of descendants = %d\n",Hpt->NumDescendants(s));
	// Note: skipping last, reject set...
	set_typ *RtnSet=MkSubSets(NodeSet,Hpt->NumSets()-1,subtree);
	set_typ	Set=RtnSet[0]; RtnSet[0]=0; *rtnSets=RtnSet; NilSet(subtree);
        a_type  AB=AlphabetCMSA(cma);
        e_type *EList; NEW(EList,CardSet(Set)+5,e_type);
        j=0;
        if(AddCsq){ j++; EList[j]=MKCsqSubCMSA(Set,cma); }
        for(sq=1; sq <= NumSeqsCMSA(cma); sq++){
                if(!MemberSet(sq,Set)) continue;
                j++; EList[j]=CopySeq(TrueSeqCMSA(sq,cma));
        }
        ss_type data=Array2SeqSet(EList,j,NameCMSA(cma),AB); // data absorbs EList!!!
        gss_typ *gss=gssCMSA(cma);
        cma_typ rcma=EmptyCMSA(nBlksCMSA(cma),LengthsCMSA(cma),data,gss->GapOpen(),
                gss->GapExtend(),PerNatsCMSA(cma),0,0);
        gsq_typ *gsq,*gsq2;
        j=0;
        if(AddCsq){
                gsq2 = new gsq_typ[1];
                Int4 *pos; NEW(pos,nBlksCMSA(rcma) + 5, Int4); pos[1]=1;
                for(i=1; i < nBlksCMSA(rcma); i++) pos[i+1]=pos[i] + LengthCMSA(i,rcma);
                gsq2->initialize(nBlksCMSA(rcma),LengthsCMSA(rcma),pos,TrueSeqCMSA(1,rcma));
                j++; ReplaceCMSA(j,gsq2,rcma); // replace sequence s in CMSA & fmodel.
                for(blk=1; blk <= nBlksCMSA(rcma); blk++) AddSiteCMSA(blk,j,pos[blk],rcma); free(pos);
        }
        for(sq=1; sq <= NumSeqsCMSA(cma); sq++){
                if(!MemberSet(sq,Set)) continue;
                gsq=gsqCMSA(sq,cma);
                gsq2 = new gsq_typ[1]; gsq2->initialize(*gsq);
                // gsq2= new gsq_typ(*gsq);     // core dumps due to array of [1] convention in gss_typ.
                // gsq->Put(stderr,AB); gsq2->Put(stderr,AB);
                j++; ReplaceCMSA(j,gsq2,rcma); // replace sequence s in CMSA & fmodel.
                sites=GetPosSitesCMSA(sq,cma);
                for(blk=1; blk <= nBlksCMSA(rcma); blk++) AddSiteCMSA(blk,j,sites[blk],rcma);
                free(sites);
        } NilSet(Set);
#if 0
	if(debug) PutHieralnDebug("b",0,"junk",Hpt->SetName(s),0,rcma);
// ExtendFakeToRealCMSA(cma);
        if(debug) PutHieralnDebug("a",0,"junk",Hpt->SetName(s),0,rcma);
#endif
	return rcma;
}


void	PutHieralnDebug(const char s[],Int4 id,char *prefix,char *suffix,e_type csqE,cma_typ cma)
// call as PutHieralnDebug(x,argv[1],Hpt->SetName(s),cma);
{
	char	str[500]; sprintf(str,"%s_%s_dbg%s%d",prefix,suffix,s,id); // RenameCMSA(str,cma);
	FILE *fp=open_file(str,".cma","w"); 
	if(csqE) PutWithCsqCMSA(fp,csqE,cma); else PutCMSA(fp,cma); fclose(fp); 
}

char *MkTemplateOperation(ssx_typ *SSX, ssx_typ *ssx, char *Oper[3])
{
    Int4	x,z;
    const Int4	p=1,s=2;
    char	*operation; NEW(operation, strlen(Oper[s]) * 2, char); 
    char	*ao,cs,cp,co; ao=operation; ao[0]='E'; ao[1]='M';
    char	State=Oper[p][1],state=Oper[s][1];
    Int4	xs,xp,xo;	// operation array indices.
    Int4	is=2,ip=2;	// cma file position indices.
    for(xs=xp=xo=2; (cs=Oper[s][xs]) != 'E'; xs++,xp++){
	cp=Oper[p][xp]; co=0;  // co=ao[z];
	switch (cs){
	  case 'm':	// match in Oper[s] relative to root.
	    if(cp == 'd'){		// m..d: Oper[p] deleted relative to root.
		co='m'; ip--;	// make this a match relative to template.
	    } else if(cp == 'I'){	// m..I: Oper[p] insert relative to root.
		co='d'; xs--; 		// Oper[s] is missing this insert...
	    } else { co='m'; }		// m..m: both a match, so retain.
	  break;
	  case 'd':	// deletion in Oper[s] relative to root.
		is--; // column not in cma alignment.
	     if(cp == 'm'){		// d..m: Oper[p] match relative to root.
		co='d'; 	// make deleted relative to template.
	     } else if(cp == 'I'){	// d..I: Oper[p] insertion relative to root.
		co='d'; xs--; 		// add an extra deletion 
	     } else { ip--; co=0;} 		// d..d: both deleted relative to root.
	     break;
	  case 'I':	// Insertion in Oper[s] relative to root.
	     if(cp == 'm'){		// I..m: Oper[s] has an insert to root.
		co='I'; xp--;	// also add an insert in s relative to template.
	     } else if(cp == 'd'){	// I..d: Oper[p] deletes next root column.
		co='I'; xp--; ip--;	// also an insert in s relatvie to template.
	     } else {			// I..I: both insertions - PROBLEMATIC!!
		//************** need some additional processing here. **
		if(state=='I' && State =='I') fprintf(stderr,"%s\n",ao);	// then 
           	assert(!(state == 'I' && State == 'I'));
           	assert(state != 'I' && State != 'I');
           	Int4	Xp,Xs,Ip,Is,LenP,LenS,LenAln,stXp,stIp,stXs,stIs,Diff;
           	Int4	StrtIs,StrtIp,StrtXs,StrtXp,endX,endI;
           	for(LenP=0,x=xp; Oper[p][x] == 'I'; x++){ LenP++; }
           	for(LenS=0,z=xs; Oper[s][z] == 'I'; z++){ LenS++; }
           	double	d,dd,DD=-999999999.9;
           	LenAln=MINIMUM(Int4,LenP,LenS); Diff = abs(LenP-LenS);
           
           	if(LenS < LenP){ 	// 
           	  StrtIs=is; StrtXs=xs;	// this never changes...
           	  for(x=0,stXp=xp,stIp=ip; x <= Diff; x++,stXp++,stIp++){
           	     Xp=stXp;Ip=stIp; Xs=StrtXs,Is=StrtIs; 
           	     for(dd=0.0,z=1; z <= LenAln; z++,Is++,Xs++,Xp++,Ip++){
           		d=SSX->ColScoreP2P(ssx,Ip,Is); dd+=d;
           		fprintf(stderr,"%d vs %d: %.3f nats\n",Ip,Is,d);
           	     } fprintf(stderr,"sum = %.3f nats\n",dd);
           	     if(dd > DD){ DD=dd; StrtIp = stIp; StrtXp=stXp; }
           	  } endX=StrtXp+LenAln-1; endI=StrtIp+LenAln-1;
           	  for(x=1; x <= LenP; x++){
           		if(xp >= StrtXp  && xp <= endX){	
           		   ao[xo]='m'; xo++; is++; ip++; xs++; xp++;  // set these to matches..
           		} else {	// these columns are missing from child alignment.
           		   ao[xo]='d'; xo++; is++; ip++; xs++; xp++;
           		}
           	  }
assert(Oper[s][xs] != 'E'); assert(Oper[p][xp] != 'E');
           	} else if(LenP < LenS){ 	// 
           	  StrtIp=ip; StrtXp=xp;	// this never changes...
           	  for(x=0,stXs=xs,stIs=is; x <= Diff; x++,stXs++,stIs++){
           	     Xs=stXs;Ip=stIs; Xp=StrtXp,Is=StrtIp; 
           	     for(dd=0.0,z=1; z <= LenAln; z++,Is++,Xs++,Xp++,Ip++){
           		d=SSX->ColScoreP2P(ssx,Ip,Is); dd+=d;
           		fprintf(stderr,"%d vs %d: %.3f nats\n",Ip,Is,d);
           	     } fprintf(stderr,"sum = %.3f nats\n",dd);
           	     if(dd > DD){ DD=dd; StrtIs = stIs; StrtXs=stXs; }
           	  } endX=StrtXs+LenAln-1; endI=StrtIs+LenAln-1;
           	  for(x=1; x <= LenP; x++){
           		if(xs >= StrtXs  && xs <= endX){	
           		   ao[xo]='m'; xo++; is++; ip++; xs++; xp++;  // set these to matches..
           		} else {	// these columns are missing from child alignment.
           		   ao[xo]='I'; xo++; is++; ip++; xs++; xp++;
           		}
           	  }
assert(Oper[s][xs] != 'E'); assert(Oper[p][xp] != 'E');
           	} else {		// else set both to matches for entire insert.
           	  for(x=1; x <= LenAln; x++){ ao[xo]='m'; xo++; is++; ip++; xs++; xp++; }
           	} co=0; xp--; xs--; ip--; is--;
		//***************************************************
	     }
	    break;
	  default: print_error("hieraln fatal error"); break;
	} if(co){ ao[xo]=co; xo++; } State=cp,state=cs; is++; ip++;
   } ao[xo]='E'; xo++; ao[xo]=0;
   return operation;
}

cma_typ	QuickOptimizeCMA(cma_typ cma,Int4 s, char **Oper)
//========  6c. Lengthen the subalignment for this subtree. ========
//========== The set includes all sequences in subtree. =============
{
	const Int4 NumIter=1;	// one iteration appears to be enough...
	cma_typ	rcma=0,xcma=cma;
	gmb_typ	*xgmb=0;
	// Int4    i,x,aa_per_io=500,aa_per_do=500,exp_ie=1,exp_de=1;
	Int4    i,x,aa_per_io=20,aa_per_do=150,exp_ie=1,exp_de=1;
	char	*ao,dms_mode='F';
	double	PriorWt=1.0;
	double	temp=0.0; 	// be conservative so don't ruin alignment...
	FILE	*efp=0; // efp=stderr;

	for(i=1; i <= NumIter; i++){
		xgmb = new gmb_typ(aa_per_io,aa_per_do,exp_ie,exp_de,xcma,dms_mode); 
		xgmb->SetPriorWt(PriorWt);
// rcma=xgmb->Sample(0,'S',0,temp,temp,30,stderr); 
xgmb->SetMaxIter(1);
xgmb->Sample(0,'S',0,temp,30,30,stderr);  
		rcma=xgmb->AddColumns(0.0,FALSE,0,0.75,FALSE); 
// return rcma;
		Oper[s]=xgmb->RtnAddOp(); // sets xgmb->AddOp=0;
		if(efp) fprintf(stderr,"Oper[%d]=%s\n",s,Oper[s]);
		for(ao=Oper[s],x=strlen(ao); x > 0; x--){
			if(ao[x] == 'd') ao[x]='I'; assert(ao[x] != 'D'); 
		} if(efp) fprintf(stderr,"Oper[%d]=%s\n",s,Oper[s]);
		delete xgmb;
		if(rcma == 0){ rcma=CopyCMSA(xcma); } // uses same gss & data...

		xgmb = new gmb_typ(aa_per_io,aa_per_do,exp_ie,exp_de,rcma,dms_mode); 
		xgmb->SetPriorWt(PriorWt);
		// if > 85% deletions then remove the column...
		xcma=xgmb->RmColumns(0.0,FALSE,0.85); 
		ao=xgmb->RtnAddOp(); // sets xgmb->AddOp=0;
		if(ao != NULL){		// then remove columns 
		   // compare Oper[s] with this...
		   if(efp) fprintf(stderr,"ao=%s\n",ao);
		   assert(strlen(ao) == strlen(Oper[s]));
		   for(x=1; ao[x] != 'E'; x++){
			if(ao[x]=='d') Oper[s][x]='d'; 
		   } if(efp) fprintf(stderr,"Oper[s]=%s\n",Oper[s]); free(ao);
		} delete xgmb;
		if(xcma != 0){ if(xcma != rcma) NilCMSA(rcma); rcma=xcma; xcma=0; }
// ExtendFakeToRealCMSA(rcma);
		xgmb = new gmb_typ(aa_per_io,aa_per_do,exp_ie,exp_de,rcma,dms_mode); 
		xgmb->SetPriorWt(PriorWt);
#if 1	// don't do sticky sampling...tends to mess up global alignment.
	{
		rcma=xgmb->Sample(0,'S',0,temp,30,30,stderr);  
		// rcma=xgmb->Sample(0,'S',0,temp); 
		// rcma=xgmb->Sample(0,'S',0,temp,-30,30,stderr); 
	}
#else
	{
		// char    Strategy[20]=" IVXFD    ";  // 'D' ruins the alignment?
		char    Strategy[20]=" IDV    ";  // does 'D' ruin the alignment?
        	double  sfrq[20]={ 0.0,0.25,0.25,0.25,0.25,0.25,0.0,0.0,0.0,0.0,0.0,0,0,0.0};
        	Int4    MinSticky=5;
        	xgmb->SetMaxIter(1);
        	// rcma=xgmb->SampleStickyTogether(stderr,MinSticky,Strategy,sfrq,temp,'R');
        	rcma=xgmb->SampleStickyTogether(stderr,MinSticky,Strategy,sfrq,temp,'O');
	} 
#endif
		delete xgmb; 
		if(xcma && xcma != cma && xcma != rcma) NilCMSA(xcma); xcma=0; 
// ExtendFakeToRealCMSA(rcma);
		if(i < NumIter){ xcma=rcma; rcma=0; }
	}
	return rcma;
}

cma_typ	*MkTemplateCMAs(hpt_typ *Hpt,cma_typ *NodeCMA,set_typ *ChildSet,e_type *CsqE,
				char **Oper, a_type AB)
//========  7. Make the template alignments using GAMBIT. ========
{
	Int4	p,i,j,s,Start,*Parent=0;
	// Int4    aa_per_io=500,aa_per_do=500,exp_ie=1,exp_de=1;
	// char    dms_mode='S';
	Int4    aa_per_io=20,aa_per_do=150,exp_ie=1,exp_de=1;
	char	dms_mode='F';
	double	PriorWt=1.0;
	BooLean	debug=FALSE;
	gmb_typ	*gmb,*xgmb;
	cma_typ	xcma,tcma,*TmplCMA; NEW(TmplCMA,Hpt->NumSets()+5, cma_typ);
	Hpt->IsTree(Parent);
	char	*operation;
	for(p=1; p < Hpt->NumSets(); p++){	// skip root and reject sets!!!
	   set_typ childset=ChildSet[p]; 	// set of child nodes of p (parent).
	   if(CardSet(childset) == 0) continue;	// this probably should not happen.
	   if(childset==0) continue;  // this is a leaf node.

	   //======= 7a. Create a sequence data set for the template alignment. ========
	   e_type *EList; NEW(EList,CardSet(childset) + 5,e_type); xcma=NodeCMA[p];
	   assert(CsqE[p] != 0); EList[1]=CopySeq(CsqE[p]); EList[2]=CopySeq(CsqE[p]);
	   for(j=2,s=2; s < Hpt->NumSets(); s++){	// skip root and reject sets!!!
		if(!MemberSet(s,childset)) continue;
	   	assert(CsqE[s] != 0);
	       	j++; EList[j]=CopySeq(CsqE[s]);
	   } ss_type data=Array2SeqSet(EList,j,NameCMSA(xcma),AB); // data absorbs EList!!!

	   //======= 7b. Create the template multiple alignment for parent (tcma). ========
	   gss_typ *gss=gssCMSA(xcma);
	   tcma=EmptyCMSA(nBlksCMSA(xcma),LengthsCMSA(xcma),data,gss->GapOpen(),
                				gss->GapExtend(),PerNatsCMSA(xcma),0,0);
	   AddSiteCMSA(1,1,1,tcma);	// add the root consensus seq to the alignment
	   AddSiteCMSA(1,2,1,tcma);	// add a 2nd root consensus seq to the alignment
	   RenameCMSA(Hpt->SetName(p),tcma);
	   xgmb = new gmb_typ(aa_per_io,aa_per_do,exp_ie,exp_de,xcma,dms_mode);
ssx_typ	*SSX=xgmb->RtnSSX();
	   fprintf(stderr,"######### Making parent node %d (\"%s\") template #########\n",p,Hpt->SetName(p));
	   for(j=2,s=2; s < Hpt->NumSets(); s++){	// skip root and reject sets!!!
		if(!MemberSet(s,childset)) continue;	  // if not in child set, then skip.
		j++; 
		if(Oper[p]) fprintf(stderr,"Oper[p=%d]=%s\n",p,Oper[p]);
		if(Oper[s]) fprintf(stderr,"Oper[s=%d]=%s\n",s,Oper[s]);
#if 1
		operation=xgmb->AlignSeq(CsqE[s],Start); // Sequence-to-profile alignment.
	        // operation=xgmb->AlignP2P_Operation(NodeCMA[s],Start); // Profile-to-profile align...
#else
		Start=1;
	        fprintf(stderr,"*********** node %d (\"%s\") versus parent %d (\"%s\") **********\n",
				s,Hpt->SetName(s),Parent[s],Hpt->SetName(Parent[s]));
		if(p==1) operation=AllocString(Oper[s]);	// then use as is.
		else {
		  gmb_typ *sgmb=new gmb_typ(aa_per_io,aa_per_do,exp_ie,exp_de,NodeCMA[s],dms_mode);
		  ssx_typ *ssx=sgmb->RtnSSX();
		  char *oper[3]; oper[1]=Oper[p]; oper[2]=Oper[s]; 
		  operation=MkTemplateOperation(SSX,ssx,oper);
		  delete sgmb;
		} 
#endif
		fprintf(stderr,"operation=%s\n",operation);
#if 1	// turn this into another routine...
		Int4	RmLeft=0,RmRight=0,limit=1;
		char	*Op=FixTmplOperation(operation,CsqE[s],Start,RmLeft,RmRight);
		// fprintf(stderr,"Operation=%s\n",Op);
		cma_typ ncma,mcma=0;	// node cma, modified cma...
#if 1		// assert only for now.
		assert(RmLeft == 0 && RmRight == 0);
#else		// this is tricky; wait to implement this...
		if(RmLeft > 0 || RmRight> 0){
			// trim node alignment.
			ncma=NodeCMA[s]; assert(nBlksCMSA(ncma) == 1);
			mcma=TrimBlkCMSA(ncma,1,RmLeft,RmRight,limit); 
			if(mcma==0) print_error("TrimBlkCMSA(): too many columns removed from alignment");
			NilCMSA(ncma); NodeCMA[s]=mcma;

			// trim consensus seq.
			e_type  tmpE=MkSubSeq(1+RmLeft,LenSeq(CsqE[s])-RmRight, CsqE[s]);
			NilSeq(CsqE[s]); CsqE[s]=tmpE;

			// trim consensus alignment.
			mcma=TrimBlkCMSA(tcma,1,RmLeft,RmRight,limit); 
			if(mcma==0) print_error("TrimBlkCMSA(): too many columns removed from alignment");
			NilCMSA(tcma); tcma=mcma;
		}
#endif
#endif
	        if(debug){ xgmb->AlignSeq(stderr,'S',CsqE[s],0.0); PutSeq(stderr,CsqE[s],AB); 
			   xgmb->AlignSeq(stderr,'S',TrueSeqCMSA(j,tcma),0.0); }
	        fprintf(stderr,"*********** node %d (\"%s\") versus parent %d (\"%s\") **********\n",
				s,Hpt->SetName(s),Parent[s],Hpt->SetName(Parent[s]));
		Int4    sites[3]; sites[1]=Start; sites[2]=0;  // don't need to set this, not used for gsq.
	        gsq_typ *gsq=new gsq_typ[1];
	        // gsq->initialize(operation,TrueSeqCMSA(j,tcma),sites); // sets sites[1] to position.
		gsq->initialize(Op,TrueSeqCMSA(j,tcma),sites); // sets sites[1] to position.
		ReplaceCMSA(j,gsq,tcma); AddSiteCMSA(1,j,sites[1],tcma); free(operation); free(Op);
	   } RenameCMSA(Hpt->SetName(p),tcma); TmplCMA[p]=tcma; 
	   // xgmb->AlignP2P(stderr,tcma); // core dump: assert(fakeE != NULL) fails in gsq_typ::AddTransitions()
	   delete xgmb;
	   if(debug) PutHieralnDebug("",2,"junk",Hpt->SetName(p),CsqE[p],NodeCMA[p]);
	} free(Parent);
	return TmplCMA;
}


cma_typ	PurgeSampleReAlign(char *name, cma_typ TrueMainCMA)
// call as rcma=PurgeSampleReAlign(argv[1], TrueMainCMA)
// replace this code; can skip this perhaps...?
{
	cma_typ	old_cma,new_cma=0,xcma,tcma,rcma,cma=0,CMA,*NodeCMA,*TmplCMA; 
	// Int4	aa_per_io=500,aa_per_do=500,exp_ie=1,exp_de=1;
	// double	PriorWt=2;
	// char	dms_mode='S',str[200];
	Int4    i,x,aa_per_io=20,aa_per_do=150,exp_ie=1,exp_de=1;
	char	dms_mode='F',str[200];
	double	PriorWt=1.0;
	a_type	AB=AlphabetCMSA(TrueMainCMA);
	gmb_typ	*gmb=0,*xgmb=0;
    
	Int4	card,similarity=0,maxiter=30;
	set_typ	InSet=MakeSet(NumSeqsCMSA(TrueMainCMA)+4); FillSet(InSet);
	Int4	min_card=(Int4) ceil(0.25*(double)NumSeqsCMSA(TrueMainCMA));
	min_card = MINIMUM(Int4,min_card,100); min_card = MAXIMUM(Int4,min_card,25);
        //=================== 3a. Purge main cma file. ====================
	Int4	percent_ident=50;  
	set_typ PurgeSet=0;
	do {
	   if(PurgeSet) NilSet(PurgeSet);
	   FILE *efp=0; // efp=stderr;
	   PurgeSet=RtnFastRepSetCMSA(efp, percent_ident,InSet,TrueMainCMA);
	   card=CardSet(PurgeSet);
	   fprintf(stderr," purge%d -> %d seqs\n",percent_ident,card);
	   percent_ident++;
	} while(card < min_card && percent_ident <= 80);
	FILE *fp=tmpfile(); PutInSetCMSA(fp,PurgeSet, TrueMainCMA); rewind(fp); 
	old_cma=ReadCMSA(fp,AB); fclose(fp);
	ExtendFakeToRealCMSA(old_cma);
	NilSet(PurgeSet);

	//============== 3b. Run gismo++ on purged set to realign. ===============
	sprintf(str,"%s_gmb",name); RenameCMSA(str,old_cma);
        Int4 MinSticky=2;
        double MaxFrctn=0.49,BildCut=0.0,prior_wt=0.0;
        gmb=new gmb_typ(aa_per_io,aa_per_do,exp_ie,exp_de,old_cma);
        if(prior_wt > 0) gmb->SetPriorWt(prior_wt); gmb->SetMaxIter(3); 
	new_cma=gmb->Sample(str,'S',similarity);
        do {
              rcma=gmb->AddColumns(BildCut,FALSE);
              if(rcma){
		delete gmb; if(new_cma) NilCMSA(new_cma); new_cma=rcma; 
		gmb=new gmb_typ(aa_per_io,aa_per_do,exp_ie,exp_de,new_cma);
		if(prior_wt > 0) gmb->SetPriorWt(prior_wt);
		gmb->SetMaxIter(3); rcma=gmb->Sample(str,'S',similarity);
		break;
              }
        } while(rcma); delete gmb;
        gmb=new gmb_typ(aa_per_io,aa_per_do,exp_ie,exp_de,new_cma);
        do {
              // rcma=gmb->SampleStickyTogether(stderr,MinSticky,MaxFrctn,str);
              rcma=gmb->SampleStickyTogether(stderr,MinSticky,MaxFrctn,str,0.0,'R');
              gmb->SetMaxIter(9); rcma=gmb->Sample(str,'S',similarity);
              break;
        } while(rcma); delete gmb;
	if(new_cma == 0) print_error("gismo++ failed to find conserved regions");
	else {
	   if(new_cma != old_cma) NilCMSA(old_cma); ExtendFakeToRealCMSA(new_cma); 
	} old_cma=0;
	NilSet(InSet);

	//============== 3d. Realign the entire main set of sequences. ===============
	gmb = new gmb_typ(aa_per_io,aa_per_do,exp_ie,exp_de,new_cma,dms_mode);
	rcma=gmb->CreateAlign(TrueDataCMSA(TrueMainCMA)); delete gmb; 
	TotalNilCMSA(new_cma); new_cma=0;	// new_cma data no longer needed.

	//============== 3e. Output new alignment. ===============
	// fp = open_file("junkOutX",".cma","w"); PutCMSA(fp,rcma); fclose(fp); fp=0;
	fprintf(stderr,"============ Creating full alignment. ===========\n");
	double temp=300;
ExtendFakeToRealCMSA(rcma);
	assert(NumSeqsCMSA(TrueMainCMA) == NumSeqsCMSA(rcma));
	return rcma;
}

//***************** hybrid routines ***************8

void	PutSinglesCMSA(FILE *fp,cma_typ cma)
{
	Int4 n,sq,N=NumSeqsCMSA(cma),MaxLen=0;
	e_type sE=0;
	BooLean *skip; NEW(skip,N+3,BooLean);
	for(sq=1; sq <= N; sq++){
	   sE=TrueSeqCMSA(sq,cma); MaxLen=MAXIMUM(Int4,LenSeq(sE),MaxLen); skip[sq]=TRUE;
	}
	for(sq=1; sq <= N; sq++){
	   	sE=TrueSeqCMSA(sq,cma);
		char str[202]; StrSeqID(str, 200, sE);
		RenameCMSA(str,cma); skip[sq]=FALSE;
		PutSelectCMSA(fp,skip,cma); skip[sq]=TRUE;
	} free(skip);
}

cma_typ	GrowViaTemplateCMSA(cma_typ Template, cma_typ SubGrp, Int4 target, e_type TargetCsq)
// expand the SuperGrp alignment to match the subgroup alignment by adding deletions in omitted regions.
#if 0	// Comments:
  a. Loosen the template requirement to allow extensions on either end of the template alignment.

  A. Find correspondence between the two IN_CMAs.
	0. Pass in two IN_CMA files: one the SuperGrp alignment the other the subgroup.
	1. Find all sequences in SuperGrp aln corresponding to a subgroup aln.
	2. 
  B. Order the SuperGrp sequences so that those at the top match the subgroup aln.
	1. 
  C. Create a fake expanded SuperGrp alignment for sampling.
	1. Follow first and second sequences in template alignment.
	2. MAke sure that first sequence has no indels at all.
	3. Put a deletion everywhere in InCMA where the 2nd Template has an insert.
	4. 
  D. Output a hybrid alignment with real subgroup sequences at top  and fake SuperGrp seqs at bottom.
	1. Order 
	2. Include both the hybrid aln and the original SuperGrp alignment in output (two CMAs).
#endif
{
	FILE	*fp=tmpfile(); // fp=stdout;
	Int4	i,j,o,s,r,L,tru,fk,x,Start,sq,Sq=0;
	e_type	tE,sE,fE;
	a_type	AB=AlphabetCMSA(Template);	
	FILE	*efp=0; // efp=stderr;

	//============ A. Confirm that the input is compatible. ===============
	// 1. Template seq #1 same length as Template cma file.
#if 1	// not necessary...
	if(LengthCMSA(1,Template) != LengthCMSA(1,SubGrp)){
		fprintf(stderr,"Len SubGrp = %d; Len Template = %d\n",
			LengthCMSA(1,Template),LengthCMSA(1,SubGrp));
		print_error("FATAL: SubGroup and Template alignments are incompatible");
	} 
#endif
	tE=TrueSeqCMSA(1,Template); fE=TargetCsq;
	// PutSeq(stderr,fE,AB);
	for(j=1,i=2; i <= NumSeqsCMSA(Template); j++,i++){	// template has extra first seq.
		sE=TrueSeqCMSA(i,Template); 
		// if(IsSubSeq(sE,fE,&Start,FALSE)){ Sq=i; assert(target == j); break; } else sE=0;
		if(IsSubSeq(sE,fE,&Start,FALSE)){ Sq=i; break; } else sE=0;
	} if(sE == 0) print_error("FATAL: Template file does not contain a match to subgroup file");
	// if(sE) AlnSeqSW(stderr,11, 1, sE, fE, AB);
	if(efp) fprintf(stderr,"Len SubGrp = %d; lenSeq(sE)=%d\n",LengthCMSA(1,SubGrp),LenSeq(sE));

	//============ B. Find subseq in Template Seq#2 with subgroup seq #1. ============
	char	rtn=IsSubSeq(sE,fE,&Start,FALSE);
	// rtn==1 --> E1 subseq of E2; ==2 --> E2 < E1; ==3 -->  E1 == E2; ==0 --> no subseq.
	assert(Start >= 0); assert(rtn==1 || rtn == 3);
	unsigned short RmLeft=Start,RmRight = LenSeq(fE)-LenSeq(sE)-RmLeft,LenTplSq=LenSeq(sE);
	if(efp) fprintf(stderr,"RmLeft=%d; RmRight=%d;Start=%d\n",RmLeft,RmRight,Start);

	//============ C. Find the template sequence for SubGrp... ===============
	char name[200],*Name;
	StrSeqID(name, 150, TargetCsq); Name=AllocString(name);
	tE=TrueSeqCMSA(Sq,Template);
	StrSeqID(name, 150, tE);
	if(strcmp(name,Name) != 0) print_error("FATAL: subgroup not in template file!");
	LenTplSq=LenSeq(TargetCsq); free(Name);

  	//============ D. Create Template aln with deletions for inserts in subgrp aln. ==============
	gss_typ *gss=gssCMSA(Template);
	char c,*operation=gss->Operation(Sq);	// This is target template alignment against ancestral node.
	if(efp) fprintf(stderr,"operation=%s\n",operation);
	Int4 ptr,Len=LengthCMSA(1,SubGrp), len=LengthCMSA(1,Template);
	if(efp) fprintf(stderr,"len subgrp = %d; len tmplt =%d\n",Len,len);
	Int4 NetLen=LenTplSq+RmLeft+RmRight;
	char	*buffer;  NEW(buffer,LenSeq(fE) + NetLen +9, char);
	if(efp) fprintf(stderr,"NetLen = %d = TplSeq (%d) + RmLeft (%d) + RmRight (%d)\n",
			NetLen,LenTplSq,RmLeft,RmRight);
	fprintf(fp,"[0_(1)=%s(%d){go=10000,gx=2000,pn=1000.0,lf=500,rf=500}:\n(%d)",
		NameCMSA(SubGrp),NumSeqsCMSA(SubGrp),NetLen);
	for(i=1; i<= NetLen;i++) fprintf(fp,"*"); fprintf(fp,"\n\n");
	Int4 start=1,end=RmRight;
	if(RmLeft > 0){ start = 2; }
	if(RmRight > 0){ end--; }
	for(sq=1; sq <= NumSeqsCMSA(SubGrp); sq++){
	   ptr=0;tru=0;fk=0;
	   if(RmLeft > 0){ tru++; fk++; buffer[ptr]='A'; ptr++; } // work around no '-' at N-terminus issue.
	   for(i=start; i <= RmLeft; i++){ fk++; buffer[ptr]='-'; ptr++; }
	   for(o=1,s=0; (c=operation[o]) != 'E'; o++){
		switch (c){
		  case 'M': case 'm': 	// target alignment matches root template.
		   {
                       s++; r=ResidueCMSA(1,sq,s,SubGrp);
                       if(r != UndefAlpha(AB)){
				fk++; tru++; buffer[ptr]=AlphaChar(r,AB); 
		       } else { 
#if 1
			  if(IsDeletedCMSA(1,sq,s,SubGrp)){ fk++; buffer[ptr]='-'; }
			  else { fk++; tru++; buffer[ptr]='X'; }
#else
			  fk++; buffer[ptr]='-'; 
#endif
		       } ptr++;
		   } break;
		  case 'I': {		// target alignment has an insert compared to root template.
			fk++; buffer[ptr]='-'; ptr++; 	// add deletions to subcma...
		  } break; 
		  case 'd': case 'D': {	// target alignment is deleted at this position.
		   	s++; 		// so delete these as well by skipping this position.
		  } break;	
		  default: break;
                } 
	   } 
	   for(i=1; i <= end; i++){ fk++; buffer[ptr]='-'; ptr++; }
	   if(RmRight > 0){ tru++; fk++; buffer[ptr]='A'; ptr++; } // work around no '-' at C-terminus. }
 	   buffer[ptr]=0;
	   fprintf(fp,"$%d=%d(%d):\n>",sq,tru,fk); PutSeqInfo(fp,TrueSeqCMSA(sq,SubGrp)); 
           fprintf(fp,"{()%s()}*\n\n",buffer);
	} fprintf(fp,"_0].\n"); 
	free(operation); free(buffer); rewind(fp);
	cma_typ rcma=ReadCMSA(fp,AB); fclose(fp);
	return rcma;
}



