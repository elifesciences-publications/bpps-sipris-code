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

double	omc_typ::ComputeCrossScore(Int4 key,Int4 site, set_typ FGSet,set_typ BGSet, cmh_typ *cmh)
// Compute FG and BG sets corresponding to cross-conserved sibling nodes...
// GetOptPttrnLPR() computes LPR using weighted counts as passed into lpr_typ!!!
{
	FILE *ofp=stderr;
	double	d=0;
	if(cmh->InHeap(FGSet,BGSet) == 0){	// if this configuration is not yet in the heap.
		// 1. Set up data structures for the analysis.
		Int4    i=site,j,y,r,x,NumCols=mcs->RtnDefaultMaxCol();
	        char TypNode='R';
		Int4	*Parent; Hpt->IsTree(Parent);
		ClearSet(TmpFG); ClearSet(TmpBG);
		assert(Set==0);
		NEW(Set,Hpt->NumSets() +3, set_typ);

		// 2. Create the sequence sets for the relevant nodes...
		set_typ *set=mcs->CopyOfSeqSets();
		e_type tmpE=GetSeqAsCsqCMSA(set[key],TrueMainCMA);
                Set[1]=CopySet(set[1]);
		for(y=2; y < Hpt->NumSets(); y++){     // find sibling nodes...in tree format.
		   if(MemberSet(y,FGSet)){ 
                	Set[y]=MakeSet(mcs->GetSetSize());  mcs->CopySubTreeSeqSet(y,Set[y]);
		   } else if(MemberSet(y,BGSet)){
                	Set[y]=MakeSet(mcs->GetSetSize());  mcs->CopySubTreeSeqSet(y,Set[y]);
		   } 
                }
		// 3. Create new FG and BG sequence sets.
	        for(y=2; y < Hpt->NumSets(); y++){
		   if(Parent[y] != Parent[key]) continue;
		   if(MemberSet(y,FGSet)){ UnionSet(TmpFG,Set[y]); }
		   else if(MemberSet(y,BGSet)){ UnionSet(TmpBG,Set[y]); }
		}
		// 4. Find the optimal pattern for the new FG and BG seqs.
		if(ofp) fprintf(ofp,"%d.     FG=%d; BG=%d: --- ",i,CardSet(TmpFG),CardSet(TmpBG));
		sst_typ *zsst=this->GetOptPttrnLPR(0,TmpFG,TmpBG,d,NumCols,TypNode,tmpE);
		if(ofp) fprintf(ofp,"LPR=%.3f ---\n", d);

	        // 5. If found something new and significant, then insert into cmh heap.
		x=cmh->Insert(FGSet,BGSet,d,zsst); // free(zsst);
		if(x > 0 && ofp) fprintf(ofp,"%d. insert = %d; key = %.2f.\n",i,x,d);
		// 6. Free up memory.
		if(x==0) free(zsst);
		for(j=1; j <= Hpt->NumSets(); j++){
			if(Set[j]) NilSet(Set[j]); 
			if(set[j]) NilSet(set[j]); 
		} free(set); free(Set); Set=0; free(Parent); NilSeq(tmpE);
	} return d;
}

#define DUMMY_MTRX "-4  -4  -4 \
                    -4   5  -4 \
                    -4  -4   5 "

double  calc_binary_bild(UInt4 *kount, double wt_factor, unsigned char key_res, a_type AB)
/*	Call this program once for each column, with a vector of independent	*/
/*	counts for the 20 amino acids, in this order:   ARNDCQEGHILKMFPSTWYV	*/
/*      as well as a vector of background frequencies in the same order.  The	*/
/*	program returns a floating point column score, with the unit of nats.	*/
{
	unsigned char r;
	double	Count;			/*  Aggregate number of observations	*/
	double	info;			/*  BILD score in nats			*/
	double	logprob;		/*  Log-likelihood of data		*/
	double	count[3]; count[1]=count[2]=0.0;

	Count=count[1]=(double)kount[key_res]/wt_factor;
	for(r=0; r <= nAlpha(AB); r++){ 
	   if(r != key_res){ count[2] += (double) kount[r]/wt_factor; }
	} // if(Count<=1.0) return(0.0);
	Count+=count[2];
#if 1
        logprob = lgamma(1 + count[1]) - lgamma(1 + count[2]);
#elif 0
        logprob = log(1 + count[1]) - log(1 + count[2]);
#else
	logprob = -lgamma(1+Count)
        logprob += lgamma(1 + count[1]) + lgamma(1 + count[2]);
        logprob -= count[1] * log(0.05) + count[2]*log(0.95);
#endif
	return(logprob);
}

double	omc_typ::CrossScoreBILD(Int4 key, Int4 site, set_typ FGSet, set_typ BGSet, dms_typ *dms,FILE *ofp)
// Compute the BILD score for this FG/BG configuration.
// PdbSeq(sE)...pdb sequences...
//  set_typ PdbSet=mcs->RtnPdbSet(Int4 g);
{
	Int4	m=key,i=site;
	Int4	n_fg,n_bg,j,n,x,y,z,*Parent=0;
	double  sum,dd,d,Dd,FG,BG,DD,dD,tmp=0.0,sumFG,sumBG,diff,Diff;
	UInt4   WtCntsFG[50],WtCntsBG[50],WtCntsBoth[50],TF,TB;
	che_typ **che=mcs->RtnChe( );   
	if(che[m]->PttrnPos(i)) return 0;
	if(che[1]->PttrnPos(i)) return 0;
	double	WtFact=(double)che[m]->GetWtFactor(),KeyFreq=0;
	unsigned char KeyRes = ResSeq(site,che[m]->KeyE());
	static Int4 calls=0;
	if(calls==0 && ofp) che[m]->PutBestPatterns(ofp,TRUE); calls++;
	//che[m]->PutPattern(stdout);	// there will be no pattern at 'site'.

#if 0	// Allow larger sets...
	bpps_typ *bpps=che[m]->BPPS();
	rst_typ *rst=new rst_typ('G');
	rst->Put(stdout,AB); exit(1);
	rst_typ *RST=bpps->RtnRST();
	sst_typ **xsst=RST->LegalResSets( );
	sst_typ *legal_sst=xsst[KeyRes];
#endif
#if 0
	Int4	pernats=1000;
	char	dms_mode='F',wt_factor=100;	// 
	a_type	ab=MkAlpha("XCN",DUMMY_MTRX);
	dms_typ	*dms0=new dms_typ(wt_factor,pernats,dms_mode,ab);
#endif
	
	UInt4	*KeyWtCntsFG=che[m]->FGResWtsBILD(i,TF),*KeyWtCntsBG=che[m]->BGResWtsBILD(i,TB);
	KeyFreq=(double)KeyWtCntsFG[KeyRes]/(double)TF;
	for(x=0; x <= nAlpha(AB); x++){ WtCntsBoth[x] = KeyWtCntsFG[x] + KeyWtCntsBG[x]; }
#if 1
	Diff = calc_binary_bild(WtCntsBoth,WtFact,KeyRes,AB);
	Diff -= calc_binary_bild(KeyWtCntsFG,WtFact,KeyRes,AB);
	Diff -= calc_binary_bild(KeyWtCntsBG,WtFact,KeyRes,AB);
#else
	Diff = dms->bild(WtCntsBoth) - (dms->bild(KeyWtCntsFG) + dms->bild(KeyWtCntsBG));
#endif
	// Diff = baseline for comparison...this should be negative...
	for(x=0; x <= nAlpha(AB); x++){ WtCntsFG[x]=KeyWtCntsFG[x]; WtCntsBG[x]=0; }
	ClearSet(FGSet); AddSet(m,FGSet); ClearSet(BGSet); Hpt->IsTree(Parent);
	free(KeyWtCntsFG); free(KeyWtCntsBG);
	// fprintf(stdout,"%d.Freq(%c) = %.3f\n",i,AlphaChar(KeyRes,AB),KeyFreq);
	if(KeyFreq < 0.50) return 0;

	for(n_fg=1,n_bg=0,sumFG=sumBG=0.0,n=2; n < Hpt->NumSets(); n++){
	        if(n == m) continue; 	// examine siblings only...
	        if(Parent[n] != Parent[m]) continue; 	// examine siblings only...
	        // if(MemberSet(n,keySet) || MemberSet(n,FGSet)) continue; 
		// see whether the sibling is similar to key node...
	        UInt4  *WtCntsTst=che[n]->FGResWtsBILD(i,TF);
		for(x=0; x <= nAlpha(AB); x++){ WtCntsBoth[x] = WtCntsFG[x] + WtCntsTst[x]; }
#if 1
		diff = calc_binary_bild(WtCntsBoth,WtFact,KeyRes,AB);
		diff -= calc_binary_bild(WtCntsFG,WtFact,KeyRes,AB);
		diff -= calc_binary_bild(WtCntsTst,WtFact,KeyRes,AB);
#else
		diff = dms->bild(WtCntsBoth) - (dms->bild(WtCntsFG) + dms->bild(WtCntsTst));
#endif
		// fprintf(stdout,"%d(%d vs %d): diff = %.3f\n",i,m,n,diff);
		if(diff <= 0.0){	// then put into background...
		   AddSet(n,BGSet); sumBG += -diff; n_bg++; 
	      	   for(x=0; x <= nAlpha(AB); x++){ WtCntsBG[x] += WtCntsTst[x]; }
		} else {		// else put  into foreground.
		   AddSet(n,FGSet); sumFG += diff; n_fg++;
	      	   for(x=0; x <= nAlpha(AB); x++){ WtCntsFG[x] += WtCntsTst[x]; }
		} free(WtCntsTst);
	} if(Parent) free(Parent);
	// fprintf(stdout,"%d: n_fg=%d; n_bg=%d\n",site,n_fg,n_bg);
	if(n_bg > 1 && n_fg > 1){
	  for(x=0; x <= nAlpha(AB); x++){ WtCntsBoth[x] = WtCntsFG[x] + WtCntsBG[x]; }
#if 1
	  diff = calc_binary_bild(WtCntsBoth,WtFact,KeyRes,AB);
	  diff -= calc_binary_bild(WtCntsFG,WtFact,KeyRes,AB);
	  diff -= calc_binary_bild(WtCntsBG,WtFact,KeyRes,AB);
#else
	  diff = dms->bild(WtCntsBoth) - (dms->bild(WtCntsFG) + dms->bild(WtCntsBG));
#endif
	  tmp = Diff - diff;	// negative difference is used...
#if 1	// come up with key pattern...
	  if(ofp){
	    fprintf(ofp,"========== %c%d: ==========\n",AlphaChar(KeyRes,AB),site);
	    fprintf(ofp,"  FG:"); PutSet(ofp,FGSet);
	    fprintf(ofp,"  BG:"); PutSet(ofp,BGSet);
	  }

	  double totalFG=0.0,totalBG=0.0;
	  for(x=0; x <= nAlpha(AB); x++){ totalFG+=(double)WtCntsFG[x]; totalBG+=(double)WtCntsBG[x]; } 
	  for(x=1; x <= nAlpha(AB); x++){
		double re,p,q;
		p = (double) WtCntsFG[x]/totalFG; q = (double) WtCntsBG[x]/totalBG; 
		if(p > 0.05 && p > q){
		   re=p*log(p/q);
		   if(ofp) fprintf(ofp,"%d.RE(%c}=%.2f; p=%.3f; q=%.3f; cntsFG=%.1f; cntsBG=%1.f.\n",
			  i,AlphaChar(x,AB),re,p,q,
			  (double)WtCntsFG[x]/WtFact,(double)WtCntsBG[x]/WtFact);
		}
	  } totalFG=totalFG/WtFact;  totalBG=totalBG/WtFact; 
#if 0
	  //  TotWtMN=TotWtMN/WtFact;
	  fprintf(stdout,"%d.FG residues=%.1f (%.1f%c), deletions=%.1f; BG residues=%.1f\n", i,
		totalFG,100.0*totalFG/TotWtMN,'%',TotWtMN-totalFG,totalBG);
#endif
#endif
	  // Also want to check deletions in original foreground set...
#if 0
	  fprintf(stdout,"%d.%d %s BILD scores: (FG=%.1f; BG =%.1f; Both=%.1f;) %.1f-%.1f=%.1f.\n",
			i,m,Hpt->SetName(m),d,dd,DD,Diff,diff,tmp);
#else
	  if(ofp) fprintf(ofp,"%d.%d %s delta BILD score: %.1f nats  (p = %.3g).\n",
			i,m,Hpt->SetName(m),tmp,exp(-tmp)/(1.0 + exp(-tmp)));
	  // logistic function used at the end to obtain a p-value (i.e., p).
#endif
	} return tmp;
} 

Int4 	omc_typ::PutCrossConserved(Int4 key, cmh_typ *cmh)
{
	// 1. Obtain information...
	assert(key > 1 && key < Hpt->NumSets());  // avoid root and reject nodes.
	char str[200]; assert(outfilename != 0); 
	Int4	i,j,n,x,r,new_node=0,*Parent; Hpt->IsTree(Parent);
	hpt_typ *Hpt=mcs->GetHpt();
	Int4	id,*chld,dpth=Hpt->NodeDepth(key); NEW(chld,Hpt->NumSets()+3,Int4);
	for(i=0, j=key; j > 0; j=Parent[j]){ i=Parent[j]; chld[i]=j; }

	// 2. Obtain hpt...
	hpt_typ *nhpt=0,*xhpt= new hpt_typ(20000);
	// nhpt = Hpt->MkLineageTree(key,new_node); // nhpt->Put(stdout); 
	// delete nhpt;
	// xhpt->Put(stderr);
	for(j=1,i=1; chld[i] > 0; j++,i=chld[i]){
	   id=Hpt->ItoSetID(chld[i]); nhpt = xhpt->AddLeaf(j,id); // nhpt->Put(stderr);
	   delete xhpt; xhpt=nhpt;
	} free(chld);

	// 3. Pop off hits...
	Int4 ID=0;
	keytyp	k;
	set_typ FGSet=0,BGSet=0;
	sst_typ	*isst;
#if 0
     UInt4 itot,jtot,idnt,jdnt,isim,jsim;
     sst_typ *xsst,*ysst,zsst,usst;
     for(i=1; i <= cmh->Size(); i++){
	xsst=cmh->RtnSST(i);
	if(xsst == 0) continue; 
        for(j=i+1; j <= cmh->Size(); j++){
	   ysst=cmh->RtnSST(j);
	   if(ysst == 0) continue; 
	   fprintf(stdout,"%d.%d: ",i,j);
	   itot=idnt=isim=0; jtot=jdnt=jsim=0;
	   for(r=1; r <= LengthCMSA(1,TrueMainCMA); r++){ 
		zsst=IntersectSset(xsst[r],ysst[r]);
		if(xsst[r] > 0) itot++; if(ysst[r] > 0) jtot++;
		if(xsst[r] > 0 && xsst[r] == ysst[r]){ idnt++; jdnt++; }
		else if(zsst > 0){ isim++; jsim++; }
		if(zsst != 0){
		  usst=UnionSset(xsst[r],ysst[r]);
		  PutSST(stdout,usst,AB); fprintf(stdout,"%d ",r); 
		}
	   } fprintf(stdout,"\n  itot=%d; jtot=%d; idnt=%d; sim=%d\n",itot,jtot,idnt,isim);
	}
     } exit(1);
#endif
     while((n=cmh->DeleteMax(FGSet,BGSet,k,isst)) != 0){
	ID++;
	sprintf(str,"%s_%d_%d",outfilename,key,ID);
	FILE *ofp=0;
	// ofp=open_file(str,".hpt","w"); xhpt->Put(ofp); fclose(ofp);
	fprintf(stdout,"key = %d: %d hits; LPR = %.3f\n",key,n,k);
	for(r=1,x=0; r <= LengthCMSA(1,TrueMainCMA); r++){ 
               if(isst[r]){ x++; PutSST(stdout,isst[r],AB); fprintf(stdout,"%d ",r); }
        } fprintf(stdout," (%d pttrn pos.)\n",x); 
	fprintf(stdout,"  FG:"); PutSet(stdout,FGSet);
	fprintf(stdout,"  BG:"); PutSet(stdout,BGSet);

	// 4. Obtain sets...
	set_typ	*xset; NEW(xset,xhpt->NumSets() + 3, set_typ);
	Int4	p,rjct=xhpt->NumSets(),Rjct=Hpt->NumSets();
        xset[rjct]=MakeSet(mcs->GetSetSize());  mcs->CopySubTreeSeqSet(Rjct,xset[rjct]);
	for(j=rjct-1,i=key; j > 0 && i > 0; j--,i=Parent[i]){
            xset[j]=MakeSet(mcs->GetSetSize());  mcs->CopySubTreeSeqSet(i,xset[j]);
	    for(Int4 z=j+1; z < rjct; z++){
		IntersectNotSet(xset[j], xset[z]); // xset[j] = xset[j] && not xset[z].
	    }
	}
	// 4b. Move selected BG sets to FG.
	set_typ iset=0;
	Int4 keyX=rjct-1;
	Int4 Prnt=keyX-1;
	for(p=Parent[key],i=2; i < Hpt->NumSets(); i++){
	    if(i == key) continue;
	    if(Parent[i] == p && MemberSet(i,FGSet)){
        	// fprintf(stdout,"moving %d to FG)\n",i); 
        	iset=MakeSet(mcs->GetSetSize());  mcs->CopySubTreeSeqSet(i,iset);
		UnionSet(xset[keyX],iset); // add iset to target set.
		IntersectNotSet(xset[Prnt], iset); // xset[Prnt] = xset[Prnt] && not iset.
		NilSet(iset);
	    }
	}
	// ofp=open_file(str,".sets","w"); WriteSets(ofp,xhpt->NumSets(),xset); fclose(ofp);

	// 5. Obtain sma...
	ofp=0; // ofp=open_file(str,".sma","w");
	for(i=1; i < rjct; i++){ 
		cma_typ	xsma=GetBestCsqAsCMSA(xset[i],TrueMainCMA);
		ReNameCMSA(xhpt->SetName(i),xsma);
		if(ofp) PutCMSA(ofp,xsma); TotalNilCMSA(xsma);
	} if(ofp) fclose(ofp);
        // in_sma[1]=MakeConsensusCMSA(TrueMainCMA); RenameCMSA("Set1",in_sma[1]);

	// 6. Free memory...
	for(i=1; i<= rjct; i++) NilSet(xset[i]); free(xset);
	NilSet(FGSet); NilSet(BGSet); free(isst);
    } free(Parent); delete xhpt;
    return ID;
}

#include "blosum62.h"

double	*omc_typ::FindCrossConserved(Int4 key, FILE *ofp)
//****************** Compare BPPS with BILD scores. *******************
// This version uses Altschul's approach: sum BILD - total BILD.
// Use bild scores & BPPS to find inconsistently conserved residue positions.
// 
{
	Int4	heapsize=XC_size;
	hpt_typ *Hpt=mcs->GetHpt();
	if(!(key > 1 && key < Hpt->NumSets())) print_error("-XC option input error");
			// ignore root & reject nodes.
	// if(Hpt->TypeOfSet(key) != '!') return 0; // ignore non-leaf nodes.

	set_typ SibSet=Hpt->MkSiblingSet(key);
	if(CardSet(SibSet) < 1){ NilSet(SibSet); return 0; }
		// print_error("\n!!!  Input error: selected node lacks siblings!!!\n");
	if(ofp) PutSet(ofp,SibSet); NilSet(SibSet);
	set_typ	keySet=Hpt->DescendantSet(key); 
	AddSet(key,keySet);
	set_typ FGSet=CopySet(keySet),BGSet=CopySet(keySet); 
	// set_typ otherSet=CopySet(keySet),TstSet=CopySet(keySet);

        if(ofp) mcs->PutHyperPartition(ofp); 
	double	Max=(double)-INT4_MAX,Min=(double)INT4_MAX; 
	double	dd,d,Dd,DD,dD;
	Int4	i,j,k,n,x,pernats=1000;
	char	dms_mode='F',wt_factor=100;	// 
	dms_typ	*dms=new dms_typ(wt_factor,pernats,dms_mode,AB);
	dh_type dH=dheap(mcs->RtnLengthMainCMSA()+3,4);
	double	*RtnValue; NEW(RtnValue,mcs->RtnLengthMainCMSA()+3,double);
	cmh_typ *cmh=new cmh_typ(heapsize);
	// 1. Find cross conserved positions in the alignment.
	for(i=1; i <= mcs->RtnLengthMainCMSA(); i++){
	    DD=this->CrossScoreBILD(key,i,FGSet,BGSet,dms,ofp); // Get FG and BG sets for key node.
	    // if(DD > 0) this->ComputeCrossScore(key,i,FGSet,BGSet,cmh);
	    if(DD > 0){ insrtHeap(i,-(keytyp)DD,dH); RtnValue[i]=DD; }
	    else RtnValue[i]=-1.0;
	    if(DD > Max) Max=DD;
	    if(DD < Min) Min=DD;
	} if(ofp) fprintf(ofp,"\n\n"); RtnValue[0]=Max;
	h_type	HG=Histogram("Relative BILD scores",-50,100,5);
	set_typ SiteSet=MakeSet(TotalLenCMSA(TrueMainCMA)+5); ClearSet(SiteSet);
	for(n=1; !emptyHeap(dH); n++){
		dd=-minkeyHeap(dH); assert((i=delminHeap(dH)) != 0);
		Dd = 100*(dd - Min)/(Max - Min); IncdHist(Dd, HG);
		if(n <= heapsize && ofp) fprintf(ofp,"%d.col %d: %.1f average dF (%.1f)\n",n,i,dd,Dd);
		if(n <= heapsize) AddSet(i,SiteSet);
	} if(ofp) fprintf(ofp,"\n"); Nildheap(dH);
#if 1
	// 1. Add current pattern sites as well...
	sst_typ *xsst=mcs->RtnCopyOfSST(key);
	set_typ LineSet=Hpt->MkLineageSet(key); DeleteSet(key,LineSet);
	sst_typ **rsst=0; 
	if(XCuseL){
	  Int4 ns=0;
	  NEWP(rsst, Hpt->NumSets()+5,sst_typ);
	  for(i=1; i < Hpt->NumSets(); i++){ 
#if 1
	    if(MemberSet(i,LineSet)) rsst[i]=mcs->RtnCopyOfSST(i);
#else	// merging only a few...
	    if(MemberSet(i,LineSet)){
		ns++;
		if(ns >= 3) rsst[i]=mcs->RtnCopyOfSST(i);
	    }
#endif
	  } if(ofp) PutSet(ofp,LineSet);
        }
	mcs->RtnCopyOfSST(1);
	che_typ **che=mcs->RtnChe( );   
	e_type KeySeq=che[key]->KeyE();
	char *SubTyp; NEW(SubTyp,TotalLenCMSA(TrueMainCMA)+5,char);
	for(j=1; j <= TotalLenCMSA(TrueMainCMA); j++){
	   if(XCuseX && xsst[j]==0 && MemberSet(j,SiteSet)){	 // 'X' = cross conserved.
		unsigned char c=ResSeq(j,KeySeq);
		SubTyp[j]='X';
		if(c > 0){ xsst[j]=UnionSset(xsst[j],SsetLet(c)); }
	   } else if(xsst[j]){				// 'F' = family-specific
		AddSet(j,SiteSet);
		if(ofp){ PutSST(ofp,xsst[j],AB); fprintf(ofp,"%d,",j); }
		SubTyp[j]='F';
	   } else if(XCuseL){		// 'L' = lineage conserved.
	      for(i=Hpt->NumSets(); i > 0; i--){
		if(rsst[i] && rsst[i][j]){
		  AddSet(j,SiteSet); xsst[j]=rsst[i][j]; SubTyp[j]='L'; break; 
		}
	      }
	   }
	} if(ofp) fprintf(ofp,"\n");
	// 2. Compute BILD scores, Relative entropy and entropy...
	if(ofp) fprintf(ofp,"======= Cross conserved(X);Family(F); Lineage(L) ======\n");
	dH=dheap(TotalLenCMSA(TrueMainCMA)+3,4);
	UInt4  **WtCntsFG=che[key]->GetResWtsFG();
	UInt4   *TotalWtCntsFG=che[key]->GetTotalResWtFG();
	Int4	start=1;
	for(j=1; j <= TotalLenCMSA(TrueMainCMA); j++){
	    if(!MemberSet(j,SiteSet)) continue;
	    if(xsst[j]==0) continue;
	    double p,q,pp,qq,totalFG=0.0;
	    UInt4 *wtCntsFG=WtCntsFG[j];
	    for(x=start; x <= nAlpha(AB); x++){ totalFG+=(double)wtCntsFG[x]; } 
	    if(TotalWtCntsFG[0] >= totalFG) continue;  // if more gaps than residues then skip.
	    for(dd=Dd=p=0,x=start;  x <= nAlpha(AB); x++){
	       if(MemSset(x,xsst[j])) p += (double)wtCntsFG[x]/totalFG; 
	       pp=(double)wtCntsFG[x]/totalFG;
	       qq=blosum62freq[x];
	       if(pp > 0){ dd += pp*log(pp/0.05); Dd += pp*log(pp/qq); }
	    } q = 1.0 - p;
	    dD = 0.0; if(p > 0) dD -= p*log(p); if(q > 0) dD -= q*log(q);  
	    d=dms->bild(WtCntsFG[j]); 
#if 0	// this doesn't work too well...use BG to calculate RE?
	    DD=-dD;	// functional vs non-functional entropy.
#elif 0	// works okay...
	    DD=d;	// BILD scores.
#elif 0	// not the best...
	    DD=dd;	// 20-residue relative entropy with flat null model.
#else	// this might be the best...
	    DD=Dd;	// 20-residue relative entropy.
#endif
	    if(SubTyp[j]=='L') DD= DD*XC_wt;	// down weight these...appears to help (anecdotally).
	    if(SubTyp[j]=='X') DD= DD*XC_wt;	//
	    insrtHeap(j,-(keytyp)DD,dH); 
	    if(ofp){ PutSST(ofp,xsst[j],AB); fprintf(stdout,"%d=%.3f %c\n",j,DD,SubTyp[j]); }
	}
	Int4 *Site; NEW(Site,CardSet(SiteSet)+5,Int4);
	for(n=0; !emptyHeap(dH); ){
	    dd=-minkeyHeap(dH); assert((j=delminHeap(dH)) != 0);
	    n++; Site[n]=j; 
	} Nildheap(dH);
	if(rsst){
	  for(i=1; i < Hpt->NumSets(); i++) if(rsst[i]) free(rsst[i]);
	  free(rsst);
	} free(xsst); delete dms; free(SubTyp);
#endif
#if 1
	Int4 N=NumSeqsCMSA(TrueMainCMA);
	set_typ pdbSet=mcs->RtnPdbSet(key);
	for(Int4 gg=2; gg < Hpt->NumSets(); gg++){
	   if(MemberSet(gg,keySet)){
		set_typ tmpSet=mcs->RtnPdbSet(gg);
		UnionSet(pdbSet,tmpSet); NilSet(tmpSet);
	   }
	}
	if(ofp){
	   PutSet(ofp,pdbSet); fprintf(ofp,"...%d pdb seqs\n",CardSet(pdbSet)); 
	   PutSet(ofp,keySet); PutSet(ofp,SiteSet);
	}
	e_type sE=0;
	for(i=1; i <= N; i++){
	   if(MemberSet(i,pdbSet)){ 
		Int4 s1,pos[4];
		if(ofp) fprintf(ofp,"%c=",XC_mode);
		for(k=1; Site[k] > 0; k++){
		   j = Site[k];
		   assert(PosSiteCMSA(1,i,pos,TrueMainCMA));
                   if(IsDeletedCMSA(i,pos[1]+j-1,TrueMainCMA)) continue;
		   s1=TruePosCMSA(i,j,TrueMainCMA) + OffSetSeq(TrueSeqCMSA(i,TrueMainCMA));
		   // fprintf(stdout,"%d(%d),",s1,pos[1]+j-1);
		   if(ofp) fprintf(ofp,"%d,",s1);
		} if(ofp) fprintf(ofp,"\n");
		sE=TrueSeqCMSA(i,TrueMainCMA); if(ofp) PutSeq(ofp,sE,AB);
	   }
	} NilSet(pdbSet); NilSet(SiteSet);
#endif
	if(ofp){ fprintf(ofp,"Max=%.3f; Min=%.3f\n",Max,Min); PutHist(ofp,60,HG); fflush(ofp); }
	NilHist(HG);
	//  this->PutCrossConserved(key,cmh);
	delete cmh; NilSet(FGSet); NilSet(BGSet);
	NilSet(keySet); NilSet(LineSet); 
	return RtnValue;
}

