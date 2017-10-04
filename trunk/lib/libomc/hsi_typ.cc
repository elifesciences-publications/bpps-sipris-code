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

#include "hsi_typ.h"

void	hsi_typ::Put(FILE *fp,set_typ *Set,Int4 num)
{
h_type	HG=Histogram("single node lprs",-10000,10000,10);
	double max_lpr=-999999999999.0,min_lpr=DBL_MAX,d;
	Int4   min_id,max_id;
	set_typ BG=MakeSet(SetN(Set[1])); ClearSet(BG);
	for(Int4 c=1; c <= num; c++){ if(MemberSet(c,SetBG)) UnionSet(BG,Set[c]); }
	for(Int4 c=1; c <= num; c++){
	   if(MemberSet(c,SetFG)){	
		d=L->CalcSetvsPttrnLPR(0,Set[c],BG,qsst,FALSE,'M');
#if 1
		double dummy,WtCntsFG,WtCntsBG;
		L->WtCardFG_BG_Sets(WtCntsFG,dummy);
		// WARNING: need to have called CalcSetvsPttrnLPR() to set FG/BG counts!!
		Int4 wt_cnt=(Int4)ceil(WtCntsFG); 
		if(wt_cnt > 0) IncdMHist(d,wt_cnt,HG);
#else
		IncdHist(d,HG);
#endif
		// fprintf(stderr," Set[%d] single node LPR = %.3f\n",c,d);
		if(d > max_lpr){ max_lpr=d; max_id=c; }
		if(d < min_lpr){ min_lpr=d; min_id=c; }
	   }
	}
	fprintf(stderr,"Single node FG LPR: %.2f(%d) ... %.2f(%d)\n",
			min_lpr,CardSet(Set[min_id]),max_lpr,CardSet(Set[max_id]));
	double sd=sqrt(VarianceHist(HG)),m=MeanHist(HG);
	fprintf(stderr,"Mean +/- 1sd = %.2f ... %.2f ... %.2f [%.2f]\n",m-sd,m,m+sd,sd/m);
NilHist(HG);
	if(CardSet(SetBG) > 1){
	      max_lpr=-999999999999.0; min_lpr=DBL_MAX;
	      for(Int4 c=1; c <= num; c++){
		if(MemberSet(c,SetBG)){	
		   d=L->CalcSetvsPttrnLPR(0,Set[c],BG,qsst,TRUE,'M');
		   // fprintf(stderr," Set[%d] single node LPR = %.3f\n",c,d);
		   if(d > max_lpr){ max_lpr=d; max_id=c; }
		   if(d < min_lpr){ min_lpr=d; min_id=c; }
		}
	      } fprintf(stderr,"Single node BG LPR: %.2f(%d) ... %.2f(%d)\n",
			min_lpr,CardSet(Set[min_id]),max_lpr,CardSet(Set[max_id]));
	} NilSet(BG);
	fprintf(fp,"--------------------------------------------------\n");
}

void	hsi_typ::Put(FILE *fp)
{
        fprintf(fp,"FG (%d)= Sets: ",CardFG); PutSet(fp,SetFG);
        fprintf(fp,"------  optimum subLPR=%.3f  ------\n",LPR);
        // fprintf(stderr,"Single node LPR: %.3f ... %.3f\n",min_lpr,max_lpr);
        fprintf(fp,"BG (%d)= Sets: ",CardBG); PutSet(fp,SetBG);
	fprintf(fp,"--------------------------------------------------\n");
}

set_typ hsi_typ::GetSetBG(Int4 num, set_typ *Set){
        set_typ X=MakeSet(SetN(Set[1])); ClearSet(X);
        for(Int4 c=1; c <= num; c++){ if(MemberSet(c,SetBG)) UnionSet(X,Set[c]); }
        return X;
}

set_typ hsi_typ::GetSetFG(Int4 num, set_typ *Set){
        set_typ X=MakeSet(SetN(Set[1])); ClearSet(X);
        for(Int4 c=1; c <= num; c++){ if(MemberSet(c,SetFG)) UnionSet(X,Set[c]); }
        return X;
}

double  hsi_typ::GetWtCntsFG(set_typ FG, set_typ BG, sst_typ *xsst,double &WtCntsFG,double &WtCntsBG){
        double d=L->CalcSetvsPttrnLPR(0,FG,BG,xsst,FALSE,'M');
        // WARNING: requires calling CalcSetvsPttrnLPR() to set FG/BG counts!!
        L->WtCardFG_BG_Sets(WtCntsFG,WtCntsBG); return d;
}

double	hsi_typ::CrossCompare(FILE *fp,hsi_typ *other,set_typ *Set,Int4 num,double &tLpr,
			double &oLpr, double &tWtCntFG, double &oWtCntFG)
{
	double	dummy;
	set_typ FG=this->GetSetFG(num,Set),BG=this->GetSetBG(num,Set);
	set_typ oFG=other->GetSetFG(num,Set),oBG=other->GetSetBG(num,Set);

	tLpr=this->GetWtCntsFG(FG,BG,other->qsst,tWtCntFG,dummy);

	oLpr=other->GetWtCntsFG(oFG,oBG,this->qsst,oWtCntFG,dummy);

	double ratio=0.0;
	if(tLpr > 0.0 && oLpr > 0.0){ 
		ratio = (tLpr*oWtCntFG)/(oLpr*tWtCntFG);
		if(ratio > 1.0) ratio = 1.0/ratio;
	}
	if(fp){
            fprintf(fp,"tFG (%d:%.1f)= Sets: ",this->CardFG,tWtCntFG); PutSet(fp,this->SetFG);
            fprintf(fp,"------  cross LPR=%.2f (%.2f) ------\n",tLpr,this->LPR);
            fprintf(fp,"oFG (%d:%.1f)= Sets: ",other->CardFG,oWtCntFG); PutSet(fp,other->SetFG);
            fprintf(fp,"------  cross LPR=%.2f (%.2f) ------\n\n",oLpr,other->LPR);
            // fprintf(fp,"------  ratio = %.2f ------\n", ratio);
	}
	NilSet(FG); NilSet(oFG); NilSet(BG); NilSet(oBG);
	return ratio;
}

#if 0
void	hsi_typ::AddChild(hsi_typ *other)
{ if(this->child == 0){ this->child = other; other->parent = this; } }
#endif

double  hsi_typ::TestSubset(hsi_typ *other,set_typ *Set,Int4 num)
{ double d; sst_typ *xsst=TestSubset(other,Set,num,d); free(xsst); return d; }

sst_typ	*hsi_typ::TestSubset(hsi_typ *other,set_typ *Set,Int4 num,double &Lpr)
// Test whether 
{
	sst_typ	*xsst;
	Lpr=0;
	set_typ FG=other->GetSetFG(num,Set),BG=this->GetSetFG(num,Set);
	IntersectNotSet(BG,FG); // modifies BG to equal BG intersect not FG */
        xsst=L->GetOptPttrnLPR(0,FG,BG,FALSE,Lpr,'M'); NilSet(FG); NilSet(BG); 
	return xsst;
}

hsi_typ	*hsi_typ::PlaceIntoTree(hsi_typ *other,set_typ *Set,Int4 num, double min_nats)
// 'this' is assumed to be the 'root' though with siblings...
// Find the best place in which to put 'other' within the tree, otherwise return 0.
{
	FILE *fp=stderr;
	double	tLpr,oLpr,Lpr,min_oLpr,tWtCnt,oWtCnt,ratio,d;
	Int4	CardIS,xC,oC;
	hsi_typ	*X=0;
	sst_typ	*xsst;
	xC=this->CardFGset(); oC=other->CardFGset();
	CardIS=this->IntersectionFGs(other);
	if(CardIS > 1){	// The sets intersect...
	   // Case A:
	   if(oC == CardIS){		// 'other' is a proper subset of 'this'...
	      // Then o < t --> o could also be a subset of a t-child; recursion will check this case.
	      ratio=this->CrossCompare(stderr,other,Set,num,tLpr,oLpr,tWtCnt,oWtCnt);
	      if(tLpr < 0.0 && oLpr > 0.0){		// then this may be a usable superset.
		xsst=this->TestSubset(other,Set,num,d);
		if(d > min_nats){			
		   if(this->child == 0){
			free(other->qsst); other->qsst=xsst; 
			IntersectNotSet(this->SetFG, other->SetFG,other->SetBG);
			// ^ modifies other->SetBG to equal this->SetFG intersect not other->SetFG.
			other->CardBG=CardSet(other->SetBG);	// need to update this as well.
			this->child = other; other->parent = this; return this; 
		  } else {	// check the children of 'this' both for sub-supersets and siblings. 
			free(xsst); 
			if(this->child->PlaceIntoTree(other,Set,num, min_nats)) return this;
			else return 0;
			// return this->child->PlaceIntoTree(other,Set,num, min_nats);
			// ^ this will put 'other' at the end of the child list if no conflicts...
			// will return NULL if conflict found...
		   }
		} else { free(xsst); return 0; }
	      } else return 0;
	   // Case B:
	   } else if(xC == CardIS){	// 'this'  a subset of 'other'... need to return 'other'.
	      // then o > t --> o not a subset of children of t (recursion will only check t-siblings).
		// then make entire tree a subset of 'other'...
	      ratio=other->CrossCompare(stderr,this,Set,num,tLpr,oLpr,tWtCnt,oWtCnt);
	      if(tLpr < 0.0 && oLpr > 0.0){		// then this may be a usable superset.
		d=other->TestSubset(this,Set,num);
		if(d > min_nats){			
		   assert(other->child == 0);	// siblings of this can't be a subset...
		   other->child = this; this->parent = other; return other;
		} else return 0;
	      } else return 0;
	   // Case C:
	   } else {	// there is a conflict!
	      // if(this->sibling != 0) return this->sibling->PlaceIntoTree(other,Set,num, min_nats);
	      return 0;
	   }
	} else if(CardIS == 1){ return 0;	// i.e., This is a conflict! Always reject.
	} else {				// i.e., This is okay, but need to check other children.
		if(this->sibling == 0){ this->sibling=other; other->parent=this->parent; return this; }
#if 0
		else return this->sibling->PlaceIntoTree(other,Set,num, min_nats);
#else
		// A sibling of 'this' will not share any nodes with it; but still return 'this' as root.
		if(this->sibling->PlaceIntoTree(other,Set,num, min_nats)) return this; else return 0; 
#endif
	}
}

void	hsi_typ::PutTree(FILE *fp,set_typ *Set,Int4 num,Int4 depth)
// Depth first search of tree rooted at this'.
{
	switch (depth){
	  case 0: fprintf(fp,"*********** Root ***********\n");  break;
	  case 1: fprintf(fp,"=========== child ===========\n");  break;
	  case 2: fprintf(fp,"=========== grandchild ===========\n");  break;
	  case 3: fprintf(fp,"=========== great grandchild ===========\n");  break;
	  default: fprintf(fp,"----------- %dth generation -----------\n",depth);  break;
	}
	this->Put(stderr); this->Put(stderr,Set,num);
	if(child){ child->PutTree(fp,Set,num,depth+1); }
	if(sibling){ sibling->PutTree(fp,Set,num,depth); }
}

void    hsi_typ::Init()
{
	SetFG=MakeSet(SetSize); SetBG=MakeSet(SetSize); qsst=0;
	child=0; parent=0; sibling=0;
}

void    hsi_typ::Free()
{ if(SetFG) NilSet(SetFG); if(SetBG) NilSet(SetBG); if(qsst) free(qsst); }

BooLean hsi_typ::IsSubSet(hsi_typ *hsi)
// Is hsi a subset of this?
{
        if(this->CardFG < hsi->CardFG) return FALSE;
        if(this->CardBG < hsi->CardBG) return FALSE;
        if(CardInterSet(this->SetFG,hsi->SetFG) != CardSet(hsi->SetFG)) return FALSE;
        if(CardInterSet(this->SetBG,hsi->SetBG) != CardSet(hsi->SetBG)) return FALSE;
	return TRUE;
}

BooLean hsi_typ::IsTheSame(hsi_typ *hsi)
{
        if(this->LPR != hsi->LPR) return FALSE;
        if(this->CardFG != hsi->CardFG) return FALSE;
        if(this->CardBG != hsi->CardBG) return FALSE;
        if(CardSet(this->SetFG) != CardSet(hsi->SetFG)) return FALSE;
        if(CardSet(this->SetBG) != CardSet(hsi->SetBG)) return FALSE;
        if(CardInterSet(this->SetFG,hsi->SetFG) != CardSet(hsi->SetFG)) return FALSE;
        if(CardInterSet(this->SetBG,hsi->SetBG) != CardSet(hsi->SetBG)) return FALSE;
        return TRUE;
}

//******************** hsi min-max heap operations ************************
void    hmh_typ::Init(Int4 hs)
{
	hpsz=hs;
	mheap=Mheap(hpsz, 3);
	NEWP(HSI,hpsz +5, hsi_typ);
}

void    hmh_typ::Free()
{
	double	key;
	Int4	item;
	hsi_typ *hsi;

	while((hsi=DelMin(&key,&item)) != NULL){ delete hsi; }
	free(HSI); NilMheap(mheap); 
}

BooLean	hmh_typ::IsSuperInHeap(hsi_typ *hsi)
// Is hsi a subset of an itme in the heap?
{
	for(Int4 i=1; i <= hpsz; i++){
		if(HSI[i] && HSI[i]->IsSubSet(hsi)) return TRUE;
	} return FALSE;
}

Int4	hmh_typ::Insert(hsi_typ *hsi, double key)
// Will insert in the heap if key is > than minimum key currently in the heap.
{
	// throw out duplicates.
	for(Int4 i=1; i <= hpsz; i++){
		if(HSI[i] && hsi->IsTheSame(HSI[i])) return 0;
	}
	Int4 item=InsertMheap(key, mheap);
	if(item==0) return 0;
	if(HSI[item]!= NULL) delete HSI[item];
	HSI[item]=hsi; 
	return item;
}

hsi_typ	*hmh_typ::DelMin(double *key, Int4 *Item)
{
	*key = MinKeyMheap(mheap);
	Int4	item = DelMinMheap(mheap); *Item = item;
	if(item==0) return NULL;
	hsi_typ *hsi = HSI[item]; HSI[item]=0;
	return hsi;
}

hsi_typ	*hmh_typ::DelMax(double *key, Int4 *Item)
{
	*key = MaxKeyMheap(mheap);
	Int4	item = DelMaxMheap(mheap); *Item = item;
	if(item==0) return NULL;
	hsi_typ *hsi= HSI[item]; HSI[item]=NULL;
	return hsi;
}


