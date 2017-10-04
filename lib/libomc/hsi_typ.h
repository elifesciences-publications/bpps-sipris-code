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

#if !defined (HSI_TYP)
#define HSI_TYP

#include "stdinc.h"
#include "afnio.h"
#include "btn_typ.h"
#include "random.h"
#include "dheap.h"
#include "mheap.h"
#include "set_typ.h"
#include "dsets.h"
#include "probability.h"
#include "lpr_typ.h"
#include "sqd_typ.h"
#include "histogram.h"

#include "sset.h"
#include "hpt_typ.h"

/*****************************************************************************
  Create a binary tree representing a N-ary tree with these objects...
 *****************************************************************************/

class hsi_typ  { 	// hierarchical subgroup information type 
public:
		hsi_typ(set_typ fg,set_typ bg,sst_typ *sst,double Lpr,Int4 card_fg,
				Int4 card_bg,lpr_typ *lpr)
		   { SetSize=SetN(fg); L=lpr; Init(); Store(fg,bg,sst,Lpr,card_fg,card_bg); }
		~hsi_typ( ){ Free(); }
	void	Store(set_typ fg,set_typ bg,sst_typ *sst,double Lpr,Int4 card_fg, Int4 card_bg){
			CopySet(SetFG,fg); CopySet(SetBG,bg);  CardFG=card_fg; CardBG=card_bg;
			if(qsst) free(qsst); qsst=sst; LPR=Lpr;  
		}
	double	Lpr(){ return LPR; }
	void	Put(FILE *fp,set_typ *Set,Int4 num);
	void	Put(FILE *fp);
	void	PutTree(FILE *fp,set_typ *Set,Int4 num){
			PutTree(fp,Set,num,0);
		}
	void	PutTree(FILE *fp,set_typ *Set,Int4 num,Int4 depth);
	double	CrossCompare(FILE *fp,hsi_typ *other,set_typ *Set,Int4 num,double &tLpr,double &oLpr,
					double &tWtCntFG, double &oWtCntFG);
	Int4	IntersectionFGs(hsi_typ *hsi){ return CardInterSet(this->SetFG,hsi->SetFG); }
	Int4	CardFGset( ){ return CardSet(SetFG); }
	sst_typ	*GetSST(){ return qsst; }
	hsi_typ	*Child(){ return child; }
	hsi_typ	*Sibling(){ return sibling; }
	hsi_typ	*GetParent(){ return parent; }
	hsi_typ *Link(hsi_typ *X){ assert(this->child==0); this->child=X; return this; }
	double	TestSubset(hsi_typ *hsi,set_typ *Set,Int4 num);
	sst_typ	*TestSubset(hsi_typ *hsi,set_typ *Set,Int4 num, double &Lpr);
	hsi_typ	*PlaceIntoTree(hsi_typ *other,set_typ *Set,Int4 num, double min_nats);
	// void	AddChild(hsi_typ *other);
	set_typ	RtnSetFG(){ return SetFG; }
	set_typ	RtnSetBG(){ return SetBG; }
	BooLean	IsTheSame(hsi_typ *hsi);
	BooLean	IsSubSet(hsi_typ *hsi);
	Int4	NumChildren(){
		    hsi_typ *X; Int4 n=0;
		    for(X=this->child; X!= 0; X=X->sibling) { n++; } return n;
		}
private:
	hsi_typ	*sibling;	// points to a sibling of this node...
				//    which may point to other siblings with same parent.
	hsi_typ	*child;		// poings to this node's child
				//    which may point to other children (i.e., follow sibling list).
	set_typ	GetSetFG(Int4 num, set_typ *Set);	// large set of seqs in foreground
	set_typ	GetSetBG(Int4 num, set_typ *Set);	
	double	GetWtCntsFG(set_typ FG, set_typ BG, sst_typ *xsst,double &WtCntsFG,double &WtCntsBG);
	void	Init();
	void	Free();
	sst_typ	*qsst;
	double	LPR;
	lpr_typ	*L;
	Int4	SetSize;
	set_typ	SetFG,SetBG;	// small set of nodes in foreground and background.
	Int4	CardFG,CardBG;
	// N-ary tree structure...implemented as a binary tree.
	hsi_typ	*parent;	// points to the parent of this node.
};

class hmh_typ {		// hmh_typ = hsi Min/Max heap type
public:
		hmh_typ(){ print_error("not allowed"); }
		hmh_typ(Int4 h){ Init(h); }
		~hmh_typ(){ Free(); }
	Int4    Insert(hsi_typ *hsi,double key);
	hsi_typ	*DelMin(double *key, Int4 *Item);
	hsi_typ	*DelMax(double *key, Int4 *Item);
	BooLean IsSuperInHeap(hsi_typ *hsi);
private:
	mh_type mheap;
	Int4	hpsz;
	hsi_typ	**HSI;
	void	Init(Int4 hs);
	void	Free();
};

#endif

