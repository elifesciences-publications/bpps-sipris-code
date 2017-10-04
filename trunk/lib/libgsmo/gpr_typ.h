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

#if !defined (_GPR_TYP_)
#define _GPR_TYP_

// #include "gsm_typ.h"
#include "gmb_typ.h"
#include "cma_gmb.h"

class gpr_typ {		// gambit parameters & routines type.
public:
	gpr_typ(char *Name,Int4 Iter, Int4 Stage){
			write=FALSE; init(); name=AllocString(Name); stage=Stage; iter=Iter; 
		}
	gpr_typ(char *Name,Int4 Iter, Int4 Stage,ssx_typ *issx){
			write=FALSE; init(); name=AllocString(Name); stage=Stage; iter=Iter; 
			issx->GetParameters(aa_per_io,aa_per_do,exp_ie,exp_de,prior_wt,dms_mode);
		}
	~gpr_typ(){ free(name); if(gmb) delete gmb; }
	void	SetParameters(Int4 eio,Int4 edo,Int4 ie,Int4 de,double pr_wt,char dm){
			aa_per_io=eio; aa_per_do=edo; exp_ie=ie ;exp_de=de; prior_wt=pr_wt;
			dms_mode=dm;
		}
	void	SetAddRmCols(double bc,double ba,double temp,double af,double rf)
		     { bild_cut=bc; bild_add=ba; Temp=temp; afrq=af; rfrq=rf; }
	void	SetSticky(Int4 min_sticky,char st_md=' ')
		     { assert(min_sticky > 0); MinSticky=min_sticky; StMd=st_md; }
	void	init(){
		aa_per_io=30; aa_per_do=200; exp_ie=1 ;exp_de=1; prior_wt=1; dms_mode='f'; 
		bild_cut=-2.0; bild_add=-2.0; Temp=300; afrq=0.60; rfrq=0.40;
		pn=1000; timeR=time(NULL); gmb=0; MinSticky=2; similarity=0;  StMd=' ';
		routine='x'; Cycle=0;
	}
	// routines:
	cma_typ DoStickySample(FILE *rfp, cma_typ &cma);
	cma_typ DoSingleSample(FILE *rfp, cma_typ &cma);
	BooLean	DoMvColumns(FILE *rfp, cma_typ &cma,Int4 printCycle=0);
	void    DoWrite(cma_typ cma);
	void    DoWriteX(cma_typ cma) { routine='x'; if(write) this->DoWrite(cma); }
	Int4	DoAddRmColumns(FILE *rfp, cma_typ &cma,Int4 MaxCycles=4);
	Int4	DoAddColumns(FILE *rfp, cma_typ &cma);
	Int4	DoPurgeSampling(FILE *rfp, cma_typ &cma);
	gmb_typ	*RtnGambit(){ return gmb; } 
	void	DoNotWrite(){ write=FALSE; }
	void	DoWrite(){ write=TRUE; }
private:
	BooLean write;	// write out intermediate files.
	Int4    similarity,aa_per_io,aa_per_do,exp_ie,exp_de;
	double	prior_wt,pn;
	char	dms_mode,str[100],*name,routine,StMd;
	double	bild_cut,bild_add,Temp,afrq,rfrq;
	Int4	timeR,iter,stage,Cycle,MinSticky;
	gmb_typ	*gmb;
	//static const char      StMd='R',Strategy[20]=" FRTWFRT           ";  // 'D' ruins the alignment!!!
                          //   F    R    T    W    F    R    T    D
        // static const double sfrq[20]={ 0.0,0.49,0.49,0.25,0.33,0.33,0.33,0.49,0.0,0.0};

};

#endif	

