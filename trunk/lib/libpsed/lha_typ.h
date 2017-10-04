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

#if !defined (LHA_TYP)
#define LHA_TYP

#include "gmb_typ.h"
#include "gsm_typ.h"
#include "hpt_typ.h"
#include "set_typ.h"
#include "blosum62.h"
#include "cma_gmb.h"
// #include "omc_typ.h"

#include "hat_typ.h"
#include "residues.h"
#include "swaln.h"
#include "dheap.h"
#include "blosum62.h"

#if 0
class fba_typ {
public:
        fba_typ(Int4 len, UInt4 twc, BooLean *del, a_type ab)
        {
                Int4    i,j,r;
                Leng=len; AB=ab; TotalWtCnts=twc;
                NEW(Deleted,Leng+3,BooLean);
                NEWP(WtCnts,Leng+3,UInt4);
                for(i=1; i <= Leng; i++){
                   NEW(WtCnts[i],nAlpha(AB)+5,UInt4);
                   if(del[i]){
                     Deleted[i]=TRUE;
                     for(r=0; r <= nAlpha(AB); r++){
                        double  d=(double)twc*blosum62freq[r];
                        WtCnts[i][r] = (UInt4) floor(d+0.5);
                     }
                   }
                }
        }
        ~fba_typ( )
        {
                free(Deleted);
                for(Int4 i=1; i <= Leng; i++) free(WtCnts[i]); free(WtCnts);
        }
private:
        Int4    Leng;
        UInt4   **WtCnts,TotalWtCnts;
        BooLean *Deleted;

        a_type  AB;

};
#endif

// Lineage hiearchical alignment type

// hieraln:
set_typ *MkSubSets(set_typ *set, Int4 NumSets, set_typ Merge);
cma_typ MkSubSetCMSA(set_typ *NodeSet,hpt_typ *Hpt,Int4 s,BooLean  AddCsq,cma_typ cma,set_typ **rtnSets);
void    PutHieralnDebug(const char s[],Int4 id,char *prefix,char *suffix,e_type csqE,cma_typ cma);
char	*MkTemplateOperation(ssx_typ *SSX, ssx_typ *ssx, char *Oper[3]);
cma_typ QuickOptimizeCMA(cma_typ cma,Int4 s, char **Oper);
cma_typ PurgeSampleReAlign(char *name,cma_typ TrueMainCMA);
cma_typ *MkTemplateCMAs(hpt_typ *Hpt,cma_typ *NodeCMA,set_typ *ChildSet,e_type *CsqE,
			char **Oper, a_type AB);

// hybrid:
void    PutSinglesCMSA(FILE *fp,cma_typ cma);
cma_typ GrowViaTemplateCMSA(cma_typ Template, cma_typ SubGrp, Int4 target, e_type TargetCsq);
int     run_hybrid(int argc,char *argv[],FILE *mfp=0,FILE *hfp=0,FILE *stfp=0,FILE *sdfp=0);
int	run_hieraln(int argc, char *argv[]);
int	run_hieraln(char sc,int argc, char *argv[]);

#endif

