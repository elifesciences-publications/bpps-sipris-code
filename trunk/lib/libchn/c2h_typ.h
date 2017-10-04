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

#if !defined (_C2H_TYP_)
#define _C2H_TYP_

#include "bpcp_typ.h"
// #include "hat_typ.h"
#include "hpt_typ.h"
#include "fdt_typ.h"
#include "set_typ.h"
#include "dheap.h"
#include <time.h>

class c2h_typ {	// compare two hierarchies type.
public:
	c2h_typ( ) { assert(!"Illegal constructor"); }
	c2h_typ(Int4 argc, char *argv[]){ Init(0,argc,argv); }
	c2h_typ(char code, Int4 argc, char *argv[]){ Init(code,argc,argv); }
	Int4	Run( );
	~c2h_typ( ) { Free(); }

	Int4	*EqualArray[2];
	Int4	*NodeMap;	// mapping of nodes from the rearranged to the original hierarchy. 
	Int4	*BackMap;	// mapping of nodes from the original to the rearranged hierarchy. 
	void	NoTree(){ print_tree=FALSE; }

	// return procedures...
	wdg_typ	RtnTreeI(){  wdg_typ T=TreeI; TreeI=0; return T; }
	set_typ	*RtnSetsI(Int4 &N){ set_typ *S=isets; N=NumISets; NEW(isets,N+3,set_typ); return S; }
	set_typ	*RtnMiscSetsI(Int4 &N){
		   set_typ *S=misc_isets; N=NumISets; NEW(misc_isets,N+3,set_typ); return S; 
		}
	// wdg_typ	RtnTreeO(){  wdg_typ T=TreeO; TreeO=0; return T; }
	Int4	*RtnMatchesI(){  Int4 *A=EqualArray[1]; EqualArray[1]=0; return A; }
	hpt_typ	*CopyOfHptI(){  return ihpt->Copy(); }

private:
#if 1	// for finding a consensus hierarchy.
	wdg_typ	TreeI,TreeO;
	Int4	LenEqOI[4];
	char	**EqualOI; // y = o; x = i;
	BooLean	print_tree;
	BooLean	print_stdout;
	// FILE	*sofp,*ofp,*afp,*efp;
#endif
	set_typ IdenticalSets(Int4 Num[2], set_typ *sets[2],char **NameRow[2])
		{ set_typ set; IdenticalSets(Num, sets,NameRow,0,set); return set; }
	Int4	IdenticalSets(Int4 Num[2], set_typ *sets[2],char **NameRow[2],hpt_typ *Hpt,set_typ &set);

	set_typ IdenticalFuzzySets(Int4 Num[2], set_typ *sets[2],set_typ *Psets[2],char **NameRow[2])
		{ set_typ set; IdenticalFuzzySets(Num,sets,Psets,NameRow,0,set); return set; }
	Int4	IdenticalFuzzySets(Int4 Num[2], set_typ *sets[2],set_typ *Psets[2],char **NameRow[2],hpt_typ *hpt,
			set_typ &RtnSet);
	Int4	*MapO2I( );
	Int4    *GetTransitiveParentO(FILE *fp,Int4 Root,Int4 *Parent);
	Int4    *SortByOverlapSetO(FILE *fp,Int4 RootO,Int4 *Po, Int4 *Pi);
	void    ReInit(Int4 *Row);

	set_typ SplitSets(Int4 Num[2], set_typ *Sets[2],char **NameRow[2],hpt_typ *hpt[2],set_typ &RtnISet);
	set_typ NewSets(Int4 Num[2], set_typ *sets[2],char **NameRow[2],hpt_typ *Hpt[2]);
	BooLean CompareCMSAs(FILE *fp);
	void	Init(char code, Int4 argc, char *argv[]);
	void	PrintError(char code,char *program_name);
	// void	GetArg(Int4 argc, char *argv[]);
	void	Free();
	char	*infile,*onfile;
	Int4	CheckInput(Int4 N_cma, Int4 N_rows, cma_typ *CMSA,char **Name_row);
	Int4    Intersection(ss_type P1, ss_type P2);
	// fdt_typ *FDT;	// FD-table type.
	//***************** input versus read/create objects. *****************
	set_typ *GetJackknifeSets(hpt_typ *hpt, cma_typ *cma,Int4 num_cma);
	a_type	AB;
	cma_typ	MainCMA,*IN_CMA,*ON_CMA;
	set_typ	*osets,*isets;
	Int4	NumISets,NumOSets;
	hpt_typ	*ihpt,*ohpt;
	hpt_typ	*iHpt,*oHpt,*HPT[2];
	Int4	NumIn,NumOn;
	char	tmpstr[200],tmp_str[200];
	char	*filename1,*filename2;
	//***************** input versus read/create objects. *****************
	Int4	AddChildrenToSet(Int4 num, char **row, set_typ Set,char **NameSet){
		  for(Int4 c=num-1; c > 1; c--){
                    if(MemberSet(c,Set)){	// previously added children will be skipped.
                      for(Int4 r=1; r < num; r++){ // columns...
                        if(row[r][c] == '+'){
			   fprintf(stdout,"   child: %d.%s\n",r,NameSet[r]);
			   AddSet(r,Set);  // row 'r' column 'c'.
			}
                      }
                    }
                  }
		}
	void	PrintBuffer(FILE *fp,Int4 card, const char *str){ PrintBuffer(fp,card, str,0,0); }
	void	PrintBuffer(FILE *fp,Int4 card, const char *str,Int4 kard, const char *str2){
			assert(buffer[50000] == 0);	// overflow check.
		        if(str2) fprintf(fp,"----------- %d %s %d %s. -----------\n",card,str,kard,str2);
			else fprintf(fp,"----------- %d %s. -----------\n",card,str);
			if(buffer[0] != 0) fprintf(fp,"%s",buffer);
		}
	char	buffer[50005],*bptr;
	Int4	CountSets(Int4 Num[2], set_typ *Sets[2], Int4 &num0, Int4 &num1)
		{
			Int4	s,n; // num0,num1;
			for(n=0,s=1; s < Num[0]; s++) if(Sets[0][s]) n++; num0=n;
			for(n=0,s=1; s < Num[1]; s++) if(Sets[1][s]) n++; num1=n;
			n=num0 - num1;
			// fprintf(fp,"%d vs %d subtrees (%d)\n",num0,num1,n);
			return (n);
		}
	Int4	Nempty,Mempty;
	set_typ	*oSets,*iSets;	// Sets for intermediate (parent) node leftover seqs.
	set_typ	*misc_osets,*misc_isets;
	set_typ	leaf_osets,leaf_isets,tmp_isets,tmp_osets,inter_osets,inter_isets;
	set_typ	match_osets,match_isets;
	set_typ	deleted_osets,deleted_isets;
	set_typ	source_osets,source_isets;
	set_typ	split_osets,split_isets;
	set_typ	fuzzy_osets,fuzzy_isets;
	set_typ	Fuzzy_osets,Fuzzy_isets;	// fuzzy with oSets or iSets.
	set_typ	Fuzzy1_osets,Fuzzy1_isets;	// fuzzy with oSets or iSets.
	set_typ	Fuzzy2_osets,Fuzzy2_isets;
	set_typ	m2m_osets,m2m_isets;
	set_typ	expanded_osets,expanded_isets;
	set_typ	condensed_osets,condensed_isets;
	set_typ	new_osets,new_isets;
	Int4	NumORows,NumIRows;
	char	**orow,**irow,**NameORow,**NameIRow;
	time_t  time1;
};

#endif

