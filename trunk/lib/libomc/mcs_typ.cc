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
#include "blosum62.h"
#include "rst_typ.h"
#include "swt_typ.h"

mcs_typ::mcs_typ(cma_typ in_cma, cma_typ in_mcma, hsw_typ hsw, Int4 NumInSets, set_typ *InSet, 
			hpt_typ *in_hpt, cma_typ *in_sma, Int4 argc, char *argv[])
{
	num_passed_in_sma=in_hpt->NumSets()-1; passed_in_sma=in_sma; passed_in_hpt=in_hpt; 
	passed_in_cma=in_cma; passed_in_hsw=hsw; passed_in_mcma=in_mcma; // own_in_sma=FALSE;
	num_passed_in_sets=NumInSets; passed_in_sets=InSet;
	// Make sure that sets are disjoint.
	for(Int4 i=1; i <= NumInSets; i++){
	   for(Int4 j=i+1; j <= NumInSets; j++){
		Int4 X=CardInterSet(InSet[i],InSet[j]);
		if(X != 0){
		   fprintf(stderr,"Set %d = %d\n",i,CardSet(InSet[i]));
		   fprintf(stderr,"Set %d = %d\n",j,CardSet(InSet[j]));
		   set_typ SetXXX=CopySet(InSet[i]); IntersectSet1(InSet[i],InSet[j],SetXXX); PutSet(stderr,SetXXX); 
		   fprintf(stderr,"Intersection of input sets %d and %d = %d seqs\n",i,j,X);
		   print_error("mcs_typ input error: sequence sets not disjoint");
		}
	   }
	} Init(argc, argv); 
}

mcs_typ::mcs_typ(cma_typ in_cma, cma_typ in_mcma, hsw_typ hsw, Int4 NumInSets, set_typ *InSet, 
			Int4 argc, char *argv[])
{
	num_passed_in_sma=0; passed_in_sma=0; passed_in_hpt=0; 
	passed_in_cma=in_cma; passed_in_hsw=hsw; passed_in_mcma=in_mcma; // own_in_sma=FALSE;
	num_passed_in_sets=NumInSets; passed_in_sets=InSet;
	// Make sure that sets are disjoint.
	for(Int4 i=1; i <= NumInSets; i++){
	   for(Int4 j=i+1; j <= NumInSets; j++){
		Int4 X=CardInterSet(InSet[i],InSet[j]);
		if(X != 0) print_error("mcs_typ input error: sequence sets not disjoint");
	   }
	} Init(argc, argv); 
}

mcs_typ::mcs_typ(cma_typ in_cma, cma_typ in_mcma, hsw_typ hsw, hpt_typ *in_hpt, cma_typ *in_sma,
					Int4 argc, char *argv[])
{
	num_passed_in_sma=in_hpt->NumSets()-1; passed_in_sma=in_sma; // own_in_sma=FALSE;
	passed_in_hpt=in_hpt; passed_in_cma=in_cma; passed_in_hsw=hsw; passed_in_mcma=in_mcma;
	num_passed_in_sets=0; passed_in_sets=0; Init(argc, argv); 
}

mcs_typ::mcs_typ(cma_typ in_cma, cma_typ in_mcma, hsw_typ hsw, Int4 argc, char *argv[])
{
	num_passed_in_sma=0; passed_in_sma=0; passed_in_hpt=0; // own_in_sma=FALSE;
	passed_in_cma=in_cma; passed_in_hsw=hsw; passed_in_mcma=in_mcma;
	num_passed_in_sets=0; passed_in_sets=0; Init(argc, argv); 
}

BooLean	mcs_typ::MvUpContribution(FILE *fp,Int4 node,Int4 parent,double &Rtn)
// put the contributions to the Map for each set and intermediate node.
{
     Int4	*Parent;
     char	x,y,**HP=Hpt->RtnHyperPartition();
     double	d,D,lprFG,lprBG;

     assert(IsTreeHpt);	// make sure this is a tree...
     assert(node > 1 && node < Hpt->NumSets()); assert(!IsFailedSet[node]);
     assert(parent > 0 && parent < node); assert(!IsFailedBPPS[parent]);
     if(HP[node][parent] != '+') return FALSE;
     if(fp) fprintf(fp,"%d.%s(%d):",node,Hpt->ElmntSetName(node),CardSet(GrpSet[node]));
     assert(Hpt->TypeOfSet(parent) == '?'); // this should correspond to an internal node.

     lprFG=CalcTotalLPR(0,FALSE);
     assert(Hpt->IsScrambledTree(Parent));
     set_typ subtree=MakeSet(Hpt->NumBPPS()+2); ReSetSubTree(subtree,node);
     Int4 p = Parent[node], gp = Parent[p]; assert(p == parent);
     if(this->MoveUp(gp,p,node,subtree)) ReSetSubTree(subtree,node);
     lprBG=CalcTotalLPR(0,FALSE);
     if(this->MoveDown(gp,p,node,subtree)) ReSetSubTree(subtree,node);
     NilSet(subtree); free(Parent);

     if(fp) fprintf(fp," (%d: '%c') %.2f",parent,HP[node][parent],lprFG-lprBG);
     Rtn=lprFG-lprBG;
     if(fp) fprintf(fp,"\n");
     return TRUE;
}

BooLean	mcs_typ::QuickMvDownContrib(FILE *fp, Int4 node,Int4 sibling,double &Rtn)
// put the contributions to the Map if node were to be attached to one of its siblings.
{
     Int4	*Parent,p,gp;
     char	x,y,**HP=Hpt->RtnHyperPartition();
     double	d,D,lprN,lprM;

     assert(IsTreeHpt);	// make sure this is a tree...
     assert(node > 1 && node < Hpt->NumSets()); assert(!IsFailedSet[node]);
     assert(sibling > 1 && sibling < Hpt->NumSets()); assert(!IsFailedSet[sibling]);
     assert(Hpt->IsScrambledTree(Parent));
     if(node == sibling){ free(Parent);  return FALSE; }
     if(Parent[node] != Parent[sibling]){ free(Parent); return FALSE; }
     assert(Parent[node] == Parent[sibling]);
     gp = Parent[node]; p = sibling; 
// if(fp) fprintf(fp,"Moving node %d from %d to %d\n",node,gp,sibling);
     assert(HP[node][sibling] == '-'); assert(HP[sibling][node] == '-');
     assert(Hpt->TypeOfSet(gp) == '?'); // this should correspond to an internal node.

     if(fp) fprintf(fp,"%d.%s(%d):",node,Hpt->ElmntSetName(node),CardSet(GrpSet[node]));
     lprN=CalcTotalLPR(0,FALSE);
     set_typ subtree=MakeSet(Hpt->NumBPPS()+2); ReSetSubTree(subtree,node);
// if(fp) fprintf(fp,"\n=============== Moving down (%.3f). ===============\n",lprN);
     if(this->MoveDown(gp,p,node,subtree)) ReSetSubTree(subtree,node);
     this->SampleColumns(); lprM=CalcTotalLPR(0,FALSE);
// Hpt->Put(stderr);
     if(fp) fprintf(fp," (%d: '%c') %.2f",sibling,HP[node][sibling],lprN-lprM);
     Rtn=lprN-lprM;
     if(fp) { if(Rtn < 0.0) fprintf(fp," <===\n"); else fprintf(fp,"\n"); }
     // if(Rtn < 0.0) Hpt->Put(stderr);
     NilSet(subtree); free(Parent);
     return TRUE;
}

BooLean	mcs_typ::MvDownContribution(FILE *fp, Int4 node,Int4 sibling,double &Rtn)
// put the contributions to the Map if node were to be attached to one of its siblings.
{
     Int4	*Parent,p,gp;
     char	x,y,**HP=Hpt->RtnHyperPartition();
     double	d,D,lprN,lprM;

     assert(IsTreeHpt);	// make sure this is a tree...
     assert(node > 1 && node < Hpt->NumSets()); assert(!IsFailedSet[node]);
     assert(sibling > 1 && sibling < Hpt->NumSets()); assert(!IsFailedSet[sibling]);
     assert(Hpt->IsScrambledTree(Parent));
     if(node == sibling){ free(Parent);  return FALSE; }
     if(Parent[node] != Parent[sibling]){ free(Parent); return FALSE; }
     assert(Parent[node] == Parent[sibling]);
     gp = Parent[node]; p = sibling; 
// if(fp) fprintf(fp,"Moving node %d from %d to %d\n",node,gp,sibling);
     assert(HP[node][sibling] == '-'); assert(HP[sibling][node] == '-');
     assert(Hpt->TypeOfSet(gp) == '?'); // this should correspond to an internal node.

     if(fp) fprintf(fp,"%d.%s(%d):",node,Hpt->ElmntSetName(node),CardSet(GrpSet[node]));
// this->SampleColumns();
     lprN=CalcTotalLPR(0,FALSE);
     set_typ subtree=MakeSet(Hpt->NumBPPS()+2); ReSetSubTree(subtree,node);
// if(fp) fprintf(fp,"\n=============== Moving down (%.3f). ===============\n",lprN);
     if(this->MoveDown(gp,p,node,subtree)) ReSetSubTree(subtree,node);
// this->SampleColumns();
// Hpt->Put(stderr);
     lprM=CalcTotalLPR(0,FALSE);
// if(fp) fprintf(fp,"\n=============== Moving up (%.3f). ===============\n",lprM);
     if(fp) fprintf(fp," (%d: '%c') %.2f",sibling,HP[node][sibling],lprN-lprM);
     Rtn=lprN-lprM;
     if(fp) { if(Rtn < 0.0) fprintf(fp," <===\n"); else fprintf(fp,"\n"); }
     // if(Rtn < 0.0) Hpt->Put(stderr);
     if(this->MoveUp(gp,p,node,subtree)) ReSetSubTree(subtree,node);
     NilSet(subtree); free(Parent);
     return TRUE;
}

// #include "blosum62.h"

#if 0	//**********************************************************************
Syntax: <gi>(begin..end) <category>[;<category>] [(begin..end) <category>[;<category>]]

Example input:
362794	(25..135) Set35=-1.34; Set3=3.4
9117264	(5..234) Set35=2.34 (8..230) Set15=2.01

Space characters may be present but are unnecessary and can be easily parsed out.

> sed 's/[\t ]//g' example | sed 's/^.*[);];Set3=/Set3=/'
362794(25..135)Set35=-1.34;Set3=3.4
9117264(5..234)Set35=2.34(8..230)Set15=2.01

We can unambiguously parse out any fields we like with a command like the following:

> sed 's/[\t ]//g' example | grep '[);]Set3=' | sed 's/(.*[);]Set3=/ /'
362794 3.4

#endif	//**********************************************************************

set_typ *mcs_typ::RtnSeqSets() { if(DidRestoreBest) return BestSet; else return GrpSet; }

set_typ	*mcs_typ::RtnSubTreeSeqSet( )
// return full sets for non-root internal nodes...
{
	Int4	r,c,*P;
	// assert(Hpt->IsScrambledTree(P));
	assert(Hpt->IsTree(P));
	set_typ *set=this->RtnSeqSets();
	set_typ *rtn_sets=0; NEW(rtn_sets,Hpt->NumSets() +2, set_typ);
	for(r=2; r < Hpt->NumSets(); r++){
	   if(Hpt->TypeOfSet(r) == '?') rtn_sets[r]=CopySet(set[r]);
	}
	for(c=2; c <= Hpt->NumBPPS(); c++){
	   if(Hpt->TypeOfSet(c) != '?') continue;
	   for(r=c+1; r < Hpt->NumSets(); r++){
              if(Hpt->Cell(r,c) == '+') UnionSet(rtn_sets[c],set[r]);
	   }
	} free(P);
	return rtn_sets;
}

BooLean	mcs_typ::PutMapContributions(FILE *fp,lpr_typ *xlpr)
// put the contributions to the Map for each set and intermediate node.
{
     assert(IsTreeHpt);	// make sure this is a tree...
     Int4	n,g,sq,num_sq = NumSeqsCMSA(TrueMainCMA);
     char	x,**HP=Hpt->RtnHyperPartition();
     BooLean	rtn=TRUE;

#if 0	// can replace with this...
     for(g=1; g <= Hpt->NumBPPS(); g++) PutMapContributions(fp,g); return TRUE;
#endif
#if 0	// Probably don't need this...
     for(g=2; g < Hpt->NumSets(); g++){	  // check to see if hpt is in tree format...
	if(g > Hpt->NumBPPS()){ rtn=FALSE; break; }
	if(IsFailedSet[g]) continue;
        for(n=1; n < g; n++){ 
	   if(IsFailedBPPS[n]) continue;  // this very rarely be true;
	   if(Hpt->TypeOfSet(n) != '?'){ rtn=FALSE; break; } // should be an internal node.
	} if(rtn==FALSE) break;
     }
     if(!rtn){
	fprintf(fp,"PutMapContributions() failed; fd-table not in tree format\n");
	return rtn;
     }
#endif
     double	lpr,lpr0,*subLLR;
     NEW(subLLR,Hpt->NumBPPS() +3,double); CalcTotalLPR(0,FALSE); 
     for(g=0; g <= Hpt->NumBPPS(); g++){ subLLR[g]=Map[g]; }
     fprintf(fp,"MAP Contributions= <node>.<set>(<size>): ...(<subnode>) <subset> [<fullset>}\n");
     fprintf(fp,"%d.%s(%d): (*) %.2f\n",1,Hpt->ElmntSetName(1),CardSet(GrpSet[1]),subLLR[1]);
     set_typ  *FullSet=this->RtnSubTreeSeqSet( );
     for(g=2; g < Hpt->NumSets(); g++){	
	assert(g <= Hpt->NumBPPS());
	if(IsFailedSet[g]) continue;
	// Int4 card=CardSet(GrpSet[g]);
	if(FullSet[g]) fprintf(fp,"%d.%s(%d|%d):",g,Hpt->ElmntSetName(g),CardSet(GrpSet[g]),CardSet(FullSet[g]));
	else fprintf(fp,"%d.%s(%d):",g,Hpt->ElmntSetName(g),CardSet(GrpSet[g]));
        for(n=1; n < g; n++){ 
	   if(IsFailedBPPS[n]) continue;  // this should never be true;
	   if(HP[g][n] == '+'){
		assert(Hpt->TypeOfSet(n) == '?'); // this should correspond to an internal node.
		x='-';  // Put in background partition.
		// x='o';  // Omit from analysis.

		for(sq=1; sq <= num_sq; sq++){
		   if(MemberSet(sq,GrpSet[g])){ che[n]->SetPartition(x,sq); } 
		} lpr0=che[n]->CalcLLR( ); lpr=subLLR[n];
		fprintf(fp," (%d) %.1f",n,lpr-lpr0);
		for(sq=1; sq <= num_sq; sq++){
		   if(MemberSet(sq,GrpSet[g])){ che[n]->SetPartition('+',sq); } 
		} 

		if(FullSet[g] == 0) continue;
		for(sq=1; sq <= num_sq; sq++){
		   if(MemberSet(sq,FullSet[g])){ che[n]->SetPartition(x,sq); } 
		} lpr0=che[n]->CalcLLR( ); lpr=subLLR[n];
		fprintf(fp," [%.1f]",lpr-lpr0);
		for(sq=1; sq <= num_sq; sq++){
		   if(MemberSet(sq,FullSet[g])){ che[n]->SetPartition('+',sq); } 
		} 
	   }
	} fprintf(fp," (*) %.1f\n",subLLR[n]);
     } 
     for(g=1; g < Hpt->NumSets(); g++) if(FullSet[g]) NilSet(FullSet[g]);
     fprintf(fp,"\n"); free(subLLR); free(FullSet);
     return rtn;
}

BooLean	mcs_typ::RtnContribLLR(Int4 row, Int4 col, double &subLpr)
// n=column; g= row; n < g.
{
	Int4	sq,num_sq=NumSeqsCMSA(TrueMainCMA);
	char	x,y,**HP=Hpt->RtnHyperPartition();
	double	lpr,lpr0;
	if(row == col) return FALSE;  // allow row same as column?
	if(row < 1 || row >= Hpt->NumSets()) return FALSE;
	if(col < 1 || col > Hpt->NumBPPS()) return FALSE;
	if(Hpt->TypeOfSet(col) != '?') return FALSE; // column should correspond to an internal node.
	if(IsFailedBPPS[col]) return FALSE;  // this should never be true;
#if 0
	if(HP[row][col] != '+') return FALSE;
	x='-'; // x='o';  // Put in background partition ('-') or omit from analysis ('o').
	lpr=che[col]->CalcLLR( ); 
	set_typ FullSet=0,TheSet=0;
	if(Hpt->TypeOfSet(row) == '?') FullSet=this->RtnSubTreeSeqSet(row); 
	if(FullSet) TheSet=FullSet; else TheSet=GrpSet[row];
	for(sq=1; sq <= num_sq; sq++){ if(MemberSet(sq,TheSet)) che[col]->SetPartition(x,sq); }
	subLpr=lpr-che[col]->CalcLLR( ); 
	for(sq=1; sq <= num_sq; sq++){ if(MemberSet(sq,TheSet)) che[col]->SetPartition('+',sq); }
	if(FullSet) NilSet(FullSet);
#else
	// x='-' or x='o';  // Put in background partition ('-') or omit from analysis ('o').
	// if(HP[row][col] == '+'){ x='-'; y='+'; }  else return FALSE;
	if(HP[row][col] == '+'){ x='-'; y='+'; }  
	else if(HP[row][col] == '-'){ x='+'; y='-'; } else return FALSE;
	lpr=che[col]->CalcLLR( ); 
	set_typ FullSet=0,TheSet=0;
	if(Hpt->TypeOfSet(row) == '?') FullSet=this->RtnSubTreeSeqSet(row); 
	if(FullSet) TheSet=FullSet; else TheSet=GrpSet[row];
	for(sq=1; sq <= num_sq; sq++){ if(MemberSet(sq,TheSet)) che[col]->SetPartition(x,sq); }
	subLpr=lpr-che[col]->CalcLLR( ); 
	for(sq=1; sq <= num_sq; sq++){ if(MemberSet(sq,TheSet)) che[col]->SetPartition(y,sq); }
	if(FullSet) NilSet(FullSet);
#endif
	return TRUE;
}

void	mcs_typ::CopySubTreeSeqSet(Int4 node, set_typ rtnSet)
// return full set for node's subtree == foreground.
{
	Int4	r,c,*P;
	assert(Hpt->IsTree(P)); free(P);
	assert(node > 0 && node <= Hpt->NumSets());
	set_typ *set=this->RtnSeqSets();
	assert(SetN(set[1]) == SetN(rtnSet));
	CopySet(rtnSet,set[node]);
	if(Hpt->TypeOfSet(node) != '?') return;
	for(r=1; r <= Hpt->NumSets(); r++){
	      if(r==node) continue;
              if(Hpt->Cell(r,node) == '+') UnionSet(rtnSet,set[r]);
	} return;
}

void	mcs_typ::CopyBkGrndSeqSet(Int4 node, set_typ rtnSet)
// return full set for node's background.
{
	Int4	r,c,*P;
	assert(Hpt->IsTree(P)); free(P);
	assert(node > 0 && node <= Hpt->NumSets());
	set_typ *set=this->RtnSeqSets();
	assert(SetN(set[1]) == SetN(rtnSet));
	ClearSet(rtnSet);
	for(r=1; r <= Hpt->NumSets(); r++){
              if(Hpt->Cell(r,node) == '-') UnionSet(rtnSet,set[r]);
	} return;
}

set_typ	mcs_typ::RtnSubTreeSeqSet(Int4 col)
// return full sets for non-root internal nodes...
{
	Int4	r,c,*P;
	// assert(Hpt->IsScrambledTree(P));
	assert(Hpt->IsTree(P)); free(P);
	set_typ *set=this->RtnSeqSets();
	set_typ rtn_set=0; 
#if 0
	if(Hpt->TypeOfSet(col) == '?') rtn_set=CopySet(set[col]); else return 0;
	for(r=col+1; r < Hpt->NumSets(); r++){
              if(Hpt->Cell(r,col) == '+') UnionSet(rtn_set,set[r]);
	} return rtn_set;
#else
	rtn_set=CopySet(set[col]);
	for(r=1; r < Hpt->NumSets(); r++){
              if(r != col && Hpt->Cell(r,col) == '+') UnionSet(rtn_set,set[r]);
	} return rtn_set;
#endif
}

void	mcs_typ::PutHyperPartition(FILE *fp)
// finding out why some GoldStd sequences are getting pushed out.
{
     Int4	n,g,p;
     char	c,**HP=Hpt->RtnHyperPartition();
if(Hpt->NumBPPS() <= 25){ 
     fprintf(fp,"\n ");
     for(n=1; n<= Hpt->NumBPPS(); n++){ fprintf(fp,"="); }
     fprintf(fp," HyperPartition: ");
     for(n=1; n<= Hpt->NumBPPS(); n++){ fprintf(fp,"="); } fprintf(fp,"\n");
     fprintf(fp,"      _Category_\nSet: ");
     for(n=1; n<= Hpt->NumBPPS(); n++){ fprintf(fp,"%2d ",n); } fprintf(fp,"\n");
     for(g=1; g<= Hpt->NumSets(); g++){	
        fprintf(fp,"%3d: ",g);
	Int4 NumZero=0;
        for(p=0,n=1; n<= Hpt->NumBPPS(); n++){
	   c = HP[g][n];
	   if(c != 0 && c != 'o'){ fprintf(fp," %c ",c); if(c=='+') p++; }
	   // else { NumZero++; fprintf(fp," o "); }
	   else { NumZero++; fprintf(fp,"   "); }
	} 
	Int4 card=CardSet(GrpSet[g]);
	char c=' '; if(IsFailedSet[g]) c='x'; else c=' '; 
	if(strcmp("Random",Hpt->ElmntSetName(g)) == 0){
		// fprintf(fp," %s (%d)\n",Hpt->ElmntSetName(g),CardSet(GrpSet[g])-NumRandom);
		fprintf(fp," Rejected (%d)\n",card-NumRandom);
	} else {
	  if(p > 2) fprintf(fp," _"); else fprintf(fp," ");
	  if(Hpt->TypeOfSet(g) == '?'){
		fprintf(fp,"%d.%s {%d}* %c\n",g,Hpt->ElmntSetName(g),card,c);
	  } else fprintf(fp,"%d.%s (%d) %c\n",g,Hpt->ElmntSetName(g),card,c);
	}
	// else fprintf(fp," %s (%d)\n",Hpt->ElmntSetName(g),CardSet(GrpSet[g]));
	if(NumZero == Hpt->NumBPPS()) print_error("Fatal: seq grp absent from all partitions");
     }
     fprintf(fp," =============================================\n"); 
 } else {
     fprintf(fp,"\nHyperPartition:\n ");
     for(n=1; n<= Hpt->NumBPPS(); n++){
	if(n % 10 == 0) fprintf(fp,"|"); else if(n % 5 == 0) fprintf(fp,"+"); else fprintf(fp,"-"); 
     } fprintf(fp,"\n");
     for(g=1; g<= Hpt->NumSets(); g++){	
        fprintf(fp," ");
	Int4 NumZero=0;
        for(p=0,n=1; n<= Hpt->NumBPPS(); n++){
	   c = HP[g][n];
	   if(c != 0 && c != 'o'){ fprintf(fp,"%c",c); if(c=='+') p++; }
	   // else { NumZero++; fprintf(fp,"o"); }
	   else { NumZero++; fprintf(fp," "); }
	} 
	Int4 card=CardSet(GrpSet[g]);
	char c=' '; if(IsFailedSet[g]) c='x'; else c=' '; 
	if(strcmp("Random",Hpt->ElmntSetName(g)) == 0){
		fprintf(fp," %d.Rejected (%d)\n",g,card-NumRandom);
	} else {
	   if(p > 2) fprintf(fp," _"); else fprintf(fp," ");
	   if(Hpt->TypeOfSet(g) == '?'){
		fprintf(fp,"%d.%s {%d}* %c\n",g,Hpt->ElmntSetName(g),card,c);
	   } else fprintf(fp,"%d.%s (%d) %c\n",g,Hpt->ElmntSetName(g),card,c);
	}
	if(NumZero == Hpt->NumBPPS()) print_error("Fatal: seq grp absent from all partitions");
     }
     fprintf(fp," ");
     for(n=1; n<= Hpt->NumBPPS(); n++){
	if(n % 10 == 0) fprintf(fp,"|"); else if(n % 5 == 0) fprintf(fp,"+"); else fprintf(fp,"-"); 
     } fprintf(fp,"\n\n Column LPRs:\n");
 }
     double d=CalcTotalLPR();
     Int4 numFailedCols=0;
     for(n=1; n<= Hpt->NumBPPS(); n++){
#if 0	// DEBUG...
		double *d=che[n]->SubMap( );
		if(d[0] < 0.0){ 
			che[n]->PutSubLPRs(stderr);
			fprintf(stderr,"WARNING: the set size computed for rho is too small.\n");
			fprintf(stderr,"some pattern residue set(s) is/are larger than expected.\n");
			fprintf(stderr,"see GetRhoCategoricalPriors( ) in cmc_init.cc file.\n");
		} // need to see why this is happening...
		// only seems to happen when reading a cmc checkpoint file.
#endif
		fprintf(fp,"%2d: ",n); 
		double lpr=che[n]->PutInfoShort(fp); 
		if(lpr <= 0.0 || che[n]->NumColumns( ) < 1) numFailedCols++;
     } fprintf(fp," ====== Total LPR = %.4f (%.1f K) (%d/%d failed) ======\n\n",
			d,temperature,numFailedCols,Hpt->NumBPPS());
}

void	mcs_typ::Free()
{
	Int4	g,i,n;
	// Note: ElmntSetCMA[i]'s freed up above...
	if(cfp) fclose(cfp); if(ifp) fclose(ifp);
	for(g=1; g < MAX_NUM_ELMENTARY_SETS; g++){
		if(SeedCMA[g]) TotalNilCMSA(SeedCMA[g]);
		if(sst_str[g]) free(sst_str[g]);
		// if(passed_in_sets==0 && InitSet[g]) NilSet(InitSet[g]);
		if(InitSet[g]) NilSet(InitSet[g]);
		if(WorstToBest[g]) free(WorstToBest[g]);
	}
	for(g=1; g <= Hpt->NumSets(); g++){
		// free(Hpt->ElmntSetName(g));
		if(GrpSet[g]) NilSet(GrpSet[g]);
		if(BestSet[g]) NilSet(BestSet[g]);
	}
	for(n=1; n <= Hpt->NumBPPS(); n++){
		assert(che && che[n]);
		Int4 length = che[n]->BPPS()->LenPattern( );
		if(che && che[n]) delete che[n]; 
		if(chn && chn[n]) delete chn[n]; 
		if(QryCMAs[n]) TotalNilCMSA(QryCMAs[n]);
		for(i=1; i <= length; i++){
			if(sst && sst[n] && sst[n][i]) free(sst[n][i]);
		} if(sst && sst[n]) free(sst[n]);
		if(RelateFGs && RelateFGs[n]) free(RelateFGs[n]);
		if(RelateBGs && RelateBGs[n]) free(RelateBGs[n]);
		if(SetFG && SetFG[n]) NilSet(SetFG[n]);
		if(SetBG && SetBG[n]) NilSet(SetBG[n]);
		if(best_sst && best_sst[n]) free(best_sst[n]);
		if(BestCsq[n]) NilSeq(BestCsq[n]);
// fprintf(stderr,"cmc->Free %d_4\n",n);
	} if(QryCMAs) free(QryCMAs);
	if(RelateFGs) free(RelateFGs);
	if(RelateBGs) free(RelateBGs);
	if(SetFG) free(SetFG);
	if(SetBG) free(SetBG);

	if(Labeled) NilSet(Labeled);
	if(RandomSet) NilSet(RandomSet);
	if(chn && chn[0]) delete chn[0]; 
	if(che) free(che); 
	if(sst) free(sst); 
	if(chn) free(chn);
   	if(MaxNumCol) free(MaxNumCol);
   	if(program_name) free(program_name);
   	if(infile) free(infile);
	if(SFBG) free(SFBG);
	if(TrueMainCMA && !passed_in_cma) TotalNilCMSA(TrueMainCMA);
	if(MainCMA && !passed_in_mcma) TotalNilCMSA(MainCMA);
	if(dummyCMA) TotalNilCMSA(dummyCMA);
#if 0	// == DisplayCMA[n]
	if(passed_in_sma && own_in_sma){
	   for(n=1; n <= num_passed_in_sma; n++){ TotalNilCMSA(passed_in_sma[n]); }
	}
#endif
	if(DisplayCMA){ 
		for(n=1; n <= NumDisplayCMA; n++){
		    if(DisplayCMA[n]) TotalNilCMSA(DisplayCMA[n]);
		} free(DisplayCMA);
	} if(Hpt) delete Hpt;
	if(BestHpt) delete BestHpt;
	if(outfp && outfp != stdout) fclose(outfp);
	if(cfp) fclose(cfp);
	if(ifp) fclose(ifp);
}

cma_typ *mcs_typ::RtnLeafNodeCMSA(Int4 &Number)
{
    FILE *fp=tmpfile();
    Int4 N,g;
    for(g=1; g < Hpt->NumSets(); g++){	
	if(Hpt->TypeOfSet(g) != '?'){
	   char *name=Hpt->ElmntSetName(g);
	   ReNameCMSA(name,MainCMA);		// don't need to worry about Main Set name.   
	   PutInSetCMSA(fp,GrpSet[g],MainCMA);	// keep labels as inputted
	}
    } rewind(fp);
    cma_typ *CMA=MultiReadCMSA(fp,&N,AB); fclose(fp); Number=N;
    return CMA;
}

double	mcs_typ::RevertToBest()
{
	if(!SaveBest || BestHpt == 0) return CalcTotalLPR(0,FALSE);  // then use current state...
	Int4	s,sq,n,i,Length=LenSeq(che[1]->KeyE()),num_sq = NumSeqsCMSA(MainCMA);
	if(!Hpt->TheSame(BestHpt)){
		delete Hpt; Hpt=BestHpt->Copy(); 
		// for(i=1; i <= Hpt->NumSets(); i++){ HyperPartition[i]=Hpt->RtnHyperPartition(i); }
	} for(i=1; i <= Hpt->NumBPPS(); i++){ che[i]->SetConsensus(BestCsq[i]); }
	char	**HP=Hpt->RtnHyperPartition();
	for(s=1; s<= Hpt->NumSets(); s++){	
		IsFailedSet[s]=IsFailedBestSet[s];
		CopySet(GrpSet[s],BestSet[s]);	// copies BestSet[s] to GrpSet[s].
	}
     	for(sq=1; sq <= num_sq; sq++){
	  if(MemberSet(sq,Labeled)) continue; // no need to reset this; never changes.
	  for(n=1; n<= Hpt->NumBPPS(); n++){
	    for(s=1; s<= Hpt->NumSets(); s++){	
	        if(MemberSet(sq,GrpSet[s])){ che[n]->SetPartition(HP[s][n],sq); }
	    }
	  }
	}
	for(n=1; n<= Hpt->NumBPPS(); n++){
	  IsFailedBPPS[n]=IsFailedBestBPPS[n];
	  bpps_typ *bpps=che[n]->BPPS();
          for(i=1; i <= Length; i++){
		// che[n]->RmColButKeepLPR(i); // WARNING: won't work if # columns <= mincol!!!
		bpps->RmColumn(i); // this will work...
		if(best_sst[n][i]) che[n]->AddColumn(i,best_sst[n][i]);
	  }
	} ReSetRelations();
	return CalcTotalLPR(0,FALSE);	// don't want to store it again.
}

double	mcs_typ::RestoreBest()
{
	if(DidRestoreBest) return CalcTotalLPR(0,FALSE); // don't want to store it again.
#if 0
	if(!SaveBest) return CalcTotalLPR(0,FALSE);	 // just use current configuration...
#else
	if(!SaveBest || BestHpt == 0){ 	// never saved; use current...
		SaveBest=TRUE; StoreBest(); return CalcTotalLPR(0,FALSE);
	}
#endif
#if 0	// DEBUG...
this->PutPttrnVsConsSeq(stderr,2);
this->PutHyperPartition(stderr);
CheckPttrnCsqMatch("Error 0 in mcs_typ->RestoreBest()");
#endif
	// CheckPttrnCsqMatch("debug 1");
	Int4	s,sq,n,i,Length;

	// New: for sampling over Hpt.
	if(!Hpt->TheSame(BestHpt)){
		delete Hpt; Hpt=BestHpt->Copy(); 
		// for(i=1; i <= Hpt->NumSets(); i++){ HyperPartition[i]=Hpt->RtnHyperPartition(i); }
	} for(i=1; i <= Hpt->NumBPPS(); i++){ che[i]->SetConsensus(BestCsq[i]); }

	char	**HP=Hpt->RtnHyperPartition();
	for(s=1; s<= Hpt->NumSets(); s++){	
		IsFailedSet[s]=IsFailedBestSet[s];
		CopySet(GrpSet[s],BestSet[s]);	// copies BestSet[s] to GrpSet[s].
	}
        Int4 num_sq = NumSeqsCMSA(MainCMA);
     	for(sq=1; sq <= num_sq; sq++){
	  if(MemberSet(sq,Labeled)) continue; // no need to reset this; never changes.
	  for(n=1; n<= Hpt->NumBPPS(); n++){
	    for(s=1; s<= Hpt->NumSets(); s++){	
	        if(MemberSet(sq,GrpSet[s])){ che[n]->SetPartition(HP[s][n],sq); }
	    }
	  }
	}
	for(n=1; n<= Hpt->NumBPPS(); n++){
	  Length=LenSeq(che[n]->KeyE());
	  IsFailedBPPS[n]=IsFailedBestBPPS[n];
	  bpps_typ *bpps=che[n]->BPPS();
          for(i=1; i <= Length; i++){
		// che[n]->RmColButKeepLPR(i);	// Warning: Won't remove if # columns < mincol!!
		bpps->RmColumn(i); // this will work...
		if(best_sst[n][i]) assert(che[n]->AddColumn(i,best_sst[n][i]));
	  }
	} DidRestoreBest=TRUE; ReSetRelations();
CheckPttrnCsqMatch("Error in mcs_typ->RestoreBest()");
	return CalcTotalLPR(0,FALSE);	// don't want to store it again.
}

void	mcs_typ::StoreBest()
{
	Int4	s,n,i,Length;
	if(TotalLPR <= 0.0) return;
	if(this->NoFailureMode && this->NumFailedNodes() > 0) return;
CheckPttrnCsqMatch("Error in mcs_typ->StoreBest()");
	SaveBest=TRUE;
	BestLPR=TotalLPR;
	for(s=1; s<= Hpt->NumSets(); s++){	
		IsFailedBestSet[s]=IsFailedSet[s];
		CopySet(BestSet[s],GrpSet[s]);	// copies GrpSet[s] to BestSet[s].
	}
     	for(n=1; n<= Hpt->NumBPPS(); n++){
	    IsFailedBestBPPS[n]=IsFailedBPPS[n];
	    Length=LenSeq(che[n]->KeyE());
	    bpps_typ *bpps=che[n]->BPPS();
            for(i=1; i <= Length; i++) best_sst[n][i]=bpps->RtnSST(i);
	}
	// new: storing csq and hpt...
	if(BestHpt ==0) BestHpt=Hpt->Copy(); 
	else if(!Hpt->TheSame(BestHpt)){ delete BestHpt; BestHpt=Hpt->Copy(); } 
	assert(Hpt->NumSets() == Hpt->NumBPPS() + 1);
	for(s=1; s < Hpt->NumSets(); s++){
	   if(BestCsq[s]!=0) NilSeq(BestCsq[s]); BestCsq[s]=CopySeq(che[s]->KeyE());
	}
}

double	mcs_typ::CalcTotalLPR(FILE *fp,BooLean StoreBestOK)
{
	Int4	n;
	BooLean	okay=TRUE;
     	for(TotalLPR=0.0, n=1; n<= Hpt->NumBPPS(); n++){
		if(che[n]->NumColumns() > MaxNumCol[n]){
			// che[n]->RtnMaxNumColumns();
#if 0
	fprintf(stderr,"che[%d]->NumColumns() = %d; MaxNumCol[%d] = %d\n",
			n,che[n]->NumColumns(), n,MaxNumCol[n]);
#endif
			StoreBestOK = FALSE;
		} 
		SubLPR[n]=che[n]->SubMap( ); Map[n]=SubLPR[n][0]; // CheckValue(Map[n]);
		if(fp){ che[n]->PutInfo(fp,n); }
		TotalLPR+=Map[n];
// if(n > 1) assert(!(Map[n-1] > 0.0 && Map[n] < -99856392.0));	// DEBUG!!!
		if(this->NoFailureMode && Map[n] <= 0.0) okay = FALSE;
	}
#if 1	// MDL: Adjust for number of n-ary trees (very conservative as only depth=5 allowed).
	{
		Int4 d=Hpt->NumSets();
     	        double D = lnfact(2*d) - (lnfact(d) + lnfact(d+1));	// Catalan number.
		// fprintf(stderr,"Catalan number (%d nodes) = %.3g\n",Hpt->NumSets(),D);
		TotalLPR -= D;
	} Map[0]=TotalLPR;
#else
	Map[0]=TotalLPR;
#endif

	CheckValue(TotalLPR);
	if(fp) fprintf(fp,"Total LPR = %.3f (%.1f K)\n",TotalLPR,temperature);
#if 0
	fprintf(stderr,"Total LPR = %.3f\n", TotalLPR);
	fprintf(stderr,"okay = %d\n", okay);
	fprintf(stderr,"StoreBestOK = %d\n", StoreBestOK);
	fprintf(stderr,"BestLPR=%.2f\n",BestLPR);
#endif
	if(okay && StoreBestOK && TotalLPR > 0.0 && TotalLPR > BestLPR){
			StoreBest();
	} // if(StoreBestOK && TotalLPR > 0.0 && TotalLPR > BestLPR) StoreBest();
	return TotalLPR;
}

BooLean	mcs_typ::CheckValue(double x)
{
        if(!isfinite(x)){       // this is not a valid floating point number.
          if(isnan(x)){
		print_error("No signficant hyperpartition found (LPR:NaN).");
	  } else if(isinf(x)){
                  if(isinf(x) == 1) print_error("No signficant hyperpartition found (LPR:+infinity).");
                  else if(isinf(x) == -1){      // negative infinity
                        print_error("No signficant hyperpartition found (LPR: neg. infinity).");
                  } else print_error("No signficant hyperpartition found (LPR: infinite).");
          } else {
	    print_error("No signficant hyperpartition found (LPR: not finite).");
            return FALSE;
	  }
        } return TRUE;
}

void	mcs_typ::MoveSeqs(Int4 from, Int4 target)
// Move all sequences in Group from to group target.
{
	Int4	cnt=0,n,sq,N=NumSeqsCMSA(TrueMainCMA); 
	char    **HP=Hpt->RtnHyperPartition();
	assert(target > 0 && target <= Hpt->NumSets());
	assert(from > 0 && from <= Hpt->NumSets());
	for(sq=1; sq <= N; sq++){
	    if(MemberSet(sq,Labeled)) continue; //  DeleteSet(sq,Labeled);;
	    if(MemberSet(sq,GrpSet[from])){
		DeleteSet(sq,GrpSet[from]); AddSet(sq,GrpSet[target]);
		for(n=1; n<= Hpt->NumBPPS(); n++){
		   che[n]->SetPartition(HP[target][n],sq);
		} cnt++;
	    }
	} if(cnt > 0) DidRestoreBest=FALSE;   // this is now different from the best.
}

