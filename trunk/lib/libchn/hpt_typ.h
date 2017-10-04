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

#if !defined (_HPT_TYP_)
#define _HPT_TYP_

#include "alphabet.h"
#include "afnio.h"
#include "dheap.h"
#include "wdigraph.h"
#include "histogram.h"
#include "set_typ.h"

#define	HPT_PUBLIC_USAGE "\n\
************************** hyperpartition file syntax: **************************\n\
HyperParTition:\n\
<partitions> 1.<name><type>\n\
<partitions> 2.<name><type>\n\
   :    :      :   :   :\n\
   :    :      :   :   :\n\
<partitions> <int1>.Random=<int2>.   (<int1> = # sets. <int2> = # random sequences to be generated; default=10000)\n\
\n\
Settings:                            (this section is optional)\n\
1.<category_name> <options>\n\
2.<category_name> <options>\n\
   :    :      :   :\n\
   :    :      :   :\n\
<int>.<category_name> <options>         (where <int> corresponds to the total # categories)\n\
\n\
(Note: add comments by placing a '#' character at the beginning of the line.)\n\
(Note: 'End.' can be used to separate multiple hpts within a single text file.)\n\
=========================== non-terminal definitions: ===========================\n\
  <partitions> = a string of '+','-' & 'o' characters representing the partitions to which the set belongs\n\
	'+'  --> the set belongs to the foreground partition\n\
	'-'  --> the set belongs to the background partition\n\
     	'o'  --> the set belongs to neither partition (i.e., set omitted)\n\
                 (Note that there must be one such character for each category)\n\
                 For example, <partitions> = \"+----++++++--oo-o\" (17 categories)\n\
  <name>  =  name of the corresponding seed alignment within the *.sma input file\n\
  <type>  =  '.'  --> set is specific (e.g., <name><type> = 'Ran.') or\n\
             '?'  --> set is generic (e.g., <name><type> = 'MiscRasLike?')\n\
             '!'  --> a specific set for which a *.rtf contrast alignment will be generated.\n\
\n\
  <category_name> = an arbitrary character string that describes the category\n\
  <options> =\n\
     -A<real>:<real> - alpha hyperparameters A0:B0 (default: A0=B0=1.0)\n\
     -col=<int>:<int>  - Specify the min and max number of columns allowed\n\
     -N=<int>      - Maximum number of significant pattern positions to highlight\n\
                   - (This sets the contrast for the output alignment.)\n\
     -Ri=<real>    - Set prior probability that a row (seq) is in the foreground (default: 0.5)\n\
\n\
\n"

#define MAX_NUM_ELMENTARY_SETS 503 
#define MAX_NUM_SUBGROUPS 503 

class hpt_typ {         // comprehensive multiple category partitioning with pattern selection
public:
        hpt_typ( ){ assert(!"Illegal constructor"); }
        hpt_typ(char *filename){ Init(filename); }
        hpt_typ(FILE *fp){ Init(fp); }
        hpt_typ(Int4 NumRandom){ 	// create a one node 'hierarchy'.
	   FILE *fp=tmpfile(); fprintf(fp,"HyperParTition:\n!\n+ 1.Set1?\n- 2.Random=%d.\n\n",NumRandom);
	   rewind(fp); Init(fp); fclose(fp);
	}
private:
	hpt_typ	*Next;
	void	Read(FILE *fp,BooLean ignore=FALSE);
public:
	char	*RtnFileName(){ return FileName; }
	void	SetFileName(char *name)
		    { if(FileName) free(FileName); FileName=AllocString(name); }
	void	SetNext(hpt_typ *hpt){ this->Next=hpt; }
	hpt_typ	*RtnNext( ){ return this->Next; }
        // ~hpt_typ( ){ Free( ); }
        ~hpt_typ( ){ if(this->Next) delete this->Next; this->Next=0; Free( ); }
	void    PrintError();
	void	ReNameGroup(Int4 x,char *new_name){
			if(x < 1 || x > NumberBPPS) print_error("hpt_typ input error 0");
			if(GroupName[x]) free(GroupName[x]);
			GroupName[x] = AllocString(new_name);
		}
	void	ReNameSet(Int4 Set,char *new_name){
			if(Set < 1 || Set >= NumElementarySets) print_error("hpt_typ input error 1");
			if(NameElementarySet[Set]) free(NameElementarySet[Set]);
			NameElementarySet[Set] = AllocString(new_name);
		}
	char	*SetName(Int4 x){ return ElmntSetName(x); }
	char	*ElmntSetName(Int4 x){
		   if(x < 1 || x > NumElementarySets){
			fprintf(stderr,"ElmntSetName(): x = %d; NumElementarySets=%d\n",x,NumElementarySets);
		   	assert(x > 0 && x <= NumElementarySets); 
		   } return NameElementarySet[x]; 
		}
	Int4	NumSets(){ return NumElementarySets; }
	BooLean	SetArgStr(Int4 b,char *NewArg){
			if(b < 1 || b > NumberBPPS) return FALSE;
			if(ArgmntStr[b]) free(ArgmntStr[b]);
			ArgmntStr[b]=AllocString(NewArg);
			return TRUE;
		}
	BooLean	SetPttrn(Int4 b,char *NewPttrn){
			if(b < 1 || b > NumberBPPS) return FALSE;
			if(PttrnStr[b]) free(PttrnStr[b]);
			PttrnStr[b]=AllocString(NewPttrn);
			return TRUE;
		}
	set_typ	DescendantSet(Int4 n);
	void	PutClean(FILE *fp){ Put(fp,FALSE,TRUE); } 
	void	Put(FILE *fp){ Put(fp,TRUE,FALSE); } 
	void	Put(FILE *fp,BooLean put_settings, BooLean ignore=TRUE,BooLean putArg=FALSE);
	void	PutHyperPartition(FILE *fp);
	char	Cell(Int4 row, Int4 col){
			if(row < 1 || row > NumElementarySets) return 0;
			if(col < 1 || col > NumberBPPS) return 0;
			return HyperPartition[row][col];
		}
	char	*RtnHyperPartition(Int4 x) {
			assert(x > 0 && x <= NumElementarySets); return HyperPartition[x]; }
	char	**RtnHyperPartition( ){ return HyperPartition; }
	char	RtnHyperPartition(Int4 r,Int4 c) {
			assert(r > 0 && r <= NumElementarySets); assert(c > 0 && c <= NumberBPPS);
			return HyperPartition[r][c]; 
		}
	BooLean	TheSame(hpt_typ *xhpt){
			assert(xhpt->NumSets() == NumElementarySets);
			assert(xhpt->NumBPPS() == NumberBPPS);
			char **HP=HyperPartition,**xHP=xhpt->RtnHyperPartition( );
			for(Int4 x=1; x <= NumElementarySets; x++){
			    if(strcmp(HP[x]+1,xHP[x]+1) != 0) return FALSE;
			     if(this->TypeOfSet(x) != xhpt->TypeOfSet(x)) return FALSE;
			} return TRUE;
		}
	Int4	NodeDepth(Int4 r){
		   assert(r > 0 && r <= NumElementarySets); 
		   Int4 num=0,n;
		   for(n = NumberBPPS; n > 0; n--){
			if(this->Cell(r,n) == '+') num++;
		   } return num;
		}
	BooLean	FixedElmntSet(Int4 x){
		   assert(x > 0 && x <= NumElementarySets); 
		   if(SetType[x] == '=') return TRUE;
		   // else if(SetType[x] == '!') return TRUE;
		   else return FALSE;
		}
	Int4    NumBPPS(){ return NumberBPPS; }
        char    *GrpName(Int4 x){ assert(x > 0 && x <= NumberBPPS);  return GroupName[x]; }
        Int4    *nGrpsFG(){ return nGroupsFG; }
        Int4    *nGrpsBG(){ return nGroupsBG; }
        Int4    **GrpsFG(){ return GroupsFG; }
        Int4    **GrpsBG(){ return GroupsBG; }
        Int4    nGrpsFG(Int4 x){ assert(x > 0 && x <= NumberBPPS);  return nGroupsFG[x]; }
        Int4    nGrpsBG(Int4 x){ assert(x > 0 && x <= NumberBPPS);  return nGroupsBG[x]; }
        Int4    *GrpsFG(Int4 x){ assert(x > 0 && x <= NumberBPPS);  return GroupsFG[x]; }
        Int4    *GrpsBG(Int4 x){ assert(x > 0 && x <= NumberBPPS);  return GroupsBG[x]; }
        Int4    GrpsFG(Int4 x,Int4 y){ assert(x > 0 && x <= NumberBPPS);  
			assert(y > 0 && y <= nGroupsFG[x]);  return GroupsFG[x][y]; }
        Int4    GrpsBG(Int4 x,Int4 y){ assert(x > 0 && x <= NumberBPPS);  
			assert(y > 0 && y <= nGroupsBG[x]);  return GroupsBG[x][y]; }

        Int4    nArg(Int4 x ){ assert(x > 0 && x <= NumberBPPS);  return nArgmnt[x]; }
	void	ReSetArgv(Int4 x,Int4 n, char **argv){
			Int4 i,y; assert(x > 0 && x <= NumberBPPS);
			y = nArgmnt[x]; 
			for(i=0; i < nArgmnt[x]; i++){
			   if(Argmntv[x][i]) free(Argmntv[x][i]);  Argmntv[x][i]=0;
			} if(Argmntv[x]) free(Argmntv[x]);
			nArgmnt[x]=n; NEWP(Argmntv[x],n+3,char);
			for(i=0; i < nArgmnt[x]; i++){ Argmntv[x][i]=AllocString(argv[i]); }
		}
        char	**Argv(Int4 x ){ assert(x > 0 && x <= NumberBPPS);  return Argmntv[x]; }
        char    *sst_str(Int4 x){ assert(x > 0 && x <= NumberBPPS);  return sst_string[x]; }
	Int4	NumberRandom(){ return NumRandom; }
	void	ChangeNumRandom(Int4 x){ assert(x > 0); NumRandom=x; }
	char	RtnMode(){ return Mode; }
	void	UseFormalMode(){ Mode='F'; }
	// ******************* hpt information *******************
	Int4	ColumnsInRow(Int4 row,Int4 &col_neg,Int4 &col_pos,Int4 &col_omit);
	Int4	RowsInColumn(Int4 col,Int4 &row_neg,Int4 &row_pos,Int4 &row_omit);
	//************************* Hpt Operations ****************************
	// ******************* non-tree hpt_operations *******************
	void	Change(char from, char to, Int4 g, Int4 n);
	void    PutSwappedRows(FILE *fp,Int4 s1, Int4 s2);
	hpt_typ	*SwapRows(Int4 s1, Int4 s2){
		   FILE *fp=tmpfile(); PutSwappedRows(fp,s1,s2); rewind(fp); 
		   hpt_typ *h=new hpt_typ(fp); fclose(fp); return h; 
		}
	void    PutSwappedColumns(FILE *fp,Int4 g1, Int4 g2){ PutSwapped(fp,g1,g2); }
	void    PutSwapped(FILE *fp,Int4 g1, Int4 g2);
	hpt_typ *SwapColumns(Int4 g1, Int4 g2){
		   FILE *fp=tmpfile(); PutSwapped(fp,g1,g2); rewind(fp); 
		   hpt_typ *h=new hpt_typ(fp); fclose(fp); return h; 
		}
	void    PutWithout(FILE *fp,Int4 col, Int4 row);

	//-------------- Obsolete operations (eventually remove) --------------
	//-------------- Obsolete operations (eventually remove) --------------
	void    PutInsertGrp(FILE *fp,Int4 gi);
	BooLean	PutAddInternal(FILE *fp,Int4 ParentColRow);
	hpt_typ	*AddInternal(Int4 ParentColRow) {
		   FILE *fp=tmpfile(); BooLean b=PutAddInternal(fp,ParentColRow); 
		   if(b==FALSE){ fclose(fp); return 0; } else rewind(fp); 
		   hpt_typ *h=new hpt_typ(fp); fclose(fp); return h; 
		}
	//-------------- Obsolete operations (eventually remove) --------------

	//************************* Hpt Operations ****************************
        BooLean	OutputThis(Int4 x){ assert(x > 0 && x <= NumberBPPS);  return OutputColumns[x]; }
	//************** New ************
        char	*SetArgV(Int4 x, Int4 a){
		   assert(x > 0 && x <= NumElementarySets);  
		   if(ArgC[x] <= 0) return 0;
		   else { assert(a >= 0 && a < ArgC[x]); return ArgV[x][a]; }
		}
        Int4	RtnArgC(Int4 x){ assert(x > 0 && x <= NumElementarySets);  return ArgC[x]; }
	//************** New ************
	Int4	SameSet(Int4 i, hpt_typ *hptO);
	Int4	PutSettings(FILE *fp,Int4 index,Int4 c);
	Int4	AddEdgeToTree(Int4 col,Int4 row, Int4 root, Int4 EdgeWt, wdg_typ Tree, wdg_typ &NewTree);
	Int4    MergeOldTreeIntoNew(Int4 Root,wdg_typ OT,wdg_typ &Tree);
	BooLean	IsTree(Int4 *&Parent, Int4 *&children);
	BooLean	IsTree(Int4 *&Parent);   // Does the Hpt correspond to a (properly formated) tree?
	BooLean	IsScrambledTree(Int4 *&Parent);	// Is this a tree when sorted?
	BooLean	PutAsTree(FILE *fp);
	wdg_typ	RtnAsTree( );
	BooLean PutAsSmartArt(FILE *fp);

	BooLean PutDFTree(FILE *fp);
	BooLean	RtnDFTree(Int4 *&rtn);
private:
	void    DFTreeSrch(Int4 *rtn, Int4 depth, Int4 v, set_typ S, wdg_typ T);
	// void    DFTreeSrch(FILE *fp, Int4 depth, Int4 v, set_typ S, wdg_typ T);
public:

	void    PrintSmartArtTree(FILE *fp,Int4 root, wdg_typ T);
	void    SmartArtDFS(FILE *fp, Int4 depth, Int4 v, set_typ S, wdg_typ T);
	char    TypeOfSet(Int4 x){ assert(x > 0 && x <= NumElementarySets);  return SetType[x]; }
	Int4	NumInternalNodes(){
		   Int4 s,n=0;
		   for(s=2; s < this->NumElementarySets; s++) if(SetType[s]=='?') n++;
		   return n;
		}
	BooLean	IsParentNode(Int4 x){
			assert(x > 0 && x <= NumElementarySets);
			if(SetType[x]=='?') return TRUE; else return FALSE;  
		}
	void	ChangeTypeOfSet(Int4 x,char X){
			assert(X == '?' || X == '!');
			assert(x > 0 && x <= NumElementarySets);  SetType[x]=X; 
		}
	char    *TypeOfSet( ){ return SetType; }
	char	TypeOfBPPS(Int4 x){ assert(x > 0 && x <= NumberBPPS); return OutputCols[x-1]; }
	void	DeleteBPPS(Int4 x){ assert(x > 0 && x <= NumberBPPS); OutputCols[x-1]='*'; }
	void	DeleteRow(Int4 x){ assert(x > 0 && x <= NumElementarySets); SetType[x]='*'; }
	BooLean	DeletedSet(Int4 x){ assert(x > 0 && x <= NumElementarySets); return (SetType[x]=='*'); }
	// void	FixChildlessInternalNodes( ){ ; }
	BooLean	IsMiscSetForCol(Int4 row, Int4 col){ // assumes correct format...
		   assert(row > 0 && row <= NumElementarySets);
		   assert(col > 0 && col <= NumberBPPS);
		   if(this->TypeOfSet(row) != '?') return FALSE;
		   for(Int4 n = NumberBPPS; n > 0; n--){
			if(this->Cell(row,n) == '+' && n == col) return TRUE;
		   } return FALSE;
		}
	//------------------ Basic operations on trees ----------------
	void    PutInsert(FILE *fp,Int4 level){ PutInsert(fp,level,0); }
	void    PutInsert(FILE *fp,Int4 level,Int4 id);
	hpt_typ *Insert(Int4 x) { return RtnPut(&hpt_typ::PutInsert,x); }
	hpt_typ *Insert(Int4 x,Int4 id) {
			Int4 pID=this->ItoSetID(x);
			hpt_typ *other=RtnPut(&hpt_typ::PutInsert,x,id);
			Int4 pI = other->SetIDtoI(pID);
			if(other->TypeOfSet(pI) != '?') other->ChangeTypeOfSet(pI,'?');
			return other;
		}
	void    PutMoveUp(FILE *fp,Int4 target);
	hpt_typ *MoveUp(Int4 x){ return RtnPut(&hpt_typ::PutMoveUp,x); }
#if 1	// not sure needed; need to specify how going down...think about it.
	void    PutMoveDown(FILE *fp,Int4 Mv, Int4 To);
	hpt_typ *MoveDown(Int4 M,Int4 T){ return RtnPut(&hpt_typ::PutMoveDown,M,T); }
#endif

	void	PutSubTree(FILE *fp, Int4 x);
	hpt_typ	*SubTree(Int4 x){ return RtnPut(&hpt_typ::PutSubTree,x); }

	void	PutAddLeaf(FILE *fp, Int4 x){ PutAddLeaf(fp,x,0); }
	void	PutAddLeaf(FILE *fp, Int4 x,Int4 id);
	hpt_typ	*AddLeaf(Int4 x){ return RtnPut(&hpt_typ::PutAddLeaf,x); }
	hpt_typ	*AddLeaf(Int4 x,Int4 id) {	// id assumes using Set<id> convention.
			Int4 pID=this->ItoSetID(x);
			hpt_typ *other= RtnPut(&hpt_typ::PutAddLeaf,x,id); 
			Int4 pI=other->SetIDtoI(pID);  
			if(other->TypeOfSet(pI) != '?') other->ChangeTypeOfSet(pI,'?');
			return other;
		}
	void	PutDelete(FILE *fp, Int4 x);
	hpt_typ	*Delete(Int4 x){ return RtnPut(&hpt_typ::PutDelete,x); }
	hpt_typ	*MkLineageTree(Int4 x,Int4 &new_x);

	hpt_typ	*Copy( ) { return RtnPut(&hpt_typ::Put); }
	hpt_typ	*CopyClean( ) { return RtnPut(&hpt_typ::PutClean); }
	// ******************* end hpt_operators.cc *******************
	//************************* hpt_sort.cc ****************************
	void	PutRandomTree(FILE *fp,Int4 x,char mode); 
	void	PutRandomize(FILE *fp,Int4 x); 
	hpt_typ	*Randomize(Int4 x){ return RtnPut(&hpt_typ::PutRandomize,x); }
	void	PutRandomFlat(FILE *fp,Int4 x); 
	hpt_typ	*RandomFlat(Int4 x){ return RtnPut(&hpt_typ::PutRandomFlat,x); }
	void    PutSorted(FILE *fp);	// uses private PutRearranged( ) below.
	hpt_typ	*Sort( ) { return RtnPut(&hpt_typ::PutSorted); }
	void	PutShuffled(FILE *fp);
	hpt_typ	*Shuffle( ){ return RtnPut(&hpt_typ::PutShuffled); }
	void    PutRearranged(FILE *fp,BooLean put_settings,Int4 *Row, Int4 *Col);
	hpt_typ *Rearrange(BooLean put_settings,Int4 *Row, Int4 *Col){
		   FILE *fp=tmpfile(); PutRearranged(fp,put_settings,Row,Col); rewind(fp); 
		   hpt_typ *h=new hpt_typ(fp); fclose(fp); return h; 
		}
	// ******************* end hpt_sort.cc *******************
	Int4	ItoSetID(Int4 r){ // WARNING: this assumes that Hpt is a tree!
		    Int4 x,id;
		    x=sscanf(this->ElmntSetName(r),"Set%d",&id);
		    assert(x == 1); return id;
		}
	Int4	SetIDtoI(Int4 id){ // WARNING: this assumes that Hpt is a tree!
		    Int4 x,r;
		    for(r=1; r < this->NumSets(); r++){
			if(sscanf(this->ElmntSetName(r),"Set%d",&x) == 1)
			   { if(x == id) return r; }
		    } return 0;
		}
	BooLean MoveNodeUp(Int4 child);
	Int4	NumDescendants(Int4 n){
		    Int4 N,r;
		    for(N=0,r=1; r < this->NumSets(); r++){
			if(this->Cell(r,n) == '+') N++;
		    } return N-1;
		}
	set_typ	MkLeafSet( ){
		    Int4 i,s,*P;
		    assert(this->IsTree(P));   // Does the Hpt correspond to a (properly formated) tree?
		    set_typ Leaves=MakeSet(this->NumSets()+3); FillSet(Leaves);
        	    for(s=2; s < this->NumSets(); s++){ i=P[s]; if(i > 0) DeleteSet(i,Leaves); }
        	    DeleteSet(1,Leaves); free(P); return Leaves;
		}
	set_typ	MkLineageSet(Int4 c){
		   assert(c > 1 && c < this->NumSets());
		   Int4	*P;
		   assert(IsTree(P));   // Does the Hpt correspond to a (properly formated) tree?
		   set_typ Lineage=MakeSet(this->NumSets()+2); ClearSet(Lineage);
		   for(  ; c != 1; c=P[c]) AddSet(c,Lineage); AddSet(1,Lineage);
		   free(P); return Lineage;
		}
	Int4	NumChildren(Int4 x){
		   assert(x >= 1 && x < this->NumSets());
		   Int4	N,n,*P;
		   assert(IsTree(P));   // Does the Hpt correspond to a (properly formated) tree?
        	   for(N=0,n=2; n < this->NumSets(); n++){ if(P[n]==x) N++; }
		   free(P); return N;
		}
	set_typ	MkSiblingSet(Int4 c){
		    assert(c > 1 && c < this->NumSets());
		    Int4 *P; assert(IsScrambledTree(P));
		    set_typ siblings=MakeSet(this->NumSets()+2); ClearSet(siblings);
		    for(Int4 r=2; r < this->NumSets(); r++){
			if(r == c) continue;
			if(P[r] == P[c]) AddSet(r,siblings); 
		    } free(P); return siblings;
		} 
	set_typ	MkSubTreeSet(Int4 c){
		    set_typ subtree=MakeSet(this->NumBPPS()+2); ClearSet(subtree);
		    for(Int4 r=1; r < this->NumSets(); r++){
			if(this->Cell(r,c) == '+') AddSet(r,subtree); 
		    } return subtree;
		} 
private:
	set_typ	*MkSubTreeSets( ){
		set_typ *subtree=0; NEW(subtree,this->NumBPPS()+2,set_typ);
		// WARNING: this assumes that Hpt is a tree!
		for(Int4 c=1; c <= this->NumBPPS(); c++){
#if 0
		    subtree[c]=MakeSet(Hpt->NumBPPS()+2); ClearSet(subtree[c]);
		    for(Int4 row=1; row < Hpt->NumSets(); row++){
			if(Hpt->Cell(row,c) == '+') AddSet(row,subtree[c]);
		    }
#else
		    subtree[c]=this->MkSubTreeSet(c);
#endif
		} return subtree;
	}
	// void	(hpt_typ::*func)(FILE *, Int4) = &hpt_typ::PutMoveUp;
	hpt_typ *RtnPut(void (hpt_typ::*func)(FILE *)){
		   FILE *fp=tmpfile(); (this->*func)(fp); rewind(fp); 
		   hpt_typ *h=new hpt_typ(fp); fclose(fp); return h; 
		}
	hpt_typ *RtnPut(void (hpt_typ::*func)(FILE *,Int4),Int4 x) {
		   FILE *fp=tmpfile(); (this->*func)(fp,x); rewind(fp); 
		   hpt_typ *h=new hpt_typ(fp); fclose(fp); return h; 
		}
	hpt_typ *RtnPut(void (hpt_typ::*func)(FILE *,Int4,Int4),Int4 x,Int4 id) {
		   FILE *fp=tmpfile(); (this->*func)(fp,x,id); rewind(fp); 
		   hpt_typ *h=new hpt_typ(fp); fclose(fp); return h; 
		}
	// routines for PutSorted( ) function.
	Int4    dfs_childHeap(Int4 &ptr, Int4 *Row, Int4 *Col, Int4 *RowToCol,dh_type *childH, dh_type dH)
		{
		    Int4 r,c;
		    while(!emptyHeap(dH)){
			assert((r=delminHeap(dH)) != 0);
			assert(ptr <= NumberBPPS); ptr++;
			c=RowToCol[r]; Col[ptr]=c; Row[ptr]=r; 
			if(childH[r]) dfs_childHeap(ptr,Row,Col,RowToCol,childH,childH[r]);
        	    }
		}
	void    PrintNewickTree(FILE *fp,Int4 root, wdg_typ T) { TreeDFS(fp,root,T); fprintf(fp,";\n"); }
	void    TreeDFS(FILE *fp, Int4 v, wdg_typ T);
	Int4	RemoveCycle(Int4 cs,wdg_typ &G);
	Int4	RootUnConnectedTree(Int4 Root, wdg_typ &Tree);
	Int4    FindChildOfMRCA(Int4 Parent,Int4 child,Int4 Root,wdg_typ Tree);
	Int4    PathToRoot(Int4 child, Int4 Root, Int4 *&path,set_typ S, wdg_typ Tree);
	wdg_typ CopyGraphWithoutEdge(Int4 E,wdg_typ G_in);
	Int4    GetLastColumn(char symbol, Int4 r){
			for(Int4 c=NumberBPPS; c > 0; c--) if(HyperPartition[r][c]==symbol) return c;
        		return 0;
		}
        //****************** Multiple category routines: ******************
        void    Init(char *);
        void    Init(FILE *);
        void    Free( );
        void    initialize( );
	void	SortSets( );	// Sort sets and remove duplicates.
	void    CreateGroups( );
	void    GetSettings(Int4 N, char *tmp_str);
	void	ValidityCheck();

	//*********************** core structure ***********************
	char	*FileName;
        BooLean	OutputColumns[MAX_NUM_ELMENTARY_SETS];
        char	*OutputCols;		// '!' = ouput; '_' = ignore. 
        BooLean	OutputRows[MAX_NUM_ELMENTARY_SETS];
        char    *ArgmntStr[MAX_NUM_ELMENTARY_SETS];
        char    *PttrnStr[MAX_NUM_ELMENTARY_SETS];
	char	Mode;		// formal ('F') or informal ('I') mode...
        Int4    QryGroup[MAX_NUM_ELMENTARY_SETS];
	Int4    NumRandom;

	//*************** Redundant with cmc_typ ******************8
	char	SetType[MAX_NUM_ELMENTARY_SETS];

	Int4    NumElementarySets;
	char    *NameElementarySet[MAX_NUM_ELMENTARY_SETS];	//

        Int4    nArgmnt[MAX_NUM_ELMENTARY_SETS];
        char    **Argmntv[MAX_NUM_ELMENTARY_SETS];
        char    *sst_string[MAX_NUM_ELMENTARY_SETS];

	Int4    NumberBPPS;
        char    *GroupName[MAX_NUM_ELMENTARY_SETS];
        Int4    nGroupsFG[MAX_NUM_ELMENTARY_SETS];
        Int4    nGroupsBG[MAX_NUM_ELMENTARY_SETS];
        Int4    *GroupsFG[MAX_NUM_ELMENTARY_SETS];
        Int4    *GroupsBG[MAX_NUM_ELMENTARY_SETS];

        // char    *HyperPartition[MAX_NUM_ELMENTARY_SETS];
        char    **HyperPartition;   // HyperPartition[grp][bpps]: FG='+'; BG='-'; RM='o'.
	//*************** End Redundant with cmc_typ ******************8
        char	**ArgV[MAX_NUM_ELMENTARY_SETS];
        Int4	ArgC[MAX_NUM_ELMENTARY_SETS];
};

hpt_typ *MultiReadHpt(Int4 &N, char *filename);
void    MultiPutHpt(FILE *fp, hpt_typ *head);

#endif


