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

#include "hpt_typ.h"

BooLean hpt_typ::IsScrambledTree(Int4 *&Parent)
// Is this a tree when sorted?
{
         Int4 *Ps=0,*Pt=0,*T2S=0,*S2T=0;
         hpt_typ *S=this->Sort();
         if(!S->IsTree(Ps)){ if(Ps) free(Ps); delete S; return FALSE; }
         else {
            Int4 s,t,n,N=this->NumSets(); NEW(T2S,N+3,Int4); NEW(S2T,N+3,Int4);
            for(n=0,s=1; s <= N; s++){
               char *NameS=S->NameElementarySet[s];
               for(t=1; t <= N; t++){
                  char *NameT=this->NameElementarySet[t];
                  if(strcmp(NameS,NameT) == 0){
			assert(T2S[t]==0); T2S[t]=s; S2T[s]=t; n++; break; 
		  }
               }
            } assert(n==N); NEW(Pt,N+3,Int4);
            for(t=1; t <= N; t++){ assert(T2S[t]!=0); s=T2S[t]; Pt[t]=S2T[Ps[s]]; }
	    Parent=Pt; free(Ps); free(T2S); free(S2T); delete S;
         } return TRUE;
}


void	hpt_typ::PutSorted(FILE *fp)
// Puts hpt in canonical tree format (assumes hpt corresponds to a tree!)...
// Assume first and last rows and first columns are unchanged!!!
#if 0
	STRATEGY:  Starting with the root node, 
  Hints:
    Columns:
	* Root column == 1. (Doesn't change).
	* Leaf columns will contain only one '+' symbol.
	* Internal columns will contain one '+' symbol for each node in its subtree.
	// * A leaf's parent column will have as many '+' symbols as that leaf's column has '+' & '-'. 
	// * A column with a '-' in row 1 is a child of the Root node. 
	// * A leaf column will contain a '-' symbol for sibling and parent rows. 
    Rows:
	* Root row == 1. (Doesn't change).
	* A node's level in the tree == the number of '+' symbols.
#endif
{
	FILE *dbfp=0; // dbfp=stderr;
	Int4	r,c,R,C,neg,plus,omit,NR=NumElementarySets-1,NC=NumberBPPS; 
	if(NC!=NR){ fprintf(stderr,"NC = %d; NR = %d\n",NC,NR); assert(NC == NR); }
	// ignore Random set in NR (i.e., row).

	// 1. Reality check: 
	assert(this->TypeOfSet(1) == '?'); 
	this->ColumnsInRow(1,neg,plus,omit); assert(plus==1);

	// 2. Gather information about rows. 
	dh_type *childRow; NEW(childRow,NR+3,dh_type);
	Int4	*Level; NEW(Level,NR+3,Int4);	// level within heirarchy...
	for(r=1; r <= NR; r++){ 	// ignore Random set (i.e., row).
	   this->ColumnsInRow(r,neg,plus,omit); Level[r]=plus; 	// Level in tree; Root == 1.
	}

	// 3. Gather information about columns. 
	Int4	*ColPlus; NEW(ColPlus,NR+3,Int4);
	for(c=1; c <= NC; c++){ this->RowsInColumn(c,neg,plus,omit); ColPlus[c]=plus; }

	// 4. Find the correspondence between Rows and Columns...
	Int4	*RowToCol,*ColToRow; NEW(RowToCol,NR+3,Int4); NEW(ColToRow,NC+3,Int4);
	// RowToCol[1]=1; ColToRow[1]=1;	// Which row map to which columns and vice versa.
	for(c=1; c <= NC; c++){
	     if(ColPlus[c] > 1){	// Internal node column...
		Int4 the_row=0,the_level=NR+3;
	        for(r=1; r <= NR; r++){		// find the corresponding row.
		    if(this->Cell(r,c) == '+'){  
			if(Level[r] < the_level){ the_level=Level[r]; the_row=r; }
		    }
		}  assert(RowToCol[the_row] == 0 && ColToRow[c]== 0);
		RowToCol[the_row]=c; ColToRow[c]=the_row;
	     } else { // The one '+' row in this column == the corresponding leaf column.
	        for(r=1; r <= NR; r++){
		   if(this->Cell(r,c) == '+'){ RowToCol[r]=c; ColToRow[c]=r; break; }
		}
	     }
	}
	if(dbfp){
	  fprintf(dbfp,"\n");
	  for(r=1; r <= NR; r++) fprintf(dbfp,"row %d%c (%d) <--> col %d\n",
				r,this->TypeOfSet(r),Level[r],RowToCol[r]);
	  fprintf(dbfp,"\n");
	  for(c=1; c <= NC; c++) fprintf(dbfp,"col %d (%d '+') <--> row %d\n",c,ColPlus[c],ColToRow[c]);
	  fprintf(dbfp,"\n\n");
	}

	// 5. Find leaf and internal rows and find the child nodes.
	for(r=1; r <= NR; r++){ 	
	   if(ColPlus[RowToCol[r]] == 1){ 
	     if(this->TypeOfSet(r) == '?'){
		fprintf(stderr,"\nr = %d ('Set%d'); RowToCol[r] = %d; set type = '%c'\n",
			r,ItoSetID(r),RowToCol[r],this->SetType[r]);
		this->Put(stderr,FALSE);
	 	print_error("hpt_typ::PutSorted( ) input error: '?' node without child nodes");
		assert(this->TypeOfSet(r) != '?'); 
	     }
	   } else {	// find the direct child nodes.
	        if(this->TypeOfSet(r) != '?'){
			fprintf(stderr,"r = %d; RowToCol[r] = %d; type = '%c'\n",
				r,RowToCol[r],this->TypeOfSet(r));
			assert(this->TypeOfSet(r) == '?'); 
		}
		childRow[r]=dheap(NR+2,4); 	// childRow[c]== 0 --> leaf row.
		if(dbfp) fprintf(dbfp,"row %d(%d): ",r,Level[r]);
		C=RowToCol[r];	
		for(c=1; c <= NC; c++){
		   R=ColToRow[c];	// ancestor or descendant node
		   if(this->Cell(R,C) == '+'){
			if(Level[r] == Level[R]-1){	// next level in the hierarchy...
				if(dbfp) fprintf(dbfp," %d(%d)",R,Level[R]);
				insrtHeap(R,-(keytyp)ColPlus[C],childRow[r]);
			}
		   }
		} if(dbfp) fprintf(dbfp,"\n");
	   }
	} if(dbfp) fprintf(dbfp,"\n\n");

	// 6. Find the mapping from the shuffled to canonical hpt format.
	Int4	ptr,*Row,*Col; NEW(Row,NR+3,Int4); NEW(Col,NC+3,Int4);
	ptr=1; Row[1]=1; Col[1]=1; Row[NR+1]=NR+1; 	// Root & Reject rows & columns don't change.
	dfs_childHeap(ptr,Row,Col,RowToCol,childRow,childRow[1]);

	if(dbfp){
	  for(r=1; r <= NR+1; r++){ fprintf(dbfp,"Row[%d] <--> row %d.\n",r,Row[r]); } fprintf(dbfp,"\n");
	  for(c=1; c <= NC; c++){ fprintf(dbfp,"Col[%d] <--> col %d.\n",c,Col[c]); } fprintf(dbfp,"\n\n");
	}

	// 7. Print canonical hpt and deallocate memory.
	free(RowToCol); free(ColToRow); free(ColPlus);
	for(r=1; r <= NR; r++) if(childRow[r]) Nildheap(childRow[r]); free(childRow);
	free(Level);
	this->PutRearranged(fp,FALSE,Row,Col);  // Put Hpt in canonical tree form.
	free(Row); free(Col);
}

void    hpt_typ::PutShuffled(FILE *fp)
// Randomly permute the rows and columns (except for Root & Reject).
{
	Int4 i,j,NR=NumElementarySets-1,NC=NumberBPPS; assert(NC==NR);
        Int4 *Row,*Col; NEW(Row,NR+3,Int4); NEW(Col,NC+3,Int4);
        Row[1]=1; Row[NR+1]=NR+1; Col[1]=1;
        dh_type dH=dheap(NR+2,4);
        for(i=2; i <= NR; i++) insrtHeap(i,(keytyp)Random(),dH);
        for(j=2; (i=delminHeap(dH)) != 0; j++) Row[j]=i;
        for(i=2; i <= NC; i++) insrtHeap(i,(keytyp)Random(),dH);
        for(j=2; (i=delminHeap(dH)) != 0; j++) Col[j]=i;
        PutRearranged(fp,FALSE,Row,Col); free(Row); free(Col);
        Nildheap(dH);
}

void	hpt_typ::PutRearranged(FILE *fp,BooLean put_settings,Int4 *Row, Int4 *Col)
// Output the Hpt using the ordering given by Row and Col arrays.
{
	Int4	R,C,s,i,j,n,col;
	if(Mode == 'I') fprintf(fp,"\nHyperpartition:\n"); else fprintf(fp,"\nHyperParTition:\n");
	for(n=0; n < NumberBPPS; n++) fprintf(fp,"%c",OutputCols[n]); fprintf(fp,"\n");

	for(s=1; s <= NumElementarySets; s++){
	   R=Row[s];
           for(n=1; n <= NumberBPPS; n++){
	      C=Col[n];
              if(HyperPartition[R][C] != 0){ fprintf(fp,"%c",HyperPartition[R][C]); }
              else { fprintf(fp,"o"); }
           } 
	   if(SetType[R] == '=') fprintf(fp," %d.%s=%d.",s,NameElementarySet[R],NumRandom);
	   else fprintf(fp," %d.%s%c",s,NameElementarySet[R],SetType[R]);
	   for(Int4 a=0; a < ArgC[R]; a++) fprintf(fp," %s",ArgV[R][a]); 
	   fprintf(fp,"\n");
	} fprintf(fp,"\n\n"); 

	if(put_settings && i > 0){
	  fprintf(fp,"\nSettings:\n");
	  for(n=1; n <= NumberBPPS; n++){
	    C=Col[n];
	    if(nArgmnt[n] > 0) fprintf(fp,"%d.%s\n",n,ArgmntStr[C]);
	    else fprintf(fp,"%d.Column%d\n",n,C);
	  } fprintf(fp,"\n\n");
	}
}

void    hpt_typ::PutRandomize(FILE *fp, Int4 x) { PutRandomTree(fp,x,'A'); }

void    hpt_typ::PutRandomFlat(FILE *fp, Int4 x) { PutRandomTree(fp,x,'F'); }

void    hpt_typ::PutRandomTree(FILE *fp, Int4 x, char mode)
// Add a new leaf to node x.
{
        Int4 *P=0,i,j,n,m,op,N,M,p;
	enum operation {addleaf=0,insrt=1,del=2,mvup=3,mvdwn=4};
	char *routine[] = {"addleaf","insert","delete","moveup","movedown"};
	hpt_typ *tmp_hpt=0,*Hpt=0;

        if(!IsTree(P)){ if(P) free(P); print_error("PutRandomize() error: Hpt not in tree format"); }
	Hpt=this->Copy();
	h_type	HG= Histogram("hpt operations",-1,10,1);
	fprintf(stderr,"PutRandomTree() mode = '%c'\n",mode);
	for(i=1; i <= x; i++){
	   switch (mode){
	     case 'A': n=(Random()%Hpt->NumberBPPS) + 1; op=0; break;
	     case 'F': n=1; op=0; break;	// flat tree.
	     default:
#if 0
	   	if(Hpt->NumSets() <= 15) op = 0; else op=Random()%7; 
	   	if(Hpt->NumSets() <= 8) n=1; else n=(Random()%Hpt->NumberBPPS) + 1;
#else
	   	op=Random()%7; n=(Random()%Hpt->NumberBPPS) + 1;
#endif
	       break;
	   } 
	   if(op >= 5) op = 0;
	   if(n == 1 && op == insrt){ i--; continue; }
	   if(n == 1 && op == mvdwn){ i--; continue; }
	   if(op == del && Hpt->NumSets() <= 3){ i--; continue; }
	   if(n == 1 && op != addleaf && op != insrt){ i--; continue; }
	   if(Hpt->SetType[n] == '!' && op == insrt){ i--; continue; }
	   if(P[n] == 1 && op == mvup){ i--; continue; }
	   if(op == mvdwn){ 	// find a sibling node to n
		set_typ siblings=Hpt->MkSiblingSet(n); N=CardSet(siblings);
		if(N == 0){ i--;  NilSet(siblings); continue; } else p = (Random()%N) + 1;
		for(m=0,M=0,j=2; j <= Hpt->NumberBPPS; j++){
			if(MemberSet(j,siblings)) M++;
			if(M >= p){ m=j; break; }
		} assert(m > 0); NilSet(siblings);
	   }
if(op == del){ i--; continue; }
if(op == mvup || op == mvdwn){ i--; }
	   fprintf(stderr," === %d: operation = %s node=%d; parent=%d",i,routine[op],n,P[n]); 
	   if(op == mvdwn) fprintf(stderr,"; move to node %d ===\n",m); 
	   else fprintf(stderr," ===\n");
	   switch ((operation) op){
		case addleaf: tmp_hpt=Hpt->AddLeaf(n); break;
		case insrt: tmp_hpt=Hpt->Insert(n); break;
		case del: tmp_hpt=Hpt->Delete(n); break;
		case mvup: tmp_hpt=Hpt->MoveUp(n); break;
		case mvdwn: tmp_hpt=Hpt->MoveDown(n,m); break;
		default: print_error("hpt_typ::PutRandomize(): this should not happen"); break;
	   } delete Hpt; Hpt=tmp_hpt; 
	   IncdHist(op, HG);
	     // if(op==mvdwn) Hpt->Put(stderr); 
	   tmp_hpt=Hpt->Sort( ); delete Hpt; Hpt=tmp_hpt; if(P) free(P); 
	   if(!Hpt->IsTree(P)){ if(P) free(P); print_error("PutRandomize() error: Hpt non-tree format"); }
	   if(Hpt->NumSets() > x) break;
	} Hpt->Put(fp); delete Hpt; if(P) free(P);
	PutHist(stderr,60,HG); NilHist(HG);
}



