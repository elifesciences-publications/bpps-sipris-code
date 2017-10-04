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
#include "blosum62.h"

static Int4    dummyID=1;

void	hpt_typ::PutMoveDown(FILE *fp,Int4 Mv, Int4 To)
{
	Int4	N,*Parent=0,x=Mv,c,p=To,gp,r;
	hpt_typ *hpt=this->Copy();
        if(!hpt->IsTree(Parent)) print_error("hpt_typ PutMoveDown() error: hpt is not a tree");
        if(x < 2 || x > hpt->NumSets()) print_error("hpt_typ PutMoveDown() error: node 1 is not in tree");
        if(p < 2 || p > hpt->NumSets()) print_error("hpt_typ PutMoveDown() error: node 2 is not in tree");
        if(Parent[x] != Parent[p]) print_error("hpt_typ PutMoveDown() error: Move Down operation not possible");
	char    **HP=this->RtnHyperPartition();
	set_typ subtree=this->MkSubTreeSet(x);
#if 0
	gp=Parent[p];
	for(N=0,c=2; c < this->NumSets(); c++){
		if(c == p) continue;
		if(Parent[c] == gp) N++;
	} if(N == 0){ hpt->ChangeTypeOfSet(gp,'!'); }	// change this to a leaf node.
#endif
	gp=Parent[p];
        if(this->TypeOfSet(p) == '!'){ hpt->ChangeTypeOfSet(p,'?'); }
	
	for(c=1; c <= hpt->NumBPPS(); c++){
	    if(c == gp) continue;	// hpt is the current parent/new grandparent.
	    if(c == p){	// hpt is the new parent...
	        for(r=1; r <= hpt->NumBPPS(); r++){	// over all rows == sets...
	           if(r==c) continue;	// diagonal does not change...
		   if(MemberSet(r,subtree)){
			assert(HP[r][c] == '-'); 
			hpt->Change('-','+',r,c);
		   }
		}
	    } else if(c == x){
		for(r=1; r <= hpt->NumBPPS(); r++){	// over all rows == sets...
	           if(r==c) continue;	// diagonal does not change...
		   if(MemberSet(r,subtree)) continue; // child's descendants don't change..
		   if(HP[r][gp] == '+' && HP[r][p] == '-'){ // r is a current sibling...
			assert(HP[r][c] == '-'); 	// so that c is in the background.
			hpt->Change('-','o',r,c);
		   }
		}
	    } else if(Parent[c] == p){
	        if(!MemberSet(c,subtree)){ // New sibling of child.
		  for(r=1; r <= hpt->NumBPPS(); r++){	// then need to remove these from BG.
		    if(r==c) continue;
		    if(MemberSet(r,subtree)){
			assert(HP[r][c] == 'o');
			// fprintf(stderr," Moving %d down from %d to %d (%d,%d)\n",ch,gp,p,r,c);
			hpt->Change('o','-',r,c);
		    }
		  }
		}
	    } // else do nothing
	} NilSet(subtree); free(Parent);
	hpt->PutSorted(fp); delete hpt;
}

void	hpt_typ::PutAddLeaf(FILE *fp, Int4 target,Int4 id)
// Add a leaf to node x.
{
	Int4	s,i,j,n,S,col,parent,*Parent;
	char    state,**HP=HyperPartition;

	if(this->IsTree(Parent) == FALSE){
                if(Parent) free(Parent);
                print_error("PutAddLeaf() input error: Hpt not in tree format ");
        }
	if(target < 1 || target > NumberBPPS) print_error("PutAddLeaf() error: 'target' node out of range.");
	parent = Parent[target];

	if(Mode == 'I') fprintf(fp,"\nHyperpartition:\n"); else fprintf(fp,"\nHyperParTition:\n");
	for(n=1; n <= NumberBPPS; n++){ fprintf(fp,"%c",OutputCols[n-1]); if(n == target) fprintf(fp,"!"); }
	fprintf(fp,"\n");
	for(S=1,s=1; s <= NumElementarySets; S++,s++){
           for(n=1; n <= NumberBPPS; n++){
	      fprintf(fp,"%c",HP[s][n]);
	      if(n == target){				// We are within target's column.
	    	if(s == NumElementarySets) fprintf(fp,"o");
		else if(s==target) fprintf(fp,"-");     // make new node a child of the target node.
		// else if(Parent[n]==target) fprintf(fp,"-");
		else if(HP[s][target] == '+') fprintf(fp,"-");
		else fprintf(fp,"o");
	      } // else fprintf(fp,"%c",HP[s][n]);
           } 
	   if(SetType[s] == '=') fprintf(fp," %d.%s=%d.",S,NameElementarySet[s],NumRandom);
	   else {
		fprintf(fp," %d.%s",S,NameElementarySet[s]);
	   	if(s==target) fprintf(fp,"?"); else fprintf(fp,"%c",SetType[s]);
	   }
	   for(Int4 a=0; a < ArgC[s]; a++) fprintf(fp," %s",ArgV[s][a]); 
	   fprintf(fp,"\n");
	   if(s == target){      // this is the new node right after target node...
                S++;
                for(n=1; n <= NumberBPPS; n++){
                    if(n == target) fprintf(fp,"++");
                    else if(Parent[n] == target) fprintf(fp,"-");
                    else if(TRUE) fprintf(fp,"%c",HP[s][n]);	// same as parent at other column positions.
                    else if(HP[s][n] == '+') fprintf(fp,"-"); else fprintf(fp,"o");
                } 
		if(id < 1) fprintf(fp," %d.dummy%d!\n",S,dummyID++); 
		else fprintf(fp," %d.Set%d!\n",S,id);
           }
	}
	if(NumberBPPS > 5){
	 for(col=i=0,n=1; n <= NumberBPPS; n++){
	   if(OutputCols[n-1] == '*') continue; else col++;
	   if(col==1) fprintf(fp,"#"); else if(col%5 ==0) fprintf(fp,"|"); else fprintf(fp," ");
	   if(nArgmnt[n] > 0) i++; 
	 }
	 fprintf(fp,"\n");
	 for(n=1; n <= col; n++){
		if(n==1) fprintf(fp,"#"); else if(n%5 ==0) fprintf(fp,"%4d ",n);
	 } fprintf(fp,"\n");
	} else fprintf(fp,"\n"); 
	if(Parent) free(Parent);
}

void    hpt_typ::PutDelete(FILE *fp, Int4 x)
// Delete target node 'x' from the tree; attach any child nodes to parent of x.
{
	Int4 *P=0,neg,plus,omit,s,n,g,p;
	hpt_typ *hpt=this->Copy();
	char    **HP=hpt->HyperPartition;
        if(!hpt->IsTree(P)){ if(P) free(P); print_error("PutDelete() error: Hpt not in tree format"); }
	if(x < 2 || x > NumberBPPS) print_error("PutDelete() error: input node out of range");

	hpt->RowsInColumn(P[x],neg,plus,omit); 
	if(plus == 2) hpt->SetType[P[x]] = '!';

	hpt->RowsInColumn(x,neg,plus,omit); 
	if(plus > 1){	// then an internal node was deleted; attach x's child nodes to Parent node.
	   p=P[x];
	   for(n=1; n <= NumberBPPS; n++){ 
	     if(P[n] == x){	// n is a child of x.
	       for(s=1; s < NumElementarySets; s++){
		 if(HP[s][n] == 'o' && HP[s][p] == '+') HP[s][n] = '-';
	       }
	     }
	   }
	}
        hpt->DeleteBPPS(x); hpt->DeleteRow(x); 
	free(P); hpt->Put(fp,FALSE); delete hpt;
}

void	hpt_typ::PutMoveUp(FILE *fp,Int4 target)
// Changes the parent of the target node to the grandparent node.
{
	Int4	c,r,s,n,N,i,j,g,*Parent=0,parent,grandpa;
	char	state,**HP=HyperPartition;

	if(this->IsTree(Parent) == FALSE){
		if(Parent) free(Parent);
		print_error("PutMoveUp() input error: Hpt not in tree format "); 
	}
	if(target < 1 || target > NumberBPPS) print_error("PutMoveUp() input error: 'target' node out of range.");
	// if(SetType[target] == '?'){ this->Put(fp);  return; }	// don't add anything...

	if(target == 1) print_error("PutMoveUp() input error 2");
	if(Parent[target] == 1) print_error("PutMoveUp() input error 3");
	parent=Parent[target]; grandpa=Parent[parent];
	set_typ Sibling=MakeSet(NumElementarySets+3); ClearSet(Sibling);
	for(s=1; s < NumElementarySets; s++) if(Parent[s] == parent) AddSet(s,Sibling);

	if(Mode == 'I') fprintf(fp,"\nHyperpartition:\n"); else fprintf(fp,"\nHyperParTition:\n");
	fprintf(fp,"%s\n",OutputCols);
	for(s=1; s <= NumElementarySets; s++){
	  if(HP[s][target] == '+') { // s is target or a current child of target.
	    for(n=1; n <= NumberBPPS; n++){
	      if(n != target && MemberSet(n,Sibling)){ fprintf(fp,"o"); }
	      else if(n == parent) fprintf(fp,"-"); // target is now a sibling to parent.
	      else if(Parent[n] == grandpa) fprintf(fp,"-");
	      else fprintf(fp,"%c",HP[s][n]);
	    }
	  } else {		// other rows.
	    for(n=1; n <= NumberBPPS; n++){
	      // if(HP[s][target] == '+') { // s is target or a current child of target.
	      if(n == target){	// target column within non-target row.
		if(n==grandpa) { fprintf(fp,"-"); 	// make new node the sole child of target node.
		// } else if(HP[s][parent] == '+') fprintf(fp,"+"); // s is a current sibling of target.
		} else if(HP[s][grandpa] == '+') {
			fprintf(fp,"-");
	        } else fprintf(fp,"%c",HP[s][n]);
	      } else fprintf(fp,"%c",HP[s][n]);
	    }
	  }
	  if(SetType[s] == '=') fprintf(fp," %d.%s=%d.",s,NameElementarySet[s],NumRandom);
	  else {
		if(s == parent){  // see whether current parent has only one child == target
		    for(c=0,r=1; r < NumElementarySets; r++){ // count the remaining children of parent.
			if(r == parent || r == target) continue;
			if(HP[r][parent] == '+' && HP[r][target] != '+') c++;
			//                         ^ make sure this node in not in target's subtree!
		    }
		    if(c == 0) fprintf(fp," %d.%s!",s,NameElementarySet[s]);
		    else fprintf(fp," %d.%s?",s,NameElementarySet[s]);
		} else fprintf(fp," %d.%s%c",s,NameElementarySet[s],SetType[s]);
	  }
	  for(Int4 a=0; a < ArgC[s]; a++) fprintf(fp," %s",ArgV[s][a]); 
	  fprintf(fp,"\n");
	} fprintf(fp,"\n");
	if(Parent) free(Parent); NilSet(Sibling);
}

BooLean	hpt_typ::MoveNodeUp(Int4 child)
// move the node 'child' from its parent to its grand parent.
{
	Int4	gp,p,ch,r,c,x,sq,*Parent=0;
	hpt_typ *oHpt=this->Copy();
	char    state,**HP=oHpt->RtnHyperPartition();
	assert(NumSets() == NumberBPPS +1); // need direct correspondence between column & rows!
	assert(this->IsScrambledTree(Parent));

	assert(child > 1 && child <= NumberBPPS);
	ch=child; p=Parent[ch]; assert(p > 1); gp = Parent[p]; assert(gp > 0);
	set_typ subtree=this->MkSubTreeSet(ch);


	// if(this->TypeOfSet(p) == '!'){ this->ChangeTypeOfSet(p,'?'); }
	state=this->Cell(ch,p);	
	assert(state == '-' || state == '+');	// must be the case if parent in a tree.

	// fprintf(stderr,"Now moving %d up from %d to %d\n",c,p,gp);
	for(c=1; c <= NumberBPPS; c++){
	    if(c == gp) continue;
	    if(c == p){				// current parent...
		// fprintf(stderr," ... %d == %d\n",c,p);
		for(r=1; r <= NumberBPPS; r++){	// over all rows == sets...
	           if(r==c) continue;	// diagonal does not change...
		   if(MemberSet(r,subtree)){
			assert(HP[r][c] == '+'); 
			// fprintf(stderr," ... moving %d up from %d to %d\n",ch,p,gp);
			this->Change('+','-',r,c);
		   }
		}
	    } else if(c == child){	// child column
		for(r=1; r <= NumberBPPS; r++){	// over all rows == sets...
	           if(r==c) continue;	// diagonal does not change...
		   if(MemberSet(r,subtree)) continue; // child & descendants don't change..
		   if(HP[r][gp] == '+' && HP[r][p] == '-'){
			assert(HP[r][c] == 'o'); 
			this->Change('o','-',r,c);
		   }
		}
	    } else if(!MemberSet(c,subtree) && Parent[c] == p){
		for(r=1; r <= NumberBPPS; r++){	// over all rows == sets...
		    if(r==c) continue;
		    if(MemberSet(r,subtree)){
			assert(HP[r][c] == '-');
			this->Change('-','o',r,c);
		    }
		}
	    } // else do nothing...
	} free(Parent); delete oHpt; NilSet(subtree);
	return TRUE;
}

void	hpt_typ::PutInsert(FILE *fp,Int4 target,Int4 id)
// attach a node to the target in tree and assign target child nodes to the new node.
// if NumElementarySets - NumberBPPS  == 1 then add a row as well.
{
	Int4	c,r,s,S,n,N,i,j,g,*Parent=0;
	char	state,**HP=HyperPartition;

	if(this->IsTree(Parent) == FALSE){
		if(Parent) free(Parent);
		print_error("PutInsert() error: Hpt not in tree format "); 
	}
	if(target < 1 || target > NumberBPPS) print_error("PutInsert() error: 'target' node out of range.");
#if 0	// allow insertions in leaf nodes.
	if(SetType[target] != '?'){ this->Put(fp);  return; }	// don't add anything...
#endif

	if(Mode == 'I') fprintf(fp,"\nHyperpartition:\n"); else fprintf(fp,"\nHyperParTition:\n");
	for(c=0,g=1; c < NumberBPPS; c++,g++){ fprintf(fp,"%c",OutputCols[c]); if(g == target) fprintf(fp,"!"); }
	fprintf(fp,"\n");
	for(s=S=1; s <= NumElementarySets; S++,s++){
	  for(n=1; n <= NumberBPPS; n++){
	    if(n == target){	// add a node to target node and attach target's child nodes to it...
		fprintf(fp,"%c",HP[s][n]);
		if(s == NumElementarySets) fprintf(fp,"o"); 
		else if(s==target){
			fprintf(fp,"-"); 	// make new node the sole child of target node.
		} else {
		  switch(HP[s][target]){
		    case '+': fprintf(fp,"+"); break;
		    case '-': fprintf(fp,"o"); break;
		    case 'o': fprintf(fp,"o"); break;
		    default: print_error("PutInsert() switch error");
		  } 
		}
	    } else if(s== target && Parent[n] == target){
			fprintf(fp,"o"); 	// make new node the sole child of target node.
	    } else fprintf(fp,"%c",HP[s][n]);
	  }
	  if(SetType[s] == '=') fprintf(fp," %d.%s=%d.",S,NameElementarySet[s],NumRandom);
	  else fprintf(fp," %d.%s%c",S,NameElementarySet[s],SetType[s]);
	  for(Int4 a=0; a < ArgC[s]; a++) fprintf(fp," %s",ArgV[s][a]); 
	  fprintf(fp,"\n");
	  if(s == target){	// this is the new node right after target node...
		S++;
		for(n=1; n <= NumberBPPS; n++){
		    if(n == target) fprintf(fp,"++");
		    else if(Parent[n] == target){
			fprintf(fp,"-");
		    } else if(TRUE){ fprintf(fp,"%c",HP[s][n]);
		    } else if(HP[s][n] == '+'){
			fprintf(fp,"-");
		    } else fprintf(fp,"o"); 
		} 
		if(id < 1) fprintf(fp," %d.dummy%d?\n",S,dummyID++); 
		else fprintf(fp," %d.Set%d?\n",S,id);
	  }
	}
#if 0
	fprintf(fp,"\nSettings:\n"); N=1;
	for(n=1; n <= NumberBPPS; n++,N++){
	    fprintf(fp,"%d.%s\n",N,ArgmntStr[n]);
	    if(n == KeySet){ N++; fprintf(fp,"%d.Dummy \n",N); }
	} fprintf(fp,"\n\n");
#else
	fprintf(fp,"\n");
#endif
	if(Parent) free(Parent); 
}

BooLean	hpt_typ::PutAddInternal(FILE *fp,Int4 gi)
// insert a dummy group directly before column gi
// if NumElementarySets - NumberBPPS  == 1 then add a row as well.
{
	Int4	c,s,S,n,N,i,j,g;
	BooLean	InsertRow=FALSE;

	assert(gi >= 0 && gi <= NumberBPPS); 
	assert(NumElementarySets - NumberBPPS  == 1);

	// 1. Check to ensure that this is being added directly AFTER an internal node.
	if(gi <= 0) return FALSE;
	if(this->TypeOfSet(gi) != '?') return FALSE;	// must be an internal node...
	for(n=0,s=1; s < NumElementarySets; s++){
	   if(HyperPartition[s][gi] == '+'){ n++; }
	} if(n < 4) return FALSE;			// ..with at least 3 child nodes.
	
	// 2. Print out first 2 lines.
	if(Mode == 'I') fprintf(fp,"\nHyperpartition:\n"); else fprintf(fp,"\nHyperParTition:\n");
	for(c=0,g=1; c < NumberBPPS; c++,g++){
	    fprintf(fp,"%c",OutputCols[c]);
	    if(g == gi){ fprintf(fp,"!"); }
	} fprintf(fp,"\n");

	// 3. Print out the hyperpartition with new row added...
	for(s=S=1; s <= NumElementarySets; S++,s++){
	  // if(gi == 0){ fprintf(fp,"+"); }
	  for(n=1; n <= NumberBPPS; n++){
	    fprintf(fp,"%c",HyperPartition[s][n]);
	    if(n == gi){	// add the new column to existing rows here...
		if(s == NumElementarySets) fprintf(fp,"o"); 
		else if(s == gi) fprintf(fp,"-"); 	// parent node...
		else if(HyperPartition[s][n] == '+') fprintf(fp,"-");
		else fprintf(fp,"o"); 
	    }
	  }
	  if(SetType[s] == '=') fprintf(fp," %d.%s=%d.",S,NameElementarySet[s],NumRandom);
	  else fprintf(fp," %d.%s%c",S,NameElementarySet[s],SetType[s]);
	  for(Int4 a=0; a < ArgC[s]; a++) fprintf(fp," %s",ArgV[s][a]); 
	  fprintf(fp,"\n");
	  if(s == gi){		// Add an additional row here...
		S++;
		// fprintf(fp,"+");	// root node always '+'.
		for(n=1; n <= NumberBPPS; n++){
		    fprintf(fp,"%c",HyperPartition[gi][n]); // same as parent node...
		    if(n == gi) fprintf(fp,"+");
		    // else if(HyperPartition[s][gi] == '+') fprintf(fp,"-");
		    // else fprintf(fp,"o"); 
		}
		fprintf(fp," %d.dummy%d?\n",S,dummyID++); 
	  }
	}

	// 4. Print out the Settings.
	fprintf(fp,"\nSettings:\n"); N=1;
	if(gi == 0){ fprintf(fp,"%d.Dummy \n",N); N++; }
	for(n=1; n <= NumberBPPS; n++,N++){
	    fprintf(fp,"%d.%s\n",N,ArgmntStr[n]);
	    if(n == gi){ N++; fprintf(fp,"%d.Dummy \n",N); }
	} fprintf(fp,"\n\n");
	return TRUE;
}


void	hpt_typ::Change(char from, char to, Int4 g, Int4 n)
{
   assert((from == '-' && to == '+') || (from == 'o'&& to == '-') || (from == '-'&& to == 'o') ||
		(from == '+'&& to == '-') || (from == '+'&& to == 'o') || (from == 'o'&& to == '+'));
   assert(g > 0 && g <= NumElementarySets); assert(n > 0 && n <= NumberBPPS);
   Int4	i,j,N,x;
   if(from == '+'){	// from '+' to '-'.
    switch(to){
     case '-':	// from '+' to '-'.
     {
       assert(HyperPartition[g][n] == '+'); HyperPartition[g][n] = '-';
       N=nGroupsFG[n];
       for(i=1; i <= N; i++){
	  if(GroupsFG[n][i] == g){
		for(j=i; j < N; j++) GroupsFG[n][j] = GroupsFG[n][j+1];
		GroupsFG[n][j]=0;  nGroupsFG[n]--; break;
	  }
       } assert(i <= N);
       N=nGroupsBG[n];
       for(i=1; GroupsBG[n][i] < g && i <= N; i++) ;
       assert(GroupsBG[n][i] != g);
       for(j=N; j >= i; j--) GroupsBG[n][j+1]=GroupsBG[n][j];  // make room for g.
       GroupsBG[n][i]=g; nGroupsBG[n]++;
     } break; 
     case 'o':	// from '+' to 'o'.
     {
       assert(HyperPartition[g][n] == '+'); HyperPartition[g][n] = 'o';
       N=nGroupsFG[n];
       for(i=1; i <= N; i++){
	  if(GroupsFG[n][i] == g){
		for(j=i; j < N; j++) GroupsFG[n][j] = GroupsFG[n][j+1];
		GroupsFG[n][j]=0;  nGroupsFG[n]--; break;
	  }
       } assert(i <= N);
     } break; 
     default: print_error("illegal input to Hpt->Change()"); break;
    }
   } else if(from == '-'){
    assert(HyperPartition[g][n] == '-'); 
    switch(to){
     case '+':	// from '-' to '+'.
       HyperPartition[g][n] = '+';
       // remove g from BG:
       N=nGroupsBG[n];
       for(i=1; i <= N; i++){
	  if(GroupsBG[n][i] == g){
		for(j=i; j < N; j++) GroupsBG[n][j]=GroupsBG[n][j+1]; 
		GroupsBG[n][j]=0; nGroupsBG[n]--; break;
	  }
       } assert(i <= N);
       // add g to FG:
       N=nGroupsFG[n];
       for(i=1; GroupsFG[n][i] < g && i <= N; i++) ;
       assert(GroupsFG[n][i] != g);
       for(j=N; j >= i; j--) GroupsFG[n][j+1]=GroupsFG[n][j];  // make room for g.
       GroupsFG[n][i]=g; nGroupsFG[n]++;
       // for(i=1; i <= nGroupsFG[n]; i++){ fprintf(stderr,"GroupFG %d = %d\n",i,GroupsFG[n][i]); }
     break;
     case 'o':		// from '-' to 'o'.
       HyperPartition[g][n] = 'o';
       // remove g from BG:
       N=nGroupsBG[n];
       for(i=1; i <= N; i++){
	  if(GroupsBG[n][i] == g){
		for(j=i; j < N; j++) GroupsBG[n][j]=GroupsBG[n][j+1]; 
		GroupsBG[n][j]=0; nGroupsBG[n]--; break;
	  }
       } assert(i <= N);
       break;
     default: print_error("illegal input to Hpt->Change()"); break;
    }
   } else {	// from == 'o'.
    assert(from == 'o'); assert(HyperPartition[g][n] == 'o'); 
    switch(to){
     case '-': 		// from 'o' to '-'.
       HyperPartition[g][n] = '-';
       N=nGroupsBG[n];
       // for(i=1; i <= N; i++){ fprintf(stderr,"GroupBG %d/%d = %d (%d)\n",i,N,GroupsBG[n][i],g); }
       for(i=1; GroupsBG[n][i] < g && i <= N; i++) ; 
       if(i <= N) assert(GroupsBG[n][i] != g);
       for(j=N; j >= i; j--) GroupsBG[n][j+1]=GroupsBG[n][j];  // make room for g.
       GroupsBG[n][i]=g; nGroupsBG[n]++;
     break;
     case '+':	// from 'o' to '+'.
       HyperPartition[g][n] = '+';
       N=nGroupsFG[n];
       for(i=1; GroupsFG[n][i] < g && i <= N; i++) ; 
       if(i <= N) assert(GroupsFG[n][i] != g);
       for(j=N; j >= i; j--) GroupsFG[n][j+1]=GroupsFG[n][j];  // make room for g.
       GroupsFG[n][i]=g; nGroupsFG[n]++;
     break;
     // break;
     default: print_error("illegal input to Hpt->Change()"); break;
    }
   }
}

Int4    hpt_typ::AddEdgeToTree(Int4 col,Int4 row, Int4 root, Int4 EdgeWt, wdg_typ Tree, wdg_typ &NewTree)
// Add an edge to the NewTree.
{
	Int4 Wt,wt,c,p,e,v;
        // fprintf(stderr,"col=%d; row=%d\n",col,row); fflush(stderr);
        Int4 cs,child = row,ps,parent = col; // 'col' corresponds to the parent node for set in row == col.
        char *tmp=this->ElmntSetName(child); if(sscanf(tmp,"Set%d",&cs) != 1) print_error("SampleHpt() graph error 1");
        tmp=this->ElmntSetName(parent); if(sscanf(tmp,"Set%d",&ps) != 1) print_error("SampleHpt() graph error 2");
        assert(cs!=ps);
#if 1	// See whether the original graph has the sets arranged oppositely...if so skip this...
	assert(NumElementarySets == NumberBPPS + 1);
	if(HyperPartition[col][row] == '+'){	// looks like original graph has this the other way around...
         	e=FirstInWdgraph(ps,Tree);     // see if already joined...
		v=TailWdgraph(e,Tree);
		if(v == cs){
			fprintf(stderr,"\nProposed child %d (Set%d) is parent of proposed parent %d (Set%d)(wt=%d); ignored.\n",
				child,cs,parent,ps,EdgeWt);
			return 0;
		}
	}
#endif
        Int4 ChildOfMRCA=FindChildOfMRCA(ps,cs,root,Tree);

	if(cs == ChildOfMRCA) fprintf(stderr,"\nChild %d (Set%d) <-- parent %d (Set%d)\n",child,cs,parent,ps); 
	else fprintf(stderr,"\nChild %d (Set%d)...Set%d <-- parent %d (Set%d)\n",child,cs,ChildOfMRCA,parent,ps); 
        // fprintf(stderr,"Child of MRCA = Set%d\n",ChildOfMRCA);

        Int4 edge=FirstInWdgraph(ChildOfMRCA,NewTree);     // see if already joined...
	if(edge == 0){ 	// No edge into ChildOfMRCA node
	    if(IsJoined(ChildOfMRCA,ps,NewTree)){	// Are the nodes joined the other way around?
		e=FirstInWdgraph(ps,NewTree);	// then find the current edge into p.
		assert(NextInWdgraph(e,NewTree) == 0);
		wt= WeightWdgraph(e,NewTree);
		if(wt < EdgeWt){        // If new edge is better than current edge; replace the edge.
			wdg_typ TmpTree=CopyGraphWithoutEdge(e,NewTree);
         	        fprintf(stderr," removing %d --> %d (wt=%d)\n",ChildOfMRCA,ps,wt); fflush(stderr);
			JoinWdgraph(ps,ChildOfMRCA,EdgeWt,TmpTree); //  
         		fprintf(stderr," joining %d --> %d (wt=%d)\n",ps,ChildOfMRCA,EdgeWt); fflush(stderr);
			NilWdgraph(NewTree); NewTree=TmpTree;
		}  else { // do nothing...
		   fprintf(stderr," %d --> %d oppositely joined (wt=%d; new wt = %d) do nothing.\n",
						ChildOfMRCA,ps,wt,EdgeWt);
		   return 0;
		}
	    } else {
		edge=JoinWdgraph(ps,ChildOfMRCA,EdgeWt,NewTree);
         	fprintf(stderr," joining %d --> %d (wt=%d)\n",ps,ChildOfMRCA,EdgeWt); 
		fflush(stderr); // return edge; // add ps --> cs.
	    }
	} else {	// ChildOfMRCA is already connected to a node; this is more complicated...
            assert(NextInWdgraph(edge,NewTree) == 0);	// Confirm there is only one 'in' edge.
            p=TailWdgraph(edge,NewTree); c=HeadWdgraph(edge,NewTree);
	    assert(c==ChildOfMRCA);
	    if(p == ps){
	        wt=WeightWdgraph(edge,NewTree);
		if(wt < EdgeWt) SetWeightWdgraph(EdgeWt,edge,NewTree);
		fprintf(stderr," %d --> %d already joined (wt=%d;new wt = %d).\n",ps,ChildOfMRCA,wt,EdgeWt);
		fflush(stderr); return 0;     // The nodes are already joined as desired, do nothing..
	    }

	    // ChildOfMRCA has another parent; see if it should be replaced.
	    Wt=WeightWdgraph(edge,NewTree);	// get the weight for this edge.
	    if(Wt < EdgeWt){		// If new edge is better than current edge; replace the edge.
		// But before replacing it, first see if this would form a cycle...
	    	if(IsJoined(ChildOfMRCA,ps,NewTree)){	// Are the nodes also connected the other way?
		  e=FirstInWdgraph(ps,NewTree); // then find the current edge into p.
         	  assert(NextInWdgraph(e,NewTree) == 0);
		  wt= WeightWdgraph(e,NewTree);
		  if(wt < EdgeWt){	// If new edge is better than current edge; replace the edge.
		    wdg_typ TmpTree=CopyGraphWithoutEdge(e,NewTree);  // remove current parent edge.
		    // JoinWdgraph(ps,ChildOfMRCA,EdgeWt,TmpTree);  // Add new edge below...
         	    fprintf(stderr," removing %d --> %d (wt=%d; new wt = %d)\n",ChildOfMRCA,ps,wt,EdgeWt); fflush(stderr);
		    NilWdgraph(NewTree); NewTree=TmpTree;
		  } else {
		    fprintf(stderr," %d --> %d oppositely joined (wt=%d; new wt = %d) do nothing.\n",ChildOfMRCA,ps,wt,EdgeWt);
		    return 0;	// don't remove current edge into ps and can't create a cycle.
		  }
		} // Now replac the parent edge...
		wdg_typ TmpTree=CopyGraphWithoutEdge(edge,NewTree);
         	fprintf(stderr," removing %d --> %d (wt=%d)\n",p,c,Wt); fflush(stderr);
		JoinWdgraph(ps,ChildOfMRCA,EdgeWt,TmpTree); 
         	fprintf(stderr," joining %d --> %d (wt=%d)\n",ps,ChildOfMRCA,EdgeWt); fflush(stderr);
		NilWdgraph(NewTree); NewTree=TmpTree;
	    } else {
		fprintf(stderr," %d --> %d already better joined (wt=%d; new wt = %d) do nothing.\n",ps,ChildOfMRCA,Wt,EdgeWt);
		fflush(stderr); return 0;   
	    }
	} 
	// PrintNewickTreeWDG(stderr,root,NewTree); PutWdgraph(stderr,NewTree); 
	RemoveCycle(cs,NewTree); // if cycle found, then remove weakest link.
        return edge;
}

Int4	hpt_typ::RemoveCycle(Int4 node,wdg_typ &G) // if cycle found, then remove weakest edge.
// Use to connect an edge from parent to direct child of MRCA along path to descendent child.
{
	Int4	v,e,n,*path,*pathWt,*pathE,i,j,worstWt,worst;

	NEW(path,WdgraphN(G)+3,Int4); NEW(pathWt,WdgraphN(G)+3,Int4); NEW(pathE,WdgraphM(G)+3,Int4); 
	set_typ S=MakeSet(WdgraphN(G)+1); ClearSet(S);

	n=1; v=path[n]=node;
	while(v != 0){	// these should end up without a parent.
	   if(MemberSet(v,S)){  // Input graph cycle found; remove weakest link...
		for(i=1; path[i] != v; i++) ;	// go to first occurance of v.
		assert(i < n && path[i]==v && path[n]==v);
		worst=0; worstWt=INT4_MAX;
		// for(j=i; j <= n; j++){	// nth node is v, same as previous InEdge; so pathE[n]=0;
		for(j=i; j < n; j++){	// nth node is v, same as previous InEdge.
		    if(worstWt > pathWt[j]){ worstWt=pathWt[j]; worst=pathE[j]; }
		}
#if 1	//DEBUG:
		fprintf(stderr," Node=%d; i=%d; path[i]=%d; n=%d; path[n]=%d\n",node,i,path[i],n,path[n]);
		for(Int4 x=1; path[x]; x++) fprintf(stderr,"%d <--",path[x]); fprintf(stderr,"\n");
		fprintf(stderr," Removing edge: Set%d --(%d)--> Set%d (wt=%d) to eliminate cycle\n",
				TailWdgraph(worst,G),worst,HeadWdgraph(worst,G),worstWt);
#endif
		wdg_typ TmpG=CopyGraphWithoutEdge(worst,G); NilWdgraph(G); G = TmpG;
		{ NilSet(S); free(path); free(pathWt); free(pathE); }
		return 0;
	   } AddSet(v,S); e=FirstInWdgraph(v,G);
	   if(e==0){	// i.e., reach a node without a parent; implies no cycle.
		fprintf(stderr," RemoveCycle( ) node %d lacks a parent; couldn't find a cycle\n",v);
		break; // free memory and return;
	   } assert(NextInWdgraph(e,G) == 0);  // should only be one 'in' edge per node.
	   pathWt[n]=WeightWdgraph(e,G); pathE[n]=e; n++; v=path[n]=TailWdgraph(e,G); 
	} NilSet(S); free(path); free(pathWt); free(pathE);
	return n;
}

#if 1	// Move this code to new tre_typ later.
wdg_typ	hpt_typ::CopyGraphWithoutEdge(Int4 E,wdg_typ G_in)
// Make a copy the graph G, but without the edge e; set e = 0 to copy entire graph.
{
	Int4    wt,p,c,e; 
	wdg_typ G=MkWdgraph(WdgraphN(G_in),WdgraphM(G_in));
	for(e=1; e <= mWdgraph(G_in); e++){ // // arrow points from parent (tail) to child (head).
		if(e == E) continue;	// skip this edge.
		p=TailWdgraph(e,G_in); c=HeadWdgraph(e,G_in); wt= WeightWdgraph(e,G_in);
		JoinWdgraph(p,c,wt,G);
	} return G;
}

Int4    hpt_typ::MergeOldTreeIntoNew(Int4 Root,wdg_typ OT,wdg_typ &Tree)
// merge the OldTree (OT) in with Tree, but don't connect children in OT that are connected in Tree.
{
	Int4    p,c,v,u,i,j,e,f,n,wt; 
#if 1	// DEBUG.
	PutWdgraph(stderr,OT); PrintNewickTreeWDG(stderr,Root,OT); PutWdgraph(stderr,Tree); 
#endif
	for(e=1; e <= mWdgraph(OT); e++){ // // arrow points from parent (tail) to child (head).
		p=TailWdgraph(e,OT); c=HeadWdgraph(e,OT); 
		f=FirstInWdgraph(c,Tree);
		if(f == 0){	// no edge into c in new tree.
			wt=WeightWdgraph(e,OT);
			JoinWdgraph(p,c,wt,Tree);
		} else if(NextInWdgraph(f,Tree)!=0){
			p=TailWdgraph(f,Tree); c=HeadWdgraph(f,Tree);
			fprintf(stderr,"edge %d: %d --> %d\n",f,p,c);
			print_error("MergeOldTreeWithNew( ) input graph is not a tree.");
		}  
	}
	if((n=RootUnConnectedTree(Root,Tree)) > 0){
		fprintf(stderr,"Rooted %d unconnected subtree(s).\n",n);
	}
}

Int4	hpt_typ::RootUnConnectedTree(Int4 Root, wdg_typ &Tree)
{
	Int4	v,e,n,*path=0,oldest,N=0,restarts=0,M=mWdgraph(Tree);
	// 1. go through leaf nodes to find one that does not lead back to the Root node.
	set_typ S=MakeSet(WdgraphN(Tree)+1);
	for(e=1; e <= mWdgraph(Tree); e++){
	   v=HeadWdgraph(e,Tree);	// ---(e)---> v is a child node.
	   assert(v != 0);
	   if(FirstOutWdgraph(v,Tree) == 0){	// then v is a leaf node...
	     n=PathToRoot(v,0,path,S,Tree); // using 0 as Root find the path to the oldest node.
	     if(n == 0){
		fprintf(stderr,"RootUnConnectedTree( ): cycle found & removed, starting over...\n");
		fprintf(stderr," Root=%d: edge =  %d --(%d)--> %d\n",Root,TailWdgraph(e,Tree),e,v);
		restarts++; 
		if(restarts > M){ 
			fprintf(stderr,"Root=%d: edge =  %d --(%d)--> %d\n",Root,TailWdgraph(e,Tree),e,v);
			PrintNewickTreeWDG(stderr,Root,Tree); PutWdgraph(stderr,Tree); 
			for(Int4 j=1; path[j]; j++) fprintf(stderr,"%d <--",path[j]); fprintf(stderr,"\n");
			assert(restarts <= M); 
		} RemoveCycle(v,Tree); free(path); e=0; // continue; 
	     } else {
	       assert(n > 2 && path[n]==0); // path must be at least of length 3 and end in 0 [v,p,0] 
	       oldest=path[n-1];  free(path);
	       if(oldest == Root) continue;
	       // 2. If one is found, connect the 'oldest' node to the Root.
	       JoinWdgraph(Root,oldest,1,Tree);	// make Root the parent of oldest.
	       N++;	// count number of unconnected subtrees that are rooted.
	     }
	   }	// else v is an internal or root node...continue;
	} NilSet(S); return N;
}

//************************* This will replace other FD-table operations...
Int4	hpt_typ::PathToRoot(Int4 node, Int4 Root, Int4 *&path, set_typ S, wdg_typ Tree)
{
	Int4	v,u,e,n;

	NEW(path,WdgraphN(Tree)+3,Int4); ClearSet(S); 
	n=1; v=path[n]=node;
	while(v != Root){
	   assert(v!=0);
	   if(MemberSet(v,S)) return 0; else AddSet(v,S);
	   e=FirstInWdgraph(v,Tree);	// u --(e)--> v
	   assert(NextInWdgraph(e,Tree) == 0);  // only one in per node; it's a tree!
	   if(e==0) u = 0; else u=TailWdgraph(e,Tree); // if e == 0 then u = 0.
	   assert(u != v);	// this should not happen as FirstOutWdgraph(v,Tree) == 0 above.
	   n++; v=path[n]=u;
	} return n;
}

Int4	hpt_typ::FindChildOfMRCA(Int4 Parent,Int4 child,Int4 Root,wdg_typ Tree)
// Find the direct child of the most recent common ancestor along path to descendent child.
// Use to connect an edge from parent to direct child of MRCA along path to descendent child.
{
	Int4    c,v,u,i,j,*pathC,*pathA,stepsC,stepsA;
	set_typ S=MakeSet(WdgraphN(Tree)+1); 
	stepsC = PathToRoot(child,Root,pathC,S,Tree); 
 	if(stepsC==0) print_error("hpt_typ::PathToRoot( ): input graph not a tree (cycle found).");
	stepsA = PathToRoot(Parent,Root,pathA,S,Tree); NilSet(S);
 	if(stepsA==0) print_error("hpt_typ::PathToRoot( ): input graph not a tree (cycle found).");
	for(i=1; i <= stepsC; i++){
	   v = pathC[i];
	   for(j=1; j <= stepsA; j++){
		u = pathA[j];
		if(v == u){
		   assert(i > 1);
		   c=pathC[i-1]; free(pathC); free(pathA); return c;  
		}
	   }
	} print_error("FindChildOfMRCA( ) input graph is not a tree.");
}
#endif

void	hpt_typ::PutInsertGrp(FILE *fp,Int4 gi)
// insert a dummy group directly before column gi
// if NumElementarySets - NumberBPPS  == 1 then add a row as well.
{
	Int4	c,s,S,n,N,i,j,g;
	BooLean	InsertRow=FALSE;

	assert(gi >= 0 && gi <= NumberBPPS); 
	if(NumElementarySets - NumberBPPS  == 1) InsertRow=TRUE;
	if(Mode == 'I') fprintf(fp,"\nHyperpartition:\n"); else fprintf(fp,"\nHyperParTition:\n");
	if(gi == 0){ fprintf(fp,"!"); }
	for(c=0,g=1; c < NumberBPPS; c++,g++){
	    fprintf(fp,"%c",OutputCols[c]);
	    if(g == gi){ fprintf(fp,"!"); }
	} fprintf(fp,"\n");
	for(s=S=1; s <= NumElementarySets; S++,s++){
	  if(gi == 0){ fprintf(fp,"+"); }
	  for(n=1; n <= NumberBPPS; n++){
	    fprintf(fp,"%c",HyperPartition[s][n]);
	    if(n == gi){	// new column...
		if(s == NumElementarySets) fprintf(fp,"o"); 
		// else if(s == gi) fprintf(fp,"+"); 
		else fprintf(fp,"-"); 
	    }
	  }
	  if(SetType[s] == '=') fprintf(fp," %d.%s=%d.",S,NameElementarySet[s],NumRandom);
	  else fprintf(fp," %d.%s%c",S,NameElementarySet[s],SetType[s]);
	  for(Int4 a=0; a < ArgC[s]; a++) fprintf(fp," %s",ArgV[s][a]); 
	  fprintf(fp,"\n");
	  if(InsertRow && s == gi){
		S++;
		fprintf(fp,"+");
		for(n=1; n <= NumberBPPS; n++){
		    if(n == gi) fprintf(fp,"+");
		    else fprintf(fp,"o"); 
		}
		fprintf(fp," %d.dummy%d?\n",S,dummyID++); 
	  }
	}
#if 0
	for(i=0,n=1; n <= NumberBPPS; n++){ if(nArgmnt[n] > 0) i++; }
	if(i > 0){
	  fprintf(fp,"\nSettings:\n"); N=1;
	  if(gi == 0){ fprintf(fp,"%d.Dummy \n",N); N++; }
	  for(n=1; n <= NumberBPPS; n++,N++){
	    if(nArgmnt[n] > 0) fprintf(fp,"%d.%s\n",N,ArgmntStr[n]);
	    if(n == gi){ N++; fprintf(fp,"%d.Dummy \n",N); }
	  } fprintf(fp,"\n\n");
	}
#else
	fprintf(fp,"\nSettings:\n"); N=1;
	if(gi == 0){ fprintf(fp,"%d.Dummy \n",N); N++; }
	for(n=1; n <= NumberBPPS; n++,N++){
	    fprintf(fp,"%d.%s\n",N,ArgmntStr[n]);
	    if(n == gi){ N++; fprintf(fp,"%d.Dummy \n",N); }
	} fprintf(fp,"\n\n");
#endif
}

void	hpt_typ::PutSwappedRows(FILE *fp,Int4 s1, Int4 s2)
{
	Int4	c,s,S,n,i,j,g,row,col;
	BooLean	put_settings=TRUE;

	assert(s1 > 0 && s1 < NumElementarySets); 
	assert(s2 > 0 && s2 < NumElementarySets);	// Don't allow Reject set to be swapped. 
	if(s1 == s2){ Put(fp); return; }
	if(Mode == 'I') fprintf(fp,"\nHyperpartition:\n"); else fprintf(fp,"\nHyperParTition:\n");
	for(n=0; n < NumberBPPS; n++) if(OutputCols[n] != '*') fprintf(fp,"%c",OutputCols[n]);
        fprintf(fp,"\n");

        for(row=0,s=1; s <= NumElementarySets; s++){
           if(SetType[s] == '*') continue; else row++;
	   if(s == s1) S=s2; else if(s == s2) S=s1; else S=s;
           for(n=1; n <= NumberBPPS; n++){
              if(OutputCols[n-1] == '*') continue;
              if(HyperPartition[S][n] != 0){ fprintf(fp,"%c",HyperPartition[S][n]); }
              else { fprintf(fp,"o"); }
           }
           if(SetType[S] == '=') fprintf(fp," %d.%s=%d.",row,NameElementarySet[S],NumRandom);
           else fprintf(fp," %d.%s%c",row,NameElementarySet[S],SetType[S]);
           for(Int4 a=0; a < ArgC[S]; a++) fprintf(fp," %s",ArgV[S][a]);
           fprintf(fp,"\n");
        } fprintf(fp,"\n\n");

        if(put_settings && i > 0){
          fprintf(fp,"\nSettings:\n");
          for(col=0,n=1; n <= NumberBPPS; n++){
           if(OutputCols[n-1] == '*') continue; else col++;
            if(nArgmnt[n] > 0) fprintf(fp,"%d.%s\n",col,ArgmntStr[n]);
            else fprintf(fp,"%d.Column%d\n",col,col);
          } fprintf(fp,"\n\n");
        }
}


void	hpt_typ::PutSwapped(FILE *fp,Int4 g1, Int4 g2)
{
	Int4	c,s,n,i,j,g;

	assert(g1 > 0 && g1 <= NumberBPPS); assert(g2 > 0 && g2 <= NumberBPPS); 
	if(g1 == g2){ Put(fp); return; }
	if(Mode == 'I') fprintf(fp,"\nHyperpartition:\n"); else fprintf(fp,"\nHyperParTition:\n");
	for(c=0,g=1; c < NumberBPPS; c++,g++){
	    if(g == g1) fprintf(fp,"%c",OutputCols[g2-1]);
	    else if(g == g2) fprintf(fp,"%c",OutputCols[g1-1]);
	    else fprintf(fp,"%c",OutputCols[c]);
	} fprintf(fp,"\n");
	for(s=1; s <= NumElementarySets; s++){
	  for(n=1; n <= NumberBPPS; n++){
	    if(n == g1) fprintf(fp,"%c",HyperPartition[s][g2]);
	    else if(n == g2) fprintf(fp,"%c",HyperPartition[s][g1]);
	    else fprintf(fp,"%c",HyperPartition[s][n]);
	  }
	  if(SetType[s] == '=') fprintf(fp," %d.%s=%d.",s,NameElementarySet[s],NumRandom);
	  else fprintf(fp," %d.%s%c",s,NameElementarySet[s],SetType[s]);
	  for(Int4 a=0; a < ArgC[s]; a++) fprintf(fp," %s",ArgV[s][a]); 
	  fprintf(fp,"\n");
	}
#if 0
	for(i=0,n=1; n <= NumberBPPS; n++){ if(nArgmnt[n] > 0) i++; }
	if(i > 0){
	  fprintf(fp,"\nSettings:\n");
	  for(n=1; n <= NumberBPPS; n++){
	    if(n == g1){ 
	    	if(nArgmnt[g2] > 0) fprintf(fp,"%d.%s\n",n,ArgmntStr[g2]);
	    } else if(n == g2){
	    	if(nArgmnt[g1] > 0) fprintf(fp,"%d.%s\n",n,ArgmntStr[g1]);
	    } else if(nArgmnt[n] > 0) fprintf(fp,"%d.%s\n",n,ArgmntStr[n]);
	  } fprintf(fp,"\n\n");
	}
#else
	fprintf(fp,"\nSettings:\n");
	for(n=1; n <= NumberBPPS; n++){
	    if(n == g1){ 
	    	if(nArgmnt[g2] > 0) fprintf(fp,"%d.%s\n",n,ArgmntStr[g2]);
	    } else if(n == g2){
	    	if(nArgmnt[g1] > 0) fprintf(fp,"%d.%s\n",n,ArgmntStr[g1]);
	    } else if(nArgmnt[n] > 0) fprintf(fp,"%d.%s\n",n,ArgmntStr[n]);
	} fprintf(fp,"\n\n");
#endif
}

void	hpt_typ::PutWithout(FILE *fp,Int4 col,Int4 row)
{
       Int4    c,r,s,m,n,i,j,g;

        assert(col > 0 && col <= NumberBPPS); 
        if(Mode == 'I') fprintf(fp,"\nHyperpartition:\n"); else fprintf(fp,"\nHyperParTition:\n");
        for(c=0; c < NumberBPPS; c++){
            if(c != col) fprintf(fp,"%c",OutputCols[c]);
        } fprintf(fp,"\n");
        for(r=0,s=1; s <= NumElementarySets; s++){
	  if(s== row) continue; else r++;
          for(n=1; n <= NumberBPPS; n++) if(n != col) fprintf(fp,"%c",HyperPartition[s][n]);
          if(SetType[s] == '=') fprintf(fp," %d.%s=%d.",r,NameElementarySet[s],NumRandom);
          else fprintf(fp," %d.%s%c",r,NameElementarySet[s],SetType[s]);
          for(Int4 a=0; a < ArgC[s]; a++) fprintf(fp," %s",ArgV[s][a]);
          fprintf(fp,"\n");
        }
        fprintf(fp,"\nSettings:\n");
        for(m=1,n=1; n <= NumberBPPS; n++){
            if(n != col){
		if(nArgmnt[n] > 0) fprintf(fp,"%d.%s\n",m,ArgmntStr[n]);
		else fprintf(fp,"%d.Column%d\n",m,m); m++; 
	    }
        } fprintf(fp,"\n\n");
}

void	hpt_typ::PutSubTree(FILE *fp, Int4 target)
// Add a leaf to node x.
{
	Int4	s,i,j,n,S,N,col,parent,*Parent;
	char    state,**HP=HyperPartition;

	if(this->IsTree(Parent) == FALSE){
                if(Parent) free(Parent);
                print_error("PutSubTree() input error: Hpt not in tree format ");
        }
	if(target < 1 || target > NumberBPPS) print_error("PutSubTree() error: 'target' node out of range.");
	set_typ	*set,Set=MkSubTreeSet(target);
	AddSet(NumElementarySets,Set); // add random set...
	NEW(set,NumElementarySets+3,set_typ);
	for(s=1; s < NumElementarySets; s++){
		if(MemberSet(s,Set)){ set[s]=MkSubTreeSet(s); }
	}

	if(Mode == 'I') fprintf(fp,"\nHyperpartition:\n"); else fprintf(fp,"\nHyperParTition:\n");
	for(n=1; n <= NumberBPPS; n++){
	     if(MemberSet(n,Set)){ fprintf(fp,"%c",OutputCols[n-1]); }
	} fprintf(fp,"\n");
	for(S=0,s=1; s <= NumElementarySets; s++){	// == Row.
	   if(!MemberSet(s,Set)) continue; else S++;
           for(n=1; n <= NumberBPPS; n++){		// == column.
	      if(!MemberSet(n,Set)) continue;
	      parent=Parent[n];
	      if(n==target && SetType[s] == '=') fprintf(fp,"-"); 
	      else if(n==target || n==s) fprintf(fp,"+"); // all s in target and its own foreground.
	      else if(set[parent]==0) fprintf(fp,"o");	// this should only for target and random.
	      else if(MemberSet(s,set[n])){		// s is in foreground (i.e., subtree) of n.
			fprintf(fp,"+");
	      } else if(MemberSet(s,set[parent])){	// s is in background of n.
			fprintf(fp,"-");
	      } else fprintf(fp,"o");			// s is in omitted set of n.
	      // else fprintf(fp,"%c",HP[s][n]);
           } 
	   if(SetType[s] == '=') fprintf(fp," %d.%s=%d.",S,NameElementarySet[s],NumRandom);
	   else fprintf(fp," %d.%s%c",S,NameElementarySet[s],SetType[s]);
	   for(Int4 a=0; a < ArgC[s]; a++) fprintf(fp," %s",ArgV[s][a]); 
	   fprintf(fp,"\n");
	} N=CardSet(Set);
	if(N > 5){
	 for(n=1; n <= N; n++){
	   if(n ==1) fprintf(fp,"#"); else if(n%5 ==0) fprintf(fp,"|"); else fprintf(fp," ");
	 } fprintf(fp,"\n");
	 for(n=1; n <= N; n++){
		if(n==1) fprintf(fp,"#"); else if(n%5 ==0) fprintf(fp,"%4d ",n);
	 } fprintf(fp,"\n");
	} else fprintf(fp,"\n"); 

	// 4. Print out the Settings.
	fprintf(fp,"\nSettings:\n");
	for(S=0,n=1; n <= NumberBPPS; n++){
	    if(!MemberSet(n,Set)) continue; else S++;
	    fprintf(fp,"%d.%s\n",S,ArgmntStr[n]);
	} fprintf(fp,"\n\n");

	for(s=1; s <= NumElementarySets; s++){ if(set[s]) NilSet(set[s]); }
	NilSet(Set); free(set);
	if(Parent) free(Parent);
}


