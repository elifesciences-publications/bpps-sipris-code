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

BooLean	mcs_typ::MvUpDown(Int4 grandpa,Int4 parent,Int4 child, set_typ subtree, char State)
// WARNING: this routine assumes direct correspondence between column and rows!!
// So that the diagonal equates rows (sets) with columns (FD-categories)!!
/**********************************************************************
This operation will not change node rows and columns, but could put a
parent '+' symbol after the node diagonal for a given row.

    *     Up      *                     *   Down      *         * = grandpa
    |             |\		        |\            |         X = toggled parent
    X     -->     X x                   X x   -->     X         x = child
   /|\           /|                    /|            /|\
  o o x         o o                   o o           o o x

 ************************** Down (State == '+') ****************************
HyperParTition:			// Sample Set5 into Set15?
+--oo--ooo----- 1.Set1?		+--oo--ooo----- 1.Set1?   (original parent)
++-oo--ooo----- 2.Set12!
   :     :         :
+--oo+-ooo----- 6.Set8!
+--oo-+-------- 7.Set15?        +--oo-+-------- 7.Set15?  (NEW, INTERMEDIATE PARENT)
+--oo-++------- 8.Set16!
+--oo-+-+------ 9.Set14!
+--oo-+--+----- 10.Set7!
+--oo--ooo+---- 11.Set6!
+--oo--ooo-+--- 12.Set5!  FROM: +--oo--ooo-+--- 12.Set5!  (original row)
                                   :        :       :
                            TO: +--oo-+----+--- 12.Set5!  (NEW ROW)
				      ^-new (intermediate) parent column
+--oo--ooo--+-- 13.Set4! 	(needs to be the same as parent except for diagonal).
   :     :          :

 ************************** Up (State == '-') ******************************
+--oo--ooo----- 1.Set1?		+--oo--ooo----- 1.Set1?   (NEW PARENT)
   :     :          :
+--oo-+-------- 7.Set15?        +--oo-+-------- 7.Set15?  (original intermediate parent)
   :     :          :
+--oo--ooo-+--- 12.Set5!  FROM: +--oo-+----+--- 12.Set5!  (original row)
                                   V     V         V
                            TO: +--oo--ooo-+--- 12.Set5!  (NEW ROW)

Bottom Line: all rows needs to be the same as parent except for parent and node positions.
 **********************************************************************/
// 3. Change that cell and recompute LPR.
// 4. If this improves the LPR then accept the change, else revert back... 
//************ Testing this for cdhBPPS routine: adding an internal node...
// WARNING: assumes a direct correspondence between columns and rows!!
{
	char	state;
	double	olpr,nlpr,d;
	Int4	row_neg,row_pos,row_omit,col_neg,col_pos,col_omit;
	Int4	n,gp,p,ch,r,c,x,sq,*Parent=0,num_sq=NumSeqsCMSA(MainCMA);
	BooLean	IsChanged=FALSE;

	assert(Hpt->NumSets() == Hpt->NumBPPS() +1); // need correspondence between column & rows!
	gp = grandpa; p = parent; ch=child;
	if(parent <= 1) return FALSE;	// don't look at root node.
	if(child < 2 || child >= Hpt->NumSets()) return FALSE;
	if(Hpt->Cell(ch,p) == State) return FALSE;	// doesn't change state...
	assert(State == '+' || State == '-'); // only toggle between parent and grandparent.

	hpt_typ	*oHpt=Hpt->Copy();
	char    **HP=oHpt->RtnHyperPartition();
	assert(Hpt->IsScrambledTree(Parent));
// for(c=1; c <= Hpt->NumBPPS(); c++){ fprintf(stderr,"Parent[%d]=%d\n",c,Parent[c]); }
	// if(Hpt->IsScrambledTree(Parent) == FALSE){ if(Parent) free(Parent); return FALSE; }
	// if(Hpt->IsTree(Parent) == FALSE){ if(Parent) free(Parent); return FALSE; }
	// if(Parent[ch] != Parent[p]) continue;  // look only at siblings of n.
	state=Hpt->Cell(ch,p);	
	assert(state == '-' || state == '+');	// must be the case if parent in a tree.
	// if(Hpt->TypeOfSet(ch) == '?') assert(subtree != 0); // then bring along all child nodes as well.
	// else assert(subtree == 0);
	if(State == '+'){	// then move child __DOWN__ from grandpa to parent.
	  if(Hpt->TypeOfSet(p) == '!'){ Hpt->ChangeTypeOfSet(p,'?'); }
	  for(c=1; c <= Hpt->NumBPPS(); c++){
	    if(c == gp) continue;	// this is the current parent/new grandparent.
	    if(c == parent){	// this is the new parent...
	        for(r=1; r <= Hpt->NumBPPS(); r++){	// over all rows == sets...
	           if(r==c) continue;	// diagonal does not change...
		   if(MemberSet(r,subtree)){
			assert(HP[r][c] == '-'); 
			Hpt->Change('-','+',r,c);
		   	for(sq=1; sq <= num_sq; sq++){
		     	   if(MemberSet(sq,GrpSet[r])){ che[c]->SetPartition('+',sq); } 
		   	}
		   }
		}
	    } else if(c == child){
		for(r=1; r <= Hpt->NumBPPS(); r++){	// over all rows == sets...
	           if(r==c) continue;	// diagonal does not change...
		   if(MemberSet(r,subtree)) continue; // child's descendants don't change..
#if 0
		   if(HP[r][gp] == '+' && HP[r][p] == '+'){ // r is a current sibling...
			assert(HP[r][c] == '-'); 	// so that c is in the background.
		   } else 
#endif
		   if(HP[r][gp] == '+' && HP[r][p] == '-'){ // r is a current sibling...
			assert(HP[r][c] == '-'); 	// so that c is in the background.
			Hpt->Change('-','o',r,c);
		        for(sq=1; sq <= num_sq; sq++){
			   if(MemberSet(sq,GrpSet[r])){ che[c]->SetPartition('o',sq); }
                        }
		   }
		}
	    } else if(Parent[c] == parent){
	        if(!MemberSet(c,subtree)){ // New sibling of child.
		  for(r=1; r <= Hpt->NumBPPS(); r++){	// then need to remove these from BG.
		    if(r==c) continue;
		    if(MemberSet(r,subtree)){
			assert(HP[r][c] == 'o');
			// fprintf(stderr," Moving %d down from %d to %d (%d,%d)\n",ch,gp,p,r,c);
			Hpt->Change('o','-',r,c);
		        for(sq=1; sq <= num_sq; sq++){
			   if(MemberSet(sq,GrpSet[r])){ che[c]->SetPartition('-',sq); }
                        }
		    }
		  }
		}
	    } // else do nothing
	   }
	} else if(State == '-'){	// then move child __UP__ from parent to grandpa. 
	  // fprintf(stderr,"Now moving %d up from %d to %d\n",c,p,gp);
	  if(Hpt->TypeOfSet(p) == '?'){
	        for(n=0,r=1; r < Hpt->NumSets(); r++){ if(Parent[r] == p) n++; }
		if(n == 1) Hpt->ChangeTypeOfSet(p,'!');  // if only one child, then make a leaf
	  }
	  for(c=1; c <= Hpt->NumBPPS(); c++){
	    if(c == gp) continue;
	    if(c == parent){		// current parent...
		// fprintf(stderr," ... %d == %d\n",c,p);
		for(r=1; r <= Hpt->NumBPPS(); r++){	// over all rows == sets...
	           if(r==c) continue;	// diagonal does not change...
		   if(MemberSet(r,subtree)){
			assert(HP[r][c] == '+'); 
			// fprintf(stderr," ... moving %d up from %d to %d\n",ch,p,gp);
			Hpt->Change('+','-',r,c);
		        for(sq=1; sq <= num_sq; sq++){
			   if(MemberSet(sq,GrpSet[r])){ che[c]->SetPartition('-',sq); }
                        }
		   }
		}
	    } else if(c == child){	// child column
		for(r=1; r <= Hpt->NumBPPS(); r++){	// over all rows == sets...
	           if(r==c) continue;	// diagonal does not change...
		   if(MemberSet(r,subtree)) continue; // child & descendants don't change..
		   if(HP[r][gp] == '+' && HP[r][p] == '-'){
			assert(HP[r][c] == 'o'); 
			Hpt->Change('o','-',r,c);
		        for(sq=1; sq <= num_sq; sq++){
			   if(MemberSet(sq,GrpSet[r])){ che[c]->SetPartition('-',sq); }
                        }
		   }
		}
	    } else if(!MemberSet(c,subtree) && Parent[c] == parent){
		for(r=1; r <= Hpt->NumBPPS(); r++){	// over all rows == sets...
		    if(r==c) continue;
		    if(MemberSet(r,subtree)){
			assert(HP[r][c] == '-');
			Hpt->Change('-','o',r,c);
		        for(sq=1; sq <= num_sq; sq++){
			   if(MemberSet(sq,GrpSet[r])){ che[c]->SetPartition('o',sq); }
                        }
		    }
		}
	    } // else do nothing...
	  }
	} else print_error("this should not happen");
	ReSetRelations(); free(Parent); delete oHpt; DidRestoreBest=FALSE; // current state may not be best...
	return TRUE;
}

