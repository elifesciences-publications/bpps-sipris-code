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

Int4    hpt_typ::ColumnsInRow(Int4 row,Int4 &col_neg,Int4 &col_pos,Int4 &col_omit)
{
         col_neg=0,col_pos=0,col_omit=0;
         for(Int4 c=1; c <= NumberBPPS; c++){
            switch (Cell(row,c)){
             case '+': col_pos++; break;
             case '-': col_neg++; break;
             case 'o': col_omit++; break;
             default: print_error("ColumnsInRow(): this should not happen"); break;
            }
         } return col_neg+col_pos+col_omit;
}

Int4    hpt_typ::RowsInColumn(Int4 col,Int4 &row_neg,Int4 &row_pos,Int4 &row_omit)
// Don't count Reject row...
{
    row_neg=0,row_pos=0,row_omit=0;
    for(Int4 r=1; r < NumElementarySets; r++){
      switch (Cell(r,col)){
         case '+': row_pos++; break;
         case '-': row_neg++; break;
         case 'o': row_omit++; break;
         default: print_error("RowsInColumn(): this should not happen"); break;
      }
    } return row_neg+row_pos+row_omit;
}

Int4    hpt_typ::SameSet(Int4 i, hpt_typ *hptO)
// return the set in hptO corresponding to set i in this.
{
    Int4    j;
    assert(this->NumElementarySets == hptO->NumElementarySets);
    assert(i > 0 && i < NumElementarySets); // don't look at Random set.
    char *name=this->ElmntSetName(i);
    for(j=1; j <= NumElementarySets; j++){
    char *nameO=hptO->ElmntSetName(j);
    if(strcmp(name,nameO) == 0) return j;
    } return 0;
}

Int4    hpt_typ::PutSettings(FILE *fp,Int4 index,Int4 c)
{
    assert(c > 0 && c <= NumberBPPS);
    if(nArgmnt[c] > 0) fprintf(fp,"%d.%s\n",index,ArgmntStr[c]);
    else fprintf(fp,"%d.Column%d\n",index,index);
}

#if 0	// OLD VERSION
void	hpt_typ::Put(FILE *fp, BooLean put_settings)
{
	Int4	s,i,j,n;
	if(Mode == 'I') fprintf(fp,"\nHyperpartition:\n"); else fprintf(fp,"\nHyperParTition:\n");
	fprintf(fp,"%s\n",OutputCols);
	for(s=1; s <= NumElementarySets; s++){
	   // HyperPartition[s][NumberBPPS +1] = 0;
	   if(SetType[s] == '='){
	   	fprintf(fp,"%s %d.%s=%d.",HyperPartition[s]+1,s,NameElementarySet[s],NumRandom);
	   	// fprintf(fp,"%s %s=%d.\n",HyperPartition[s]+1,NameElementarySet[s],NumRandom);
	   } // else fprintf(fp,"%s %s%c\n",HyperPartition[s]+1,NameElementarySet[s],SetType[s]);
	     else fprintf(fp,"%s %d.%s%c",HyperPartition[s]+1,s,NameElementarySet[s],SetType[s]);
	   for(Int4 a=0; a < ArgC[s]; a++) fprintf(fp," %s",ArgV[s][a]); 
	   fprintf(fp,"\n");
	}
	for(i=0,n=1; n <= NumberBPPS; n++){
		if(n==1) fprintf(fp,"#"); else if(n%5 ==0) fprintf(fp,"|"); else fprintf(fp," ");
		if(nArgmnt[n] > 0) i++; 
	} fprintf(fp,"\n");
	for(n=1; n <= NumberBPPS; n++){
		if(n==1) fprintf(fp,"#"); else if(n%5 ==0) fprintf(fp,"%4d ",n);
	} fprintf(fp,"\n");
	if(put_settings && i > 0){
	  fprintf(fp,"\nSettings:\n");
	  for(n=1; n <= NumberBPPS; n++){
#if 0
	    fprintf(fp,"%d.%s",n,GroupName[n]);
	    for(i=0; i < nArgmnt[n]; i++){
	        fprintf(fp," %s",Argmntv[n][i]);
	    } fprintf(fp,"\n");
#else
	    if(nArgmnt[n] > 0) fprintf(fp,"%d.%s\n",n,ArgmntStr[n]);
	    else fprintf(fp,"%d.Column%d\n",n,n);
#endif
	  } fprintf(fp,"\n\n");
	}
}
#endif


void	hpt_typ::Put(FILE *fp, BooLean put_settings,BooLean ignore, BooLean putArg)
//************** New version to deal with deleted rows and columns. ****************
{
	Int4	s,i,j,n,row,col;
#if 1
	Int4	*children=0;
	if(ignore){	// then find childless (internal) nodes so can be changed to leaves.
	  // Int4 *p; assert(IsTree(p)); free(p);
	  assert(NumberBPPS == NumElementarySets-1);
	  NEW(children,NumElementarySets+3,Int4);
	  for(col=1; col<= NumberBPPS; col++){
	     assert(HyperPartition[col][col] == '+');
	     for(n=0,row=1; row < NumElementarySets; row++){
		if(col == row) continue;
		if(SetType[row] != '*' && HyperPartition[row][col] == '+'){ children[col]++; }
	     }
	  }
	}
#endif
	if(Mode == 'I'){
	    if(FileName) fprintf(fp,"\nHyperpartition(%s):\n",FileName); 
	    else fprintf(fp,"\nHyperpartition:\n"); 
	} else {
	    if(FileName) fprintf(fp,"\nHyperParTition(%s):\n",FileName); 
	    else fprintf(fp,"\nHyperParTition:\n");
	}
	for(n=0; n < NumberBPPS; n++) if(!ignore || OutputCols[n] != '*') fprintf(fp,"%c",OutputCols[n]);
	fprintf(fp,"\n");
	for(row=0,s=1; s <= NumElementarySets; s++){
	   if(ignore && SetType[s] == '*') continue; else row++;
           for(n=1; n <= NumberBPPS; n++){
	      if(ignore && OutputCols[n-1] == '*') continue;
              if(HyperPartition[s][n] != 0){ fprintf(fp,"%c",HyperPartition[s][n]); }
              else { fprintf(fp,"o"); }
           } 
	   if(SetType[s] == '=') fprintf(fp," %d.%s=%d.",row,NameElementarySet[s],NumRandom);
	   else {
		if(ignore && SetType[s] == '?' && children[s] < 1){
		    fprintf(fp," %d.%s%c",row,NameElementarySet[s],'!');
		} else fprintf(fp," %d.%s%c",row,NameElementarySet[s],SetType[s]);
	   }
	   for(Int4 a=0; a < ArgC[s]; a++) fprintf(fp," %s",ArgV[s][a]); 
	   fprintf(fp,"\n");
	}
	if(NumberBPPS > 5){
	 for(col=i=0,n=1; n <= NumberBPPS; n++){
	   if(ignore && OutputCols[n-1] == '*') continue; else col++;
	   if(col==1) fprintf(fp,"#"); else if(col%5 ==0) fprintf(fp,"|"); else fprintf(fp," ");
	   if(nArgmnt[n] > 0) i++; 
	 }
	 fprintf(fp,"\n");
	 for(n=1; n <= col; n++){
		if(n==1) fprintf(fp,"#"); else if(n%5 ==0) fprintf(fp,"%4d ",n);
	 } fprintf(fp,"\n");
	} else { fprintf(fp,"\n\n");  i=0; }
	// if(put_settings && i > 0)
	if(put_settings)
	{
	  fprintf(fp,"\nSettings:\n");
	  for(col=0,n=1; n <= NumberBPPS; n++){
	   if(ignore && OutputCols[n-1] == '*') continue; else col++;
	   if(putArg){
	    fprintf(fp,"%d.%s",n,GroupName[n]);
	    for(i=0; i < nArgmnt[n]; i++){
	        fprintf(fp," %s",Argmntv[n][i]);
	    } fprintf(fp,"\n");
	   } else {
	    if(nArgmnt[n] > 0) fprintf(fp,"%d.%s\n",col,ArgmntStr[n]);
	    else fprintf(fp,"%d.Column%d\n",col,col);
	   }
	  } fprintf(fp,"\n\n");
	} if(children) free(children);
}

void    hpt_typ::PutHyperPartition(FILE *fp)
{
     Int4       n,g;
     fprintf(fp,"\n ");
     for(n=1; n<= NumberBPPS; n++){ fprintf(fp,"="); }
     fprintf(fp," HyperPartition: ");
     for(n=1; n<= NumberBPPS; n++){ fprintf(fp,"="); } fprintf(fp,"\n");
     fprintf(fp,"      _Category_\nSet: ");
     for(n=1; n<= NumberBPPS; n++){ fprintf(fp,"%2d ",n); } fprintf(fp,"\n");
     for(g=1; g<= NumElementarySets; g++){
        fprintf(fp,"%3d: ",g);
        Int4 NumZero=0;
        for(n=1; n<= NumberBPPS; n++){
           if(HyperPartition[g][n] != 0){ fprintf(fp," %c ",HyperPartition[g][n]); }
           else { NumZero++; fprintf(fp," o "); }
        } 
	// fprintf(stdout," %s (%d)\n",NameSubGroupCMA[g],CardSet(GroupSet[g]));
	fprintf(fp," %s\n",NameElementarySet[g]);
        if(NumZero == NumberBPPS) print_error("Fatal: seq group absent from all partitions");
     } fprintf(fp,"\n");
}

void	hpt_typ::Free()
{
	Int4	n,s,i;
	for(s=1; s<= NumElementarySets; s++){
		free(HyperPartition[s]); free(NameElementarySet[s]);
	} free(HyperPartition);
        for(n=1; n<= NumberBPPS; n++){
	   free(GroupsFG[n]); free(GroupsBG[n]);
	   if(GroupName[n]) free(GroupName[n]);
	   if(ArgmntStr[n]) free(ArgmntStr[n]);
	   for(i=0; i < nArgmnt[n]; i++){
		if(Argmntv[n][i]) free(Argmntv[n][i]);
	   } free(Argmntv[n]);
	}
	if(FileName) free(FileName);
	if(OutputCols) free(OutputCols);
}


