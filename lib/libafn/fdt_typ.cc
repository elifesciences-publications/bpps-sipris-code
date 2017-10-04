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

#include "fdt_typ.h"

void	fdt_typ::Put(FILE *fp)
{
	FILE	*efp=stderr;
	Int4	o,i,j,n,x,total,*TotalI;
	Int4	NumCorrespond=0,NumIdent=0,NumLeaves=0,NumNew=0,Num50=0,Num80=0;
	NEW(TotalI,NumIRows+3,Int4);
	dh_type   dH=dheap(NumIRows,3);
	set_typ Set=MakeSet(SetN(ParentNodeI[1]));
	set_typ MiscSet=MakeSet(SetN(ParentNodeI[1]));
	for(o=1; o<=NumORows; o++){
#if 1	// eliminate 'o' characters.
	   char *tmp=AllocString(orow[o]);
	   for(n=0,i=0; tmp[i] != 0; i++){ if(tmp[i]=='o') tmp[i]=' '; else if(tmp[i] == '+') n++; }
	   if(TypeOfSetO[o] == '=') fprintf(fp," %s %3d.Reject (%d)",tmp,o,Hits[o][0]);
	   // else if(TypeOfSetO[o] == '?' && n <= 2){
	   else if(n <= 2){
	      fprintf(fp," %s %3d.%s%c (%d)",tmp,o,NameORow[o],TypeOfSetO[o],Hits[o][0]);
	   } else fprintf(fp," %s  %3d.%s%c (%d)",tmp,o,NameORow[o],TypeOfSetO[o],Hits[o][0]);
	   free(tmp);
#else
	   if(TypeOfSetO[o] == '=') fprintf(fp," %s %3d.Reject (%d)",orow[o],o,Hits[o][0]);
	   else fprintf(fp," %s %3d.%s%c (%d)",orow[o],o,NameORow[o],TypeOfSetO[o],Hits[o][0]);
#endif
	   if(TypeOfSetO[o] == '!' && Hits[o][0] > 0) NumLeaves++;
	   Int4 misc_total=0,other_misc_total=0; ClearSet(Set); ClearSet(MiscSet);
	   for(i=1; i<=NumIRows; i++){
		if(Hits[o][i] > 0){ 
		    if(TypeOfSetI[i] == '!') {
			insrtHeap(i,(keytyp) -Hits[o][i],dH); 
		        if(Hits[o][i] > 2) AddSet(i,Set);
		    } else if(TypeOfSetI[i] == '?') AddSet(i,MiscSet);
		    TotalI[i]+=Hits[o][i];
		}
	   }
	   if(Hits[o][0] > 0) fprintf(fp,"->");
	   for(total=0,x=0; ((i=delminHeap(dH)) != 0); ){
		x++; if(x < 5 && Hits[o][i] > 1 || x < 2){
		   fprintf(fp," %s(%d)",NameIRow[i],Hits[o][i]); total+=Hits[o][i]; 
		}
		if(x == 1 && Hits[o][0] > 0 && TypeOfSetO[o] == '!'){
		   double d=(double)Hits[o][i]/(double)Hits[o][0];
		   if(d > 0.90) NumCorrespond++;
		   if(d > 0.80) Num80++;
		   if(d > 0.50) Num50++;
		   if(Hits[o][i]==Hits[o][0]) NumIdent++;
		}
	   } 
	   for(i=1; i<=NumIRows; i++){
		if(MemberSet(i,MiscSet)){	// is this a detected Misc Set?
	   	  for(j=1; j <= NumIRows; j++){  // ...that is a parent of a detected leaf set?
		     if(MemberSet(j,Set)){
			if(MemberSet(i,ParentNodeI[j])){
			 misc_total+=Hits[o][i]; break; // don't count more than once...
		     	} 
		     }
		  } if(j > NumIRows) other_misc_total+=Hits[o][i]; // not found above.
		}
	   }
	   Int4 rand=Hits[o][NumIRows];
	   Int4 leftover=Hits[o][0] - total - misc_total -other_misc_total -rand;

	   if(total == 0 && TypeOfSetO[o] == '!' && Hits[o][0] > 0){
		 fprintf(fp," New!"); NumNew++;
	   }
	   if(misc_total > 0 && rand > 0) fprintf(fp," [%d parent; %d reject]",misc_total,rand);
	   else if(misc_total > 0) fprintf(fp," [%d parent]",misc_total);
	   else if(rand > 0) fprintf(fp," [%d reject]",rand);
	   if(other_misc_total > 0) fprintf(fp," (%d other internal)",other_misc_total);
	   if(leftover > 0) fprintf(fp," + %d.\n",leftover); else fprintf(fp,".\n");
	   // PutSet(fp,ParentNodeO[o]);
	} Nildheap(dH); NilSet(Set); NilSet(MiscSet); free(TotalI);  
#if 0
	fprintf(efp,"  Report50N --> %d/%d (%.1f%c) > 50%c identical (non-new) leaf sets.\n",
			Num50,NumLeaves-NumNew,
			100.0*(double)Num50/(double)(NumLeaves-NumNew),'%','%');
	fprintf(efp,"  Report50I --> %d/%d (%.1f%c) > 50%c identical sets.\n\n",Num50,NumORows-1,
		 	100.0*(double)Num50/(double)(NumORows-1),'%','%');

	fprintf(efp,"  ReportA --> %d/%d (%.1f%c) > 80%c identical (non-new) leaf sets.\n",
			Num80,NumLeaves-NumNew,
			100.0*(double)Num80/(double)(NumLeaves-NumNew),'%','%');
	fprintf(efp,"  ReportB --> %d/%d (%.1f%c) > 80%c identical sets.\n\n",Num80,NumORows-1,
		 	100.0*(double)Num80/(double)(NumORows-1),'%','%');

	fprintf(efp,"  ReportC --> %d/%d (%.1f%c) > 90%c identical (non-new) leaf sets.\n",
			NumCorrespond,NumLeaves-NumNew,
			100.0*(double)NumCorrespond/(double)(NumLeaves-NumNew),'%','%');
	fprintf(efp,"  ReportD --> %d/%d (%.1f%c) > 90%c identical sets.\n\n",NumCorrespond,NumORows-1,
		 	100.0*(double)NumCorrespond/(double)(NumORows-1),'%','%');

	fprintf(efp,"  ReportE --> %d/%d (%.1f%c) identical (non-new) leaf sets.\n",
			NumIdent,NumLeaves-NumNew,
			100.0*(double)NumIdent/(double)(NumLeaves-NumNew),'%');
	fprintf(efp,"  ReportF --> %d/%d (%.1f%c) identical sets.\n",NumIdent,(NumORows-1),
			100.0*(double)NumIdent/(double)(NumORows-1),'%');
#endif
}

void	fdt_typ::FindMisc() 
{
   Int4 i,o,c,r;
   NEW(ParentNodeI,NumIRows+3,set_typ);
   NEW(ParentNodeO,NumORows+3,set_typ);
   for(i=1; i <= NumIRows; i++){
      ParentNodeI[i]=MakeSet(NumIRows+3); ClearSet(ParentNodeI[i]);
      for(c=1; c < i; c++){
	 if(irow[i][c]=='+') AddSet(c,ParentNodeI[i]);
      }
   }
   for(o=1; o <= NumORows; o++){
      ParentNodeO[o]=MakeSet(NumORows+3); ClearSet(ParentNodeO[o]);
      for(c=1; c < o; c++){
	 if(orow[o][c]=='+') AddSet(c,ParentNodeO[o]);
      }
   }
}

void	fdt_typ::Free()
{
	for(Int4 i=1; i<=NumIRows; i++){ NilSet(ParentNodeI[i]); } 
	for(Int4 o=1; o<=NumORows; o++){ NilSet(ParentNodeO[o]); free(Hits[o]); }
	free(Hits); free(ParentNodeI); free(ParentNodeO);
}


