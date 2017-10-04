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

#include "scl_typ.h"

sch_typ	*sch_typ::GetSets(Int4 &II,double *&pval, Int4 *&II2ii, set_typ *&SetII)
{
	Int4	i,s,JJ,ii,bestm,MaxStr,*String;
	double	prob;

    sch_typ *sch=new sch_typ(size,MaxNumPttrns,resStart,resEnd,Clss,SetX,ResALL,AB);
    NEW(SetII, size+5, set_typ); NEW(II2ii, size+5,Int4); NEW(pval, size+5, double);
    for(II=0; this->DeleteMin(ii, bestm,MaxStr,String,prob) != 0; ){
	if(MaxStr <= 1){ if(String) free(String); continue; }  // with -K= option can return MaxStr==1.
	II++; pval[II]=prob;
	SetII[II]=MakeSet(resEnd+5); ClearSet(SetII[II]);
	for(i=0,s=1; s < MaxStr; s++){
		if(!MemberSet(String[s],SetX)) continue; else i++;
		JJ=ResidueID(ResALL[String[s]]); AddSet(JJ,SetII[II]);
	} JJ=ResidueID(ResALL[String[s]]); AddSet(JJ,SetII[II]);
	II2ii[II]=sch->Insert(ii,bestm,MaxStr,String,prob); 
    } return sch;
}

Int4	sch_typ::PrintSite(FILE *fptr, Int4 ii, char mode,FILE *tfp)
{
        if(ii < 1 || ii > size || SCI[ii]==0) return 0;
        Int4	i,s,JJ,*Str=SCI[ii]->String,nput=0;
	char	rc='X';
// if(mode=='m') PutSet(stdout,SetX);
// if(mode=='W') PutSet(stderr,SetX);
	if(Str[1] > 0){
           rc = GetCharResidue(ResALL[Str[1]],AB); JJ=ResidueID(ResALL[Str[1]]);
	} else { rc='X'; JJ=0; }  // Using KeyAtmDist: -K= option.
#if 1
        if(mode=='m') fprintf(fptr,"%d: start=%c%d; p=%.3g(x%d=%.2g): ", SCI[ii]->ID,
			rc,JJ,SCI[ii]->pval,Bonferroni,SCI[ii]->pval*(double)Bonferroni);
#else
        if(mode=='m') fprintf(fptr,"%d: start=%c%d; p = %.3g: ",SCI[ii]->ID,rc,JJ,SCI[ii]->pval);
#endif
        for(i=0,s=1; s <= SCI[ii]->MaxStr; s++){
            if(!MemberSet(Str[s],SetX) && mode != 'W') continue; 
	    else if(mode == 'W' && s==1 && Str[s] < 1) continue; else i++;
	    rc = GetCharResidue(ResALL[Str[s]],AB); JJ=ResidueID(ResALL[Str[s]]);
	    if(mode=='W'){
		if(!MemberSet(Str[s],SetX)){
		  if(nput > 0) fprintf(fptr,",");
		  fprintf(fptr,"%c%d.W",rc,JJ); nput++; 
	        }
	    } else if(i <= SCI[ii]->bestm){
              // if(mode=='m') fprintf(fptr,"%c%d[%d:%d],",rc,JJ,s,Str[s]);
              if(mode=='m'){ 
		if(nput > 0) fprintf(fptr,",");
		fprintf(fptr,"%c%d",rc,JJ); nput++; 
	      } else {	// mode == 'v'??
		if(nput > 0) fprintf(fptr,",");
		char cl=Clss; 
		if(strchr(" AROYGCBMPWXDLTSK",cl)==NULL) cl='X';
		fprintf(fptr,"%c%d.%c",rc,JJ,cl); nput++;
	      }
	    }
        } JJ=SCI[ii]->MaxStr;
	if(Str[1]==0) { JJ--; } // then don't count first position in String.
	if(mode == 'm') fprintf(fptr," (%d/%d)\n",i,JJ); else fprintf(fptr,"\n");
#if 1
	if(tfp){
	    if(Bonferroni == 1){
		fprintf(tfp,"-\t1\t%.3g\t%d\t%d\t%d\t%d\n",
		  SCI[ii]->pval,MaxNumPttrns,JJ,i,resEnd-resStart+1);
	    } else {
		fprintf(tfp,"%.3g\t%d\t%.3g\t%d\t%d\t%d\t%d\n",
		  SCI[ii]->pval,Bonferroni,SCI[ii]->pval*(double)Bonferroni,
		  MaxNumPttrns,JJ,i,resEnd-resStart+1);
	    }
	}
#endif
// fprintf(stderr,"-->%d: maxstr=%d; s=%d.\n", ii,SCI[ii]->MaxStr,s);
	return s;
}

set_typ	sch_typ::Put(FILE *fptr, FILE *vsifp,FILE *efp,FILE *tfp,char FullCluster, char *aStr, 
			double cutoff, BooLean PutAll, BooLean TakeUnion)
{
	Int4	i,j,s,II,JJ,ii,bestm,MaxStr,*String,*II2ii;
	double	prob,*pval; 
	BooLean notdone=TRUE;
	set_typ	*SetII,rtnSet=0;

    sch_typ *sch=this->GetSets(II,pval,II2ii,SetII);
// for(i=1; i <= II; i++) fprintf(stderr,"%d: ii=%d.\n",i,II2ii[i]);
    sch->Bonferroni=this->Bonferroni;
    if(II >= 1){
	Int4 N,si,sj,xi,xj,xx,nu,ni,red,black,out,outred,all;
	double D,d,pp;
	// 1. Count the number of Sets...eliminate empty sets.
	for(N=0,i=1; SetII[i]; i++){
	   N++; xi=CardSet(SetII[i]); 
	   if(xi == 0){ NilSet(SetII[i]); SetII[i]=0; } 
	} assert(II == N); 
	// 2. Eliminate  worst overlapping sets...
	BooLean	*Skip=0; NEW(Skip,N+5,BooLean);
	for(i=1; i < N; i++){
	   if(Skip[i]) continue;
	   if(SetII[i]==0) continue; 
	   xi=CardSet(SetII[i]); 
	   for(j=i+1; j <= N; j++){
	        if(Skip[j]) continue;
	        if(SetII[j]==0) continue; 
		xj=CardSet(SetII[j]);
		nu=CardUnionSet(SetII[i],SetII[j]);
		ni=CardInterSet(SetII[i],SetII[j]);
		if(ni==nu){	// i && j are identical sets; j > i so use i;
		   sch->Bonferroni--;
#if 0
// PutSet(stderr,SetII[i]); // PutSet(stderr,SetII[j]);
fprintf(stderr,"%d:%d.Bonferroni = %d; id=%d;pval=%.3g.\n",
	i,j,sch->Bonferroni,sch->SCI[II2ii[i]]->ID,pval[j]);
#endif
		   Skip[j]=TRUE;
		   if(PutAll) continue;	// show identical sets also.
		   NilSet(SetII[j]); SetII[j]=0; 
		} // else if(PutAll) continue;
#if 0
	        else if(ni==xi){	// i is a subset of j;
		   // this should always be true...pval[j] should be largest.
		   if(pval[i] <= pval[j]){ NilSet(SetII[j]); SetII[j]=0; }
		}
		if(SetII[i]==0) continue;
		d=(double)ni/(double)xi; D=(double)ni/(double)xj;
		if(d > 0.66 || D > 0.66){	// sizable overlap...
		   // if i is smaller and more significant then delete j cluster.
		   if(pval[i] <= pval[j]){ NilSet(SetII[j]); SetII[j]=0; }
		}
#endif
	   }
	} free(Skip);
	for(i=1; i <= N; i++){
	    double pv=(double)sch->Bonferroni * pval[i];
	    if(SetII[i] != 0 && pv > cutoff){ NilSet(SetII[i]); SetII[i]=0; }
	}
	for(i=1; i <= N; i++){
	    if(SetII[i] == 0) sch->Remove(II2ii[i]); 
	    // else if(pval[i] < cutoff){ NilSet(SetII[i]); SetII[i]=0; }
	}
	// 3. Combine surviving clusters into disjoint sets...
	ds_type sets = DSets(II);
	for(i=1; i < N; i++){
	   if(SetII[i]==0) continue;
	   xi=CardSet(SetII[i]); si = findDSets(i,sets);
	   for(j=i+1; j <= N; j++){
	        if(SetII[j]==0) continue; 
		xj=CardSet(SetII[j]);
		nu=CardUnionSet(SetII[i],SetII[j]);
		ni=CardInterSet(SetII[i],SetII[j]);
		d=(double)ni/(double)nu;
#if 0
		if(xi < xj){ red=xi; out=xj; } else { red=xj; out=xi; }
		all=CardSet(SetX); black=all-red; outred=ni; 
		pp=CumHypGeomProb(red,black,out,outred);
		if(pp < 0.01){
		   if(0) fprintf(fptr,"%d_%d: R=%d; B=%d; out=%d/%d; outred=%d; p = %.3g\n",
				i,j,red,black,out,all,outred,pp);
#else
		d=(double)ni/(double)xi; D=(double)ni/(double)xj;
		if(d > 0.50 || D > 0.50){	// sizable overlap...
		   if(0) fprintf(fptr,"%d.%d: d=%.3f; D=%.3f.\n", i,j,d,D);
#endif
		   sj = findDSets(j,sets);
		   if(si != sj){ si=sj=linkDSets(si,sj,sets); }
		}
	   }
	}
	Int4	NG,*grpID; NEW(grpID, N+5, Int4);
	for(NG=0,i=1; i <= N; i++){
	  if(SetII[i] && findDSets(i,sets) == i){
	    for(NG++,j=1; j <= N; j++){
	        if(SetII[j] && findDSets(j,sets) == i){ grpID[j]=NG; }
	    }
	  }
	} // PutDSets(stdout,sets); 
	NilDSets(sets);

     if(NG > 1) fprintf(fptr,"----- Class %c (%d subgroups) ----\n",Clss,NG);
     else fprintf(fptr,"----- Class %c ----\n",Clss);
     Int4	g,*NIG=0; NEW(NIG,N+5,Int4); 
     for(g=1; g <= NG; g++){ for(i=1; i <= N; i++){ if(grpID[i]==g) NIG[g]++; } }
#if 0	// DEBUG...
     for(g=1; g <= NG; g++){
	if(NIG[g] > 0) fprintf(stderr,"Grp %d: ",g);
	for(i=1; i <= N; i++){ if(grpID[i]==g) fprintf(stderr,"%d ",i);
	} fprintf(stderr,"(%d)\n",NIG[g]);
     }
#endif
     for(g=1; g <= NG; g++){
        if(NG > 1) fprintf(fptr,"----- subgroup %c%d ----\n",Clss,g);
	// 4. Print results...
	rtnSet=CopySet(SetX);
// PutSet(stderr,rtnSet);
	Int4	*cnsHits=0; 	// Consensus set.
	if(NIG[g] > 1){ NEW(cnsHits,SetN(SetX)+5,Int4);  }
	for(j=0,i=1; i <= N; i++){ 
	        if(SetII[i]==0) continue; 
		if(grpID[i] != g) continue; else j++;
		ii=II2ii[i];	// PutSet(stderr,rtnSet);
		Int4 rtn=0;
		if(j==1 && tfp != 0){
		   Int4 cl=10;
		   const char *category[]={ 0, "superfam","family","subfamily","catalytic","ligand",
			"structure","trace","classless","subtype","unknown" };
		   switch(Clss){
			case 'M': cl=1; break;
			case 'R': cl=2; break;
			case 'O': cl=3; break;
			case 'C': cl=4; break;
			case 'L': cl=5; break;
			case 'S': cl=6; break;
			case 'T': cl=7; break;
			case 'E': cl=8; break;
			default: cl=10; break;
		   } fprintf(tfp,"%s\tnone\t",category[cl]); 
		   rtn=sch->PrintSite(fptr,ii,'m',tfp);
		} else rtn=sch->PrintSite(fptr,ii,'m',0);
		if(rtn > 0){
		   // fprintf(stderr,"-->%d: Set %d; ii=%d.\n", i,grpID[i],II2ii[i]);
		   for(s=1; s <= sch->SCI[ii]->MaxStr; s++) DeleteSet(sch->SCI[ii]->String[s],rtnSet);
		   if(cnsHits) for(s=1; s <= sch->SCI[ii]->MaxStr; s++){
			 Int4 x=sch->SCI[ii]->String[s]; 
			 if(MemberSet(x,SetX)) cnsHits[x]++;
		   }
#if 1	// print out categories for Aravind: 'R' 'I' 'T' ... root one two...
		   if(aStr){
		     for(s=1; s < SetN(SetX); s++){
			if(MemberSet(s,SetX) && !MemberSet(s,rtnSet)) aStr[s]=Clss; 
		     }
		   } // PutSet(stdout,rtnSet); PutSet(stdout,SetX); // Invisible residues skipped.
#endif
		}
	        if(cnsHits == 0 && vsifp){
	          sch->PrintVSI(vsifp,II2ii[i]); 
	          if(notdone && Clss == FullCluster){
		     if(sch->PrintWhiteVSI(vsifp,II2ii[i]) > 0){
			// fprintf(stderr,"Done outputting full cluster (%d)\n",II2ii[i]);
			notdone=FALSE; 
		     }
	          }
		}
	} 
	if(cnsHits){	// consensus set for current group...
	    Int4 nput=0; fprintf(fptr,"  Consensus: ");
	    for(s=1; s < SetN(SetX); s++){ 
		if(cnsHits[s] > 0){
	           if(!TakeUnion){
		     double dd=(double)cnsHits[s]/(double)NIG[g];
		     if(dd < 0.50) continue;
		   }
           	   char rc = GetCharResidue(ResALL[s],AB); Int4 JJ=ResidueID(ResALL[s]);
		   if(nput > 0){ fprintf(fptr,","); if(vsifp) fprintf(vsifp,","); }
		   if(vsifp) fprintf(vsifp,"%c%d.%c",rc,JJ,Clss);
		   fprintf(fptr,"%c%d",rc,JJ); nput++;
		   if(0) fprintf(fptr," --> %d: %d\n",s,cnsHits[s]);
		}
	    } if(vsifp) fprintf(vsifp,"\n"); fprintf(fptr," (%d)\n",nput); 
	    free(cnsHits); cnsHits=0;
	}
      } free(grpID); free(NIG);
    }
    fflush(fptr); fflush(vsifp);
    for(i=1; i <= size; i++){ if(SetII[i]) NilSet(SetII[i]); }
    free(SetII);  free(II2ii); free(pval); delete sch;
    return rtnSet;
}

