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

#include "tax_typ.h"

Int4	PutPhylaSeqSet(FILE *fptr, FILE *sfp, ss_type data)
{ char k; Int4 n; return PutPhylaSeqSet(fptr, sfp, data, k,n,0); }

Int4	PutPhylaSeqSet(FILE *fptr, FILE *sfp, ss_type data, char &NumKing, Int4 &NumNull,set_typ Set)
// create tax_typ for this at some point?
{ 
	Int4	i;
	char	c='P',Kingdom=' ',*Phylum=0,phylum_str[100],k;
	a_type	A=SeqSetA(data);
	e_type	Seq;
	char	*phyla[1002];
	char	kingdom[1002],king;
	Int4	seqs_in_phylum[1002];
	UInt4	j,num_phyla=0;
    	char	*phylum;

	phyla[0]=AllocString("null"); kingdom[0]='U';
        for(i=0; i <= 1000; i++) seqs_in_phylum[i]=0;
        for(i=1; i <= NSeqsSeqSet(data); i++){
	        if(Set && !MemberSet(i,Set)) continue;
		Seq = SeqSetE(i,data);
    		phylum=PhylumSeq(Seq);
		king=KingdomSeq(Seq);
		if(king==0) king='U';
		if(phylum==0){
		      seqs_in_phylum[0]++;
		} else if(num_phyla > 0){
		   for(j=1; j <= num_phyla; j++){
			if(strcmp(phylum,phyla[j]) == 0){
				seqs_in_phylum[j]++; break;
			}
		   } 
		   if(j > num_phyla){
			num_phyla++; phyla[num_phyla]=phylum; 
			kingdom[num_phyla]=king; 
			seqs_in_phylum[num_phyla]=1;
			if(sfp) PutSeq(sfp,Seq,A);
		   }
		} else {
			num_phyla++;
			phyla[1]=phylum; 
			kingdom[1]=king; 
			seqs_in_phylum[1]=1;
			if(sfp) PutSeq(sfp,Seq,A);
		}
		assert(num_phyla < 1000);
	}
        if(fptr){
	   if(Set) fprintf(fptr,"\n%s(%d seqs): num_phyla = %d \n",
			NameSeqSet(data),CardSet(Set),num_phyla);
	   else fprintf(fptr,"\n%s(%d seqs): num_phyla = %d \n",
			NameSeqSet(data),NSeqsSeqSet(data),num_phyla);
	   fprintf(fptr,"=======================\n");
	}
	Int4 num_of_kingdoms=0;
	for(king='A'; king <= 'Z'; king++){
	 Int4 num_in_kingdom=0;
	 for(char let1='A'; let1<= 'Z'; let1++){
	  for(char let2='A'; let2<= 'Z'; let2++){
	    for(j=0; j <= num_phyla; j++){
		if(kingdom[j] != king) continue;
		char let= toupper(phyla[j][0]);
		if(let == let1){
		   let=toupper(phyla[j][1]);
		   if(let == let2){
		     num_in_kingdom+=seqs_in_phylum[j];
        	     if(fptr) fprintf(fptr,"  %s (%d)\n",phyla[j],seqs_in_phylum[j]);
		   }
		}
	    }
	  }
	 }
	 if(num_in_kingdom > 0){
	  num_of_kingdoms++;
	  switch(king){
	    case 'V': 
	      if(fptr) fprintf(fptr,"===== %d plants (%c) =====\n\n",
		num_in_kingdom,king); break;
	    case 'A': 
	      if(fptr) fprintf(fptr,"===== %d archaea (%c) =====\n\n",
		num_in_kingdom,king); break;
	    case 'E': 
	      if(fptr) fprintf(fptr,"===== %d protozoa (%c) =====\n\n",
		num_in_kingdom,king); break;
	    case 'B': 
	      if(fptr) fprintf(fptr,"===== %d bacteria (%c) =====\n\n",
		num_in_kingdom,king); break;
	    case 'M': 
	      if(fptr) fprintf(fptr,"===== %d metazoa (%c) =====\n\n",
		num_in_kingdom,king); break;
	    case 'F': 
	      if(fptr) fprintf(fptr,"===== %d fungi (%c) =====\n\n",
		num_in_kingdom,king); break;
	    default:
	      num_of_kingdoms--; // don't count these...
	      if(fptr) fprintf(fptr,"===== %d unknown(%c) =====\n\n",
		num_in_kingdom,king); break;
	  }
	 }
	}
#if 1
        if(fptr){
	   if(Set) fprintf(fptr,"\n%s(%d seqs): %d phyla; %d kingdoms.\n\n",
			NameSeqSet(data),CardSet(Set),num_phyla,num_of_kingdoms);
	   else fprintf(fptr,"\n%s(%d seqs): %d phyla; %d kingdoms.\n\n",
			NameSeqSet(data),NSeqsSeqSet(data),num_phyla,num_of_kingdoms);
	}
#else
	if(fptr) fprintf(fptr,"\n number of kingdoms = %d; phyla = %d\n\n",
		num_of_kingdoms,num_phyla);
#endif
	free(phyla[0]);
	NumKing=num_of_kingdoms;
	NumNull=seqs_in_phylum[0];
	return num_phyla;
}


Int4	*GetPhylaSeqSet(FILE *fptr, Int4  *NumPhyla, ss_type data)
// create tax_typ for this at some point?
// Merge this with routine above at some point....
{ 
	Int4	i;
	char	c='P',Kingdom=' ',*Phylum=0,phylum_str[100],k;
	a_type	A=SeqSetA(data);
	e_type	Seq;
	char	*phyla[1002];
	char	kingdom[1002],king;
	Int4	seqs_in_phylum[1002];
	UInt4	j,num_phyla=0;
	Int4	*phylum_rtn=0;
    	char	*phylum;

        for(i=0; i <= 1000; i++){ kingdom[i]=0; phyla[i]=0; seqs_in_phylum[i]=0; }
	phyla[0]=AllocString("null");
	NEW(phylum_rtn,NSeqsSeqSet(data) +2,Int4);
        for(i=1; i <= NSeqsSeqSet(data); i++){
		Seq = SeqSetE(i,data);
    		phylum=PhylumSeq(Seq);
		king=KingdomSeq(Seq);
		if(king==0) king='U';
		if(phylum==0){
		      seqs_in_phylum[0]++;
		      kingdom[0]='U';
		} else if(num_phyla > 0){
		   for(j=1; j <= num_phyla; j++){
			if(strcmp(phylum,phyla[j]) == 0){
		   		phylum_rtn[i]=j;
				seqs_in_phylum[j]++; break;
			}
		   } 
		   if(j > num_phyla){
			num_phyla++; phyla[num_phyla]=phylum; 
			kingdom[num_phyla]=king; 
			seqs_in_phylum[num_phyla]=1;
			phylum_rtn[i]=num_phyla;
		   } 
		} else {
			num_phyla++; assert(num_phyla == 1);
			phyla[1]=phylum; 
			kingdom[1]=king; 
			seqs_in_phylum[1]=1;
			phylum_rtn[i]=num_phyla;
		}
		assert(num_phyla < 1000);
	}
	if(fptr){
          fprintf(fptr,"\n%s(%d seqs): num_phyla = %d \n",
			NameSeqSet(data),NSeqsSeqSet(data),num_phyla);
	  fprintf(fptr,"=======================\n");
	}
	Int4 num_of_kingdoms=0;
	for(king='A'; king <= 'Z'; king++){
	 Int4 num_in_kingdom=0;
	 for(char let1='A'; let1<= 'Z'; let1++){
	  for(char let2='A'; let2<= 'Z'; let2++){
	    for(j=0; j <= num_phyla; j++){
		if(kingdom[j] != king) continue;
		char let= toupper(phyla[j][0]);
		if(let == let1){
		   let=toupper(phyla[j][1]);
		   if(let == let2){
		     num_in_kingdom+=seqs_in_phylum[j];
		     if(fptr) fprintf(fptr,"  %s (%d)\n",phyla[j],seqs_in_phylum[j]);
		   }
		}
	    }
	  }
	 }
	 if(num_in_kingdom > 0){
	  num_of_kingdoms++;
	  switch(king){
	    case 'V': 
	     if(fptr) fprintf(fptr,"===== %d plants (%c) =====\n\n",
		num_in_kingdom,king); break;
	    case 'A': 
	     if(fptr) fprintf(fptr,"===== %d archaea (%c) =====\n\n",
		num_in_kingdom,king); break;
	    case 'E': 
	     if(fptr) fprintf(fptr,"===== %d protozoa (%c) =====\n\n",
		num_in_kingdom,king); break;
	    case 'B': 
	     if(fptr) fprintf(fptr,"===== %d bacteria (%c) =====\n\n",
		num_in_kingdom,king); break;
	    case 'M': 
	     if(fptr) fprintf(fptr,"===== %d metazoa (%c) =====\n\n",
		num_in_kingdom,king); break;
	    case 'F': 
	     if(fptr) fprintf(fptr,"===== %d fungi (%c) =====\n\n",
		num_in_kingdom,king); break;
	    default:
	      num_of_kingdoms--; // don't count these...
	     if(fptr) fprintf(fptr,"===== %d unknown(%c) =====\n\n",
		num_in_kingdom,king); break;
	  }
	 }
	} if(fptr) fprintf(fptr,"\n number of kingdoms = %d; phyla = %d\n\n",
		num_of_kingdoms,num_phyla);
	*NumPhyla=num_phyla; free(phyla[0]);
	return phylum_rtn;
}

