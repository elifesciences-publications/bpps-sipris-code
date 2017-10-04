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

#include "afnio.h"
#include "sequence.h"
#include "my_ncbi.h"
#include "residues.h"

#if 0
sap_typ	fhp_typ::ViterbiToGSAP(...,double &Overall_Evalue)
// 
{
	:	:	:	:
	hmm_path = Viterbi(...,...);
	:	:	:	:
	sap_typ sap=PathHMMToGSAP(qE,sE,hmm_path,num_rpts,Evalue,bit_score,score);
	:	:	:	:
	return sap;
}

sap_typ MakeGSeqAlign(Int2 numseg, UInt4 query_id,
        UInt4 subject_id, e_type qE, e_type sE, Int4 *starts,
        Int4 *lens)
// not used in psiblast code yet...
{
        sap_typ sap = GSeqAlignNew();
        sap->dim =2; // only two dimention alignment.
        dsp_typ dsp = GDenseSegNew(); // Denseg Object is only type for GSeqAlign.
        dsp->dim = 2; dsp->numseg = numseg;
        dsp->subject_id = SeqI(sE); dsp->sE=sE;
        dsp->query_id=0; dsp->qE=qE;
        dsp->starts = starts; dsp->lens = lens;
        sap->segs = dsp; sap->next = NULL;
        return sap;
}

#endif

static Int4	CountNumSegs(register char *hmm_path)
{
        register Int4	i,numseg = 0;
	register char	c,state;

	if(hmm_path[0] != 'B') print_error("CountNumSegs() input error");
        for(state='B',i=1; (c=hmm_path[i]); i++){
             switch(c){
		case 'm': if(state!='m') numseg++; break;
		case 'd': if(state!='d') numseg++; break;
		case 'i': if(state!='i') numseg++; break;
		case 'E': return numseg; break;
		case 'n':	// N-terminal insertions...ignore.
		case 'c':	// C-terminal insertions...ignore.
		default: print_error("CountNumSegs() input error");
		   break;
	     } state=c;
		
	} print_error("CountNumSegs() input error");
}

#if 0
print_error(char *msg, char *hmm_path)
{
}
#endif

sap_typ PathHMMToGSAP(e_type hmmE, e_type sE, char *hmm_path, double *Evalue,
	double *bit_score,Int4 score,Int4 *scores)
// Create a sap_typ from a path through an HMM.
// Uses Eddy's Plan7 architecture
// State types: 'S','n','B','m','d','i','E','c','j','T'.
// Uses a finite state machine to parse the states.
// Input string ends in 0.
{
        Int4	i,j,seg,Nrpts,rpt,numseg = 0;
	char	state,c;
        Int4	*starts, *lens;
	Int4	qs,ss,L;
	sap_typ	head=0,sap=0;

	if(hmm_path[0] != 'S') print_error("PathHMMToGSAP() parse error 1");
        for(Nrpts=0,i=1; (c=hmm_path[i]); i++){
		if(c=='B') Nrpts++;
	} if(Nrpts==0) return NULL;
	starts=lens=0;
	
        for(j=seg=qs=ss=rpt=0,state='S',i=1; (c=hmm_path[i]); i++){
             switch(c){
		case 'n':	// N-terminal insertions..
		   if(!strchr(" nS",state)) print_error("PathHMMToGSAP() parse error 2");
		   ss++; 
		   break;
		case 'B':	// Begin state; encountered a repeat...
		   if(!strchr(" njES",state)){
	     		fprintf(stderr,"state=%c; c = %c\npath=%s\n",state,c,hmm_path);
			print_error("PathHMMToGSAP() parse error 3");
		   }
		   rpt++;
		   numseg = CountNumSegs(hmm_path+i);
		   // fprintf(stderr,"rpt(%d) = %d segs\n",rpt,numseg);
		   if(numseg==0) print_error("PathHMMToGSAP() parse error 4");
		   NEW(starts,2*numseg+2,Int4); NEW(lens,numseg+2,Int4);
		   qs=0; seg=j=L=0;
		   break;
		case 'm':	// match state; only reachable from 'B','m','i' or 'd'.
		   if(!strchr(" dimB",state)) print_error("PathHMMToGSAP() parse error 5");
		   if(state!='m'){
			if(seg > 0){ lens[j]=L; j++; }
			starts[seg] = ss; seg++; starts[seg] = qs; seg++;
			L=0;
		   } qs++; ss++; L++;
		   break;
		case 'd':	// delete state; only reachable from 'd', 'm' or 'B'.
		   if(!strchr(" dmB",state)) { 
			fprintf(stderr,"current state = %c; last state = %c\n",
				c,state);
			print_error("PathHMMToGSAP() parse error 6");
		   }
		   if(state!='d'){
			if(seg > 0){ lens[j]=L; j++; }
			starts[seg] = -1; seg++; starts[seg] = qs; seg++;
			L=0;
		   } qs++; L++;
		   break;
		case 'i':	// insert state; only reachable from 'i' or 'm'.
		   if(!strchr(" im",state)) print_error("PathHMMToGSAP() parse error 7");
		   if(state!='i'){
			if(seg > 0){ lens[j]=L; j++; }
			starts[seg] = ss; seg++; starts[seg] = -1; seg++;
			L=0;
		   } ss++; L++;
		   break;
		case 'E':	// End state; only reachable from 'd',or 'm'.
		   if(qs != LenSeq(hmmE)) print_error("PathHMMToGSAP() seq length error 8");
		   if(!strchr(" dm",state)) print_error("PathHMMToGSAP() parse error 9");
		   if(seg > 0){ lens[j]=L; }
		   if(head){
			// sap->next= MakeGSeqAlign(numseg,SeqI(hmmE),SeqI(sE),hmmE,sE,starts,lens);
			sap->next= MakeGSeqAlign(numseg,SeqI(sE),SeqI(hmmE),sE,hmmE,starts,lens);
			sap=sap->next;
			sap->score = scores[rpt];
			sap->evalue = Evalue[rpt];
			sap->bit_score = bit_score[rpt];
		   } else {
			// head=sap=MakeGSeqAlign(numseg,SeqI(hmmE),SeqI(sE),hmmE,sE,starts,lens);
			head=sap=MakeGSeqAlign(numseg,SeqI(sE),SeqI(hmmE),sE,hmmE,starts,lens);
			sap->score = scores[rpt];
			sap->evalue = Evalue[rpt];
			sap->bit_score = bit_score[rpt];
		   } qs=0;
		   break;
		case 'j':	// repeat insert state; only reachable from 'j' or 'E'.
		   if(!strchr(" jE",state)) print_error("PathHMMToGSAP() parse error 10");
		   ss++;
		   break;
		case 'T':
		   if(!strchr(" cE",state)) print_error("PathHMMToGSAP() parse error 11");
        	   return head;
		   break;
		case 'c':	// Stop here and return head.
		   if(!strchr(" cE",state)) print_error("PathHMMToGSAP() parse error 12");
		   ss++;
        	   return head;
		   break;
		default: print_error("PathHMMToGSAP() parse error 13");
		   break;
	     }
	     if(ss > LenSeq(sE)) print_error("PathHMMToGSAP() subject seq error 14");
	     if(qs > LenSeq(hmmE)){
	     	fprintf(stderr,
		  "qs = %d > LenSeq(hmmE) = %d; state=%c; c = %c\n",
						qs,LenSeq(hmmE),state,c);
		print_error("PathHMMToGSAP() hmm seq error 15");
	     }
	     // fprintf(stderr,"state=%c; c = %c\n",state,c);
	     state=c;
        } print_error("PathHMMToGSAP() parse error 16");
}

#if 0
int	main(Int4 argc,char *argv[])
{
	e_type	hmmE=0,sE=0;
	Int4	length;
	a_type          AB=MkAlpha(AMINO_ACIDS,PROT_BLOSUM62);
	hmmE=MkEmptySeq(1,"hmm_seq hmm concensus seq",100);
	sE=MkEmptySeq(1,"sub_seq subject seq",300);
	char hmm_path[] = "SnnnBmmmiiimmmdddmmEjjjjBmmmiiimmmdddmmEcccT";
	sap_typ	sap=PathHMMToGSAP(hmmE,sE,hmm_path, 0,0,0);
	// PutGSeqAlign(stdout,sap,60,AB);
	PutGSeqAlignList(stdout,sap,60,AB);
	FreeGSeqAlignList(sap);
	NilAlpha(AB);
	NilSeq(hmmE); NilSeq(sE);
}
#endif

