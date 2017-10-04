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

#include "hspheap.h"

hhp_typ MkHSPHeap(Int4 hpsz)
{
	hhp_typ	H;

	NEW(H,1,hspheap_type);
	NEW(H->hsp_array,hpsz+2,hsp_typ);
	H->mH=Mheap(hpsz, 3);
	H->size = hpsz;
	return H;
}

void    NilHSPHeap(hhp_typ H)
{
	hsp_typ	hsp;
	double	key;

	while((hsp=DelMinHSPHeap(&key,H)) != NULL){ NilHSPTyp(hsp); }
        NilMheap(H->mH); free(H->hsp_array);
	free(H);
}

Int4    InsertHSPHeap(hsp_typ hsp, hhp_typ H)
{
        Int4    item;
	double key = hsp->score;
        item=InsertMheap(-key, H->mH);
        if(item==NULL) return NULL;
        if(H->hsp_array[item]!= NULL){ NilHSPTyp(H->hsp_array[item]); }
        H->hsp_array[item]=hsp; 
        return item;
}

hsp_typ SeeHSPHeap(Int4 item, hhp_typ H) { return H->hsp_array[item]; }

hsp_typ DelMinHSPHeap(double *key, hhp_typ H)
{
        *key= -MinKeyMheap(H->mH); 
	Int4 item = DelMinMheap(H->mH);
	hsp_typ M = H->hsp_array[item]; H->hsp_array[item]=NULL;
        return M; 
}


hsp_typ DelMaxHSPHeap(double *key, hhp_typ H)
{
        Int4    item;
	hsp_typ M;

        *key= -MaxKeyMheap(H->mH);
        item = DelMaxMheap(H->mH);
	M = H->hsp_array[item]; H->hsp_array[item]=NULL;
        return M; 
}

Int4	PurgeHSPHeap(hhp_typ H)
// Remove overlapping HSPs from heap
// WARNING: Assumes that all HSPs are from the same query & subject sequences.
{
    hsp_typ	hsp1,hsp2,*hsp_array;
    double	key;
    Int4 	i,j,n;

    n=nHSPHeap(H);
    hsp_array = new hsp_typ [n+2];
    // 1. store HSPs from highest to lowest score.
    for(i=1;(hsp1=DelMinHSPHeap(&key, H)) != 0; i++) hsp_array[i]=hsp1;
    // 2. discard HSPs contained within higher scoring HSPs
    for(i=1; i <= n; i++){
	if((hsp1=hsp_array[i]) != 0){
          for(j=i+1; j <= n; j++){
	     if((hsp2=hsp_array[j]) != 0){
    		// 2a. Discard hsp2 if contained within hsp1:
		if(WithinHSPTyp(hsp2,hsp1) || WithinHSPTyp(hsp1,hsp2)){
			NilHSPTyp(hsp2); hsp_array[j]=0;
		} // 2b. Discard hsp2 if substantial overlap with hsp1:
#if 1		// Adapted from Tom Madden's code: CheckGappedAlignmentsForOverlap( )
                else if((hsp2->query_start== hsp1->query_start &&
                          hsp2->subject_start == hsp1->subject_start) ||
                            (hsp1->query_end == hsp2->query_end &&
                              hsp1->subject_end == hsp2->subject_end)) {
// fprintf(stderr,"Overlap = %g\n",OverlapHSPTyp(hsp1,hsp2));
// PutHSPTyp(stderr, hsp1); PutHSPTyp(stderr,hsp2);
			NilHSPTyp(hsp2); hsp_array[j]=0;
		} 
#endif
#if 0		// AFN code...
		else if(OverlapHSPTyp(hsp1,hsp2) >= 0.50){
			NilHSPTyp(hsp2); hsp_array[j]=0;
		}
#endif
	     }
	  } // 3. Reinsert retained HSPs into heap
	  InsertHSPHeap(hsp1, H);
	}
    } 
    delete [ ] hsp_array;
    return nHSPHeap(H);
}

BooLean	WithinHSPHeap(Int4 qs, Int4 ss, Int4 qe, Int4 se, hhp_typ hhp)
// See whether ungapped HSP is contained within any existing HSP.
// WARNING: Assumes that all HSPs are from the same query & subject sequences.
{
    hsp_typ	hsp;

    if(EmptyHSPHeap(hhp)) return FALSE;
    for(Int4 i=1; i <= hhp->size; i++){
	hsp=SeeHSPHeap(i,hhp);
	if(hsp != 0){ // 2a. return TRUE if contained within hsp:
	   if(WithinHSP(qs,ss,qe,se,hsp)) return TRUE;
	}
    } return FALSE;  // i.e., not in existing hsp
}

//**************************** HSP_TYP ROUTINES ***********************

static BooLean CONTAINED_IN_HSP(Int4 qs1,Int4 qe1,Int4 qi,Int4 ss1,Int4 se1,Int4 si)
//   TRUE if qi is between qs1 and qe1; si between ss1 and se1.  Determines if the
//   coordinates are already in an HSP that has been evaluated. 
{ return (((qs1 <= qi && qe1 >= qi) && (ss1 <= si && se1 >= si))?TRUE : FALSE);}

BooLean WithinHSP(Int4 query_start, Int4 subject_start,
	Int4 query_end, Int4 subject_end, hsp_typ hsp)
// Is sq/ss..eq/ew within hsp?
{
    BooLean     start_is_contained=FALSE, end_is_contained=FALSE;

    if(CONTAINED_IN_HSP(hsp->query_start, hsp->query_end,
         query_start, hsp->subject_start,
         hsp->subject_end, subject_start) == TRUE) start_is_contained = TRUE;
    if(CONTAINED_IN_HSP(hsp->query_start, hsp->query_end,
         query_end, hsp->subject_start,
         hsp->subject_end, subject_end) == TRUE) end_is_contained = TRUE;
    return (start_is_contained && end_is_contained);
}

double	OverlapHSPTyp(hsp_typ hsp1,hsp_typ hsp2)
{
	double overlap;
	Int4	ss_min,ss_max,se_min,se_max;
	Int4	qs_min,qs_max,qe_min,qe_max;
	Int4	q_max,q_min,s_max,s_min;

	ss_min= MINIMUM(Int4,hsp1->subject_start,hsp2->subject_start);
	ss_max= MAXIMUM(Int4,hsp1->subject_start,hsp2->subject_start);

	se_min= MINIMUM(Int4,hsp1->subject_end,hsp2->subject_end);
	se_max= MAXIMUM(Int4,hsp1->subject_end,hsp2->subject_end);

	qs_min= MINIMUM(Int4,hsp1->query_start,hsp2->query_start);
	qs_max= MAXIMUM(Int4,hsp1->query_start,hsp2->query_start);

	qe_min= MINIMUM(Int4,hsp1->query_end,hsp2->query_end);
	qe_max= MAXIMUM(Int4,hsp1->query_end,hsp2->query_end);

	q_max = qe_max - qs_min;  q_min = qe_min - qs_max;
	s_max = se_max - ss_min;  s_min = se_min - ss_max;
	overlap = (double)(s_min + q_min)/(double)(s_max + q_max);
	// overlap = (double)(s_max + q_max)/(double)(s_min + q_min);
	return overlap;
}

BooLean	WithinHSPTyp(hsp_typ hsp, hsp_typ hsp1)
// Is hsp within hsp1?
{
    BooLean	start_is_contained=FALSE, end_is_contained=FALSE;

    if(CONTAINED_IN_HSP(hsp1->query_start, hsp1->query_end,
         hsp->query_start, hsp1->subject_start,
         hsp1->subject_end, hsp->subject_start) == TRUE) {
         start_is_contained = TRUE;
    }
    if(CONTAINED_IN_HSP(hsp1->query_start, hsp1->query_end,
         hsp->query_end, hsp1->subject_start,
         hsp1->subject_end, hsp->subject_end) == TRUE) {
         end_is_contained = TRUE;
    }
    return (start_is_contained && end_is_contained);
}

hsp_typ MakeHSPTyp(gab_typ gap_align, double e_value,double bit_score)
// call as: MakeHSPTyp(gab_typ gap_align,GappedScoreToEvalueSBP(hsp->score, sbp),
// GappedBitScoreSBP(hsp->score,sbp));
{
	hsp_typ hsp;
	NEW(hsp,1,hsp_type); 
        hsp->query_gapped_start=gap_align->q_start;
        hsp->subject_gapped_start=gap_align->s_start;
        hsp->query_start = gap_align->query_start;
        hsp->subject_start = gap_align->subject_start;
        // The end is one further for BLAST than for the gapped align.
        hsp->query_end = gap_align->query_stop + 1;
        hsp->subject_end = gap_align->subject_stop + 1;
        hsp->query_length = hsp->query_end - hsp->query_start;
        hsp->subject_length = hsp->subject_end - hsp->subject_start;
        hsp->score = gap_align->score;
        hsp->gap_info = gap_align->edit_block;
	hsp->e_value = e_value;
        hsp->p_value = GBlastKarlinEtoP(hsp->e_value);
	hsp->bit_score = bit_score;
#if 0
fprintf(stderr,"\n************ Eval = %g; Pval = %g; score = %d\n",
                        hsp->e_value,hsp->p_value,hsp->score);
#endif
        // only one alignment considered for blast[np].
        // This may be changed by LinkHsps for blastx or tblastn.
        hsp->number = 1; // hsp->num = 1;  // ?? in blastcore.c ??
	return hsp;
}

// Delete this eventually.
hsp_typ MakeHSPTyp(gab_typ gap_align, sbp_typ sbp)
{
	hsp_typ hsp;
	NEW(hsp,1,hsp_type); 
        hsp->query_gapped_start=gap_align->q_start;
        hsp->subject_gapped_start=gap_align->s_start;
        hsp->query_start = gap_align->query_start;
        hsp->subject_start = gap_align->subject_start;
        // The end is one further for BLAST than for the gapped align.
        hsp->query_end = gap_align->query_stop + 1;
        hsp->subject_end = gap_align->subject_stop + 1;
        hsp->query_length = hsp->query_end - hsp->query_start;
        hsp->subject_length = hsp->subject_end - hsp->subject_start;
        hsp->score = gap_align->score;
        hsp->gap_info = gap_align->edit_block;
	hsp->e_value = GappedScoreToEvalueSBP(hsp->score, sbp);
        hsp->p_value = GBlastKarlinEtoP(hsp->e_value);
	hsp->bit_score = GappedBitScoreSBP(hsp->score,sbp);
#if 0
fprintf(stderr,"\n************ Eval = %g; Pval = %g; score = %d\n",
                        hsp->e_value,hsp->p_value,hsp->score);
#endif
        // only one alignment considered for blast[np].
        // This may be changed by LinkHsps for blastx or tblastn.
        hsp->number = 1; // hsp->num = 1;  // ?? in blastcore.c ??
	return hsp;
}

void	PutHSPTyp(FILE *fp, hsp_typ hsp)
{
	fprintf(fp,"number = %d\n",hsp->number);
	fprintf(fp,"score = %d\n",hsp->score);
	fprintf(fp,"bit score = %.0f\n",hsp->bit_score);
	fprintf(fp,"e_value = %g; p_value = %g\n",hsp->e_value,hsp->p_value);
	fprintf(fp,"query: start = %d; length = %d; end = %d\n",
		hsp->query_start,hsp->query_length,hsp->query_end);
	fprintf(fp,"subject: start = %d; length = %d; end = %d\n",
		hsp->subject_start,hsp->subject_length,hsp->subject_end);
	fprintf(fp,"query: gapstart = %d; subject: gapstart = %d\n",
			hsp->query_gapped_start,hsp->subject_gapped_start);
	fprintf(fp,"\n");
}

void	NilHSPTyp(hsp_typ hsp)
{
	if(hsp->gap_info) GXEBDelete(hsp->gap_info); 
	if(hsp) free(hsp); 
}

