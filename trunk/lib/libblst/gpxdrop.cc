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

#include "gpxdrop.h"
// From gapxdrop.h by Gennadiy Savchuk, Jinqhui Zhang, Tom Madden
// [From code originally supplied by Webb Miller's lab].
// Contents: Functions to perform a gapped alignment on two sequences.

/* A PACKAGE FOR LOCALLY ALIGNING TWO SEQUENCES WITHIN A BAND:

   To invoke, call SEMI_G_ALIGN(A,B,M,N,W,G,H,S,X,&EI,&EJ, score_only, positionBased).
   The parameters are explained as follows:
	A, B : two sequences to be aligned
	M : the length of sequence A
	N : the length of sequence B
	W : scoring table for matches and mismatches
	G : gap-opening penalty
	H : gap-extension penalty
	X : maximum drop off
	S : script for DISPLAY routine
	*EI : ending position of sequence A in the optimal local alignment
	*EJ : ending position of sequence B in the optimal local alignment
	score_only: indicate if ony score is needed. 1--score only 0--alignment
	            is also needed.

	Use best_score to cut unfeasible points. quadratic space and row-wise
*/

/* Append "Delete k" op */
#define DEL_(k) \
data.last = (data.last < 0) ? (data.sapp[-1] -= (k)) : (*data.sapp++ = -(k));

/* Append "Insert k" op */
#define INS_(k) \
data.last = (data.last > 0) ? (data.sapp[-1] += (k)) : (*data.sapp++ = (k));

/* Append "Replace" op */
#define REP_ \
{ data.last = *data.sapp++ = 0; }

/* Divide by two to prevent underflows. */
// #define MININT (INT4_MIN-1)/2	// ?? original form...
#define MININT INT4_MIN/2
#define REPP_ \
{ *data.sapp++ = MININT; data.last = 0; }

Int4	AlnEndFlagGXEB( ) {	return MININT; }

/* Dynamic Programming structure. */
typedef struct DP {
  Int4 CC, DD, FF;	/* Values for gap opening and extensions (?). */
} *dp_ptr, dp_node;

typedef struct {
  Int4 * sapp;			/* Current script append ptr */
  Int4  last;	
} data_t;

/* #define	CHUNKSIZE	1048576 */
#define	CHUNKSIZE	2097152

// static gxdsa_typ GapXDropGetState(gxdsa_typ *head, Int4 length)

gxdsa_typ GapXDropGetState(gxdsa_typ *head, Int4 length)
{
	gxdsa_typ	retval, var, last;
	Int4				chunksize = MAX(CHUNKSIZE, length + length/3);

	length += length/3;	/* Add on about 30% so the end will get reused. */
	retval = NULL;
	if (*head == NULL) {
		retval = (gxdsa_typ) GMemNew(sizeof(gxdsa_type));
		retval->state_array = (unsigned char*) GMemNew(chunksize*sizeof(unsigned char));
		retval->length = chunksize;
		retval->used = 0;
		*head = retval;
	} else {
		var = *head;
		last = *head;
		while (var) {
			if (length < (var->length - var->used)) {
				retval = var; break; }
			else if (var->used == 0)
			{ /* If it's empty and too small, replace. */
				var->state_array = (unsigned char*) GMemFree(var->state_array);
				var->state_array = (unsigned char*) GMemNew(chunksize*sizeof(unsigned char));
				var->length = chunksize;
				retval = var;
				break;
			}
			last = var; var = var->next;
		}
		if (var == NULL) {
			retval = (gxdsa_typ) GMemNew(sizeof(gxdsa_type));
			retval->state_array = (unsigned char*) GMemNew(chunksize*sizeof(unsigned char));
			retval->length = chunksize;
			retval->used = 0;
			last->next = retval;
		}
	}
	if (retval->state_array == NULL) print_error("state array is NULL");
	return retval;
}

// static Boolean GapXDropPurgeState(gxdsa_typ state_struct)
Boolean GapXDropPurgeState(gxdsa_typ state_struct)
{
	while (state_struct) {
		memset(state_struct->state_array, 0, state_struct->used);
		state_struct->used = 0;
		state_struct = state_struct->next;
	} return TRUE;
}

gxdsa_typ GapXDropStateDestroy(gxdsa_typ state_struct)
{
	gxdsa_typ next;
	while (state_struct) {
		next = state_struct->next;
		GMemFree(state_struct->state_array);
		GMemFree(state_struct);
		state_struct = next;
	} return NULL;
}

Int4 ALIGN(unsigned char * A, unsigned char * B, Int4 M, Int4 N,
		Int4 * S, Int4 * pei, Int4 * pej, Int4 **sapp, 
		gab_typ gap_align, Int4 query_offset, Boolean reversed)
{ 
  data_t data;
  Int4 i, j, cb,  j_r, s, k;
  unsigned char st, std, ste;
  Int4 gap_open, gap_extend, decline_penalty;
  register Int4 c, d, e, m,t, tt, f, tt_start;
  Int4 best_score = 0;
  register Int4 *wa;
  register dp_ptr dp, dyn_prog;
  unsigned char **state, *stp, *tmp;
  unsigned char *state_array;
  Int4	**matrix;
  register Int4 X;
  gxdsa_typ state_struct;
  
  matrix = gap_align->matrix;
  *pei = *pej = 0;
  data.sapp = *sapp = S;
  data.last= 0;
  m = gap_align->gap_open + gap_align->gap_extend;
  decline_penalty = gap_align->decline_align;
  gap_open = gap_align->gap_open;
  gap_extend = gap_align->gap_extend;
  X = gap_align->x_parameter;

  if (X < m) X = m;
  if(N <= 0 || M <= 0) { *pei = *pej; return 0; }
  GapXDropPurgeState(gap_align->state_struct);

  j = (N + 2) * sizeof(dp_node);
  dyn_prog = (dp_ptr)malloc(j);

  state = (unsigned char**) malloc(sizeof(unsigned char *)*(M+1));
  dyn_prog[0].CC = 0;
  dyn_prog[0].FF = -m;
  c = dyn_prog[0].DD = -m;

  /* Protection against divide by zero. */
  if (gap_extend > 0)
  	state_struct = GapXDropGetState(&gap_align->state_struct, X/gap_extend+5);
  else state_struct = GapXDropGetState(&gap_align->state_struct, N+3);

  state_array = state_struct->state_array;
  state[0] = state_array;
  stp  = state[0];
  for(i = 1; i <= N; i++) {
    if(c < -X) break;
    dyn_prog[i].CC = c;
    dyn_prog[i].DD = c - m; 
    dyn_prog[i].FF = c-m;
    c -= gap_extend;
    stp[i] = 1;
  }
  state_struct->used = i+1;
  
  tt = 0;  j = i;
  for(j_r = 1; j_r <= M; j_r++) {
     /* Protection against divide by zero. */
    if (gap_extend > 0)
    	state_struct = GapXDropGetState(&gap_align->state_struct, j-tt+5+X/gap_extend);
    else
	state_struct = GapXDropGetState(&gap_align->state_struct, N-tt+3);
    state_array = state_struct->state_array + state_struct->used + 1;
    state[j_r] = state_array - tt + 1;
    stp = state[j_r];
    tt_start = tt; 
    if (!(gap_align->positionBased)) /*AAS*/
      wa = matrix[A[j_r]]; 
    else {
      if(reversed) wa = gap_align->posMatrix[M - j_r];
      else wa = gap_align->posMatrix[j_r + query_offset];
    } e = c = f= MININT;
    
    for (cb = i = tt, dp = &dyn_prog[i-1]; i < j; i++) {
	d = (++dp)->DD;
	if (d < f) { d=f; std=83;} else std=42;
	if (e < f) {e=f; ste = 83;} else ste = 41;
	if (c < d || c < e) {
	    if (d < e) { c = e; st = ste; } else { c = d; st = std; }
	    if (best_score - c > X) {
		c = dp->CC+wa[B[i+1]]; f = dp->FF;
		if (tt == i) tt++; else { dp->CC = MININT; }
	    } else {
		cb = i;
		e -= gap_extend; dp->DD = d-gap_extend;
		d = dp->FF; dp->FF = f-decline_penalty; f = d;
		d = dp->CC+wa[B[i+1]]; dp->CC = c; c=d;
	    }
	    st += 5;
	} else {
	    st = 0;
	    if (best_score - c > X){
		c = dp->CC+wa[B[i+1]]; f= dp->FF;
		if (tt == i) tt++; else { dp->CC =  MININT; }
	    } else {
		cb = i;
		if (c > best_score) { best_score = c; *pei = j_r; *pej = i; }
		if ((c-=m) > (d-=gap_extend)) { dp->DD = c; } 
		else { dp->DD = d; if (std == 83) st+=60; else st+=30; } 
		if (c > (e-=gap_extend)) { e = c; }
		else { if (ste> 50) st+=20; else st +=10; }
		c+=m; d = dp->FF;
		if (f < c-gap_open) { dp->FF = c-gap_open-decline_penalty; }
		else { dp->FF = f-decline_penalty; st+= 5; }
		f = d;
		d = dp->CC+wa[B[i+1]]; dp->CC = c; c = d;      
	    }
	}
	stp[i] = st;
    }
    if (tt == j) { j_r++; break;}
    if (cb < j-1) { j = cb+1;}
    else {
	while (e >= best_score-X && j <= N) {
	    dyn_prog[j].CC = e; 
	    dyn_prog[j].DD = e-m; dyn_prog[j].FF = e-gap_open-decline_penalty;
	    e -= gap_extend; stp[j] = 1;
	    j++;
	}
    }
    if (j <= N) {
	dyn_prog[j].DD = dyn_prog[j].CC= dyn_prog[j].FF = MININT;
	j++; 
    }
    state_struct->used += (MAX(i, j) - tt_start + 1);
  }
  i = *pei; j = *pej;
  tmp = (unsigned char*) malloc(i+j);
  for (s=0, c = 0; i> 0 || j > 0; c++) {
      t = state[i][j]; k  = t %5;
      if (s == 1) 
	  if ((t/10)%3 == 1) k = 1; else if ((t/10)%3 == 2) k = 3;
      if (s == 2)
	  if ((t/30) == 1) k = 2; else if ((t/30) == 2) k = 3;
      if (s == 3 && ((t/5)%2) == 1) k = 3;
      if (k == 1) { j--;} else if (k == 2) {i--;} else {j--; i--;}
      tmp[c] = s = k;
  }
  c--;
  while (c >= 0) {
      if (tmp[c] == 0) REP_
      else if (tmp[c] == 1) INS_(1)
      else if (tmp[c] == 3) REPP_			  
      else DEL_(1)     
      c--;
  }
  GMemFree(tmp); GMemFree(state); GMemFree(dyn_prog);
  *sapp = data.sapp;
  return best_score;
}

Int4 SEMI_G_ALIGN(unsigned char * A, unsigned char * B, Int4 M, Int4 N,
		Int4 * S, Int4 * pei, Int4 * pej, 
		Boolean score_only, Int4 **sapp, gab_typ gap_align,
		Int4 query_offset, Boolean reversed)
		
{ 
  dp_ptr dyn_prog;
  Int4 i, j, cb, j_r, g, decline_penalty;
  register Int4 c, d, e, m, tt, h, X, f;
  Int4 best_score = 0;
  Int4 * *matrix;
  register Int4 * wa;
  register dp_ptr dp;
  
  if(!score_only)
    return ALIGN(A, B, M, N, S, pei, pej, sapp, gap_align, query_offset, reversed);
  
  matrix = gap_align->matrix;
  *pei = *pej = 0;
  m = (g=gap_align->gap_open) + gap_align->gap_extend;
  h = gap_align->gap_extend;
  decline_penalty = gap_align->decline_align;
  X = gap_align->x_parameter;

  if (X < m) X = m;
  if(N <= 0 || M <= 0) return 0;

  j = (N + 2) * sizeof(dp_node);
  dyn_prog = (dp_ptr)malloc(j);

  dyn_prog[0].CC = 0; c = dyn_prog[0].DD = -m;
  dyn_prog[0].FF = -m;
  for(i = 1; i <= N; i++) {
    if(c < -X) break;
    dyn_prog[i].CC = c;
    dyn_prog[i].DD = c - m; 
    dyn_prog[i].FF = c-m;
    c -= h;
  }

  tt = 0;  j = i;
  for (j_r = 1; j_r <= M; j_r++) {
      if (!(gap_align->positionBased)) /*AAS*/
	  wa = matrix[A[j_r]]; 
      else {
	  if(reversed) wa = gap_align->posMatrix[M - j_r];
	  else wa = gap_align->posMatrix[j_r + query_offset];
      }
      e = c =f = MININT;
      for (cb = i = tt, dp = &dyn_prog[i]; i < j; i++) {
	  d = dp->DD;
	  if (e < f) e = f;
	  if (d < f) d = f;
	  if (c < d || c < e) {
	      if (d < e) { c = e; } else { c = d; }
	      if (best_score - c > X) {
		  c = dp->CC+wa[B[i+1]]; f = dp->FF;
		  if (tt == i) tt++; else { dp->CC = MININT;}
	      } else {
		  cb = i; e -= h; dp->DD =  d- h;
		  d = dp->CC+wa[B[i+1]]; dp->CC = c; c=d;
		  d = dp->FF; dp->FF = f-decline_penalty; f = d;
	      }
	  } else {
	      if (best_score - c > X){
		  c = dp->CC+wa[B[i+1]]; f= dp->FF;
		  if (tt == i) tt++; else { dp->CC = MININT;}
	      } else {
		  cb = i; 
		  if (c > best_score) { best_score = c; *pei = j_r; *pej = i; } 
		  if ((c-=m) > (d-=h)) { dp->DD = c; } else { dp->DD = d; } 
		  if (c > (e-=h)) { e = c; } 
		  c+=m;
		  d = dp->FF;
		  if (c-g>f) dp->FF = c-g-decline_penalty; else dp->FF = f-decline_penalty;
		  f = d;
		  d = dp->CC+wa[B[i+1]]; dp->CC = c; c = d;
	      }
	  }
	  dp++;
      }
      if (tt == j) break;
      if (cb < j-1) { j = cb+1;}
      else while (e >= best_score-X && j <= N) {
	  dyn_prog[j].CC = e; dyn_prog[j].DD = e-m;dyn_prog[j].FF = e-g-decline_penalty;
	  e -= h; j++;
      }
      if (j <= N) {
	  dyn_prog[j].DD = dyn_prog[j].CC = dyn_prog[j].FF = MININT; j++;
      }
  }
  GMemFree(dyn_prog);
  return best_score;
}

/* Allocates an EditScriptPtr and puts it on the end of
	the chain of EditScriptPtr's.  Returns a pointer to the
	new one.  */
static gxes_typ GXESNew(gxes_typ old)
{
	gxes_typ New;
	New = (gxes_typ) GMemNew(sizeof(gxes_type));
	if (old == NULL) return New;
	while (old->next != NULL) { old = old->next; }
	old->next = New;
	return New;
}

gxeb_typ GXEBNew(Int4 start1, Int4 start2)
{
	gxeb_typ edit_block;
	edit_block = (gxeb_typ) GMemNew(sizeof(gxeb_type));
	edit_block->start1 = start1; edit_block->start2 = start2;
	return edit_block;
}

gxeb_typ GXEBDelete(gxeb_typ edit_block)
{
	if (edit_block == NULL) return NULL;
	gxes_typ next,old=edit_block->esp;
 	while(old){ next = old->next; free(old); old=next; } 
	free(edit_block); return NULL;
}

gxeb_typ TracebackToGXEB(unsigned char * A, unsigned char * B, 
	Int4 M, Int4 N, Int4 * S, Int4 start1, Int4 start2)
/* Alignment display routine */
{

  register Int4 i, j, op, number_of_subs, number_of_decline;
  gxes_typ e_script;
  gxeb_typ edit_block;

  if (S == NULL) return NULL;

  i = j = op = number_of_subs = number_of_decline = 0;
  edit_block = GXEBNew(start1, start2);
  edit_block->esp = e_script = GXESNew(NULL);

  while(i < M || j < N) {
	op = *S;
#if 0   // AFN: diagnostic
        fprintf(stderr,"i = %d; j = %d; op = %d; s1 = %d; s2 = %d\n",i,j,op,start1,start2);
	assert(start1==start2);
#endif
	if (op != MININT && number_of_decline > 0) {
               e_script->op_type = GAPXALIGN_DECLINE;
               e_script->num = number_of_decline;
               e_script = GXESNew(e_script);
		number_of_decline = 0;
	}
	if (op != 0 && number_of_subs > 0) {
                        e_script->op_type = GAPXALIGN_SUB;
                        e_script->num = number_of_subs;
                        e_script = GXESNew(e_script);
                        number_of_subs = 0;
        }      
	if (op == 0) { i++; j++; number_of_subs++; }
	else if (op == MININT) { i++; j++; number_of_decline++; }
	else {
		if(op > 0) {
			e_script->op_type = GAPXALIGN_DEL;
			e_script->num = op;
			j += op;
			if (i < M || j < N)
                                e_script = GXESNew(e_script);
		} else {
			e_script->op_type = GAPXALIGN_INS;
			e_script->num = ABS(op);
			i += ABS(op);
			if (i < M || j < N)
                                e_script = GXESNew(e_script);
		}
    	} S++;
  }

  if (number_of_subs > 0) {
	e_script->op_type = GAPXALIGN_SUB;
	e_script->num = number_of_subs;
  } else if (number_of_decline > 0) {
        e_script->op_type = GAPXALIGN_DECLINE;
        e_script->num = number_of_decline;
  }
  return edit_block;
}

// static Int4 reverse_seq(unsigned char * seq, unsigned char * pos, unsigned char * target) 
Int4 reverse_seq(unsigned char * seq, unsigned char * pos, unsigned char * target) 
{
    unsigned char * c;

    for (c = pos; c >= seq; c--) *target++ = *c;
    *target = '\0';
    return (Int4) (pos-c);
}

gab_typ GABDelete(gab_typ gap_align)
// Destruct Function for gab_type.  If "state" != NULL, it's deallocated. 
{
	if (gap_align == NULL) return NULL;
	gap_align->state_struct = GapXDropStateDestroy(gap_align->state_struct);
	gap_align = (gab_typ) GMemFree(gap_align);
	return gap_align;
}

gab_typ GABNew(e_type qE, Int4 open, Int4 extend, Int4 **matrix, Int4 decl_aln,
	double x_parameter, Int4 min_query_align)
{
	gab_typ gap_align = GABNew(open, extend);
	gap_align->include_query = min_query_align;
	// query length from q_start that MUST be aligned.
	gap_align->query = SeqPtr(qE) + 1;
	gap_align->query_length = LenSeq(qE);
	gap_align->decline_align = decl_aln;
	gap_align->matrix=matrix;  // NOTE: may need to change sbp->matrix!!!
	gap_align->posMatrix=NULL; gap_align->positionBased = FALSE;
	gap_align->x_parameter=x_parameter;
	return gap_align;
}

gab_typ GABNew(Int4 open, Int4 extend)
// NOTE: ALIGN allocates and deallocates state memory. 
{
	gab_typ gap_align;
	gap_align = (gab_typ) GMemNew(sizeof(gab_type));
	if (gap_align == NULL) return NULL;
	gap_align->decline_align = INT2_MAX;
	gap_align->gap_open = open;
	gap_align->gap_extend = extend;
	return gap_align;
}

// AFN: the next two routines are ones that I made.
Int4    ScoreGAB(gab_typ gap_align){ return gap_align->score; }

sap_typ ExtractAlnGAB(gab_typ gap_align, double e_value, e_type sE, e_type qE,
        double PerNat)
{
        sap_typ sap;

        assert(gap_align->edit_block != NULL);
        assert(gap_align->edit_block->esp->num > 0);
        sap=GXEBToGSeqAlign(gap_align->edit_block,sE, qE);
        GXEBDelete(gap_align->edit_block); gap_align->edit_block=NULL;
        sap->evalue = e_value;
        sap->score = gap_align->score;
        sap->bit_score = (Int4) floor((1.442695*(gap_align->score/PerNat))+0.5);
        sap->next=NULL;
        return sap;
}

void    SetSubjectGAB(Int4 q_start, Int4 s_start, unsigned char *pssq, 
			Int4 len_ssq, gab_typ gab)
{
	gab->subject = pssq + 1; // starts at 0 not 1.
	gab->subject_length = len_ssq;
	gab->q_start=q_start; gab->s_start=s_start;
}

void    SetSubjectGAB(Int4 q_start, Int4 s_start, e_type sE, gab_typ gab)
{ SetSubjectGAB(q_start,s_start,SeqPtr(sE),LenSeq(sE),gab); }

Boolean PerformGappedAlignment(gab_typ gap_align)
{
	Boolean		found_start, found_end;
	Int4		q_length=0, s_length=0, middle_score;
	Int4		score_right, score_left, private_q_start, private_s_start;
	Int4		include_query, index;
	unsigned char	*q_left=NULL, *s_left=NULL;
	unsigned char	*query, *subject, *query_var, *subject_var;

	if (gap_align == NULL) return FALSE;
	found_start = FALSE; found_end = FALSE;

	query = gap_align->query;
	subject = gap_align->subject;
	include_query = gap_align->include_query;
	score_left = 0;
	if (gap_align->q_start != 0 && gap_align->s_start != 0) {
		found_start = TRUE;
#if 0		// afn: 1/22/08 edit to fix core dump
		q_left = (unsigned char *) malloc((gap_align->q_start+2)*sizeof(unsigned char));
		s_left = (unsigned char *) malloc((gap_align->s_start+2)*sizeof(unsigned char));
#else
		MEW(q_left,(gap_align->q_start+2),unsigned char);
		MEW(s_left,(gap_align->s_start+2),unsigned char);
#endif

		q_length = reverse_seq(query, query+gap_align->q_start, q_left);
		s_length = reverse_seq(subject, subject+gap_align->s_start, s_left);

		score_left = SEMI_G_ALIGN(q_left-1, s_left-1, q_length, s_length, NULL, &private_q_start, &private_s_start, TRUE, NULL, gap_align, gap_align->q_start, TRUE);
		gap_align->query_start = q_length - private_q_start;
		gap_align->subject_start = s_length - private_s_start;
	} else {
		q_length = gap_align->q_start;
		s_length = gap_align->s_start;
	}
	middle_score = 0;
	query_var = query+gap_align->q_start;
	subject_var = subject+gap_align->s_start;
	for (index=0; index<include_query; index++) {
		query_var++;
		subject_var++;
		if (!(gap_align->positionBased))  /*AAS*/
		  middle_score += gap_align->matrix[*query_var][*subject_var];
		else 
		  middle_score += MtrxGScoreGapAlign(gap_align,
				gap_align->q_start+1 + index,*subject_var);
	}
	score_right = 0;
	if (gap_align->q_start+include_query < gap_align->query_length 
		&& gap_align->s_start+include_query < gap_align->subject_length)
	{
		found_end = TRUE;
		score_right = SEMI_G_ALIGN(query+gap_align->q_start+include_query,
			subject+gap_align->s_start+include_query,
			gap_align->query_length-q_length-include_query,
			gap_align->subject_length-s_length-include_query,
			NULL, &(gap_align->query_stop),
			&(gap_align->subject_stop), TRUE,
			NULL, gap_align, gap_align->q_start+include_query, FALSE);
		gap_align->query_stop += gap_align->q_start+include_query;
		gap_align->subject_stop += gap_align->s_start+include_query;
	}

	if (found_start == FALSE) {	/* Start never found */
		gap_align->query_start = gap_align->q_start;
		gap_align->subject_start = gap_align->s_start;
	}
	if (found_end == FALSE) {
		gap_align->query_stop = gap_align->q_start + include_query - 1;
		gap_align->subject_stop = gap_align->s_start + include_query - 1;
	}
	q_left = (unsigned char*) GMemFree(q_left);
	s_left = (unsigned char*) GMemFree(s_left);
	gap_align->score = score_right+score_left+middle_score;
	return TRUE;
}

Boolean PerformGappedAlignmentWithTraceback(gab_typ gap_align)
/* Perform a gapped alignment and return the score.  A sap_typ with
	the traceback is also returned.  */
{
	Boolean found_start, found_end;
	Int4	q_length=0, s_length=0, score_right, middle_score, score_left,
		private_q_length, private_s_length, tmp;
	Int4 include_query, index;
	Int4		*tback, *tback1, *p, *q;
	unsigned char	*q_left=NULL, *s_left=NULL;
	unsigned char	*query, *subject, *query_var, *subject_var;

	if (gap_align == NULL) return FALSE;
	found_start = FALSE;
	found_end = FALSE;
	query = gap_align->query;
	subject = gap_align->subject;
#if 0	// afn: 1/22/08 edit to fix core dump
	tback = tback1 = (Int4*) 
		malloc((gap_align->subject_length+gap_align->query_length)*sizeof(Int4));
#else
	MEW(tback,(gap_align->subject_length+gap_align->query_length),Int4);
	tback1 = tback;
#endif
	include_query = gap_align->include_query;

	score_left = 0;
	if (gap_align->q_start != 0 && gap_align->s_start != 0) {
		found_start = TRUE;
#if 0		// afn: 1/22/08 edit to fix core dump
		q_left = (unsigned char *) 
			malloc((gap_align->q_start+2)*sizeof(unsigned char));
		s_left = (unsigned char *) 
			malloc((gap_align->s_start+2)*sizeof(unsigned char));
#else
		MEW(q_left,(gap_align->q_start+2),unsigned char);
		MEW(s_left,(gap_align->s_start+2),unsigned char);
#endif
		q_length = reverse_seq(query, query+gap_align->q_start, q_left);
		s_length = reverse_seq(subject, subject+gap_align->s_start, s_left);
		score_left = SEMI_G_ALIGN(q_left-1, s_left-1, q_length, s_length, 
			tback, &private_q_length, &private_s_length, 
			FALSE, &tback1, gap_align, gap_align->q_start, TRUE);
	        for(p = tback, q = tback1 - 1; p < q; p++, q--)
			{ tmp = *p; *p = *q; *q = tmp; }
		gap_align->query_start = q_length - private_q_length;
		gap_align->subject_start = s_length - private_s_length;
	} else { q_length = gap_align->q_start; s_length = gap_align->s_start; }

	middle_score = 0;
	query_var = query+gap_align->q_start;
	subject_var = subject+gap_align->s_start;
	for (index=0; index<include_query; index++) {
		query_var++;
		subject_var++;
		if (!(gap_align->positionBased))  /*AAS*/
		  middle_score += gap_align->matrix[*query_var][*subject_var];
		else 
		  middle_score += MtrxGScoreGapAlign(gap_align,
			gap_align->q_start+1 + index,*subject_var);
		*tback1 = 0;
		tback1++;
	}
	
	score_right = 0;
	if ((gap_align->q_start+include_query) < gap_align->query_length 
		&& (gap_align->s_start+include_query) < gap_align->subject_length)
	{
		found_end = TRUE;
		score_right = SEMI_G_ALIGN(query+gap_align->q_start+include_query,
			subject+gap_align->s_start+include_query,
			gap_align->query_length-q_length-include_query,
			gap_align->subject_length-s_length-include_query,
			tback1, &private_q_length, &private_s_length,
			FALSE, &tback1, gap_align, gap_align->q_start+include_query,
			FALSE);
		gap_align->query_stop = gap_align->q_start
				+ private_q_length+include_query;
		gap_align->subject_stop = gap_align->s_start
				+ private_s_length+include_query;
	}
	if (found_start == FALSE) {	/* Start never found */
		gap_align->query_start = gap_align->q_start;
		gap_align->subject_start = gap_align->s_start;
	}
	if (found_end == FALSE) {
		gap_align->query_stop = gap_align->q_start + include_query - 1;
		gap_align->subject_stop = gap_align->s_start + include_query - 1;
	}
	gap_align->edit_block = TracebackToGXEB(query, subject, 
		gap_align->query_stop - gap_align->query_start + 1, 
		gap_align->subject_stop-gap_align->subject_start + 1,
		tback, gap_align->query_start, gap_align->subject_start);

	gap_align->edit_block->length1 = gap_align->query_length;
	gap_align->edit_block->length2 = gap_align->subject_length;

        assert(gap_align->edit_block != NULL);
        assert(gap_align->edit_block->esp->num > 0);

	tback = (Int4*) GMemFree(tback);
	q_left = (unsigned char*) GMemFree(q_left);
	s_left = (unsigned char*) GMemFree(s_left);
	gap_align->score = score_right+score_left+middle_score;
	return TRUE;
}

// static Int4 get_current_pos(Int4 * pos, Int4 length)
// Int4 get_current_pos(Int4 * pos, Int4 length)
Int4	get_current_pos(Int4 * pos, Int4 length)
/* Get the current position. */
{
        Int4 val;
        if(*pos < 0) val = -(*pos + length -1); else val = *pos;
        *pos += length;
        return val;
}

/*
	SimpleIntervalToGXEB(Int4 start1, Int4 start2, Int4 length)
	Int4 start1, start2: offset of this interval in sequence 1 and 2.
	Int4 length: length of the interval.
	May be used to produce a gapXEditBlock when there are no subs. or deletions
	in the interval (e.g., ungapped BLAST HSP).
*/
gxeb_typ SimpleIntervalToGXEB(Int4 start1, Int4 start2, Int4 length)
{
	gxeb_typ edit_block;
  	gxes_typ e_script;
	edit_block = GXEBNew(start1, start2);
	edit_block->esp = e_script = GXESNew(NULL);
	e_script->op_type = GAPXALIGN_SUB;
	e_script->num = length;
	return edit_block;
}

sap_typ GXEBToGSeqAlign(gxeb_typ edit_block, e_type sE, e_type qE)
// 	Convert an EditScript chain to a GSeqAlign of type GDenseSeg.
// 	Used for a non-simple interval (i.e., one without subs. or deletions).  
// 	The first GSeqIdPtr in the chain of subject_id and query_id is duplicated for
// 	the GSeqAlign.
{
	gxes_typ	curr, esp;
	Int2		numseg;
	Int4		begin1,begin2,index,start1,start2,length1,length2;
	Int4Ptr		length, start;
	dsp_typ	dsp;
	sap_typ	sap;

	numseg=0;
	start1 = edit_block->start1; start2 = edit_block->start2;
	length1 = edit_block->length1; length2 = edit_block->length2;
	esp = edit_block->esp;
	for(curr=esp; curr; curr=curr->next) { numseg++; }
	start = (Int4*) GMemNew((2*numseg+1)*sizeof(Int4));
	length = (Int4*) GMemNew((numseg+1)*sizeof(Int4));
	for(index=0,curr=esp; curr; curr=curr->next) {
	   switch(curr->op_type) {
		case GAPXALIGN_SUB: case GAPXALIGN_DECLINE:
		  begin1 = get_current_pos(&start1, curr->num);
		  begin2 = get_current_pos(&start2, curr->num);
		  start[2*index] = begin1; start[2*index+1] = begin2;
		  break;
		case GAPXALIGN_DEL:
		  begin1 = -1; begin2 = get_current_pos(&start2, curr->num);
		  start[2*index] = begin1; start[2*index+1] = begin2;
		  break;
		case GAPXALIGN_INS:
		  begin1 = get_current_pos(&start1, curr->num); begin2 = -1;
		  start[2*index] = begin1; start[2*index+1] = begin2;
		  break;
		default: break;
	   } length[index] = curr->num; index++;
	}
	sap = GSeqAlignNew();
	sap->dim =2; // only two dimention alignment.
	dsp = GDenseSegNew(); // Denseg Object is only type for GSeqAlign.
	dsp->dim = 2; dsp->numseg = numseg;
	dsp->subject_id = SeqI(sE); dsp->sE=sE;
	dsp->query_id=0; dsp->qE=qE;
	dsp->starts = start; 
	dsp->lens = length;
	sap->segs = dsp; sap->next = NULL;
	if(IsSplitGSeqAlign(sap)) FixSplitsGSeqAlign(sap); // AFN: fixes gapxdrop bug.
	return sap;
}

void	PutGapAlign(FILE *fp,gab_typ gab)
{
	fprintf(fp,"query: %d..%d; subject %d..%d; score = %d\n",gab->query_start,gab->query_stop,
		gab->subject_start,gab->subject_stop,gab->score);
	fprintf(fp,"...query: %d(%d); subject %d\n",gab->q_start,gab->include_query,gab->s_start);
}


#if 0	//********************** PSG ALIGNMENT ROUTINES *************************
	//************************* WORK IN PROGRESS ****************************
Boolean PerformPSGappedAlignment(gab_typ gap_align, idp_typ *idp)
{
	Boolean		found_start, found_end;
	Int4		q_length=0, s_length=0, middle_score;
	Int4		score_right, score_left, private_q_start, private_s_start;
	Int4		include_query, index;
	unsigned char	*q_left=NULL, *s_left=NULL;
	unsigned char	*query, *subject, *query_var, *subject_var;

	if (gap_align == NULL) return FALSE;
	found_start = FALSE; found_end = FALSE;

	query = gap_align->query;
	subject = gap_align->subject;
	include_query = gap_align->include_query;
	score_left = 0;
	if (gap_align->q_start != 0 && gap_align->s_start != 0) {
		found_start = TRUE;
		q_left = (unsigned char *) malloc((gap_align->q_start+2)*sizeof(unsigned char));
		s_left = (unsigned char *) malloc((gap_align->s_start+2)*sizeof(unsigned char));

		q_length = reverse_seq(query, query+gap_align->q_start, q_left);
		s_length = reverse_seq(subject, subject+gap_align->s_start, s_left);

		score_left = SEMI_G_PSGALIGN(q_left-1, s_left-1, q_length, s_length, 
			NULL, &private_q_start, &private_s_start, TRUE, NULL, 
			gap_align, gap_align->q_start, TRUE,idp);
		gap_align->query_start = q_length - private_q_start;
		gap_align->subject_start = s_length - private_s_start;
	} else {
		q_length = gap_align->q_start;
		s_length = gap_align->s_start;
	}
	middle_score = 0;
	query_var = query+gap_align->q_start;
	subject_var = subject+gap_align->s_start;
	for (index=0; index<include_query; index++) {
		query_var++;
		subject_var++;
		if (!(gap_align->positionBased))  /*AAS*/
		  middle_score += gap_align->matrix[*query_var][*subject_var];
		else 
		  middle_score += MtrxGScoreGapAlign(gap_align,
				gap_align->q_start+1 + index,*subject_var);
	}
	score_right = 0;
	if (gap_align->q_start+include_query < gap_align->query_length 
		&& gap_align->s_start+include_query < gap_align->subject_length)
	{
		found_end = TRUE;
		score_right = SEMI_G_PSGALIGN(query+gap_align->q_start+include_query,
			subject+gap_align->s_start+include_query,
			gap_align->query_length-q_length-include_query,
			gap_align->subject_length-s_length-include_query,
			NULL, &(gap_align->query_stop),
			&(gap_align->subject_stop), TRUE,
			NULL, gap_align, gap_align->q_start+include_query, FALSE,idp);
		gap_align->query_stop += gap_align->q_start+include_query;
		gap_align->subject_stop += gap_align->s_start+include_query;
	}

	if (found_start == FALSE) {	/* Start never found */
		gap_align->query_start = gap_align->q_start;
		gap_align->subject_start = gap_align->s_start;
	}
	if (found_end == FALSE) {
		gap_align->query_stop = gap_align->q_start + include_query - 1;
		gap_align->subject_stop = gap_align->s_start + include_query - 1;
	}
	q_left = (unsigned char*) GMemFree(q_left);
	s_left = (unsigned char*) GMemFree(s_left);
	gap_align->score = score_right+score_left+middle_score;
	return TRUE;
}

Boolean PerformPSGappedAlignmentWithTraceback(gab_typ gap_align, idp_typ *idp)
/* Perform a gapped alignment and return the score.  A sap_typ with
	the traceback is also returned.  */
{
	Boolean found_start, found_end;
	Int4	q_length=0, s_length=0, score_right, middle_score, score_left,
		private_q_length, private_s_length, tmp;
	Int4 include_query, index;
	Int4		*tback, *tback1, *p, *q;
	unsigned char	*q_left=NULL, *s_left=NULL;
	unsigned char	*query, *subject, *query_var, *subject_var;

	if (gap_align == NULL) return FALSE;
	found_start = FALSE;
	found_end = FALSE;
	query = gap_align->query;
	subject = gap_align->subject;
	tback = tback1 = (Int4*) 
		malloc((gap_align->subject_length+gap_align->query_length)*sizeof(Int4));
	include_query = gap_align->include_query;

	score_left = 0;
	if (gap_align->q_start != 0 && gap_align->s_start != 0) {
		found_start = TRUE;
		q_left = (unsigned char *) 
			malloc((gap_align->q_start+2)*sizeof(unsigned char));
		s_left = (unsigned char *) 
			malloc((gap_align->s_start+2)*sizeof(unsigned char));
		q_length = reverse_seq(query, query+gap_align->q_start, q_left);
		s_length = reverse_seq(subject, subject+gap_align->s_start, s_left);
		score_left = SEMI_G_PSGALIGN(q_left-1, s_left-1, q_length, s_length, 
			tback, &private_q_length, &private_s_length, 
			FALSE, &tback1, gap_align, gap_align->q_start, TRUE,idp);
	        for(p = tback, q = tback1 - 1; p < q; p++, q--)
			{ tmp = *p; *p = *q; *q = tmp; }
		gap_align->query_start = q_length - private_q_length;
		gap_align->subject_start = s_length - private_s_length;
	} else { q_length = gap_align->q_start; s_length = gap_align->s_start; }

	middle_score = 0;
	query_var = query+gap_align->q_start;
	subject_var = subject+gap_align->s_start;
	for (index=0; index<include_query; index++) {
		query_var++;
		subject_var++;
		if (!(gap_align->positionBased))  /*AAS*/
		  middle_score += gap_align->matrix[*query_var][*subject_var];
		else 
		  middle_score += MtrxGScoreGapAlign(gap_align,
			gap_align->q_start+1 + index,*subject_var);
		*tback1 = 0;
		tback1++;
	}
	
	score_right = 0;
	if ((gap_align->q_start+include_query) < gap_align->query_length 
		&& (gap_align->s_start+include_query) < gap_align->subject_length)
	{
		found_end = TRUE;
		score_right = SEMI_G_PSGALIGN(query+gap_align->q_start+include_query,
			subject+gap_align->s_start+include_query,
			gap_align->query_length-q_length-include_query,
			gap_align->subject_length-s_length-include_query,
			tback1, &private_q_length, &private_s_length,
			FALSE, &tback1, gap_align, gap_align->q_start+include_query,
			FALSE,idp);
		gap_align->query_stop = gap_align->q_start
				+ private_q_length+include_query;
		gap_align->subject_stop = gap_align->s_start
				+ private_s_length+include_query;
	}
	if (found_start == FALSE) {	/* Start never found */
		gap_align->query_start = gap_align->q_start;
		gap_align->subject_start = gap_align->s_start;
	}
	if (found_end == FALSE) {
		gap_align->query_stop = gap_align->q_start + include_query - 1;
		gap_align->subject_stop = gap_align->s_start + include_query - 1;
	}
	gap_align->edit_block = TracebackToGXEB(query, subject, 
		gap_align->query_stop - gap_align->query_start + 1, 
		gap_align->subject_stop-gap_align->subject_start + 1,
		tback, gap_align->query_start, gap_align->subject_start);

	gap_align->edit_block->length1 = gap_align->query_length;
	gap_align->edit_block->length2 = gap_align->subject_length;

        assert(gap_align->edit_block != NULL);
        assert(gap_align->edit_block->esp->num > 0);

	tback = (Int4*) GMemFree(tback);
	q_left = (unsigned char*) GMemFree(q_left);
	s_left = (unsigned char*) GMemFree(s_left);
	gap_align->score = score_right+score_left+middle_score;
	return TRUE;
}

Int4 PSGALIGN(unsigned char * A, unsigned char * B, Int4 M, Int4 N,
		Int4 * S, Int4 * pei, Int4 * pej, Int4 **sapp, 
		gab_typ gap_align, Int4 query_offset, Boolean reversed,
		idp_typ *idp)
{ 
  data_t data;
  Int4 i, j, cb,  j_r, s, k;
  unsigned char st, std, ste;
  Int4 gap_open, gap_extend, decline_penalty;
  register Int4 c, d, e, gap_ox,t, tt, f, tt_start;
  Int4 best_score = 0;
  register Int4 *wa;
  register dp_ptr dp, dyn_prog;
  unsigned char **state, *stp, *tmp;
  unsigned char *state_array;
  Int4	**matrix;
  register Int4 X;
  gxdsa_typ state_struct;
#if 1
  Int4	del_ox, ins_ox, del_open, ins_open, del_extend, ins_extend;
  Int4	seq_inso, seq_insx;
#endif
  
  matrix = gap_align->matrix;
  *pei = *pej = 0;
  data.sapp = *sapp = S;
  data.last= 0;
  gap_ox = gap_align->gap_open + gap_align->gap_extend;
  decline_penalty = gap_align->decline_align;
  gap_open = gap_align->gap_open;
  gap_extend = gap_align->gap_extend;
  X = gap_align->x_parameter;

  if (X < gap_ox) X = gap_ox;
  if(N <= 0 || M <= 0) { *pei = *pej; return 0; }
  GapXDropPurgeState(gap_align->state_struct);

  j = (N + 2) * sizeof(dp_node);
  dyn_prog = (dp_ptr)malloc(j);

// Initialize dyn_prog matrix...
  state = (unsigned char**) malloc(sizeof(unsigned char *)*(M+1));
  dyn_prog[0].CC = 0;
  dyn_prog[0].FF = -gap_ox;
  c = dyn_prog[0].DD = -gap_ox;

  /* Protection against divide by zero. */
  if (gap_extend > 0)
  	state_struct = GapXDropGetState(&gap_align->state_struct, X/gap_extend+5);
  else state_struct = GapXDropGetState(&gap_align->state_struct, N+3);

  state_array = state_struct->state_array;
  state[0] = state_array;
  stp  = state[0];
  for(i = 1; i <= N; i++) {
    if(c < -X) break;
    dyn_prog[i].CC = c;
    dyn_prog[i].DD = c - gap_ox; 
    dyn_prog[i].FF = c - gap_ox;
    c -= gap_extend;
    stp[i] = 1;
  }
  state_struct->used = i+1;
  
  tt = 0;  j = i;
  for(j_r = 1; j_r <= M; j_r++) {
     /* Protection against divide by zero. */
    if (gap_extend > 0)
    	state_struct = GapXDropGetState(&gap_align->state_struct, j-tt+5+X/gap_extend);
    else
	state_struct = GapXDropGetState(&gap_align->state_struct, N-tt+3);
    state_array = state_struct->state_array + state_struct->used + 1;
    state[j_r] = state_array - tt + 1;
    stp = state[j_r];
    tt_start = tt; 
    if (!(gap_align->positionBased)) /*AAS*/
      wa = matrix[A[j_r]]; 
    else {
      if(reversed) wa = gap_align->posMatrix[M - j_r];
      else wa = gap_align->posMatrix[j_r + query_offset];
    } e = c = f= MININT;
    
    for (cb = i = tt, dp = &dyn_prog[i-1]; i < j; i++) {
	d = (++dp)->DD;
	if (d < f) { d=f; std=83;} else std=42;
	if (e < f) {e=f; ste = 83;} else ste = 41;
	if (c < d || c < e) {
	    if (d < e) { c = e; st = ste; }	// AFN: extension???
	    else { c = d; st = std; }		// AFN: deletion???
	    if (best_score - c > X) {
		c = dp->CC+wa[B[i+1]]; f = dp->FF;		// AFN: match???
		if (tt == i) tt++; else { dp->CC = MININT; }
	    } else {
		cb = i;
		e -= gap_extend; dp->DD = d-gap_extend;
		d = dp->FF; dp->FF = f-decline_penalty; f = d;
		d = dp->CC+wa[B[i+1]]; dp->CC = c; c=d;
	    } st += 5;
	} else {
	    st = 0;
	    if (best_score - c > X){
		c = dp->CC+wa[B[i+1]]; f= dp->FF;
		if (tt == i) tt++; else { dp->CC =  MININT; }
	    } else {
		cb = i;
		if (c > best_score) { best_score = c; *pei = j_r; *pej = i; }
		if ((c-=gap_ox) > (d-=gap_extend)) { dp->DD = c; } 
		else { dp->DD = d; if (std == 83) st+=60; else st+=30; } 
		if (c > (e-=gap_extend)) { e = c; }
		else { if (ste> 50) st+=20; else st +=10; }
		c+=gap_ox; d = dp->FF;
		if (f < c-gap_open) { dp->FF = c-gap_open-decline_penalty; }
		else { dp->FF = f-decline_penalty; st+= 5; }
		f = d;
		d = dp->CC+wa[B[i+1]]; dp->CC = c; c = d;      
	    }
	}
	stp[i] = st;
    }
    if (tt == j) { j_r++; break;}
    if (cb < j-1) { j = cb+1;}
    else {
	while (e >= best_score-X && j <= N) {
	    dyn_prog[j].CC = e; 
	    dyn_prog[j].DD = e-gap_ox; dyn_prog[j].FF = e-gap_open-decline_penalty;
	    e -= gap_extend; stp[j] = 1;
	    j++;
	}
    }
    if (j <= N) {
	dyn_prog[j].DD = dyn_prog[j].CC= dyn_prog[j].FF = MININT;
	j++; 
    }
    state_struct->used += (MAX(i, j) - tt_start + 1);
  }
  i = *pei; j = *pej;
  tmp = (unsigned char*) malloc(i+j);
  for (s=0, c = 0; i> 0 || j > 0; c++) {
      t = state[i][j]; k  = t %5;
      if (s == 1) 
	  if ((t/10)%3 == 1) k = 1; else if ((t/10)%3 == 2) k = 3;
      if (s == 2)
	  if ((t/30) == 1) k = 2; else if ((t/30) == 2) k = 3;
      if (s == 3 && ((t/5)%2) == 1) k = 3;
      if (k == 1) { j--;} else if (k == 2) {i--;} else {j--; i--;}
      tmp[c] = s = k;
  }
  c--;
  while (c >= 0) {
      if (tmp[c] == 0) REP_
      else if (tmp[c] == 1) INS_(1)
      else if (tmp[c] == 3) REPP_			  
      else DEL_(1)     
      c--;
  }
  GMemFree(tmp); GMemFree(state); GMemFree(dyn_prog);
  *sapp = data.sapp;
  return best_score;
}

Int4 SEMI_G_PSGALIGN(unsigned char * A, unsigned char * B, Int4 M, Int4 N,
		Int4 * S, Int4 * pei, Int4 * pej, 
		Boolean score_only, Int4 **sapp, gab_typ gap_align,
		Int4 query_offset, Boolean reversed,idp_typ *idp)
		
{ 
  dp_ptr dyn_prog;
  Int4 i, j, cb, j_r, gap_open, decline_penalty;
  register Int4 c, d, e, gap_ox, tt, gap_extend, X, f;
  Int4 best_score = 0;
  Int4 * *matrix;
  register Int4 * wa;
  register dp_ptr dp;
  
  if(!score_only)
    return PSGALIGN(A,B,M,N,S,pei,pej,sapp,gap_align,query_offset,reversed,idp);
  
  matrix = gap_align->matrix;
  *pei = *pej = 0;
  gap_ox = (gap_open=gap_align->gap_open) + gap_align->gap_extend;
  gap_extend = gap_align->gap_extend;
  decline_penalty = gap_align->decline_align;
  X = gap_align->x_parameter;

  if (X < gap_ox) X = gap_ox;
  if(N <= 0 || M <= 0) return 0;

  j = (N + 2) * sizeof(dp_node);
  dyn_prog = (dp_ptr)malloc(j);

  dyn_prog[0].CC = 0; c = dyn_prog[0].DD = -gap_ox;
  dyn_prog[0].FF = -gap_ox;
  for(i = 1; i <= N; i++) {
    if(c < -X) break;
    dyn_prog[i].CC = c;
    dyn_prog[i].DD = c - gap_ox; 
    dyn_prog[i].FF = c-gap_ox;
    c -= gap_extend;
  }

  tt = 0;  j = i;
  for (j_r = 1; j_r <= M; j_r++) {
      if (!(gap_align->positionBased)) /*AAS*/
	  wa = matrix[A[j_r]]; 
      else {
	  if(reversed) wa = gap_align->posMatrix[M - j_r];
	  else wa = gap_align->posMatrix[j_r + query_offset];
      }
      e = c =f = MININT;
      for (cb = i = tt, dp = &dyn_prog[i]; i < j; i++) {
	  d = dp->DD;
	  if (e < f) e = f;
	  if (d < f) d = f;
	  if (c < d || c < e) {
	      if (d < e) { c = e; } else { c = d; }
	      if (best_score - c > X) {
		  c = dp->CC+wa[B[i+1]]; f = dp->FF;
		  if (tt == i) tt++; else { dp->CC = MININT;}
	      } else {
		  cb = i; e -= gap_extend; dp->DD =  d- gap_extend;
		  d = dp->CC+wa[B[i+1]]; dp->CC = c; c=d;
		  d = dp->FF; dp->FF = f-decline_penalty; f = d;
	      }
	  } else {
	      if (best_score - c > X){
		  c = dp->CC+wa[B[i+1]]; f= dp->FF;
		  if (tt == i) tt++; else { dp->CC = MININT;}
	      } else {
		  cb = i; 
		  if (c > best_score) { best_score = c; *pei = j_r; *pej = i; } 
		  if ((c-=gap_ox) > (d-=gap_extend)) { dp->DD = c; } else { dp->DD = d; } 
		  if (c > (e-=gap_extend)) { e = c; } 
		  c+=gap_ox;
		  d = dp->FF;
		  if (c-gap_open>f) dp->FF = c-gap_open-decline_penalty;
		  else dp->FF = f-decline_penalty;
		  f = d;
		  d = dp->CC+wa[B[i+1]]; dp->CC = c; c = d;
	      }
	  }
	  dp++;
      }
      if (tt == j) break;
      if (cb < j-1) { j = cb+1;}
      else while (e >= best_score-X && j <= N) {
	  dyn_prog[j].CC = e; dyn_prog[j].DD = e-gap_ox;
	  dyn_prog[j].FF = e-gap_open-decline_penalty;
	  e -= gap_extend; j++;
      }
      if (j <= N) {
	  dyn_prog[j].DD = dyn_prog[j].CC = dyn_prog[j].FF = MININT; j++;
      }
  }
  GMemFree(dyn_prog);
  return best_score;
}


#endif


