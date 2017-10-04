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

#ifndef __GPXDROP__
#define __GPXDROP__

// Adapted from gapxdrop.h by Gennadiy Savchuk, Jinqhui Zhang, Tom Madden
// which they adapted from Webb Miller.
/* $Revision: 6.5 $ * gapxdrop.h,v $ * Revision 1.1  1996/12/12  14:02:51  madden * */

#include "stdinc.h"
#include "my_ncbimath.h"
#include "sequence.h"

#define CODON_LENGTH 3

#define GAPXALIGN_SUB ((unsigned char)0)  // op types within the edit script
#define GAPXALIGN_INS ((unsigned char)1)
#define GAPXALIGN_DEL ((unsigned char)2)
#define GAPXALIGN_DECLINE ((unsigned char)3)

typedef struct gap_x_edit_script {
        unsigned char op_type;  // GAPXALIGN_SUB, GAPXALIGN_INS, or GAPXALIGN_DEL.
        Int4 num;       	// Number of operations.
        struct gap_x_edit_script *next;
} gxes_type, *gxes_typ;

typedef struct gap_x_edit_block {
	Int4	start1,		// starts of alignments.
		start2,	
		length1,	// total lengths of the sequences.
		length2;
	gxes_typ esp;
} gxeb_type, * gxeb_typ;

// Structure to keep memory for state structure.
typedef struct gxd_state_array_struct {
	Int4 	length,		/* length of the state_array. */
		used;		/* how much of length is used. */
	unsigned char	*state_array;	/* array to be used. */
	struct gxd_state_array_struct * next;
} gxdsa_type, *gxdsa_typ;

#define  MtrxGScoreGapAlign(S,x,y)      ((S)->posMatrix[(x)][(y)])

// Structure used to pass arguments for gapped alignment functions.
typedef struct _afn_gapalign_blk {
	// The following must be supplied by caller. 
	unsigned char *query,	/* The query sequence. */
		 *subject;	/* The subject sequence. */
	Int4	query_length,	/* the length of the query. */
		subject_length,	/* The subject length. */
		q_start,	/* query letter to start the gapped align. */
		s_start,	/* subject letter to start the gapped align.*/ 
		include_query,	/* length of query (starting from q_start) that 
				   MUST be included in an alignment. */
		gap_open,	/* Cost to open a gap. */
		gap_extend,	/* cost to extend a gap. */
		decline_align,  /* decline to align penalty */
		x_parameter;	/* values of X-dropoff parameter. */	
	Int4	**matrix;	/* Matrix for the alignment. */
	Int4	**posMatrix;	/* Matrix for position-based searches. */

	// The state, state_column_length, and state_row_length are used by ALIGN.
	// If state is NULL, then ALIGN allocates a memory block and frees it 
	// before returning (for "call and forget" applications).
	gxdsa_typ state_struct;
	// The score, start, and stop of alignments are filled in by the
	// functions PerformGappedAlignment and PerformGappedAlignmentWithTraceback
	Int4	score, 		/* score of alignment. */
		query_start,	/* start of alignment on query. */
		query_stop,	/* end of alignment on query. */
		subject_start,	/* start of alignment on subject. */
		subject_stop;	/* end of alignment on subject. */
		// gxeb_typ filled in by PerformGappedAlignmentWithTraceback,
		// used to make a sap_typ. 
	gxeb_typ edit_block;	//  Another TLM kludge for display. 
        Int4	*tback;
	Boolean	positionBased;		// Is the search position based?
	Boolean posConverged;
} gab_type, * gab_typ;

gxdsa_typ GapXDropStateDestroy(gxdsa_typ state_struct);
gxeb_typ SimpleIntervalToGXEB(Int4 start1, Int4 start2, Int4 length);
gxeb_typ GXEBNew (Int4 start1, Int4 start2);
gxeb_typ GXEBDelete (gxeb_typ edit_block);

Boolean PerformGappedAlignment (gab_typ);
Boolean PerformGappedAlignmentWithTraceback (gab_typ);

gxeb_typ TracebackToGXEB (unsigned char *A, unsigned char *B,
	Int4 M, Int4 N, Int4 * S, Int4 start1, Int4 start2);
Int4 	AlnEndFlagGXEB( );

/*************************************************************************
        Allocates gab_typ and "state", if state_column_length and 
        state_row_length are not NULL. 
        For "call and forget" applications, state_column_length and
        state_row_length should both be set to zero.  ALIGN will
        then allocate and deallocate this memory.
 *************************************************************************/
gab_typ GABNew(Int4 open, Int4 extend);
gab_typ GABNew(e_type qE, Int4 open, Int4 extend, Int4 **matrix, Int4 decl_aln,
        double x_parameter, Int4 min_query_align);
void	SetSubjectGAB(Int4 q_start, Int4 s_start, e_type sE, gab_typ gab);
void    SetSubjectGAB(Int4 q_start, Int4 s_start, unsigned char *pssq, 
                        Int4 len_ssq, gab_typ gab);
void    PutGapAlign(FILE *fp,gab_typ gab);
Int4    ScoreGAB(gab_typ gap_align);
sap_typ ExtractAlnGAB(gab_typ gap_align, double e_value, e_type sE, e_type qE,
        double PerNat);

// Destruct Function for gap_type.  If "state" is not NULL, then it's deallocated.
gab_typ GABDelete(gab_typ gap_align);

Int4 SEMI_G_ALIGN(unsigned char *A, unsigned char *B, Int4 M, Int4 N,
                Int4 * S, Int4 * pei, Int4 * pej,
                Boolean score_only, Int4 **sapp, gab_typ gap_align,
                Int4 query_offset, Boolean reversed);

Int4 ALIGN(unsigned char *A, unsigned char *B, Int4 M, Int4 N,
                Int4 * S, Int4 * pei, Int4 * pej, Int4 **sapp,
                gab_typ gap_align, Int4 query_offset, Boolean reversed);

sap_typ GXEBToGSeqAlign(gxeb_typ edit_block, e_type sE, e_type qE);

//********************** PSG ALIGNMENT ROUTINES *************************
//************************* WORK IN PROGRESS ****************************
#include "idp_typ.h"

Boolean PerformPSGappedAlignment(gab_typ gap_align, idp_typ *idp);
Boolean PerformPSGappedAlignmentWithTraceback(gab_typ gap_align, idp_typ *idp);
Int4 PSGALIGN(unsigned char * A, unsigned char * B, Int4 M, Int4 N,
                Int4 * S, Int4 * pei, Int4 * pej, Int4 **sapp,
                gab_typ gap_align, Int4 query_offset, Boolean reversed,
                idp_typ *idp);
Int4 SEMI_G_PSGALIGN(unsigned char * A, unsigned char * B, Int4 M, Int4 N,
                Int4 * S, Int4 * pei, Int4 * pej,
                Boolean score_only, Int4 **sapp, gab_typ gap_align,
                Int4 query_offset, Boolean reversed,idp_typ *idp);

//**** WARNING THESE ARE PRIVATE ROUTINES THAT SHOULD BE STATIC ******
gxdsa_typ GapXDropGetState(gxdsa_typ *head, Int4 length);
Boolean GapXDropPurgeState(gxdsa_typ state_struct);
Int4 reverse_seq(unsigned char *seq,unsigned char *pos,unsigned char *target);
Int4	get_current_pos(Int4 * pos, Int4 length);

#endif /* !__GPXDROP__ */

