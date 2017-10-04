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

/* fsm.h - finite state machine program. */
#if !defined (_FSM_)
#define _FSM_
#include "stdinc.h"
#include "afnio.h"
#include "dheap.h"
#include "mheap.h"
#include "mlist.h"

/*********************** Finite State Machine ************************
FSM = (Q,q0,A,Sigma,d)
	Q = States
	q0 e Q = start state
	A subset Q = accepting states
	Sigma = input alphabet
	d = function from Q x Sigma -> Q (transition function)

if in an accepting state then execute appropriate action.
	- go to positions on list and try extending alignment

input text T[1..n];  pattern to be found P[1..m]
 n = input string length; m = pattern length 

	(note: lex is a FSM)

  input tokens = A-Y and '-' 'x' '*' '\0'

 if q = a then  go to query sequence:
   pos[q][1..n][A...Y] = list of positions matching pattern in accepting 
	state = NULL if not an accepting state.

  blast method:
	1. compile list of high scoring words and make fsm.
	2. scan database for hits.
	3. extend hits.
(for purge extend only until find that score >= cutoff.)

	QWL
	..:  -> S = R(Q,N) + R(W,Y) + R(L,L).
	NYL	
		if(S > Threshold) then extend hit to find MSP.

		need drop-score.
 *********************************************************************/

/*************************** generic fsm type **************************/
typedef struct {
	Int4	nQ;		/** number of States **/
	Int4	**d;		/** d[q][r] = function from Q x A -> Q **/
	ml_type	mlist;		/** lists for accept **/
	Int4	*tmp;		/** temporary list **/
	Int4	T;
	/******************************************************/
	Int4	nhits;		/** number of hits **/
	/******************* New: matrix **********************/
	Int4	**matrix;	/** score[pos][res] ***/
	Int4	nlet;		/** number of letters in alphabet **/
	Int4	len_mtrx;	/** leng of matrix **/
} fsm_type;
typedef fsm_type *fsm_typ;
/*********************************************************************/
#define MAX_SEQ_LENG_GB 30000
#define MAX_ID_LENG_GB  10000
#define MAX_OVERLAP_GB	9 
#define MAX_FRACTION_GB	7 

/******************************* private *******************************/
void	fsm_error(char *s);
/****************************** macros ********************************/
#define StateGB(q,t)		(((t)<<9)|(q))
#define MaxStateGB(n)		(((n)<<9)|0777)

/******************************* PUBLIC *******************************/
fsm_typ	MakeFSM(Int4 T, Int4 numlet, Int4 leng, Int4 **matrix);
Int4    MatcherFSM(Int4 len, register unsigned char *seq, register fsm_typ B);
void    NilFSM(fsm_typ B);

Int4    *Scan4MatchesFSM(register Int4 len, register unsigned char *seq,
        register Int4 i, Int4 *pos, Int4 *number, register fsm_typ B);
/******************************** MACROS *****************************/



#endif

