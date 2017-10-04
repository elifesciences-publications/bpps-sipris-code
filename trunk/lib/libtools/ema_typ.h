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

/* ema_typ.h - codes and constants for edit cma_typ. */
#if !defined (_EMA_TYP_)
#define _EMA_TYP_
#include "sequence.h"

typedef struct editcma_type {
	Int4		fsq,ssq;	// fullseq, subseq.
	BooLean		new_seq;	// new subseq between ssq && ssq + 1.
	Int4		start,end;	// start and end of potential repeat.
	Int4		score;		// score for detection.
	double		Evalue;		// 
	struct editcma_type *next;
} editcma_type_type;
typedef editcma_type_type       *ema_typ;

/********************************* private ********************************/
ema_typ initEMA(Int4 fsq, Int4 ssq, BooLean new_seq, Int4 start,
        Int4 end, Int4 score, double Evalue);

/********************************* PUBLIC ********************************/
ema_typ MakeEMA(Int4 fsq, Int4 ssq, Int4 start, Int4 end);
ema_typ MakeEMA(Int4 fsq, Int4 ssq, Int4 start, Int4 end, Int4 score,
        double Evalue);
ema_typ AppendEMA(ema_typ ema, ema_typ ema2);
void    PutEMA(FILE *fp, ema_typ ema);
ema_typ	PurgeEMA(ema_typ ema);
Int4	LengthEMA(ema_typ ema);
void    NilEMA(ema_typ ema);

/********************************* MACROS ********************************/
#define	NextEMA(ema)	((ema)->next)
#define	FullSeqEMA(ema)	((ema)->fsq)
#define	SubSeqEMA(ema)	((ema)->ssq)
#define	ScoreEMA(ema)	((ema)->score)
#define	EvalueEMA(ema)	((ema)->Evalue)
#define	StartEMA(ema)	((ema)->start)
#define	EndEMA(ema)	((ema)->end)
#define	NewEMA(ema)	((ema)->new_seq)

#endif

