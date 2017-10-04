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

/************************************************************************/
/***  pseg.h                                                          ***/
/*** #include precomputed ln(fact(n)) array in "lnfac.h"              ***/
/************************************************************************/
#if !defined(PSEG)
#define PSEG
#include "stdinc.h"
#include "afnio.h"
#include "alphabet.h"
#include "sequence.h"
#include <fcntl.h>

#define LENCORRLIM 120
#define ALPHA_SIZE 20
#define ALPHA_SIZE_PLUS 21

struct Sequence
  {
   char *header;
   char *seq;
   int start;                     /* for windows */
   int length;
   struct Sequence *parent;       /* for windows */
   struct Sequence *root;         /* for subwindows */
   int *state;
   double entropy;
   int *composition;
};

struct Segment {
   int begin;
   int end;
   struct Segment *next;
};

/******************************  private **********************************/
struct Segment *per_mergesegs(struct Sequence *seq, struct Segment **persegs);
struct Sequence *seq_phase(struct Sequence *seq, int phase, int period);
void    per_seqprep(struct Sequence *seq, struct Segment **persegs);
int	findlo(int i, int limit, double *H);
int	findhi(int i, int limit, double *H);
void	trim(struct Sequence *seq, int *leftend, int *rightend);
double	lnperm(int *sv, int tot);
double	lnass(register int *sv);
void    appendseg(struct Segment *segs, struct Segment *seg);
void    freesegs(struct Segment *segs);

void    segseq(struct Sequence *seq, struct Segment **segs,int  offset);
double *seqent(struct Sequence *seq);
double *seqent(struct Sequence *seq);
void    mergesegs(struct Sequence *seq, struct Segment *segs);
void    singreport(struct Sequence *seq, struct Segment *segs);
/******************************  private **********************************/
void    genwininit(void);
void	closeseq(struct Sequence *seq);
struct Sequence *openwin(struct Sequence *parent, int start,
        int length);
BooLean	shiftwin1(struct Sequence *win);
void	closewin(struct Sequence *win);
void	compon(struct Sequence *win);
static int state_cmp(int *s1,int *s2);
void	stateon(struct Sequence  *win);
void	enton(struct Sequence *win);
void	entropy_init(int window);
double	entropy(register int *sv);
void	decrementsv(register int *sv, register int Class);
void	incrementsv(register int *sv, register int Class);
void	upper(register char *string, size_t len);
BooLean	process_pseg(struct Segment *segs, e_type E);
void    per_process_pseg(struct Sequence *seq, e_type E);
struct Sequence *CreatePSegSeq(e_type E, a_type A);
/****************************** PUBLIC ***********************************/
BooLean	ProcessSeqPSeg(int window, double lowcut, double highcut,
        int maxtrim, e_type E, a_type A);

#endif

