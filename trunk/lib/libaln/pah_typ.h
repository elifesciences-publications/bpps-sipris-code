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

/* psaln_heap.h - scan hit min/max heap. */
#if !defined (PSALN_HEAP)
#define PSALN_HEAP
#include "afnio.h"
#include "alphabet.h"
#include "histogram.h"
#include "sequence.h"
#include "mheap.h"
#include "set_typ.h"

/************************** Scan Heap Type ***************************

/*********************************************************************/
typedef struct psaln_info_link {
        e_type                  E;      /** database sequence **/
	char			*operation; // operation string
	Int4			Score;	// alignment score
	Int4			start;	// alignment start
        unsigned char           rpts;   /** number of rpts **/
} psaln_info_type;

class pah_typ {
public:
		pah_typ(Int4 );
		~pah_typ( ){ Free(); }
	// Int4	Insert(pai_typ I);
	Int4	Insert(e_type , char *, Int4 , unsigned char ,Int4 );
	Int4	Insert(e_type , char *, Int4 , Int4 );
	Int4    DelMin(e_type *, char **, unsigned char *,Int4 *,Int4*);
	Int4    DelMax(e_type *, char **, unsigned char *,Int4 *,Int4*);
	void	Put(FILE *);
private:
	void	Free();		// free memory.
	void	NilInfo(Int4);
	void	MkInfo(Int4 ,e_type , char *, Int4 , unsigned char ,Int4);
	psaln_info_type *info;	// information on hits.
	mh_type	mH;		/** heap for best hits **/
	Int4    hpsz;		/** heap size **/
	h_type  HG;		/** histogram for key scores **/
};

#endif

