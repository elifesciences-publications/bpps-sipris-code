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

#if !defined(_GHM_TYP_)
#define _GHM_TYP_

#if 0
	labels:	i O o M;
	other States: S B E T
	Dimensions (for dynamic programming): ?? 4...

	Applications: 
		TMHMM.
		coiled coils.
		beta barrel repeats.
		other repeats.
#endif

class ghs_typ {		// general HMM state type.
public:
	
private:
	Int4	id;		// identifier integer for this state.
	Int4    *out;           // transition in list for each state.
	Int4    *in;            // transition out list for each state.
	Int4	Nin,Nout;	// 
};

class ghm_typ { 	// general HMM type (for an arbitrary HMM).
public:
	ghm_typ();		//
private:
	// States
	sta_typ	*state;
	Int4	num_states;	// 
	// Transition probabilities
	float	**trns_prb;	// transition probabilities for states.
	Int4	*out;		// transition in list for each state.
	Int4	*in;		// transition out list for each state.
	char	**trace;	// trace back list for each state.
	char	*label;		// label for path; max_num = 255.
	float	**emit_prb;	// emission probabilities.
	Int4	max_num_trans;	// maximum number of transitions.
};

#endif

