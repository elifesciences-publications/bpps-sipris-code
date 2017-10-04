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

#if !defined(_MEMORY_)
#define _MEMORY_
#include "stdinc.h"
#include "set_typ.h"

/*************************** Memory Type ******************************

action 356 (p=0.472):
S: .....
L: ........*****************

'.' = failed; '*' = succeeded.


**********************************************************************/
class mem_typ {
public:
		mem_typ(Int4,char);
		mem_typ(Int4,char,unsigned short,unsigned short,double);
		~mem_typ( ){ Free(); }
	void	Clear(Int4);	// Set all events to failures.
	void    Failure(unsigned short);
	void    Success(unsigned short);
	void	Put(FILE *fptr, unsigned short);
	double  Chance(unsigned short);
	BooLean	Empty(Int4);
	BooLean	Full(Int4);
private:
	void	init(unsigned short ,unsigned short , Int4, double ,char);
	void	Free();
	char	ActionType;
	set_typ		*longterm;	// long term memory; one for each action.
	unsigned short	longLen;	// length of long term memory
	unsigned short	*eventL;	// eventL[i] number of successful long term actions.
	set_typ		*shorterm;	// short term memory.
	unsigned short	shortLen;	// length of short term memory
	unsigned short	*eventS;	// N short term events.
	double		Ps;		// pseudo counts
	Int4		N;		// number of actions that can be performed.
};


#endif

