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

#if !defined(_LKP_TYP_)
#define _LKP_TYP_

#include "shm_typ.h"

class lkp_typ {
public:
	lkp_typ(shm_typ **Shm, Int4 Data_len);
	lkp_typ(Int4 Data_len,Int4 Total,Int4 Max_cnt,Int4 *Start,Int4 *Look);
	lkp_typ(FILE *fp);
	~lkp_typ() { Free(); }
	Int4	 	*GetLook() { return look; }
	Int4   		GetDataLen() { return data_len; }
	Int4 		GetMaxCnt() { return max_cnt; }
	Int4 		GetTotal() { return total; }
	Int4		*GetStart() { return start; }
	void		Write(FILE *fp);
private:
	Int4 		*start;
	Int4 		data_len,max_cnt,total;
	Int4		*look;
	void		Free();
	
};
#endif


