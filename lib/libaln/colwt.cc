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

#include "colwt.h"

double	mv_col_weight(Int4 c, Int4 w, Int4 n, Int4 lemon, BooLean *null)
/********************************************************************

   If a fmodel has c on columns and a length of w.  Then there are 
   w-2 choose c-2  configurations for that model.

          n     1 lemon       w_old
          |     |  |           |
        ........*??*???????????*.............
                 |--- (w-2) --|
               .....
              -10234
      n = -inf..0 or +2..+inf.	 (i.e., n != 1)

       lemon = 1..w and n = -f..end.

 ********************************************************************/
{
	Int4 start=1, end=w, w_new,w_old;
#if 0
	double	b0,b2,P;
	Int4	i,s,e;
#endif

	w_old = w;
	if(lemon == 1) {     /** case 1: **/
	   if(n < 1){          /** case 1a: increase width by |n| **/
		w_new = w_old - n + 1;  /** note: n is <= 0 **/
	   } else if(n > 1){   /** case 1b: **/
		for(start=2; null[start]; start++) ;
		if(n < w_old){
		    start = MINIMUM(Int4,start,n);
		    w_new = w_old - start + 1;
		} else if(n > w_old){
		    w_new = n - start + 1;   
		} else print_error("input error in mv_col_weight( )");
	   } else print_error("input error in mv_col_weight( )");
	} else if(lemon == w_old){	/** case 2: **/
	   for(end = w_old-1; null[end]; end--) ;
	   if(n < 1){
		w_new = end - n + 1;  /** note: n is <= 0 **/
	   } else if(n > 1){
		if(n < end) w_new = end;
		else if(n > end) w_new = n; 
		else print_error("input error in mv_col_weight( )");
	   } else print_error("input error in mv_col_weight( )");
	} else { 		/** case 3: 1 < lemon < w_old **/
	   if(n < 1){
		w_new = w_old - n + 1;  /** note: n is <= 0 **/
	   } else if(n > 1){
		if(n < w_old) w_new = w_old; 
		else if(n > w_old) w_new = n;
		else print_error("input error in mv_col_weight( )");
	   } else print_error("input error in mv_col_weight( )");
	}
#if 0	/** DEBUG **/
	b0 = bico(w_old-2,c-2);
	b2 = bico(w_new-2,c-2);
	P = b0/b2;
if(P > 1000000){
	fprintf(stderr,"w_old: ");
	s=MINIMUM(Int4,n,start);
	e=MAXIMUM(Int4,n,w_old);
	for(i=s; i<=w_old; i++){
	   if(i >= 1 && i <= w_old){
		if(null[i]) fprintf(stderr,".");
		else fprintf(stderr,"*");
	   }else fprintf(stderr," ");
	}
	fprintf(stderr,"\nw_new: ");
	s=MINIMUM(Int4,n,start);
	e=MAXIMUM(Int4,n,end);
	for(i=s; i<=e; i++){
	   if(i==lemon) fprintf(stderr,".");
	   else if(i==n) fprintf(stderr,"*");
	   else if(i >= s && i <= e){
		if(i >= start && i <= end && !null[i]) fprintf(stderr,"*");
		else fprintf(stderr,".");
	   }else fprintf(stderr," ");
	}
	s=MINIMUM(Int4,n,start);
	e=MAXIMUM(Int4,n,end);
	fprintf(stderr," (%d-%d)",s,e);
	fprintf(stderr,"\nw_old = %d w_new = %d; c= %d\n",
			w_old,w_new,c);
	fprintf(stderr,"P = %g/%g = %g; n=%d; lemon = %d\n",b0,b2,P,n,lemon);
}
#endif
	return (bico(w_old-2,c-2)/bico(w_new-2,c-2));
}

double	add_col_weight(Int4 c, Int4 w_old, Int4 New)
/********************************************************************

   If a fmodel has c on columns and a length of w then there are 
   w-2 choose c-2  configurations for that model.

          n     1             w_old
          |     |              |
        ........*??*???????????*.............
                 |--- (w-2) --|
               .....
              -10234
      n = -inf..0 or +2..+inf.	 (i.e., n != 1)

      n = -f..end.	(treat as above with lemon always in middle.)

 ********************************************************************/
{
	Int4	w_new;

	/** case 3: 1 < lemon < w_old **/
	if(New < 1){
		w_new = w_old - New + 1;  /** note: New is <= 0 **/
	} else if(New > 1){
		if(New < w_old) w_new = w_old; 
		else if(New > w_old) w_new = New;
		else print_error("input error in add_col_weight( )");
	} else print_error("input error in add_col_weight( )");
	return (bico(w_old-2,c-2)/bico(w_new-2,c-1));	/** c-1->added column **/
}


