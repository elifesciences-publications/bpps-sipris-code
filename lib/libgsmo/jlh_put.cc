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

#include "jlh_typ.h"

void    jlh_typ::put_dp_matrixSW(Int4 **D, Int4** T, Int4 len, smx_typ M,
        unsigned char *profile, unsigned char *seq, a_type A,
        Int4 startM,Int4 RelstartM,Int4 start, Int4 end)
{
        Int4    r,m,i,lenM,v;

	assert(startM > 0 && RelstartM > 0);
        end = MINIMUM(Int4,len,end);
        start = MAXIMUM(Int4,start,0);
        end = MAXIMUM(Int4,start,end);
        start = MINIMUM(Int4,start,end);
        lenM=LenSMatrix(M);
        fprintf(stderr,"     |");
        for(i=startM-1,m=RelstartM-1; m<=lenM; m++,i++) fprintf(stderr,"  %c  |",
                                AlphaChar(profile[i],A));
        fprintf(stderr,"\n     |");
        for(m=RelstartM-1; m<=lenM; m++) fprintf(stderr,"%3d  |",m);
        fprintf(stderr,"\n");
        for(r=start; r<=end; r++) {
           fprintf(stderr,"%c%4d|",AlphaChar(seq[r],A),r);
           for(i=startM-1,m=RelstartM-1; m<=lenM; m++,i++) {
		v=D[i][r];
                if(v <= -infinity) fprintf(stderr," -inf");
                else {
		   if(D != T) v = v/100;
		   fprintf(stderr,"%5d",v);
		}
                if(T[i][r] < 0) fprintf(stderr,"^");
                else if(T[i][r] == 0) fprintf(stderr,"\\");
                else fprintf(stderr,"<");
           }
           fprintf(stderr,"\n");
        }
        fprintf(stderr,"\n");
}

void    jlh_typ::put_dp_matrixSMX(Int4 len, smx_typ M, unsigned char *profile,
        unsigned char *seq, a_type A, Int4 startM,Int4 RelstartM,Int4 start, Int4 end)
{
        Int4    r,m,i,lenM;

	assert(startM > 0 && RelstartM > 0);
        end = MINIMUM(Int4,len,end);
        start = MAXIMUM(Int4,start,0);
        end = MAXIMUM(Int4,start,end);
        start = MINIMUM(Int4,start,end);
        lenM=LenSMatrix(M);
        fprintf(stderr,"     |");
        for(i=startM,m=RelstartM; m<=lenM; m++,i++) fprintf(stderr,"  %c  |",
                                AlphaChar(profile[i],A));
        fprintf(stderr,"\n     |");
        for(m=RelstartM; m<=lenM; m++) fprintf(stderr,"%3d  |",m);
        fprintf(stderr,"\n");
        for(r=start; r<=end; r++) {
           fprintf(stderr,"%c%4d|",AlphaChar(seq[r],A),r);
           for(m=RelstartM; m<=lenM; m++) {
                fprintf(stderr,"%5d|",ValSMatrix(m,seq[r],M)/100);
           }
           fprintf(stderr,"\n");
        }
        fprintf(stderr,"\n");
}

