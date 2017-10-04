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

Int4    ExtendGPsiBlstStr(Int4 lenMaster, Int4 i1, Int4 len2, unsigned char *seq2,
        Int4 i2, Int4 *left, Int4 *right, Int4 **matrix, register Int4 x_dropoff)
{
        register Int4   i,s,max,end,**mtx;
        register unsigned char  *p2;

        mtx=matrix+i1; p2=seq2+i2;
        if((max=len2-i2) > (s=lenMaster-i1)){ end=s; } else { end=max; }
        *right = 1;
        s = mtx[1][p2[1]];
        for(max=s, i=2; i <= end; i++) {
                if((s += mtx[i][p2[i]]) > max) { *right=i; max=s; }
                else if((s <= (max-x_dropoff)) || s < 0) break;
        } *right += i2; *left = 1;
        if(i1<=i2) end=-i1; else end=-i2;
        for(s=max, i=0; i > end; i--){
                if((s += mtx[i][p2[i]]) > max) { *left=i; max=s; }
                else if((s <= (max-x_dropoff)) || s < 0) break;
        } *left += i2;
        return max;
}


