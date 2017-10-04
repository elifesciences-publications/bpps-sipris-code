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

#include "segment.h"

s_type	Segment(Int4 I, Int4 offset)
{
	s_type  S;
	Int4	*s;

	if(offset > INT4_MAX || offset < INT4_MIN) seg_error("offset out of range");
	if(I < 1) seg_error("illegal segment identifier");
	s = (Int4 *) &S; s[0] = I; s[1] = offset;
	return S;
}

// Int4     SegmentI(s_type S){ short *s; s = (short*) &S; return (Int4) (s[0]>>1); }
// FIX BUG...
Int4     SegmentI(s_type S) { Int4 *s; s = (Int4 *) &S; return s[0]; }

Int4     SegmentStart(s_type S)
{ Int4 *s; s = (Int4 *) &S; return s[1]; }

Int4     BubbleSortSegments(s_type *S)
{
        UInt4 i,j,n;
	s_type	t;

        for(n=i=0; S[i]!=0; i++) n++;
        for(i=1; i < n; i++){
             for(j=i;j > 0 && (S[j-1] > S[j]);j--){ t=S[j]; S[j]=S[j-1]; S[j-1]=t; }
        } return n;
}

Int4     CopySegments(register s_type *SC, register s_type *S)
{
        register Int4     i;
        for(i=0; S[i]!=0;i++){ SC[i] = S[i];}
        SC[i] = 0;
        return i;
}
	
Int4     OffsetSegments(register Int4 os,register s_type *S1,register s_type *S2)
/***************************************************************
  set S2 to be S1 with offset os.

        offsetSegments(2,[1:2,1:7,2:9,2:12,NULL],S2) 
		-> S2 = [1:4,1:9,2:11,2:14,NULL]
  or
        offsetSegments(-2,[1:2,1:7,2:9,2:12,NULL],S2) 
		-> S2 = [0,5,7,10,-1]

***************************************************************/
{
        register UInt4     i;
	register Int4 *s1,*s2;
	register Int4	L1,L2;

        for(i=0; S1[i]!=0;i++){ 
#if 1
	        s1 = (Int4 *) &S1[i]; 
		s2 = (Int4 *) &S2[i];
		s2[1] = s1[1] + os; s2[0] = s1[0];
#elif 0	// old, problematic code:
	        s1 = (unsigned short*) &S1[i]; 
		s2 = (unsigned short*) &S2[i];
		s2[1] = s1[1] + os; s2[0] = s1[0];
		if((s2[1] < 0) != (s1[1] < 0)) s2[0]=s2[0] ^ 1;
#else	// modify to guarrantee correct bit configureation.
	        s1 = (unsigned short*) &S1[i]; 
		s2 = (unsigned short*) &S2[i];
		if((s1[0] & 1) == 0){ L1 = -1 * (Int4) s1[1]; }
		else { L1 = (Int4) s1[1]; }
		L2 = L1 + os;
		s2[0] = s1[0]; 
		if(L2 < 0) s2[0] = s2[0] & 0xFFFE;	// set last bit to 0.
		else s2[0] = s2[0] | 1;	// set last bit to 1.
		s2[1] = (unsigned short) abs(L2);
#endif
	} S2[i] = 0;
        return i;
}

#if 0
BooLean	MemberSegments(s_type *L, s_type S)
/* return TRUE if S is a member of L; else return FALSE */
/* assumes that segment list is arranged in increasing order */
{
        Int4 i;

        for(i=0; L[i]!=0; i++){ 
		if(S < L[i]) return FALSE;
		else if(S == L[i]) return TRUE; 
	}
        return FALSE;
}
#endif

Int4     IntersectSegments(register s_type *SC, register s_type *S1,
	register s_type *S2)
/* set SC to be the intersection of S1 and S2; returns cardinality of SC */
{
        register Int4 i,j,c;

        for(i=j=c=0; S1[i]!=0 && S2[j]!=0; ){
                if(S1[i] < S2[j]) i++; 
                else if(S1[i] > S2[j]) j++; 
                else { SC[c++] = S1[i]; i++; j++; }
        } SC[c]=0; 
        return c;
}

Int4     UnionSegments(s_type *SC,s_type *S1,s_type *S2)
/* set SC to be the union of S1 and S2; returns cardinality of SC */
{
        Int4 i,j,c;

        for(i=j=c=0; S1[i]!=0 && S2[j]!=0; ){
           if(S1[i] < S2[j]) SC[c++] = S1[i++];
           else if(S1[i] > S2[j]) SC[c++] = S2[j++];
           else { SC[c++] = S1[i]; i++; j++; }
        }
        while(S1[i]!=0) SC[c++]=S1[i++];
        while(S2[j]!=0) SC[c++]=S2[j++];
        SC[c]=0; 
        return c;
}

#if 0
Int4     DiffSegments(s_type *SD,s_type *S1,s_type *S2)
/* set SD to be S2 - S1, that is, SD contains items in list S2 that 
   are not in list S1; returns cardinality of SD */
{
        Int4 i,j,c;

	BubbleSortSegments(S1);
	BubbleSortSegments(S2);
        for(i=j=c=0; S1[i]!=0 && S2[j]!=0; ){
           if(S1[i] < S2[j]) i++;
           else if(S1[i] > S2[j]) SD[c++] = S2[j++];
           else { i++; j++; }
        }
        SD[c]=0; 
        return c;
}
#endif

void	PutSegment(s_type S)
{
	Int4 *s = (Int4 *) &S;
	fprintf(stderr,"S.I = %d, S.o = %d\n",s[0],s[1]);
}

void	seg_error(char *s) { fprintf(stderr,"Segment: %s\n",s); exit(1); }

Int4	IntersectOSegments(Int4 *C, Int4 r, Int4 *os, Int4 k1, Int4 k2, register s_type *S1,
	register s_type *S2)
/************************************************************************
  set S to be S1 with offset os.

        offsetSegments(2,[1:2,1:7,2:9,2:12,NULL],S2) 
		-> S = [1:4,1:9,2:11,2:14,NULL]
  or
        offsetSegments(-2,[1:2,1:7,2:9,2:12,NULL],S2) 
		-> S = [1:0,1:5,2:7,2:10,NULL]
 
  then find intersection of S and S2.
  check all offsets and set *C and *os to the best of these 
  return the number of offsets with c > 1.

	for r = 3

        S1 (k1=11)
       xxxxxxxxAxx				xxxxxxxxAxx
       xxxxxxxxAxx                              xxxxxxxxAxx
       xxxxxxxxAxx	o = 3-k2..k1-3		xxxxxxxxAxx
  Axxxxxxx		  = -5..8		        Axxxxxxx
  Axxxxxxx                                              Axxxxxxx
  Axxxxxxx                                              Axxxxxxx
  Axxxxxxx                                              Axxxxxxx
  S2 (k2=8)

	NOTE: s = (unsigned short*) &S; cannot be placed outside of loop!!!
	returns the number of offsets with an intersect > 0
 ************************************************************************/
{
        register Int4	i,j,c,o;
	register Int4 	*s,*s1;
	s_type		S;
	Int4		n;
	Int4		L;

        for(*C=n=0,o=r-k2,k1-=r; o <= k1; o++){
            for(i=j=c=0; S1[i]!=0 && S2[j]!=0; ){
#if 1
		S = S1[i];  s = (Int4 *) &S; s[1]+=o;
#elif 0
		S = S1[i];  s = (unsigned short*) &S; 
		if((s[1] < 0) != ((s[1]+o)<0)) s[0]=s[0] ^ 1;
		s[1]+=o;
#else 	// fix problems due to (presumably) twos-complement issues...
		S = S1[i];  s = (unsigned short*) &S; 

		if((s[0] & 1) == 0){ L = -1 * (Int4) s[1]; }
		else { L = (Int4) s[1]; }
		L += o;

		if(L < 0){
			s[0] = s[0] & 0xFFFE; // set last bit to 0.
			s[1] = (unsigned short) abs(L);
		} else {
			s[0] = s[0] | 1;      // set last bit to 1.
			s[1] = (unsigned short) L;
		}
#endif
                if(S < S2[j]) i++; 
                else if(S > S2[j]) j++; 
                else { c++; i++; j++; }
            }
	    if(c > 0){  /*** or use 1 ***/ 
		if(c > *C) { *C=c; *os=o; }
		else if(c == *C){ if(abs(o) < abs(*os)) *os=o; }
		n++; 
	    }
	} return n;
}


