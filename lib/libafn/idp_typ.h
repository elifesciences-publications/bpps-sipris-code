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

#if !defined(_IDP_TYP_)
#define _IDP_TYP_

#include "stdinc.h"
#include "afnio.h"

/********************** Insertion & Deletion Penalty Class *************************

   length  = length of motif block.
   prob[i] = probability of 

 *******************************************************************************/
class idp_typ {
public: 
		idp_typ( ){ assert(!"Illegal constructor"); }
		idp_typ(idp_typ *, double *, double );
		idp_typ(char*,char,double*,unsigned short,unsigned short,Int4,float*);
		idp_typ(char* IdpArg, char Mode,double *Prob,unsigned short Length,
			unsigned short Rpts, Int4 NumRes,float *SeqProb,
			unsigned short nBlk, unsigned short *BlkLen,Int4 PerNats);

		idp_typ(char*,char,double*,unsigned short,unsigned short,Int4,float*,
		  unsigned short,unsigned short*);
		idp_typ& operator=(const idp_typ&);	// assignment operator.
		~idp_typ( ){ Free(); }
	void	Put(FILE *);
	BooLean	IsSeqPenalty(){ return (ins_seq_high!=0 || ins_seq_low!=0);}
	// OLD function names
	Int4	*InsOpen( ){ return m2i; }
	Int4	*InsExtend( ){ return i2i; }
	Int4	*DelOpen( ){ return m2d; }
	Int4	*DelExtend( ){ return d2d; }
	Int4	InsOpen(Int4 i){ return m2i[i]; }
	Int4	InsExtend(Int4 i){ return i2i[i]; }
	Int4	DelOpen(Int4 i){ return m2d[i]; }
	Int4	DelExtend(Int4 i){ return d2d[i]; }
	// END of OLD function names
	Int4	Length( ){ return length; }
	Int4	Length(Int4 b){ assert(b > 0 && b <= nblk); return blklen[b]; }
	Int4	nBlks( ){ return nblk; }
	Int4    StartBlk(Int4);
	Int4    EndBlk(Int4);
	Int4	MaxRpts( ){ return rpts; }
	Int4	**GapFunct( ){ return gpen; }
	Int4    *SeqPenalties(unsigned char *seq, Int4 len)
		{ if(seq_ins_prob) return ComputeGREASE(seq,len,3); else return 0; }
	Int4	*MatToIns( ){ return m2i; }
	Int4	*MatToDel( ){ return m2d; }
	Int4	*InsToIns( ){ return i2i; }
	Int4	*DelToDel( ){ return d2d; }
	Int4	*MatToMat( ){ return m2m; }
	Int4	*InsToMat( ){ return i2m; }
	Int4	*DelToMat( ){ return d2m; }
	Int4	MatToIns(Int4 i){ return m2i[i]; }
	Int4	MatToDel(Int4 i){ return m2d[i]; }
	Int4	InsToIns(Int4 i){ return i2i[i]; }
	Int4	DelToDel(Int4 i){ return d2d[i]; }
	Int4	MatToMat(Int4 i){ return m2m[i]; }
	Int4	InsToMat(Int4 i){ return i2m[i]; }
	Int4	DelToMat(Int4 i){ return d2m[i]; }
private:
	Int4	*ComputeGREASE(unsigned char*,Int4,Int4);
	void	Init(char*,char,double *,unsigned short,unsigned short, unsigned short,
        	unsigned short *,Int4 ,float *);
	void	Init(char*,char,double *,unsigned short,unsigned short, unsigned short,
        	unsigned short *,Int4 ,float *,Int4);
	void	init( );
	void	Free();		// free memory...
	Int4	inso_low,inso_high;  // range for insertion opening penalties.
	Int4	insx_low,insx_high;  // range for profile insertion extension penalties.
	Int4	ins_seq_low,ins_seq_high; // range for seq. insertion extension penalties.
	Int4	delo;			// deletion opening penalty.
	Int4	delo_low,delo_high;	// deletion opening penalty.
	Int4	delx_low,delx_high; 	// range for deletion extension penalties.
	unsigned short	rpts, length, nblk,*blklen;
	Int4	*i2i,*i2m,*i2d;
	Int4	*m2i,*m2m,*m2d;
	Int4	*d2d,*d2m,*d2i;
	Int4	pernats;
	double	*prob;
	Int4	num_res;	// number of residue types (21)
	float	*seq_ins_prob; // corresponds to residues...
	char	mode;
	char	*idp_arg;
	Int4	**gpen;		// gap function between blocks...
	Int4	interblkL,interblkH;
};

#endif

