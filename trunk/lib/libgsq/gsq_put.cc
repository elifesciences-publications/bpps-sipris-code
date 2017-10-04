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

#include "gsq_typ.h"

void	gsq_typ::PutFake(FILE *fptr, a_type A) { PutSeq(fptr, fakeE, A); }

void    gsq_typ::Put(FILE *fp, a_type A)
{
	unsigned char *rseq,*fseq;
	Int4	i,j,k,length,TrueLen;
	char	c1,c2,operation;

	assert(fakeE != NULL);
	TrueLen=LenSeq(realE);
	length=LenSeq(fakeE); rseq=SeqPtr(realE); fseq=SeqPtr(fakeE);
	fprintf(fp,"     r f   ins del f2r oper\n");
        for(i=0; i<=length; i++){
	    if(!deletion(i))
	    { // no deletion here
                c2 = AlphaChar(rseq[f2r[i]],A);
                c1 = AlphaChar(fseq[i],A);
		if(f2r[i] > 0 && f2r[i] <= TrueLen) operation='m';
		else operation='d';
	        fprintf(fp,"%3d: %c %c  %3d %3d %3d %c\n",i,c2,c1,
			insert(i),deletion(i),f2r[i],operation);
	    } else {				// deletion ='-'.
                c1 = AlphaChar(fseq[i],A);
		if(f2r[i] > 0 && f2r[i] <= TrueLen) operation='d';
		else operation='d';
		fprintf(fp,"%3d: . %c  %3d %3d %3d %c\n",i,c1,
			insert(i), deletion(i), f2r[i],operation);
	    }
	    k = insert(i);
	    for(j=f2r[i]; k > 0; k--){
                j++; c2 = tolower(AlphaChar(rseq[j],A));
		operation='i';
		fprintf(fp,"%3d: %c .    :   - %3d %c\n",i,c2,j,operation);
	    } 
        }
	fprintf(fp,"num_insert=%d; num_del=%d; num_open=%d.\n",
		num_insert,num_del,num_open);
	fprintf(fp,"#### offset real = %d; offset fake = %d\n",
                OffSetSeq(realE),OffSetSeq(fakeE));
	PutSeq(fp,realE,A); PutSeq(fp,fakeE,A);
}

void    gsq_typ::PutX(FILE *fp, a_type A)
{
        unsigned char *rseq,*fseq;
        Int4    i,j,k,length,TrueLen,real;
        char    c1,c2,operation;

        assert(fakeE != NULL);
        TrueLen=LenSeq(realE);
        length=LenSeq(fakeE); rseq=SeqPtr(realE); fseq=SeqPtr(fakeE);
        fprintf(fp,"     r f   ins del f2r state r2f \n");
        for(real=0,i=0; i<=length; i++){
            if(!deletion(i))
	    { // no deletion here?
                c2 = AlphaChar(rseq[f2r[i]],A);
                c1 = AlphaChar(fseq[i],A);
                if(f2r[i] > 0 && f2r[i] <= TrueLen) operation='m';
                else operation='d';
                fprintf(fp,"%3d: %c %c  %3d %3d %3d  %c   %3d\n",i,c2,c1,
                        insert(i),deletion(i),f2r[i],operation,RealToFake(f2r[i]));
		real=f2r[i];
            } else {                            // deletion ='-'.
                c1 = AlphaChar(fseq[i],A);
                if(f2r[i] > 0 && f2r[i] <= TrueLen) operation='d';
                else operation='d';
                fprintf(fp,"%3d: . %c  %3d %3d %3d  %c     -\n",i,c1,
                         insert(i), deletion(i), f2r[i],operation);
            }
            k = insert(i);
            for(j=f2r[i]; k > 0; k--){
                j++; c2 = tolower(AlphaChar(rseq[j],A));
                operation='i'; real++;
                fprintf(fp,"%3d: %c .    :   -(%3d) %c\n",i,c2,real,operation);
            }
        }
        // fprintf(fp,"num_insert=%d; num_del=%d; num_open=%d.\n", num_insert,num_del,num_open);
        // fprintf(fp,"#### offset real = %d; offset fake = %d\n", OffSetSeq(realE),OffSetSeq(fakeE));
        // PutSeq(fp,realE,A);
}

void    gsq_typ::Put(FILE *fp, Int4 line_len, a_type A)
/********************************************************************
'-' = deletion; lowercase = insertion; uppercase = retained residue.
>gsqid1 {|N-term(C-term)[tax_group]|}description
plwlbfPLWLGVRPFFRALFRIEVIG--KEnipEGACIVASNHRSHLDPPVL
PLVflaKEELFK----GILKsHMRAIPLRRGRSEDISTLEECvslLKLGCKI 
GIFPEGTRANPGEFKRAKPGAVGFLAINSGFQPVLPVYIDGTDRAFPRGKga
 ********************************************************************/
{
	Int4	strlen,i;
	char	*buffer=new char[LenSeq(fakeE)+num_insert+5];

	assert(fakeE != NULL); 
	fprintf(fp,">"); PutSeqInfo(fp,fakeE); 
	// strlen=Region(buffer,0,LenSeq(fakeE)+1,A);
	strlen=Region(buffer,1,LenSeq(fakeE),A);
	for(i=0; i < strlen; ){
		fprintf(fp,"%c",buffer[i]); i++;
		if(i%line_len == 0) fprintf(fp,"\n");
	} fprintf(fp,"\n\n");
	delete []buffer;
}

void	gsq_typ::Put_cma_format(FILE *fp,Int4 J0,Int4 nblk,Int4 *start, Int4 *blk_len,a_type A)
// blk_len[nblk+1]==0;
{
	Int4	len_real,i,j,k,b,site,s,true_len=LenSeq(realE),fake_len=LenSeq(fakeE);
	unsigned char *rseq=SeqPtr(realE),*fseq=SeqPtr(fakeE);
	char    c;
	// assert(blk_len[nblk+1]==0);

	fprintf(fp,"$%d=%d(%d):\n",J0,true_len,fake_len);
	// if(Labeling == 'L') LabelSeq(realE); else if(Labeling == 'U') UnLabelSeq(realE);
	fprintf(fp,">"); PutSeqInfo(fp,realE);
	// PutSeqRFF(fp,realE);     // output rff if and only if defined.
	for(len_real=0,s=1; s <= this->OverHangN(); s++){		// print out N-terminal extension
		fprintf(fp,"%c",AlphaChar(ResSeq(s,realE),A));  len_real++;
	} fprintf(fp,"{(");
	for(site=LenSeq(fakeE),b=0,i=1; i<=fake_len; i++,site++){
	   if(b < nblk && i == start[b+1]){ b++; site=1; fprintf(fp,")"); }
	   if(!deletion(i))		// no deletion here
		{ c = AlphaChar(rseq[f2r[i]],A); len_real++; } else c='-';
	   fprintf(fp,"%c",c);
	   if(site == blk_len[b]){ fprintf(fp,"("); assert(b <= nblk);  }
	   if(i == fake_len) fprintf(fp,")}");
	   k =insert(i);
           for(j=f2r[i]; k > 0; k--){			// insertion here.
		j++; c=AlphaChar(rseq[j],A); 
		if(i < fake_len) c = tolower(c); 
           	fprintf(fp,"%c",c); len_real++;
           } 
        }
#if 1	// Bug fix ...
	while(len_real < true_len){		// print out C-terminal extension
		len_real++; fprintf(fp,"%c",AlphaChar(ResSeq(len_real,realE),A));  
	} 
#endif
	fprintf(fp,"*\n\n");
// this->Put(stderr,A);
	assert(len_real == true_len);
}

