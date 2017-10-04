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

#include "cmsa.h"
#include "blosum62.h"

static unsigned char RndRes()
{
        unsigned char a;
        double rand_num;
        rand_num = SampleUniformProb();
        for(Int4 i=1;i<=20;i++){
                if((rand_num-=blosum62freq[i]) <= 0.0) { a=i; break; } 
        }
        return a;
}

static Int4	**SimulateBlock(cma_typ cma, Int4 block)
{
	Int4		**freq;
	unsigned char	res;
	a_type 		A = AlphabetCMSA(cma);
	Int4 		n_seq = NumSeqsCMSA(cma);
	gss_typ 	*gss=gssCMSA(cma);
        Int4 		i,j,l,s,pos[5];
	Int4		len = LengthCMSA(block,cma);
	e_type 		fakeE;

	NEWP(freq,len+1,Int4);
	for(j=0;j<=len;j++)
		NEW(freq[j],nAlpha(A)+1,Int4);
	for(i=1;i<=n_seq;i++){
              	PosSiteCMSA(block,i,pos,cma);
               	fakeE=gss->FakeSeq(i);
		for(j=1;j<=len;j++){
	              	s=pos[1]+j-1;
			res = ResSeq(s,fakeE);
			if(res==0) res=RndRes();
			freq[j][res]+=1;
		}
        }
	return freq;
}

static Int4 **SimulateLeftGap(cma_typ cma, Int4 block, Int4 length)
{
	Int4            **freq;
        unsigned char   res,*ptr;
        a_type          A = AlphabetCMSA(cma);
        Int4            n_seq = NumSeqsCMSA(cma);
        Int4            nBl = nBlksCMSA(cma);
        gss_typ         *gss=gssCMSA(cma);
        Int4            i,j,len_bef,s,pos[5],pos_bef[5],gap;
        e_type          fakeE;
        double          rand_num;

        NEWP(freq,length+1,Int4);
        for(j=0;j<=length;j++)
                NEW(freq[j],nAlpha(A)+1,Int4);
        if(block>1){
	        len_bef = LengthCMSA(block-1,cma);
                for(i=1;i<=n_seq;i++){
                        fakeE=gss->FakeSeq(i);
                        ptr=SeqPtr(fakeE);
                        PosSiteCMSA(block,i,pos,cma);
                        PosSiteCMSA(block-1,i,pos_bef,cma);
			s=pos[1];
                        gap=s-pos_bef[1]-len_bef;
                        for(j=1;j<=length;j++){
                                if(j<=gap) res=ptr[s-j];
                                else res=RndRes();
                                freq[length-j+1][res]+=1;
                        }
                }
        }
        else{
                for(i=1;i<=n_seq;i++){
                        fakeE=gss->FakeSeq(i);
                        ptr=SeqPtr(fakeE);
                        PosSiteCMSA(block,i,pos,cma);
                        s=pos[1];
                        gap=s-1;
                        for(j=1;j<=length;j++){
                                if(j<=gap) res=ptr[s-j];
                                else res=RndRes();
                                freq[length-j+1][res]+=1;
                        }
                }
        }
	return freq;
}

static Int4 **SimulateRightGap(cma_typ cma, Int4 block, Int4 length)
{
	Int4 		**freq;
        unsigned char   res,*ptr;
        a_type          A = AlphabetCMSA(cma);
        Int4            n_seq = NumSeqsCMSA(cma);
	Int4 		nBl = nBlksCMSA(cma);	
        gss_typ         *gss=gssCMSA(cma);
        Int4            i,j,s,pos[5],pos_next[5],gap;
	Int4		len = LengthCMSA(block,cma);
        e_type          fakeE;

        NEWP(freq,length+1,Int4);
        for(j=0;j<=length;j++)
                NEW(freq[j],nAlpha(A)+1,Int4);
	if(block<nBl){
		for(i=1;i<=n_seq;i++){
			fakeE=gss->FakeSeq(i);
			ptr=SeqPtr(fakeE);
			PosSiteCMSA(block,i,pos,cma);
			PosSiteCMSA(block+1,i,pos_next,cma);
			s=pos[1]+len-1;
			gap=pos_next[1]-(s+1);
			for(j=1;j<=length;j++){
				if(j<=gap) res=ptr[s+j];
				else res=RndRes();
				freq[j][res]+=1;					
			}
		}
	}
	else{
		for(i=1;i<=n_seq;i++){
			fakeE=gss->FakeSeq(i);
			ptr=SeqPtr(fakeE);
			PosSiteCMSA(block,i,pos,cma);
			s=pos[1]+len-1;
			gap=LenSeq(fakeE)-s;
			for(j=1;j<=length;j++){
				if(j<=gap) res=ptr[s+j];
				else res=RndRes();
				freq[j][res]+=1;
			}
		}
	}
	return freq;
}

static unsigned char *ConvFreqToSeq(Int4 **freq, Int4 nalpha, Int4 length, Int4 nSeq)
{
	unsigned char 	*a;
	double 		rand_num;
	Int4 		*tfreq,j,i;

	NEW(a,length+1,unsigned char);

        for(j=1;j<=length;j++){
                tfreq=freq[j];
                rand_num=nSeq*SampleUniformProb();
                for(i=1;i<=nalpha;i++){
                        if((rand_num -= tfreq[i]) <= 0.0) { a[j]=i; break; }
                }
        }
	return a;
}


e_type *SimulatedSeqsCMSA(cma_typ cma, Int4 nSimSeq, Int4 rpts, Int4 *Gap_Len)
{
	e_type	*E;
	assert(nBlksCMSA(cma)==1);
	Int4 		i,j,k=1,l,gap,rem_gap,tot_len,maxgap=0;
	unsigned char 	*ptr,*sim_bl,*sim_left,*sim_right;
	Int4 		bl_len = LengthCMSA(1,cma);
	char		str[20];
	Int4 		n_seq=NumSeqsCMSA(cma);
	NEW(E,nSimSeq+1,e_type);
	a_type A=AlphabetCMSA(cma);
	for(tot_len=0,i=1;i<rpts;i++){
		tot_len+=Gap_Len[i];
		if(Gap_Len[i]>maxgap) maxgap=Gap_Len[i];
	}
	tot_len+=rpts*bl_len;
	maxgap /= 2; maxgap++;
	NEW(ptr,tot_len+1,unsigned char);

	Int4 **bl_freq=SimulateBlock(cma,1);
	Int4 **left_freq=SimulateLeftGap(cma,1,maxgap);
	Int4 **right_freq=SimulateRightGap(cma,1,maxgap);
	for(i=1;i<=nSimSeq;i++){
		k=1;
		sprintf(str,"simulated_%d",i);
		for(j=1;j<rpts;j++){
			l=1;
			sim_bl = ConvFreqToSeq(bl_freq,nAlpha(A),bl_len,n_seq);
			while(l<=bl_len){
				ptr[k++]=sim_bl[l++];
			}
			free(sim_bl);
			l=1;gap=Gap_Len[j];
			sim_right = ConvFreqToSeq(right_freq,nAlpha(A),gap/2,n_seq);
			while(l<=(gap/2)){
				ptr[k++]=sim_right[l++];
			}
			free(sim_right);
			l=1;rem_gap=gap-(gap/2);
			sim_left = ConvFreqToSeq(left_freq,nAlpha(A),rem_gap,n_seq);
			while(l<=rem_gap){
				ptr[k++]=sim_left[l++];
			}
			free(sim_left);
		}
		l=1;
		sim_bl = ConvFreqToSeq(bl_freq,nAlpha(A),bl_len,n_seq);
		while(l<=bl_len){
			ptr[k++]=sim_bl[l++];
		}
		free(sim_bl);
		E[i] = MkSeq(str,tot_len,ptr);
		// PutSeq(stdout,E[i],A);
	}
	for(i=0;i<=bl_len;i++) free(bl_freq[i]); free(bl_freq);
	for(i=1;i<=maxgap;i++) free(left_freq[i]); free(left_freq);
	for(i=1;i<=maxgap;i++) free(right_freq[i]); free(right_freq);
	return E;
}

#if 0
#define USAGE_START "USAGE: simulate cmafile <nseq> <gap1>,<gap2>,...<gapn> [options]\n\
options:\n\
   -n<int>			- number of simulated sequences\n\
   -r<int>:<gap1>,<gap2>..	-number of rpts and gaps between them\n\
   -s				- seed for the random number generator\n\
   -x         		  	- dummy\n\n"

main(Int4 argc,char *argv[])
{
	Int4		i,t,arg,*Gap_Len,nSimSeq=10,rpts;
	char		*Arg;
        UInt4 	seed=18364592;
	e_type 		*E;
	a_type          A=MkAlpha(AMINO_ACIDS,PROT_BLOSUM62);
        if(argc < 4) print_error(USAGE_START);
        for(arg = 4; arg < argc; arg++){
           if(argv[arg][0] != '-') print_error(USAGE_START);
           switch(argv[arg][1]) {
	     case 'n': if(sscanf(argv[arg],"-n%d",&nSimSeq)!=1)
		       print_error(USAGE_START); break;
             case 's': if(sscanf(argv[arg],"-s%d",&seed)!=1)
                       print_error(USAGE_START); break;
             default:  print_error(USAGE_START);
           }
	}
        if(seed == 18364592)  seed = (UInt4) time(NULL)/2;
        sRandom(seed);
	cma_typ cma = ReadCMSA2(argv[1],A);
	rpts = atoi(argv[2]);
	Arg = argv[3];
	i=0;t=1;
	while(Arg[i]!=NULL){
		if(Arg[i++] == ',') t++;
	}
	if(t+1!=rpts) {
		printf("There should be %d gaps!\n",rpts-1); 
		exit(1);
	} 

	NEW(Gap_Len,rpts,Int4);
	ParseIntegers(Arg,Gap_Len,"Check command line arguments!");
	E = SimulatedSeqsCMSA(cma,nSimSeq,rpts,Gap_Len);
}
#endif

