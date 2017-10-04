#include "psc_typ.h"

void	p2c_typ::DebugMapSqAln2Strct(FILE *fp, Int4 S,Int4 C, Int4 R, Int4 I, Int4 *TmpColToSeq,
					e_type pdbS, e_type pdbIC, e_type pdbE)
// print_error("p2c_typ::MapSqAlnToStruct( ) error: Redundant internal repeats?");
{
#if 0	// for debugging...can all be turned off.
	PutSeq(stderr,pdbIC,AB);
	AlnSeqSW(stderr,11,1,pdbIC,FullSeq[S],AB);
	PutSeq(stderr,pdbS,AB);
	AlnSeqSW(stderr,11,1,pdbS,FullSeq[S],AB);
	PutSeq(stderr,pdbE,AB);
	AlnSeqSW(stderr,11,1,pdbE,FullSeq[S],AB);
	PutSeq(stderr,FullSeq[S],AB);
#elif 1
	Int4	r,site,RR,col,os,os_cma,os_pdb,NumX;
	char	c2;
	fprintf(stderr,"\nFullSeq[%d]:\n",S);
	for(col=1; col <= NumCol(); col++){
	   fprintf(stderr,"col %d: ",col);
	   for(RR=1; RR <= NumFullRpts[S]; RR++){
		site=Col2FullSeq[S][RR][col];
		r=ResSeq(site,FullSeq[S]); c2=AlphaChar(r,AB); 
		fprintf(stderr," %c%d",c2,site);
	   }
	   site=TmpColToSeq[col];
	   r=ResSeq(site,pdbIC); c2=AlphaChar(r,AB);
	   fprintf(stderr," --> %c%d\n",c2,site);
	}
#if 0
	fprintf(stderr,"\npdbSeq[%d][%d]:\n",I,C);
	for(col=1; col <= NumCol(); col++){
	   fprintf(stderr,"col %d: ",col);
#if 0
	   for(RR=1; RR <= NumRpts[I][C]; RR++){
		Int4 site=Col2pdbSeq[I][C][RR][col];
		Int4 r=ResSeq(site,pdbIC); c2=AlphaChar(r,AB);
		fprintf(stderr," %c%d",c2,site);
	   } fprintf(stderr,"\n");
#else
	   Int4 site=TmpColToSeq[col];
	   Int4 r=ResSeq(site,pdbIC); c2=AlphaChar(r,AB);
	   fprintf(stderr," %c%d\n",c2,site);
#endif
	}
#endif
	os_pdb=OffSetSeq(pdbS); os_cma=OffSetSeq(pdbIC);
	
	char rtn=IsSameSeqFastX(pdbS,pdbIC,&os,&NumX,esc->MinSeqOverlap); // ignoring 'X' residues...
	Int4 d_os;
	switch (rtn){
		case 1: d_os=-os_pdb; break; 	// seqI N-terminus starts within seqJ.
		case 2: d_os=os_cma-os_pdb; break; 	// seqJ N-terminus starts within seqI.
	        default: assert(rtn < 3 && rtn > 0);
	} os = os-d_os;	// adjust for inherent offset in sequence.
	if(rtn==2) PutDiagonalSeq(stderr, os+d_os, pdbS,pdbIC,AB);
	else PutDiagonalSeq(stderr, os+d_os, pdbIC,pdbS,AB);
	PutSeq(stderr,pdbIC,AB);
	AlnSeqSW(stderr,11,1,pdbIC,FullSeq[S],AB);
	PutSeq(stderr,pdbS,AB);
	AlnSeqSW(stderr,11,1,FullSeq[S],FullSeq[S],AB);
	PutSeq(stderr,pdbE,AB);
#endif
#if 0
	for(col = 1; col <= NumCol; col++){
				c2=AlphaChar(ResSeq(TmpColToSeq[col],pdbE),AB);
				fprintf(stderr,"col %d: %c%d\n",col,c2,TmpColToSeq[col]);
	}
#endif
}

