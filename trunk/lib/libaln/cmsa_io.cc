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

cma_typ AddRandomCMSA(cma_typ cma, Int4 num_random)
{ cma_typ rcma; MkMainFileCMSA(cma, num_random,rcma); TotalNilCMSA(rcma); }

cma_typ MkMainFileCMSA(cma_typ cma, Int4 num_random,cma_typ &rcma)
{
	a_type AB=AlphabetCMSA(cma);
        cma_typ TmpCMA[3]; TmpCMA[1]=cma;
        FILE *tfp = tmpfile();
        PutRandomCMA(tfp,blosum62freq,LengthCMSA(1,cma),num_random,AB);
        rewind(tfp); TmpCMA[2]=ReadCMSA(tfp,AB); fclose(tfp);
        tfp = tmpfile(); PutMergedCMSA(tfp,2,TmpCMA); rewind(tfp);
        cma_typ mcma=ReadCMSA(tfp,AB); fclose(tfp);
        // TotalNilCMSA(TmpCMA[2]); // destroy the temporary Random sequence alignment.
        rcma=TmpCMA[2];
        return mcma;
}

cma_typ	MakeConsensusCMSA(cma_typ cma)
{
	FILE *fp=tmpfile(); PutConsensusCMSA(fp,cma); rewind(fp);
	cma_typ CSQ=ReadCMSA(fp,AlphabetCMSA(cma)); fclose(fp); return CSQ;
}

void    PutConsensusCMSA(FILE *fp, cma_typ cma)
// print out a concensus sequence cma file
{
        e_type  cE=MkConsensusCMSA(cma);
        a_type A=AlphabetCMSA(cma);
        fprintf(fp, "[0_(1)=%s(1){go=10000,gx=2000,pn=1000.0,lf=0,rf=0}:\n",
                        NameCMSA(cma));
        fprintf(fp,"(%d)",LenSeq(cE));
        for(Int4 i=1; i <= LenSeq(cE); i++) fprintf(fp,"*");
        fprintf(fp,"\n\n$1=%d(%d):\n",LenSeq(cE),LenSeq(cE));
	 
	if(strchr(NameCMSA(cma),'/') != 0 || strlen(NameCMSA(cma)) > 50) fprintf(fp,">consensus seq\n");
        else fprintf(fp,">%s_consensus {<Consensus(C)>}\n",NameCMSA(cma));
        char *seq; NEW(seq,LenSeq(cE) +3, char);
        SeqToString(seq, cE, A);
        fprintf(fp,"{()%s()}*\n\n_0].\n",seq);
        free(seq);
        NilSeq(cE);
}

void PutRandomCMA(FILE *fp, double *freq, Int4 Length, Int4 NumSeqs, a_type AB)
// output a cma file of 'NumSeqs" random sequences based on the amino acid freq array.
{
	Int4	i,j,s,r,R,n;
	double	*Freq,total;
	static Int4 NumRandEqualOne=0;
	for(r=1,total=0.0; r <= nAlpha(AB); r++){ total+=freq[r]; }
	NEW(Freq,nAlpha(AB) +3, double);
	for(r=1; r <= nAlpha(AB); r++){ Freq[r]=freq[r]/total; }
	fprintf(fp,"[0_(1)=Random(%d){go=19,gx=2,pn=5.0,lf=0,rf=0}:\n",NumSeqs);
	fprintf(fp,"(%d)",Length);
	for(n=1; n<=Length; n++) fprintf(fp,"*"); fprintf(fp,"\n\n");
	// h_type HG=Histogram("amino acid residues",0,25,1.0);
	for(s=1; s <= NumSeqs; s++){
	   fprintf(fp,"$%d=%d(%d):\n>",s,Length,Length);
	   fprintf(fp,"Random%d {<simulated(S)>} \n{()",s);
	   for(i=1; i <= Length; i++){
		double rand = (double) Random()/(double) RANDOM_MAX;
		assert(rand >= 0.0 && rand <= 1.0);
		for(R=0,total=0.0,r=1; r <= nAlpha(AB); r++){
			total += Freq[r]; if(total >= rand){ R=r; break; } 
		} 
		if(rand == 1){ NumRandEqualOne++; assert(NumRandEqualOne < 3); }
		// 1 out of every 2,147,483,647 calls to Random will return 1!
		// total may be less than 1.0 due to rounding errors (e.g. got 0.99999999999999967).
		else assert(R != 0);
		fprintf(fp,"%c",AlphaChar(R,AB));
	        // IncdHist((double)R,HG);
	   } fprintf(fp,"()}*\n\n");
	} fprintf(fp,"_0].\n\n");
        // PutHist(stderr,50,HG); NilHist(HG); 
	free(Freq);
}

void	WriteCMSA(char *filename,cma_typ cma)
{ FILE *fp = open_file(filename,"","w"); PutCMSA(fp,cma); fclose(fp); }

#if 0
/***************************************************************
========================= ALIGNED SEQUENCES =================
 A way to represent all information for a gapped multiple 
 sequence alignment in order to read and write *.msa type
 file information more efficiently.

 alignment_template (can be used to create sites data structure)
 '*' = use always.  '^' = ignore always.
 '.' = ignore by default.  '!' = key residue.

 ')' = start of a motif region.  '(' = end of a motif region.
 )( = motif entirely absent in sequence (i.e., location unimportant).
 uppercase = motif residue (or region between motifs).
 lowercase = insertion (in motif) relative to model.
 '-' = deletion relative to model.
 (no deletions outside of motif; ignored).
 '{' & '}' denotes start and end of fakeSeq.

 Can convert to alignment operations:

 LTPVGKYLLL {(PL)WLGVrPFF(RALFR)IEVIGKEN--()--AcIVASN-RS-(LDPP)} STL
 .......... B ii MmmmImmm iiiii Mmmmmmmmdd  DdmImmmmmdmmd iiii E ...

 FOR EXAMPLE: --------------------------------------------
[0_(6)=barth(32){go=10,gx=5,pn=5.0,lf=9,rf=9}:
(13)*..*.*.*.****
(22)********.*.****.*..***
(32)***...*..**.*.***.**.*****.....*
(21)**.*..*.*******....**
(13)******.***.**
(10)********.*

 :   :   :

_0].                    // end of alignment file.

 If Fake==Real then 

$1=RealLength:
>id_TrueSeq1
(LGRKLLKLYFKL)YHRIDVRGLENLP...(AYEAWSVYD)*

  or:

$1=4{(5)"EMdmmmmmMmmmmmmmmmiMmImmmmmmmmmmmmiiiiiiiiMmmmmmmiMmmmmDmmmmmmmmiiii...E"(4)}7*

 for another sequence derived from $1 realsq (for storing files with multiple cma alignments
 of the same input sequences).

 output file == *.cma  (can convert this to *.msa, *.rtf, *.psi etc.
 ***************************************************************/
#endif
void	PutCMSA(FILE *fp,cma_typ cma) {	PutSelectCMSA(fp,0,cma); }

void	WriteGoodCMSA(char *filename,double cutoff, cma_typ cma)
{
	FILE *fp = open_file(filename,".cma","w");
	PutGoodCMSA(fp, cutoff ,cma); fclose(fp);
}

void	PutGoodCMSA(FILE *fp,double cutoff, cma_typ cma)
{
        Int4 	N = NumSeqsCMSA(cma);
        ss_type data = TrueDataCMSA(cma);
   	BooLean *skip; NEW(skip,N+3,BooLean);
	double prob; 
	for(Int4 J=1; J <= N; J++){ 
	     for(Int4 m=1; m <= nBlksCMSA(cma); m++){
	       // prob = GetProbCMSA(m,J,cma);
	       prob = GetGappedProbCMSA(m,J,cma);
	       if(prob < cutoff) {
			PutSeqSetE(stdout,J,data);
			skip[J]=TRUE; break;
	       }
	     }
	} PutSelectCMSA(fp,skip,cma); free(skip);
}

void	PutBestCMSA(FILE *fp,Int4 Num, BooLean KeepFirst, cma_typ cma)
// Print out the Num best sequences in the alignment.
{
        Int4 	i,j,N = NumSeqsCMSA(cma),start;
	dh_type dH=dheap(N+2,4);
        ss_type data = TrueDataCMSA(cma);
   	BooLean *skip; NEW(skip,N+3,BooLean);
	if(KeepFirst){ start=2; Num--; } else start=1;
	assert(Num > 0);
	double prob; 
	for(Int4 J=start; J <= N; J++){ 
	     prob=0.0;
	     for(Int4 m=1; m <= nBlksCMSA(cma); m++){
	       prob += GetProbCMSA(m,J,cma);
	       // prob = GetGappedProbCMSA(m,J,cma);
	     }
	     if(!isfinite(prob)) prob=-99999999999999.9;
	     insrtHeap(J,(keytyp) prob,dH);
	} for(i=0; ((j=delminHeap(dH)) != 0); ){
		assert(j <= N); skip[j]=TRUE;
		if(ItemsInHeap(dH) <= Num) break;
	} PutSelectCMSA(fp,skip,cma); free(skip);
	Nildheap(dH);
}

Int4 PutFullSeqCMSA(FILE *fp,BooLean *skip, cma_typ cma)
{
        Int4 fullN,FullN,s,J,i;
        unsigned short  *rpts=0;
        ss_type data = cma->FullSeq;
        fullN=FullN=NSeqsSeqSet(data);
	gss_typ& gss=*gssCMSA(cma);
	Int4 number=gss.NumSeq();

        BooLean  *Skip=0;
	Int4	totalRpts=0;
        if(skip){
             NEW(Skip,FullN+3,BooLean);
             NEW(rpts,FullN+3,unsigned short);
             for(J=1; J <= number; J++){
                s=cma->SubToFull[J];
                assert(s > 0 && s <= FullN);
                if(!skip[J]) rpts[s]++;
             }
             for(s=1; s <= FullN; s++){
                if(rpts[s] == 0){ fullN--; Skip[s]=TRUE; }
                else totalRpts+=rpts[s];
             }
        } else {
             rpts=cma->FullRpts;
             for(s=1; s <= FullN; s++){
                if(rpts[s] == 0){ fullN--; } else totalRpts+=rpts[s];
             }
        } assert(fullN > 0);
        fprintf(fp,"<FullSeqs(%d)=",fullN);
        for(i=0,s=1; s <= FullN; s++){
             if(rpts[s] != 0){
                   if(i) fprintf(fp,",%d(%d)",SqLenSeqSet(s,data),rpts[s]);
                   else fprintf(fp,"{%d(%d)",SqLenSeqSet(s,data),rpts[s]);
                   i++;
             }
        } fprintf(fp,"};\n");
        for(s=1; s <= FullN; s++){
             if(rpts[s] != 0){ PutSeqSetE(fp,s, data); }
        } fprintf(fp,">.\n");
        if(cma->Domains) cma->Domains->Put(fp,Skip);
        if(Skip) { free(rpts); free(Skip); }
	return totalRpts;
}

void	PutFastaAlnCMSA(FILE *fp,cma_typ cma)
// output a single block cmsa file in fasta format...
/// WARNING: assumes that seqset for ma is derived from gss_typ.
// eventually change so that uses gss_typ& gss = gssCMSA(cma);
{
	Int4	lenM,LenStr,AlignLen,J,x,*s;
	char	**alignment;
	st_type	S=SitesCMSA(cma);
	ss_type	data=DataCMSA(cma);
	a_type	A=SeqSetA(data);
	e_type	sgE,E;
	fm_type	*M=ModelsCMSA(cma);
	Int4	site,number,buffer_size;
	char	*ps;
	BooLean	**null;
	gss_typ     *gss=gssCMSA(cma);

	// assert(gss->Gapped());
	number=gss->NumSeq();
	assert(nBlksCMSA(cma) == 1);
	x = nColsFModel(M[1]);
	buffer_size = gss->MaxTrueSqLen();
	for(J=1; J <= number; J++) buffer_size+=gss->NumIns(J)+gss->NumDel(J);
	s = new Int4[number+3];
	null = NullCMSA(cma);

	//// 1. Get locations of sites in gapped sequences from cmsa.
	lenM = SiteLen(1,S);
	NEWP(alignment,number+3,char);
	for(J=1; J <= number; J++){ 
		sgE = gss->FakeSeq(J);
		s[J] =SitePos(1,J,1,S);
// fprintf(stderr,"s[%d] = %d\n",J,s[J]);
		NEW(alignment[J],buffer_size+3,char); 
	   	//// 2. Get corresponding gapped alignable segments.
		LenStr=gss->Region(J,alignment[J],s[J],lenM);
	}
	//// 3. Align gapped sequences.
	AlignLen=cmsa_fix_alignment(lenM,number,alignment,null[1]);

	//// 4. Print gapped alignment.
	for(J=1; J <= number; J++){ 
		sgE = gss->FakeSeq(J);
		site = gss->TrueSite(J, s[J]);
		ps=alignment[J];
		E = gss->TrueSeq(J);
        	// fprintf(fp,">"); PutSeqInfo2(fp,E);
        	fprintf(fp,">"); PutSeqInfo(fp,E);
		for(Int4 i=0; alignment[J][i]; i++){
		   char r=alignment[J][i];
		   if(r == '.') r='-';
		   fprintf(fp,"%c",r);
		}
		fprintf(fp,"\n\n");
		for(ps++; *ps != 0; ) { if(isalpha(*ps)) site++; ps++; }
	} fflush(fp);
	for(J=0; J <= number; J++){ free(alignment[J]); } free(alignment);
	// alignment[0] allocated above....
	free(null[1]); free(null);
	delete []s;
}

void	PutStockholmCMSA(FILE *fp,cma_typ cma,set_typ set)
// output a single block cmsa file in sto format...
/// WARNING: assumes that seqset for ma is derived from gss_typ.
// eventually change so that uses gss_typ& gss = gssCMSA(cma);
{
	Int4	lenM,LenStr,AlignLen,i,j,J,x,*s;
	char	**alignment;
	st_type	S=SitesCMSA(cma);
	ss_type	data=DataCMSA(cma);
	a_type	A=SeqSetA(data);
	e_type	E;
	fm_type	*M=ModelsCMSA(cma);
	Int4	site,number,buffer_size;
	char	*ps;
	BooLean	**null;
	gss_typ     *gss=gssCMSA(cma);

	// assert(gss->Gapped());
	number=gss->NumSeq();
	if(number >= 1000000) 
		print_error("FATAL: PutStockholmCMSA(): >= 1 million input sequences!");
	if(set && SetN(set) <= number) 
		print_error("FATAL: PutStockholmCMSA() set size too small");
	assert(nBlksCMSA(cma) == 1);
	x = nColsFModel(M[1]);
	buffer_size = gss->MaxTrueSqLen();
	for(J=1; J <= number; J++) buffer_size+=gss->NumIns(J)+gss->NumDel(J);
	s = new Int4[number+3];
	null = NullCMSA(cma);

	//// 1. Get locations of sites in gapped sequences from cmsa.
	lenM = SiteLen(1,S);
	NEWP(alignment,number+3,char);
	for(J=1; J <= number; J++){ 
		// if(set && !MemberSet(J,set)) continue;
		s[J] =SitePos(1,J,1,S);
// fprintf(stderr,"s[%d] = %d\n",J,s[J]);
		NEW(alignment[J],buffer_size+3,char); 
	   	//// 2. Get corresponding gapped alignable segments.
		LenStr=gss->Region(J,alignment[J],s[J],lenM);
	}
	//// 3. Align gapped sequences.
	AlignLen=cmsa_fix_alignment(lenM,number,alignment,null[1]);

	//// 4. Print gapped alignment.
	fprintf(fp,"# STOCKHOLM 1.0\n#=GF ID %s\n\n",NameCMSA(cma));
	char	str[60],str2[100];
	Int4	start,end;
	for(start=0,end=50; start < AlignLen; ){
	   for(J=1; J <= number; J++){ 
		if(alignment[J] == 0) continue;
		E = gss->TrueSeq(J); StrSeqID(str, 40,E);
		sprintf(str2,"%d.%s",J,str); fprintf(fp,"%-48s",str2);
		for(i=start; i < end && alignment[J][i]; i++){
		   char r=alignment[J][i];
		   if(r == '-') r='.';
		   fprintf(fp,"%c",r);
		} fprintf(fp,"\n");
	   } start=end; end=MINIMUM(Int4, AlignLen, end+50);
	   if(start==AlignLen) fprintf(fp,"//\n"); else fprintf(fp,"\n");
	} fflush(fp);
	for(J=0; J <= number; J++){ free(alignment[J]); } free(alignment);
	// alignment[0] allocated above....
	free(null[1]); free(null);
	delete []s;
}

Int4	PutFastaCMSA(FILE *fp,cma_typ cma, BooLean add_consensus)
// output cma file in fasta alignment format (so can be used by other routines
// == common language for alignments.
{
	Int4	d,i,j,lenM,LenStr,AlignLen,J,m,s,t,end;
	char	***alignment,*ps;
	st_type	S=SitesCMSA(cma);
	a_type	A=SeqSetA(DataCMSA(cma));
	e_type	E;
	BooLean	**null;
	gss_typ& gss=*gssCMSA(cma);
	Int4	site,number=gss.NumSeq(),buffer_size,*MaxInsert=0;

	//******************************************************************
	//*********** Get alignment and alignment template information: ****
	//******************************************************************
	AlignLen = 0;
	buffer_size = gss.MaxTrueSqLen();
	// length, which is needed to ELIMINATE GAPS OUTSIDE OF BLOCKS.
	for(J=1; J <= number; J++){
		E = gss.TrueSeq(J); 
		buffer_size+=gss.NumIns(J)+gss.NumDel(J); 
	} null = NullCMSA(cma);
	NEWPP(alignment,nBlksCMSA(cma)+3,char);
	for(m=1; m <= nBlksCMSA(cma); m++){
	   lenM = SiteLen(m,S);
	   NEWP(alignment[m],number+3,char);
	   for(J=1; J <= number; J++){ 
		NEW(alignment[m][J],buffer_size+3,char); 
	   	//*** 2. Get corresponding gapped alignable segments.
		LenStr=gss.Region(J,alignment[m][J],SitePos(m,J,1,S),lenM);
	   }
	   //// 3. Align gapped sequences.
	   AlignLen+=cmsa_fix_alignment(lenM,number,alignment[m],null[m]);
#if 0
	   for(J=1; J <= number; J++){ 
		char c;
		for(s=0; (c=alignment[m][J][s])!=0; s++){
		   if(c=='.') alignment[m][J][s]='-';
		}
	   } // alignment[0] allocated from cmsa_fix_alignment...
#endif
	} 
	//**********************************************************************
	//*** Find the distance from end of blk m and beginning of blk m+1. ****
	//**********************************************************************
	NEW(MaxInsert,nBlksCMSA(cma)+3,Int4);
	for(J=1; J <= number; J++){ 
	   Int4 fakesite = SitePos(1,J,1,S);
	   Int4 truesite = gss.TruePos(J,fakesite);
	   s=truesite;
	   if(gss.IsDeleted(J,fakesite)) s++; 
	   for(m=1; m < nBlksCMSA(cma); m++){
		for(ps=alignment[m][J]; *ps != 0; ps++) if(isalpha(*ps)) s++; 
		site = SitePos(m+1,J,1,S); 
		if(IsDeletedCMSA(J,site,cma)){ site = gss.TruePos(J,site) + 1; }
		else site = gss.TruePos(J,site);

		Int4 gap_insert=0;
	   	while(s < site){ s++; gap_insert++; }
		MaxInsert[m]= MAXIMUM(Int4,MaxInsert[m],gap_insert);
	   }
	} 
	//**********************************************************************
	//********************* Print Consensus ********************************
	//**********************************************************************
#if 0	// need to add gaps!!
	if(add_consensus){
	   e_type  cE=MkConsensusCMSA(cma);
	   fprintf(fp,">"); PutSeqInfo(fp,cE);
	   for(s=m=1; m <= nBlksCMSA(cma); m++){
	     // fprintf(stderr,"MaxInsert[%d]=%d\n",m,MaxInsert[m]);
	     for(i=1; i <= LengthCMSA(m,cma); i++){
		char res=AlphaChar(ResSeq(s,cE),A);
		fprintf(fp,"%c", res); s++;
	     }
	     for(i=1; i <= MaxInsert[m]; i++) fprintf(fp,"-");
	   } NilSeq(cE);
	   fprintf(fp,"\n\n");
	}
#endif
	//*************************************************************************
	//************************* Print input sequences *************************
	//*************************************************************************
	Int4 fakesite,truesite;
	for(J=1; J <= number; J++){ 
	   E = gss.TrueSeq(J); 
	   fakesite = SitePos(1,J,1,S);
	   truesite = gss.TruePos(J,fakesite);
	   if(gss.IsDeleted(J,fakesite)){
		// if start site in FakeSeq is deleted then corresponding position 
		// in TrueSeq is just before the deletion. So increment s.
	   	s=truesite + 1;
	   } else s=truesite;

	   if((OffSetSeq(E)+s-1) >= 0) AdjustOffsetSeq(s-1,E);	// add the difference to the offset.
	   fprintf(fp,">"); PutSeqInfo(fp,E);	// Print sequence id.

	   for(m=1; m <= nBlksCMSA(cma); m++){
	        fprintf(fp,"%s",alignment[m][J]);
	   	if(m < nBlksCMSA(cma)){
		    for(ps=alignment[m][J]; *ps != 0; ps++) if(isalpha(*ps)) s++;  
		    site = SitePos(m+1,J,1,S);	// site of next match.
		    // if start site in FakeSeq is deleted then corresponding position 
		    // in TrueSeq is just before the deletion and needs to be printed...
		    if(IsDeletedCMSA(J,site,cma)){ site = gss.TruePos(J,site) + 1; }
		    else site = gss.TruePos(J,site);
		    Int4 gap_insert=0;
	   	    while(s < site){
			char res=AlphaChar(ResSeq(s,E),A); res=tolower(res);
			fprintf(fp,"%c",res); s++; gap_insert++;
		    } 
	   	    while(gap_insert < MaxInsert[m]){ fprintf(fp,"-"); gap_insert++; }
		}
	   } fprintf(fp,"\n\n");
	} 
	// ************************************************************
	// ******************* Deallocate memory **********************
	free(MaxInsert);
	for(t=1; t <= nBlksCMSA(cma); t++) {
	   for(J=0; J <= number; J++){ free(alignment[t][J]); } free(alignment[t]);
	   free(null[t]); 
	} free(null); free(alignment);
	return AlignLen;
}

cma_typ	ReadCMSA2(char *filename,a_type A)
{
	cma_typ	cma;
	FILE *fp = open_file(filename,"","r"); 
	cma=ReadCMSA(fp,A); fclose(fp); return cma;
}

typedef struct cma_list_type {
	cma_typ	cma;
	char	*status;
        cma_list_type *next;          /* next[i] is successor of i in list */
};

cma_typ	*MultiReadCMSA(FILE *fp,Int4 *Number,a_type A)
{ return MultiReadCMSA(fp,Number,0,A); }

cma_typ	*MultiReadCMSA(FILE *fp,Int4 *Number,char ***Status,a_type A)
{
	cma_list_type *cma_list=NULL,*tmp;
	cma_typ	cma,*CMA; 
	Int4	i=0;
	char	*status;
	if(Status){
	  while((cma=ReadCMSA(fp,&status,A)) != NULL){
		NEW(tmp,1,cma_list_type); 
		tmp->status = status;
		tmp->next=cma_list; tmp->cma=cma; cma_list=tmp; i++;
	  } *Number=i;
	  NEW(CMA,i+3,cma_typ);
	  char **str; NEWP(str,i+3,char); *Status=str;
	  for(Int4 j=i; j > 0; j--){ 
		CMA[j]=cma_list->cma; cma_list->cma=NULL;
		str[j] = cma_list->status;
		tmp=cma_list; cma_list=tmp->next; free(tmp);
	  }
	} else {
	  while((cma=ReadCMSA(fp,A)) != NULL){
		NEW(tmp,1,cma_list_type); 
		tmp->next=cma_list; tmp->cma=cma; cma_list=tmp; i++;
	  } *Number=i;
	  NEW(CMA,i+3,cma_typ);
	  for(Int4 j=i; j > 0; j--){ 
		CMA[j]=cma_list->cma; cma_list->cma=NULL;
		tmp=cma_list; cma_list=tmp->next; free(tmp);
	  }
	} return CMA;
}

cma_typ	ReadCMSA(FILE *fp,a_type A) { return ReadCMSA(fp,0,A); }

cma_typ	ReadCMSA(FILE *fp,char **column_status,a_type A)
// read in a cma_typ object.
{
	cma_typ cma;
	char	c,**null,name[200];
	gss_typ *gss;
	Int4	i,J,m,s,t,**pos,number,depth=0,nblks;
	Int4	N,open,extend,left,right,err,*len_elem;
        double	Pernats;

	assert(A != NULL);
#if 0
<Fullseq(200)={234(2),544(4), ... };
>b_adapt_drome g434902|emb|CAA53509| (X75910) beta-adaptin 1 [Drosophi
VGKDVSALFPDVVNCMQTDNLELKKLVYLYLMNYAKSQPDMAIM

>xyz_seq

>.
#endif
        unsigned short  *FullRpts=0;      // number of repeats in each full seq.
	e_type	*FullE=0;
	Int4	FullN=0;
	dom_typ *dom=0;

	do { c=fgetc(fp); } while(isspace(c)); ungetc(c,fp);
	err=fscanf(fp,"<FullSeqs(%d)%c",&N,&c);
	if(err == 2 && c == '='){ // if full data set seqs are provided...
	   FullN=N;
	   NEW(FullE, N+3, e_type);
	   NEW(FullRpts, N+3, unsigned short);
	   unsigned short  *size;
	   NEW(size, N+3, unsigned short);
	   UInt4 rpts,x;
	   err=fscanf(fp,"{%u(%u)",&x,&rpts); assert(err == 2);
	   assert(rpts > 0 && x > 0);
	   FullRpts[1] = rpts; size[1] = x;
	   for(s=2; s <= N; s++){
		err=fscanf(fp,",%u(%u)",&x,&rpts); assert(err == 2);
	   	assert(rpts > 0 && x > 0);
	   	FullRpts[s] = rpts; size[s] = x;
	   } err=fscanf(fp,"};\n"); assert(err == 0);
	   for(s=1; s <= N; s++){
		FullE[s] = ReadSeq(fp, s, size[s], A);
	   	assert(FullE[s]); 
	   } free(size);
	   err=fscanf(fp,">.\n"); assert(err == 0);
	   // Check for domains information.
	   if((c=getc(fp)) == '<'){ ungetc(c,fp); dom = new dom_typ(fp); }
	   else ungetc(c,fp);
	}
#if 0
        err=fscanf(fp,"[%d_(%d)=%[^(](%d){go=%d,gx=%d,pn=%lf,lf=%d,rf=%d}:",
		&depth,&nblks,name,&number,&open,&extend,&Pernats,&left,&right);
	// if(err!=9) { char str[200]; fscanf(fp,"%[^ ]",str); fprintf(stderr,"%s\n",str); }
	if(err!=9) return NULL;
	assert(err==9);
	fscanf(fp,"\n");
#else
	char	str[1003],set_mode=' ';
	double	alpha=-1;
	Int4	A0=-1,B0=-1;
	if(fgets(str,1000,fp) == NULL) return NULL;	// reads up to and including '\n'.
        err=sscanf(str,"[%d_(%d)=%[^(](%d){go=%d,gx=%d,pn=%lf,lf=%d,rf=%d};BPPS=(%c,%d,%d:%lf):",
		&depth,&nblks,name,&number,&open,&extend,&Pernats,&left,&right,
			&set_mode,&A0,&B0,&alpha);
	if(err!=13){
	   alpha=-1; A0=B0=-1; set_mode=' ';
           err=sscanf(str,"[%d_(%d)=%[^(](%d){go=%d,gx=%d,pn=%lf,lf=%d,rf=%d}:",
		&depth,&nblks,name,&number,&open,&extend,&Pernats,&left,&right);
	   if(err!=9) return NULL;
	} 
#endif

	NEW(len_elem,nblks+2,Int4);
	NEWP(null,nblks+2,char);
    if(column_status){
	if(nblks != 1) print_error("Fatal error: input alignment is block based!");;
	assert(fscanf(fp,"(%d)",&len_elem[1])==1);
	NEW(null[1],len_elem[1]+2,char); null[1]++;

	char *status; 
	NEW(status,len_elem[1]+4,char); 
	*column_status = status; status++;

	assert(fscanf(fp,"%[!^.*]",status)==1);
	// printf("%s\n",status);
	assert(len_elem[1]==strlen(status));
	for(i=0; i < len_elem[1]; i++){
		if(status[i]=='.') null[1][i]=TRUE;
		else if(status[i]=='^') null[1][i]=TRUE;
		else null[1][i]=FALSE;
	} null[1]--;
	fscanf(fp,"\n"); fscanf(fp,"\n");
    } else {
	for(m=1; m <= nblks; m++){
	   assert(fscanf(fp,"(%d)",&len_elem[m])==1);
	   NEW(null[m],len_elem[m]+2,char); null[m]++;
	   assert(fscanf(fp,"%[!^.*]",null[m])==1);
	   assert(len_elem[m]==strlen(null[m]));
	   for(i=0; i < len_elem[m]; i++){
		if(null[m][i]=='.') null[m][i]=TRUE;
		else if(null[m][i]=='^') null[m][i]=TRUE;
		else null[m][i]=FALSE;
	   } null[m]--;
	   fscanf(fp,"\n");
	} fscanf(fp,"\n");
    }
	
	// 2. scan aligned sequences:
	gsq_typ	**gsq = new gsq_typ*[number+2];
	NEWP(pos,number+2,Int4);
	for(J=1; J <= number; J++){ 
	   rff_typ *rff=0;
	   gsq[J] = new gsq_typ[1];
	   NEW(pos[J],nblks+2,Int4);
	   gsq[J]->read(fp,nblks,len_elem,pos[J],A);
#if 0	// call this function recursively here to read tree...
	   cma->next[J] = ReadCMSA(fp,A,depth+1);
#endif
// if(J==164) { gsq[J]->Put(stderr,A); gsq[J]->Put(stderr,60,A); }
	}
	Int4	depth2=0;
        if(fscanf(fp,"\n_%d].\n",&depth2) != 1) print_error("ReadCMSA(): cma input error"); 
	if(depth != depth2){
	     fprintf(stderr,"Warning: CMA file level inconsistency! (%d != %d)\n",
		depth,depth2);
	}
	// assert(depth == depth2);
	// 3. create sequence set.
	gss = new gss_typ[1];
	gss->FromArray(number,gsq,open,extend,Pernats,left,right,name,A);
	cma=MakeCMSA(MakeSites(nblks, len_elem, *gss),null);
	cma->Level=depth;
	if(alpha > 0.0 && alpha <= 1.0) { 
		cma->alpha=alpha; cma->A0 = A0; cma->B0=B0; cma->set_mode = set_mode;
	}
#if 1	// NEW: FullSeq. 
	if(FullE){ 
	  ss_type FullSeq=Array2SeqSet(FullE,FullN,"full_seqs",A);
	  AddFullCountsCMSA(FullSeq, FullRpts,cma);
	  if(dom) AddDomainsCMSA(dom, cma);
	}
#endif
	delete [] gss;// A copy was made.
	for(J=1; J <= number; J++){ 
		for(m=1; m <= nblks; m++) AddSiteCMSA(m,J,pos[J][m],cma);
		free(pos[J]);
	} free(pos);
	for(t=1; t <= nblks; t++) { free(null[t]); }
	free(null); free(len_elem);
	return cma;
}

void    PutGoodRptsCMSA(FILE *fp,Int4 min_rpt,Int4 min_spacing,cma_typ cma)
{
        gss_typ *gss=gssCMSA(cma);
        e_type  E;
	BooLean	*GoodSpace;
	unsigned short *SqNum;
        Int4    Sq,J,pos[5],start,end,repeats;
        Int4    number=gss->NumSeq(),N=NumSeqsCMSA(cma);
        char    id0[202],id[202];

        if(nBlksCMSA(cma) != 1) print_error("PutGoodRptsCMSA( ) input error");
	NEW(SqNum,N+3,unsigned short);
	NEW(GoodSpace,N+3,BooLean);
	// 1. Find numbers of repeats within min_spacing from another repeat.
        for(id0[0]=0,Sq=end=repeats=0,J=1; J <= number; J++) {
           E = gss->TrueSeq(J);
           StrSeqID(id,200,E);
           if(strncmp(id,id0,200) != 0){ repeats=0; Sq++; strcpy(id0,id); }
	   SqNum[J]=Sq; 
           PosSiteCMSA(1,J,pos,cma);
           start = gss->TrueSite(J,pos[1]); // offset is added here...
           if(repeats!=0) {
		if((start-end) <= min_spacing) GoodSpace[J]=TRUE; // else GoodSpace[J]=FALSE;
	   }
           end = pos[1]+LengthCMSA(1,cma)-1; end = gss->TrueSite(J,end);
           repeats++;
        }
	// 2. Find the number of good repeats in each Full sequence.
	BooLean	good;
	unsigned short *GoodRpts;
	NEW(GoodRpts, Sq+3, unsigned short);
        for(J=1; J <= number; J++) {
	   good=FALSE;
	   if(SqNum[J-1]==SqNum[J] && GoodSpace[J-1]) good=TRUE;
	   if(SqNum[J]==SqNum[J+1] && GoodSpace[J]) good=TRUE;
	   if(good) GoodRpts[SqNum[J]]++; 
	} free(GoodSpace);
        BooLean *skip; NEW(skip,N+3,BooLean);
        for(J=1; J <= number; J++) {
	   Sq = SqNum[J];
	   if(GoodRpts[SqNum[J]] < min_rpt) skip[J]=TRUE; else skip[J]=FALSE;
	} PutSelectCMSA(fp,skip,cma); free(skip); free(SqNum); free(GoodRpts);
}

Int4	*SortByQueryOneCMSA(cma_typ cma)
// return a list of aligned sequences sorted by score to query (first sequence).
{
	Int4	b,i,j,N=NumSeqsCMSA(cma);
	Int4	score,nblk=nBlksCMSA(cma);
	Int4	s,si,sj,pos[3];
        a_type  AB=AlphabetCMSA(cma);
	unsigned char	*isq,*jsq;

	dh_type dH=dheap(N+2,4);
	isq = SeqPtrCMSA(1,cma);
	for(j=2; j <= N; j++){
	      jsq = SeqPtrCMSA(j,cma);
	      for(score=0,b=1; b <= nblk; b++){
		PosSiteCMSA(b,1,pos,cma); si=pos[1];
		PosSiteCMSA(b,j,pos,cma); sj=pos[1];
		for(s=1; s <= LengthCMSA(b,cma); s++,si++,sj++){
			if(isq[si] > 0) score+=valAlphaR(isq[si],jsq[sj],AB);
		}
	      } insrtHeap(j,(keytyp)-score,dH);
	}
	Int4	*list;
	NEW(list, NumSeqsCMSA(cma)+2,Int4);
	list[1]=1;
	for(i=2; !emptyHeap(dH); ){
		assert((j=delminHeap(dH)) != 0);
		list[i]=j; i++;
	} return list;
}

Int4	PutSelectCMSA(FILE *fp,BooLean *skip, cma_typ cma)
{
	if(nBlksCMSA(cma) == 1) return PutSelectOneCMSA(fp,skip,cma);
	else return PutSelectCMSA0(fp,skip,cma);
}

Int4	PutSelectCMSA0(FILE *fp,BooLean *skip, cma_typ cma)
{ return PutSelectCMSA0(fp,skip, FALSE, cma); }

void	PutCsqWithCMSA(FILE *fp,cma_typ cma) { PutSelectCMSA0(fp,0, TRUE, cma); }

Int4	PutSelectCMSA0(FILE *fp,BooLean *skip, BooLean put_csq, cma_typ cma)
// put in aligned sequence format (so can be used by other routines
// == common language for alignments.
// == Psiblast too where master sequence defines a single block template.
{
	Int4	d,i,j,lenM,LenStr,AlignLen,J,J0,m,s,t,end;
	// Int4	netlen;
	Int4	insertions,deletions,matches,interblks,endseqs;
	st_type	S=SitesCMSA(cma);
	a_type	A=SeqSetA(DataCMSA(cma));
	e_type	fakeE,E;
	BooLean	**null;
	gss_typ& gss=*gssCMSA(cma);
	Int4	site,number=gss.NumSeq(),buffer_size,depth=0;
	Int4	NetNum;

	Int4	Num;
	if(skip == 0) Num = number;
	else for(Num=0,J=1; J <= number; J++) if(!skip[J]) Num++;
	if(cma->FullSeq) assert(1 == PutFullSeqCMSA(fp,skip, cma));
	NetNum=Num;
	if(put_csq) NetNum++;
        fprintf(fp,"[%d_(%d)=%s(%d){go=%d,gx=%d,pn=%.1f,lf=%d,rf=%d}",
		cma->Level,nBlksCMSA(cma),NameCMSA(cma),NetNum,
		gss.GapOpen(),gss.GapExtend(),gss.PerNats(),
        	gss.LeftFlank(),gss.RightFlank());
	if(cma->alpha > 0.0 && cma->alpha <= 1.0){
            fprintf(fp,";BPPS=(%c,%d,%d:%.4f):\n",cma->set_mode,cma->A0,cma->B0,cma->alpha);
	} else { fprintf(fp,":\n"); }
	// 2a. Print column indicators:
#if 0	// BPPS pattern positions here....
	if(sstP != 0){
	   fprintf(fp,"(%d)",LengthCMSA(1,cma));
	   for(m=1; m <= LengthCMSA(1,cma); m++){
		if(sstP[m]) fprintf(fp,"!"); else fprintf(fp,"*");
	   } fprintf(fp,"\n");
	} 
#endif
	// null = NullCMSA(cma); null[m][s] == ?
	for(m=1; m <= nBlksCMSA(cma); m++){
	   fprintf(fp,"(%d)",LengthCMSA(m,cma));
	   fm_type fm=ModelCMSA(m,cma);
	   for(s=1; s <= LengthCMSA(m,cma); s++){
		if(NullSiteFModel(s,fm)) fprintf(fp,".");
		else fprintf(fp,"*");
	   } fprintf(fp,"\n");
	} fprintf(fp,"\n");
	fprintf(fp,"\n");
	//******  Print consensus sequence here:
	if(put_csq){
	  fprintf(fp,"$0=%d(%d):\n>%s consensus seq\n",
			TotalLenCMSA(cma),TotalLenCMSA(cma),NameCMSA(cma));
	  fprintf(fp,"{()");
	  if(nBlksCMSA(cma) == 1){
	     char    *csq=ConsensusSeqCMSA(cma); fprintf(fp,"%s()",csq+1); free(csq);
	  } else {
	    for(m=1; m <= nBlksCMSA(cma); m++){
	     char    *csq=ConsensusSeqCMSA(m, cma);
	     fprintf(fp,"%s()",csq); free(csq);
	    }
	  }fprintf(fp,"}*\n\n");
	}
	create_gsq_cmsa(cma);	// create gsq objects for cma if not present...
	for(J0=J=1; J <= number; J++){ 
	   if(skip && skip[J]) continue;
	   gsq_typ *gsq=gsqCMSA(J,cma);
	   assert(gsq != 0);	// this should not happen...
	   Int4    *sites=GetPosSitesCMSA(J,cma);
	   gsq->Put_cma_format(fp,J0,nBlksCMSA(cma),sites,LengthsCMSA(cma),A);
	   free(sites);
	   J0++;  // always number consecutively 
	} fprintf(fp,"_%d].\n",cma->Level); 
	return NetNum;
}

Int4    create_gsq_cmsa(cma_typ cma)
// for each gsq == 0 create a gsq.
{
        // a_type AB=AlphabetCMSA(cma);
        Int4    N=0,sq,blk,*sites,len,s;
        for(sq=1; sq <= NumSeqsCMSA(cma); sq++){
           gsq_typ *gsq=gsqCMSA(sq,cma);
           if(gsq==0){
#if 1	// seq in alingment?
	      if(!SeqIsInCMSA(sq,cma)){
		gsq = new gsq_typ[1]; 
		gsq->initialize(TrueSeqCMSA(sq,cma));
	      } else 
#else
	      assert(SeqIsInCMSA(sq,cma));
#endif
	      {
                gsq = new gsq_typ[1]; N++;
                sites=GetPosSitesCMSA(sq,cma);
                // fprintf(stderr,"%d",sq);
                gsq->initialize(nBlksCMSA(cma),LengthsCMSA(cma),sites,TrueSeqCMSA(sq,cma));
                // gsq->Put_cma_format(stderr,sq,nBlksCMSA(cma),sites,LengthsCMSA(cma),AB);
                ReplaceCMSA(sq,gsq,cma); // replace sequence sq in CMSA; vacates sites.
                for(blk=1; blk <= nBlksCMSA(cma); blk++) AddSiteCMSA(blk,sq,sites[blk],cma);
	        free(sites);
	      }
           }
        } return N;
}

Int4	PutSelectOneCMSA(FILE *fp,BooLean *skip, cma_typ cma)
{ return PutSelectOneCMSA(fp,skip,0,cma); }

Int4 PutSelectOneCMSA(FILE *fp,BooLean *skip, Int4 *sortedlist, cma_typ cma)
// put in aligned sequence format (so can be used by other routines
// == common language for alignments.
// == Psiblast too where master sequence defines a single block template.
// check that only one block before calling.  Low memory requirement.
// IMPORTANT: PutSelectCMSA0( ) SHOULD ALSO BE MODIFIED IN THE SAME WAY.
{
	st_type	S=SitesCMSA(cma);
	a_type	A=SeqSetA(DataCMSA(cma));
	e_type	fakeE,E;
	gss_typ& gss=*gssCMSA(cma);
	Int4	site,number=gss.NumSeq(),depth=0;
	Int4	m,i,d,j,lenM,LenStr,J,J0,s,t,end,netlen,tmp;
	// netlen Used to store fakeE effective length, which is 
	// needed to ELIMINATE GAPS OUTSIDE OF BLOCKS.
	Int4	Num;

	// assert(nBlksCMSA(cma) == 1);
	if(skip == 0){
	    if(sortedlist){
		for(Num=0, i=1; i <= number; i++) if(sortedlist[i] > 0) Num++;
	    } else Num = number;
	} else if(sortedlist){
	    for(Num=0,i=1; i <= number; i++){
		J=sortedlist[i]; 
		if(J == 0 || skip[J]) continue; else Num++;
	    }
	} else { for(Num=0,J=1; J <= number; J++) if(!skip[J]) Num++; }
	if(cma->FullSeq) assert(1 == PutFullSeqCMSA(fp,skip,cma));
	depth=cma->Level;
        fprintf(fp,"[%d_(%d)=%s(%d){go=%d,gx=%d,pn=%.1f,lf=%d,rf=%d}",
		depth,nBlksCMSA(cma),NameCMSA(cma),Num,
		gss.GapOpen(),gss.GapExtend(),gss.PerNats(),
        	gss.LeftFlank(),gss.RightFlank());
	if(cma->alpha > 0.0 && cma->alpha <= 1.0){
            fprintf(fp,";BPPS=(%c,%d,%d:%.4f):\n",cma->set_mode,cma->A0,cma->B0,cma->alpha);
	} else { fprintf(fp,":\n"); }
	for(m=1; m <= nBlksCMSA(cma); m++){
	   fprintf(fp,"(%d)",LengthCMSA(m,cma));
	   fm_type fm=ModelCMSA(m,cma);
	   for(s=1; s <= LengthCMSA(m,cma); s++){
		if(NullSiteFModel(s,fm)) fprintf(fp,".");
		else fprintf(fp,"*");
	   } fprintf(fp,"\n");
	} fprintf(fp,"\n");
	// 2b. Print aligned sequences:
	create_gsq_cmsa(cma);	// create gsq_typ for cma if not present...
	Int4	I0,II;
	for(J0=I0=II=1; II <= number; II++){ // replaces: for(J0=J=1; J <= number; J++)
	   if(sortedlist){
		if(cma->FullSeq) print_error("sortedlist option not available w/ FullSeq!");
		if(sortedlist[II] == 0) continue;
		assert(sortedlist[II] <= number && sortedlist[II] > 0);
		J=sortedlist[II];
	   	if(skip && skip[J]) continue;
	   } else {
		J=II;
	   	if(skip && skip[J]) continue;
	   }
	   gsq_typ *gsq=gsqCMSA(J,cma);
	   Int4    *sites=GetPosSitesCMSA(J,cma);
	   gsq->Put_cma_format(fp,J0,nBlksCMSA(cma),sites,LengthsCMSA(cma),A);
	   free(sites);
	   J0++;  // always number consecutively 
	} fprintf(fp,"_%d].\n",depth); 
	// Deallocate memory... may want to do this above for recursion
	return Num;
}

double  AveRelEntropyCMSA(BooLean *skip, cma_typ cma)
{
        Int4    t,k,sq,r1,n;
        double  sum,re,cnt[50],tcnt;
        a_type A = AlphabetCMSA(cma);
	double	*freq=tFreqSeqSet(TrueDataCMSA(cma));

     for(sum=0.0,n=0,t=1; t <= nBlksCMSA(cma); t++){
        for(k=1; k <= LengthCMSA(t,cma); k++){
	  n++;
          for(tcnt=0.0,r1=0; r1 <= nAlpha(A); r1++) cnt[r1]=0.0;
          for(sq=1; sq <= NumSeqsCMSA(cma); sq++){
              if(skip!=0 && skip[sq]) continue;
              r1=ResidueCMSA(t,sq,k,cma);
	      if(r1){ cnt[r1]+=1.0; tcnt +=1.0; }
          }
	  double p,q;
          for(re=0.0,r1=1; r1 <= nAlpha(A); r1++){
		p = cnt[r1]/tcnt; q = freq[r1];
		if(p > 0.0) re += p * log(p/q);
	  } sum += re;
        }
     } return (sum/(double)n);
}

double	AveRelEntropyCMSA(set_typ set, cma_typ cma)
{
	BooLean *skipseq;
	NEW(skipseq,NumSeqsCMSA(cma)+3,BooLean);
	for(Int4 sq=1; sq <=NumSeqsCMSA(cma); sq++) if(!MemberSet(sq,set)) skipseq[sq]=TRUE;
	double d=AveRelEntropyCMSA(skipseq,cma); free(skipseq);
	return d;
}

double  AvePercentIdentityCMSA(register set_typ set, register unsigned char *qsq,
	register cma_typ cma)
// ignores ambiguity residues
{
     register UInt4	i,sq;
     register double		sum,total;
     register unsigned char	*ssq;
     register char		**R=AlphaR(AlphabetCMSA(cma));

     assert(nBlksCMSA(cma) == 1);
     for(sum=total=0.0,sq=NumSeqsCMSA(cma); sq > 0; sq--){
	  if(!MemberSet(sq,set)) continue;
          ssq=SeqSeqSet(sq,DataCMSA(cma));
          for(i=LengthCMSA(1,cma); i > 0; i--){
		if(ssq[i] > 0){
			total+=1.0;
			sum+=R[qsq[i]][ssq[i]];
			// if(qsq[i] == ssq[i]) sum++;
		}
          }
     } return (100.0*sum)/total;
}

double  ResidueDiversityCMSA(register Int4 *sqid, register cma_typ cma)
// ignores ambiguity residues
{
     register sst_typ   set;
     register Int4      k,i,sum;
     register unsigned char *seq,r;

     assert(sqid); assert(nBlksCMSA(cma) == 1);
     for(sum=0,k=LengthCMSA(1,cma); k > 0; k--){
          for(set=0,i=0; sqid[i]!=0; i++){
                seq=SeqSeqSet(sqid[i],DataCMSA(cma));
		r=seq[SitePos(1,sqid[i],1,SitesCMSA(cma))+k-1];
		set = UnionSset(set,SsetLet(r));
                // set = UnionSset(set,SsetLet(seq[SitePos(1,sqid[i],1,SitesCMSA(cma))+k-1]));
          } for(i=nAlpha(AlphabetCMSA(cma)); i > 0; i--) if(MemSset(i,set)) sum++;
     } 
     // double diversity = (double)sum/(double) LengthCMSA(1,cma);
     // fprintf(stderr,"diversity = %g\n",diversity);
     // return diversity;
     return (double)sum/(double) LengthCMSA(1,cma);
}

double	ResidueDiversityNotCMSA2(set_typ set, cma_typ cma)
{
	if(NumSeqsCMSA(cma) >= 10000) print_error("ResidueDiversityCMSA() limit exceeeded");
	Int4	sq,sqid[10003],i;
	for(i=0,sq=1; sq <=NumSeqsCMSA(cma); sq++){
		if(!MemberSet(sq,set)){ sqid[i]=sq; i++; } 
	} sqid[i]=0;
	return ResidueDiversityCMSA(sqid,cma);
}

double	ResidueDiversityCMSA2(set_typ set, cma_typ cma)
{
#if 1
	Int4	sq,*sqid,i;
	NEW(sqid,NumSeqsCMSA(cma)+3,Int4);
	for(i=0,sq=1; sq <=NumSeqsCMSA(cma); sq++){
		if(MemberSet(sq,set)){ sqid[i]=sq; i++; } 
	} sqid[i]=0;
	double d=ResidueDiversityCMSA(sqid,cma); free(sqid);
	return d;
#else
	if(NumSeqsCMSA(cma) >= 10000) print_error("ResidueDiversityCMSA() limit exceeeded");
	Int4	sq,sqid[10003],i;
	for(i=0,sq=1; sq <=NumSeqsCMSA(cma); sq++){
		if(MemberSet(sq,set)){ sqid[i]=sq; i++; } 
	} sqid[i]=0;
	return ResidueDiversityCMSA(sqid,cma);
#endif
}

double	ResidueDiversityCMSA(set_typ set, cma_typ cma)
// WARNING: fix 
{
	BooLean skipseq[10003];
	if(NumSeqsCMSA(cma) >= 10000) print_error("ResidueDiversityCMSA() limit exceeeded");
	for(Int4 sq=1; sq <=NumSeqsCMSA(cma); sq++){
	     if(MemberSet(sq,set)) skipseq[sq]=FALSE; else skipseq[sq]=TRUE;
	} return ResidueDiversityCMSA(skipseq,cma);
}

void	PutDiffSetsCMSA(char *filename, Int4 set_id,set_typ SubSet,set_typ SuperSet,
		cma_typ cma)
{
	char str[100];
	a_type	AB = AlphabetCMSA(cma);
	if(set_id < 0){ sprintf(str,".set%di",-set_id); }
	else { sprintf(str,".set%d",set_id); }
        FILE *fp = open_file(filename,str,"w");
	for(Int4 sq=1; sq <=NumSeqsCMSA(cma); sq++){
	    if(SubSet == 0){
	       if(MemberSet(sq,SuperSet)) PutSeq(fp,TrueSeqCMSA(sq,cma),AB);
	    } else if(MemberSet(sq,SuperSet) && !MemberSet(sq,SubSet)){
		PutSeq(fp,TrueSeqCMSA(sq,cma),AB);
	    }
	} fprintf(fp,"\n"); fclose(fp);
}

void    PutSelectCMSA(FILE *fp,FILE *fp2, const char *name1,const char *name2,
	set_typ Set, cma_typ cma)
{
	BooLean *skip,*notskip;
	char	old_name[200];
	char	new_name[200];
	if(name1 || name2) strcpy(old_name,NameCMSA(cma));
	NEW(skip,NumSeqsCMSA(cma)+3,BooLean);
	NEW(notskip,NumSeqsCMSA(cma)+3,BooLean);
	for(Int4 sq=1; sq <=NumSeqsCMSA(cma); sq++){
	   	if(MemberSet(sq,Set)){ skip[sq]=FALSE; notskip[sq]=TRUE; }
		else { skip[sq]=TRUE; notskip[sq]=FALSE; }
	}
	if(fp){
	    if(name1){
		 strcpy(new_name,name1);
		 ReNameCMSA(new_name,cma);
	    } PutSelectCMSA(fp,skip,cma); 
	}
	if(fp2){ 
	    if(name2){
		 strcpy(new_name,name2);
		 ReNameCMSA(new_name,cma);
	    } PutSelectCMSA(fp2,notskip,cma);
	} free(skip); free(notskip);
	if(name1 || name2) ReNameCMSA(old_name,cma);
}

void    PutSelectCMSA(FILE *fp,FILE *fp2, set_typ Set, cma_typ cma)
{ PutSelectCMSA(fp,fp2,0,0,Set,cma); }

void	PutSeqsInSetCMSA(FILE *fp,set_typ set, cma_typ cma)
{
	gss_typ& gss=*gssCMSA(cma);
	a_type	A=SeqSetA(DataCMSA(cma));
	for(Int4 J=1; J <= gss.NumSeq(); J++){
	  if(MemberSet(J,set)){ PutSeq(fp,gss.TrueSeq(J),A); }
	}
}

void	PutMergedCMSA(FILE *fp,unsigned short nsets,cma_typ *cma)
{
	set_typ	*set; NEW(set,nsets+3,set_typ);
	PutMergedCMSA(fp,nsets,set, cma, NULL); free(set);
}

void	PutMergedCMSA(FILE *fp,unsigned short nsets,set_typ *set, cma_typ *cma, sst_typ	*sstP)
{ PutMergedCMSA(fp,0, nsets,set, cma,sstP,0); }

void	PutMergedCMSA(FILE *fp,char *name, unsigned short nsets,set_typ *set, cma_typ *cma,
		sst_typ	*sstP)
{ PutMergedCMSA(fp,name,nsets,set,cma,sstP,0); }

void	PutMergedCMSA(FILE *fp,char *name, unsigned short nsets,set_typ *set, cma_typ *cma,
		sst_typ	*sstP,char Labeling)
// == Psiblast too where master sequence defines a single block template.
// check that only one block before calling.  Low memory requirement.
{
	e_type	fakeE,E;
	Int4	site,number,I,i,d,j,lenM,LenStr,J,J0,s,end;

	assert(nsets > 0);
	a_type	A=SeqSetA(DataCMSA(cma[1]));
	st_type	S;

	Int4	Num=0,depth=0;
	depth=cma[1]->Level;
	for(i = 1; i <= nsets; i++){
	    assert(nBlksCMSA(cma[i]) == 1);
	    if(LengthCMSA(1,cma[1]) != LengthCMSA(1,cma[i])){
		fprintf(stderr,"FATAL: lengthCMA inconsistency: \n\t cma[1](%s) = %d; cma[%d](%s) = %d\n",
			NameCMSA(cma[1]),LengthCMSA(1,cma[1]),i,
			NameCMSA(cma[i]),LengthCMSA(1,cma[i])); 
			assert(LengthCMSA(1,cma[1]) == LengthCMSA(1,cma[i])); 
			exit(1);
	    }
	    if(set[i]) assert(NumSeqsCMSA(cma[i]) < SetN(set[i]));
	    if(set[i]) Num += CardSet(set[i]);
	    else Num += NumSeqsCMSA(cma[i]);
	} number=Num;

	// if(cma->FullSeq) assert(J0 == PutFullSeqCMSA(fp,skip,cma));
	{
	  gss_typ& gss=*gssCMSA(cma[1]);
	  if(name){
            fprintf(fp,"[%d_(1)=%s(%d){go=%d,gx=%d,pn=%.1f,lf=%d,rf=%d}",
		depth,name,Num,gss.GapOpen(),gss.GapExtend(),
		gss.PerNats(), gss.LeftFlank(),gss.RightFlank());
	  } else {
	    assert(NameCMSA(cma[1]));
            fprintf(fp,"[%d_(1)=%s(%d){go=%d,gx=%d,pn=%.1f,lf=%d,rf=%d}",
		depth,NameCMSA(cma[1]),Num,gss.GapOpen(),gss.GapExtend(),
		gss.PerNats(), gss.LeftFlank(),gss.RightFlank());
	  }
	  if(cma[1]->alpha > 0.0 && cma[1]->alpha <= 1.0){
            fprintf(fp,";BPPS=(%c,%d,%d:%.4f):\n",
			cma[1]->set_mode,cma[1]->A0,cma[1]->B0,cma[1]->alpha);
	  } else { fprintf(fp,":\n"); }
	  if(sstP){	// if BPPS pattern positions passed in...
	     fprintf(fp,"(%d)",LengthCMSA(1,cma[1]));
	     for(s=1; s <= LengthCMSA(1,cma[1]); s++){
		if(sstP[s]) fprintf(fp,"!"); else fprintf(fp,"*");
	     } fprintf(fp,"\n");
	  } else {	// original code.
	    fprintf(fp,"(%d)",LengthCMSA(1,cma[1]));
	    fm_type fm=ModelCMSA(1,cma[1]);
	    for(s=1; s <= LengthCMSA(1,cma[1]); s++){
		if(NullSiteFModel(s,fm)) fprintf(fp,".");
		else fprintf(fp,"*");
	    } fprintf(fp,"\n\n");
	  }
	}

	for(J0=1,I = 1; I <= nsets; I++){
	 S=SitesCMSA(cma[I]); 
	 create_gsq_cmsa(cma[I]);	// create gsq objects for cma if not present...
	 for(J=1; J <= NumSeqsCMSA(cma[I]); J++){ 
	   if(nSites(1,J,S) == 0) continue;
	   if(set[I] && !MemberSet(J,set[I])) continue;
	   gsq_typ *gsq=gsqCMSA(J,cma[I]);
	   if(Labeling == 'L') LabelSeq(gsq->TrueSeq( ));
	   else if(Labeling == 'U') UnLabelSeq(gsq->TrueSeq( ));
	   Int4    *sites=GetPosSitesCMSA(J,cma[I]);
	   gsq->Put_cma_format(fp,J0,nBlksCMSA(cma[I]),sites,LengthsCMSA(cma[I]),A);
	   free(sites);
	   J0++;  // always number consecutively 
	 }
	} fprintf(fp,"_%d].\n",depth); 
}

void	PutInSetCMSA(FILE *fp,set_typ set, cma_typ cma)
{ set_typ s[3]; s[1]=set; cma_typ c[3]; c[1]=cma; PutMergedCMSA(fp,1,s,c,0); }

void	PutInSetCMSA(FILE *fp,set_typ set, sst_typ *sstP, cma_typ cma)
{ set_typ s[3]; s[1]=set; cma_typ c[3]; c[1]=cma; PutMergedCMSA(fp,1,s,c,sstP); }

void	LabelPutInSetCMSA(FILE *fp,set_typ set, cma_typ cma)
{ set_typ s[3]; s[1]=set; cma_typ c[3]; c[1]=cma; PutMergedCMSA(fp,0,1,s,c,0,'L'); }

void	UnLabelPutInSetCMSA(FILE *fp,set_typ set, cma_typ cma)
{ set_typ s[3]; s[1]=set; cma_typ c[3]; c[1]=cma; PutMergedCMSA(fp,0,1,s,c,0,'U'); }

