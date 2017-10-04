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

//double	num_pseudo_msa(double N) { return pow(N,(0.71 - 0.0036*N)); }
double	num_pseudo_msa(double N) { return sqrt(N); }

cma_typ MakeColinearMSA(st_type S, BooLean **null)
{cma_typ L=MkMSA(S);create_models_cmsa(L,null);fill_models_cmsa(L);return L;}

cma_typ MakeCMSA(e_type *ListE, Int4 N, char **operation, Int4 *start, Int4 nblks,
	Int4 *lengths, Int4 gapo, Int4 gapx, double pernats, Int4 leftflank, 
	Int4 rightflank, char *name, a_type A, e_type *FullE, unsigned short *FullR)
// e.g.,  cma=MakeCMSA(ListE, NumHits, operation, start, nblks,PrtnModelLengths(PM),
//            gapo,gapx,5,left,right,argv[1],PrtnModelA(PM),FullE,FullR);
{

	cma_typ cma;
        cma=MakeCMSA(ListE,N,operation,start,nblks,lengths,gapo,gapx,
				pernats,leftflank,rightflank,name,A);
        Int4 fullN,i;
        for(fullN=0,i=1; FullR[i] > 0; i++,fullN++) assert(FullE[i]);
        ss_type FullSeq = Array2SeqSet(FullE,fullN,"full_seqs",A);
        AddFullCountsCMSA(FullSeq,FullR,cma);
        for(Int4 s=1; s <= N; s++) free(operation[s]);
	return cma;
}

cma_typ MakeCMSA(e_type *ListE, Int4 N, char **operation, Int4 *start, cma_typ cma)
// create cma_typ from output from a database search 
// ListE = list of (sub)sequences found (numbered from 1..N);  
// operation = list of operation strings
// start = list of start sites
// cma = cma_typ used for search
{
	gss_typ *gss=gssCMSA(cma);
	return MakeCMSA(ListE, N, operation, start, nBlksCMSA(cma),
        	LengthsCMSA(cma), gss->GapOpen(), gss->GapExtend(), 
		gss->PerNats(), gss->LeftFlank(),gss->RightFlank(),
        	NameCMSA(cma), AlphabetCMSA(cma));
}

cma_typ MakeCMSA(e_type *ListE, Int4 N, char **operation, Int4 *start, Int4 nblks,
	Int4 *lengths, Int4 open, Int4 extend, double pernats, Int4 leftflank,
	Int4 rightflank, char *name, a_type A)
// nblks = number of blocks 
// 
{
	cma_typ	cma;
        ss_type data=Array2SeqSet(ListE,N,name,A);
        cma=EmptyCMSA(nblks,lengths,data,open,extend,pernats,leftflank,rightflank);
        // split up operations into substrings with start sites for each...
        gss_typ *gss=gssCMSA(cma);
	Int4    *newpos;
	NEW(newpos,nBlksCMSA(cma)+2,Int4);
	for(Int4 s=1; s <= N; s++){
	  e_type E = gss->TrueSeq(s);
          gsq_typ *gsq; gsq = new gsq_typ[1];
          Int4 trace_length=strlen(operation[s]);
          gsq->initialize(leftflank,rightflank,operation[s],trace_length,start[s],E,newpos);
          ReplaceCMSA(s,gsq,cma);
          for(Int4 t=1; t<=nBlksCMSA(cma); t++) AddSiteCMSA(t,s,newpos[t], cma);
	} free(newpos);
	return cma;
}

cma_typ MakeCMSA(st_type S, BooLean **null)
{ 	
	cma_typ cmsa;
	BooLean **null0;
	if(null==NULL) {
		NEWP(null0,nTypeSites(S) +2, BooLean);
		cmsa=make_cmsa(MakeColinearMSA(S,null0)); free(null0);
	} else cmsa=make_cmsa(MakeColinearMSA(S,null));
	return cmsa;
}

#if 0	// WARNING: Not yet tested or used in any way...
BooLean GetColResFreqCMSA(Int4 t, BooLean right, Int4 *pos,char *blocked, cma_typ cma)
// Int4    *AlwaysGetSiteFreq(st_type S,Int4 t,BooLean right,Int4 *pos,char *blocked)
// blocked and *pos are arrays the length of the number of sequences (plus 2).
{
	Int4	site_freq;

	site_freq[1]=AlwaysGetSiteFreq(S,t,right,location,blocked);
	return site_freq;
}
#endif

void    PutMinimalSeqCMA(FILE *fp, Int4 sq, cma_typ cma)
// print out a concensus sequence cma file
{
 {      
	if(sq < 1 || sq > NumSeqsCMSA(cma)) return;
	if(nBlksCMSA(cma) > 1) print_error("PutMinimalSeqCMA( ) requires only 1 blk");

        gss_typ *gss = gssCMSA(cma);
        char	*aln;
        Int4    len=LengthCMSA(1,cma),pos[5],max_id_len=30;
	Int4	x = len + gss->NumIns(sq) + gss->NumDel(sq);
	NEW(aln,x+55,char);
        e_type sqE= TrueSeqCMSA(sq,cma);
        StrSeqID(aln, max_id_len, sqE);
        fprintf(fp, "[0_(1)=%s(1){go=10000,gx=2000,pn=1000.0,lf=0,rf=0}:\n",aln);
        fprintf(fp,"(%d)",len);
        for(Int4 i=1; i <= len; i++) fprintf(fp,"*");
        fprintf(fp,"\n\n$1=%d(%d):\n",len + gss->NumIns(sq) - gss->NumDel(sq),len);
        StrSeqID(aln, 30, sqE);
        fprintf(fp,">%s ",aln);
        StrSeqDescript(aln,50, sqE);
        fprintf(fp,"%s\n",aln);
        PosSiteCMSA(1, sq, pos, cma);
        Int4 LenStr=gss->Region(sq,aln,pos[1],len);
        fprintf(fp,"{()%s()}*\n\n_0].\n",aln);
	free(aln);
 }
}

cma_typ MinimizeFirstSeqCMSA(cma_typ cma)
{
        Int4    Number,n;
        a_type  A=AlphabetCMSA(cma);
        if(nBlksCMSA(cma) > 1) print_error("MinimizeFirstSeqCMSA( ) input error");

        // create concensus cma
        FILE *fp=tmpfile();
	PutMinimalSeqCMA(fp,1,cma);
	BooLean	*skip; NEW(skip, NumSeqsCMSA(cma) + 3, BooLean); skip[1]=TRUE;
	PutSelectOneCMSA(fp,skip,cma); free(skip);
        rewind(fp);
        cma_typ *IN_CMA=MultiReadCMSA(fp,&Number,A);
        fclose(fp);

        // merge csq.cma and original cma.
        fp=tmpfile(); PutMergedCMSA(fp,Number,IN_CMA); rewind(fp);
        cma_typ CMA=ReadCMSA(fp,A); fclose(fp);
        for(n = 1; n<= Number; n++) TotalNilCMSA(IN_CMA[n]); free(IN_CMA);
        return CMA;
}

void    AddDomainsCMSA(dom_typ *domains, cma_typ cma)
// Add domains structure to cma.
{
	if(cma->Domains) RmDomainsCMSA(cma);
	assert(cma->FullSeq);
	Int4 F = NSeqsSeqSet(cma->FullSeq);
	assert(domains->NumSeq() == F);
	for(Int4 s=1; s <= F; s++){
	   for(Int4 d=1; d <= domains->NumDomains(s); d++){
		assert(domains->SqEnd(s,d) <= SqLenSeqSet(s,cma->FullSeq));
	   }
	}
	cma->Domains = domains;
}

void    CopyDomainsCMSA(cma_typ cma2, cma_typ cma)
{
	if(cma2->Domains) RmDomainsCMSA(cma2);
	cma2->Domains = new dom_typ;
	cma2->Domains[0] = *(cma->Domains);
}

void    RmDomainsCMSA(cma_typ cma)
{ if(cma->Domains) delete cma->Domains; cma->Domains=0; }

void	fill_models_cmsa(cma_typ L)
{
    st_type	S=SitesCMSA(L);
    ss_type	D=DataCMSA(L); 
    mdl_typ	*mdl=mdlCMSA(L);
    for(Int4 t=1; t<= nBlksCMSA(L); t++){
	for(Int4 n=NSeqsSeqSet(D); n > 0; n--){
	   PosTSites(t,n,L->pos,S);
	   for(Int4 k=1; k <= nSites(t,n,S); k++){
		mdl->Add(t,SeqPtrCMSA(n,L), L->pos[k]);
	   }
	}
    }
}

void	delete_models_cmsa(cma_typ cma)
{ delete cma->mdl; cma->mdl=0; }
 
void	create_models_cmsa(cma_typ cma,BooLean **null)
{
	cma->mdl=new mdl_typ(nBlksCMSA(cma),LengthsCMSA(cma),
	 	MaxSeqCMSA(cma),CountsCMSA(cma),
		num_pseudo_msa(NumSeqsCMSA(cma)),
		null,AlphabetCMSA(cma));
}

void	InitMAPCMSA(cma_typ L)
{
     if(L->best){
	if(L->sites != NULL)  NilSites(L->sites);
	L->sites=ExtractSites(L->best); 
	delete_models_cmsa(L); create_models_cmsa(L,L->bestnull);
	fill_models_cmsa(L);
     }
}

void	InitCMSA(cma_typ L) // Initialize (shuffle) alignment.
{
   assert(L->sites != NULL); ShuffleCLSites(L->sites); 
   delete_models_cmsa(L);
   BooLean **null; NEWP(null,nBlksCMSA(L)+3,BooLean);
   create_models_cmsa(L,null); free(null); fill_models_cmsa(L);
}

void    SaveBestCMSA(cma_typ cma)
{
	if(cma->best != NULL) NilArchiveSites(cma->best);
	cma->best = ArchiveSites(cma->sites);
	for(Int4 t = 1; t <= nTypeSites(cma->sites); t++){
	   assert(cma->maxlen[t] >= LengthCMSA(t,cma));
	   NullSitesFModel(cma->bestnull[t],ModelCMSA(t,cma));
	}
	// PutConfigCMSA(stderr,cma); // testing...
}

BooLean	SamePerSeqCMSA(cma_typ L)
/** return TRUE if each sequence has the number of motifs **/
{
	Int4	N,n,t,k,ntyps;

	N = NumSeqsCMSA(L); ntyps=nBlksCMSA(L);
	for(t=1; t<= ntyps; t++) {
                k=nSites(t,1,L->sites); 
                for(n = 2; n <= N; n++){
                   if(k!=nSites(t,n,L->sites)) return FALSE;
                }
	} return TRUE;
}

cma_typ EmptyCMSA(Int4 N, Int4 *lengths, ss_type data,Int4 open,
        Int4 extend, double pernats, Int4 leftflank, Int4 rightflank)
// create a cmsa without initialized sites.
{
  cma_typ cma=make_cmsa(MkMSA(MkSites(N,lengths,data,open,extend,
					pernats,leftflank,rightflank)));
  BooLean **null; NEWP(null,nBlksCMSA(cma)+3,BooLean);
  create_models_cmsa(cma,null); free(null); 
  return cma;
}

cma_typ RandomCMSA(Int4 N, Int4 *lengths, gss_typ& gss)
// create a cmsa that is randomly aligned (i.e., BooLean aligned == FALSE).
{ return make_cmsa(MkMSA(StartSites(N,lengths, gss))); }

cma_typ make_cmsa(cma_typ C)
{ C->lngamma=NULL; calc_lngamma_cmsa(C); return C; }

cma_typ MkMSA(st_type S)
{
        cma_typ L;
        // ss_type data=SitesSeqSet(S);
        ss_type data=SitesTrueSqSet(S);	// longer???
        Int4    t,max,ntyps;

        max= MaxSeqSeqSet(data);
	L = new colinearmaln_type [1]; // replaces NEW(L,1,colinearmaln_type);
        L->sites = S; ntyps=nTypeSites(S);
	L->Level=0;
	L->alpha=0; L->A0=0; L->B0=0; L->set_mode=0;
	L->FirstRpt=0; L->FullRpts=0; L->FullSeq=0; L->best = 0; L->SubToFull=0;
	L->Domains=0; // domain information...
        NEW(L->pos,(2*max)+2,Int4);
        NEW(L->null, (2*max)+2,BooLean);
        NEW(L->maxlen,ntyps +2,Int4);;
        NEWP(L->bestnull, ntyps+2, BooLean);
        for(t = 1; t <= ntyps; t++){
                MEW(L->bestnull[t], (2*max)+2, BooLean);
                L->maxlen[t] = MaxSeqSeqSet(data);
        }
	L->mdl=0;
        return L;
}

BooLean NullSiteCMSA(Int4 t,Int4 s,cma_typ cma)
{ assert(t > 0 && t <= nBlksCMSA(cma)); return NullSiteFModel(s,ModelCMSA(t,cma)); }

BooLean	**NullCMSA(cma_typ ma)
// return the null column configuration for ma
// to free add: for(t=1; t <= nBlksCMSA(ma); t++) free(null[t]); free(null);
{
	BooLean	**null;
	Int4	len,t,nblks = nBlksCMSA(ma);
	fm_type	M;
	
	NEWP(null,nblks+2,BooLean);
        for(t=1; t <= nblks; t++){
	    M = ModelCMSA(t,ma); len = LenFModel(M);
	    NEW(null[t],len+2,BooLean); NullSitesFModel(null[t], M);
	}
	return null;
}

cma_typ	CopyCMSA(cma_typ msa)
// creates and returns an exact copy of msa 
{
	Int4	i,t,n,nblks,*len_elem;
	fm_type	M;
	BooLean	**null;
	cma_typ	ma;

	assert(msa);
	st_type S1=SitesCMSA(msa); nblks=nBlksCMSA(msa); 
	NEWP(null,nblks+2,BooLean); NEW(len_elem,nblks+2,Int4);
	for(t=1; t<= nblks; t++){ 
		M = ModelCMSA(t,msa); len_elem[t] = LenFModel(M);
		NEW(null[t],len_elem[t]+5,BooLean);
	    	NullSitesFModel(null[t], M);
	}
	ma = MakeCMSA(CopySites(S1), null);
        for(t=1; t <= nblks; t++) free(null[t]); free(null);
	free(len_elem); 
	if(msa->best !=NULL){
	    ma->best=CopyArchiveSites(msa->best);
	    for(t=1; t <= nTypeSites(msa->sites); t++){
		for(i=1; i<=ArchivedSiteLen(t,msa->best); i++)
			ma->bestnull[t][i]=msa->bestnull[t][i];
	    }
	}
	if(msa->FullSeq) { CopyFullCountsCMSA(ma, msa); }
	if(msa->Domains) { CopyDomainsCMSA(ma, msa); }
	ma->Level=msa->Level;
	ma->alpha=msa->alpha;
	ma->A0=msa->A0;
	ma->B0=msa->B0;
	ma->set_mode=msa->set_mode;
        for(t=1; t <= nblks; t++) ma->maxlen[t]=msa->maxlen[t];
        return ma;
}

void	TotalNilCMSA(cma_typ cmsa)
{
	ss_type data = TrueDataCMSA(cmsa);
        gss_typ *gssX=gssCMSA(cmsa);    // not owned by cmsa_in!!!
        if(gssX) gssX->~gss_typ();	// if sequences are gapped
        NilCMSA(cmsa);
	NilSeqSet(data);
}

void	NilCMSA(cma_typ L)
{
	for(Int4 i=0; i<=NumSeqsCMSA(L); i++) free(L->lngamma[i]); 
	free(L->lngamma);
        for(Int4 t=1; t <= nBlksCMSA(L); t++){
           if(L->bestnull[t]!=NULL) free(L->bestnull[t]);
        }
        delete_models_cmsa(L);
        free(L->pos); free(L->maxlen); free(L->bestnull);
        free(L->null); 
        if(L->best != NULL) NilArchiveSites(L->best);
        if(L->sites != NULL) NilSites(L->sites); 
	if(L->FullSeq) RmFullCountsCMSA(L);
	if(L->Domains) RmDomainsCMSA(L);
	delete []L;
}

void	PutConfigCMSA(FILE *fp, cma_typ cma)
{
	Int4	min,n,i,b,B=nBlksCMSA(cma);
	UInt8   total;
	for(min=MaxSeqCMSA(cma),n=NumSeqsCMSA(cma); n > 0; n--){
	   i=GapBetweenSites(n,0,SitesCMSA(cma));
	   min = MINIMUM(Int4,i,min);
	}
	fprintf(fp,"%2d blocks [ (%d)",B,min);
	for(b=1; b<=B; b++){
	   for(total=0,min=MaxSeqCMSA(cma),n=NumSeqsCMSA(cma); n > 0; n--){
	   	i=GapBetweenSites(n,b,SitesCMSA(cma));
		total+=i;
		min = MINIMUM(Int4,i,min);
	   }
	   // fprintf(fp," %d (%d)",LengthCMSA(b,cma),min);
	   fprintf(fp," %d (%.0f)",LengthCMSA(b,cma),(double)total/(double)NumSeqsCMSA(cma));
	} fprintf(fp,"]\n");
}

double	GetBPPS_CMA(Int4 *A0, Int4 *B0, char *set_mode, cma_typ cma)
{ *A0=cma->A0; *B0=cma->B0; *set_mode=cma->set_mode; return cma->alpha; }

void    SetBPPS_CMA(double alpha,Int4 A0, Int4 B0,char set_mode, cma_typ cma)
{ cma->alpha=alpha; cma->A0=A0; cma->set_mode = set_mode, cma->B0=B0; }

double ExpectedGapLengthCMSA(cma_typ cma, Int4 blk)
//returns expected insert length between blk and blk+1 including ends
// From Aleksandar's code
{
        double          exp_length=0;
        unsigned short  *BlkLen;
        Int4            i, nBl=nBlksCMSA(cma);

        assert((blk >=0) && (blk <= nBl));

        Int4            n_seq = NumSeqsCMSA(cma);
        Int4            pos[5],posp1[5];
        gss_typ         *gss = gssCMSA(cma);
        e_type          fakeE;

        NEW(BlkLen,nBl+2,unsigned short);
        for(i=1;i<=nBl;i++){ BlkLen[i]=LengthCMSA(i,cma); }
        if(blk>0 && blk<nBl){
                for(i=1;i<=n_seq;i++){
                        PosSiteCMSA(blk+1,i,posp1,cma);
                        PosSiteCMSA(blk,i,pos,cma);
                        exp_length += posp1[1] - (pos[1] + BlkLen[blk]);
                }
        } else if(blk == 0){
                for(i=1;i<=n_seq;i++){
                        PosSiteCMSA(1,i,posp1,cma);
                        exp_length += posp1[1] - 1;
                }
        } else{
                for(i=1;i<=n_seq;i++){
                        fakeE = gss->FakeSeq(i);
                        PosSiteCMSA(nBl,i,pos,cma);
                        exp_length += LenSeq(fakeE) - (pos[1] + BlkLen[nBl]) + 1;
                }
        } exp_length /= n_seq;
        free(BlkLen);
        return exp_length;
}

Int4	ConfigCMSA2(Int4 *width, cma_typ L)
/** returns # blocks and sets width == widths of blocks **/
{
    Int4	t,k;
    st_type	S=SitesCMSA(L);
    fm_type	M;

    k = nTypeSites(S);
    for(t=1; t <= k; t++){ M=ModelCMSA(t,L); width[t]=LenFModel(M); }
    return k;
}

Int4	ConfigCMSA(Int4 *ncols, cma_typ L)
/** returns # blocks and sets ncols == number of columns **/
{
    Int4	t,k;
    st_type	S=SitesCMSA(L);
    fm_type	M;

    k = nTypeSites(S);
    for(t=1; t <= k; t++){ M=ModelCMSA(t,L); ncols[t]=nColsFModel(M); }
    return k;
}

void	UnLabelSeqsCMSA(cma_typ cma)
{
	for(Int4 i=1; i <= NumSeqsCMSA(cma); i++){
		e_type tE=TrueSeqCMSA(i,cma),fE=FakeSeqCMSA(i,cma);
		UnLabelSeq(tE); UnLabelSeq(fE);
	}
}

void	LabelSeqsCMSA(cma_typ cma)
{
	for(Int4 i=1; i <= NumSeqsCMSA(cma); i++){
		e_type tE=TrueSeqCMSA(i,cma),fE=FakeSeqCMSA(i,cma);
		LabelSeq(tE); LabelSeq(fE);
	}
}

unsigned short InsertionCMSA(UInt4 blk, UInt4 n, UInt4 i, cma_typ cma)
{
	Int4	pos[4];
	assert(blk > 0 && blk <= nBlksCMSA(cma));
	assert(n > 0 && n <= NumSeqsCMSA(cma));
	assert(i > 0 && i <= LengthCMSA(blk,cma));
	assert(PosSiteCMSA(blk,n,pos,cma));
	gss_typ	*gss=gssCMSA(cma);
	if(gss->Gapped()) return (gss->Insertion(n,pos[1]+i-1)); else return 0; 
}

BooLean IsDeletedCMSA(UInt4 blk, UInt4 n, UInt4 r, cma_typ cma)
// returns true if residue r in sequence n is deleted.
{
	Int4	pos[4];
	assert(blk > 0 && blk <= nBlksCMSA(cma));
	assert(n > 0 && n <= NumSeqsCMSA(cma));
	assert(r > 0 && r <= LengthCMSA(blk,cma));
	assert(PosSiteCMSA(blk,n,pos,cma));
	gss_typ	*gss=gssCMSA(cma);
	if(gss->Gapped()) return (gss->IsDeleted(n,pos[1]+r-1)); else return FALSE; 
}


unsigned short InsertionCMSA(UInt4 n, UInt4 i, cma_typ cma)
// WARNING: THIS FUNCTION HAS PROBLEMS; IT FAILS TO REPORT DELETIONS AT ENDS;  NEED TO FIX.
// IT APPEARS THAT i NEEDS TO BE THE POSITION IN THE FAKE SEQ RATHER THEN TRUE SEQ; NEED TO CHECK!!!!
// returns the number of residues inserted immediately following 
// position i in sequence n.
{
	gss_typ	*gss=gssCMSA(cma);
	if(gss->Gapped()) return (gss->Insertion(n,i)); else return 0; 
}

BooLean IsDeletedCMSA(UInt4 n, UInt4 r, cma_typ cma)
// returns true if residue r in sequence n is deleted.
{
	gss_typ	*gss=gssCMSA(cma);
	if(gss->Gapped()) return (gss->IsDeleted(n,r)); else return FALSE; 
}

#if 0
double	GetGappedSqLPR_CMSA(Int4 s, cma_typ cma)
// get the probability ratio of sequence s being in the alignment 
// WARNING: assumes gap penalties are meaningful!!!
{
	Int4	*oldpos;
	double	oldmap,newmap;

	oldmap=RelMapCMSA(cma);
	// 1. Remove sequence s from alignment.
	NEW(oldpos,nBlksCMSA(cma)+2,Int4);  // save old sites.
	for(Int4 t=1; t<=nBlksCMSA(cma); t++){
		PosSiteCMSA(t, s, cma->pos, cma); oldpos[t]=cma->pos[1];
	} VacateSitesCMSA(s,cma);
	newmap=RelMapCMSA(cma);
	for(t=1 ; t <= nBlksCMSA(cma); t++) AddSiteCMSA(t,s,oldpos[t], cma);
	free(oldpos); 
	return ((oldmap-newmap)*0.43429448); // convert to log10
}
#endif

double	GetGappedProbCMSA(UInt4 t,UInt4 n, cma_typ cma)
{
	gss_typ     *gss=gssCMSA(cma);
	double	prob=GetProbCMSA(t,n,cma);

	if(gss->Gapped()){
		Int4	inso,insx,del,e,pos[3];
		double	p;
		PosSiteCMSA(t,n,pos,cma);
		e = pos[1] + LengthCMSA(t,cma) - 1;
		del = gss->InDels(n,pos[1],e,&inso,&insx);
// fprintf(stderr,"#### del = %d; inso = %d; insx = %d\n", del,inso,insx);
		if(insx || del) p = gss->IndelPenalty(del,inso,insx); else p = 0;
// fprintf(stderr,"#### prob = %g; p = %f; pos[1] = %d\n", prob,p,pos[1]);
		return (prob - p);
	} return prob;
}

double  GetTotalGappedProbCMSA(Int4 n, cma_typ cma)
{
	assert(n > 0 && n <= NumSeqsCMSA(cma));
	double	prob=0.0;
	for(Int4 b=1; b<=nBlksCMSA(cma); b++){
		prob += GetGappedProbCMSA(b,n,cma);
	} return prob;
}

double  GetTotalProbCMSA(Int4 n, cma_typ cma)
{
	assert(n > 0 && n <= NumSeqsCMSA(cma));
        assert(SeqIsInCMSA(n,cma));
	double	prob=0.0;
	for(Int4 b=1; b<=nBlksCMSA(cma); b++){
		prob += GetProbCMSA(b,n,cma);
	} return prob;
}

e_type	GetBestConSqCMSA(cma_typ mcma)
{
    a_type  AB=AlphabetCMSA(mcma);
    FILE    *fp=tmpfile(); PutBestCMSA(fp,1,FALSE,mcma); rewind(fp);
    cma_typ cma=ReadCMSA(fp,AB); fclose(fp);
    e_type  rtn=MkConsensusCMSA(cma); TotalNilCMSA(cma);
    return rtn;
}

cma_typ	GetBestCsqAsCMSA(set_typ set, cma_typ cma)
// return the ith sequence in the alignment as a consensus (no indels).
{
        Int4    N = NumSeqsCMSA(cma);
        a_type  AB=AlphabetCMSA(cma);
        BooLean *skip; NEW(skip,N+3,BooLean);
        for(Int4 j=1; j <= N; j++){ if(!MemberSet(j,set)) skip[j]=TRUE; }
        FILE    *fp=tmpfile(); PutSelectCMSA(fp,skip,cma); free(skip);
        rewind(fp); cma_typ tcma=ReadCMSA(fp,AB); fclose(fp);
        cma_typ rcma=GetBestCsqCMSA(tcma); TotalNilCMSA(tcma);
        return rcma;
}


cma_typ RtnConSqAsCMSA(char *name,sst_typ *xsst,e_type  keyE, a_type AB)
// return the consensus sequence corresponding to the FG for column n.
// if(xsst != 0) then redefine incompatible sets at pattern positions to match keyE.
{
        unsigned char r;
        FILE *fp=tmpfile();
        fprintf(fp,"[0_(1)=%s(1){go=10000,gx=2000,pn=1000.0,lf=0,rf=0}:\n",name);
        fprintf(fp,"(%d)",LenSeq(keyE));
        for(Int4 i=1; i <= LenSeq(keyE); i++) fprintf(fp,"*");
        fprintf(fp,"\n\n$1=%d(%d):\n",LenSeq(keyE),LenSeq(keyE));
        fprintf(fp,">%s consensus\n{()",name);
        for(Int4 i=1; i <= LenSeq(keyE); i++){
            r=ResSeq(i,keyE);
            if(xsst && xsst[i]){ if(!MemSset(r,xsst[i])){ xsst[i] = SsetLet(r); } }
            fprintf(fp,"%c",AlphaChar(r,AB));
        } fprintf(fp,"()}*\n\n_0].\n"); rewind(fp);
        cma_typ cma=ReadCMSA(fp,AB); fclose(fp); 
        return cma;
}

e_type  GetSeqAsCsqCMSA(Int4 i, cma_typ cma)
// return the ith sequence in the alignment as a consensus (no indels).
{
        Int4    N = NumSeqsCMSA(cma);
        a_type  AB=AlphabetCMSA(cma);
        assert(i > 0 & i <= N);
        BooLean *skip; NEW(skip,N+3,BooLean);
        for(Int4 j=1; j <= N; j++){ skip[j]=TRUE; } skip[i]=FALSE;
        FILE    *fp=tmpfile(); PutSelectCMSA(fp,skip,cma); free(skip);
        rewind(fp); cma_typ tcma=ReadCMSA(fp,AB); fclose(fp);
        e_type  rtnE=MkConsensusCMSA(tcma); TotalNilCMSA(tcma);
        return rtnE;
}

e_type  GetSeqAsCsqCMSA(set_typ set, cma_typ cma)
// return the ith sequence in the alignment as a consensus (no indels).
{
        Int4    N = NumSeqsCMSA(cma);
        a_type  AB=AlphabetCMSA(cma);
        BooLean *skip; NEW(skip,N+3,BooLean);
        for(Int4 j=1; j <= N; j++){ if(!MemberSet(j,set)) skip[j]=TRUE; }
        FILE    *fp=tmpfile(); PutSelectCMSA(fp,skip,cma); free(skip);
        rewind(fp); cma_typ tcma=ReadCMSA(fp,AB); fclose(fp);
        e_type  rtnE=MkConsensusCMSA(tcma); TotalNilCMSA(tcma);
        return rtnE;
}

cma_typ GetBestCsqCMSA(cma_typ mcma)
{
    a_type  AB=AlphabetCMSA(mcma);
    FILE    *fp=tmpfile(); PutBestCMSA(fp,1,FALSE,mcma); rewind(fp);
    cma_typ cma=ReadCMSA(fp,AB); fclose(fp);
    cma_typ rtn=MakeConsensusCMSA(cma); TotalNilCMSA(cma);
    return rtn;
}

cma_typ GetInSetCMSA(set_typ set, cma_typ mcma)
{
    a_type  AB=AlphabetCMSA(mcma);
    FILE    *fp=tmpfile(); PutInSetCMSA(fp,set,mcma); rewind(fp);
    cma_typ cma=ReadCMSA(fp,AB); fclose(fp);
    return cma;
}

double	GetProbCMSA(Int4 t, Int4 n, cma_typ L)
{
    ss_type		data=DataCMSA(L);
    st_type		sites = SitesCMSA(L);
    fm_type		*model=ModelsCMSA(L);
    Int4		s,end;
    double		prob;
    unsigned char	*seq;
    mdl_typ		*mdl=mdlCMSA(L);

    if(n < 1 || n > NSeqsSeqSet(data) || t < 1 || t > nTypeSites(sites))
			print_error("GetProbCMSA( ) input error");
    end = SqLenSeqSet(n,data) + 1;
    end += 1 - LenFModel(model[t]);
    seq = SeqSeqSet(n,data);
    s = SitePos(t,n,1,sites);
    mdl->Remove(t,seq,s);
    prob = (double) LikelihoodFModel(seq, s, model[t]);
    mdl->Add(t,seq,s);
    if(prob > 0.0) prob=log10(prob);
    else {
	fprintf(stderr,"Seq %d: prob = %g!\n",n,prob);
	// PutSeq(stderr,TrueSeqCMSA(n,L),AlphabetCMSA(L));
	prob=(double) INT4_MIN;
	// assert(prob > 0.0);
    }
    return prob;
}

Int4	FastDangerousResCMSA(register Int4 t, register Int4 n, register Int4 s, 
	register cma_typ cma)
// assumes input is okay...may be dangerous.
{
    	register unsigned char *seq = SeqSeqSet(n,DataCMSA(cma));
	return seq[SitePos(t,n,1,SitesCMSA(cma))+s-1];
}

Int4	SetResidueCMSA(Int4 t, Int4 n, Int4 s, unsigned char r, cma_typ cma)
// return the residue at position s in sequence n of block t
{
	if(t < 1 || t > nBlksCMSA(cma)) return 0;
	if(n < 1 || n > NumSeqsCMSA(cma) || s < 1 || s > LengthCMSA(t,cma)) return 0;
    	unsigned char *seq = SeqSeqSet(n,DataCMSA(cma)); // fake seq...
	Int4 site=SitePos(t,n,1,SitesCMSA(cma))+s-1; seq[site]=r;
	return site;
}

Int4	ResidueCMSA(register Int4 t, register Int4 n, register Int4 s, cma_typ cma)
// return the residue at position s in sequence n of block t
{
	if(t < 1 || t > nBlksCMSA(cma)) return 0;
	if(n < 1 || n > NumSeqsCMSA(cma) || s < 1 || s > LengthCMSA(t,cma)) return 0;
    	register unsigned char *seq = SeqSeqSet(n,DataCMSA(cma));
	return seq[SitePos(t,n,1,SitesCMSA(cma))+s-1];
}

Int4	ResBetweenCMSA(register Int4 t, register Int4 n, register cma_typ cma)
// return the number of residues between block t and t+1.
{
	if(t < 1 || t >= nBlksCMSA(cma)) return 0;
	assert(n > 0 && n <= NumSeqsCMSA(cma));
	return GapBetweenSites(n,t,SitesCMSA(cma));
}

Int4	ResBetweenCMSA(register Int4 t, register cma_typ cma)
// return the number of residues between block t and t+1.
{
	register Int4	n,res;
	if(t < 1 || t >= nBlksCMSA(cma)) return 0;
	else for(res=0,n=1; n <= NumSeqsCMSA(cma); n++){
		res += GapBetweenSites(n,t,SitesCMSA(cma));
	}
	return res;
}

Int4	TotalLenCMSA(cma_typ msa)
{
    Int4	t,n;
    for(n=0,t=1; t <= nBlksCMSA(msa); t++) n+=LengthCMSA(t,msa);
    return n;
}

Int4	NumColumnsCMSA(cma_typ msa)
{
    Int4	t,k,n;
    k = nBlksCMSA(msa);
    for(n=0, t=1; t <= k; t++) n += nColsFModel(ModelCMSA(t,msa));
    return n;
}

Int4    RealToFakeCMSA(Int4 sq, Int4 s, cma_typ cma)
// return the fake seq position for residue s in sq (with offset).
// returns 0 if no corresponding position in fake seq.
{
	Int4	pos[5];
	assert(nBlksCMSA(cma)==1);
	assert(sq > 0 && sq <= NumSeqsCMSA(cma));
        gss_typ *gss=gssCMSA(cma);
        assert(gss != NULL);
        return gss->RealToFake(sq,s);
}

Int4	FakeToRealCMA(Int4 sq,Int4 s, cma_typ cma)
// return the real seq position for residue s in sq (with offset)
// returns 0 if no corresponding position in real seq.
// NOTE: this is different from TruePosCMSA(), which returns the position w/o offset.
{
	Int4	pos[5];
	assert(nBlksCMSA(cma)==1);
	assert(sq > 0 && sq <= NumSeqsCMSA(cma));
        gss_typ *gss=gssCMSA(cma);
        assert(gss != NULL);
        return gss->TrueSite(sq,s);
}

Int4    TruePosCMSA(Int4 sq, Int4 s, cma_typ cma)
// return the true position (possibly in subsequence) for sq at column s in cma.
{
	Int4	pos[5];
	assert(nBlksCMSA(cma)==1);
        assert(s > 0 && s <= LengthCMSA(1,cma));
        gss_typ *gss=gssCMSA(cma);
        assert(gss != NULL);
	assert(sq > 0 && sq <= NumSeqsCMSA(cma));
	// assert(PosSiteCMSA(1,sq,pos,cma));
	PosSiteCMSA(1,sq,pos,cma);
        // return gss->TruePos(sq,pos[1]); fixed ...
        return gss->TruePos(sq,pos[1]+s-1);
}

Int4    TruePosCMSA(Int4 sq, Int4 blk, Int4 s, cma_typ cma)
// return the true position (possibly in subsequence) for sq at column s in block='blk' of cma.
{
        Int4    pos[5];
        assert(blk > 0 && blk <= nBlksCMSA(cma));
        assert(s > 0 && s <= LengthCMSA(blk,cma));
        gss_typ *gss=gssCMSA(cma);
        assert(gss != NULL);
        assert(sq > 0 && sq <= NumSeqsCMSA(cma));
        // assert(PosSiteCMSA(1,sq,pos,cma));
        PosSiteCMSA(blk,sq,pos,cma);
        // return gss->TruePos(sq,pos[1]); fixed ...
        return gss->TruePos(sq,pos[1]+s-1);
}

Int4	PosSiteCMSA(Int4 t, Int4 n, Int4 *pos, cma_typ cmsa)
// return the position of the t'th sites in sequence n
// returns the number of sites found.
{ return PosTSites(t,n,pos,SitesCMSA(cmsa)); }

void    PutAlnSeqsCMSA(FILE *fp,cma_typ cma)
{
        Int4    sq,NN=NumSeqsCMSA(cma);
        set_typ Set=MakeSet(NN+5);  ClearSet(Set);
        for(sq=1; sq <= NN; sq++) if(SeqIsInCMSA(sq,cma)) AddSet(sq,Set);
        PutInSetCMSA(fp,Set,cma);  NilSet(Set);
}

BooLean SeqIsInCMSA(Int4 sq,cma_typ cmsa)
{ Int4 pos[5]; if(PosTSites(1,sq,pos,SitesCMSA(cmsa))) return TRUE; else return FALSE; }

Int4	*CntsFieldCMSA(cma_typ L, Int4 t, Int4 *len)
/**********************************************************************
  Return the residue counts for the region surrounding element t.
  If element t does not exist then return NULL.
/**********************************************************************/
{
   Int4		n,N,i,r,start,end,ntyp;
   ss_type	data = DataCMSA(L);
   st_type	S=SitesCMSA(L);
   Int4	*counts;
   e_type	E;
   a_type	A;

   N = NSeqsSeqSet(data); A = SeqSetA(data);
   ntyp = nTypeSites(S);
   NEW(counts,nAlpha(A)+2,Int4);
   if(ntyp==1){
	for(r=0; r<=nAlpha(A); r++) counts[r]=CountsSeqSet(r,data); 
	for(n=1; n <= N; n++) len[n]=LenSeq(SeqSetE(n,data));
   } else {
	for(n=1; n <= N; n++){
	  E = SeqSetE(n,data);
	  if(nSites(t,n,S)!=1) print_error("CntsFieldCMSA( ) error");
	  if(t==1){
	    if(nSites(t+1,n,S)!=1) print_error("CntsFieldCMSA( ) error");
	    start=1; end=SitePos(t+1,n,1,S);
	  } else if(t==ntyp){
	    if(nSites(t-1,n,S)!=1) print_error("CntsFieldCMSA( ) error");
	    start = SitePos(t-1,n,1,S) + SiteLen(t-1,S);
	    end = LenSeq(E);
	  } else {
	    if(nSites(t+1,n,S)!=1) print_error("CntsFieldCMSA( ) error");
	    if(nSites(t-1,n,S)!=1) print_error("CntsFieldCMSA( ) error");
	    start = SitePos(t-1,n,1,S) + SiteLen(t-1,S);
	    end=SitePos(t+1,n,1,S);

	  }
	  len[n] = end - start + 1;
	  for(i=start; i < end; i++){ r=ResSeq(i,E); counts[r]++; }
	}
   }
   return counts;
}

e_type	MkConsensusCMSA(cma_typ cma)
//  Return an e_type consensus sequence for cma;
{
	a_type	A = AlphabetCMSA(cma);
	char	*buffer = ConsensusSeqCMSA(cma);
	e_type E=StringToSeq(buffer+1, "consensus seq", 1, A);
	free(buffer);
	return E;
}

e_type	MkConsensusCMSA(set_typ Set, cma_typ cma)
//  Return an e_type consensus sequence for cma;
{
	a_type	A = AlphabetCMSA(cma);
	assert(nBlksCMSA(cma) == 1);
	char	*buffer = ConsensusSeqCMSA(1,Set,cma);
	e_type E=StringToSeq(buffer, "consensus seq", 1, A);
	free(buffer); return E;
}

char	*ConsensusSeqCMSA(cma_typ cma)
//  Return a consensus sequence for entire sequence;
{
	ss_type data=DataCMSA(cma);
	st_type	S=SitesCMSA(cma);
	Int4	best,s,r,j,*cnts,d,lenM;
	char	*rtnseq;
	a_type	A = AlphabetCMSA(cma);

   NEW(rtnseq,TotalLenCMSA(cma)+3,char); s=1;
   for(Int4 t=1; t <= nBlksCMSA(cma); t++){
    lenM=LengthCMSA(t,cma);
    // 2. generate consensus sequence.
    Int4 N=NSeqsSeqSet(data);
    for(d=0; d < lenM; d++){
	cnts = GetSiteFreq(S,t,d);
	for(r=best=0,j=1; j <= nAlpha(A); j++){
	   if(cnts[j] > best) { best = cnts[j]; r = j; }
	} free(cnts);
	//fprintf(stderr,"best = %d; r = %c\n", best, AlphaChar(r,A));
	rtnseq[s] = AlphaChar(r,A); s++;
    }
   } return rtnseq;
}

char	*ConsensusSeqCMSA(Int4 t,cma_typ cma){ return ConsensusSeqCMSA(t,0,cma); }

char	*ConsensusSeqCMSA(Int4 t, set_typ Set, cma_typ cma)
/*************************************************************************
 Return a consensus sequence for block t;
 WARNING: assumes that "InitMAPMSA(cma); " has been invoked!!!
 WARNING: assumes alignment is COLINEAR!!
 *************************************************************************/
// AlwaysGetSiteFreq(st_type S,Int4 t,BooLean right,Int4 *pos,char *blocked)
// case 1: -inf...-1 is left of site
// case 2:  0...+inf is within or beyond site
// if right == TRUE then ...
{
	ss_type data=DataCMSA(cma);
	st_type	S=SitesCMSA(cma);
	double  *freq = tFreqSeqSet(data),likelihood,best;
	Int4	s,f,r,i,j,*cnts,d,lenM,lenSeq;
	char	*rtnseq;
	a_type	A = AlphabetCMSA(cma);

    assert(t > 0 && t <= nBlksCMSA(cma));
    lenM=LengthCMSA(t,cma);
    NEW(rtnseq,lenM+3,char);
    // 2. generate consensus sequence.
    Int4 N,N0; N=N0=NSeqsSeqSet(data);
    if(Set){ for(Int4 n=1; n <= N0; n++) if(!MemberSet(n,Set)) N--; }
    for(s=d=0; d < lenM; d++){
	if(Set) cnts = GetSiteFreq(S,Set,t,d);
	else cnts = GetSiteFreq(S,t,d);
	for(r=0,best=0.0,j=1; j <= nAlpha(A); j++){
	   if(freq[j] > 0.0) likelihood = ((double) cnts[j]/(double) N)/freq[j];
	   else likelihood = 0.0;
           if(best < likelihood){ best=likelihood; r=j; }
	} free(cnts);
	//fprintf(stderr,"best = %g; r = %c\n", best, AlphaChar(r,A));
	rtnseq[s] = AlphaChar(r,A); s++;
    } rtnseq[s] = 0;
    if(lenM != s) print_error("Consensus length error");
    return rtnseq;
}

double	**ColResFreqsCMSA(Int4 t, cma_typ cma)
{ BooLean *skip=0; return ColResFreqsCMSA(t,skip,cma); }

double	**ColResFreqsCMSA(Int4 t, BooLean *skip, cma_typ cma)
{
	Int4	**observed;
	double	**rtnfreq=ColResFreqsCMSA(t,skip, &observed, cma);
	for(Int4 s=1; s <= LengthCMSA(t,cma); s++) free(observed[s]); 
	free(observed);
	return rtnfreq;
}

double	**ColResFreqsCMSA(Int4 t, double ***observed, cma_typ cma)
{ return ColResFreqsCMSA(t,0,observed,cma); }

double	**ColResFreqsCMSA(Int4 t, BooLean *skip, double ***observed, cma_typ cma)
{
	Int4	**obs; 
	double **result= ColResFreqsCMSA(t,skip,&obs,cma); 
	double 	**Obs;
	a_type  A = AlphabetCMSA(cma);
	NEWP(Obs,LengthCMSA(t,cma)+3,double);
	for(Int4 d=1; d <= LengthCMSA(t,cma); d++){
	   NEW(Obs[d],nAlpha(A)+3,double);
	   for(Int4 r=0; r <= nAlpha(A); r++){
		 Obs[d][r]=(double) obs[d][r];
	   } free(obs[d]);
	} *observed=Obs; free(obs);
	return result;
}


double	**ColResFreqsCMSA(Int4 t, Int4 ***observed, cma_typ cma)
{ return ColResFreqsCMSA(t,(BooLean *) 0,observed,cma); }

double	**ColResFreqsCMSA(Int4 t, set_typ Set, Int4 ***observed, cma_typ cma)
/*************************************************************************
 Return a consensus sequence for block t;
 WARNING: assumes that "InitMAPMSA(cma); " has been invoked!!!
 WARNING: assumes alignment is COLINEAR!!
 *************************************************************************/
// AlwaysGetSiteFreq(st_type S,Int4 t,BooLean right,Int4 *pos,char *blocked)
// case 1: -inf...-1 is left of site
// case 2:  0...+inf is within or beyond site
// if right == TRUE then ...
{
	ss_type data=DataCMSA(cma);
	st_type	S=SitesCMSA(cma);
	double  *freq = tFreqSeqSet(data),*tmpfreq;
	Int4	s,j,*cnts,d,lenM,**obs;
	double	**rtnfreq;
	a_type	A = AlphabetCMSA(cma);

    assert(t > 0 && t <= nBlksCMSA(cma));
    lenM=LengthCMSA(t,cma);
    NEWP(rtnfreq,lenM+3,double);
    NEWP(obs,lenM+3,Int4);
    // 2. generate consensus sequence.
    Int4 N,N0; N=N0=NSeqsSeqSet(data);
    if(Set){ for(Int4 n=1; n <= N0; n++) if(!MemberSet(n,Set)) N--; }
    for(s=d=0; d < lenM; d++){
	if(Set) cnts = GetSiteFreq(S,Set,t,d);
	else cnts = GetSiteFreq(S,t,d);
	NEW(tmpfreq, nAlpha(A)+3, double);
	for(j=0; j <= nAlpha(A); j++) tmpfreq[j] = ((double)cnts[j]/(double) N);
	s++; rtnfreq[s] = tmpfreq;  obs[s]=cnts;
    } if(lenM != s) print_error("Consensus length error");
    *observed=obs;
    return rtnfreq;
}

double	**ColResFreqsCMSA(Int4 t, BooLean *skip, Int4 ***observed, cma_typ cma)
/*************************************************************************
 Return a consensus sequence for block t;
 WARNING: assumes that "InitMAPMSA(cma); " has been invoked!!!
 WARNING: assumes alignment is COLINEAR!!
 *************************************************************************/
// AlwaysGetSiteFreq(st_type S,Int4 t,BooLean right,Int4 *pos,char *blocked)
// case 1: -inf...-1 is left of site
// case 2:  0...+inf is within or beyond site
// if right == TRUE then ...
{
	ss_type data=DataCMSA(cma);
	st_type	S=SitesCMSA(cma);
	double  *freq = tFreqSeqSet(data),*tmpfreq;
	Int4	s,j,*cnts,d,lenM,**obs;
	double	**rtnfreq;
	a_type	A = AlphabetCMSA(cma);

    assert(t > 0 && t <= nBlksCMSA(cma));
    lenM=LengthCMSA(t,cma);
    NEWP(rtnfreq,lenM+3,double);
    NEWP(obs,lenM+3,Int4);
    // 2. generate consensus sequence.
    Int4 N,N0; N=N0=NSeqsSeqSet(data);
    if(skip){ for(Int4 n=1; n <= N0; n++) if(skip[n]) N--; }
    for(s=d=0; d < lenM; d++){
	if(skip) cnts = GetSiteFreq(S,skip,t,d);
	else cnts = GetSiteFreq(S,t,d);
	NEW(tmpfreq, nAlpha(A)+3, double);
	for(j=0; j <= nAlpha(A); j++) tmpfreq[j] = ((double)cnts[j]/(double) N);
	s++; rtnfreq[s] = tmpfreq;  obs[s]=cnts;
    } if(lenM != s) print_error("Consensus length error");
    *observed=obs;
    return rtnfreq;
}

float	*RelEntropyCMA(Int4 blk, cma_typ cma)
{
	fm_type fm = ModelCMSA(blk,cma);
	return InfoFModel(fm);
}

Int4	RepeatsInfoCMSA(Int4 *start,Int4 *end,Int4 *gap0,Int4 *gap1,Int4 sq, cma_typ cma)
{
	Int4	s,e,e0,lenM,N,x;
	ss_type	data = DataCMSA(cma);
	st_type	sites=SitesCMSA(cma);
	gss_typ	*gss=gssCMSA(cma);
	char	str[108],str2[108];
	e_type	E;

	N = NSeqsSeqSet(data);

	assert(sq > 0 && sq <= N);
	assert(N > 1); assert(gss);
	assert(nBlksCMSA(cma) == 1);

	E = SeqSetE(sq,data); lenM = LengthCMSA(1,cma); StrSeqID(str,100,E);
	Int4 qst = SitePos(1,sq,1,sites);
	s = gss->TrueSite(sq, qst); *start = s;
	*gap0 = s-1;  // default for first gap.
	e = gss->TrueSite(sq,qst+lenM-1); *end = e;
	x = gss->TrueSite(sq,LenSeq(E))+CtermExtendSeq(E);
	*gap1 = x-e;  // default for second gap.
	if(sq > 1){  // i.e., room for previous repeat..
	   StrSeqID(str2,100,SeqSetE(sq-1,data));
	   if(strncmp(str2,str,100) == 0){
		x = SitePos(1,sq-1,1,sites)+lenM-1;
		x = gss->TrueSite(sq-1,x);
		*gap0 = s-x-1;
	   } // otherwise use default setting.
	} 
	if(sq < N){	// room for next repeat
	   StrSeqID(str2,100,SeqSetE(sq+1,data));
	   if(strncmp(str2,str,100) == 0){  // i.e., the same sequence
		x = gss->TrueSite(sq+1,SitePos(1,sq+1,1,sites));
		*gap1 = x-e-1;
	   } // otherwise use default setting.
	} return qst;
}

e_type	MaxWordAndSubSeqCMSA(Int4 t, Int4 qsq, Int4 ssq, Int4 *Word_score,
		Int4 *Q_start, Int4 *S_start, Int4 *Ssq_start, cma_typ cma)
{
	Int4	os,N;
	ss_type	data = DataCMSA(cma);
	st_type	sites=SitesCMSA(cma);
	gss_typ	*gss=gssCMSA(cma);
	ss_type	fulldata = FullSeqCMSA(cma);

	N = NSeqsSeqSet(data);

	assert(N > 1); assert(gss); assert(fulldata);
	assert(qsq > 0 && qsq <= N);
	assert(ssq > 0 && ssq <= N);
	assert(nBlksCMSA(cma) == 1);  // Not needed???

	Int4 qst = SitePos(t,qsq,1,sites);
	Int4 sst = SitePos(t,ssq,1,sites);

	e_type qE = SeqSetE(qsq,data);
	e_type sE = SeqSetE(ssq,data);
	Int4 lenM = LengthCMSA(t,cma);
	os=FindMaxWordSeq(qst,sst,lenM, qE, sE, Word_score, AlphabetCMSA(cma));
	Int4 q_start = gss->TrueSite(qsq, qst+os);
	Int4 s_start = gss->TrueSite(ssq, sst+os);

	e_type fsE = SeqSetE(SubToFullCMSA(ssq,cma),fulldata);

	Int4 ssq_start = gss->TrueSite(ssq,sst)-lenM;
	ssq_start = MAXIMUM(Int4,ssq_start,1);

	Int4 ssq_end = gss->TrueSite(ssq,sst+lenM-1)+lenM;
	ssq_end = MINIMUM(Int4,ssq_end,LenSeq(fsE));

	*Q_start=q_start; *S_start=s_start; *Ssq_start=ssq_start;
	return MkSubSeq(ssq_start, ssq_end, fsE);
}

smx_typ SMXforSeqCMSA(Int4 t,Int4 sq,cma_typ cma)
{
	Int4	s,j,c;
	ss_type	data = DataCMSA(cma);
	a_type	A=AlphabetCMSA(cma);
	Int4	lenM = LengthCMSA(t,cma);
	smx_typ smx = MkSMatrix(2.0,lenM,tFreqSeqSet(data),A);
	e_type	E = SeqSetE(sq,data);
	unsigned char *psq = SeqPtr(E);
	Int4	st = SitePos(t,sq,1,SitesCMSA(cma));
	for(s=st,j=1; j <= lenM; j++,s++) {
	  for(c=0; c <= nAlpha(A); c++){
		SetSMatrix(c,j,valAlphaR(psq[s],c,A),smx);
	  }
	}
	return smx;
}

BooLean	*FindCloseCMA(double cutoff, e_type E,Int4 *start, cma_typ cma)
// return list of sequences closely related to sequence E starting at start
// cutoff == percent identity
{
	Int4	totlen,s,t,N,score,s1,s2,item1,hits,num_seqs=0;
	BooLean	*itemList;
	a_type	A=AlphabetCMSA(cma);
	unsigned char	*seq1,*seq2;
	double	prcntID;
	ss_type	data=DataCMSA(cma);
	Int4	exists,pos[4];

	// assert(nBlksCMSA(cma) == 1);
	assert(cutoff > 0.0);
	N=NumSeqsCMSA(cma);
	NEW(itemList,N+3,BooLean);
	if(cutoff > 100.) return itemList;
	for(totlen=0,t=1; t <= nBlksCMSA(cma); t++) totlen+=LengthCMSA(t,cma);
	seq1 = SeqPtr(E);
	for(hits=0,item1 =1 ;item1 <= N;item1++) {
	   seq2 =  SeqPtrCMSA(item1,cma);
	   for(score=0,t=1; t <= nBlksCMSA(cma); t++){
 		s1 = start[t];
		exists=PosSiteCMSA(t,item1,pos,cma); 
		if(exists) {
                   for(s2=pos[1],s=1; s <= LengthCMSA(t,cma); s++,s1++,s2++){
			if(seq1[s1] == seq2[s2]) score++;
		   }
                }
	   } prcntID = 100.*(double)score/(double)totlen;
	   if(prcntID >= cutoff){ itemList[item1]=TRUE; hits++; }
	} 
	fprintf(stderr,"%d sequences removed\n",hits);
	return itemList;
}

unsigned char	*GetAlnResInSiteCMSA(Int4 t, Int4 sq, cma_typ cma)
// Return an array of the residues aligned with (the first) site t
// in sequence sq. Return null if no sites present
{
	Int4	pos[9];
	unsigned char	*seq;
	if(t < 1 || t > nBlksCMSA(cma)) return 0;
	if(sq < 1 || sq > NumSeqsCMSA(cma)) return 0;
	if(PosSiteCMSA(t,sq,pos,cma) > 0){
		seq= SeqPtrCMSA(sq,cma);
		return (seq+pos[1]-1);
	} else return 0;
}

BooLean	*RepSetCMA(double cutoff, cma_typ cma)
/**************************************************************************
  Purge closely related sequences from a list of sequences.
  WARNING: not yet tested...
 **************************************************************************/
{
	Int4	totlen,i,k,s,t,entry,N,score,max,s1,s2;
	Int4	item,item1,*hits,num_seqs=0;
	b_type	*edges;
	dh_type	H;
	keytyp	key;
	BooLean	*itemList;
	a_type	A=AlphabetCMSA(cma);
	unsigned char	*seq1,*seq2;
	double	prcntID;
	ss_type	data=DataCMSA(cma);
	Int4	pos[4];

	assert(cutoff > 0.0 && cutoff <= 100);
	N=NumSeqsCMSA(cma);
	for(totlen=0,t=1; t <= nBlksCMSA(cma); t++) totlen+=LengthCMSA(t,cma);
	H = dheap(N+2,3);
	NEW(edges,N+2,b_type); NEW(hits,N+2,Int4);
	for(entry =1 ;entry <= N;entry++) edges[entry]=Block(N+1);
	entry = 1;
	for(item =1 ;item <= N;item++) {
	    seq1 = SeqPtrCMSA(item,cma);
	    for(item1 =item+1 ;item1 <= N;item1++) {
		seq2 =  SeqPtrCMSA(item1,cma);
		for(score=0,t=1; t <= nBlksCMSA(cma); t++){
		     PosSiteCMSA(t,item,pos,cma); s1=pos[1];
		     PosSiteCMSA(t,item1,pos,cma); s2=pos[1];
                     for(s=1; s <= LengthCMSA(t,cma); s++,s1++,s2++){
			if(seq1[s1] == seq2[s2]) score++;
                     }
                } 
		prcntID = 100.*(double)score/(double)totlen;
		if(prcntID >= cutoff){
			AddBlock(item,edges[item1]);
			AddBlock(item1,edges[item]);
			hits[item]++; hits[item1]++;
		}
	    } entry++;
	}
	fprintf(stderr,"\n%d items compared; cutoff %.2f\n", entry-1,cutoff); 
	for(entry=1; entry <= N; entry++){
		key = -(keytyp) hits[entry];
		insrtHeap(entry,key,H);
	}
	while(TRUE){
	    if(emptyHeap(H)) print_error("error in RepSetCMA( )");
	    else if(minkeyHeap(H) >= 0) break;
	    else {
	 	item=delminHeap(H);
		for(i=1;i <= N;i++){
		  if(MemberBlock(i,edges[item]) && memHeap(i,H)){
			hits[i]--;
			key = - (keytyp) hits[i];
			chkeyHeap(i,key,H);
		  }
		}
	    }
	}
	num_seqs = ItemsInHeap(H);
	NEW(itemList,N+3,BooLean);
	for(k=0,entry=1; entry <= N; entry++){
	    if(memHeap(entry,H)){ itemList[entry]=TRUE; k++; }
	    else itemList[entry] = FALSE;
	}
	fprintf(stderr,"%d sequences left after purging\n",k);
	for(entry=1; entry <= N; entry++) { NilBlock(edges[entry]); }
	free(hits); free(edges); Nildheap(H);
	return itemList;
}

char	*NewConservedCMA(Int4 *observed, double cutoff, double *freq, a_type A)
/* return the conserved residues in column "observed". */
{
	Int4	i,r,r2,nres,n;
	double	total,p,q,f,low_cut,high_cut;
	char	*conserved;
	BooLean	flag,debug=0;

	assert(cutoff <= 0.0);
	low_cut=cutoff + 1.5;
	high_cut=cutoff - 1.5;
	if(observed != NULL){
	    nres=0;
	    NEW(conserved,nAlpha(A)+2,char); // conserved[r]=0; ---> 'u';
	    for(total=0.0, r=1; r<=nAlpha(A); r++){
		total += (double)observed[r];
	    }
	    if(debug){
		fprintf(stderr,"cutoff = %.2f; low_cut = %.2f\n",cutoff,low_cut);
		for(r = 1; r <= nAlpha(A); r++){
		   if(observed[r] > 0){
			q = freq[r];
			fprintf(stderr,"%c: %d/%d (p=%.1f)\n", 
				AlphaChar(r,A), observed[r], (Int4)total,
				Log10CBP(observed[r], total, q));
   		   }
		} fprintf(stderr,"\n");
	    }
	    /** 1. find clearly conserved residues **/
	    for(f=0.0, n=0, r=1; r <= nAlpha(A); r++){
		if(observed[r] > 2){
		   q = freq[r];
		   p = Log10CBP(observed[r], total, q);
		   if(p <= low_cut) {
			if(debug) fprintf(stderr," %c: (%.1f)[%d/%d]", 
				AlphaChar(r,A), 
				p,observed[r], (Int4)total);
			n += observed[r]; f += q; 
			if(p <= cutoff){
		   		if(p <= high_cut) conserved[r] = 's'; 
				else  conserved[r] = 'm';
			} else conserved[r] = 'w';
			nres++;
		   }
		}
	    }
            /** 2. find weakly conserved related residue pairs **/
            for(n=0, r=1; r < nAlpha(A); r++){
                if(!conserved[r] && observed[r] > 2){
                  for(r2 = 1; r2 <= nAlpha(A); r2++){
                        if(!conserved[r2] && r2 != r && observed[r2] > 2 
                                        && valAlphaR(r,r2,A) > 0){
                           q = freq[r] + freq[r2];
		  	   p = Log10CBP(observed[r] + observed[r2], total, q);
                           if(p <= low_cut) {
			     if(debug) fprintf(stderr," %c %c (%.1f)", 
					AlphaChar(r,A), AlphaChar(r2,A), p);
                             nres+=2;
			     n += observed[r] + observed[r2];
			     f += freq[r] + freq[r2];
			     if(p <= cutoff){
		   		if(p <= high_cut){
                                  conserved[r]='s'; conserved[r2]='s';
				} else { conserved[r]='m'; conserved[r2]='m'; }
			     } else { conserved[r]='w'; conserved[r2]='w'; }
                           }
                        }
                  }
                }
            } 
	    /** 3. find weakly conserved residues related to clear ones **/
            do {
	      flag = FALSE;
	      for(r=1; r <= nAlpha(A); r++){
		if(conserved[r]){
	          for(r2 = 1; r2 <= nAlpha(A); r2++){
		   if(!conserved[r2] && observed[r2] >= 2){
			if(valAlphaR(r,r2,A) >= 0){
		   	   q = freq[r2];
		           p = Log10CBP(observed[r2], (total-n), q/(1.0 - f));
		   	   if(p <= low_cut){ 
				if(debug) fprintf(stderr," %c (%.1f)", 
						AlphaChar(r2,A), p);
				flag = TRUE; nres++;
				conserved[r2] = 'w';
			   }
			}
		   }
		  }
		}
	      }
	    } while(flag);
	    if(nres==0){ free(conserved); return NULL; }
#if 1
	    // 3. reset to total conserved residues if more highly conserved
	    for(q=0.0,r2=0, r = 1; r <= nAlpha(A); r++){
	      if(conserved[r]){ r2 += observed[r]; q += freq[r]; }
	    }
	    p = Log10CBP(r2, total, q);
	    if(p <= low_cut) {
	      for(r=1; r <= nAlpha(A); r++){
		if(conserved[r]){
	   	  if(p <= cutoff){
		     if(p <= high_cut) conserved[r] = 's';
		     else if(conserved[r] != 's') conserved[r] = 'm';
		  } // else if(conserved[r] != 'm' && conserved[r] != 's') conserved[r] = 'm';
		  // else conserved[r] = 'w';
		}
	      }
	    }
#endif
	    return conserved;
	} else return NULL;
}

double *HenikoffWeightsCMSA(cma_typ cma)
{
	Int4	**rawCnts=0,i;
	double *wt=HenikoffWeightsCMSA(cma,&rawCnts);
        for(i=TotalLenCMSA(cma);i > 0;i--){ free(rawCnts[i]); } free(rawCnts); 
	return wt;
}

double *HenikoffWeightsCMSA(cma_typ cma,Int4 ***RawCnts)
{
        double  *wt;
        Int4    i,n,t,p,l,res,pos[5],col,nCol,posit,**rawCnts,*num_aln;
        Int4    nSeq = NumSeqsCMSA(cma);
        Int4    nBlks = nBlksCMSA(cma);
        a_type  A = AlphabetCMSA(cma);
        ss_type SeqSet = DataCMSA(cma);

        NEW(wt,nSeq+2,double); NEW(num_aln,nSeq+2,Int4);
        for(nCol=0,i=1;i<=nBlks;i++) { nCol += LengthCMSA(i,cma); }
        NEWP(rawCnts,nCol+2,Int4);
        for(t=1;t<=nCol;t++){ NEW(rawCnts[t],nAlpha(A)+2,Int4); }
        for(n=1;n<=nSeq;n++){
                col = 0;
                for(t=1;t<=nBlks;t++){
                        assert(PosSiteCMSA(t,n,pos,cma)==1);
                        p = pos[1];
                        for(l=1;l<=LengthCMSA(t,cma);l++){
                            col++; posit = p+l-1;
                            res=SeqP(n,posit,SeqSet);
                            rawCnts[col][res]++;
                        }
                }
        }
        for(col=0,t=1;t<=nBlks;t++){
            for(l=1;l<=LengthCMSA(t,cma);l++){
                 col++;
		 Int4 nResTyp=0;
             	 for(i=1;i<=nAlpha(A);i++) if(rawCnts[t][i] > 0) nResTyp++;
                 for(n=1;n<=nSeq;n++){
                        assert(PosSiteCMSA(t,n,pos,cma)==1);
                        p = pos[1];
                        posit = p+l-1;
                        res = SeqP(n,posit,SeqSet);
			// Aleksandar's original code:
                        // wt[n] += 1./(double)(rawCnts[col][res]*nResTyp);
			// AFN modifications:
                        if(res > 0){
			   wt[n] += (double)nResTyp/(double)(rawCnts[col][res]);
			   num_aln[n]++;  // AFN: ignore gaps.
			}
                 }
            }
        }
        for(n=1;n<=nSeq;n++) wt[n] /= (double)num_aln[n];  //AFN: avg. over aligned res.
        double  max;
        for(max=0.,n=1;n<=nSeq;n++){ if(wt[n] > max) { max = wt[n]; } }
	for(n=1;n<=nSeq;n++) wt[n] /= max; 
	*RawCnts=rawCnts;
        // for(t=0;t<=nCol;t++){ free(rawCnts[t]); } free(rawCnts); 
	free(num_aln);
        return wt;
}

void    ReNameCMSA(char *newname, cma_typ cma)
{
	RenameSeqSet(newname,DataCMSA(cma));
	RenameSeqSet(newname,TrueDataCMSA(cma));
	if(cma->FullSeq) RenameSeqSet(newname,cma->FullSeq);
}

Int4    ConsensusScoreCMSA(cma_typ cma, e_type E, Int4 gapo, Int4 gapx, Int4 sq)
// Find gapped score of sequence E against sequence sq in cma
{
        assert(LenSeq(E) == LengthCMSA(1,cma) && nBlksCMSA(cma)==1);
        if(sq < 1 || sq > NumSeqsCMSA(cma)) return 0;

        Int4    score = 0;
        char    state=' ';
        a_type  A=AlphabetCMSA(cma);
        for(Int4 s=1; s <= LenSeq(E); s++){
                if(IsDeletedCMSA(sq,s,cma)){
                        if(state=='d') score -= gapx;
                        else if(state!= ' ') score -= gapo+gapx;
                        state='d';
                } else {
                    Int4 ins= InsertionCMSA(sq,s,cma);
                    if(ins > 0){
                        score -= gapo; score -= gapx*ins;
                    }
                    Int4 r1 = ResSeq(s,E);
                    Int4 r2 = ResidueCMSA(1,sq,s,cma);
                    score += valAlphaR(r1,r2,A);
                    state='m';
                }
        }
        if(state=='d'){
                score += gapo;
                for(Int4 s=LenSeq(E); s > 0; s--){
                   if(IsDeletedCMSA(sq,s,cma)) score += gapx;
                   else break;
                }
        } return score;
}

double  ResidueDiversityCMSA(FILE *fp,cma_typ cma)
{ return ResidueDiversityCMSA(fp,0,cma); }

double  ResidueDiversityCMSA(register BooLean *skip, register cma_typ cma)
// ignores ambiguity residues
{
     register sst_typ	set;
     register Int4	k,i,sum,s;
     register unsigned char *seq;

     assert(skip); assert(nBlksCMSA(cma) == 1);
     for(sum=0,k=LengthCMSA(1,cma); k > 0; k--){
          for(set=0,i=NumSeqsCMSA(cma); i > 0; i--){
	      if(!skip[i]){
                seq=SeqPtrCMSA(i,cma);
                s = SitePos(1,i,1,SitesCMSA(cma))+k-1;
                set=UnionSset(set,SsetLet(seq[s]));
	      }
          } for(i=nAlpha(AlphabetCMSA(cma)); i > 0; i--) if(MemSset(i,set)) sum++;
     } return (double)sum/(double) TotalLenCMSA(cma);
}

double  ResidueDiversityCMSA(FILE *fp,BooLean *skip, cma_typ cma)
{
        Int4    t,k,sq,r1,r2,ntypes;
        double  sum;
        BooLean found[50];
        a_type A = AlphabetCMSA(cma);
        h_type HG=0;

        if(fp) HG=Histogram("diversity of residue types",0,nAlpha(A)+3,1.0);
     for(sum=0.0,t=1; t <= nBlksCMSA(cma); t++){
        for(k=1; k <= LengthCMSA(t,cma); k++){
          ntypes=0; found[0]=TRUE;      // ignore ambiguity and insert residues.
          for(r1=1; r1 <= nAlpha(A); r1++) found[r1]=FALSE;
          for(sq=1; sq <= NumSeqsCMSA(cma); sq++){
	      if(skip!=0 && skip[sq]) continue;
              r1=ResidueCMSA(t,sq,k,cma);
              if(!found[r1]) { found[r1]=TRUE; ntypes++; }
          }
          if(HG) IncdHist(ntypes,HG);
          sum += ntypes;
        } 
     }
	if(HG){ PutHist(stdout,60,HG); NilHist(HG); }
        sum /= (double) TotalLenCMSA(cma);
        // fprintf(stdout,"Average diversity = %.2f\n\n",sum);
        return sum;
}

static Int4 pseudo_aln_score_cmsa(register Int4 Length, register unsigned char *seq1,
        register unsigned char *seq2, register char **R)
{
        // register Int4 k,score=0;
        // for(k = 0; k < Length; k++){ score += R[seq1[k]][seq2[k]]; }
#if 1
        register Int4 score;
        for(score=0; Length > 0; Length--,seq1++,seq2++){
                // Skip over gap residues...
                if(*seq1 && *seq2) score += R[*seq1][*seq2];
                // score += R[*seq1][*seq2];
        }
#else	// Find local optimal pseudo alignment score.
        register Int4 k=Length,score=0;
	while(R[seq1[k]][seq2[k]] < 0) k--; 
	Length=k; k=1; 
	while(R[seq1[k]][seq2[k]] < 0){ k++; }
	while(k <= Length){ score += R[seq1[k]][seq2[k]]; k++; }
#endif
	return score;
}

Int4    PseudoAlnScoreCMSA(Int4 sq1, Int4 sq2, cma_typ cma)
// use cma alignment to obtain pseudo pairwise scores for two aligned sequences.
{
	// assert(!"THIS ROUTINE NEEDS TO BE DEBUGGED!\n");
	// NOTE: OFF BY ONE!!!???
        assert(nBlksCMSA(cma) ==1);
        Int4 N=NumSeqsCMSA(cma);
        assert(sq1 >0 && sq1 <= N && sq2 > 0 && sq2 <= N);
        a_type  A=AlphabetCMSA(cma);
        unsigned char *seq1 = SeqSeqSet(sq1,DataCMSA(cma));
        Int4 s1 = SitePos(1,sq1,1,SitesCMSA(cma));
        unsigned char *seq2 = SeqSeqSet(sq2,DataCMSA(cma));
        Int4 s2 = SitePos(1,sq2,1,SitesCMSA(cma));
        return pseudo_aln_score_cmsa(LengthCMSA(1,cma),seq1+s1,seq2+s2,AlphaR(A));
}

Int4	PseudoAlnScoreTwoCMSA(Int4 sq1, cma_typ cma1, Int4 sq2, cma_typ cma2)
// use cma alignment to obtain pseudo pairwise scores for two aligned sequences.
{
        Int4 score;
        assert(nBlksCMSA(cma1) ==1); assert(nBlksCMSA(cma2) ==1);
	assert(LengthCMSA(1,cma1) == LengthCMSA(1,cma2));
        Int4 N1=NumSeqsCMSA(cma1),N2=NumSeqsCMSA(cma2);
        assert(sq1 >0 && sq1 <= N1 && sq2 > 0 && sq2 <= N2);
        a_type  A=AlphabetCMSA(cma1);
        unsigned char *seq1 = SeqSeqSet(sq1,DataCMSA(cma1));
        Int4 s1 = SitePos(1,sq1,1,SitesCMSA(cma1));
        unsigned char *seq2 = SeqSeqSet(sq2,DataCMSA(cma2));
        Int4 s2 = SitePos(1,sq2,1,SitesCMSA(cma2));
        // return pseudo_aln_score_cmsaX(LengthCMSA(1,cma),seq1+s1,seq2+s2,AlphaR(A));
	Int4 Length = LengthCMSA(1,cma1);
	seq1 += s1; seq2 += s2;
	char **R = AlphaR(A);
        for(score=0; Length > 0; Length--,seq1++,seq2++){
		// Skips over gap residues...
                if(*seq1 && *seq2) score += R[*seq1][*seq2];
        } return score;
}

//*******************************************************************************
static Int4 pseudo_aln_score_sq2cma(register Int4 Length, register unsigned char *seq1,
        	register unsigned char *seq2, register char **R)
{
        register Int4 score;
        for(score=0; Length > 0; Length--,seq1++,seq2++){
                // Skip over gap residues...
                if(*seq1 && *seq2) score += R[*seq1][*seq2];
        }
        return score;
}

Int4    PseudoAlnScoreSqToCMSA(e_type E, Int4 sq2, cma_typ cma)
// use cma alignment to obtain pseudo pairwise scores for two aligned sequences.
{
        assert(nBlksCMSA(cma) ==1);
        assert(LenSeq(E) == LengthCMSA(1,cma));
        Int4 N=NumSeqsCMSA(cma);
        assert(sq2 > 0 && sq2 <= N);
        a_type  A=AlphabetCMSA(cma);
        unsigned char *seq1 = SeqPtr(E);
        unsigned char *seq2 = SeqSeqSet(sq2,DataCMSA(cma));
        Int4 s2 = SitePos(1,sq2,1,SitesCMSA(cma));
        return pseudo_aln_score_sq2cma(LengthCMSA(1,cma),seq1+1,seq2+s2,AlphaR(A));
}
//*******************************************************************************

static double pseudo_aln_score_cmsa2(register Int4 Length, register unsigned char *seq1,
        register unsigned char *seq2, register char **R,register Int4 *Rank)
// modified from cmsa.c file...
{
        // register Int4 k,score=0;
        // for(k = 0; k < Length; k++){ score += R[seq1[k]][seq2[k]]; }
        register Int4 s;
        register double score;
        for(score=0; Length > 0; Length--,seq1++,seq2++,Rank++){
                // Skip over gap residues...
                s = R[*seq1][*seq2];
                if(s > 0) score+=10*s*(1.0/(double)*Rank);
                // if(*seq1 && *seq2) score += R[*seq1][*seq2];
                // score += R[*seq1][*seq2];
        }
        return score;
}

double  PseudoAlnScoreCMSA(unsigned char *qseq, Int4 sq, Int4 *Rank, cma_typ cma)
// use cma alignment to obtain pseudo pairwise scores for two aligned sequences.
{
        // assert(!"THIS ROUTINE NEEDS TO BE DEBUGGED!\n");
        // NOTE: OFF BY ONE!!!???
        assert(nBlksCMSA(cma) ==1);
        Int4 N=NumSeqsCMSA(cma);
        assert(sq > 0 && sq <= N);
        a_type  A=AlphabetCMSA(cma);
        unsigned char *sseq = SeqSeqSet(sq,DataCMSA(cma));
        Int4 s2 = SitePos(1,sq,1,SitesCMSA(cma));
        return pseudo_aln_score_cmsa2(LengthCMSA(1,cma),qseq,sseq+s2,AlphaR(A),Rank);
}

e_type	**GetRepSetCMSA(FILE *fp, Int4 percent_ident,Int4 *Nset,cma_typ cma)
// return a representative set of sequences from cma.
{
	ds_type sets;
	Int4	b,i,j,k,N=NumSeqsCMSA(cma);
	Int4	score,nblk=nBlksCMSA(cma);
	Int4	s,si,sj,pos[3],seti,setj;
	ss_type	data=TrueDataCMSA(cma);
        a_type  A=AlphabetCMSA(cma);
	unsigned char	*isq,*jsq;

	assert(percent_ident > 0 && percent_ident <= 100);
	// 1. Cluster sequences into sets at the percent identity cutoff.
	Int4 total = TotalLenCMSA(cma);
	sets = DSets(N);
	for(i = 1; i < N; i++){
	  isq = SeqPtrCMSA(i,cma);
	  seti=findDSets(i,sets);
	  for(j=i+1; j <= N; j++){
	     setj=findDSets(j,sets);
	     if(seti != setj){
	      jsq = SeqPtrCMSA(j,cma);
	      for(score=0,b=1; b <= nblk; b++){
		PosSiteCMSA(b,i,pos,cma); si=pos[1];
		PosSiteCMSA(b,j,pos,cma); sj=pos[1];
		for(s=1; s <= LengthCMSA(b,cma); s++,si++,sj++){
			if(isq[si] == jsq[sj]) score++;
		}
	      } score = (Int4) floor(((double)score*100.0/(double)total) +0.5);
	      if(score >= percent_ident) seti=linkDSets(seti,setj,sets);
	     }
	  }
	}
	// 2. Within each set pick the sequence with the highest profile score.
	e_type	**ListE;
	NEWP(ListE, NumSeqsCMSA(cma)+2, e_type);
        for(s=1,i=1; i <= N; i++){
	  seti=findDSets(i,sets);
	  if(i==seti){	// is this the canonical sequence?
	     NEW(ListE[s], NumSeqsCMSA(cma)+2, e_type);
	     Int4 sq=0;
	     ListE[s][sq] = SeqSetE(i,data); sq++;
	     for(j=1; j <= N; j++){
		if(j!=i){
		  setj=findDSets(j,sets);
		  if(seti == setj){
			ListE[s][sq] = SeqSetE(j,data); sq++;
		  }
		}
	     } s++;
	  }
	} *Nset=s-1;
	// 3. return the sequence sets 
	NilDSets(sets);
	return ListE;
}

Int4	*AlignSeqCMSA(FILE *ftpr,Int4 *Score, e_type  E,cma_typ cmsa)
// aligns a sequence against an fmodel of the input cmsa object.
{
	Int4	n,t,ntyps,score;
	Int4	*pos;
	smx_typ	*smx;
	double	pernats=1000.0;

	ntyps=nBlksCMSA(cmsa);
	NEW(smx,ntyps+1,smx_typ);
	NEW(pos,ntyps+1,Int4);
	for(t=1; t<=ntyps; t++){
		smx[t]=GetSmatrixFModel(pernats,ModelCMSA(t,cmsa));
	}
	score=AlnSeqSMatrix(LenSeq(E),SeqPtr(E),ntyps,smx,pos);
	// score=AlnSeqSMatrix(LenSeq(E),XSeqPtr(E),ntyps,smx,pos);
	if(ftpr) PutSeqInfo2(ftpr,E);
	for(t=1; t<=ntyps; t++){
		if(ftpr){
		       	PutSeqRegion(ftpr,pos[t],LengthCMSA(t,cmsa),E,AlphabetCMSA(cmsa));
			fprintf(ftpr,"\n");
		} NilSMatrix(smx[t]);
	} if(ftpr) fprintf(ftpr,"score = %d\n",score);
	free(smx);
	*Score=score;
	return pos;
}

Int4    NumInsertsCMSA(Int4 blk, Int4 site, cma_typ cma, Int4 &ins_res)
// return number of insertions following column 'site'.
{
        Int4    sq,n,nins=0;
        for(ins_res=0,sq=1; sq <= NumSeqsCMSA(cma); sq++){
#if 1
                if(!SeqIsInCMSA(sq,cma)) continue;
#else
                assert(SeqIsInCMSA(sq,cma));
#endif
                n = InsertionCMSA(blk,sq,site,cma);
                if(n > 0){ nins++; ins_res += n; }
        } return nins;
}

