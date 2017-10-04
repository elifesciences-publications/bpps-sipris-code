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

#include "hat_typ.h"

hat_typ::hat_typ(FILE *fp)
{ 
	Init(fp);
}

void	hat_typ::Init(FILE *fp)
{

}

void	hat_typ::Free( )
{

}

void    hat_typ::hat_srch(FILE *fp)
// search a set of sequences 
{
}

void    hat_typ::Put(FILE *fp)
{
}

void	TrimToTemplateCMSA(cma_typ *IN_CMA, Int4 num_cma_files)
// Check for truncations at either end and modify cma files accordingly.
{
   cma_typ tpl_cma=IN_CMA[0];
   for(Int4 s=1; s <= num_cma_files; s++){ 
        fprintf(stderr,"================================ %d: %s.\n", s,NameCMSA(IN_CMA[s])); 
	e_type tplSq=TrueSeqCMSA(s+1,tpl_cma);
	if(LenSeq(tplSq) > LengthCMSA(1,IN_CMA[s])) print_error("TrimToTemplateCMSA() input error 1");
	if(LenSeq(tplSq) != LengthCMSA(1,IN_CMA[s])){	// then need to change input cma file.
		Int4 Start,i;
		cma_typ cmaX=IN_CMA[s];
		e_type csqSq=TrueSeqCMSA(1,cmaX);
		char rtn=IsSubSeq(tplSq,csqSq,&Start,FALSE);
		// rtn = 1 if tplSq is a subseq of csqSq.
		if(rtn != 1) print_error("Template and cma files are incompatible");
		if(Start > 0){					// remove N-terminal columns.
			for(i=1; i<=Start; i++){
                           if(LengthCMSA(1,cmaX) <= 3) print_error("TrimToTemplateCMSA() input error 2");
                           RmColumnMSA(1,1,cmaX); // block 1, first column removed.
                        }
		}
		if(LenSeq(tplSq) < LengthCMSA(1,cmaX)) {	// remove C-terminal columns.
			Int4 lenrm = LengthCMSA(1,cmaX) - LenSeq(tplSq);
                        for(i=1; i<=lenrm; i++){
                           Int4 lemon = LengthCMSA(1,cmaX);
                           if(lemon <= 3) print_error("TrimToTemplateCMSA() input error 3");
                           RmColumnMSA(1, lemon, cmaX);
                        }
		}
		IN_CMA[s] = MinimizeFirstSeqCMSA(cmaX); TotalNilCMSA(cmaX);
	}
   }
}

#if 0	// pasted this here for testing...can delete (afn: 10-3-12).
st_type SplitElementSites3(Int4 x, Int4 len_left, Int4 minlen, st_type S)
// From sites.cc (need to put it back there....
/*****************************************************************
 if element x is >= minlen + len_left then split it into two blocks;
 create and return the resultant object S2.
 *****************************************************************/
{
    Int4        t,n,s,N,T,T2,*len_elem,len_right,leftleng,rightleng;
    st_type     S2;
    double      r;

    // 0. Test for input errors
    T = nTypeSites(S); T2=T+1;
    if(x < 1 || x > T) return NULL;
    len_right = SiteLen(x,S) - len_left;
    if(len_right < minlen || len_right < 1 || len_left < 1){
	fprintf(stderr,"len = %d < minlen = %d || len = %d < 1 || len_left = %d < 1?\n",
		len_right,minlen,len_right,len_left);
	fprintf(stderr,"SiteLen(%d,S) = %d.\n",x,SiteLen(x,S));
	return NULL;
    } leftleng = len_left; rightleng = len_right;

    // 1. Create new sites object  (same as below; eventually merge...)
    NEW(len_elem, T2+3, Int4);
    for(t=1; t <= T2; t++) {
        if(t < x) len_elem[t]=SiteLen(t,S);
        else if(t == x) len_elem[t]=leftleng;
        else if(t == (x+1)) len_elem[t]=rightleng;
        else len_elem[t]=SiteLen(t-1,S);
    }
    S2=MakeSites(T2,len_elem, S->gss);
    N=NSeqsSites(S2);
    free(len_elem);

    // 2. Add sites to new object
    for(n=1; n<=N; n++){
        for(t=1; t <= T2; t++) {
          if(t <= x) AddSite(t,n,SitePos(t,n,1,S),S2);
          else if(t == (x+1)) AddSite(t,n, leftleng+SitePos(x,n,1,S),S2);
          else AddSite(t,n,SitePos(t-1,n,1,S),S2);
        }
    }
    return S2;
}

cma_typ SplitBlkCMSA3(Int4 x, Int4 left_leng, Int4 minlen, cma_typ L)
// From cmsa_recombine.cc; need to put it back eventually...
/***********************************************************************
 Splits a block into left_leng making sure that right_leng >= minlen.
 ***********************************************************************/
{
    cma_typ     ma;
    Int4        j,t2,t,T,T2;
    st_type     S,S2=NULL;
    fm_type     *model;
    BooLean     **null;

    S=SitesCMSA(L); T = nTypeSites(S); T2=T+1;
    if(x < 1 || x > T) return NULL;
    S2 = SplitElementSites3(x, left_leng, minlen, S);
    // rest same as below... eventually merge...
    if(S2 == NULL) return NULL;
    model=ModelsCMSA(L);
    NEWP(null,T2+2,BooLean);
    if(x==0) t2=2; else t2=1;
    for(t=1,j=1; j<= T; j++,t2++,t++){
        if(j==x){ t2++; continue; } // use all columns in split block
        NEW(null[t2],LenFModel(model[t])+5,BooLean);
        NullSitesFModel(null[t2], model[t]);
    }
    ma=MakeCMSA(S2,null);
    for(t=1; t<=T2; t++) if(null[t]!=NULL)free(null[t]);
    free(null);
    return ma;
}

cma_typ InsertColumnsCMA3(Int4 start, Int4 N, cma_typ cma)
// add n new columns just to the right of start ...
{
        Int4 min_len=1;
        assert(nBlksCMSA(cma) == 1);

        // 1. Split block right after start point for inserting residues.
        // fprintf(stderr,"1. Split block right after start point for inserting residues.\n");
        cma_typ cma0=SplitBlkCMSA3(1,start,min_len,cma);
        if(cma0 == 0){
          FILE *fp = open_file("debug_error",".cma","w"); PutCMSA(fp,cma); fclose(fp);
          assert(cma0 != 0);
          print_error("Error in InsertColumnsCMA3()");
        } // WriteCMSA("debug0.cma", cma0);

        // 2. Append n columns to first block right after start of insert
        // fprintf(stderr,"2. Append %d columns to first block right after start of insert\n",N);
        BooLean right=TRUE;
        for(Int4 col=1; col<= N; col++){
                if(InsertColCMSA(1,right,cma0) == FALSE)
                        fprintf(stderr,"InsertColCMSA() inserted gaps\n");
        } // WriteCMSA("debug1.cma", cma0);

        // 3. Fuse blocks again.
        // fprintf(stderr,"3. Fuse blocks again\n");
        cma_typ cma1=SimpleFuseBlksCMSA(1, 10000, cma0);
        if(cma1 == 0) print_error("Fuse error in InsertColumnsCMA3()");
        // WriteCMSA("debug2.cma", cma1);
        // NilCMSA(cma);
        NilCMSA(cma0);
        return cma1;
}
#endif	// end of inserted code...can delete.

/********************************************************************** 

{()MVIRPATPELLRDGELDAVLMVKLLL()}*	template root consensus.
{()---NVTgnevlAPA--VVAP--NTLFYKV--()}*	template subgroup consensus.

   {()NVTGNEVLAPAVVAPNTLFYKV()}*	subgroup consensus.

  Need to add deletion ('-') columns at ends to  match template...
 **********************************************************************/
cma_typ	GrowEdgesCMSA(BooLean add2Nterm, Int4 NumGaps, cma_typ cma)
// Trim the cma input file at ends 
{
    Int4        i,j;

// fprintf(stderr,"st = %d\n",NumGaps); PutCMSA(stderr,cma);
//     ExtendFakeToRealCMSA(cma);	// else need to fix GetEdgeBlocksSite() routine!!!
// fprintf(stderr,"st = %d\n",NumGaps); PutCMSA(stderr,cma);
    for(i=1; i <= NumGaps; i++){
      if(add2Nterm){ if(!InsertColCMSA(1,FALSE,cma)) return 0; }
      else { if(!InsertColCMSA(1,TRUE,cma)) return 0; }
// PutCMSA(stderr,cma);
    } return cma;
}

BooLean	IsOkayTemplateCMA(cma_typ TemplateCMA)
{
	// fprintf(stderr,"len(cma %d)=%d\n",II,LengthCMSA(1,cma));
	for(Int4 sq=1; sq <= NumSeqsCMSA(TemplateCMA); sq++){
	   gsq_typ *gsq=gsqCMSA(sq,TemplateCMA);
	   // check for insertions at N-terminus.
	   Int4 i;
	   if((i=gsq->Insertion(0)) > 0){
		// gsq->Put(stderr,AlphabetCMSA(TemplateCMA));
		fprintf(stderr,
			"**** Fatal: %d residue insert at N-terminus of template sequence %d: ****\n",
			i,sq);	
		gsq->Put(stderr,60,AlphabetCMSA(TemplateCMA));
		print_error("********** Template file formating error *********");
		return FALSE;
	   }
	   // check for insertions at C-terminus.
	   if((i=gsq->CheckForInsertsAtEnds( )) > 0){
		// gsq->Put(stderr,AlphabetCMSA(TemplateCMA));
		fprintf(stderr,
			"**** Fatal: %d residue insert at C-terminus of template sequence %d: ****\n",
			i,sq);	
		gsq->Put(stderr,60,AlphabetCMSA(TemplateCMA));
		print_error("********** Template file formating error *********");
		return FALSE;
	   }
	}
	return TRUE;
}

cma_typ	*ConversionViaTemplateCMSA(cma_typ TemplateCMA, cma_typ *IN_CMA)
//************************************************************************************
//            NEW ROUTINE TO CONVERT INPUT CMAFILES GIVEN A TEMPLATE.
//************************************************************************************
{ 
#if 1	// replacing old routines; delete below after testing (afn: 8-13-2015).
 	return ConversionViaTemplateCMSA3(TemplateCMA, IN_CMA);
#else
	Int4	II,num_cma_files=NumSeqsCMSA(TemplateCMA)-1;
	e_type	trueE,fakeE;
	cma_typ	*OUT_CMA,cma;

   NEW(OUT_CMA,num_cma_files+3,cma_typ);
   for(II=1; II <= num_cma_files; II++){	
	if(IN_CMA[II]==0) continue;
	cma=0;		// CMA[1] is template cmafile.
	// fprintf(stderr,"\n%s:\n",NameCMSA(IN_CMA[II]));
	Int4 len_template=LengthCMSA(1,TemplateCMA);
	Int4 len_cma=LengthCMSA(1,IN_CMA[II]);
	Int4 st,sc,seq,ins,ndel=0;
	trueE = TrueSeqCMSA(II+1,TemplateCMA);
	fakeE = FakeSeqCMSA(1,IN_CMA[II]);
	if(!IdentSeqs(trueE,fakeE)){
	    fprintf(stderr,"\n%s:\n",NameCMSA(IN_CMA[II]));
	    if(!IdentSeqs(trueE,fakeE)){
		PutSeq(stderr,trueE,AlphabetCMSA(TemplateCMA));
		PutSeq(stderr,fakeE,AlphabetCMSA(TemplateCMA));
		fprintf(stderr,"Fatal: trueE and fakeE don't match in cma file %d\n",II+1);
		exit(1);
	    }
	}
	assert(IdentSeqs(TrueSeqCMSA(1,IN_CMA[II]),fakeE));	
	assert(LenSeq(trueE) == len_cma);	// sequence 
	if(IsDeletedCMSA(1,II+1,1,TemplateCMA)){	// no deletions allowed on ends.
		fprintf(stderr,"Fatal: deletion at first position in cma file %d\n",II+1);
		exit(1);
	}
	if(IsDeletedCMSA(1,II+1,len_template,TemplateCMA)){
		fprintf(stderr,"Fatal: deletion at last position in cma file %d\n",II+1);
		// exit(1);
	}
	if(InsertionCMSA(1,II+1,len_template,TemplateCMA) > 0){
		fprintf(stderr,"Fatal: insertion within cma file %d\n",II+1);
		exit(1);
	}
	cma=CopyCMSA(IN_CMA[II]);
	// sc = site in cma
	// st = site in template
	for(st=len_template, sc=len_cma; st > 1; ){
	   if(sc < 1){
	      fprintf(stderr,"sc < 1: template %d, site %d(%d)\n", II,st,sc);
	      assert(sc > 0);
	   }
	   seq=1;
	   ins=InsertionCMSA(1,II+1,st,TemplateCMA);

	   //********  CASE 1: Deletion & insertion in template alignment.
	   if(IsDeletedCMSA(1,II+1,st,TemplateCMA) && ins > 0){
	     ndel=0;
             do {
                if(st <= 1) break;
                ndel++; st--;
             } while(IsDeletedCMSA(1,II+1,st,TemplateCMA));

#if 1		// Print out changes in cma at each step...
		fprintf(stderr,"***** template %d: site %d(%d) *****\n", II,st,sc);
		gsq_typ *gsq=gsqCMSA(II+1,TemplateCMA);
		gsq->Put(stderr,60,AlphabetCMSA(TemplateCMA));

		fprintf(stderr,"************* subject seq. ******************\n");
		gsq=gsqCMSA(1,cma);
		gsq->Put(stderr,60,AlphabetCMSA(cma));
#endif

		// this is the tricky part...
		ins= ins - ndel;
		if(ins > 0){
		 sc-=ins; 
		 ColumnsToInsertCMSA(cma,sc+1,sc+ins); // start with column right after sc.
		 sc--; 
		 st--;
		} else if(ins < 0){
		 st = st - ins; // move up to start.
		 sc = sc - ins; // move up to start.
		}
	   	assert(sc > 0);
#if 1
		gsq=gsqCMSA(1,cma);
		// gsq->Put(stderr,AlphabetCMSA(cma));
		gsq->Put(stderr,60,AlphabetCMSA(cma));
		// exit(1);
#endif
	   } 
	   //******** CASE 2: Deletion only in template alignment.
	   else if(IsDeletedCMSA(1,II+1,st,TemplateCMA)){
		st--;
		// fprintf(stderr,"template %d: ndel= %d at site %d(%d)\n", II,ndel,st,sc);
	        cma_typ cma0=InsertColumnsCMA(sc,1,cma); // add '-'s to right of sc.
		if(cma0==0) print_error("InsertColumnsCMA() error");
		cma=cma0;
	   } 

	   //******** CASE 3: Insertion only in template alignment.
	   else if(ins > 0){			// Insertion in template alignment.
		// fprintf(stderr,"template %d: insert = %d at site %d(%d)\n", II,ins,st,sc);
		sc-=ins; 
		ColumnsToInsertCMSA(cma,sc+1,sc+ins); // start with column right after sc.
		st--; sc--; 
	   	assert(sc > 0);
	   } 
	   //******** CASE 4: Match in template alignment.
	   else { st--; sc--; }	
	}
	OUT_CMA[II]=cma;
   } return OUT_CMA;
#endif
}

#if 0
BooLean IsDeletedCMSA(UInt4 n, UInt4 r, cma_typ cma);
BooLean IsDeletedCMSA(UInt4 blk, UInt4 n, UInt4 r,
        cma_typ cma);
unsigned short InsertionCMSA(UInt4 blk,UInt4 n,UInt4 i,
        cma_typ cma);
unsigned short InsertionCMSA(UInt4 n,UInt4 i,cma_typ cma);
Int4    NumColumnsCMSA(cma_typ msa);
Int4    TotalLenCMSA(cma_typ msa);
Int4    FakeToRealCMA(Int4 sq,Int4 s, cma_typ cma);
Int4    ResidueCMSA(register Int4 t, register Int4 n, register Int4 s, cma_typ cma);
#endif

#if 1	// Add routine for outputting regions in vsi format.
void	PutVSIregionsCMSA(FILE *fp,Int4 sq,char *color,Int4 *start,Int4 *end,
		Int4 num_regions,cma_typ cma)
{
	Int4	r,i,c,n;
	ss_type	data=DataCMSA(cma);
	a_type	AB=SeqSetA(data);

	n=1;
	assert(nBlksCMSA(cma) ==1);
	c=TotalLenCMSA(cma);
	r = FakeToRealCMA(sq,c,cma);
	fprintf(fp,"1-%d.W20\n",r);
	for(c=1; c <= TotalLenCMSA(cma); c++){
	   if(n <= num_regions){
	     r = FakeToRealCMA(sq,c,cma);
	     if(c == start[n]){
		if(IsDeletedCMSA(sq,c,cma)){
			fprintf(fp,"%d-",r+1);
		} else fprintf(fp,"%d-",r);
	     } else if(c == end[n]){
		if(IsDeletedCMSA(sq,c,cma)){
			fprintf(fp,"%d.%c\n",r,color[n]);
		} else fprintf(fp,"%d.%c\n",r,color[n]);
		n++;
	     } else if(c > start[n] && c < end[n]){
		if((i=InsertionCMSA(sq,c,cma))){
		     printf("%d.%c,",r,color[n]);
	     	     r = FakeToRealCMA(sq,c+1,cma);
		     printf("%d-",r);
		}
	     } else if(c > end[n]) n++;
	   }
	} fprintf(fp,"\n\n");
#if 0	// DEBUG.
	for(c=1; c <= TotalLenCMSA(cma); c++){
		if(IsDeletedCMSA(sq,c,cma)){
			printf("%d(%d): -\n",c,r);
		} else {
			Int4 aa=ResidueCMSA(1, sq, c, cma);
			printf("%d(%d): %c\n",c,r,AlphaChar(aa,AB));
		} 
		if((i=InsertionCMSA(sq,c,cma))){
			printf("%d(%d): (%d)\n",c,r,i);
		}
	} fprintf(fp,"\n\n");
#endif
}
#endif

Int4	ParseInputRegionsVSI(char *input_string,Int4 **Start,Int4 **End,
		char **Colors,const char *usage)
	// -m=R:6..20;O:32..45;Y55..68.
// input_string = "R:6..20,O:32..45,Y:55..68."
{
	Int4	n,i,start,end,*S,*E,NumRegions=0;
	char	color,*C,*str;

	for(i=0; input_string[i]; i++){
		if(input_string[i] == ':') NumRegions++;
	}
	
	NEW(S,NumRegions+3,Int4); NEW(E,NumRegions+3,Int4);
	NEW(C,NumRegions+3,char);
	n=0;
	str=input_string;
	for(i=1; i < NumRegions; i++){
	  if(sscanf(str,"%c:%d..%d,",&color,&start,&end) != 3){
		print_error(usage);
	  } else { S[i]=start; E[i]=end; C[i]=color; }
	  while(str[0]!=','){ str++; } str++;
	}
	if(sscanf(str,"%c:%d..%d.",&color,&start,&end) != 3){
		print_error(usage);
	} else { S[i]=start; E[i]=end; C[i]=color; }
	// Check start and end...
	for(i=1; i < NumRegions; i++){
		if(E[i] < S[i]) print_error(usage);
		if(S[i+1] < E[i]) print_error(usage);
	} if(E[i] < S[i]) print_error(usage);
	*Start=S; *End=E; *Colors=C;
	return NumRegions;
}

#if 0	// not used...
char	*AddInsertToOperationArray3(Int4 start, Int4 end, char *operation)
// ========== Add an insertion to an operational array. ===========
{
	
	char state,*new_operation=0;
	Int4 j,o,no,column;
	Int4 trace_length=strlen(operation);
	NEW(new_operation,trace_length+5,char);
	new_operation[0]='E'; 

	for(no=1,o=j=1,column=1,state='E'; operation[o] != 'E'; o++){
          switch(operation[o]){
            case 'M': 
            case 'm': 
		if(column >= start && column <=end){
			   new_operation[no]='I'; 
		} else new_operation[no]=operation[o];
		no++; j++; column++; 
		break;
            case 'D': 
            case 'd': // deletion in sequence relative to profile.
		if(column >= start && column <=end){
			   // do nothing in new_operation; 
		} else { new_operation[no]=operation[o]; no++; }
                column++; break;
            case 'i': // insert is between profile blocks;
		new_operation[no]=operation[o]; no++; 
		j++; break;
            case 'I': // Insert ('-') within a profile block; delete from seq.
		new_operation[no]=operation[o]; no++; 
		j++; break;
            default:
            print_error("operation( ): input error"); break;
          }  state=operation[o];
       	}
	new_operation[no]='E'; no++; new_operation[no]=0;
	return new_operation;
}
#endif

cma_typ	*ConversionViaTemplateCMSA3(cma_typ TemplateCMA, cma_typ *IN_CMA)
//************************************************************************************
//            NEW ROUTINE TO CONVERT INPUT CMAFILES GIVEN A TEMPLATE.
// any changes to this code should be copied to the corresponding files in hat_typ.cc
//************************************************************************************
{ 
	Int4	II,num_cma_files=NumSeqsCMSA(TemplateCMA)-1;
	e_type	trueE,fakeE;
	cma_typ	*OUT_CMA,cma,cma0;
	a_type	AB=AlphabetCMSA(TemplateCMA);

   NEW(OUT_CMA,num_cma_files+3,cma_typ);
   for(II=1; II <= num_cma_files; II++){	
	if(IN_CMA[II]==0) continue; // else fprintf(stderr,"\n%d.%s:\n",II,NameCMSA(IN_CMA[II]));
#if 0	// DEBUG
	{
	  char    str[100]; sprintf(str,"_aln%d.cma",II);
	  FILE *fptr=open_file("debug",str,"w"); PutCMSA(fptr,IN_CMA[II]); fclose(fptr);
	}
#endif
	cma=0;		// CMA[1] is template cmafile.

	//************* 1. Confirm consistency between the template and input alignments. ***************
	Int4 len_template=LengthCMSA(1,TemplateCMA);
	Int4 len_cma=LengthCMSA(1,IN_CMA[II]);
	Int4 st,sc,ins,ndel=0;
	trueE = TrueSeqCMSA(II+1,TemplateCMA); fakeE = FakeSeqCMSA(1,IN_CMA[II]);
	if(!IdentSeqs(trueE,fakeE)){
	    fprintf(stderr,"\n%s:\n",NameCMSA(IN_CMA[II]));
	    PutSeq(stderr,trueE,AB); PutSeq(stderr,fakeE,AB);
	    fprintf(stderr,"Fatal: trueE and fakeE don't match in cma file %d\n",II+1); exit(1);
	}
	assert(IdentSeqs(TrueSeqCMSA(1,IN_CMA[II]),fakeE));	
	assert(LenSeq(trueE) == len_cma);	// sequence 
	if(InsertionCMSA(1,II+1,len_template,TemplateCMA) > 0){
		fprintf(stderr,"Fatal: insertion within cma file %d\n",II+1); exit(1);
	}

	//********** 2. If deletions on either end of template, find C- and N-terminal adjustments. ***********
	cma=CopyCMSA(IN_CMA[II]);
// PutCMSA(stderr,cma);	// DEBUG...
	Int4 Nt_adj=0,Ct_adj=0;
	// This solves issues with deletions on the ends. Need to convert GK{()---  ->  {(GK)---- .
	if(IsDeletedCMSA(1,II+1,1,TemplateCMA) || IsDeletedCMSA(1,II+1,len_template,TemplateCMA)){
        	ExtendFakeToRealCMSA(cma);	// else need to fix GetEdgeBlocksSite() routine!!!
	}
	if(IsDeletedCMSA(1,II+1,1,TemplateCMA)){	// no deletions allowed on ends.
		fprintf(stderr,"Warning: first position deleted in template alignment %d (File %d)\n",II+1,II);
		for(st=0; IsDeletedCMSA(1,II+1,st+1,TemplateCMA); st++) ;  // move in until reach a residue.
		if(GrowEdgesCMSA(TRUE,st,cma) == 0){
			WriteCMSA("debug1.cma", cma);
			print_error("hat_typ GrowEdgesCMSA( ) error");
		} Nt_adj=st;
	}
	if(IsDeletedCMSA(1,II+1,len_template,TemplateCMA)){
		fprintf(stderr,"Warning: last position deleted in template alignment %d (File %d)\n",II+1,II);
		for(st=0; IsDeletedCMSA(1,II+1,len_template-st,TemplateCMA); st++) ;
		if(GrowEdgesCMSA(FALSE,st,cma) == 0){
			WriteCMSA("debug2.cma", cma);
			print_error("hat_typ GrowEdgesCMSA( ) error");
		} Ct_adj=st;
	}
	len_cma=LengthCMSA(1,cma);  // Ct_adj=Nt_adj=0;
// PutCMSA(stderr,cma);	// DEBUG...
	//*********** 3. Iterate through full lengths of both template and input alignments. ***********
	// sc = site in cma;  st = site in template
	// Start from the C-terminal end of both template and cma alignment.
	for(st=len_template-Ct_adj,sc=len_cma-Ct_adj; st > Nt_adj; )
	{
	   assert(sc > Nt_adj);
	   //******************** 3a. Disallow insertions directly following deletions. ***************
	   if(IsDeletedCMSA(1,II+1,st,TemplateCMA) && InsertionCMSA(1,II+1,st,TemplateCMA)){
		fprintf(stderr,"FATAL! -> %d = '%s':",II,NameCMSA(IN_CMA[II]));
		print_error("Input alignment contains a deletion next to an insertion");
	   }
	   ins=InsertionCMSA(1,II+1,st,TemplateCMA);  // Get # insertions at site in template.
	   //******************** 3b. Insertion in template found. ***************
	   if(ins > 0){			// Insertion in template alignment.
		// Disallow insertions next to deletions. (REDUNDANT CHECK)
		assert(!(IsDeletedCMSA(1,II+1,st,TemplateCMA) 
			|| IsDeletedCMSA(1,II+1,st+1,TemplateCMA)));
		sc-=ins; 	// decrement site in cma by # inserted residues...
		cma0=ConvertColsToInsertsCMSA(cma,1,sc+1,sc+ins); NilCMSA(cma); cma=cma0;
// fprintf(stderr,"ConvertColsToInsertsCMSA(cma,%d,%d)\n",sc+1,sc+ins);
// PutCMSA(stderr,cma);	// DEBUG...
		// Convert aligned columns at sc+1 to sc+ins into an insertion.
		st--; sc--; 
	   	assert(sc >= 0);
	   //******************** 3c. One or more deletions in template found. ***************
	   } else if(IsDeletedCMSA(1,II+1,st,TemplateCMA)){
	     // Deletion in template alignment.
	     ndel=0;
	     //******************** 3ci. Count the number of deletions to right of site in input alignment. ***************
	     do {
		if(st <= 1) break;
#if 1		// Disallow insertions next to deletions. (REDUNDANT CHECK)
		if(!(InsertionCMSA(1,II+1,st,TemplateCMA) == 0 && 
				InsertionCMSA(1,II+1,st-1,TemplateCMA) == 0)){
		    fprintf(stderr,"II = %d; st = %d\n",II,st);
		    fprintf(stderr,"InsertionCMSA(st) = %d\n",
				InsertionCMSA(1,II+1,st,TemplateCMA));
		    e_type tmpE=TrueSeqCMSA(II+1,TemplateCMA);
		    PutSeqID(stderr,tmpE);
		    fprintf(stderr,"\n");
		    PutSeq(stderr,tmpE,AlphabetCMSA(TemplateCMA));
		    gsq_typ *gsq=gsqCMSA(II+1,TemplateCMA);
		    gsq->Put(stderr,AlphabetCMSA(TemplateCMA));
		    print_error("Insertions next to deletions disallowed in template");
		    fprintf(stderr,"InsertionCMSA(st-1) = %d\n",
				InsertionCMSA(1,II+1,st-1,TemplateCMA));
			assert(InsertionCMSA(1,II+1,st,TemplateCMA) == 0 && 
				InsertionCMSA(1,II+1,st-1,TemplateCMA) == 0);
		}
#endif
		ndel++; st--;
	     } while(IsDeletedCMSA(1,II+1,st,TemplateCMA));
	     // fprintf(stderr,"end of loop ...6-%d.c\n",II);

	   //******************** 3cii. Add deletions to right of site in input alignment. ***************
	     if(ndel > 0){	// REDUNDANT IF STATEMENT...
	        cma0=InsertColumnsCMSA(cma,1,sc,ndel); // add '-'s to right of sc.
		if(cma0==0) print_error("InsertColumnsCMSA() error");
		NilCMSA(cma); cma=cma0;
// fprintf(stderr,"InsertColumnsCMSA(%d,%d,cma)\n",sc,ndel);
// PutCMSA(stderr,cma);
	     } 
	   //******************** 3d. Input and template match at this position. ***************
	   } else { st--; sc--; }	// Match in template alignment.
	}
	OUT_CMA[II]=cma;
   } return OUT_CMA;
}

