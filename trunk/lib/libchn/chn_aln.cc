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

#include "chn_aln.h"

void PutMasterSlaveCMSA(FILE *fp, e_type qE, cma_typ cma)
// eventually move this to cmsa_io.cc file.
{
	// 1. First find alignment to qE
	Int4	BlkScore,i,j,t,a=11,b=2,s;
	Int4    *mtfpos=AlignSeqCMSA(stderr,&BlkScore,qE,cma);
	a_type	AB=AlphabetCMSA(cma);

#if 0	// Use ungapped for now...worry about this later...
	GapOperationsSMatrix ????
	char    *operation=GapAlnTraceSMatrix2(Int4 a, Int4 b, Int4 len, unsigned char *seq2,
        	Int4 nmod, smx_typ *M, Int4 **gapscore, Int4 *start,Int4 *score);
	Int4	len,nmod;
	unsigned char *seq2;
	smx_typ	*M=
	char    *operation=GapAlnTraceSMatrix2(a,b,len, unsigned char *seq2,
        	Int4 nmod, smx_typ *M, Int4 **gapscore, Int4 *start,Int4 *score);
        double  pernats = SitesGSS(SitesCMSA(cma))->PerNats();
        double  AdjScore = 1.0/exp((double) BlkScore/pernats);
        Int4    netleng = LenSeq(qE) - TotalLenCMSA(cma);
        if(netleng > 0) AdjScore = AdjScore*bico(netleng+nBlksCMSA(cma),nBlksCMSA(cma));
        fprintf(stderr,"Adjusted probability = %.3g\n",AdjScore);
        if(AdjScore > 0.05){   // ~ 0.05 significance w/o adjusting for length.
                print_error("Input cma file fails to find a significant match to query");
        }
#endif
#if 1
        fprintf(fp,"[0_(1)=%s(%d){go=%d,gx=%d,pn=%.1f,lf=%d,rf=%d}:\n",
                NameCMSA(cma),NumSeqsCMSA(cma)+1,30,5,5.0,0,0);
	fprintf(fp,"(%d)",LenSeq(qE));
	for(s=1; s<=LenSeq(qE); s++){ fprintf(fp,"*"); } fprintf(fp,"\n\n");
	fprintf(fp,"$%d=%d(%d):\n>",SeqI(qE),LenSeq(qE),LenSeq(qE));
	PutSeqInfo(fp,qE);
	fprintf(fp,"{()");
	for(s=1; s<=LenSeq(qE); s++){ 
		fprintf(fp,"%c",AlphaChar(ResSeq(s,qE),AB)); 
	        if(s%60==0) fprintf(fp,"\n");
	} fprintf(fp,"()}*\n\n");
        st_type S=SitesCMSA(cma);
	for(s=1; s<=NumSeqsCMSA(cma); s++){
	   fprintf(fp,"$%d=%d(%d):\n>",s+1,TotalLenCMSA(cma),LenSeq(qE));
	   e_type E = TrueSeqCMSA(s,cma);
	   PutSeqInfo(fp,E);
	   fprintf(fp,"{()");
	   for(j=1,t=1; t <= nBlksCMSA(cma); t++){
              Int4 lenM = SiteLen(t,S);
              Int4 site = mtfpos[t];
              Int4 end = site+lenM-1;
	      for( ; j < site; j++) { fprintf(fp,"-"); if(j%60==0) fprintf(fp,"\n"); }
	      unsigned char *seq=GetAlnResInSiteCMSA(t,s,cma);
	      for(i=1 ; j <= end; i++,j++) { 
		fprintf(fp,"%c",AlphaChar(seq[i],AB)); if(j%60==0) fprintf(fp,"\n");
	      }
	      
	   }
	   while(j <=LenSeq(qE)){ fprintf(fp,"-"); if(j%60==0) fprintf(fp,"\n"); }
	   fprintf(fp,"()}*\n\n");
	}
	fprintf(fp,"_0].\n");
	// exit(1);
#endif
}

e_type	PsiBlast(e_type qE,ss_type data,int argc,char *argv[],const char *USAGE,
	BooLean	silent)
{
	Int4    	T=11,width=60;
	time_t		time1=time(NULL);
        double		x_parameter=25.0;
	UInt4	hpsz=DEFAULT_PSIBLST_HPSZ,printopt=4;
	double		Ethresh=0.001,Ecutoff=0.001;
	e_type		sE=0;
	char		*checkin=0;

	if(argc < 2) print_error(USAGE);
        for(Int4 arg = 3; arg < argc; arg++){
           if(argv[arg][0] != '-') continue; // filenames... print_error(USAGE);
           switch(argv[arg][1]) {
             case 'e': Ecutoff=RealOption(argv[arg],'e',0.0,10000,USAGE); break;
	     case 'H': hpsz=IntOption(argv[arg],'H',1,1000000,USAGE); break;
             case 'h': Ethresh=RealOption(argv[arg],'h',0.0,10000,USAGE); break;
	     case 'm': printopt=IntOption(argv[arg],'m',0,6,USAGE); break;
	     case 'T': T=IntOption(argv[arg],'T',1,100,USAGE); break;
	     case 'X': x_parameter=IntOption(argv[arg],'X',1,1000,USAGE); break;
	     case 'r': checkin = AllocString(argv[arg]+2); break;
	     case ' ': break; // ignore
	     case 'C': break; // ignore
	     case 'n': 
		if(strcmp("-nochk",argv[arg]) != 0) print_error(USAGE);
		break; // ignore -nochk option intended for PsiBLAST( ).
	     case 'D': break; // ignore
	     case 'R': break; // ignore
	     case 'S': break; // ignore
	     case 'c': break; // ignore
	     case 't': break; // ignore
	     case 'u': break; // ignore
             default: print_error(USAGE); break;
           }
        }
	// gpsi_type gpsi(qE,data,Ethresh,Ecutoff,x_parameter,hpsz,1,0);
	gpsi_type gpsi(qE,data,Ethresh,Ecutoff,x_parameter,hpsz,1,checkin);
	if(silent) gpsi.KeepQuite();
	sap_typ sap = gpsi.Search(T,' ');
	if(sap) sE=SubjectSeqGSAP(sap);
	gpsi.ShowAlign(stdout,width,printopt);
	if(!silent) fprintf(stderr,"time = %0.1f seconds\n",difftime(time(NULL),time1));
	if(checkin) free(checkin);
	return sE; 
}

void    pbcSeqAlignToCMA(FILE *fp,sap_typ head,Int4 min_sq_id, Int4 max_sq_id,
        Int4 leftflank, Int4 rightflank, a_type AB)
// from: void    PutGSeqAlignList(FILE *fp, sap_typ head, Int4 width, a_type AB)
// Add leftflank & rightflank.
// [0_(1)=gold_pcna.full(201){go=18,gx=2,pn=5.0,lf=9,rf=9}:
// (46)**********************************************
// go=19,gx=2 is about 11,1 in half bits...
{
        if(head==0) return;
        Int4            i,n,ins,del,mat,min_qs,max_qe;
        UInt4   id;
        // n = NumHSPsListGSAP(head, &min_qs, &max_qe);
        n = NumHSPsListGSAP(head,min_sq_id,max_sq_id, &min_qs, &max_qe);
	Int4 len=LenSeq(head->segs->qE);
        min_qs=0; max_qe=len-1;
        fprintf(fp,"[0_(1)=Family(%d){go=19,gx=2,pn=5.0,lf=0,rf=0}:\n",n);
        fprintf(fp,"(%d)",len);
        for(i=0; i < len; i++) fprintf(fp,"*");
        fprintf(fp,"\n"); i=0;
        for(sap_typ gsap=head; gsap!=NULL; gsap=gsap->next){
          id = SeqI(SubjectSeqGSAP(gsap));
          if(id < min_sq_id || id > max_sq_id) continue;
          dsp_typ dsp=gsap->segs; i++;
          // fprintf(fp,"\t\tlength = %d\n",LenSeq(dsp->sE));
          xPutSeqAlignToCMA(fp,gsap,leftflank,rightflank,min_qs,max_qe,id,AB);
        } fprintf(fp,"\n");
        fprintf(fp,"_0].\n");
}

void    tpbcSeqAlignToCMA(FILE *fp, sap_typ head, Int4 min_sq_id, Int4 max_sq_id,
        Int4 leftflank, Int4 rightflank, a_type AB)
// from: void    PutGSeqAlignList(FILE *fp, sap_typ head, Int4 width, a_type AB)
// Add leftflank & rightflank.
// [0_(1)=gold_pcna.full(201){go=18,gx=2,pn=5.0,lf=9,rf=9}:
// (46)**********************************************
// go=19,gx=2 is about 11,1 in half bits...
{
        if(head==0) return;
        UInt4   last_id=0,id;
        Int4            i,n,ins,del,mat,min_qs,max_qe;

        n = NumSeqsListGSAP(head, min_sq_id,max_sq_id, &min_qs, &max_qe);
        fprintf(fp,"[0_(1)=Family(%d){go=19,gx=2,pn=5.0,lf=0,rf=0}:\n",n);
	Int4 len=LenSeq(head->segs->qE);
        min_qs=0; max_qe=len-1;
        fprintf(fp,"(%d)",len);
        for(i=0; i < len; i++) fprintf(fp,"*");
        fprintf(fp,"\n"); i=0;
        for(sap_typ gsap=head; gsap!=NULL; gsap=gsap->next){
          id = SeqI(SubjectSeqGSAP(gsap));
          if(id < min_sq_id || id > max_sq_id) continue;
          dsp_typ dsp=gsap->segs;
          if(dsp->subject_id != last_id){
                last_id= dsp->subject_id; i++;
                mat = InDelsSeqAlign(gsap,&ins,&del,min_qs,max_qe);
                Int4 len = LenSeq(dsp->sE);
                fprintf(fp,"\n$%d=%d(%d):\n",i,len,mat+del);
                // fprintf(fp,">"); PutSeqInfo(fp,dsp->sE);
                // fprintf(fp,"\t\tlength = %d\n",LenSeq(dsp->sE));
                tPutSeqAlignToCMA(fp,gsap,leftflank,rightflank,min_qs,max_qe,AB);
          } fprintf(fp,"\n");
        } fprintf(fp,"\n_0].\n");
}


#if 0	// TODO:
	1. PSI-BLAST2CMA( ) - fix start positions when including more than the 
		highest scoring gapped HSP as another locally aligned subsequence...
	2. Option to save only those sequences that align over x% of the
		full length of the query sequence.  This avoids poorly
		related sequences from being included in the alignment.
	3. Dealing with gaps and 'X' characters in dsc?
	4. Remove sequences closely related to homologs...
	5. Statistics for clusters of conserved residues.
	     Add pdb_typ into pbc_typ.
	6. Statistical highlighting for input of a single sequence against 
		a data file.
	7. Convert probe cma files into PsiBC format (1 block) to allow output 
		in these routines.
	10. Remove sequences with < x% alignment to orthologs from main cma!!!
		Or option in gblastpgp??
	11. Print out LinearToLog scores below histogram.
	13. Get orthologs program:
		seed seq -> (fooblast) -> ids -> (getrefs.pl) -> refs. ->
		   (RefsToSeqs) -> seqs in tax groups -> (gblastpgp) ->
		   best hits in each group -> ?? consistency check??
	14. Program to select meaningful paralog set (key to analysis).
	15. Allow option to create alignments for two interacting subunits
	    at one time (used only with rasmol and pdb).
	17. Add option to print out superfamily alignment only (limit to
	     100 hits). E.g., if number of sequences == 1.
	18. Gapxdrop parameter set higher?? 100?? to extend into other
		regions after finding match??
IMPORTANT!!:
	21. Fix situation where not all of query gets aligned (truncated cma).
	23. Show all subfamilies if small enough.
#endif
const char USAGE_START[]="USAGE: chn_aln query_family database [options]\n\
USAGE(telescoping): chn_aln query_family <family> <superfamily> <supersuperfamily> ... [options]\n\
   options:\n\
        -C<file> write checkpoint file\n\
        -c<file> write sequences detected to file\n\
        -D       Do SEG query sequence\n\
        -d       don't remove sub- or identical sequences in a subset from a superset\n\
        -e<real> E-value cutoff for showing scores and alignments\n\
        -h<real> E-value cutoff for inclusion in next profile\n\
        -H<int>  heapsize (default 20000)\n\
        -j<int>  maximum number iterations (default: infinite)\n\
                  Input value must be > 1.\n\
        -L       do low complexity masking of database sequences\n\
        -M<int>  minfix score (default 0)\n\
        -M=<file.cma> Use the specified motif cma file as prior for search\n\
        -m<int>  print option m=0..6 (default 4)\n\
        -nochk   don't produce default checkpoint (.chk) files.\n\
        -o       use only the single best HSP\n\
        -R<file> Read checkpoint file\n\
        -r<file> Read checkpoint file for query family (?)\n\
        -s       remove subsequences within a superset that are present in a subset\n\
        -T<int>  blast word hit threshold (default 11)\n\
        -t<real> trim alignments at ends where marg. prob. is > <real>\n\
        -u       use the first sequence as the consensus query sequence\n\
        -v       turn off verbose output\n\
        -X<int>  X dropoff for gapxdrop\n\
\n";

/**************************** Global Variables ******************************/
int	ChainAlign(Int4 argc,char *argv[])
{ 
	Int4	i,arg,minfix=0,iterations=0,arg_start;
	cma_typ cma,mtfcma=0;
	char	*main_cma=0;
	FILE	*fp;
	BooLean	silent=FALSE,UseAllHSPs=TRUE;
	BooLean	first_as_consensus=FALSE;
	char	rm_seq='I',*filename[50];
	Int4	num_sets;
	BooLean	r_option=FALSE;
	BooLean	use_query_checkpoint=FALSE;
	Int4	iter_ticks=0;

	if(argc < 3) print_error(USAGE_START);
	time_t	time1=time(NULL); 
#if 0
        fp = open_file(argv[1],".arg","w");
        for(i = 0; i < argc; i++) { if(argv[i][1] != ' ') fprintf(fp,"%s ",argv[i]); }
        fprintf(fp,"\n"); fclose(fp);
#endif
	filename[0]=argv[1];  // subfamily of interest.
	filename[1]=argv[2];  // first superset.
	num_sets=1;
	for(arg=arg_start=3; arg < argc; arg++){
		if(argv[arg][0] != '-'){
		   if(arg_start >= 50) { print_error("Too many input sets"); }
		   num_sets++;
		   filename[num_sets] = argv[arg]; 
		   arg_start++; 
		} else break;
	}
	a_type	  A=MkAlpha(AMINO_ACIDS,GBLAST_BLOSUM62); // THIS BLOSUM IMPORTANT!!
	for(arg = arg_start; arg < argc; arg++) {
	   if(argv[arg][0] != '-') print_error(USAGE_START);
	   switch(argv[arg][1]) {
	     case 'M': 
		if(argv[arg][2] == '='){ // then this is a BG_CMA file
		    if(argv[arg][3]==0) print_error(USAGE_START);
		    mtfcma = ReadCMSA2(argv[arg]+3,A); 
		    ReNameCMSA("Motifs",mtfcma);
		    argv[arg][1]=' '; 
	   	} else {
			minfix=IntOption(argv[arg],'M',1,10000,USAGE_START); 
			argv[arg][1]=' '; 
		} break;
	     case 'd': rm_seq=0; argv[arg][1]=' '; break;
	     case 'j': iterations=IntOption(argv[arg],'j',1,10000,USAGE_START); 
			argv[arg][1]=' '; break;
	     case 'u': first_as_consensus=TRUE; break;
	     case 'v': silent=TRUE; argv[arg][1]=' '; break;
	     case 'o': UseAllHSPs=FALSE; argv[arg][1]=' '; break;
	     case 's': rm_seq='S'; argv[arg][1]=' '; break;
	     case 'r': r_option=TRUE; break;
	     case 'R': use_query_checkpoint=TRUE; break;
	     case 'i':
		if(strcmp("-iter_ticks",argv[arg]) == 0){
			iter_ticks=arg; 
			argv[arg][1]=' ';  // hide for now but turn on later.
		}
		break;
	     // -nochk == pass on to PsiBLAST( );
	     // -iter_ticks == pass on to PsiBLAST( );
	     default: ; // print_error(USAGE_START);
	   }
	}
#if 0
	// if two families (A.txs and B.txs) and a database (DBS).
	// then do family to family comparison.
	// 1. Align A.txs concensus against DBS saving a checkpoint file.
	// 2. Align B.txs concensus against previous checkpoint file.
	// 3. Compare A versus DBS alignment with B versus DBS alignment.
#endif
// fprintf(stderr,"DEBUG 1\n");

	// 1. gblastpgp keyfa against orthogs to get *.oma file.
	Int4	query_start,*Intervals=0;
	
	ss_type subfam_data=0,data=0;
	cma_typ	*tcma=0;
	e_type	qE,keyE,*ListE=0;
	Int4	g,G;
	txs_typ *txs = new txs_typ(argv[1],A);	// taxonomic sequence type...
	if(!txs) print_error("The first file must be generated by the chn_txs program");
	if(txs->NumTaxGroups( ) > 0){  // then create representative set.
// fprintf(stderr,"DEBUG 1\n");
		qE=txs->QuerySeq();
		for(Int4 r=1; r <= LenSeq(qE); r++){
		  if(!ResSeq(r,qE)) print_error("query input error: 'X' residues!");
		}
		NEW(ListE,txs->NumTaxGroups( )+3,e_type);
		for(G=0,g=1; g <= txs->NumTaxGroups( ); g++){
		   if(r_option && g==1){
		   	subfam_data=SeqSet1("query sequence",qE,A); // makes a copy of qE!
		   } else subfam_data=txs->Group(g);
		   keyE=PsiBlast(qE,subfam_data,argc,argv,USAGE_START,silent); // get best hit.
		   if(keyE){ G++; ListE[G]=CopySeq(keyE); EqSeqI(G,ListE[G]); }
		   else {
			PutSeqSetPIDs(stderr, subfam_data);
			print_error("input txs_typ error");
		   }
		   if(r_option && g==1){ NilSeqSet(subfam_data); }
		} subfam_data=Array2SeqSet(ListE,G,argv[1],A);
	} else { 
	   delete txs; txs=0;
	   subfam_data=MakeSeqSet(argv[1],200000,A); 
	}
// fprintf(stderr,"DEBUG 1\n");
	if(NSeqsSeqSet(subfam_data) < 1) print_error("input error.");

	// 2. Get Consensus sequence from the ortholog alignment.
	fprintf(stderr," Aligning query sequence set.\n");
	if(use_query_checkpoint){
	  // warning: uses -Rcheckin file not -rcheckin!
	  cma=PsiBLAST(&query_start,qE,subfam_data, argc, argv, USAGE_START,1,A,FALSE,0,FALSE,
			(Int4)INT4_MAX,(cma_typ)0,(BooLean)TRUE);
	}
	else cma=PsiBLAST(&query_start,0,subfam_data, argc, argv, USAGE_START,2,A,FALSE,0,FALSE);
	// WARNING: reorders tax groups by significance to query!
	if(cma==0) print_error("possible memory failure");
	// fprintf(stderr,"query_start = %d\n",query_start);
	if(first_as_consensus) qE=TrueSeqCMSA(1,cma);
	else qE=MkConsensusCMSA(cma);
	if(!silent){
	   fp=open_file(argv[1],".csq","w"); PutSeq(fp,qE,A); fclose(fp);
	}
	// PutSeq(stdout,qE,A); exit(1);

// fprintf(stderr,"DEBUG 1\n");
	if(txs){  // then create additional cma files.
		assert(NumSeqsCMSA(cma) == txs->NumTaxGroups( ));
		NEW(tcma,txs->NumTaxGroups( )+3,cma_typ);
		for(g=1; g <= txs->NumTaxGroups( ); g++){
		   keyE=TrueSeqCMSA(g,cma);
		   Int4 g0 = SeqI(keyE);
		   assert(g0 <= txs->NumTaxGroups( ) && g0 > 0);
		   data=txs->Group(g0);
		   tcma[g]=PsiBLAST(&query_start,keyE,data,argc,argv,USAGE_START,2,A,
				   FALSE,0,FALSE);
		   if(tcma[g]){
#if 1	// Order tcma significance of hit...
			char str[200];
			sprintf(str,"%s{%c}",txs->GroupName(g0),
				txs->Kingdom(g0));
			ReNameCMSA(str,tcma[g]); 
			// Need to get the right txs file!!!
#else
			ReNameCMSA(txs->GroupName(g),tcma[g]); 
#endif
		   } 
		}
	}

// fprintf(stderr,"DEBUG 1\n");
	// 3. gblastpgp consensus against superfamily to get *.pma file.
#if 0
	// if query absent, then add consensus sequence to set.
	// This eliminates NULL pointers in the WtFreq array...
	// This may be a problem when searching a very large database...
	// Ultimately fix this program to allow shorter alignments...
	if(!IsInSeqSet(keyE,data)) AddSeqSet(qE, data); 
#endif
	cma_typ mcma;
   if(num_sets == 1){
	// fprintf(stderr,"get superfamily alignment\n");
	data=MakeSeqSet(argv[2],200000,A);
#if 1	// if number of blocks in mtfcma is > 1 then...
     if(mtfcma && nBlksCMSA(mtfcma) > 1){
	// 1. print out master-slave cma file.
        fp = open_file(argv[1],"_ms.cma","w");
	PutMasterSlaveCMSA(fp,qE,mtfcma);
	// PutMasterSlaveCMSA(stdout,qE,mtfcma);
	fclose(fp);
	// exit(1);
	// 2. replace mtfcma with master-slave cma.
	TotalNilCMSA(mtfcma);
        fp = open_file(argv[1],"_ms.cma","r");
	mtfcma = ReadCMSA(fp,A);
	fclose(fp);
     }
#endif
	fprintf(stderr," Aligning query against main sequence set....\n");
	if(iter_ticks) argv[iter_ticks][1]='i';  // turn on for main set.
	mcma=PsiBLAST(&query_start,qE,data,argc,argv,USAGE_START,iterations,A,
			TRUE,minfix,UseAllHSPs,mtfcma);
	if(iter_ticks) argv[iter_ticks][1]=' ';  // turn off again.
	NilSeqSet(data);
	if(mcma==0) { fprintf(stdout,"No hits!\n"); return 0; }
	fp = open_file(argv[1],".chn","w");
	// fprintf(stderr,"Printing orthologs\n");
	ReNameCMSA("QuerySubFamily",cma);
	PutCMSA(fp,cma); 
	// fprintf(stderr,"Printing database hits...");
	ReNameCMSA("SuperFamily",mcma);
	PutCMSA(fp,mcma); 
	// fprintf(stderr,"done.\n");
	if(txs){
		for(g=1; g <= txs->NumTaxGroups( ); g++){
		   if(tcma[g]){ PutCMSA(fp,tcma[g]); TotalNilCMSA(tcma[g]); }
		} free(tcma); delete txs;
	} fclose(fp); 
	TotalNilCMSA(mcma); 
    } else {	//**** create telescoping sets based on multiple input files...

   ss_type	Data[50];
   Int4		MaxSqID[50];
   e_type	superE,subE;
   BooLean	found;
   Int4		N,n,m,set,subset,superset;
   UInt4	sq;

   // 1. Open up database files.
   for(n=N=txs->TotalSeqs(),set=1; set <= num_sets; set++){
	fprintf(stderr,"reading file '%s'\n",filename[set]);
	Data[set] = SeqSet(filename[set],A);
	if(n >= NSeqsSeqSet(Data[set])){
	   fprintf(stderr,"subset = %d; superset = %d\n",n,NSeqsSeqSet(Data[set]));
	   print_error("'subset' >= its 'superset'");
	}
	n = NSeqsSeqSet(Data[set]); N += n;
   }

   // 2. Create subfamily sequence set.
   NEW(ListE,N+3,e_type);
   sq=0;
   // sq=1; ListE[sq]=CopySeq(qE);
   for(g=1; g <= txs->NumTaxGroups( ); g++){
	data=txs->Group(g);
	for(n=1; n <= NSeqsSeqSet(data); n++){
              sq++; ListE[sq] = SeqSetE(n,data); 
	}
   } MaxSqID[0]=sq;

   // 3. Construct hierarchical based sequence numbering and database file...
   // Remove redundant sequences...
   for(superset=1; superset <= num_sets; superset++){
	subset = superset-1;
	fprintf(stderr,"removing %s sequences (%d) from '%s' (%d)\n",
				filename[subset],sq,filename[superset],
				NSeqsSeqSet(Data[superset]));
        for(n=1; n <= NSeqsSeqSet(Data[superset]); n++){
              superE = SeqSetE(n,Data[superset]);
	      for(found=FALSE,m=1; m <= MaxSqID[subset]; m++){
		  subE = ListE[m];
                  if(rm_seq=='S'){ if(IsSubSeq(superE,subE)){ found=TRUE; break; } }
                  else if(rm_seq=='I'){
		    if(FastIdentSeqs(superE,subE)){ found=TRUE; break; }
	          } else break;	// Leave all of the sequences in subsets
	      } if(!found){ sq++; ListE[sq] = superE; }
        } MaxSqID[superset] = sq;
	fprintf(stderr,"  %d sequences in full set\n",sq);
   } // data=Array2SeqSet(ListE,sq,argv[2],A);
   fp=tmpfile();
   for(UInt4 s=1; s <= sq; s++){ PutSeq(fp,ListE[s],A); }
   rewind(fp); data=SeqSet_fptr(fp,A); fclose(fp);

   // 4. Run PsiBC on this data file.
	// This call sorts the sap_typ within RtnSAP_PsiBLAST based on sequence IDs.
	sap_typ sap=RtnSAP_PsiBLAST(qE,data,argc,argv,USAGE_START,iterations,A,minfix,FALSE,0,FALSE);
	if(sap==0) { fprintf(stdout,"No hits!\n"); return 0; }
	else {
	  fp = open_file(argv[1],".chn","w");
	  fprintf(stderr,"Printing orthologs\n");
	  ReNameCMSA("QuerySubFamily",cma);
	  PutCMSA(fp,cma); 
	  fprintf(stderr,"Printing database hits...");

	  Int4 left_flank=0,right_flank=0;
	  query_start=QueryStartGSeqAlnList(sap);
	  // always include query in database search...
	  Int4 *starts,*lens;
	  NEW(starts, 5,Int4); NEW(lens,5,Int4);
	  starts[0]=starts[1]=0; lens[0]=LenSeq(qE);
	  // if(SortBySeqID) sap = SortBySeqIDGSAP(sap); // don't need to sort these...
	  sap_typ head=MakeGSeqAlign(1,SeqI(qE),SeqI(qE),qE,qE,starts,lens); 
	  head->next = sap;
	  for(set=1; set <= num_sets; set++){
	      if(UseAllHSPs){
	         // save all alignments as the length of the query regardless...
	         // this will result in marginal probability errors though...
		 pbcSeqAlignToCMA(fp,head,MaxSqID[set-1],MaxSqID[set],
			left_flank,right_flank,A);
	      } else {
		 tpbcSeqAlignToCMA(fp,head,MaxSqID[set-1],MaxSqID[set],
			left_flank,right_flank,A);
	      }
	      if(set==1 && txs){
		for(g=1; g<=txs->NumTaxGroups( ); g++){ if(tcma[g]) PutCMSA(fp,tcma[g]); }
	      } 
	  } head->next = NULL; GSeqAlignFree(head);
	  fclose(fp); 
	  FreeGSeqAlignList(sap);
	}
   // 5. Free up memory.
	for(set=1; set <= num_sets; set++) NilSeqSet(Data[set]); 
	// ListE=NilSeqSetRtnSeqs(data); 
	free(ListE); ListE=0;
	NilSeqSet(data);
	NilSeqSet(subfam_data);
	if(txs){
	   for(g=1; g <= txs->NumTaxGroups( ); g++){ if(tcma[g]) TotalNilCMSA(tcma[g]);}
	   free(tcma); delete txs;
	}
    } //********************** end telescoping sets...
	TotalNilCMSA(cma); 
#if 0	// need to see how to deal with subfam_data, which takes in ListE seqs!!!
	if(ListE){
	   for(i=1; ListE[i]; i++) NilSeq(ListE[i]);
	}
#endif
	if(!first_as_consensus) NilSeq(qE); 
	NilAlpha(A); 
	if(!silent){
		double runtime=difftime(time(NULL),time1);
		fprintf(stderr,"\ttime: %0.1f seconds (%0.2f minutes)\n",runtime,runtime/60.0);
	}
	return 0;
}


