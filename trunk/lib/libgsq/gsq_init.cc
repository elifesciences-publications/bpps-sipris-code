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

void	gsq_typ::init( )
{
	fakeE=NULL; indels=NULL; f2r=NULL; 
	num_open=num_insert=num_del=0;
}

Int4    gsq_typ::read(FILE *fp,Int4 *pos,a_type A)
{
        Int4    len_blk[10],nblks,id,real_len,fake_len;
        char    c;
	fpos_t  ptr;
        assert(fgetpos(fp,&ptr) == 0);

        while((c=fgetc(fp)) != '$');
        if(c==EOF) print_error("gsq_typ::read() input error");
        else ungetc(c,fp);
        if(fscanf(fp,"$%d=%d(%d):",&id,&real_len,&fake_len) != 3){
                print_error("gsq_typ::read() input error");
        } len_blk[1]=fake_len; nblks=1;
	assert(fsetpos(fp,&ptr)==0);
        read(fp,nblks,len_blk,pos,A);
}

Int4	gsq_typ::read(FILE *fp,Int4 nblks,Int4 *len_blk,Int4 *pos,a_type A)
/**********************************************************************
Read in gsq object from a file.
-------------------- File format ------------------------
$1=151(145):
>gi|2984277 {|13(50)|}(AE000770) 2-acylglycerophosphoethanolamine acyltransferas
PLWL{(GVRPF)F-RALFR()IEVIGKENIP(E)GAcIVASNHRSHLDPPVL(NAVFPEPL)VFLAKEE(L)FKPPF-GGILKHMRAIPLRR(GSEDISTL)EECVSLL()KLGCKIGIFPEGTRANPGE(FK)RAKPGVGFLAINSG()FPVLPVYIDGTDRAFP(RGK)}GTE*

state:
  'B' = begin.
  'L' = Definition line
  'N' = N-terminal extension.
  'C' = C-terminal extension.
  'D' = Deletion.
  'M' = Match.
  'I' = Insertion.
  'G' = gap.
 **********************************************************************/
{
	Int4		N,o,i,j,J0,t,r,f,b,str_len;
	Int4		real_len,fake_len,I,offsetX=0,extendX=0;
        char   		c,*realsq,state,id[MAX_SEQ_DEFLINE_LENG+100];
	char		info[MAX_SEQ_DEFLINE_LENG +100],descript[MAX_SEQ_DEFLINE_LENG+100];
	unsigned char	*fakesq;
	rff_typ		*rff=0;

	Free();
	num_insert=num_del=num_open=0;
	for(t=0,state='B'; (c=getc(fp)) != '*'; ){
	    if(isspace(c)) continue;
            switch(c){
              case '$': 	// SeqId and lengths...
		assert(state == 'B'); state = 'L';
		assert(fscanf(fp,"%d=%d(%d):",&I,&real_len,&fake_len)==3);
		realsq= new char [real_len+3]; r=0;
		fakesq= new unsigned char [fake_len+3]; f=0;
		indels = new indel_typ[fake_len+3]; indels[0]=0;
		f2r = new f2r_typ[fake_len+3];  f2r[0]=0;
		break;
	      case '>': 	// fasta defline. Must be less than MAX_SEQ_DEFLINE_LENG!!!
		assert(state=='L'); state='N'; 
		assert(fgets(info,MAX_SEQ_DEFLINE_LENG+2,fp)!=NULL); 
		str_len = strlen(info);
		if(str_len==(MAX_SEQ_DEFLINE_LENG+1)) {	// maximum length allowed...
		  if((c=info[str_len-1]) != '\n'){
			while(c != '\n' && c != EOF) c=getc(fp);
			if(c == EOF) print_error("gsq_typ::read( ): EOF encountered!");
		  }
		}
#if 0	// This can't happen...
		if(strlen(info) > MAX_SEQ_DEFLINE_LENG+1){
			fprintf(stderr,"strlen(info) = %d\ninfo=%s\n",strlen(info),info);
		  assert(strlen(info) <= MAX_SEQ_DEFLINE_LENG); // make sure line is not too long.
		}
#endif
		break;
#if 1	// read rff information...
	      case '+': 	// // read rff information...
                assert(rff == 0); ungetc(c,fp); // put back '+';
                rff = new rff_typ(fp); assert(state=='N'); state='N'; 
		break;
#endif
	      case '{': assert(state=='N'); offsetX=r; indels[0]=r; state='n'; break;
	      case '}': extendX=0; state = 'C'; break;    // C-terminal extension
	      case ')': t++;
		if(t > nblks) state='e'; else { state='M'; b=0; pos[t]=f+1; }
		break;
	      case '(': state = 'G'; break; // assert(t < 1 || b==len_blk[t]); break; 
	      case '-':
		if(state !='D'){ num_open++; 
// if(state=='I'){	} // for "...MATg-YQIL..." convert to "...MATGYQIL..."
// assert(state=='M'); 
#if 0
if(state!='M'){		// for "...MATg-YQIL..." convert to "...MATGYQIL..." ???
	fprintf(stderr,"state = '%c'; id = %d\n",state,I);
	print_error("gsq_typ read error");
}
#endif
		}
		f++; state='D'; fakesq[f]=AlphaCode('X',A); f2r[f]=r;
		indels[f]=gseq_del_mask; b++; num_del++; break;
              default: 
		if(!isalpha(c)){
		    fprintf(stderr,"seq %d: f = %d; r = %d; c = '%c'\n",I,f,r,c);
		    fprintf(stderr,"info=%s\n",info);
		    fprintf(stderr,"strlen=%d\n",str_len);
		    while((c=getc(fp)) != '\n' && c != EOF) fprintf(stderr,"%c",c);
		    fprintf(stderr,"\n");
		    if(c == EOF){
			print_error("gsq_typ::read( ): EOF encountered!");
		    }
		    assert(!"gsq_typ::read( ): input error 1"); 
		}
		switch(state){
		  case 'N': realsq[r++]=c; break;	// N-terminal extension.
		  case 'M': case 'D': case 'I':		// Match, Deletion, Insertion.
	            if(isupper(c)){
			realsq[r++]=c; f++; fakesq[f]=AlphaCode(c,A); 
			indels[f]=0; f2r[f]=r; state='M'; b++;
	            } else if(islower(c)){	// insertions...
			if(state != 'I') num_open++; 
			num_insert++; state = 'I'; realsq[r++]=c; indels[f]++;
			assert(insert(f) < gseq_ins_mask);  
			// avoid overflow into delete bit...
	            } break;
		  case 'G': realsq[r++]=c; f++; indels[f]=0; 	// Gap between blocks?
		    fakesq[f]=AlphaCode(c,A); f2r[f]=r; break;
		  case 'C': extendX++; realsq[r++]=c; break;	// C-terminal extension.
                  default: 
		    fprintf(stderr,"seq %d: f = %d; r = %d; state = '%c'\n",I,f,r,state);
		    print_error("gsq_typ::read( ): input error 2"); 
                }
	   }
	   if(f > fake_len || r > real_len){
		   fprintf(stderr,"%d: f = %d (%d); r = %d (%d)\n",
			I,f,fake_len,r,real_len);
		   std::cerr << c; std::cerr << std::endl; realsq[r]=0;
		   std::cerr << realsq; std::cerr << std::endl;
		assert(f <= fake_len && r <= real_len);
		print_error("gsq_typ::read( ): input error");
	   }
	} 
	if(!(f==fake_len && r==real_len)){
	    fprintf(stderr,"%d: f = %d (%d); r = %d (%d)\n",
			I,f,fake_len,r,real_len);
	    realsq[r]=0;
	    fprintf(stderr,"seq=%s\n",realsq);
	    print_error("f==fake_len && r==real_len failed!"); 
	} realsq[r]=0;
	
	indels[f] = (gseq_del_mask & indels[f]) | (gseq_ins_mask & extendX); 
	realE=StringToSeq(realsq,info,I,A);
	if(rff) AddRFF2Seq(realE,rff);
// std::cerr << realsq; PutSeq(stdout,realE,A);
	offsetX += OffSetSeq(realE); extendX += CtermExtendSeq(realE);
	StrSeqID(id, DEFAULT_INFO_LINE_LEN_SEQ, realE); 
	StrSeqDescript(descript,DEFAULT_INFO_LINE_LEN_SEQ,realE);
#if 0	// Don't use this!!!
fprintf(stderr,"phylum='%s'; kingdom='%c'\n",PhylumSeq(realE),kingdomSeq(realE));
	fakeE=MakeSeq(id,descript,offsetX,extendX,fake_len,fakesq,
					PhylumSeq(realE),KingdomSeq(realE));
#else
	// if(islower(kingdomSeq(realE)))fprintf(stderr,"phylum='%s'; kingdom='%c'\n",
	//	PhylumSeq(realE),kingdomSeq(realE));
	fakeE=MakeSeq(id,descript,offsetX,extendX,fake_len,fakesq);
#endif
	EqSeqI(I,fakeE);
	delete []realsq; delete []fakesq;
// if(I==136) { Put(stderr, 60, A); Put(stderr, A); }
	return 0;
}

void	gsq_typ::initialize(Int4 leftflank, Int4 rightflank, char *operation, 
	Int4 trace_length, Int4 start0, e_type E)
{ initialize(leftflank,rightflank,operation,trace_length,start0,E,0); }

void	gsq_typ::initialize(Int4 leftflank, Int4 rightflank, char *operation, 
	Int4 trace_length, Int4 start0, e_type E,Int4 *pos)
{
	e_type	oldE;
	if(fakeE != NULL) Free();
	realE = E;
	oldE=initialize(leftflank, rightflank, operation,trace_length,start0,pos);
	if(oldE != NULL) NilSeq(oldE);
}

e_type	gsq_typ::initialize(Int4 leftflank, Int4 rightflank, char *operation, 
	Int4 trace_length, Int4 start0)
{ return initialize(leftflank, rightflank, operation,trace_length,start0,(Int4*)0); }

// #include "residues.h"
e_type	gsq_typ::initialize(Int4 leftflank, Int4 rightflank, char *operation, 
	Int4 trace_length, Int4 start0, Int4 *pos)
/****************************************************************************
  create 'fake' subsequence with gaps and flanking 
  trace = "EDdmmmmmmmmmmmiiiMmmmmmmmmmmmmmmmmmmmmiiMmmmmmmmmmmE"
 ****************************************************************************/
{
	e_type		oldE;
	Int4		N,o,j,k,t,J0,flankL,flankR;
        unsigned char   *rseq,*gapseq;
	char		state,id[MAX_SEQ_DEFLINE_LENG + 100];
	char		descript[MAX_SEQ_DEFLINE_LENG + 100];
	Int4		offsetX=0,extendX=0,start;

	assert(leftflank >= 0 && rightflank >= 0); assert(realE != NULL);
	assert(operation[1]!='i'); 
	assert(operation[0]=='E'); 
	assert(operation[trace_length-2]!='i');
	if(operation[trace_length-1] != 'E'){
		fprintf(stderr,"operation = %s\n",operation);
		assert(operation[trace_length-1]=='E');
	}
	oldE=fakeE; fakeE=NULL; Free(); 
	flankL=MINIMUM(Int4,leftflank,start0-1); 
	start=start0-flankL-1;	// start relative to real subsequence.
	offsetX=OffSetSeq(realE)+start; // start relative to real full sequence.

	N = flankL + trace_length + rightflank;
	gapseq = new unsigned char[N+2];
	f2r = new f2r_typ[N+3];
	indels = new indel_typ[N+3];
	num_open=num_del=num_insert=0; indels[0]=start; f2r[0]=0; 

	// 1. Add leftflanking sequence.
	for(rseq=SeqPtr(realE),J0=0,j=start,k=1; k<= flankL; k++){
	    J0++; j++; gapseq[J0]=rseq[j]; indels[J0]=0; f2r[J0]=j;
	}
	// 2. Add aligned-region sequence.
        for(t=0,o=1,state='E'; operation[o] != 'E'; o++){
            switch(operation[o]){
               case 'M': if(pos){ t++; pos[t] = J0+1; }
	       case 'm': J0++; j++; gapseq[J0] = rseq[j];
			indels[J0] = 0; f2r[J0] = j; break;
               case 'D': if(pos){ t++; pos[t] = J0+1; }
	       case 'd': // deletion in sequence relative to profile.
                        J0++; gapseq[J0] = 0; num_del++;
			if(state != 'D' && state != 'd') num_open++;
			indels[J0] = gseq_del_mask; f2r[J0] = j; break;
               case 'i': // insert is between profile blocks;
                        J0++; j++; gapseq[J0] = rseq[j]; 
			indels[J0] = 0; f2r[J0] = j; break;
               case 'I': // Insert within a profile block; delete from seq.
		      {
			indel_typ i;
			assert(state != 'I'); num_open++; 
			for(i=0; operation[o]=='I'; ){ j++; i++; num_insert++; o++; }
			o--; assert(i < gseq_ins_mask); 
			if(state != 'D' && state != 'd') indels[J0] = i; 
			else indels[J0] += i; // keep delete bit.
		      } break;
               default:
                        fprintf(stderr,"operations = %s\n",operation);
                        print_error("gsq_typ::initialize( ): input error"); break;
            }  state=operation[o];
        }
	// 3. Add rightflank region.
	flankR = MINIMUM(Int4,rightflank,LenSeq(realE)-j);
	for(Int4 i=1; i<= flankR; i++){
	    J0++; j++; gapseq[J0]=rseq[j]; indels[J0]=0; f2r[J0]=j;
	}
	Int4 x=(LenSeq(realE)-j);
#if 1	// Bug Fix!!	Fake seq extends beyond real seq so x is negative.
	// gseq_ins_mask & negative == a large number!!!!
	while(x < 0){ J0--; x++; }
#endif
	indels[J0] = deletion(J0) | (gseq_ins_mask & (indel_typ) x); 
	extendX = CtermExtendSeq(realE) + x;
	StrSeqID(id, DEFAULT_INFO_LINE_LEN_SEQ, realE); 
	StrSeqDescript(descript,DEFAULT_INFO_LINE_LEN_SEQ, realE);
	fakeE=MakeSeq(id,descript,offsetX,extendX,J0,gapseq);
#if 0
	// NEED TO FIX THIS AS SEQUENCES WITH N-terminal flank are 
	// not working properly...[e.g., {(ADE)...]
	if(FakeToReal(LenSeq(fakeE)) > LenSeq(realE)){
		PutSeqInfo(stderr,fakeE);
		fprintf(stderr,"fake len = %d\n",LenSeq(fakeE));
		PutSeqInfo(stderr,realE);
		fprintf(stderr,"real len = %d\n",LenSeq(realE));
		assert(FakeToReal(LenSeq(fakeE)) <= LenSeq(realE));
	}
#endif
// assert(!(operation[1] == 'D' && operation[2] == 'I'));
// if(flag){ printf("%s\n",operation); Put(stderr,A); fprintf(stderr,"%s\n",operation); NilAlpha(A);}
	delete []gapseq;
	return oldE;
}

void    gsq_typ::initialize(Int4 nBlks,Int4 *len,Int4 *pos,e_type E)
{
        char    *operation=FromSitesOperation(nBlks,len,pos,E);
        this->initialize(operation,E,pos); // copies new sites to sites...
        // fprintf(stderr," = %s\n",operation);
        free(operation);
}

void	gsq_typ::initialize(char *operation, e_type E, Int4 *pos)
/****************************************************************************
  create 'fake' subsequence with gaps; no flanking regions.
  trace = "EiiiiDdmmmmmmmmmmmiiiMmmmmmmmmmmmImmmmmmmmmiiMmmmmmmmmmmiiiiE"
 ****************************************************************************/
{
	Int4		N,o,j,k,t,J0;
        unsigned char   *rseq,*gapseq;
	char		state,id[MAX_SEQ_DEFLINE_LENG + 100];
	char		descript[MAX_SEQ_DEFLINE_LENG + 100];
	assert(operation[0]=='E'); 
	Int4		offsetX=0,extendX=0,trace_length=strlen(operation);
	assert(operation[trace_length-1]=='E');
	assert(realE == NULL); assert(fakeE == NULL);
	realE=E; 
	offsetX=OffSetSeq(realE); // start relative to real full sequence.

	N = trace_length;
	gapseq = new unsigned char[N+2];
	f2r = new f2r_typ[N+3];
	indels = new indel_typ[N+3];
	num_open=num_del=num_insert=0; indels[0]=0; f2r[0]=0; 

	// 2. Add aligned-region sequence.
        for(rseq=SeqPtr(realE),t=0,J0=0,j=0,o=1,state='E'; operation[o] != 'E'; o++){
            switch(operation[o]){
               case 'M': if(pos){ t++; pos[t] = J0+1; }
	       case 'm': J0++; j++; gapseq[J0] = rseq[j];
			indels[J0] = 0; f2r[J0] = j; break;
               case 'D': if(pos){ t++; pos[t] = J0+1; }
	       case 'd': // deletion in sequence relative to profile.
                        J0++; gapseq[J0] = 0; num_del++;
			if(state != 'D' && state != 'd') num_open++;
			indels[J0] = gseq_del_mask; f2r[J0] = j; break;
               case 'i': // insert is before, between or after profile blocks;
                        J0++; j++; gapseq[J0] = rseq[j]; 
			indels[J0] = 0; f2r[J0] = j; break;
               case 'I': // Insert within a profile block; delete from seq.
		    {
			indel_typ i;
			assert(state != 'I'); num_open++; 
			for(i=0; operation[o]=='I'; ){ j++; i++; num_insert++; o++; }
			o--; assert(i < gseq_ins_mask); 
			if(state != 'D' && state != 'd') indels[J0] = i; 
			else indels[J0] += i; // keep delete bit.
		    } break;
               default:
                        fprintf(stderr,"operations = %s\n",operation);
                        print_error("gsq_typ::initialize( ): input error"); break;
            }  state=operation[o];
	    assert(j <= LenSeq(realE));
        }
#if 1	// Bug Fix: Dealing with extensions (C-terminal flanks) in MAPGAPS... AFN: 5/04/2017.
	if(j < LenSeq(realE)){
	   // fprintf(stderr,"operation=%s\n",operation);
	   Int4 X0,x=(LenSeq(realE)-j);
	   f2r_typ *nf2r = new f2r_typ[N+x+3];
	   unsigned char   *ngapseq=0;
	   ngapseq = new unsigned char[N+x+2];
	   indel_typ *nindels = new indel_typ[N+x+3];
	   for(x=0; x <=J0; x++){ nf2r[x]=f2r[x]; ngapseq[x]=gapseq[x]; nindels[x]=indels[x]; }
	   delete f2r; f2r=nf2r; delete gapseq; gapseq=ngapseq; delete indels; indels=nindels;
	   while(j < LenSeq(realE)){ J0++; j++; gapseq[J0]=rseq[j]; indels[J0]=0; f2r[J0]=j; }
	}
#endif
	if(j != LenSeq(realE)){
	        fprintf(stderr,"operation=%s\n",operation);
	        fprintf(stderr,"offset realE=%d; ",OffSetSeq(realE));
	        fprintf(stderr,"C-terminal Extend=%d\n",CtermExtendSeq(realE));
		fprintf(stderr,"j=%d; LenSeq=%d\n",j,LenSeq(realE));
		assert(j == LenSeq(realE));
	}
	Int4 x=(LenSeq(realE)-j);
	assert(x >= 0);
#if 1	// Bug Fix!!	Fake seq extends beyond real seq so x is negative.
	// gseq_ins_mask & negative == a large number!!!!
	while(x < 0){ J0--; x++; }
#endif
	indels[J0] = (gseq_del_mask & indels[J0]) | (gseq_ins_mask & (indel_typ) x); 
	extendX = CtermExtendSeq(realE) + x;
	StrSeqID(id, DEFAULT_INFO_LINE_LEN_SEQ, realE); 
	StrSeqDescript(descript,DEFAULT_INFO_LINE_LEN_SEQ, realE);
	fakeE=MakeSeq(id,descript,offsetX,extendX,J0,gapseq);
	delete []gapseq;
}

void    gsq_typ::initialize(e_type E)
/****************************************************************************
  create 'fake' subsequence without indels.
 ****************************************************************************/
{
        Int4            N,j,k;
        unsigned char   *rseq;
        char            id[MAX_SEQ_DEFLINE_LENG + 100];
        char            descript[MAX_SEQ_DEFLINE_LENG + 100];
        Int4            offsetX=0,extendX=0;
        assert(realE == NULL); assert(fakeE == NULL);
        realE=E;
        offsetX=OffSetSeq(realE); // start relative to real full sequence.

        N = LenSeq(realE)+2;
        f2r = new f2r_typ[N+3];
        indels = new indel_typ[N+3];
        num_open=num_del=num_insert=0; indels[0]=0; f2r[0]=0;

        // 2. Add aligned-region sequence.
        for(rseq=SeqPtr(realE),j=0; j < LenSeq(realE); ){
            j++; indels[j] = 0; f2r[j] = j;
        }
        if(j != LenSeq(realE)){
                fprintf(stderr,"j=%d; LenSeq=%d\n",j,LenSeq(realE));
                assert(j == LenSeq(realE));
        }
        extendX = CtermExtendSeq(realE);
        StrSeqID(id, DEFAULT_INFO_LINE_LEN_SEQ, realE);
        StrSeqDescript(descript,DEFAULT_INFO_LINE_LEN_SEQ, realE);
        fakeE=MakeSeq(id,descript,offsetX,extendX,j,rseq);
}

#if 0	// comments....
	1. Allow deletions next to insertions.
	2. Allow overhangs (now removing with tweakcma -rmflank).
	3. Allow operations in both directions...Reverse as well.
	4. Allow deletions on the ends.
#endif

void	gsq_typ::initializeX(Int4 leftflank, Int4 rightflank, char *operation, 
	Int4 trace_length, Int4 start0, e_type E,Int4 *pos)
{
	e_type	oldE;
	if(fakeE != NULL) Free();
	realE = E;
	oldE=initializeX(leftflank, rightflank, operation,trace_length,start0,pos);
	if(oldE != NULL) NilSeq(oldE);
}

e_type	gsq_typ::initializeX(Int4 leftflank, Int4 rightflank, char *operation, 
	Int4 trace_length, Int4 start0, Int4 *pos)
/****************************************************************************
  create 'fake' subsequence with gaps and flanking 
  trace = "EDdmmmmmmmmmmmiiiMmmmmmmmmmmmmmmmmmmmmiiMmmmmmmmmmmE"
 ****************************************************************************/
{
	e_type		oldE;
	Int4		N,o,j,k,t,J0,flankL,flankR;
        unsigned char   *rseq,*gapseq;
	char		state,id[MAX_SEQ_DEFLINE_LENG + 100];
	char		descript[MAX_SEQ_DEFLINE_LENG + 100];
	Int4		offsetX=0,extendX=0,start;

	assert(leftflank >= 0 && rightflank >= 0); assert(realE != NULL);
	assert(operation[1]!='i'); 
	assert(operation[0]=='E'); 
	assert(operation[trace_length-2]!='i');
	if(operation[trace_length-1] != 'E'){
		fprintf(stderr,"operation = %s\n",operation);
		assert(operation[trace_length-1]=='E');
	}
	oldE=fakeE; fakeE=NULL; Free(); 
	flankL=MINIMUM(Int4,leftflank,start0-1); 
	start=start0-flankL-1;	// start relative to real subsequence.
	offsetX=OffSetSeq(realE)+start; // start relative to real full sequence.

	N = flankL + trace_length + rightflank;
	gapseq = new unsigned char[N+2];
	f2r = new f2r_typ[N+3];
	indels = new indel_typ[N+3];
	num_open=num_del=num_insert=0; indels[0]=start; f2r[0]=0; 

	// 1. Add leftflanking sequence.
	for(rseq=SeqPtr(realE),J0=0,j=start,k=1; k<= flankL; k++){
	    J0++; j++; gapseq[J0]=rseq[j]; indels[J0]=0; f2r[J0]=j;
	}
	// 2. Add aligned-region sequence.
        for(t=0,o=1,state='E'; operation[o] != 'E'; o++){
            switch(operation[o]){
               case 'M': if(pos){ t++; pos[t] = J0+1; }
	       case 'm': J0++; j++; gapseq[J0] = rseq[j];
			indels[J0] = 0; f2r[J0] = j; break;
               case 'D': if(pos){ t++; pos[t] = J0+1; }
	       case 'd': // deletion in sequence relative to profile.
                        J0++; gapseq[J0] = 0; num_del++;
			if(state != 'D' && state != 'd') num_open++;
			indels[J0] = gseq_del_mask; f2r[J0] = j; break;
               case 'i': // insert is between profile blocks;
                        J0++; j++; gapseq[J0] = rseq[j]; 
			indels[J0] = 0; f2r[J0] = j; break;
               case 'I': // Insert within a profile block; delete from seq.
		     { 
			indel_typ  i;
			assert(state != 'I'); num_open++; 
			for(i=0; operation[o]=='I'; ){ j++; i++; num_insert++; o++; }
			o--; assert(i < gseq_ins_mask); 
			if(state != 'D' && state != 'd') indels[J0] = i; 
			else indels[J0] += i; // keep delete bit.
		     } break;
               default:
                        fprintf(stderr,"operations = %s\n",operation);
                        print_error("gsq_typ::initialize( ): input error"); break;
            }  state=operation[o];
// fprintf(stderr,"j=%d; indels[%d]=%d\n",j,J0,gseq_ins_mask & indels[J0]);
        }
	// 3. Add rightflank region.
	flankR = MINIMUM(Int4,rightflank,LenSeq(realE)-j);
	for(Int4 i=1; i<= flankR; i++){
	    J0++; j++; gapseq[J0]=rseq[j]; indels[J0]=0; f2r[J0]=j;
	}
	Int4 x=(LenSeq(realE)-j);
#if 1	// BUG FIX!!!
	while(x < 0){ J0--; x++; }
	indels[J0] = (gseq_del_mask & indels[J0]) | (gseq_ins_mask & (indel_typ) x); 
// fprintf(stderr,"x=%d; j=%d; indels[%d]=%d\n",x,j,J0,insert(J0));
#endif
	extendX = CtermExtendSeq(realE) + x;
	StrSeqID(id, DEFAULT_INFO_LINE_LEN_SEQ, realE); 
	StrSeqDescript(descript,DEFAULT_INFO_LINE_LEN_SEQ, realE);
	fakeE=MakeSeq(id,descript,offsetX,extendX,J0,gapseq);
#if 0
	// NEED TO FIX THIS AS SEQUENCES WITH N-terminal flank are 
	// not working properly...[e.g., {(ADE)...]
	if(FakeToReal(LenSeq(fakeE)) > LenSeq(realE)){
		PutSeqInfo(stderr,fakeE);
		fprintf(stderr,"fake len = %d\n",LenSeq(fakeE));
		PutSeqInfo(stderr,realE);
		fprintf(stderr,"real len = %d\n",LenSeq(realE));
		assert(FakeToReal(LenSeq(fakeE)) <= LenSeq(realE));
	}
#endif
// assert(!(operation[1] == 'D' && operation[2] == 'I'));
// if(flag){ printf("%s\n",operation); Put(stderr,A); fprintf(stderr,"%s\n",operation); NilAlpha(A);}
	delete []gapseq;
	return oldE;
}

