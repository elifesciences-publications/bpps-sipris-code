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

#include "sma.h"
#include "blosum62.h"

sma_typ *MultiReadSMA(char *msafile, Int4 *n)
{
        sma_typ	*MA;
	char	str[1000],*s2;
	Int4	i,N;
        FILE	*fp,*fp2;

        fp=open_file(msafile,"","r");
        for(i=0; (s2=fgets(str,999,fp)) != NULL; ){
           if(strncmp(s2,"//",2)==0) i++;
        } fclose(fp);
        if(i==0) return NULL;
	N=i;

        NEW(MA,N+3,sma_typ);
        fp=open_file(msafile,"","r");
	fp2=tmpfile();
        s2=fgets(str,999,fp); fprintf(fp2,"%s",s2);
        for(i=1; (s2=fgets(str,999,fp)) != NULL; ){
           if(strncmp(s2,"//",2)==0){
		rewind(fp2);
		MA[i]=ReadSMA(fp2);
		fclose(fp2); fp2=tmpfile(); i++;
           }
	   fprintf(fp2,"%s",s2);
        }
	rewind(fp2); fclose(fp);
        MA[i]=ReadSMA(fp2);
        fclose(fp2);
	*n = i;
        return MA;
}

sma_typ	ReadSMA(char *msafile)
{
	FILE *fptr=open_file(msafile,"","r");
	sma_typ sma=ReadSMA(fptr); fclose(fptr);
	return sma;
}

sma_typ	ReadSMA(FILE *fptr)
{ 
    Int4		n,s,t,t0,N,ntyps;
    Int4 		i,j,e,glen,len,length;
    char		r,c,str[1000],id[100];
    char		*sq=NULL,x;
    fpos_t		fpos;
    Int4    		nseqs;          
    float   		**prob;         /** seg probability **/
    Int4    		**start,**end;  /** start[n][t]... **/
    sma_typ		MA;
    a_type		A;
    
    NEW(MA,1,comultaln_type);
    MA->A=A=MkAlpha(AMINO_ACIDS,PROT_BLOSUM62);
    NEW(MA->freq,nAlpha(A)+2,double);
    for(i=0; i<=nAlpha(A); i++) MA->freq[i]=(double)blosum62freq[i];

    // if(fscanf(fptr,"//%*[\n]ID   %s%*[\n]AC   %s%*[\n]DE   %s%*[\n]", 
	//	MA->msaid,MA->access,MA->de) != 3) return NilSMA(MA);
// fprintf(stderr,"DEBUG1\n");
    if(fscanf(fptr,"//%*[\n]ID   %s%*[\n]AC   %s%*[\n]",
		MA->msaid,MA->access) != 2) return NilSMA(MA);
// fprintf(stderr,"DEBUG2\n");
#if 1	// afn: 1/21/08 bug fix...
    if(fgets(str,90,fptr)==NULL) return NilSMA(MA);
// fprintf(stderr,"DEBUG3\n");
    if(str[0] != 'D' || str[1] != 'E') return NilSMA(MA);
// fprintf(stderr,"DEBUG4\n");
    for(i=2; isspace(str[i]); i++) if(str[i] == '\n') return NilSMA(MA);
// fprintf(stderr,"DEBUG5\n");
    if(!isprint(str[i])) return NilSMA(MA);
// fprintf(stderr,"DEBUG6\n");
    for(j=0; isprint(str[i]); i++,j++) {
	if(str[i] == '\n') break;
	MA->de[j] = str[i];
    } MA->de[j]=0;
    if(str[i+1] != 0) return NilSMA(MA);
// fprintf(stderr,"DEBUG7\n");
#endif
    // fprintf(stderr,"DE: %s\n",MA->de);
    if(fscanf(fptr,"LPR  %f.%*[\n]", &MA->lpr) !=1){ 
	fprintf(stderr,"READSMA error: LPR=%f\n",MA->lpr);
	print_error("fatal input error");
	return NilSMA(MA);
    }
//fprintf(stderr,"DEBUG8\n");
    if(fscanf(fptr,"NU   %d blocks; %d columns.%*[\n]",&ntyps,&MA->totcol)!=2){
		fprintf(stderr,"error: blks=%d; col=%d\n",ntyps,MA->totcol);
		 return NilSMA(MA);
    }
    fscanf(fptr,"SQ   %d sequences.%*[\n]",&nseqs);
    MA->ntyps=ntyps; MA->nseqs=nseqs;

    NEWP(MA->null,ntyps+2,char);
    NEW(MA->ignore,ntyps+2,BooLean);
    NEW(MA->length,ntyps+2,Int4);

    NEWP(MA->gnull,ntyps+2,char);
    NEW(MA->glength,ntyps+2,Int4);

    NEW(MA->cols,ntyps+2,Int4);
    NEW(MA->fieldmap,ntyps+2,float);
    NEW(MA->single_map,ntyps+2,float);
    NEW(MA->map_minus,ntyps+2,float);

    NEWP(MA->info,nseqs+2,char);
    NEWPP(MA->gseq,nseqs+2,char);
    NEWPP(MA->seq,nseqs+2,unsigned char);
    NEWP(MA->id,nseqs+2,char);
    NEWP(start,nseqs+2,Int4); NEWP(end,nseqs+2,Int4); 
    NEWP(prob,nseqs+2,float);
    for(n=1; n <= nseqs; n++) {
    	NEWP(MA->seq[n],ntyps+2,unsigned char);
    	NEWP(MA->gseq[n],ntyps+2,char);
    	NEW(start[n],ntyps+2,Int4); NEW(end[n],ntyps+2,Int4);
        NEW(prob[n],ntyps+2,float);
    }
    MA->start=start; MA->end=end; MA->prob=prob;

// fprintf(stderr,"DBBUG4\n");
// fprintf(stderr,"DEBUG9\n");
    for(n = 1; n <= nseqs; n++) {
	if(fgets(str,990,fptr)==NULL) return NilSMA(MA);
	NEW(MA->info[n],strlen(str)+2,char);
	strcpy(MA->info[n],str);
    }   
// fprintf(stderr,"DEBUG5\n");
// fprintf(stderr,"DEBUG10\n");
    for(t=1; t <= ntyps; t++) {
	MA->ignore[t] = FALSE;
	if(fgetpos(fptr,&fpos) != 0) print_error("fgetpos( ) failed!\n");
	if(fscanf(fptr,"BLK  %d: %d residues (%d with gaps; %d columns).\n",
		&t0,&MA->length[t],&MA->glength[t],&MA->cols[t]) != 4){
	// fprintf(stderr,"t = %d; t0 = %d; length = %d; glength = %d; cols = %d\n",
	// 			    t,t0,MA->length[t],MA->glength[t],MA->cols[t]);
	    if(fsetpos(fptr,&fpos) != 0) print_error("fsetpos( ) failed!\n");
    	    if(fscanf(fptr,"BLK  %d: %d residues (%d columns).\n",
		&t0,&MA->length[t],&MA->cols[t]) != 3) {
		    fprintf(stderr,"t = %d; t0 = %d; length = %d; cols = %d\n",
				    t,t0,MA->length[t],MA->cols[t]);
	    	    if(fsetpos(fptr,&fpos) != 0) print_error("fsetpos( ) failed!\n");
	    	    char buffer[1000];
	    	    fgets(buffer,990,fptr);
	    	    fprintf(stderr,"%s\n",buffer);
		    return NilSMA(MA);
	    } else MA->glength[t] = MA->length[t];
	} 
//fprintf(stderr,"DEBUG6\n");
// fprintf(stderr,"DEBUG11.%d\n",t);
	if(t != t0) return NilSMA(MA);
	if(MA->glength[t] < MA->length[t]) print_error("ReadSMA( ) fatal error");
#if 0
fprintf(stderr,"t = %d; glength[t]=%d; length[t] = %d\n",
	t,MA->glength[t],MA->length[t]);
#endif
// fprintf(stderr,"DEBUG12\n");
	if(fscanf(fptr,"MAP  %f (%f only; %f w/o).%*[\n]",
		&MA->fieldmap[t], &MA->single_map[t], &MA->map_minus[t]) != 3) 
			return NilSMA(MA);
	NEW(MA->null[t], MA->length[t]+2, char);
	NEW(MA->gnull[t], MA->glength[t]+5, char);
// fprintf(stderr,"DEBUG13\n");
	if(fscanf(fptr,"AL   %[.*^!_]%*[\n]",MA->gnull[t]) != 1) return NilSMA(MA);
#if 1	// Detection of array bounds read error with exit.
	// This appears to occur when files contain insertions next to deletions.
	Int4 should_be_null = MA->glength[t]+1; // just beyond the end of the array.
	if(MA->gnull[t][should_be_null] != 0){
		print_error("SMA file input error (check for adjacent deletions-insertions)\n");
	}
#endif
	for(j=i=0; (x=MA->gnull[t][i]) != 0; i++){
		if(x != '_') { MA->null[t][j] = x; j++; }
	} MA->null[t][j] = 0; len = j;
// if(MA->gnull[t] != NULL) fprintf(stderr,"%s\n",MA->gnull[t]); 
// fprintf(stderr,"DEBUG8\n");
// fprintf(stderr,"DEBUG14\n");
	if(len != MA->length[t]) return NilSMA(MA);
	MA->ignore[t]=FALSE;
	for(s=0; MA->null[t][s] != 0; s++){
	    if(MA->null[t][s] == '^') { MA->ignore[t]=TRUE; break; }
        }
        for(n=1; n<= nseqs; n++){
	    NEW(sq, MA->glength[t]+2, char); MA->gseq[n][t] = sq;
	    NEW(MA->seq[n][t], len+2 , unsigned char);
// fprintf(stderr,"DEBUG9\n");
// fprintf(stderr,"DEBUG14.%d\n",n);
            if(fscanf(fptr,"     %[a-zA-Z-.] %d-%d %s %f%*[\n]",sq,&MA->start[n][t],
			&MA->end[n][t],id, &MA->prob[n][t]) !=5){
				fprintf(stderr,"ReadSMA() input error: %s\n",sq);
				print_error("FATAL: input sequence defline not in msa format");
		    		return NilSMA(MA);
	    }
// fprintf(stderr,"%s\n",sq);
// fprintf(stderr,"DEBUG10\n");
	    
	    if(t==1) { NEW(MA->id[n],strlen(id)+2,char); strcpy(MA->id[n],id); }
	    for(j=i=0; (x=MA->gnull[t][i]) != 0; i++){
		if(x != '_') {
		    j++;
		    if(sq[i] == '-') MA->seq[n][t][j] = 0;
		    else MA->seq[n][t][j] = AlphaCode(sq[i],A);
		} 
	    }
	    if(j != len) print_error("ReadSMA( ) input error 2");
	}
    }
    MA->use = NULL;
    return MA;
}

sma_typ	NilSMA(sma_typ MA)
{
	Int4	n,t;

	NilAlpha(MA->A);
	free(MA->freq);
	if(MA->start!= NULL){
    	   for(t = 1; t <= MA->ntyps; t++){
// if(MA->gnull[t] != NULL) fprintf(stderr,"%s\n",MA->gnull[t]); 
		if(MA->null[t] != NULL) free(MA->null[t]); 
		if(MA->gnull[t] != NULL) free(MA->gnull[t]); 
	   }
	   for(n=1; n<=MA->nseqs; n++){ 
    	        for(t = 1; t <= MA->ntyps; t++){
			if(MA->seq[n][t] !=NULL) free(MA->seq[n][t]);
			if(MA->gseq[n][t] !=NULL) free(MA->gseq[n][t]);
		}
		free(MA->seq[n]);
		free(MA->gseq[n]);
		if(MA->id[n] != NULL) free(MA->id[n]);
		free(MA->start[n]); free(MA->end[n]); free(MA->prob[n]); 
		if(MA->info[n] != NULL) free(MA->info[n]);
	   }

	   free(MA->gnull); free(MA->glength); free(MA->gseq);
	   free(MA->null); free(MA->seq); free(MA->ignore);
	   free(MA->start); free(MA->end); free(MA->prob); free(MA->info);
    	   free(MA->length); free(MA->cols); free(MA->id);
    	   free(MA->fieldmap); free(MA->single_map); free(MA->map_minus);
	}
	if(MA->use != NULL) free(MA->use);
	free(MA);
	return NULL;
}

Int4	*gnullInsrtLenSMA(Int4 blk, sma_typ MA)
// This returns an array indicating how many residues are within each insert column.
{
    assert(blk > 0 && blk <= ntypSMA(MA));
    Int4        insrt_len,j,k;
    Int4	*gnull_insrt_len;
    char        state=0;
    char	*gnull=gnullSMA(blk,MA);

    NEW(gnull_insrt_len,glengthSMA(1,MA) + 5,Int4);
    for(state=0,insrt_len=j=0; j < glengthSMA(1,MA); j++){
        if(gnull[j]  == '_'){
           if(state == 'i') insrt_len++;
           else if(state == 'm'){ insrt_len=1; }
           else insrt_len=1;
           state = 'i';
        } else {
           if(state == 'm'){    // last state was a matching column.
                insrt_len=gnull_insrt_len[j]=0;
           } else if(state == 'i'){     // last state was a insert column.
              for(k = j-1; k >= 0 && gnull[k] == '_'; k--){
                gnull_insrt_len[k]=insrt_len;
              } insrt_len = 0;
           } else {     // last state was null column (this is the first positon)
              insrt_len=0;
           }
           state = 'm';
        }
    }
#if 0	// Debug...
    for(j=0; j < glengthSMA(1,MA); j++){
        fprintf(stderr,"gnull[%d] = %c %d\n",j,gnull[j],gnull_insrt_len[j]);
    }
#endif
    return gnull_insrt_len;
}

BooLean	*PutPurgeSMA(FILE *fptr, Int4 cutoff, sma_typ MA)
{ 
    Int4	n,t,m,s,N;
    BooLean	*use;
    
    use=purge_sma(cutoff, MA);
    for(N=0,n = 1; n <= MA->nseqs; n++) if(use[n]) N++;
    fprintf(fptr,"//\nID   %s\nAC   %s\nDE   %s\n",MA->msaid,MA->access,MA->de);
    fprintf(fptr,"LPR  %.2f.\n", MA->lpr);
    fprintf(fptr,"NU   %d blocks; %d columns.\n",MA->ntyps,MA->totcol);
    fprintf(fptr,"SQ   %d sequences.\n",N);
    for(n = 1; n <= MA->nseqs; n++) if(use[n]) fprintf(fptr,"%s",MA->info[n]);
    for(t=1; t <= MA->ntyps; t++) {
	fprintf(fptr,"BLK  %d: %d residues (%d columns).\n",
				t,MA->length[t],MA->cols[t]);
	fprintf(fptr,"MAP  %.1f (%.1f only; %.1f w/o).\n",
		MA->fieldmap[t], MA->single_map[t], MA->map_minus[t]);

	fprintf(fptr,"AL   %s\n",MA->null[t]);

        for(n=1; n<= MA->nseqs; n++){
	   if(use[n]){
	        fprintf(fptr,"     ");
		for(s=1; s<= MA->length[t]; s++) 
			fprintf(fptr,"%c",AlphaChar(MA->seq[n][t][s],MA->A));
	        fprintf(fptr," %d-%d %s %.2f\n",
			MA->start[n][t],MA->end[n][t],MA->id[n],
			MA->prob[n][t]);
	   }
	}
    }
    return use;
}

void	RePutSeqsSMA(FILE *fptr, e_type	*ListE, Int4 left, Int4 right, sma_typ MA)
{
	e_type	E;
	Int4	i,j,m,n,nblks,start,end,r,offset;
	char	c;

	nblks=ntypSMA(MA);
	for(i=1; i<=MA->nseqs; i++){
	  E = ListE[i];
	  offset = OffSetSeq(E);
	  start = MA->start[i][1] - offset;
	  end = MA->end[i][nblks] - offset;
	  PutSubSeq(fptr,start-left,end+right,E,MA->A);
	  fprintf(fptr,"\n");
	}
	fprintf(fptr,"\n\n");
}

void	PutClustalwSMA(FILE *fptr, Int4 blk, sma_typ MA)
/** output a Clustalw alignment of block blk **/
{
	Int4	i,j,m,n,maxid,nblks,s;
	char	c,*str,cmd[200];
        
	nblks=ntypSMA(MA);
	if(blk < 1 || blk > nblks) print_error("error in PutClustalwSMA( )");
	/** find max length of **/
	for(maxid=0,i=1; i<=MA->nseqs; i++){
	  str = MA->id[i];
	  n = strlen(MA->id[i]);
	  if(n > maxid) maxid = n;
	}
	fprintf(fptr,"CLUSTAL W (1.7) multiple sequence alignment\n\n\n");
	sprintf(cmd,"%%%ds      ",-maxid);
	for(i=1; i<=MA->nseqs; i++){
		fprintf(fptr,cmd,MA->id[i],MA->seq[i][blk]);
		for(s=1; s<= MA->length[blk]; s++) 
			fprintf(fptr,"%c",AlphaChar(MA->seq[i][blk][s],MA->A));
		fprintf(fptr,"\n");
	}
	fprintf(fptr,cmd," ");
	fprintf(fptr,"\n\n");
}

void	PutSelexSMA(FILE *fptr, e_type	*ListE, sma_typ MA)
{
	e_type	E;
	Int4	i,j,m,n,nblks,gap,start,end,r,offset,*maxgap;
	Int4	divide,top;
	char	c;

	nblks=ntypSMA(MA);
	NEW(maxgap,nblks+3,Int4);
	for(i=1; i<=MA->nseqs; i++){
	    for(m=1; m < nblks; m++){
		end=MA->end[i][m];
		start=MA->start[i][m+1];
		gap = start - end - 1;
		if(gap > maxgap[m]) maxgap[m]=gap;
	    }
	}
	for(i=1; i<=MA->nseqs; i++){
	  E = ListE[i];
	  offset = OffSetSeq(E);
	  PutSeqID2(fptr,E,20); 
	  for(m=1; m<=nblks; m++){
		start = MA->start[i][m] - offset;
		end = MA->end[i][m] - offset;
#if 1  // new gapped version
		fprintf(fptr,"%s", MA->gseq[i][m]);
#endif
#if 0  // old pregap version...
		for(j = start; j <= end; j++){
		   r = ResSeq(j,E); 
		   fprintf(fptr,"%c",AlphaChar(r,MA->A));
		}
#endif
		/** gaps here... **/
		if(m < nblks){
		   start = MA->end[i][m] - offset + 1;
		   end = MA->start[i][m+1] - offset;
#if 1		/*** NEW ***/
		   n = end - start;
		   divide = n/2;
		   top = maxgap[m]-(n - divide);
		   for(j=start,gap=0; gap < maxgap[m]; gap++){
			if(gap <= divide){
				r = ResSeq(j,E);  j++;
				c = AlphaChar(r,MA->A);
				c = tolower(c);
			} else if(gap > top){ 
				r = ResSeq(j,E);  j++;
				c = AlphaChar(r,MA->A);
				c = tolower(c);
			} else c = '-';
			fprintf(fptr,"%c",c);
		   }
#endif
#if 0
		   for(gap=0,j = start; j < end; j++){
			r = ResSeq(j,E); 
			c = AlphaChar(r,MA->A);
			c = tolower(c);
			fprintf(fptr,"%c",c); gap++;
		   }
		   for(j=gap; j < maxgap[m]; j++) {
			fprintf(fptr,"-"); 
		   }
#endif
		} 
	  }
	  fprintf(fptr,"\n");
	}
	fprintf(fptr,"\n\n");
	free(maxgap);
}

Int4    max_int_sma_array(Int4 *I, Int4 *L)
/** sets *I to integer repeated most often; returns number of repeats **/
{
        Int4    i,j,N,maxI,n,x,t;

        /** BubbleSortIArray(Int4 *L); **/
        for(i=1; i <= L[0]; i++){
                for(j=i;j > 1 && (L[j-1] > L[j]);j--){
                        t=L[j]; L[j]=L[j-1]; L[j-1]=t;
                }
        }
#if 0	/** Debug: output array **/
        fprintf(stderr,"\n\t[");
        for(i=1; i <= L[0]; i++){ 
                fprintf(stderr,"%4d",L[i]);
                if(i%10==0) fprintf(stderr,"\n\t");
        }
        fprintf(stderr," ]\n");

#endif
        i=1; maxI=-1; N = 0;
        while(i <= L[0]){
           for(x=L[i], n=0; i <= L[0] && L[i] == x; ){
                n++; i++;
           }
           if(n > N) { N= n; maxI = x; }
        } 
        *I = maxI;
        return N;
}

Int4	GapLengthSMA(Int4 m, Int4 n, sma_typ MA)
// Return gap length between motif m and m+1 for seq n
{
	if(m < 1 || m > MA->ntyps || n < 1 || n > MA->nseqs)
	   print_error("GapLengthSMA( ) input error");
	return (MA->start[n][m+1] - MA->end[n][m] - 1);
}

void	PutGapsSMA(FILE *fp, sma_typ MA)
// get totGaps using nseqSMA(MA) macro 
{
	Int4	t,n,g,t1,t2;
	h_type	H;
	char	str[100];

   for(t1 = 1; t1 < MA->ntyps; t1++){
	t2 = t1 +1;
	sprintf(str,"gap lengths between block %d and %d",t1,t2);
	H =Histogram(str,0,100,1.0);
        for(n=1; n <= MA->nseqs; n++) {
           g = MA->start[n][t2] - MA->end[n][t1] - 1;
	   IncdHist(g,H);
	}
	PutHist(fp,60,H); NilHist(H);
   }
}


BooLean	DiffSMA(FILE *fptr, sma_typ MA1, sma_typ MA2)
/*********************************************************************
  (print the difference between MA1 and MA2)
  Optimally align MA1 and MA2 
  Want to add statistics, recombination, and subalignment routines later.
  Will allow speed up of Gibbs alignment algorithm.
  Want to incorporate these into cmsa datatype. 
  return TRUE if a outfile was created.
 *********************************************************************/
{ 
    Int4	n,t1,t2,os,max,maxos,max_t,max_n;
    Int4	*offset,Total;
    sma_typ	ma1,ma2,MA;
#if 1
    Int4	*rm_front,*rm_end;
    FILE	*fp;
#endif
    
    /** 1. Check to see that alignments use the same sequences. **/
    if(MA1->nseqs != MA2->nseqs){
    	fprintf(stderr,"inconsistent alignments; can't be compared\n");
	return FALSE;
    } else {
      for(n = 1; n <= MA1->nseqs; n++) {
	if(strncmp(MA1->info[n],MA2->info[n],50) != 0) {
    		fprintf(stderr,"inconsistent alignments; can't be compared\n");
		return FALSE;
	}
      }
    }
    if(MA1->ntyps <= MA2->ntyps) { ma1 = MA1; ma2 = MA2; }
    else { ma2 = MA1; ma1 = MA2; }
#if 1
    NEW(rm_front,MA2->ntyps +3,Int4);
    NEW(rm_end,MA2->ntyps +3,Int4);
#endif
    /** 2. For each block in MA1 find best block in MA2. **/
    NEW(offset,ma1->nseqs +3, Int4);
    for(Total=0,t1=1; t1 <= ma1->ntyps; t1++) {
	max_n = 0;
        for(t2=1; t2 <= ma2->ntyps; t2++) {
            for(n=1; n<= ma1->nseqs; n++){
		offset[n] = ma2->start[n][t2] - ma1->start[n][t1];
	    } offset[0]=ma1->nseqs;  /** store length of array **/
#if 0
	    fprintf(stderr,"motif %d versus %d: \n", t1,t2);
#endif
	    n = max_int_sma_array(&os,offset);
	    if(n > max_n){ max_n = n; maxos = os; max_t = t2; }
#if 0
	    fprintf(stderr,"  %d matches (offset = %d)\n---------\n\n", n,os);
#endif
	}
	Total += (ma1->nseqs - max_n);
	fprintf(stderr,"BLKS %d & %d: alignments differ by %d/%d sites (os = %d)\n",
		t1,max_t,ma1->nseqs - max_n,ma1->nseqs, maxos);
	fprintf(stderr," %d versus %d columns\n", ma1->length[t1], 
		ma2->length[max_t]);
	if(TRUE || strcmp(ma1->null[t1],ma2->null[max_t]) != 0){
		fprintf(stderr," column configuration difference: \n\t");
		if(maxos < 0){
			for(os = maxos; os < 0; os++) fprintf(stderr," ");
			fprintf(stderr,"%s\n\t%s\n", ma1->null[t1],ma2->null[max_t]);
		} else {
			fprintf(stderr,"%s\n\t", ma1->null[t1]);
			for(os = maxos; os > 0; os--) fprintf(stderr," ");
			fprintf(stderr,"%s\n", ma2->null[max_t]);
		}
	}
#if 0	/** report gaps here **/
	if(t2 > 1){
            for(n=1; n<= ma1->nseqs; n++){
		ma2->start[n][t2] - ma1->start[n][t1];
	    }
	}
#endif
#if 1
	if(maxos < 0){ rm_front[t1] =  0; } else { rm_front[t1] =  maxos; }
	n  = ((ma1->length[t1] - maxos) - ma2->length[max_t]);
	rm_end[t1] = MAXIMUM(Int4,n,0);
#endif
    }
    /** 3. Use dynamic programming to find the best alignment of ma1 with ma2. **/
    free(offset);
    fprintf(stderr,"Alignments differ by a total of %d sites\n\n",Total);
#if 1	/*** TEST ***/
    fp= open_file("junk_diffmsa",".msa","w");
    PutSMA(fp, ma1); fclose(fp);
    MA = ReadSMA("junk_diffmsa.msa");
    if(MA == NULL){ free(rm_front); free(rm_end); return FALSE; }
    for(t1=1; t1 <= ma1->ntyps; t1++) {
      if(rm_front[t1] > 0){
	if(t1 < 1 || t1 > MA->ntyps || (MA->length[t1] - abs(rm_front[t1])) < 4){
		NilSMA(MA); free(rm_front); free(rm_end); return FALSE;
	} 
        fp= open_file("junk_diffmsa",".msa","w");
	PutTruncateBlkSMA(fp, t1, rm_front[t1], MA);
        fclose(fp); NilSMA(MA);
        MA = ReadSMA("junk_diffmsa.msa");
	if(MA == NULL){ free(rm_front); free(rm_end); return FALSE; }
      }
      if(rm_end[t1] > 0){
	if(t1 < 1 || t1 > MA->ntyps || (MA->length[t1] - abs(rm_end[t1])) < 4){
		NilSMA(MA); free(rm_front); free(rm_end); return FALSE;
	} 
        fp= open_file("junk_diffmsa",".msa","w");
        PutTruncateBlkSMA(fp, t1, -rm_end[t1], MA);
        fclose(fp); NilSMA(MA);
        MA = ReadSMA("junk_diffmsa.msa");
	if(MA == NULL){ free(rm_front); free(rm_end); return FALSE; }
      }
    }
    NilSMA(MA);
#endif
#if 1
    free(rm_front); free(rm_end);
#endif
    return TRUE;
}

void	ReOrderSMA(FILE *fptr, char *file, sma_typ MA)
/** Put MA using the order of sequences in data **/ 
{ 
    Int4	i,n,m,t,s,len,N,*map,T,start,end;
    char	*info,str0[100],*str;
    e_type	E;
    a_type	A;
    ss_type	data;
    BooLean	flag = TRUE;
    
    A = MkAlpha(AMINO_ACIDS,PROT_BLOSUM62);
    data = SeqSet(file,A);
    N = NSeqsSeqSet(data);
    T = MA->ntyps;
    if(MA->nseqs != N) {
    	  fprintf(stderr,"number sequences inconsistent with alignment; can't be reordered\n");
	  flag = FALSE;
    } else {
      NEW(map,N+2,Int4);
      for(m=1; m<= N; m++){
	   E = SeqSetE(m,data);
	   info = SeqKey(E);
	   for(i=0;!isspace(info[i]); i++) str0[i]=info[i]; str0[i]=0;
	   len = strlen(str0); 
	   start = OffSetSeq(E) + 1;
	   end = OffSetSeq(E) + LenSeq(E);
fprintf(stderr,"sequence %d: %s\n",m,str0);
           for(n=1; n <= N; n++){
		str = MA->info[n] + 5;
		s = strlen(str);
		s = MINIMUM(Int4,len,s);
		if(strncmp(info,str,s) == 0){
		    if(MA->start[n][1] >= start && MA->end[n][T] <= end){
fprintf(stderr,"sequence %d: %s (%d <= %d <=%d <= %d)\n",
	n,str,start,MA->start[n][1],MA->end[n][T],end);
		     if(map[m] == 0) map[m] = n;
		     else flag = FALSE;
		    }
		}
	   }
      }
      for(n=1; n <= N; n++){
	if(map[n] == 0) {
	  fprintf(stderr,"sequence %d: ",n);
    	  fprintf(stderr,"sequences inconsistent with alignment; can't be reordered\n");
	  flag = FALSE;
	}
      }
    }
 if(flag){
    fprintf(fptr,"//\nID   %s\nAC   %s\nDE   %s\n",MA->msaid,MA->access,MA->de);
    fprintf(fptr,"LPR  %.2f.\n", MA->lpr);
    fprintf(fptr,"NU   %d blocks; %d columns.\n",MA->ntyps,MA->totcol);
    fprintf(fptr,"SQ   %d sequences.\n",MA->nseqs);
    for(n = 1; n <= MA->nseqs; n++){
	m = map[n];
	fprintf(fptr,"%s",MA->info[m]);
    }
    for(t=1; t <= MA->ntyps; t++) {
	fprintf(fptr,"BLK  %d: %d residues (%d columns).\n",
				t,MA->length[t],MA->cols[t]);
	fprintf(fptr,"MAP  %.1f (%.1f only; %.1f w/o).\n",
		MA->fieldmap[t], MA->single_map[t], MA->map_minus[t]);

	fprintf(fptr,"AL   %s\n",MA->null[t]);

        for(n=1; n<= MA->nseqs; n++){
		m = map[n];
	        fprintf(fptr,"     ");
		for(s=1; s<= MA->length[t]; s++) 
			fprintf(fptr,"%c",AlphaChar(MA->seq[m][t][s],MA->A));
	        fprintf(fptr," %d-%d %s %.2f\n",
			MA->start[m][t],MA->end[m][t],MA->id[m],
			MA->prob[m][t]);
	}
    }
  }
    NilSeqSet(data); NilAlpha(A); free(map);
}

void	PutTruncateBlkSMA(FILE *fptr, Int4 blk, Int4 remove, sma_typ MA)
{ 
    Int4	i,j,n,t,m,s,start,end,sos,eos;
    char	c,*str;

    if(blk < 1 || blk > MA->ntyps || (MA->length[blk] - abs(remove)) < 2)
	print_error("PutTruncateBlkSMA( ) input error");
    fprintf(fptr,"//\nID   %s\nAC   %s\nDE   %s\n",MA->msaid,MA->access,MA->de);
    fprintf(fptr,"LPR  %.2f.\n", MA->lpr);
    fprintf(fptr,"NU   %d blocks; %d columns.\n",MA->ntyps,MA->totcol);
    fprintf(fptr,"SQ   %d sequences.\n",MA->nseqs);
    for(n = 1; n <= MA->nseqs; n++) fprintf(fptr,"%s",MA->info[n]);
    for(t=1; t <= MA->ntyps; t++) {
	if(t==blk){
    	   NEW(str,MA->length[blk] +3,char);
	   if(remove < 0){ sos = 0; eos = -remove; }
	   else { sos = remove; eos = 0; }
	   start=1+sos;
	   end = MA->length[t] - eos;
	   for(n=0,m=0,s=sos; s < end; s++) {
		c = MA->null[t][s]; str[m++] = c;
		if(c != '.') n++;
	   }
	   str[m]=0;
	   fprintf(fptr,"BLK  %d: %d residues (%d columns).\n",
				t,MA->length[t]-abs(remove),n);
	} else {
	   fprintf(fptr,"BLK  %d: %d residues (%d columns).\n",
				t,MA->length[t],MA->cols[t]);
	   sos=eos=0;
	   start=1; end=MA->length[t];
	}
	fprintf(fptr,"MAP  %.1f (%.1f only; %.1f w/o).\n",
		MA->fieldmap[t], MA->single_map[t], MA->map_minus[t]);

	if(t==blk){ fprintf(fptr,"AL   %s\n",str); free(str); }
	else fprintf(fptr,"AL   %s\n",MA->null[t]);

        for(n=1; n<= MA->nseqs; n++){
	    fprintf(fptr,"     ");
	    if(MA->gseq != NULL){
#if 1 // temporary fix... Will not work right for all input...
	        sos=eos=0;
		for(s=1; s < start; s++){
		   c = MA->gseq[n][t][s-1];
		   if(c != '-' && c != '.') sos++;
		}
		for(i=1,s = start; s<= end; i++,s++)
			fprintf(fptr,"%c",MA->gseq[n][t][s-1]);
		while(s <= MA->length[t]) {
		   c = MA->gseq[n][t][s-1];
		   if(c != '-' && c != '.') eos++;
		   s++;
		} 
	        fprintf(fptr," %d-%d %s %.2f\n",
			MA->start[n][t]+sos,MA->end[n][t]-eos,
			MA->id[n], MA->prob[n][t]);
#endif
	    } else {
		for(s = start; s<= end; s++) 
			fprintf(fptr,"%c",AlphaChar(MA->seq[n][t][s],MA->A));
	        fprintf(fptr," %d-%d %s %.2f\n",
			MA->start[n][t]+sos,MA->end[n][t]-eos,
			MA->id[n], MA->prob[n][t]);
	    }
	}
    }
}

void	PutSMA(FILE *fptr, sma_typ MA)
{ 
    Int4	n,t,m,s;
    
    fprintf(fptr,"//\nID   %s\nAC   %s\nDE   %s\n",MA->msaid,MA->access,MA->de);
    fprintf(fptr,"LPR  %.2f.\n", MA->lpr);
    fprintf(fptr,"NU   %d blocks; %d columns.\n",MA->ntyps,MA->totcol);
    fprintf(fptr,"SQ   %d sequences.\n",MA->nseqs);
    for(n = 1; n <= MA->nseqs; n++) fprintf(fptr,"%s",MA->info[n]);
    for(t=1; t <= MA->ntyps; t++) {
	fprintf(fptr,"BLK  %d: %d residues (%d columns).\n",
				t,MA->length[t],MA->cols[t]);
	fprintf(fptr,"MAP  %.1f (%.1f only; %.1f w/o).\n",
		MA->fieldmap[t], MA->single_map[t], MA->map_minus[t]);

	fprintf(fptr,"AL   %s\n",MA->null[t]);

        for(n=1; n<= MA->nseqs; n++){
	        fprintf(fptr,"     ");
		for(s=1; s<= MA->length[t]; s++) 
			fprintf(fptr,"%c",AlphaChar(MA->seq[n][t][s],MA->A));
	        fprintf(fptr," %d-%d %s %.2f\n",
			MA->start[n][t],MA->end[n][t],MA->id[n],
			MA->prob[n][t]);
	}
    }
}

#if 1
	// TEST routine to find out how to improve alignments
void    PutDiffSMA(FILE *fptr, sma_typ MA1, sma_typ MA2)
{
    Int4        n,t,m,s,x,len;
	char	c,d;
   
    fprintf(fptr,"WARNING: this routine is full of bugs!!!\n");
    for(n=1; n<= MA1->nseqs; n++){
        fprintf(fptr,"%s = %s:\n",MA1->id[n],MA2->id[n]);
        for(t=1; t <= MA1->ntyps; t++) {
	    x = labs(MA1->start[n][t] -MA2->start[n][t]);
	    len = MINIMUM(Int4,MA1->length[t],MA2->length[t]);
	    if(x >= len) d = '*'; else d = ' ';
	    if(MA1->start[n][t] != MA2->start[n][t]) c = '!'; else c = ' ';
                fprintf(fptr,"   %c%c",c,d);
                for(s=1; s<= MA1->length[t]; s++)
                        fprintf(fptr,"%c",AlphaChar(MA1->seq[n][t][s],MA1->A));
                fprintf(fptr," %d-%d %.2f   \n",
                        MA1->start[n][t],MA1->end[n][t],MA1->prob[n][t]);
	    if(MA1->length[t] == MA2->length[t] && 
		MA1->start[n][t] == MA2->start[n][t])
                fprintf(fptr,"     (same)  \n\n");
	    else {
                fprintf(fptr,"     ");
                for(s=1; s<= MA2->length[t]; s++)
                        fprintf(fptr,"%c",AlphaChar(MA2->seq[n][t][s],MA2->A));
                fprintf(fptr," %d-%d %.2f   \n\n",
                        MA2->start[n][t],MA2->end[n][t], MA2->prob[n][t]);
	    }
        }
        fprintf(fptr,"\n\n");
    }
}

#endif

BooLean *ConservedSMA(Int4 t, Int4 s, sma_typ MA)
{
	Int4	i,j,n,e,end,N,start,length,*observed,gap,ntyp;
	Int4	h,p,aliphatic,polar,consH,consP,consT,consMax,obs;
	float	percentH,percentP,percentT,percentDiff,percentA,percentMax;
	double	info,total,d,f,q,cbp_cut=-3.0;
	a_type	A=MA->A;
	char	r,c,str[50];
	char	*cons_state,state,laststate;  
	BooLean	*conserved;

    ntyp = MA->ntyps;
    if(t < 0 || t > ntyp) return NULL;
    length = MA->length[t];
    if(s < 0 || s > length) return NULL;
    if(MA->use == NULL) MA->use=purge_sma(250, MA);
    NEW(observed,nAlpha(A)+2,Int4);
    for(n=1; n<=MA->nseqs; n++){
	if(MA->use[n]){ r=MA->seq[n][t][s]; observed[r]++; }
    }
    conserved=conserved_sma(observed,cbp_cut,MA);
    free(observed);
    return conserved;
}

BooLean	*conserved_sma(Int4 *observed, double cutoff, sma_typ MA)
/* return the conserved residues in column "observed". */
{
	Int4	i,r,r2,nres,n;
	a_type	A=MA->A;
	double	total,p,q,f,*freq=MA->freq,low_cut;
	BooLean	*conserved,flag,debug=FALSE;

	assert(cutoff <= 0.0);
	low_cut=cutoff - 0.7;
	if(observed != NULL){
	    nres=0;
	    NEW(conserved,nAlpha(A)+2,BooLean);
	    for(total=0.0, r=1; r<=nAlpha(A); r++){
		total += (double)observed[r];
	    }
	    if(debug){
		fprintf(stderr,"\n");
		for(r = 1; r <= nAlpha(A); r++){
		   if(observed[r] > 0){
			q = freq[r];
			fprintf(stderr,"%c: %d/%d (p=%.1f)\n", 
				AlphaChar(r,A),
				observed[r], (Int4)total,
				Log10CBP(observed[r], total, q));
   		   }
		}
		fprintf(stderr,"\n");
	    }
	    /** 1. find clearly conserved residues **/
	    for(f=0.0, n=0, r=1; r <= nAlpha(A); r++){
		if(observed[r] > 2){
		   q = freq[r];
		   p = Log10CBP(observed[r], total, q);
		   if(p <= cutoff) {
			if(debug) fprintf(stderr," %c: (%.1f)[%d/%d]", 
				AlphaChar(r,A), 
				p,observed[r], (Int4)total);
			n += observed[r]; f += q; 
			conserved[r] = TRUE; nres++;
		    }
		}
	    }
            /** 2. find moderately conserved related residue pairs **/
            if(nres==0){
              for(n=0, f=0.0, r=1; r < nAlpha(A); r++){
                if(!conserved[r] && observed[r] > 2){
                  q = freq[r];
		  p = Log10CBP(observed[r], total, q);
                  if(p <= low_cut){
                    for(r2 = 1; r2 <= nAlpha(A); r2++){
                        if(!conserved[r2] && r2 != r && observed[r2] > 2 
                                        && valAlphaR(r,r2,A) > 0){
                           q = freq[r2];
		  	   p = Log10CBP(observed[r2], total, q);
                           if(p <= low_cut) {
				if(debug) fprintf(stderr," %c %c (%.1f)", 
					AlphaChar(r,A), AlphaChar(r2,A), p);
                                nres+=2;
				n += observed[r] + observed[r2];
				f += freq[r] + freq[r2];
                                conserved[r]=conserved[r2]=TRUE;
                           }
                        }
                    }
                  }
                }
              }
            } 
	    if(nres==0){ free(conserved); return NULL; }
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
				conserved[r2] = TRUE;
			   }
			}
		   }
		  }
		}
	      }
	    } while(flag);
	    /** 3. reset to off if total conserved residues too low **/
	    for(q=0.0,r2=0, r = 1; r <= nAlpha(A); r++){
	      if(conserved[r]){ r2 += observed[r]; q += freq[r]; }
	    }
	    p = Log10CBP(r2, total, q);
	    if(debug) fprintf(stderr,"  (%.0f%c conserved)(p=%.1f)\n",
		100.0*(double)r2/total,'%',p);
	    if(p > -3.0) { free(conserved); return NULL; }
	    p = (double)r2/total;
	    if(p < 0.20) { free(conserved); return NULL; }
	    return conserved;
	} else return NULL;
}

Int4    CountsSMA(Int4 *counts, sma_typ MA)
{
	Int4	n,t,s,r,total=0;
	a_type	A=MA->A;

    for(r=0; r <= nAlpha(A); r++) counts[r]=0;
    for(t=1; t <= MA->ntyps; t++) {
    	for(n = 1; n <= MA->nseqs; n++) {
	      for(s=1; s<= MA->length[t]; s++) {
		r = MA->seq[n][t][s];
		counts[r]++;
	      }
	      total+=MA->length[t];
	}
    }
    return total; 
}

Int4    *GetGapsSMA(Int4 maxgap, float minprob, Int4 t1, Int4 t2, sma_typ MA)
/** get totGaps using nseqSMA(MA) macro ***/
{
	Int4	*scoreGap,t,n,g;

	if(t1 >= t2 || t1 < 1 || t2 > MA->ntyps) 
				print_error("GetGapsSMA( ) error 1");
        NEW(scoreGap, maxgap+2, Int4);;
        for(n=1; n <= MA->nseqs; n++) {
           g = MA->start[n][t2] - MA->end[n][t1] - 1;
	   if(g < 0) g = 0;
#if 0
           if(g <= maxgap && MA->prob[n][t1] >= minprob 
				&& MA->prob[n][t2] >=minprob) scoreGap[g]++;
#endif
           if(g <= maxgap && MA->prob[n][t2] > minprob) scoreGap[g]++;
	}
	return scoreGap;
}

void    SMA2RTF(FILE *fptr, double cbp_cut, double infoLO, double infoHI, sma_typ MA)
{ SMA2RTF(fptr,cbp_cut,infoLO,infoHI,0,MA); }

void    SMA2RTF(FILE *fptr, double cbp_cut, double infoLO, double infoHI, 
	ss_type key_seq, sma_typ MA)
/** PutSites to create a figure for publication **/
/** 1 black;2 blue;3 cyan; 4 green; 5 magenta; 6 red; 7 yellow; 8 white **/
/** 9 lt blue;10 lt cyan; 11 lt green;12 ; 13; 14 lt magenta; 15 dk grey; 16 lt grey **/
{
	Int4	t,i,j,n,s,e,end,N,start,length,*observed,gap,ntyp;
	Int4	h,p,aliphatic,polar,consH,consP,consT,consMax,obs;
	float	percentH,percentP,percentT,percentDiff,percentA,percentMax;
	double	info,total,d,f,q;
	a_type	A=MA->A;
	char	r,c,str[50];
	char	*cons_state,state,laststate;  
	BooLean	**conserved,flag,*use;

    assert(cbp_cut < 0.0 && infoLO >= 0.1 && infoLO < infoHI);
    ntyp = MA->ntyps;
    fprintf(fptr,"{\\rtf1\\ansi \\deff4\\deflang1033");
    fprintf(fptr,"{\\fonttbl{\\f0\\froman\\fcharset0\\fprq2 Tms Rmn;}\n");
    fprintf(fptr,"{\\f1\\fswiss\\fcharset0\\fprq2 Arial;}\n");
    fprintf(fptr,"{\\f2\\fswiss\\fcharset0\\fprq2 Arial Narrow;}\n");
    fprintf(fptr,"{\\f3\\fmodern\\fcharset0\\fprq1 Courier New;}\n");
    fprintf(fptr,"{\\f4\\froman\\fcharset0\\fprq2 Times New Roman;}\n");
    fprintf(fptr,"{\\f5\\fswiss\\fcharset0\\fprq2 System;}\n");
    fprintf(fptr,"{\\f6\\fmodern\\fcharset0\\fprq1 Courier New;}}\n");
    fprintf(fptr,"{\\colortbl;\\red0\\green0\\blue0;\\red0\\green0\\blue255;");
    fprintf(fptr,"\\red0\\green255\\blue255;\\red0\\green255\\blue0;");
    fprintf(fptr,"\\red255\\green0\\blue255;\\red255\\green0\\blue0;");
    fprintf(fptr,"\\red255\\green255\\blue0;\\red255\\green255\\blue255;");
    fprintf(fptr,"\\red0\\green0\\blue128;\\red0\\green128\\blue128;");
    fprintf(fptr,"\\red0\\green128\\blue0;\\red128\\green0\\blue128;");
    fprintf(fptr,"\\red128\\green0\\blue0;\\red128\\green128\\blue0;");
    /** 15 dk grey; 16 lt grey (original) **/
    fprintf(fptr,"\\red128\\green128\\blue128;\\red192\\green192\\blue192;}");
    fprintf(fptr,"{\\stylesheet{\\widctlpar \\f4\\fs18 \\snext0 Normal;}");
    fprintf(fptr,"{\\*\\cs10 \\additive Default Paragraph Font;}}\n");

#if 1
    /*** 11 x 17 portrait ***/
    fprintf(fptr,"\\paperw15840\\paperh24480\\margl720\\margr720");
    fprintf(fptr,"\\margt1080\\margb1080 ");
    fprintf(fptr,"\\widowctrl\\ftnbj\\aenddoc\\formshade \\fet0\\sectd ");
    fprintf(fptr,"\\linex0\\endnhere\\sectdefaultcl \n");
#endif
#if 0
    /*** landscape ***/
    fprintf(fptr,"\\paperw15840\\paperh12240\\margl720\\margr720");
    fprintf(fptr,"\\margt1080\\margb1080 ");
    fprintf(fptr,"\\widowctrl\\ftnbj\\aenddoc\\formshade \\fet0\\sectd ");
    fprintf(fptr,"\\lndscpsxn\\psz1\\linex0\\endnhere \n");
#endif

    fprintf(fptr,"\\pard\\plain \\widctlpar \\f3\\fs18 \n");
    // fprintf(fptr,"\\pard\\plain \\widctlpar \\f3\\fs18 \\sl180\\slmult1 ");
    fprintf(fptr,"\\pard \\widctlpar \\sl180\\slmult1 ");
    use=purge_sma(250, MA);
    NEW(observed,nAlpha(A)+2,Int4);
    h_type HG=Histogram("information content in bits",0,100,0.05);
    for(t=1; t <= ntyp; t++){
        h_type HG2=Histogram("Kyte-Doolittle hydrophobicity",0,100,0.2);
	length = MA->length[t];
	NEWP(conserved, length+2, BooLean);
	NEW(cons_state, length+2, char);
	for(s=1; s <= length; s++){
	    for(r=0; r<=nAlpha(A); r++) observed[r]=0;
	    for(n=1; n<=MA->nseqs; n++){
		if(use[n]){ r=MA->seq[n][t][s]; observed[r]++; }
	    }
	    conserved[s] = conserved_sma(observed,cbp_cut,MA);
	    obs=consMax=consT=consH=consP=aliphatic=polar=0;
	    for(r=1; r<= nAlpha(A); r++){
		n = observed[r];
		sprintf(str,"%c",AlphaChar(r,A));
		if(strstr("AILVFYWMCP",str) != NULL) h=n; else h=0;
		if(strstr("DEKRHSTQN",str) != NULL) p=n; else p = 0;
		if(conserved[s] != NULL && conserved[s][r]){
		    consMax = MAXIMUM(Int4,consMax,n);
		    consT+=n; consH+=h; consP+=p;
		}
		aliphatic+=h; polar+=p; obs+=n;
	    }
/******* compute relative entropy of postions *******/
            for(total=0.0,r=1; r <= nAlpha(A); r++){ total += (double) observed[r]; }
	    double kd_score = 0.0;
            for(info=0.0,r=1; r <= nAlpha(A); r++){
		if(observed[r] > 0){
			kd_score += observed[r]*blsm62kyte_doolittle[r];
                        f = (double) observed[r]/total;
                        q = MA->freq[r];
                        info += f*log(f/q)/log(2.0);
                } 
            }
// fprintf(stderr,"info[%d][%d] = %.2f\n",t,s,info);
	    kd_score /= total;
            IncdHist(kd_score,HG2); IncdHist(info,HG);
/******* relative entropy of postions *******/
	    percentT = ((float)consT/(float)obs);
	    // percentH = ((float)consH/(float)obs);
	    percentH = ((float)consH/(float)consT);
	    percentP = ((float)consP/(float)obs);
	    percentMax = ((float)consMax/(float)obs);
	    percentA = ((float)aliphatic/(float)obs);
	    percentDiff = (double) fabs(percentH-percentP);
	    if(percentA >= 0.80 && percentH >= 0.70){
	           for(r=1; r<= nAlpha(A); r++){
			sprintf(str,"%c",AlphaChar(r,A));
			if(strstr(" AILVFYWMCP",str) != NULL &&
			    			!conserved[s][r]){
				conserved[s][r]=TRUE;
			}
		   }
	    }
	    if(info >= infoHI){ 
		if(kd_score > 5.0) cons_state[s] = 'H'; else cons_state[s] = 'I'; 
	    } else if(info >= infoLO){ 
		if(kd_score > 5.0) cons_state[s] = 'h'; else cons_state[s] = 'i'; 
	    } else if(percentH >= 0.98 && percentA >= 0.70){ 
		cons_state[s]='h'; 
	    } else if(percentA >= 0.90 && percentH >= 0.70 ){
		cons_state[s] = 'a';
	    } else if(percentT < 0.50){
		if(percentDiff < 0.33) cons_state[s] = 'w';
		else cons_state[s] = 'W';
	    } else if(percentT < 0.90){
		if(percentDiff < 0.66) cons_state[s] = 'm';
		else cons_state[s] = 'M';
	    } else {	// i.e., if(percentT >= 0.90)...
		if(kd_score > 5.0) cons_state[s] = 'a';
		else cons_state[s] = 's'; 
	    }
	}
    	PutHist(stderr,60,HG2);  NilHist(HG2);
        fprintf(fptr,"\\pard \\widctlpar \\sl180\\slmult1 ");
	BooLean	*show_position=0;
	if(key_seq){
	  NEW(show_position,MA->glength[t] + 2, BooLean);
	  for(i=j=0; j < MA->length[t]; i++){
	    if(MA->gnull[t][i] != '_'){ j++; show_position[i]=TRUE; }
	    else for(n=1; n<= MA->nseqs; n++){
    	      if(IsInSeqSet(key_seq,MA->id[n])){
		c = MA->gseq[n][t][i];
		if(c != '.') show_position[i] = TRUE;
	      }
	    }
	  }
	}
	for(N=0, n=1; n<= MA->nseqs; n++){
    		if(key_seq && !IsInSeqSet(key_seq,MA->id[n])) continue;
	   	if(t==1){
		   fprintf(fptr,"{\\f2\\fs16 ");
		   fprintf(fptr,"%s",MA->id[n]);
		   fprintf(fptr,"\\tab %d\\tab }",MA->start[n][1]);
	   	}
		if(t < MA->ntyps) {
			gap = MA->start[n][t+1] - MA->end[n][t] - 1;
		}
		if(MA->prob != NULL){
		  if(MA->prob[n][t] < 0.0) strcpy(str,"\\b\\i\\strike");
		  else if(MA->prob[n][t] < 1.3) strcpy(str,"\\b\\i");
		  else strcpy(str,"\\b");
		} else strcpy(str,"\\b");
// new 
		for(laststate=0,i=j=0; j < MA->length[t]; i++){
		   if(MA->gnull[t][i] == '_'){
			if(!show_position || show_position[i]){
			  state='g'; c = MA->gseq[n][t][i];
			  if(state == laststate) fprintf(fptr,"%c",c);
			  else {
			    if(laststate != 0) fprintf(fptr,"}\n");
			    fprintf(fptr,"{%s\\f3\\fs14\\cf16 %c",str,c);  // gap color
			  }
			  laststate = state;
			}
		   } else {
			j++;
			r = MA->seq[n][t][j]; 
			if(r==0 && MA->gseq[n][t][i] == '-') c = '-';
			else c = AlphaChar(r,A);
			if(conserved[j]==NULL || !conserved[j][r]) state='u';
			else state = cons_state[j];
			if(state == laststate) fprintf(fptr,"%c",c);
			else {
			   if(laststate != 0) fprintf(fptr,"}\n");
			   switch(state){
			     case 'H':
				fprintf(fptr,"{%s\\f3\\cf6\\highlight7 %c",str, c);
				 break;  /** red on yellow **/
			     case 'h':
				fprintf(fptr,"{%s\\f3\\cf2\\highlight7 %c",str, c);
				 break;  /** blue on yellow **/
			     case 'a':
				fprintf(fptr,"{%s\\f3\\cf1\\highlight7 %c",str, c);
				 break;  /** black on yellow **/
			     case 'I': case 'S':
				fprintf(fptr,"{%s\\f3\\cf8\\highlight6 %c",str, c);
				 break;  /** white on red **/
			     case 'i': case 's':
				fprintf(fptr,"{%s\\f3\\cf8\\highlight5 %c",str, c);
				 break;  /** white on magenta **/
			     case 'M':
				fprintf(fptr,"{%s\\f3\\cf8\\highlight1 %c",str, c);
				 break;  /** white on black **/
			     case 'm':
				fprintf(fptr,"{%s\\f3\\cf8\\highlight15 %c",str, c);
				 break;  /** white on dark gray **/
			     case 'W':
				  fprintf(fptr,"{%s\\f3\\cf1 %c",str,c);
				break;  /** black on white **/
			     case 'w':
				  fprintf(fptr,"{%s\\f3\\cf15 %c",str,c);
				break;  /** dark gray on white **/
			     case 'u':
			     default : 
				  fprintf(fptr,"{%s\\f3\\cf16 %c",str,c);
				 break;  /** light grey on white **/
			   } 
			}
			laststate = state;
		   }
		}
		fprintf(fptr,"}\n{\\f2\\fs16 ");
		if(t==ntyp) fprintf(fptr,"\\tab %d", MA->end[n][t]);
		else if(gap <=0) fprintf(fptr,"\\tab   ");
		else fprintf(fptr,"\\tab (%d)", gap);
#if 1
		if(MA->prob != NULL) {
			fprintf(fptr,"\\tab [%.1f]",MA->prob[n][t]);
		}
#endif
		fprintf(fptr,"\\tab }{\\f3 \n\\par }");
	}  /** end sequences **/
	if(show_position) free(show_position);
        fprintf(fptr,"\\pard \\widctlpar \\sl180\\slmult1 {\\f3\n\\par\n\\par }\n");
	for(s=1; s <= length; s++)
		if(conserved[s] != NULL) free(conserved[s]);
	free(conserved);
	free(cons_state);
    }
    PutHist(stderr,60,HG);  NilHist(HG);
    free(observed); free(use);
    fprintf(fptr,"\\pard \\widctlpar \\sl180\\slmult1  {\\f3\n\\par}}\n");
}

BooLean	*purge_sma(Int4 cutoff, sma_typ MA)
/**************************************************************************
  Purge from a list of sequences.
/**************************************************************************/
{
	Int4	totlen,i,k,s,t,entry,N,score,max;
	Int4	item,item1,*hits,num_seqs=0;
	b_type	*edges;
	dh_type	H;
	keytyp	key;
	BooLean	*itemList;
	a_type	A=MA->A;
	unsigned char	**seq1,**seq2;

	N=nseqSMA(MA);
	for(totlen=0,t=1; t <= ntypSMA(MA); t++) totlen+=MA->length[t];
	H = dheap(N+2,3);
	NEW(edges,N+2,b_type); NEW(hits,N+2,Int4);
	for(entry =1 ;entry <= N;entry++) edges[entry]=Block(N+1);
	entry = 1;
	for(item =1 ;item <= N;item++) {
	    seq1 = MA->seq[item];
	    for(item1 =item+1 ;item1 <= N;item1++) {
		seq2 =  MA->seq[item1];
		for(score=max=0,t=1; t <= ntypSMA(MA); t++){
                     for(s=1; s <= MA->length[t]; s++){
                        score+=valAlphaR(seq1[t][s],seq2[t][s],A);
                        if(score < 0) score = 0;
                        else { max = MAXIMUM(Int4,score, max); }
                     }
                }
		if(max >= cutoff){
			AddBlock(item,edges[item1]);
			AddBlock(item1,edges[item]);
			hits[item]++; hits[item1]++;
		}
	    } entry++;
	}
	fprintf(stderr,"\n%d items compared; cutoff %d\n", entry-1,cutoff); 
	for(entry=1; entry <= N; entry++){
		key = -(keytyp) hits[entry];
		insrtHeap(entry,key,H);
	}
	while(TRUE){
	    if(emptyHeap(H)) print_error("error in RmHomologs( )");
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

Int4	HspPurgeSMA(Int4 cutoff, sma_typ MA)
/*** find the right cutoff score for purging sequences ***/
{
	Int4		score,n,m,s,t,N,i,j,k,max;
	unsigned char 	**seq1,**seq2,*sq1,*sq2;
	a_type		A=MA->A;
	h_type  	H;
	BooLean		*use;

	/** 1. find the lowest maximum pairwise score sequence **/
	N = nseqSMA(MA);
	H=Histogram("pairwise average scores",0,500,5.0);
  
	for(n = 1; n < N; n++){
	    seq1 =  MA->seq[n]; 
	    for(m = n+1; m <= N; m++){
	        seq2 =  MA->seq[m]; 
		for(score=max=0,t=j=1; t <= ntypSMA(MA); t++){
        	     for(s=1; s <= MA->length[t]; s++){
			score+=valAlphaR(seq1[t][s],seq2[t][s],A);
                	if(score < 0) score = 0;
			else { max = MAXIMUM(Int4,score, max); }
        	     }
		}
		IncdHist(max,H);
	    }
	}
	PutHist(stderr,60,H); NilHist(H);
	use=purge_sma(cutoff, MA);
/*** DEBUG ****/
	H=Histogram("pairwise average scores",0,500,5.0);
/*** DEBUG ****/
	for(i=n=1; n <= N; n++){
	    if(use[n]) {
			fprintf(stderr,"%d: %s\n",i++,MA->id[n]);
/*** DEBUG ****/
	      seq1 =  MA->seq[n]; 
	      for(m = n+1; m <= N; m++){
	       if(use[m]) {
	        seq2 =  MA->seq[m]; 
		for(score=max=0,t=j=1; t <= ntypSMA(MA); t++){
        	     for(s=1; s <= MA->length[t]; s++){
			score+=valAlphaR(seq1[t][s],seq2[t][s],A);
                	if(score < 0) score = 0;
			else { max = MAXIMUM(Int4,score, max); }
        	     }
		}
		IncdHist(max,H);
	       }
	      }
/*** DEBUG ****/
	    }
	}
/*** DEBUG ****/
	PutHist(stderr,60,H); NilHist(H);
/*** DEBUG ****/
	return cutoff;
}

Int4	ClusterSMA(sma_typ MA)
{
	Int4		cutoff,score,n,m,s,totlen,t,N;
	Int4		*max;
	unsigned char 	**seq1,**seq2;
	a_type		A=MA->A;
	h_type  	H;

	/** 1. find the lowest maximum pairwise score sequence **/
	N = nseqSMA(MA);
	H=Histogram("pairwise average scores",-500,500,5.0);
	NEW(max,N+2,Int4);
	for(n = 1; n < N; n++) max[n]=-99999;
	for(totlen=0,t=1; t <= ntypSMA(MA); t++) totlen+=MA->length[t];
	for(n = 1; n < N; n++){
	    seq1 =  MA->seq[n]; 
	    for(m = n+1; m <= N; m++){
	        seq2 =  MA->seq[m]; 
		for(score=0,t=1; t <= ntypSMA(MA); t++){
		     for(s=1; s <= MA->length[t]; s++){
			score+=valAlphaR(seq1[t][s],seq2[t][s],A);
		     }
		}
		score = (100*score)/totlen;
		IncdHist(score,H);
		max[m] = MAXIMUM(Int4,score,max[m]);
		max[n] = MAXIMUM(Int4,score,max[n]);
	    }
	}
	PutHist(stderr,60,H); NilHist(H);
	for(cutoff=+99999, n=1; n <= N; n++) cutoff=MINIMUM(Int4,cutoff,max[n]);
	free(max);
	/** 2. call a recursive procedure to cluster sequences **/
	return cluster_sma(cutoff+2, MA);
}

Int4	cluster_sma(Int4 cutoff, sma_typ MA)
/** cluster the sequences in MA and return as a new MA2 **/
{ 
	Int4		i,j,k,n,m,s,t,N,score,scr,*max,s1,s2,numclust;
	Int4		totlen;
	unsigned char 	**seq1,**seq2;
	a_type		A=MA->A;
	ds_type		sets;
	float		**prob;
	char		**info,**id;
	Int4		**start,**end;
	unsigned char	***seq;
	sma_typ		DUMMY;
	Int4		s0;

fprintf(stderr,"cutoff = %d:\n",cutoff);
	if(nseqSMA(MA) <= 1) return 1;
	N = nseqSMA(MA);
	for(totlen=0,t=1; t <= ntypSMA(MA); t++) totlen+=MA->length[t];
	/*** PutDSets(stdout,sets); /*****/
	sets = DSets(N);
	for(n = 1; n < N; n++){
	    seq1 =  MA->seq[n]; 
	    s1 = findDSets(n,sets);	
	    for(m = n+1; m <= N; m++){
	        seq2 =  MA->seq[m]; 
		for(score=0,t=1; t <= ntypSMA(MA); t++){
		     for(s=1; s <= MA->length[t]; s++){
			score+=valAlphaR(seq1[t][s],seq2[t][s],A);
		     }
		}
		score = (100*score)/totlen;
		if(score >= cutoff){ /** then merge sets **/
	           s2 = findDSets(m,sets);	
		   if(s1 != s2){ s1 = linkDSets(s1,s2,sets); }
		}
	    }
	}
	/*** PutDSets(stdout,sets); /***/
	NEWP(info,N+2,char); NEWP(id,N+2,char);
	NEWP(start,N+2,Int4); NEWP(end,N+2,Int4); 
	NEWP(prob,N+2,float);
	NEWPP(seq,N+2,unsigned char);
	for(s=1; s <= N; s++){
		seq[s] = MA->seq[s];
		info[s] = MA->info[s];
		id[s] = MA->id[s];
		start[s] = MA->start[s];
		end[s] = MA->end[s];
		prob[s] = MA->prob[s];
	}
	numclust=0;
	for(s=0,i=1; i <= N; i++){
	  s1 = findDSets(i,sets);	
	  if(s1 == i){
	     numclust++;
	     s0 = s;
fprintf(stderr,"cluster %d:\n",numclust);
	     for(n=1; n <= N; n++){
		s2 = findDSets(n,sets);	
		if(s2==i){
		    if(seq[n] != NULL){ /** this check shouldn't be necessary..**/
			s++;
			MA->seq[s] = seq[n];
			MA->info[s] = info[n];
			MA->id[s] = id[n];
			MA->start[s] = start[n];
			MA->end[s] = end[n];
			MA->prob[s] = prob[n];
			seq[n]=NULL;
fprintf(stderr,"\t%s\n",id[n]);
		    }
		}	
	     }
	     if((s-s0) > 1){	/*** recursive call ***/
		NEW(DUMMY,1,comultaln_type);
	     	DUMMY->A=A;
		DUMMY->nseqs=(s-s0);
		DUMMY->ntyps=ntypSMA(MA);
		DUMMY->length=MA->length;
		DUMMY->seq = MA->seq + s0;
		DUMMY->id = MA->id + s0;
		DUMMY->info = MA->info + s0;
		DUMMY->start = MA->start + s0;
		DUMMY->end = MA->end + s0;
		DUMMY->prob = MA->prob + s0;
		cluster_sma(cutoff + 2, DUMMY);
		free(DUMMY);
	     }
	  }
	}
	free(seq); free(info);  free(id); 
	free(start); free(end); free(prob); 
	NilDSets(sets);
	return numclust;
}

BooLean *null_col_sma(Int4 t, double cbp_cut, sma_typ MA)
{
        Int4    r,n,s,length,*observed;
        BooLean *conserved,*null;

        length = MA->length[t];
        NEW(null, length+2, BooLean);
        null[1] = FALSE; null[length] = FALSE;
        for(s=2; s < length; s++){
	    NEW(observed,length+2,Int4);
	    for(n=1; n<=MA->nseqs; n++){
		r=MA->seq[n][t][s]; observed[r]++;
	    }
            conserved = conserved_sma(observed,cbp_cut,MA);
            if(conserved == NULL) null[s] = TRUE;
            else null[s] = FALSE;
            free(observed); free(conserved);
        }
        return null;
}

double	*WeightsSMA(sma_typ MA)
{
        Int4    i,j,k,r,s,m,n,totlen,N;
        double  w,max,min,*wt,N2,d;
        short   **nres,*ntyp;
        unsigned char    *seq;
	h_type	HG;
	a_type	A=MA->A;

	N=nseqSMA(MA);
        for(totlen=0,m=1; m<=ntypSMA(MA); m++){ totlen+=MA->length[m]; }

        NEWP(nres,totlen+2,short); NEW(ntyp,totlen+2,short);
        for(i=0; i<=totlen; i++) { NEW(nres[i],nAlpha(A)+2,short); }
        NEW(wt,N+2,double);

        /*** 1. determine the number of residues at each position ***/
        for(s=1; s<=N; s++){
           for(j=1,m=1; m<=ntypSMA(MA); m++){
		seq=MA->seq[s][m];
                for(i=1; i<=MA->length[m]; i++){
                    r=seq[i]; 
		    if(r==0) continue;
		    if(nres[j][r] == 0) { ntyp[j]++; } 
		    nres[j][r]++; j++;
                }
           }
        }

        /*** 2. determine the sequence weights. ***/
        for(max=0.,min=DBL_MAX,s=1; s<=N; s++){
	   w = 0.0; 
           for(j=1,m=1; m<=ntypSMA(MA); m++){
		seq=MA->seq[s][m];
                for(i=1; i<=MA->length[m]; i++){
                    r=seq[i]; 
		    if(r==0) continue;
                    w+= 1.0/ (double) (ntyp[j]*nres[j][r]); j++; 
                }
           }
	   w /= (double) totlen; 
	   wt[s]=w;
	   max = MAXIMUM(double,w,max);
	   min = MINIMUM(double,w,min);
        }
	// E = RandomSeq(totlen, 1, MA->freq, A);
	// HG = Histogram("number residue types", 0, 30,1.0);
        // for(j=1; j<=totlen; j++){ IncdHist(ntyp[j],HG); }
	// PutHist(stdout,60,HG);  NilHist(HG);

        /*** 3. normalize the weights. ***/
	// HG = Histogram("sequence weights", 0, 1,0.02); 
	// HG = Histogram("sequence weights", 0, 100,1.0); 
        for(N2=0.0,s=1; s<=N; s++){
	   wt[s] /= max; 
	   // IncdHist(wt[s],HG);
	   N2 += wt[s];
           // if(wt[s] > 0.9) fprintf(stderr,">%s (%g)\n",infoSMA(s,MA),wt[s]);
	   // fprintf(stderr,"seq %d: weight = %g\n",s,wt[s]); 
        }
	// PutHist(stdout,60,HG);  NilHist(HG);

// fprintf(stderr,"%d seqs: effective number of sequences N = %f\n",N,N2);

        for(i=0; i<=totlen; i++) { free(nres[i]); }
        free(nres); free(ntyp); 
	return wt;
}

BooLean	*FixSMA(double cutoff, sma_typ *NMA, sma_typ OMA)
{
	BooLean	*good;
	Int4	nseqs;

	good = GoodSMA(cutoff, &nseqs, OMA);
	fprintf(stderr,"N = %d; new = %d\n",OMA->nseqs,nseqs);
	*NMA = RmSMA(good, nseqs, OMA);
	return good;
}

sma_typ	RmSMA(BooLean *keep, Int4 nseqs, sma_typ OMA)
{ 
    Int4		m,n,s,t,t0,N,ntyps;          
    Int4 		i,e,length;
    char		r,c,str[200],id[100];
    float   		**prob;         /** seg probability **/
    Int4    		**start,**end;  /** start[n][t]... **/
    sma_typ		MA;
    a_type		A;
    
    NEW(MA,1,comultaln_type);
    MA->A=A=MkAlpha(AMINO_ACIDS,PROT_BLOSUM62);
    NEW(MA->freq,nAlpha(A)+2,double);
    for(i=0; i<=nAlpha(A); i++) MA->freq[i]=(double)blosum62freq[i];

    strcpy(MA->msaid,OMA->msaid);
    strcpy(MA->access,OMA->access);
    strcpy(MA->de,OMA->de);
    MA->lpr = OMA->lpr;
    ntyps = MA->ntyps = OMA->ntyps;
    MA->nseqs = nseqs;
    MA->totcol = OMA->totcol;
    
    NEWP(MA->null,ntyps+2,char);
    NEWP(MA->gnull,ntyps+2,char);
    NEW(MA->ignore,ntyps+2,BooLean);
    NEW(MA->length,ntyps+2,Int4);
    NEW(MA->glength,ntyps+2,Int4);
    NEW(MA->cols,ntyps+2,Int4);
    NEW(MA->fieldmap,ntyps+2,float);
    NEW(MA->single_map,ntyps+2,float);
    NEW(MA->map_minus,ntyps+2,float);

    NEWP(MA->info,nseqs+2,char);
    NEWPP(MA->seq,nseqs+2,unsigned char);
    NEWPP(MA->gseq,nseqs+2,char);
    NEWP(MA->id,nseqs+2,char);
    NEWP(start,nseqs+2,Int4); NEWP(end,nseqs+2,Int4); 
    NEWP(prob,nseqs+2,float);
    MA->start=start; MA->end=end; MA->prob=prob;

    for(t=1; t <= ntyps; t++) {
	MA->ignore[t] = OMA->ignore[t];
	MA->length[t] = OMA->length[t];
	MA->glength[t] = OMA->glength[t];
	MA->cols[t] = OMA->cols[t];
	MA->fieldmap[t] = OMA->fieldmap[t];
	MA->single_map[t] = OMA->single_map[t];
	MA->map_minus[t] = OMA->map_minus[t];
	NEW(MA->null[t], MA->length[t]+2, char);
	NEW(MA->gnull[t], MA->glength[t]+2, char);
	strcpy(MA->null[t],OMA->null[t]);
	strcpy(MA->gnull[t],OMA->gnull[t]);
	MA->ignore[t]=OMA->ignore[t];
    }
    for(m=1,n=0; m <= OMA->nseqs; m++) {
      if(keep[m]){
	n++;
    	NEWP(MA->seq[n],ntyps+2,unsigned char);
    	NEWP(MA->gseq[n],ntyps+2,char);
    	NEW(start[n],ntyps+2,Int4); NEW(end[n],ntyps+2,Int4);
        NEW(prob[n],ntyps+2,float);
	NEW(MA->info[n],strlen(OMA->info[m])+2,char);
	strcpy(MA->info[n],OMA->info[m]);
        for(t=1; t <= ntyps; t++) {
	    NEW(MA->seq[n][t], MA->length[t]+2 , unsigned char);
    	    NEW(MA->gseq[n][t],MA->glength[t]+2,char);
	    MA->start[n][t] = OMA->start[m][t];
	    MA->end[n][t] = OMA->end[m][t];
	    MA->prob[n][t] = OMA->prob[m][t];
	    if(t==1) { 
		NEW(MA->id[n],strlen(OMA->id[m])+2,char); 
	    	strcpy(MA->id[n],OMA->id[m]);
	    }
	    for(s=0; s <= MA->length[t]; s++) MA->seq[n][t][s] = OMA->seq[m][t][s];
	    for(s=0; s <= MA->glength[t]; s++) MA->gseq[n][t][s] = OMA->gseq[m][t][s];
	}
      }
    } assert(n == nseqs);
    return MA;
}

BooLean	*GoodSMA(double cutoff, Int4 *N, sma_typ MA)
/** return a list of all sequences that have good matches to model **/
{ 
    Int4	n,t,N0;
    BooLean	*good;
    
    NEW(good,MA->nseqs+2,BooLean);
    N0 = MA->nseqs;
    for(n = 1; n <= MA->nseqs; n++){
      good[n] = TRUE;
      for(t=1; t <= MA->ntyps; t++) {
	if(MA->prob[n][t] < cutoff){ N0--; good[n]=FALSE; break;}
      }
    }
    *N = N0;
    return good;
}

float	**InfoSMA(Int4 purge_cutoff, sma_typ MA)
// Get the information content at each position in the alignment
{
	Int4	t,i,j,n,s,e,end,N,start,length,*observed,gap,ntyp;
	Int4	h,p,aliphatic,polar,consH,consP,consT,consMax,obs;
	double	info,total,d,f,q;
	a_type	A=MA->A;
	char	r,c,str[50];
	char	*cons_state,state,laststate;  
	BooLean	**conserved,*use;
	double	cbp_cut=-3.0;
	float	**Information;

    ntyp = MA->ntyps;
    use=purge_sma(purge_cutoff, MA);  // try purge cutoff of 250???
    NEW(observed,nAlpha(A)+2,Int4);
    NEWP(Information,ntyp+2,float);
    for(t=1; t <= ntyp; t++){
	length = MA->length[t];
        NEW(Information[t],length+2,float);
	NEWP(conserved, length+2, BooLean);
	NEW(cons_state, length+2, char);
	for(s=1; s <= length; s++){
	    for(r=0; r<=nAlpha(A); r++) observed[r]=0;
	    for(n=1; n<=MA->nseqs; n++){
		if(use[n]){ r=MA->seq[n][t][s]; observed[r]++; }
	    }
	    conserved[s] = conserved_sma(observed,cbp_cut,MA);
	    obs=consMax=consT=consH=consP=aliphatic=polar=0;
	    for(r=1; r<= nAlpha(A); r++){
		n = observed[r];
		sprintf(str,"%c",AlphaChar(r,A));
		if(strstr("AILVFYWMC",str) != NULL) h=n;
		else h=0;
		if(strstr("DEKRHSTQN",str) != NULL) p=n;
		else p = 0;
		if(conserved[s] != NULL && conserved[s][r]){
		    consMax = MAXIMUM(Int4,consMax,n);
		    consT+=n; consH+=h; consP+=p;
		}
		aliphatic+=h; polar+=p; obs+=n;
	    }
#if 1	    // Newcompute relative entropy of positions:
            for(total=0.0,r=1; r <= nAlpha(A); r++){
                total += (double) observed[r]+nAlpha(A)*freqSMA(r,MA);
	    }
            for(info=0.0,r=1; r <= nAlpha(A); r++){
                f = (double) (observed[r] + nAlpha(A)*freqSMA(r,MA))/total;
                if(f > 0.0){
                        q = MA->freq[r];
                        info += f*log(f/q)/log(2.0);
                } 
            } // End relative entropy of postions.
#endif
#if 0	    // Old compute relative entropy of positions:
            for(total=0.0,r=1; r <= nAlpha(A); r++){
                total += (double) observed[r];
	    }
            for(info=0.0,r=1; r <= nAlpha(A); r++){
                f = (double) observed[r]/total;
                if(f > 0.0){
                        q = MA->freq[r];
                        info += f*log(f/q)/log(2.0);
                } 
            } // End relative entropy of postions.
#endif
	    Information[t][s] = info;
	}  /** end sequences **/
	for(s=1; s <= length; s++) if(conserved[s] != NULL) free(conserved[s]);
	free(conserved);
	free(cons_state);
    }
    free(observed);
    free(use);
    return Information;
}

void    NetChargeSMA(FILE *fptr, sma_typ MA)
{
}

float	**ExcessInfoSMA(char *string, Int4 purge_cutoff, sma_typ MA)
{
	Int4	t,i,j,n,s,e,end,N,start,length,*observed,gap,ntyp;
	Int4	other,pos;
	Int4	totpos;
	a_type	A=MA->A;
	char	r,c,str[50];
	BooLean	*use;
	double	info,f,q,total,freqP,freqO;
	float	**Information;
	FILE	*fptr=stdout;

    ntyp = MA->ntyps;
    use=purge_sma(250, MA);
    NEW(observed,nAlpha(A)+2,Int4);
    n=strlen(string);
    for(i=0; i<n; i++){
	if(!isupper(string[i])) print_error("ExcessInfoSMA( ) input error");
    }
    fprintf(fptr,"motif(site):   pos    net  info(tenth bits)\n");
    totpos=0;
    freqP=freqO=0.0;
    for(r=1; r<= nAlpha(A); r++){
	sprintf(str,"%c",AlphaChar(r,A));
	if(strstr(string,str) != NULL) freqP+=MA->freq[r];
	else freqO+=MA->freq[r];
    }
    NEWP(Information,ntyp+2,float);
    for(t=1; t <= ntyp; t++){
	length = MA->length[t];
        NEW(Information[t],length+2,float);
	for(s=1; s <= length; s++){
	    for(r=0; r<=nAlpha(A); r++) observed[r]=0;
	    for(n=1; n<=MA->nseqs; n++){
		if(use[n]){ r=MA->seq[n][t][s]; observed[r]++; }
	    }
	    for(pos=other=0,r=1; r<= nAlpha(A); r++){
		n = observed[r]; sprintf(str,"%c",AlphaChar(r,A));
		if(strstr(string,str) != NULL) pos+=n; else other+=n;
	    }
   /******* compute relative entropy of postions *******/
            for(total=0.0,r=1; r <= nAlpha(A); r++) total += (double) observed[r];
	    info=0.0;
	    f = pos/total; if(f > 0.0) info += f*log(f/freqP)/log(2.0); 
	    f = other/total; if(f > 0.0) info += f*log(f/freqO)/log(2.0); 
	    // Information[t][s] = info;
	    Information[t][s] = 5*info;  // scale up the information...
   /******* compute relative entropy of postions *******/
	    totpos += pos;
	    fprintf(fptr,"%5d(%3d):   %3d    (%.0f)\n", t,s,pos,10*info);
	}
    }
    fprintf(fptr,"total:     %3d\n", totpos);
    free(observed); free(use);
    return Information;
}

