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

#include "selexCMSA.h"

cma_typ SelexToCMSA(FILE *fp, a_type A)
{
        e_type		*E;
	Int4 		width,*start,*bl_lens,j,k,n,s,t,deflen,c,cons_len,nseq,lenstr;
	char    	*def,str[10000],**toper,*toper_0,*toper_j;
	unsigned char 	*ptr;

	// print_error("This needs to be fixed!!");
	for(cons_len=nseq=0; fgets(str,9995,fp); ){
	   if(isspace(str[0])) continue;	// ignore lines starting with spaces
	   else nseq++;
	   if(nseq==1){
		lenstr=strlen(str);
		width=strlen(str);
// printf("width=%d\n",width); 
		for(c=0; !isspace(str[c]); c++){ assert(c < lenstr); }
		while(isspace(str[c])){ c++; assert(c < lenstr); } deflen=c-1;
	   } else {
		if(strlen(str) != lenstr){
		   printf("lenstr=%d\n",lenstr); 
		   printf("str=%s\n",str); 
		   print_error("Selex file input error 1");
		}
		for(c=0; !isspace(str[c]); c++){ assert(c < lenstr); }
		while(isspace(str[c])){ c++; assert(c < lenstr); } 
		if(c-1 != deflen) print_error("Selex file input error 2");
	   }
        }
	NEW(def,deflen+3,char);
	printf("nseq=%d len=%d; deflen=%d\n",nseq,width,deflen);  
	NEW(ptr,width,unsigned char);
	NEW(E,nseq+2,e_type);
	NEWP(toper,nseq+2,char);
	for(j=0;j<=nseq;j++) { MEW(toper[j],width,char); }
	NEW(toper_0,width,char);
	rewind(fp);

	fgets(str,9995,fp);  
	Int4 end=width-deflen;
	for(k=1;k<=width-deflen;k++) {
		if(isspace(str[deflen+k])){ end=k; break; }
		else if(str[deflen+k] == '-') { toper_0[k] = 'd'; }
		else { toper_0[k] = 'm'; cons_len++; }
		fprintf(stderr,"%d: %c\n",k,str[deflen+k]);
	} printf("cons_len = %d\n",cons_len); 
	rewind(fp);
	for(j=1;j<=nseq;j++){
		toper_j = toper[j]; toper_j[0] = 'E';
		fgets(str,9995,fp);
		// do { fgets(str,9995,fp); } while(isspace(str[0]));
		strncpy(def,str,deflen); def[deflen+1] = '\0';
		n=t=1;
		if(str[deflen+1] != '-') {
			ptr[n++] = AlphaCode(str[deflen+1],A);
			if (toper_0[1] != 'd') toper_j[t++] = 'M';
			else toper_j[t++] = 'I';
		} else { if (toper_0[1] != 'd') toper_j[t++] = 'D'; }
		for(s=deflen+2,k=2; k < end; s++,k++) {
			if(str[s] != '-') {
				ptr[n++] = AlphaCode(str[s],A); 
				if (toper_0[k] != 'd') toper_j[t++] = 'm';
				else toper_j[t++] = 'I';
			} else {
				if (toper_0[k] != 'd')	toper_j[t++] = 'd';
			}
		} toper_j[t] = 'E'; E[j] = MkSeq(def,n-1,ptr);
	}
#if 0
	for(j=1;j<=nseq;j++){ PutSeq(stdout,E[j],A); }
	for(j=1;j<=nseq;j++){ printf("%s\n",toper[j]); }
#endif
	NEW(start,nseq+2,Int4);
	for(j=1; j<=nseq; j++) start[j]=1;
	NEW(bl_lens,2,Int4); bl_lens[1] = cons_len;
	char *name = (char *)"selex2cma";
	return MakeCMSA(E,nseq,toper,start,1,bl_lens,10000,2000,1000.,0,0,name,A);
}

void CMSAToSelex(FILE *ofp, cma_typ cma)
{
	gss_typ		*gss = gssCMSA(cma);
	Int4 		i, j, k, m, p, t, pos[5], nins, cons_len;
	Int4		nseq = NumSeqsCMSA(cma), len = TotalLenCMSA(cma);
	e_type 		E,Etrue,Efake;
	Int4		*ins;
	unsigned char 	*ptr,*ptrf,*ptrt;
	a_type 		A = AlphabetCMSA(cma);
	BooLean		is_a2m=TRUE;

	// print_error("This needs to be checked!!");
	NEW(ins,len+1,Int4);

	E = MkConsensusCMSA(cma);
	cons_len = LenSeq(E); 

//	BooLean IsDeletedCMSA(UInt4 n, UInt4 r, cma_typ cma)
//	returns true if residue r in sequence n is deleted.

//	PosSiteCMSA(nBl,i,pos,cma);

	for(i=1;i<=nseq;i++){
		PosSiteCMSA(1,i,pos,cma);
		p = pos[1];
		for(j=1;j<=cons_len;j++){
			nins = InsertionCMSA(i,p+j-1,cma);		
			if(nins > ins[j]) ins[j] = nins;
		}
	}

	char s2[] = " ";
#if 0
	fprintf(ofp,"%-30s ",strtok(SeqKey(E),s2));

	ptr = SeqPtr(E);
	for(j=1;j<=cons_len;j++){
		fprintf(ofp,"%c",AlphaChar(ptr[j],A));
		for(k=1;k<=ins[j];k++) fprintf(ofp,"-");
	}
	fprintf(ofp,"\n");
#endif

//Int4    gss_typ::TruePos(Int4 s, Int4 site)
// get the position in true sequence corresponding to site in fake seq.
	for(i=1;i<=nseq;i++){
		Efake = gss->FakeSeq(i);
		Etrue = gss->TrueSeq(i);
		ptrf = SeqPtr(Efake);
		ptrt = SeqPtr(Etrue);
//		fprintf(ofp,"%-21.20s",SeqKey(Etrue));		
//	fprintf(ofp,"example ");
//	fprintf(ofp,"%-30s_%d ",strtok(SeqKey(Etrue),s2),i);
	if(is_a2m){ 		// a2m format
		fprintf(ofp,"\n>%s\n",SeqKey(Etrue));
	} else { // my 'selex' format
		fprintf(ofp,"%-30s ",strtok(SeqKey(Etrue),s2));
	}
		PosSiteCMSA(1,i,pos,cma);
		p=pos[1];
		for(k=p,j=1;j<=cons_len;j++){
			if(ptrf[k] == 0) { fprintf(ofp,"-"); k++; }
			else fprintf(ofp,"%c",AlphaChar(ptrf[k++],A));
			nins = InsertionCMSA(i,p+j-1,cma);
			m = gss->TruePos(i,k-1) + 1;
			for(t=1;t<=nins;t++) {
				fprintf(ofp,"%c",AlphaChar(ptrt[m++],A));
			}
			for(t=nins+1;t<=ins[j];t++) {
				fprintf(ofp,"-");
			}
		}
		fprintf(ofp,"\n");
	}
	fprintf(ofp,"\n"); fprintf(ofp,"\n");
	free(ins);
}
