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

#include "alphabet.h"

static float default_blsm62_ins_emit_nats[21] = { 0.0,  // X
//    C       G       A       S       T  
-0.3466, 0.2766,-0.1033, 0.2488, 0.0811,
//    N       D       E       Q       K 
 0.1906, 0.1615, 0.0298, 0.0312, 0.1456,
//    R       H       W       Y       F
 0.0665, 0.0735,-0.2038,-0.1726,-0.2641,
//    V       I       L       M       P
-0.2558,-0.4339,-0.3230,-0.4991, 0.2731};

#if 0
static float default_norm_kyte_doolittle[21] = {0.50, // X
// C    G    A    S    T    N    D    E    Q    K 
 0.78,0.46,0.70,0.40,0.42,0.11,0.11,0.11,0.11,0.07,
// R    H    W    Y    F    V    I    L    M    P
 0.00,0.14,0.40,0.36,0.80,0.97,1.00,0.91,0.71,0.32,};
#endif

static void  default_surface_alpha(a_type A)
{
	char	c;
	Int4	r,r2;
	if(nAlpha(A) == 20){
	  a_type tmpA=DefaultMkAlpha("XCGASTNDEQKRHWYFVILMP",0);
	  for(r=0; r <= nAlpha(A); r++){
		c=AlphaChar(r,A);
		r2 = AlphaCode(c,tmpA);
		// A->surface[r]=default_norm_kyte_doolittle[r2];
		A->surface[r]=default_blsm62_ins_emit_nats[r2];
	  } NilAlpha(tmpA);
	}
}

a_type	MkAlphabet(const char *map_s,const char *R,const char *prs)
/* plus relatedness matrix input as ordered pairs */
{
	a_type	A;
	Int4 i;

	A=MkAlpha(map_s,R);
	if(prs==NULL) return A;
	NEWP(A->pairs,nAlpha(A)+2,char);
	A->matrix = NULL;
	for(i=0;i<=nAlpha(A);i++) NEW(A->pairs[i],nAlpha(A)+2,char);
	DefAlphaP(prs,A);
	return A;
}

BooLean	MemberAlpha(char c, a_type A)
{
	for(Int4 i=0; i<=nAlpha(A); i++) if(AlphaChar(i,A)==c) return TRUE;
	return FALSE;
}

Int4	Str2SeqAlpha(char *str, unsigned char *seq, a_type A)
{
	Int4	i=0;

	seq++;
	while(str[i]!=0){ seq[i] = AlphaCode(str[i],A); i++; }
	return i;
}

a_type	DefAlphaP(const char *prs,a_type A)
{
	Int4 i,j;
	char c,d,n;

	for(i=0;i<=nAlpha(A);i++) {
	   A->paired[i]=FALSE;
	   for(j=0;j<=nAlpha(A);j++) { 
		if(i==j) A->pairs[i][j] = 0; 
		else A->pairs[i][j] = 99; 
	   }
	}
	for(n=1,A->npairs=i=0;isalpha(prs[i]);i+=2){
		if(!isalpha(prs[i+1])) break;
		else {
			c = AlphaCode(prs[i],A);
			d = AlphaCode(prs[i+1],A);
			A->pairs[c][d]=A->pairs[d][c]=n++;
			A->paired[c]=A->paired[d]=TRUE;
			A->npairs++;
		}
	}
	if(A->prs!=NULL) free(A->prs);
	NEW(A->prs,strlen(prs)+3,char);
	for(j=0,i=1;i<=A->npairs;i++) {
		A->prs[j]=prs[j];j++;
		A->prs[j]=prs[j];j++;
	} A->prs[j]=0;
	return A;
}

a_type	CopyAlphabet(a_type A0)
{
	a_type	A;
	Int4 i,j;

	if(A0 == NULL) return NULL;
	A=MkAlpha(A0->alphabet,NULL);
	if(A0->R != NULL) {
	   NEWP(A->R,A->n+2,char);			
	   for(i=0;i<=nAlpha(A);i++){
	     NEW(A->R[i],A->n+2,char);
	     for(j=0;j<=nAlpha(A);j++) A->R[i][j] = A0->R[i][j];
	   }
	   A->loR = A0->loR; A->hiR = A0->hiR;
	}
	A->matrix = NULL;
	if(A0->pairs != NULL) {
	   NEWP(A->pairs,nAlpha(A)+2,char);
	   for(i=0;i<=nAlpha(A);i++) NEW(A->pairs[i],nAlpha(A)+2,char);
	   if(A0->prs!=NULL) DefAlphaP(A0->prs,A);
	}
	return A;
}

a_type  MkAlpha(const char *letters,const char *score_matrix,float *surface)
{
	a_type  A=DefaultMkAlpha(letters,score_matrix);
	for(Int4 r=0; r <= nAlpha(A); r++) A->surface[r]=surface[r];
	return A;
}

a_type  MkAlpha(const char *letters,const char *score_matrix)
{
	a_type  A=DefaultMkAlpha(letters,score_matrix);
	default_surface_alpha(A);
	return A;
}

a_type  DefaultMkAlpha(const char *letters,const char *score_matrix)
{
	a_type	A;
	Int4	n;
	char	c,i;

	IsLicenseExpired();
	NEW(A,1,alphabet_type);
	for(n=0;(c=letters[n])!= '\0';n++) 		/* Get # letters */
		if(!isalpha(c)) alpha_error("Illegal alphabet string",A);
	A->n = n-1; A->N = n;
	NEW(A->alphabet,A->n+3,char);			/* ALPHABET */
	strncpy(A->alphabet,letters,nAlpha(A)+1);	
	NEW(A->code2let,(strlen(letters)+4),char);	/* CODE2LETTER */
	strcpy(A->code2let,letters); 
	strcat(A->code2let,"-"); 			// gap character = nAlpha(A)+1
	NEW(A->code2lower,(strlen(letters)+4),char); 	/* LOWER*/
	strcpy(A->code2lower,letters);
	strcat(A->code2lower,"-"); 			// gap character = nAlpha(A)+1
	NEW(A->let2code,ALPHA_NUM_SYMBOLS,char);	      /* LETTER2CODE */
	for(i=0;i<ALPHA_NUM_SYMBOLS;i++) A->let2code[i]= 0;	/* =error */
	for(i=0;letters[i]!=0;i++) {
		c = letters[i]; A->let2code[c] = i; 
		if(isupper(c)) { c = tolower(c); A->let2code[c] = i; }
		else if(islower(c)) { c = toupper(c); A->let2code[c] = i; }
	}
	A->let2code['-'] = (char) (nAlpha(A)+1);
	for(i=0;letters[i]!=0;i++) 
		if(isupper(letters[i])) A->code2lower[i]=tolower(letters[i]); 
	if(score_matrix==NULL) { A->R = NULL; }			/* RELATION */
	else {
	   NEWP(A->R,A->n+2,char);			
	   for(i=0;i<=nAlpha(A);i++) NEW(A->R[i],A->n+2,char);
	   A = DefAlphaR(score_matrix,A);
	}
	NEW(A->paired,nAlpha(A)+2,BooLean);
	for(i=0;i<=nAlpha(A);i++) A->paired[i]=FALSE;
	A->pairs = NULL;
	A->prs = NULL;
	A->matrix = NULL;
	NEW(A->surface,A->n+2,float);			// surface access...
	return (A);
}

a_type	DefAlphaR(const char *R,a_type A)
/* relatedness matrix */
{
	int	i,j,lo,hi,value;

	if(strlen(R)==0) { alpha_error("Illegal input in DefAlphaR()",A); }
	for(lo=9999,hi=(-9999),i=0;i<=nAlpha(A);i++) {
	   for(j=0;j<=nAlpha(A);j++) {
		while(R[0] != '-' && !isdigit(R[0])) {
		   if(R[0]== 0) alpha_error("Illegal input in DefAlphaR()",A);
		   else R++;
		}
		if(sscanf(R,"%d",&value) != 1)
			alpha_error("Illegal input in DefAlphaR()",A);
		A->R[i][j] = (char) value;
		if(value != (int) A->R[i][j]) /** value too large **/
			alpha_error("Illegal input in DefAlphaR()",A);
		lo = MINIMUM(int,value,lo);
		hi = MAXIMUM(int,value,hi);
		while(R[0] == '-' || isdigit(R[0])) R++;
	   }
	}
	A->loR = lo; A->hiR = hi;
	return (A);
}

Int4	**GetBlastMatrixAlpha(a_type A)
/*** return a blast type scoring matrix (with large negative ends) 
   for use with gapxdrop.c functions. ***/
{
     Int4	**matrix,i,j;

     if(A->R == NULL) alpha_error("scoring matrix undefined", A);
     if(A->matrix == NULL){
        NEWP(matrix,nAlpha(A) + 3, Int4);
        for(i = 0; i <= nAlpha(A); i++){
            NEW(matrix[i],nAlpha(A) + 3, Int4);
            for(j = 0; j <= nAlpha(A); j++){
                matrix[i][j] = valAlphaR(i,j,A);
            }
            matrix[i][j] = SHRT_MIN;
        }
        NEW(matrix[i],nAlpha(A) + 3, Int4);
        for(j = 0; j <= nAlpha(A) + 1; j++) {
            matrix[i][j] = -SHRT_MIN;
        }
	A->matrix = matrix;
     }
     return A->matrix;
}

a_type	PutAlpha(FILE *fptr,a_type A)
{
	Int4	i;

	fprintf(fptr,"Alphabet: %s (n=%d)\n", A->alphabet,nAlpha(A));
	fprintf(fptr," Code:\n");
	for(i=0;i<=nAlpha(A)+1;i++) fprintf(fptr,"  %c",A->code2let[i]);
	fprintf(fptr,"\n");
	for(i=0;i<=nAlpha(A)+1;i++) fprintf(fptr,"%3d",i);
	fprintf(fptr,"\n\n");
	PutAlphaR(fptr,A);
	fprintf(fptr,"\n");
	PutAlphaPairs(fptr,A);
	fprintf(fptr,"\n");
	PutAlphaSurface(fptr,A);
	fprintf(fptr,"\n");
	return A;
}

void	PutAlphaSurface(FILE *fptr,a_type A)
{
	fprintf(fptr, "  surface:\n     ");
	for(Int4 r=1;r<=nAlpha(A);r++){
		fprintf(fptr,"%c  %4.2f    ",AlphaChar(r,A),A->surface[r]);
		if(r%5 == 0) fprintf(fptr,"\n     ");
	} fprintf(fptr,"\n");
}

void	PutAlphaPairs(FILE *fptr, a_type A)
{
	Int4	i,j;
	if(A->pairs != NULL) {
		fprintf(fptr, "  ordered pairs:\n     ");
		for(j=0,i=1; i<=A->npairs; i++){
			fprintf(fptr,"%c", A->prs[j]); j++;
			fprintf(fptr,"%c ", A->prs[j]); j++;
			if(i%10==0 && i < A->npairs){
				fprintf(fptr, "\n     ");
			}
		}
		fprintf(fptr,"\n");
	} else fprintf(fptr,"  ordered pairs:   NULL\n");
	fprintf(fptr,"\n");
}

void	PutAlphaR(FILE *fptr, a_type A)
{
	Int4	i,j;
	if(A->R != NULL) {
		fprintf(fptr, "  relatedness matrix:\n     ");
		for(i=0; i<=nAlpha(A); i++)
			fprintf(fptr,"%3c", AlphaChar(i,A));
		for(i=0; i<=nAlpha(A); i++) {
			fprintf(fptr,"\n   %c:",AlphaChar(i,A));
			for(j=0; j<=nAlpha(A); j++)
				fprintf(fptr,"%3d", valAlphaR(i,j,A));
		}
		fprintf(fptr,"\n");
	} else fprintf(fptr,"  R:   NULL\n");
	fprintf(fptr,"\n");
}


void	NilAlpha(a_type A)
{
	Int4	i;
	if(A != NULL) {
	   free(A->alphabet);
	   free(A->code2let);
	   free(A->code2lower);
	   free(A->let2code);
	   if(A->R != NULL) {
		for(i=0;i<=nAlpha(A);i++) free(A->R[i]);
		free(A->R);
	   }
	   if(A->matrix != NULL) {
		for(i = 0; i <= nAlpha(A); i++) free(A->matrix[i]);
		free(A->matrix[i]); /** free extra array **/
		free(A->matrix); 
	   }
	   if(A->pairs != NULL) {
		for(i=0;i<=nAlpha(A);i++) free(A->pairs[i]);
		free(A->pairs);
	   }
	   if(A->prs!=NULL) free(A->prs);
	   free(A->paired);
	   free(A->surface);
	   free(A);
	}
}

char    *GetPatternFromSST(sst_typ sst,a_type AB)
{
   char tmp[30];
   Int4 r,c=0;
   if(sst == 0) return 0;
   for(r=1; r <= nAlpha(AB); r++){
        if(MemSset(r,sst)){ tmp[c]=AlphaChar(r,AB); c++; }
   } tmp[c]=0;
   return AllocString(tmp);
}


double	*GetJthRhoCategoricalPriors(FILE *fp, double rho, sst_typ *sstJ, a_type AB)
// categorical distribution for rho given # sets...using geometric weighting...1/2 , 1/4, 1/8,...
// weighted geometric distribution to accomodate multiple residue sets.
{
	double *RhoJ;
	assert(sstJ);
	assert(rho > 0.0 && rho <= 0.5);

	assert(sstJ);
	// Count the number of sets with 1, 2, 3, ... etc. residues.
	Int4 x,*NumSets; NEW(NumSets,nAlpha(AB) + 3, Int4);
        for(Int4 s=1; sstJ[s]; s++) NumSets[CardSST(sstJ[s],AB)]++;
	NEW(RhoJ,nAlpha(AB) + 3, double);
	if(sstJ[1] == 0){ RhoJ[0] = 1.0; free(NumSets); return RhoJ; } // insertion at this point in query...
	double p,d,total=0.0;
	for(x=1; x <= nAlpha(AB); x++){
		if(x == 1) d=rho; else d = pow(rho,(double)x);
		if(NumSets[x] > 0){
		   p = d/(double)NumSets[x]; RhoJ[x] = p; total += d;	// total for all sets...
		} else RhoJ[x] = 0;
	} RhoJ[0] = 2.0*RhoJ[1]; total += RhoJ[0];
	for(x=0; x <= nAlpha(AB); x++){ p = RhoJ[x]/total; RhoJ[x] = p; }
	for(total=0,x=0; x <= nAlpha(AB); x++){ if(RhoJ[x] > 0.0) total += RhoJ[x]; }
	assert(total <= 1.0);
        if(fp){
	    fprintf(fp," {}(0,1: %.3g)",RhoJ[0]);
	    total=RhoJ[0];
            for(Int4 s=1; sstJ[s]; s++){
		fprintf(fp,"{");
  		PutSST(fp,sstJ[s],AB);
		Int4 card=CardSST(sstJ[s],AB);
		fprintf(fp,"}(%d,%d: %.3g) ",card,NumSets[card],RhoJ[card]);
		total+=RhoJ[card];
        } 
	fprintf(fp," total:%.3f\n",total);
	} free(NumSets);
	return RhoJ;
}

double	**GetRhoCategoricalPriors(FILE *fp, Int4 Length, double rho, sst_typ **sst, a_type AB)
// categorical distribution for rho given # sets...using geometric weighting...1/2 , 1/4, 1/8,...
// for A_j == null; 
{
	double **Rho;
	assert(Length > 0);
	assert(sst);
        NEWP(Rho,Length+3, double);
	// weighted geometric distribution to accomodate multiple residue sets.
	assert(rho > 0.0 && rho <= 0.5);
        for(Int4 j=1; j <= Length; j++){
	   assert(sst[j]);
	   // Count the number of sets with 1, 2, 3, ... etc. residues.
	   Int4 x,*NumSets; NEW(NumSets,nAlpha(AB) + 3, Int4);
           for(Int4 s=1; sst[j][s]; s++) NumSets[CardSST(sst[j][s],AB)]++;
	   NEW(Rho[j],nAlpha(AB) + 3, double);
	   if(sst[j][1] == 0){
		Rho[j][0] = 1.0; free(NumSets); continue;  // insertion at this point in query...
	   }
	   double p,d,total=0.0;
	   for(x=1; x <= nAlpha(AB); x++){
		if(x == 1) d=rho; else d = pow(rho,(double)x);
		if(NumSets[x] > 0){
		   p = d/(double)NumSets[x];
		   Rho[j][x] = p;
		   total += d;	// total for all sets...
		} else Rho[j][x] = 0;
	   } Rho[j][0] = 2.0*Rho[j][1]; total += Rho[j][0];
	   for(x=0; x <= nAlpha(AB); x++){ p = Rho[j][x]/total; Rho[j][x] = p; }
	   for(total=0,x=0; x <= nAlpha(AB); x++){ if(Rho[j][x] > 0.0) total += Rho[j][x]; }
	   // if(total > 1.0){ fprintf(stderr,"total - 1 = %g\n",total-1.0); }
	   assert(total <= 1.0001);	// good enough; small rounding errors.
           if(fp){
	    fprintf(fp,"%d: {}(0,1: %.3g)",j,log(Rho[j][0]));
	    total=Rho[j][0];
            for(Int4 s=1; sst[j][s]; s++){
		fprintf(fp,"{");
  		PutSST(fp,sst[j][s],AB);
		Int4 card=CardSST(sst[j][s],AB);
		fprintf(fp,"}(%d,%d: %.1f) ",card,NumSets[card],log(Rho[j][card]));
		total+=Rho[j][card];
            } 
	    fprintf(fp," total:%.3f\n",total);
	   } free(NumSets);
	} 
	return Rho;
}

Int4	CardSST(sst_typ sst,a_type AB)
// return cardinality of set sst
{
	Int4	N=0;
        for(Int4 i=1; i <= nAlpha(AB); i++) if(MemSset(i,sst)) N++;
	return N;
}

void    PutSST(FILE *fp,sst_typ sst,a_type AB)
{
        for(Int4 i=1; i <= nAlpha(AB); i++)
          if(MemSset(i,sst)) fprintf(fp,"%c",AlphaChar(i,AB));
}

void    PutSST(char *str,sst_typ sst,a_type AB)
{
	Int4 i,j;
        for(j=0,i=1; i <= nAlpha(AB); i++){
          if(MemSset(i,sst)){ str[j] =AlphaChar(i,AB); j++; }
	} str[j]=0;
}

sst_typ	*StringToSmallSet(char *residue_str,Int4 max_pattern_length,
		const char *Usage,Int4 *pattern_length,a_type AB)
// 
{
	Int4	i,j;
	char	c,state=' ';
	sst_typ	*Sets,tmp_set;
	NEW(Sets,max_pattern_length+3,sst_typ);

	for(j=i=0,state='B'; (c=residue_str[i]) != 0; i++){
	  switch(c){
	     case '{':	// start of a set.
		if(!strchr(" }RB.",state)) print_error(Usage);
		tmp_set=0;	// empty set
		state = '{';
		break;
	     case '}':	// end of a set.
		if(!strchr(" S",state)) print_error(Usage); // must be inside a set.
		Sets[j]=tmp_set; j++;
		state = '}';
		break;
	     case '.':	// Universal set.
		if(!strchr(" }R.",state)) print_error(Usage); // must be inside a set.
		Sets[j]=UnivSset(nAlpha(AB)); j++;
		state = '.';
		break;
	     default:
		if(!isalpha(c)) print_error(Usage);
	  	if(islower(c)) c = toupper(c);
		if(!strchr(" CGASTNDEQKRHWYFVILMP",c)) print_error(Usage);
		switch (state){
		  case '{':
			tmp_set=SsetLet(AlphaCode(c,AB));
			state = 'S'; // now inside of a set.
			break;
		  case 'S': // now inside of a set.
			tmp_set=UnionSset(tmp_set,SsetLet(AlphaCode(c,AB)));
			break;
		  case 'B': // single residue set.
		  case '.': // single residue set.
		  case 'R': // single residue set.
		  case '}': // single residue set.
			Sets[j]=SsetLet(AlphaCode(c,AB)); j++;
			state = 'R'; // just saw a single residue.
			break;
		  default: print_error(Usage);
			break;
		}
		break;
	  } if(j > max_pattern_length) print_error(Usage);
	  printf("%d%c: state = %c\n",j,c, state);
	}
	*pattern_length=j;
#if 0
	for(j=0; j < *pattern_length; j++){
	  fprintf(stderr,"set %d: ",j);
	  PutSST(stderr,Sets[j],AB); fprintf(stderr,"\n");
	}
#endif
	return Sets;
}


void    PutSST2Alpha(FILE *fp,sst_typ sst,a_type AB)
{
        for(Int4 i=0; i <= nAlpha(AB); i++)
          if(MemSset(i,sst)) fprintf(fp,"%c",AlphaChar(i,AB));
}

Int4    CardSstAlpha(sst_typ sst, a_type AB)
{
	Int4	r,N;
	for(N=0,r=1; r <= nAlpha(AB); r++){ if(MemSset(r,sst)) N++; }
	return N;
}

char    *SST2ArgStrAlpha(sst_typ *sst, Int4 len, a_type AB)
// convert a small set type (sst_typ) array to a pattern argumnet string.
// "-P=M11,M14,D15,S17,N25,I28,L31,VI53..."
{
   char	tmp[20],*str;
   Int4 p,i,j,r,n,Size;
   sst_typ X=SsetLet(AlphaCode('X',AB));
   // 1. Determine length of argument string needed
   assert(sst != 0); assert(len < 1000000);
   for(Size=3,i=1; i <= len; i++){
      if(sst[i] == 0 || sst[i] == X) continue;
      Size += CardSstAlpha(sst[i],AB);
      if(i < 10) Size++;
      else if(i < 100) Size += 2; else if(i < 1000) Size += 3;
      else if(i < 10000) Size += 4; else if(i < 100000) Size += 5;
      else Size += 6;
      Size++;	// for ',' or '=' before pattern.
   } // 2. Create the argument string...
   NEW(str,Size+9,char); str[0]='-'; str[1]='P'; str[2]='='; p=3;
   for(n=0,i=1; i <= len; i++){
      if(sst[i] == 0 || sst[i] == X) continue;
      // if(sst[i] == 0) continue;
      if(n > 0){ str[p]=','; p++; }
      for(r=1; r <= nAlpha(AB); r++){
        if(MemSset(r,sst[i])){ str[p]=AlphaChar(r,AB); p++; }
      } sprintf(tmp,"%d",i); 
      for(j=0; tmp[j] != 0; j++){ str[p]=tmp[j]; p++; } n++;
   } str[p]=0; assert(p <= Size);
   return str;
}

sst_typ	*StringToSmallSetsAlpha(char *residue_str,Int4 max_pattern_length,
		const char *Usage,Int4 *pattern_length,a_type AB)
// Convert a character array to small sets of letters in alphabet AB.
// e.g., "{WFYH}.P"
{
	Int4	i,j;
	char	c,state=' ';
	sst_typ	*Sets,tmp_set;
	NEW(Sets,max_pattern_length+3,sst_typ);

	for(j=i=0,state='B'; (c=residue_str[i]) != 0; i++){
	  switch(c){
	     case '{':	// start of a set.
		if(!strchr(" }RB.",state)) print_error(Usage);
		tmp_set=0;	// empty set
		state = '{';
		break;
	     case '}':	// end of a set.
		if(!strchr(" S",state)) print_error(Usage); // must be inside a set.
		Sets[j]=tmp_set; j++;
		state = '}';
		break;
	     case '.':	// Universal set.
		if(!strchr(" }R.",state)) print_error(Usage); // must be inside a set.
		Sets[j]=UnivSset(nAlpha(AB)); j++;
		state = '.';
		break;
	     default:
		if(!isalpha(c)) print_error(Usage);
	  	if(islower(c)) c = toupper(c);
		if(!strchr(" CGASTNDEQKRHWYFVILMP",c)) print_error(Usage);
		switch (state){
		  case '{':
			tmp_set=SsetLet(AlphaCode(c,AB));
			state = 'S'; // now inside of a set.
			break;
		  case 'S': // now inside of a set.
			tmp_set=UnionSset(tmp_set,SsetLet(AlphaCode(c,AB)));
			break;
		  case 'B': // single residue set.
		  case '.': // single residue set.
		  case 'R': // single residue set.
		  case '}': // single residue set.
			Sets[j]=SsetLet(AlphaCode(c,AB)); j++;
			state = 'R'; // just saw a single residue.
			break;
		  default: print_error(Usage);
			break;
		}
		break;
	  } if(j > max_pattern_length) print_error(Usage);
	  printf("%d%c: state = %c\n",j,c, state);
	}
	*pattern_length=j;
#if 0
	for(j=0; j < *pattern_length; j++){
	  fprintf(stderr,"set %d: ",j);
	  PutSST(stderr,Sets[j],AB); fprintf(stderr,"\n");
	}
#endif
	return Sets;
}

void	alpha_error(const char *s, a_type A) 
{ fprintf(stderr,"%s\n",s); PutAlpha(stderr,A); exit(1); }

Int4     ParseResidueSets(char *str, Int4 *pos, sst_typ *sst,a_type AB,const char *msg)
// input string: "YF90,ST97,YF98"
// returns:      3 and sets pos=[3,90,97,98], sst={0,{YF},{ST},{YF}};
{
        Int4    n,v,len,i;
        sst_typ tmp_sst;
        char    tmpstr[200],c;

        if(!isalpha(str[0])) print_error(msg);
        for(n=0; str[0] != 0; ){
           tmp_sst=0;
           if(str[0] == ',') { str++; }
           else if(isalpha(str[0])){
                if(sscanf(str,"%[A-Z]%d",tmpstr,&v) != 2){ print_error(msg); }
                else {
                   n++; pos[n] = v;
                   while(isalpha(str[0])) str++;
                   while(isdigit(str[0])) str++;
                }
                len = strlen(tmpstr);
                for(i=0; i < len; i++){
                   c=tmpstr[i];
                   tmp_sst=UnionSset(tmp_sst,SsetLet(AlphaCode(c,AB)));
                }
                sst[n] = tmp_sst;
           } else print_error(msg);
        } return n;
}


