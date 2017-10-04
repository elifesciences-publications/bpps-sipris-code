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

#include "table.h"

t_type  Table(Int4 i, sst_typ up, sst_typ down, a_type A)
{
	t_type T; Int4 j,k;

	NEW(T,1,tables_type); 
	T->A = A; T->i[0] = i; 
	T->bin[0][0] = up; T->bin[0][1] = down;
	NEWP(T->cell,2,double); 
	for(j=0;j<2; j++) {
		NEW(T->cell[j],2,double); 
		for(k=0;k<2; k++) T->cell[j][k] = 0.0;
	}
	return T;
}

t_type	AddColumnsTab(Int4 i, sst_typ left, sst_typ right, t_type T)
{
	T->i[1] = i; T->bin[1][0] = left; T->bin[1][1] = right;
	return T;
}

t_type  NilTab(t_type T) 
{ 
	Int4 j;
	for(j=0;j<2; j++) free(T->cell[j]); 
	free(T->cell); free(T); return (t_type)NULL; 
}

t_type  ClearTab(t_type T)
{ Int4 i,j; for(i=0;i<2;i++) for(j=0;j<2;j++) T->cell[i][j]=0.0; return T; }

double	CnTabEntropy(t_type T)
/* returns the symmetrical dependency uxy for the input table */
{
	Int4 i,j;
	double sum=0.0,p,*sumi,*sumj,**nn;
	double h,hx,hy,hygx,hxgy,uygx,uxgy;

	nn=T->cell;
	NEW(sumi,2,double); NEW(sumj,2,double);
	for(i=0;i<2;i++) {
		sumi[i]=0.0;
		for(j=0;j<2;j++){sumi[i] += nn[i][j];sum += (Int4)nn[i][j];}
	}
	for(j=0;j<2;j++) {
		sumj[j]=0.0;
		for(i=0;i<2;i++)sumj[j] += nn[i][j];
	}
	for(hx=0.0,i=0;i<2;i++)
		if (sumi[i]) { p=sumi[i]/sum; hx -= p*log(p); }
	for(hy=0.0,j=0;j<2;j++)
		if (sumj[j]) { p=sumj[j]/sum; hy -= p*log(p); }
	for(h=0.0,i=0;i<2;i++)
		for (j=0;j<2;j++)
			if ((Int4) nn[i][j]) { p=nn[i][j]/sum; h -= p*log(p); }
	hygx=(h)-(hx); hxgy=(h)-(hy);
	uygx=(hy-hygx)/(hy+TABSTINY); uxgy=(hx-hxgy)/(hx+TABSTINY);
	free(sumj); free(sumi);
	return (2.0*(hx+hy-h)/(hx+hy+TABSTINY));
}

/*************** Fisher's Exact Test for a 2x2 contingency table *****
 *  observed:
 *
 *		c0     c1
 *	----+-------+-------+-------
 *	r0  |   a   |   b   | a + b  
 *	----+-------+-------+-------
 *	r1  |   c   |   d   | c + d  
 *	----+-------+-------+-------
 *	    | a + c | b + d |a+b+c+d
 *
 *	see p. 81 of thesis.
 *
 **********************************************************************/
double	ExactTab(t_type T)
{
	double a,b,c,d,i,end,K,P,Q,SUM=0.0;

	a = T->cell[0][0]; b = T->cell[0][1];
	c = T->cell[1][0]; d = T->cell[1][1];
	end = MINIMUM(double,a,d);
	K = lngamma(a+b+1.0) + lngamma(a+c+1.0) 
		+ lngamma(b+d+1.0) + lngamma(c+d+1.0) - lngamma(a+b+c+d+1.0);
	P = lngamma(a+1.0) + lngamma(b+1.0) + lngamma(c+1.0) + lngamma(d+1.0);
	for(i=MAXIMUM(double,-b,-c); i <= end; i = i + 1.0) {
	    Q = lngamma(a-i+1.0) + lngamma(b+i+1.0) 
			+ lngamma(c+i+1.0) + lngamma(d-i+1.0);
	    /* printf(" [%d : %d | %d : %d]\n", (Int4)(a-i), (Int4) (b+i),
				(Int4) (c+i), (Int4) (d-i));/**/
	    if(P <= Q) SUM += exp(K-Q);
	}
	return SUM;
}

Int4	expct_tab(double E[2][2], t_type T)
/*	E[i][j] = M[0][i]*M[1][j]/N*N  */
{
	double	M[2][2],N;
	Int4	i,j;

	for(i=0,N=0.0;i<2;i++) {		/* Get the row totals */
		M[0][i]=0.0;
		for(j=0;j<2;j++) M[0][i] += T->cell[i][j];
	   	N += M[0][i];			/* Get total number */
	}				
	for(j=0;j<2;j++) {
		M[1][j]=0.0;
		for(i=0;i<2;i++)M[1][j]+= T->cell[i][j]; 
	}	/*** Calculate frequencies for each cell */
	if(N)for(i=0;i<2;i++)for(j=0;j<2;j++)E[i][j]=((M[0][i]*M[1][j])/N);
	return	(Int4) N;
}

void    PutTable(FILE *fptr,t_type T)
/* T->cell[r][c] */
{
	Int4	N,i,j,sz,max[5];
	double	prob,E[2][2],depend;
	char	s[3][3][50],str[2][50],name[50];

	N=expct_tab(E,T);
	prob=ExactTab(T);
	depend=CnTabEntropy(T);

	sprintf(s[0][0],"(%d,%d)",T->i[0],T->i[1]);
	get_name_tab(name,1,0,T); sprintf(s[0][1],"%s   ",name);
	get_name_tab(name,1,1,T); sprintf(s[0][2],"%s   ",name);
	for(i=0;i<2;i++) {
	   get_name_tab(name,0,i,T); sprintf(s[i+1][0],"%s   ",name);
	   for(j=0;j<2;j++)
		sprintf(s[i+1][j+1],"%d(%.1f)",(Int4)T->cell[i][j],E[i][j]);
	}
	sprintf(str[0],"%6.2g",prob); sprintf(str[1],"%6.3g",depend); 
	for(i=0;i<3;i++) {
	   for(max[i]=0,j=0;j<3;j++)
		{sz=strlen(s[j][i]);max[i]=MAXIMUM(Int4,max[i],sz);}
	}
	for(i=0;i<2;i++) { sz=strlen(str[i]); max[i+3]=MAXIMUM(Int4,sz,6); }

	fprintf(fptr,"\n%*s   %*s    %*s  ", max[0],s[0][0], max[1],s[0][1], 
			max[2],s[0][2]);
	fprintf(fptr,"  %*s  %*s\n", max[3],"PROB ",max[4],"DEPEND");
	fprintf(fptr,  "%*s [ %*s    %*s ]", max[0],s[1][0], max[1],s[1][1], 
			max[2],s[1][2]);
	fprintf(fptr,"  %*s  %*s\n", max[3],str[0],max[4],str[1]);
	fprintf(fptr,  "%*s [ %*s    %*s ] (total:%d)\n\n",
		 max[0],s[2][0], max[1],s[2][1], max[2],s[2][2],N);
}

Int4	get_name_tab(char *name, Int4 r, Int4 c, t_type T)
{
	Int4 i,j;

	for(j=0,i=1;i<=nAlpha(T->A);i++) {
	   if(MemSset(i,T->bin[r][c])) {
		name[j] = AlphaChar(i,T->A); j++;
	   }
	}
	name[j] = 0;
	return j;
}


//********************** Move below to table.c *************************

Int4    ExpectTable(double E[2][2], t_type T)
/*      E[i][j] = M[0][i]*M[1][j]/N*N  */
// Not used...
{
        double  M[2][2],N;
        Int4    i,j;

        for(i=0,N=0.0;i<2;i++) {                /* Get the row totals */
                M[0][i]=0.0;
                for(j=0;j<2;j++) M[0][i] += T->cell[i][j];
                N += M[0][i];                   /* Get total number */
        }
        for(j=0;j<2;j++) {
                M[1][j]=0.0;
                for(i=0;i<2;i++)M[1][j]+= T->cell[i][j];
        }       /*** Calculate frequencies for each cell */
        if(N)for(i=0;i<2;i++)for(j=0;j<2;j++)E[i][j]=((M[0][i]*M[1][j])/N);
        return  (Int4) N;
}

double	RelEntropyTable(t_type T)
{
	double total=0,No,Ne,I;
	double p,q,Exp[2][2];
	Int4	i,j;
	
	assert(T);
	double **Obs=T->cell;
	Int4 N=expct_tab(Exp,T);
	assert(N > 0);
	for(No=Ne=0.0,i=0;i<2;i++){
	   for(j=0;j<2;j++){ No+=Obs[i][j]; Ne+=Exp[i][j]; }
	}
	for(I=0.0,i=0;i<2;i++){
	   for(j=0;j<2;j++){
		p = Obs[i][j]/No; q = Exp[i][j]/Ne; I+=p*log10(p/q);
	   }
	} return I;
}

double	RelEntropyTableMod(double **Obs, double Exp[2][2])
{
	double total=0,No,Ne,I;
	double p,q;
	Int4	i,j;
	for(No=Ne=0.0,i=0;i<2;i++){
	   for(j=0;j<2;j++){ No+=Obs[i][j]; Ne+=Exp[i][j]; }
	}
	for(I=0.0,i=0;i<2;i++){
	   for(j=0;j<2;j++){
		p = Obs[i][j]/No;
		q = Exp[i][j]/Ne;
		I+=p*log10(p/q);
	   }
	} return I;
}

void    PutTableMod(FILE *fptr,t_type T) { PutTableMod(fptr,0,T); }

void    PutTableMod(FILE *fptr,char *Title,t_type T)
/* T->cell[r][c] */
{
	Int4	N,i,j,sz,max[5];
	double	prob,E[2][2],depend;
	char	s[3][3][50],str[2][50],name[50],Set[2][50];

	N=expct_tab(E,T);	
	if(N==0) E[0][0]=E[0][1]=E[1][0]=E[1][1]=0;
	prob=ExactTab(T);
	depend=CnTabEntropy(T);

	sprintf(s[0][0],"(%d,%d)",T->i[0],T->i[1]);
	get_name_tab(Set[0],1,0,T); get_name_tab(Set[1],1,1,T);
	// sprintf(s[0][1],"%s   ",name); sprintf(s[0][2],"%s   ",name);
	sprintf(s[0][1],"%s   ",Set[0]); 
	sprintf(s[0][2],"^%s  ",Set[0]); 
	sprintf(s[1][0]," %s  ",Set[1]);
	sprintf(s[2][0],"^%s  ",Set[1]);
	for(i=0;i<2;i++) {
	   // get_name_tab(name,0,i,T); sprintf(s[i+1][0],"%s   ",name);
	   for(j=0;j<2;j++)
		sprintf(s[i+1][j+1],"%d(%.1f)",(Int4)T->cell[i][j],E[i][j]);
	}
	sprintf(str[0],"%6.2g",prob); sprintf(str[1],"%6.3g",depend); 
	for(i=0;i<3;i++) {
	   for(max[i]=0,j=0;j<3;j++)
		{sz=strlen(s[j][i]);max[i]=MAXIMUM(Int4,max[i],sz);}
	}
	for(i=0;i<2;i++) { sz=strlen(str[i]); max[i+3]=MAXIMUM(Int4,sz,6); }

	if(Title) fprintf(fptr,"\n%s\n",Title); else fprintf(fptr,"\n");
	fprintf(fptr,"%*s   %*s    %*s  ", max[0],s[0][0], max[1],s[0][1], 
			max[2],s[0][2]);
	fprintf(fptr,"  %*s  %*s\n", max[3],"PROB ",max[4],"DEPEND");
	fprintf(fptr,  "%*s [ %*s    %*s ]", max[0],s[1][0], max[1],s[1][1], 
			max[2],s[1][2]);
	fprintf(fptr,"  %*s  %*s\n", max[3],str[0],max[4],str[1]);
	fprintf(fptr,  "%*s [ %*s    %*s ] (total:%d)\n",
		 max[0],s[2][0], max[1],s[2][1], max[2],s[2][2],N);
#if 0
	double I=RelEntropyTableMod(T->cell,E);
	fprintf(fptr,  "Relative Entropy = %.3f\n",I);
#endif
}

#if 0	// OLD resurrected CODE
/**********************************************************************
CnTable() - Given a 2-dimensional contingency table in the form 
of a double array bin[1..ni][1..nj], this routine returns the 
chi-square 'chsq', the number of degrees of freedom 'df' and the 
significance level 'prob' (small values indicating a significant association).
    Assuming no associations between i and j the expected value 
can be derived by the formula for a two dimensional table...

        N(..) = N(i.)*N(.j)/N....                 for 2 dimensions      (1)

    where 
                ---                       ---
                \                         \
        N(.j) = /   N(ij),        N(i.) = /    N(ij), 
                ---                       ---
                 i                         j

    The number of degrees of freedom is equal to the number of entries
in the table - 1 (e.g., dim=2: ni*nj=4-1) minus the number of probabilities
estimated from the data for the particular hypothesis being tested.
Since we are here considering the hypothesis of mutual independence
there are (ni-1)+(nj-1) probabilites estimated from the
data.  Therefore 
        df = (ni*nj) - (ni+nj) + 1.  
**********************************************************************/
double  CnTable(t_type T,double *chsq, double *df)
{
        int     ni,nj,i,j,p,sg,sh;
        double  expctd,denom=1.0,N=0.0, *Ni,*Nj,temp;
        double  small=0.0,total=0.0,**bin;
        boolean okay=TRUE;

        bin = T->bin; ni=2; nj=2;
        NEW(Ni,(ni),double);            /* = Ni. values */
        NEW(Nj,(nj),double);            /* = N.j values */
        for(i=0;i<ni && okay;i++) {     /* Get the row totals */
                Ni[i]=0.0;
                for(j=0;j<nj;j++) Ni[i] += bin[i][j];
                N += Ni[i];                     /* Get total number */
                if(Ni[i] == 0.0) okay=FALSE;
        }
        denom *= N;
        for(j=0;j<nj && okay;j++) {
                Nj[j]=0.0;
                for(i=0;i<ni;i++)Nj[j]+=bin[i][j];
                if(Nj[j] == 0.0) okay=FALSE;
        }
        *chsq=0.0;
        for(i=0;i<ni;i++) {
                for(j=0;j<nj;j++) {
                        expctd=(Ni[i]*Nj[j])/(denom);
                        /****** TEST VALIDITY OF CHI SQUARED VALUE ******
                        total += 1.0;
                        if(expctd < 2.0) { okay = FALSE; break; }
                        if (expctd < 5.0) small += 1.0;
                        /****************** END OF TEST *****************/
                        temp=bin[i][j]-expctd;
                        *chsq += temp*temp/(expctd+TABSTINY); /* Here TINY */
                }       /* guarantees that any eliminated row or column */
        }               /* will not contribute to the N. */
        /****** TEST VALIDITY OF CHI SQUARED VALUE ******
        if(okay && small/total > 0.20) okay = FALSE;
        /****************** END OF TEST *****************/
        if(okay) {
           *df=(double)(ni*nj-(ni+nj)+1);
           free(Ni); free(Nj); 
           return gammq((double)(0.5*(*df)),(double)(0.5*(*chsq)));
        } else { free(Ni);free(Nj);return ILLEGAL; }
}
#endif

