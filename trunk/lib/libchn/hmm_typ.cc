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

#include "hmm_typ.h"

void    hmm_typ::Init(Int4 len, Int4 **matemit, Int4 **insemit, Int4 *mm, Int4 *mi, Int4 *md,
                        Int4 *ii, Int4 *im, Int4 *dd, Int4 *dm, Int4 *bm, Int4 *me,char mode)
{
	const char ch[] = "XACDEFGHIKLMNPQRSTVWY ";
	Int4	i,j,r;
	length=len; 
	InitAsNull();
	for(i=1;i<=len;i++){
		if(mode =='N'){	// already negative..
		   if(mm) m2m[i]=mm[i]; 
		   if(mi) m2i[i]=mi[i]; 
		   if(md) m2d[i]=md[i];
		   if(ii) i2i[i]=ii[i]; 
		   if(im) i2m[i]=im[i]; 
		   if(dd) d2d[i]=dd[i]; 
		   if(dm) d2m[i]=dm[i]; 
		   if(bm) b2m[i]=bm[i]; 
		   if(me) m2e[i]=me[i];
	        } else {	// turn negative...
		   if(mm) m2m[i]=-mm[i]; 
		   if(mi) m2i[i]=-mi[i]; 
		   if(md) m2d[i]=-md[i];
		   if(ii) i2i[i]=-ii[i]; 
		   if(im) i2m[i]=-im[i]; 
		   if(dd) d2d[i]=-dd[i]; 
		   if(dm) d2m[i]=-dm[i]; 
		   if(bm) b2m[i]=-bm[i]; 
		   if(me) m2e[i]=-me[i];
		}
		for(j=1; isalpha(ch[j]); j++){
	   		r= AlphaCode(ch[j],AB);
			mat_emit[i][r]=matemit[i][r];
			if(insemit) ins_emit[i][r]=insemit[i][r];
			else ins_emit[i][r]=DefaultInsEmit[r];
		}
	}
}

void	hmm_typ::Free()
{
}

void	hmm_typ::InitAsNull()
{
	Int4	r,i,len=length;
	const char ch[] = "XACDEFGHIKLMNPQRSTVWY ";
	NEW(m2m,len+3,Int4);NEW(m2i,len+3,Int4);NEW(m2d,len+3,Int4);
	NEW(i2i,len+3,Int4);NEW(i2m,len+3,Int4);NEW(d2d,len+3,Int4);
	NEW(d2m,len+3,Int4);NEW(b2m,len+3,Int4);NEW(m2e,len+3,Int4);	
	NEWP(mat_emit,len+3,Int4);NEWP(ins_emit,len+3,Int4);
        for(i=0;i<=len+1;i++) {
		NEW(mat_emit[i],nAlpha(AB)+2,Int4); NEW(ins_emit[i],nAlpha(AB)+2,Int4); 
		m2m[i]=-2; m2i[i]=-10254; m2d[i]=-894; i2m[i]=-1115; i2i[i]=-701; d2m[i]=-1378; d2d[i]=-42;
		b2m[i] = BigNegative; m2e[i] = BigNegative; 
	}
	NEW(nule,nAlpha(AB)+3,Int4);
	NEW(DefaultInsEmit,nAlpha(AB)+3,Int4);
	for(i=1; isalpha(ch[i]); i++){
	   r= AlphaCode(ch[i],AB);
           switch(ch[i]){
                   case 'A': nule[r] = 595; DefaultInsEmit[r]=-149; break;
                   case 'C': nule[r] = -1558; DefaultInsEmit[r]=-500; break;
                   case 'D': nule[r] = 85; DefaultInsEmit[r]=233; break;
                   case 'E': nule[r] = 338; DefaultInsEmit[r]=43; break;
                   case 'F': nule[r] = -294; DefaultInsEmit[r]=-381; break;
                   case 'G': nule[r] = 453; DefaultInsEmit[r]=399; break;
                   case 'H': nule[r] = -1158; DefaultInsEmit[r]=106; break;
                   case 'I': nule[r] = 197; DefaultInsEmit[r]=-626; break;
                   case 'K': nule[r] = 249; DefaultInsEmit[r]=210;break;
                   case 'L': nule[r] = 902; DefaultInsEmit[r]=-466; break;
                   case 'M': nule[r] = -1085; DefaultInsEmit[r]=-720; break;
                   case 'N': nule[r] = -142; DefaultInsEmit[r]=275; break;
                   case 'P': nule[r] = -21; DefaultInsEmit[r]=394; break;
                   case 'Q': nule[r] = -313; DefaultInsEmit[r]=45; break;
                   case 'R': nule[r] = 45; DefaultInsEmit[r]=96; break;
                   case 'S': nule[r] = 531; DefaultInsEmit[r]=359;break;
                   case 'T': nule[r] = 201; DefaultInsEmit[r]=117;break;
                   case 'V': nule[r] = 384; DefaultInsEmit[r]=-369;break;
                   case 'W': nule[r] = -1998; DefaultInsEmit[r]=-294;break;
                   case 'Y': nule[r] = -644; DefaultInsEmit[r]=-249;break;
                   case 'X': nule[r] = 0; DefaultInsEmit[r]=0;break;
                   default: print_error("alphabet error"); break;
           }
	}
	nb=-8455; nn=-4; ec=-1000; ej=-1000; ct=-8455; cc=-4; jb=-8455; jj=-4; gg=-4; gf=-8455; 
#if 0
	// b2m, m2e; // These are all BigNegative for global models.
	Int4 Bits=1000;
	double	FreqOfMatchAtPositionOne=0.80;	// by default.
	BM = (Int4) floor(Bits*log2(FreqOfMatchAtPositionOne));
	BD=BigNegative;    // for global models... as this approaches zero the models becomes more & more local
#else
	BM=-17;
	BD=-6412;    // for global models... as this approaches zero the models becomes more & more local
#endif
}

hmm_typ::hmm_typ(FILE *fp,a_type A)
{
	Int4	ind;
	char	str[600],*t;
	double U=0.,lambda=0.;
	double SU=0.,Slambda=0.;
	Int4 wthresh,uxparam,uthresh,xparam,gthresh,vthresh;
	Int4 i,j,n,l,len,tmpr,x;
	const char ch[] = "XACDEFGHIKLMNPQRSTVWY";

	AB=A;
	unsigned char uch[42];
	for(n=0;n<=nAlpha(AB);n++) { uch[n] = AlphaCode(ch[n],AB); }
	NEW(name,52,char); 
        do { fgets(str,550,fp); } while (str[0] != 'N');
	t=str;
	while(isalpha(t[0])) t++; while(isspace(t[0]) && t[0] != '\n') t++;
	for(l=0; t[0] != '\n' && l < 50; l++) { name[l] = t[0]; t++; }
	name[l] = 0; // fprintf(stderr,"name = %s\n",name);
        do { fgets(str,550,fp); } while (str[0] != 'D' && str[0] != 'L');
	if (str[0] == 'D') {
		NEW(desc,52,char);
		t=str; while(!isspace(t[0]))t++; while(isspace(t[0]))t++; l = 0;
		while (t[0] != '\n' && l < 50){ desc[l++] = t[0]; t++; }
		desc[l] = 0;
		do fgets(str,550,fp); while (str[0] != 'L');
	} else desc=0;
	t = str; while(!isspace(t[0]))t++; while(isspace(t[0]))t++; sscanf(t,"%d",&len);
	length=len; InitAsNull();

	do { fgets(str,550,fp); } while (str[0] != 'X');
	t = str; while(!isspace(t[0]))t++; while(isspace(t[0]))t++; sscanf(t,"%d",&nb);
	while(!isspace(t[0]))t++; while(isspace(t[0]))t++; sscanf(t,"%d",&nn);
	while(!isspace(t[0]))t++; while(isspace(t[0]))t++; sscanf(t,"%d",&ec);
	while(!isspace(t[0]))t++; while(isspace(t[0]))t++; sscanf(t,"%d",&ej);
	while(!isspace(t[0]))t++; while(isspace(t[0]))t++; sscanf(t,"%d",&ct);
	while(!isspace(t[0]))t++; while(isspace(t[0]))t++; sscanf(t,"%d",&cc);
	while(!isspace(t[0]))t++; while(isspace(t[0]))t++; sscanf(t,"%d",&jb);
	while(!isspace(t[0]))t++; while(isspace(t[0]))t++; sscanf(t,"%d",&jj);
	fgets(str,550,fp); t = str;
	while(!isspace(t[0]))t++; while(isspace(t[0]))t++; sscanf(t,"%d",&gg);
	while(!isspace(t[0]))t++; while(isspace(t[0]))t++; sscanf(t,"%d",&gf);
	fgets(str,550,fp); t = str;
	for(i=1;i<=nAlpha(AB);i++) {
		while(!isspace(t[0]))t++; while(isspace(t[0]))t++; sscanf(t,"%d",&(nule[uch[i]]));
	}
	do { fgets(str,550,fp); } while (str[0] != 'E' && str[0] != 'H');
	if (str[0] == 'E') {
		t = str;
		while(!isspace(t[0]))t++; while(isspace(t[0]))t++; sscanf(t,"%lf",&SU);
		while(!isspace(t[0]))t++; while(isspace(t[0]))t++; sscanf(t,"%lf",&Slambda);
		fgets(str,550,fp);
	}
	if (str[0] == 'E') {
		t = str;
		while(!isspace(t[0]))t++; while(isspace(t[0]))t++; sscanf(t,"%lf",&U);
		while(!isspace(t[0]))t++; while(isspace(t[0]))t++; sscanf(t,"%lf",&lambda);
		fgets(str,550,fp);
	} else { U = 0.; lambda = 0.; }
	if(str[0] == 'W') {
		t = str;
		while(!isspace(t[0]))t++; while(isspace(t[0]))t++; sscanf(t,"%d",&wthresh);
		while(!isspace(t[0]))t++; while(isspace(t[0]))t++; sscanf(t,"%d",&uxparam);
		while(!isspace(t[0]))t++; while(isspace(t[0]))t++; sscanf(t,"%d",&uthresh);
		while(!isspace(t[0]))t++; while(isspace(t[0]))t++; sscanf(t,"%d",&xparam);
		while(!isspace(t[0]))t++; while(isspace(t[0]))t++; sscanf(t,"%d",&gthresh);
		while(!isspace(t[0]))t++; while(isspace(t[0]))t++; sscanf(t,"%d",&vthresh);
	}
	do { fgets(str,550,fp); } while (str[0] != ' ');
	fgets(str,550,fp); t=str;
	while(isspace(t[0]))t++; sscanf(t,"%d",&x); m2m[0] = x;
	while(!isspace(t[0]))t++; while(isspace(t[0]))t++;
	while(!isspace(t[0]))t++; while(isspace(t[0]))t++;
	sscanf(t,"%d",&x); m2d[0] = x;

	for(i=1;i<=len;i++){
		fgets(str,550,fp); t = str;
               	while(isspace(t[0]))t++; while(!isspace(t[0]))t++; while(isspace(t[0]))t++;
               	for(j=1;j<=nAlpha(AB);j++) {
                       	sscanf(t,"%d",&tmpr); mat_emit[i][uch[j]] = tmpr;
                       	while(!isspace(t[0]))t++; while(isspace(t[0]))t++;
                }
		fgets(str,550,fp); t = str;
	        while(isspace(t[0]))t++; while(!isspace(t[0]))t++; while(isspace(t[0]))t++;
	        for(j=1;j<=nAlpha(AB);j++) {
        	        sscanf(t,"%d",&tmpr); ins_emit[i][uch[j]] = tmpr;
                        while(!isspace(t[0]))t++; while(isspace(t[0]))t++;
                } fgets(str,550,fp); t = str;
        	while(isspace(t[0]))t++; while(!isspace(t[0]))t++; while(isspace(t[0]))t++;
		if(t[0] == '*') { m2m[i] = 0; } else { sscanf(t,"%d",&(m2m[i])); }
		while(!isspace(t[0]))t++;while(isspace(t[0]))t++;	
		if(t[0] == '*') { m2i[i] = 0; } else { sscanf(t,"%d",&(m2i[i])); }
		if(m2i[i] == 0) m2i[i] = -1;
	        while(!isspace(t[0]))t++;while(isspace(t[0]))t++;
		if(t[0] == '*') { m2d[i] = 0; } else { sscanf(t,"%d",&(m2d[i])); }
		if(m2d[i] == 0) m2d[i] = -1;
	        while(!isspace(t[0]))t++;while(isspace(t[0]))t++;
		if(t[0] == '*') { i2m[i] = 0; } else { sscanf(t,"%d",&(i2m[i])); }
	        while(!isspace(t[0]))t++;while(isspace(t[0]))t++;
		if(t[0] == '*') { i2i[i] = 0; } else { sscanf(t,"%d",&(i2i[i])); }
		if(i2i[i] == 0) i2i[i] = -1;
	  	while(!isspace(t[0]))t++;while(isspace(t[0]))t++;
		if(t[0] == '*') { d2m[i] = 0; } else { sscanf(t,"%d",&(d2m[i])); }
		while(!isspace(t[0]))t++;while(isspace(t[0]))t++;
		if(t[0] == '*') { d2d[i] = 0; } else { sscanf(t,"%d",&(d2d[i])); }
		if(d2d[i] == 0) d2d[i] = -1;
               	while(!isspace(t[0]))t++;while(isspace(t[0]))t++;
		if(t[0] == '*') { b2m[i] = BigNegative; } else { sscanf(t,"%d",&b2m[i]); }
                while(!isspace(t[0]))t++;while(isspace(t[0]))t++;
		if(t[0] == '*') { m2e[i] = BigNegative; } else { sscanf(t,"%d",&m2e[i]); }
		// b2m & m2e are all really log(prob) not likelihoods! So they always need to be negative!
        }
}

void hmm_typ::Put(FILE *fp)
{
        Int4    i,j,*mat_emit_i,*ins_emit_i;
        Int4    XT[9];
		XT[1] = nb; XT[2] = nn; XT[3] = ec; XT[4] = ej; 
		XT[5] = ct; XT[6] = cc; XT[7] = jb; XT[8] = jj; 
        Int4    NULT[3];
		NULT[1] = gg; NULT[2] = gf;
        Int4    *Nle = nule;
	Int4 	Nule[50];
	const char ch[] = "XACDEFGHIKLMNPQRSTVWY";
        double  EVD[3] = {0.,0.000000,0.000000};
        unsigned char uch[50];
        for(i=0;i<=nAlpha(AB);i++) { uch[i] = AlphaCode(ch[i],AB); }
	for(i=1;i<=nAlpha(AB);i++) Nule[i] = Nle[uch[i]];
        fprintf(fp,"HMMER2.0\n");
        fprintf(fp,"NAME  %s\n",name);
	if(desc) fprintf(fp,"DESC  %s\n",desc);
        fprintf(fp,"LENG  %d\n",length);
        fprintf(fp,"ALPH  Amino\n");
        fprintf(fp,"XT    ");
        for(i=1;i<=8;i++) { fprintf(fp,"%7d",XT[i]); }
        fprintf(fp,"\n");
        fprintf(fp,"NULT ");
        for(i=1;i<=2;i++) { fprintf(fp,"%7d",NULT[i]); }
        fprintf(fp,"\n");
        fprintf(fp,"NULE ");
        for(i=1;i<=20;i++) { fprintf(fp,"%7d",Nule[i]); }
        fprintf(fp,"\n");
        fprintf(fp,"EVD  ");
        for(i=1;i<=2;i++) { fprintf(fp,"%11.6lf",EVD[i]); }
        fprintf(fp,"\n");
        fprintf(fp,"HMM        ");
        for(i=1;i<=20;i++) { fprintf(fp,"%c      ",ch[i]); }
        fprintf(fp,"\n");
        fprintf(fp,"         m->m   m->i   m->d   i->m   i->i   d->m   d->d   b->m   m->e\n");
        fprintf(fp,"         %d      *%7d\n",BM,BD);

        for(i=1;i<=length-1;i++) {
                mat_emit_i = mat_emit[i]; 
                ins_emit_i = ins_emit[i];

                fprintf(fp,"%6d",i);
                for(j=1;j<=20;j++) fprintf(fp,"%7d",mat_emit_i[uch[j]]);
                fprintf(fp,"\n");

                fprintf(fp,"     -");
                for(j=1;j<=20;j++) fprintf(fp,"%7d",ins_emit_i[uch[j]]);
                fprintf(fp,"\n");

                fprintf(fp,"     -");
                fprintf(fp,"%7d",m2m[i]); fprintf(fp,"%7d",m2i[i]); fprintf(fp,"%7d",m2d[i]); 
                fprintf(fp,"%7d",i2m[i]); fprintf(fp,"%7d",i2i[i]); fprintf(fp,"%7d",d2m[i]); 
                fprintf(fp,"%7d",d2d[i]); 
                if(b2m[i]==BigNegative) fprintf(fp,"      *"); else fprintf(fp,"%7d",b2m[i]); 
                if(m2e[i]==BigNegative) fprintf(fp,"      *"); else fprintf(fp,"%7d",m2e[i]); 

                fprintf(fp,"\n");
        }

        fprintf(fp,"%6d",length);
        for(j=1;j<=20;j++) fprintf(fp,"%7d",mat_emit[length][uch[j]]);
        fprintf(fp,"\n");
        fprintf(fp,"     -");
        for(j=1;j<=20;j++) fprintf(fp,"      *");
        fprintf(fp,"\n");
        fprintf(fp,"     -");
        for(j=1;j<=7;j++) fprintf(fp,"      *");
        if(b2m[length]==BigNegative) fprintf(fp,"      *");
	else fprintf(fp,"%7d",b2m[length]);
        fprintf(fp,"      0\n");
        fprintf(fp,"//\n");
}

