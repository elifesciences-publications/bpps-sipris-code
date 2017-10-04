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

#include "shm_typ.h"

static Int4 **Concat(Int4 **smth, Int4 smth_len, Int4 smth_width, Int4 rpts);
static Int4 Sample2(double a1);
static Int4 Sample3(double a1, double a2);
static Int4 Sample(double *a,Int4 size_a);

static Int4 **Concat(Int4 **smth, Int4 smth_len, Int4 smth_width, Int4 rpts)  {
        Int4 **rpt_smth, *rpt_smth_tmpr, *smth_i; 
        Int4 i,j,k; 
        NEWP(rpt_smth,rpts*smth_len+4,Int4); 
        for(i=0;i<=rpts*smth_len;i++) NEW(rpt_smth[i],smth_width+3,Int4); 
        for(j=0;j<rpts;j++) {
                for(i=1;i<=smth_len;i++) {
                        rpt_smth_tmpr = rpt_smth[j*smth_len+i]; 
                        smth_i = smth[i]; 
                        for(k=0;k<=smth_width;k++) {
                                rpt_smth_tmpr[k] = smth_i[k]; 
                        }
                }
        }
        for(i=0;i<=smth_len;i++) free(smth[i]); 
                free(smth); 
        return rpt_smth; 
}   

static Int4 Sample2(double a1)  {
        double rand_num; 
        rand_num=SampleUniformProb(); 
        if((rand_num -= a1) <= 0.0) return 1; 
        else return 2; 
}
 
static Int4 Sample3(double a1, double a2)  {
        double rand_num; 
        rand_num=SampleUniformProb(); 
        if((rand_num -= a1) <= 0.0) return 1; 
        else if((rand_num -= a2) <= 0.0) return 2; 
        else return 3;  
}

static Int4 Sample4(double a1, double a2, double a3)  {
        double rand_num; 
        rand_num=SampleUniformProb(); 
        if((rand_num -= a1) <= 0.0) return 1; 
        else if((rand_num -= a2) <= 0.0) return 2; 
        else if((rand_num -= a3) <= 0.0) return 3; 
        else return 4;  
}
 
static Int4 Sample(double *a,Int4 size_a)  {
        double rand_num; 
        rand_num=SampleUniformProb(); 
        for(Int4 i=1;i<=size_a-1;i++){
                if((rand_num -= a[i]) <= 0.0) return i; 
        }
	return size_a;
}                    

shm_typ **ReadShm(FILE *fp, Int4 bits, Int4 *nmbr, a_type A)
{
	shm_typ **shm;
	Int4 br=0,ind;
	char str[600],*t;

	while(fgets(str,550,fp) != NULL){
		if(str[0] == '/') br++;
	}
	*nmbr = br;
	rewind(fp);
	stpb_typ *stpb;
	double U=0.,lambda=0.;
	double SU=0.,Slambda=0.;
	Int4 wthresh,uxparam,uthresh,xparam,gthresh,vthresh;
	Int4 **mat_emit,**ins_emit;
	Int4 *nule;
	char *name,*desc;
	Int4 i,j,n,l,len;
	Int4 *m2m,*m2i,*m2d,*i2i,*i2m,*d2d,*d2m,*b2m,*m2e;
	Int4 nb,nn,ec,ej,ct,cc,jb,jj;
	Int4 bm,bd;
	Int4 gg,gf;

	NEWP(shm,br+3,shm_typ);

	Int4 tmpr;

	char ch[] = "XACDEFGHIKLMNPQRSTVWY";
	unsigned char uch[22];
	for(n=0;n<=nAlpha(A);n++) { uch[n] = AlphaCode(ch[n],A); }
	for(n=1;n<=br;n++){
		NEW(name,3*sizeof(double),char);
		NEW(desc,6*sizeof(double),char);
                do { fgets(str,550,fp); }
                while (str[0] != 'N');
		t=str;
                while(!isspace(t[0]))t++;
                while(isspace(t[0]))t++;
		l = 0;
		while (t[0] != '\n' && l < 3*sizeof(double)-2) { name[l++] = t[0]; t++; }
		name[l++] = ' ';
		name[l] = '\0';
                do fgets(str,550,fp);
                while (str[0] != 'D' && str[0] != 'L');
		if (str[0] == 'D') {
			t=str;
			while(!isspace(t[0]))t++;
			while(isspace(t[0]))t++;
			l = 0;
			while (t[0] != '\n' && l < 6*sizeof(double)-1){
				desc[l++] = t[0]; t++; 
			} desc[l] = '\0';
			do fgets(str,550,fp); 
			while (str[0] != 'L');
		}
		t = str;
		while(!isspace(t[0]))t++;
		while(isspace(t[0]))t++;
		sscanf(t,"%d",&len);
		NEW(m2m,len+3,Int4);NEW(m2i,len+3,Int4);NEW(m2d,len+3,Int4);
		NEW(i2i,len+3,Int4);NEW(i2m,len+3,Int4);NEW(d2d,len+3,Int4);
		NEW(d2m,len+3,Int4);NEW(b2m,len+3,Int4);NEW(m2e,len+3,Int4);	

	        NEWP(mat_emit,len+3,Int4);NEWP(ins_emit,len+3,Int4);
        	for(i=0;i<=len+1;i++) { 
			NEW(mat_emit[i],nAlpha(A)+2,Int4); 
			NEW(ins_emit[i],nAlpha(A)+2,Int4);
		}

	        do { fgets(str,550,fp); }
		while (str[0] != 'X');
		t = str;
		while(!isspace(t[0]))t++;
		while(isspace(t[0]))t++;
		sscanf(t,"%d",&nb);

		while(!isspace(t[0]))t++;
		while(isspace(t[0]))t++;
		sscanf(t,"%d",&nn);

		while(!isspace(t[0]))t++;
		while(isspace(t[0]))t++;
		sscanf(t,"%d",&ec);

		while(!isspace(t[0]))t++;
		while(isspace(t[0]))t++;
		sscanf(t,"%d",&ej);

		while(!isspace(t[0]))t++;
		while(isspace(t[0]))t++;
		sscanf(t,"%d",&ct);

		while(!isspace(t[0]))t++;
		while(isspace(t[0]))t++;
		sscanf(t,"%d",&cc);

		while(!isspace(t[0]))t++;
		while(isspace(t[0]))t++;
		sscanf(t,"%d",&jb);

		while(!isspace(t[0]))t++;
		while(isspace(t[0]))t++;
		sscanf(t,"%d",&jj);

		fgets(str,550,fp);
		t = str;
		while(!isspace(t[0]))t++;
		while(isspace(t[0]))t++;
		sscanf(t,"%d",&gg);

		while(!isspace(t[0]))t++;
		while(isspace(t[0]))t++;
		sscanf(t,"%d",&gf);

		NEW(nule,nAlpha(A)+3,Int4);
		fgets(str,550,fp);
		t = str;
		for(i=1;i<=nAlpha(A);i++) {
			while(!isspace(t[0]))t++;
			while(isspace(t[0]))t++;
			sscanf(t,"%d",&(nule[uch[i]]));
		}

	        do { fgets(str,550,fp); }
		while (str[0] != 'E' && str[0] != 'H');


		if (str[0] == 'E') {
			t = str;
			while(!isspace(t[0]))t++;
			while(isspace(t[0]))t++;
			sscanf(t,"%lf",&SU);
			while(!isspace(t[0]))t++;
			while(isspace(t[0]))t++;
			sscanf(t,"%lf",&Slambda);
			fgets(str,550,fp);
		}
		if (str[0] == 'E') {
			t = str;
			while(!isspace(t[0]))t++;
			while(isspace(t[0]))t++;
			sscanf(t,"%lf",&U);
			while(!isspace(t[0]))t++;
			while(isspace(t[0]))t++;
			sscanf(t,"%lf",&lambda);
			fgets(str,550,fp);
		} else { U = 0.; lambda = 0.; }
		if(str[0] == 'W') {
			t = str;
			while(!isspace(t[0]))t++;
			while(isspace(t[0]))t++;
			sscanf(t,"%d",&wthresh);
			while(!isspace(t[0]))t++;
			while(isspace(t[0]))t++;
			sscanf(t,"%d",&uxparam);
			while(!isspace(t[0]))t++;
			while(isspace(t[0]))t++;
			sscanf(t,"%d",&uthresh);
			while(!isspace(t[0]))t++;
			while(isspace(t[0]))t++;
			sscanf(t,"%d",&xparam);
			while(!isspace(t[0]))t++;
			while(isspace(t[0]))t++;
			sscanf(t,"%d",&gthresh);
			while(!isspace(t[0]))t++;
			while(isspace(t[0]))t++;
			sscanf(t,"%d",&vthresh);
		}
	        do { fgets(str,550,fp); }
		while (str[0] != ' ');
		fgets(str,550,fp);

		t=str;
		while(isspace(t[0]))t++;
		sscanf(t,"%d",&bm);
		m2m[0] = bm;
		while(!isspace(t[0]))t++;while(isspace(t[0]))t++;
		while(!isspace(t[0]))t++;while(isspace(t[0]))t++;
		sscanf(t,"%d",&bd);
		m2d[0] = bd;

	        for(i=1;i<=len;i++){
			fgets(str,550,fp);
                	t = str;
                	while(isspace(t[0]))t++;
                	while(!isspace(t[0]))t++;
                	while(isspace(t[0]))t++;
                	for(j=1;j<=nAlpha(A);j++) {
                        	sscanf(t,"%d",&tmpr);
                        	mat_emit[i][uch[j]] = tmpr;
                        	while(!isspace(t[0]))t++;
                        	while(isspace(t[0]))t++;
                	}

			fgets(str,550,fp);
			t = str;
	                while(isspace(t[0]))t++;
        	        while(!isspace(t[0]))t++;
                	while(isspace(t[0]))t++;
	                for(j=1;j<=nAlpha(A);j++) {
        	                sscanf(t,"%d",&tmpr);
                	        ins_emit[i][uch[j]] = tmpr;
                        	while(!isspace(t[0]))t++;
                	        while(isspace(t[0]))t++;
                	}

			fgets(str,550,fp);
	                t = str;
        	        while(isspace(t[0]))t++;
	                while(!isspace(t[0]))t++;
	                while(isspace(t[0]))t++;

			if(t[0] == '*') { m2m[i] = 0; }
			else {
				sscanf(t,"%d",&(m2m[i]));
			}
			while(!isspace(t[0]))t++;while(isspace(t[0]))t++;	

			if(t[0] == '*') { m2i[i] = 0; }
			else {
		                sscanf(t,"%d",&(m2i[i]));
			}
			if(m2i[i] == 0) m2i[i] = -1;
	                while(!isspace(t[0]))t++;while(isspace(t[0]))t++;

			if(t[0] == '*') { m2d[i] = 0; }
			else {
		                sscanf(t,"%d",&(m2d[i]));
			}
			if(m2d[i] == 0) m2d[i] = -1;
	                while(!isspace(t[0]))t++;while(isspace(t[0]))t++;

			if(t[0] == '*') { i2m[i] = 0; }
			else {
	        	        sscanf(t,"%d",&(i2m[i]));
			}
	                while(!isspace(t[0]))t++;while(isspace(t[0]))t++;

			if(t[0] == '*') { i2i[i] = 0; }
			else {
	        	        sscanf(t,"%d",&(i2i[i]));
      			}
			if(i2i[i] == 0) i2i[i] = -1;
	  	        while(!isspace(t[0]))t++;while(isspace(t[0]))t++;

			if(t[0] == '*') { d2m[i] = 0; }
			else {
	                	sscanf(t,"%d",&(d2m[i]));
			}
			while(!isspace(t[0]))t++;while(isspace(t[0]))t++;

			if(t[0] == '*') { d2d[i] = 0; }
			else {
		                sscanf(t,"%d",&(d2d[i]));
			}
			if(d2d[i] == 0) d2d[i] = -1;
                	while(!isspace(t[0]))t++;while(isspace(t[0]))t++;

			if(t[0] == '*') { b2m[i] = BigNegative; }
			else {
	                	sscanf(t,"%d",&b2m[i]);
			}
                	while(!isspace(t[0]))t++;while(isspace(t[0]))t++;

			if(t[0] == '*') { m2e[i] = BigNegative; }
			else {
	                	sscanf(t,"%d",&m2e[i]);
			}
        	}

		stpb = new stpb_typ(len,m2m,m2i,m2d,i2i,i2m,d2d,d2m,b2m,m2e,nb,nn,ec,ej,
			ct,cc,jb,jj,bm,bd,gg,gf,bits);

		double x = 0.;
		double *nl;
		NEW(nl,nAlpha(A)+2,double);
		for(j=1;j<=nAlpha(A);j++) nl[j] = BitsToProb(nule[j],1./nAlpha(A),bits);

		for(i=1;i<=len;i++) {
			for(x=0.,j=1;j<=nAlpha(A);j++) {
				x += nl[j]*mat_emit[i][j];
			}
			mat_emit[i][0] = (Int4) floor(0.5 + x);

			for(x=0.,j=1;j<=nAlpha(A);j++) {
				x += nl[j]*ins_emit[i][j];
			}
			ins_emit[i][0] = (Int4) floor(0.5 + x);
		}		

		shm[n] = new shm_typ(name,desc,nule,mat_emit,ins_emit,stpb,A);
		shm[n]->SetLambda(lambda);
		shm[n]->SetU(U);
		shm[n]->SetSLambda(Slambda);
		shm[n]->SetSU(SU);
		shm[n]->SetWthreshold(wthresh);
		shm[n]->SetUxparameter(uxparam);
		shm[n]->SetUthreshold(uthresh);
		shm[n]->SetXparameter(xparam);
		shm[n]->SetGthreshold(gthresh);
		shm[n]->SetVthreshold(vthresh);
		fgets(str,550,fp);
	}

	return shm;
}

shm_typ *ReadProbShm(FILE *fp, Int4 bits, a_type A)
{
	shm_typ *shm;
	char str[600],*t;

	stpb_typ *stpb;
	double U=0.,lambda=0.;
	double SU=0.,Slambda=0.;
	double **mat_emit,**ins_emit;
	double *nule;
	char *name,*desc;
	Int4 i,j,n,l,len;
	double *m2m,*m2i,*m2d,*i2i,*i2m,*d2d,*d2m,*b2m,*m2e;
	double nb,nn,ec,ej,ct,cc,jb,jj;
	double bm,bd;
	double gg,gf;

	double tmpr;

	char ch[] = "XACDEFGHIKLMNPQRSTVWY";
	unsigned char uch[22];
	for(n=0;n<=nAlpha(A);n++) { uch[n] = AlphaCode(ch[n],A); }
	NEW(name,3*sizeof(double),char);
	NEW(desc,6*sizeof(double),char);
        do { fgets(str,550,fp); }
        while (str[0] != 'N');
	t=str;
        while(!isspace(t[0]))t++;
        while(isspace(t[0]))t++;
	l = 0;
	while (t[0] != '\n' && l < 3*sizeof(double)-2) { name[l++] = t[0]; t++; }
	name[l++] = ' ';
	name[l] = '\0';
        do fgets(str,550,fp);
        while (str[0] != 'D' && str[0] != 'L');
	if (str[0] == 'D') {
		t=str;
		while(!isspace(t[0]))t++;
		while(isspace(t[0]))t++;
		l = 0;
		while (t[0] != '\n' && l < 6*sizeof(double)-1) { desc[l++] = t[0]; t++; }
		desc[l] = '\0';
		do fgets(str,550,fp); 
		while (str[0] != 'L');
	}
	t = str;
	while(!isspace(t[0]))t++;
	while(isspace(t[0]))t++;
	sscanf(t,"%d",&len);
	NEW(m2m,len+3,double);NEW(m2i,len+3,double);NEW(m2d,len+3,double);
	NEW(i2i,len+3,double);NEW(i2m,len+3,double);NEW(d2d,len+3,double);
	NEW(d2m,len+3,double);NEW(b2m,len+3,double);NEW(m2e,len+3,double);	
        NEWP(mat_emit,len+2,double);NEWP(ins_emit,len+2,double);
        for(i=0;i<=len+1;i++) { 
		NEW(mat_emit[i],nAlpha(A)+2,double); 
		NEW(ins_emit[i],nAlpha(A)+2,double);
	}

	do { fgets(str,550,fp); }
	while (str[0] != 'X');
	t = str;
	while(!isspace(t[0]))t++;
	while(isspace(t[0]))t++;
	sscanf(t,"%lf",&nb);

	while(!isspace(t[0]))t++;
	while(isspace(t[0]))t++;
	sscanf(t,"%lf",&nn);

	while(!isspace(t[0]))t++;
	while(isspace(t[0]))t++;
	sscanf(t,"%lf",&ec);

	while(!isspace(t[0]))t++;
	while(isspace(t[0]))t++;
	sscanf(t,"%lf",&ej);

	while(!isspace(t[0]))t++;
	while(isspace(t[0]))t++;
	sscanf(t,"%lf",&ct);

	while(!isspace(t[0]))t++;
	while(isspace(t[0]))t++;
	sscanf(t,"%lf",&cc);

	while(!isspace(t[0]))t++;
	while(isspace(t[0]))t++;
	sscanf(t,"%lf",&jb);

	while(!isspace(t[0]))t++;
	while(isspace(t[0]))t++;
	sscanf(t,"%lf",&jj);

	fgets(str,550,fp);
	t = str;
	while(!isspace(t[0]))t++;
	while(isspace(t[0]))t++;
	sscanf(t,"%lf",&gg);

	while(!isspace(t[0]))t++;
	while(isspace(t[0]))t++;
	sscanf(t,"%lf",&gf);

	NEW(nule,nAlpha(A)+3,double);
	fgets(str,550,fp);
	t = str;
	for(i=1;i<=nAlpha(A);i++) {
		while(!isspace(t[0]))t++;
		while(isspace(t[0]))t++;
		sscanf(t,"%lf",&(nule[uch[i]]));
	}

        do { fgets(str,550,fp); }
	while (str[0] != 'E' && str[0] != 'H');


	if (str[0] == 'E') {
		t = str;
		while(!isspace(t[0]))t++;
		while(isspace(t[0]))t++;
		sscanf(t,"%lf",&SU);
		while(!isspace(t[0]))t++;
		while(isspace(t[0]))t++;
		sscanf(t,"%lf",&Slambda);
		fgets(str,550,fp);
	}
	if (str[0] == 'E') {
		t = str;
		while(!isspace(t[0]))t++;
		while(isspace(t[0]))t++;
		sscanf(t,"%lf",&U);
		while(!isspace(t[0]))t++;
		while(isspace(t[0]))t++;
		sscanf(t,"%lf",&lambda);
		fgets(str,550,fp);
	} else { U = SU; lambda = Slambda; SU = 0.; Slambda = 0.; }
        do { fgets(str,550,fp); }
	while (str[0] != ' ');
	fgets(str,550,fp);

	t=str;
	while(isspace(t[0]))t++;
	sscanf(t,"%lf",&bm);
	while(!isspace(t[0]))t++;while(isspace(t[0]))t++;
	while(!isspace(t[0]))t++;while(isspace(t[0]))t++;
	sscanf(t,"%lf",&bd);

        for(i=1;i<=len;i++){
		fgets(str,550,fp);
               	t = str;
               	while(isspace(t[0]))t++;
               	while(!isspace(t[0]))t++;
               	while(isspace(t[0]))t++;
               	for(j=1;j<=nAlpha(A);j++) {
                       	sscanf(t,"%lf",&tmpr);
                       	mat_emit[i][uch[j]] = tmpr;
                       	while(!isspace(t[0]))t++;
                       	while(isspace(t[0]))t++;
               	}

		fgets(str,550,fp);
		t = str;
                while(isspace(t[0]))t++;
       	        while(!isspace(t[0]))t++;
               	while(isspace(t[0]))t++;
                for(j=1;j<=nAlpha(A);j++) {
       	                sscanf(t,"%lf",&tmpr);
               	        ins_emit[i][uch[j]] = tmpr;
                       	while(!isspace(t[0]))t++;
               	        while(isspace(t[0]))t++;
               	}

		fgets(str,550,fp);
                t = str;
       	        while(isspace(t[0]))t++;
                while(!isspace(t[0]))t++;
                while(isspace(t[0]))t++;

		if(t[0] == '*') { m2m[i] = 0.; }
		else {
			sscanf(t,"%lf",&(m2m[i]));
		}
		while(!isspace(t[0]))t++;while(isspace(t[0]))t++;	

		if(t[0] == '*') { m2i[i] = 0.; }
		else {
	                sscanf(t,"%lf",&(m2i[i]));
		}
                while(!isspace(t[0]))t++;while(isspace(t[0]))t++;

		if(t[0] == '*') { m2d[i] = 0.; }
		else {
	                sscanf(t,"%lf",&(m2d[i]));
		}
                while(!isspace(t[0]))t++;while(isspace(t[0]))t++;

		if(t[0] == '*') { i2m[i] = 0.; }
		else {
        	        sscanf(t,"%lf",&(i2m[i]));
		}
                while(!isspace(t[0]))t++;while(isspace(t[0]))t++;

		if(t[0] == '*') { i2i[i] = 0.; }
		else {
        	        sscanf(t,"%lf",&(i2i[i]));
     		}
	 	while(!isspace(t[0]))t++;while(isspace(t[0]))t++;

		if(t[0] == '*') { d2m[i] = 0.; }
		else {
	               	sscanf(t,"%lf",&(d2m[i]));
		}
		while(!isspace(t[0]))t++;while(isspace(t[0]))t++;

		if(t[0] == '*') { d2d[i] = 0.; }
		else {
	                sscanf(t,"%lf",&(d2d[i]));
		}
               	while(!isspace(t[0]))t++;while(isspace(t[0]))t++;

		if(t[0] == '*') { b2m[i] = 0.; }
		else {
                	sscanf(t,"%lf",&b2m[i]);
		}
               	while(!isspace(t[0]))t++;while(isspace(t[0]))t++;

		if(t[0] == '*') { m2e[i] = 0.; }
		else {
                	sscanf(t,"%lf",&m2e[i]);
		}
        }


	stpb = new stpb_typ(len,m2m,m2i,m2d,i2i,i2m,d2d,d2m,b2m,m2e,nb,nn,ec,ej,
		ct,cc,jb,jj,bm,bd,gg,gf,bits);

	double x = 0.;
	double *nl;
	NEW(nl,nAlpha(A)+2,double);
	for(j=1;j<=nAlpha(A);j++) nl[j] = BitsToProb(nule[j],1./nAlpha(A),bits);

	for(i=1;i<=len;i++) {
		for(x=0.,j=1;j<=nAlpha(A);j++) {
			x += nl[j]*mat_emit[i][j];
		}
		mat_emit[i][0] = x;

		for(x=0.,j=1;j<=nAlpha(A);j++) {
			x += nl[j]*ins_emit[i][j];
		}
		ins_emit[i][0] = x;
	}		

	shm = new shm_typ(name,desc,nule,mat_emit,ins_emit,stpb,A);
	shm->SetLambda(lambda);
	shm->SetU(U);
	shm->SetSLambda(Slambda);
	shm->SetSU(SU);
	fgets(str,550,fp);

	return shm;
}

void shm_typ::Write(FILE *fp)
{
        Int4    bits = Bits();
        Int4    **mat_emit = Matrix(); 
        Int4    **ins_emit = InsMatrix();        
        Int4    i,j,*mat_emit_i,*ins_emit_i;
        Int4    *m2m = MatToMat( );
        Int4    *m2i = MatToIns( );
        Int4    *m2d = MatToDel( );
        Int4    *i2m = InsToMat( );
        Int4    *i2i = InsToIns( );
        Int4    *d2m = DelToMat( );
        Int4    *d2d = DelToDel( );
        Int4    *b2m = BegToMat( );
        Int4    *m2e = MatToEnd( );
        Int4    length = GetLength();
        Int4    XT[9];
		XT[1] = NB(); XT[2] = NN(); XT[3] = EC(); XT[4] = EJ(); 
		XT[5] = CT(); XT[6] = CC(); XT[7] = JB(); XT[8] = JJ(); 
        Int4    NULT[3];
		NULT[1] = GG(); NULT[2] = GF();
        Int4    *Nle = NULE();
	Int4 	Nule[21];
        double  EVD[3] = {0.,0.000000,0.000000};
        a_type  A = GetAlphabet();
        char    ch[] = "XACDEFGHIKLMNPQRSTVWY";
        unsigned char uch[22];
        for(i=0;i<=nAlpha(A);i++) { uch[i] = AlphaCode(ch[i],A); }
	for(i=1;i<=nAlpha(A);i++) Nule[i] = Nle[uch[i]];
        fprintf(fp,"HMMER2.0\n");
        fprintf(fp,"NAME  %s\n",GetName());
        fprintf(fp,"DESC  %s\n",GetDesc());
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
        fprintf(fp,"         %d      *%7d\n",BM(),BD());

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
                fprintf(fp,"%7d",m2m[i]); 
                fprintf(fp,"%7d",m2i[i]); 
                fprintf(fp,"%7d",m2d[i]); 
                fprintf(fp,"%7d",i2m[i]); 
                fprintf(fp,"%7d",i2i[i]); 
                fprintf(fp,"%7d",d2m[i]); 
                fprintf(fp,"%7d",d2d[i]); 
                if(b2m[i]==BigNegative) fprintf(fp,"      *");
		else fprintf(fp,"%7d",b2m[i]); 
                if(m2e[i]==BigNegative) fprintf(fp,"      *");
 		else fprintf(fp,"%7d",m2e[i]); 

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

void WriteSHM(FILE *fp, shm_typ *shm, shm_typ *rshm)
{
	Int4 	i;
	char 	*name = shm->GetName();
	char 	*desc = shm->GetDesc();
	Int4	shm_len = shm->GetLength();
	Int4	bits = shm->Bits();
	double	lambda = shm->GetLambda();
	double	U = shm->GetU();
	double	Slambda = shm->GetSLambda();
	double	SU = shm->GetSU();
	Int4	wthreshold = shm->GetWthreshold();
	Int4	uxparameter = shm->GetUxparameter();
	Int4	uthreshold = shm->GetUthreshold();
	Int4	xparameter = shm->GetXparameter();
	Int4	gthreshold = shm->GetGthreshold();
	Int4	vthreshold = shm->GetVthreshold();
	a_type 	AB = shm->GetAlphabet();

        fwrite(&lambda,sizeof(double),1,fp);
        fwrite(&U,sizeof(double),1,fp);
        fwrite(&Slambda,sizeof(double),1,fp);
        fwrite(&SU,sizeof(double),1,fp);
        fwrite(&shm_len,sizeof(Int4),1,fp);
        fwrite(&bits,sizeof(Int4),1,fp);

        fwrite(&wthreshold,sizeof(Int4),1,fp);
        fwrite(&uxparameter,sizeof(Int4),1,fp);
        fwrite(&uthreshold,sizeof(Int4),1,fp);
        fwrite(&xparameter,sizeof(Int4),1,fp);
        fwrite(&gthreshold,sizeof(Int4),1,fp);
        fwrite(&vthreshold,sizeof(Int4),1,fp);   

        Int4 *m2m,*m2i,*m2d,*i2i,*i2m,*d2d,*d2m,*b2m,*m2e;
        Int4 *rm2m,*rm2i,*rm2d,*ri2i,*ri2m,*rd2d,*rd2m;
 
        m2m = shm->MatToMat();
        m2i = shm->MatToIns();
        m2d = shm->MatToDel();
        i2i = shm->InsToIns();
        i2m = shm->InsToMat();
        d2d = shm->DelToDel();
        d2m = shm->DelToMat();
        b2m = shm->BegToMat();
        m2e = shm->MatToEnd();
        
        fwrite(m2m,sizeof(Int4),shm_len+1,fp);
        fwrite(m2i,sizeof(Int4),shm_len+1,fp);
        fwrite(m2d,sizeof(Int4),shm_len+1,fp);
        fwrite(i2i,sizeof(Int4),shm_len+1,fp);
        fwrite(i2m,sizeof(Int4),shm_len+1,fp);
        fwrite(d2d,sizeof(Int4),shm_len+1,fp);
        fwrite(d2m,sizeof(Int4),shm_len+1,fp);
        fwrite(b2m,sizeof(Int4),shm_len+1,fp);
        fwrite(m2e,sizeof(Int4),shm_len+1,fp);

	Int4 nb = shm->NB();
	Int4 nn = shm->NN();
	Int4 ec = shm->EC();
	Int4 ej = shm->EJ();
	Int4 ct = shm->CT();
	Int4 cc = shm->CC();
	Int4 jb = shm->JB();
	Int4 jj = shm->JJ();
	fwrite(&nb,sizeof(Int4),1,fp);
	fwrite(&nn,sizeof(Int4),1,fp);
	fwrite(&ec,sizeof(Int4),1,fp);
	fwrite(&ej,sizeof(Int4),1,fp);
	fwrite(&ct,sizeof(Int4),1,fp);
	fwrite(&cc,sizeof(Int4),1,fp);
	fwrite(&jb,sizeof(Int4),1,fp);
	fwrite(&jj,sizeof(Int4),1,fp);

	Int4 bm = shm->BM();
	Int4 bd = shm->BD();
	fwrite(&bm,sizeof(Int4),1,fp);
	fwrite(&bd,sizeof(Int4),1,fp);

	Int4 gg = shm->GG();
	Int4 gf = shm->GF();
	fwrite(&gg,sizeof(Int4),1,fp);
	fwrite(&gf,sizeof(Int4),1,fp);

	Int4 *nule = shm->NULE();
	fwrite(nule,sizeof(Int4),nAlpha(AB)+2,fp);

        rm2m = rshm->MatToMat();
        rm2i = rshm->MatToIns();
        rm2d = rshm->MatToDel();
        ri2i = rshm->InsToIns();
        ri2m = rshm->InsToMat();
        rd2d = rshm->DelToDel();
        rd2m = rshm->DelToMat();
        fwrite(rm2m,sizeof(Int4),shm_len+1,fp);
        fwrite(rm2i,sizeof(Int4),shm_len+1,fp);
        fwrite(rm2d,sizeof(Int4),shm_len+1,fp);
        fwrite(ri2i,sizeof(Int4),shm_len+1,fp);
        fwrite(ri2m,sizeof(Int4),shm_len+1,fp);
        fwrite(rd2d,sizeof(Int4),shm_len+1,fp);
        fwrite(rd2m,sizeof(Int4),shm_len+1,fp);

	Int4 **mat_emit = shm->Matrix();
	Int4 **ins_emit = shm->InsMatrix();
        
        for(i=1;i<=shm_len;i++){
                fwrite(mat_emit[i],sizeof(Int4),nAlpha(AB)+1,fp);
        }
        for(i=1;i<=shm_len;i++){
                fwrite(ins_emit[i],sizeof(Int4),nAlpha(AB)+1,fp);
        }
        for(i=1;i<=shm_len;i++){
                fwrite(mat_emit[i],sizeof(Int4),nAlpha(AB)+1,fp);
        }
        for(i=1;i<=shm_len;i++){
                fwrite(ins_emit[i],sizeof(Int4),nAlpha(AB)+1,fp);
        }

	fwrite(name,sizeof(char),3*sizeof(double),fp);
	fwrite(desc,sizeof(char),6*sizeof(double),fp);
	e_type cons = shm->GetConsensus();
        unsigned char *Ptr_cons = SeqPtr(cons);
        fwrite(Ptr_cons,sizeof(unsigned char),shm_len+1,fp);
	Int4 k = (shm_len+1)%sizeof(double);
	if (k != 0) {
		for(i=1;i<=sizeof(double) - k;i++) {
			char t = '0';
			fwrite(&t,sizeof(char),1,fp);
		}
	}
}

//construct an hmm from probabilities
shm_typ::shm_typ(char *Name, char *Desc, Int4 *Nule, Int4 **Matemit, Int4 **Insemit,
                        Int4 len, idp_typ *idp, double FreqOfMatchAtPositionOne,a_type A)
// This is the only constructor used other than in Read mode...
{
	stpb_typ *Stpb = new stpb_typ(len, idp->MatToMat( ), idp->MatToIns( ), idp->MatToDel( ),
                idp->InsToIns( ), idp->InsToMat( ), idp->DelToDel( ),
                idp->DelToMat( ), FreqOfMatchAtPositionOne);
	Int4 i,j,**Mat_emit, **Ins_emit,*NULE;
	NEWP(Mat_emit,len+2,Int4); NEWP(Ins_emit,len+2,Int4);
	for(i=0;i<=len+1;i++) { NEW(Mat_emit[i],nAlpha(A)+2,Int4); NEW(Ins_emit[i],nAlpha(A)+2,Int4); }
	NEW(NULE,nAlpha(A)+2,Int4);
	for(j=1;j<=nAlpha(A);j++) {
		NULE[j] = ProbToBits(Nule[j],1./nAlpha(A),Stpb->Bits());
	}
	for(i=1;i<=len;i++) {
		for(j=1;j<=nAlpha(A);j++) {
			Mat_emit[i][j] = ProbToBits(Matemit[i][j],Nule[j],Stpb->Bits());
			Ins_emit[i][j] = ProbToBits(Insemit[i][j], Nule[j], Stpb->Bits());
		}
	} init(Name, Desc, NULE, Mat_emit, Ins_emit, Stpb, A);
}

shm_typ::shm_typ(char *Name, char *Desc, double *Nule, double **Matemit, double **Insemit, 
		stpb_typ *Stpb, a_type A)
{
	Int4 i,j,len = Stpb->Length(),**Mat_emit, **Ins_emit,*NULE;
	NEWP(Mat_emit,len+2,Int4); NEWP(Ins_emit,len+2,Int4);
	for(i=0;i<=len+1;i++) { NEW(Mat_emit[i],nAlpha(A)+2,Int4); NEW(Ins_emit[i],nAlpha(A)+2,Int4); }
	NEW(NULE,nAlpha(A)+2,Int4);
	for(j=1;j<=nAlpha(A);j++) {
		NULE[j] = ProbToBits(Nule[j],1./nAlpha(A),Stpb->Bits());
	}
	for(i=1;i<=len;i++) {
		for(j=1;j<=nAlpha(A);j++) {
			Mat_emit[i][j] = ProbToBits(Matemit[i][j],Nule[j],Stpb->Bits());
			Ins_emit[i][j] = ProbToBits(Insemit[i][j],Nule[j],Stpb->Bits());
		}
	} init(Name, Desc, NULE, Mat_emit, Ins_emit, Stpb, A);
}

shm_typ::shm_typ(char *Name, char *Desc, Int4 *NULE, Int4 **Mat_emit, Int4 **Ins_emit, stpb_typ *Stpb, a_type A)
{ init(Name, Desc, NULE, Mat_emit, Ins_emit, Stpb, A); }

shm_typ::shm_typ(char *Name, char *Desc, Int4 *NULE, Int4 **Mat_emit, Int4 **Ins_emit, stpb_typ *Stpb, a_type A, 
	double Lambda, double U, double SLambda, double SU, Int4 Wthreshold, Int4 Uxparameter, Int4 Uthreshold,
	Int4 Xparameter, Int4 Gthreshold, Int4 Vthreshold, e_type Cons)

{ init(Name, Desc, NULE, Mat_emit, Ins_emit, Stpb, A, Lambda, U, SLambda, SU, Wthreshold, Uxparameter,
	Uthreshold, Xparameter, Gthreshold, Vthreshold, Cons); }


shm_typ::shm_typ(char *Name, char *Desc, Int4 *NULE, Int4 **Mat_emit, Int4 **Ins_emit, stpb_typ *Stpb, 
	a_type A, Int4 Rpts)
{ 
	Int4 **tmat_emit = Concat(Mat_emit,Stpb->Length(),nAlpha(A),Rpts);
	Int4 **tins_emit = Concat(Ins_emit,Stpb->Length(),nAlpha(A),Rpts);
	stpb_typ *tstpb = new stpb_typ(Stpb,Rpts);
	init(Name, Desc, NULE, tmat_emit, tins_emit, tstpb, A);
	rpts = Rpts;
}

shm_typ::shm_typ(char file_type, FILE *fp, a_type A)
{
	if(file_type == 'B'){ readShmB(fp,A); }
        f_open = 5*bits;
        f_ox = (Int4) (floor) (5.5*bits + 0.5);
	desc = NULL;
}

unsigned char *shm_typ::ConsensusSeqPtr()
{	
	Int4 		i,j,*mat_i,max,best_j;
	unsigned char 	*ptr_cons;
	NEW(ptr_cons,shm_len+2,unsigned char);
	for(i=1;i<=shm_len;i++){
		mat_i = mat_emit[i];
		best_j = 1; 
		max = mat_i[1];
        	for(j=1;j<=nAlpha(AB);j++) {
                	if (mat_i[j] > max) { max = mat_i[j]; best_j = j; }
                }
        	ptr_cons[i] = best_j;   
        }
	return ptr_cons;
}

void shm_typ::init(char *Name, char *Desc, Int4 *NULE, Int4 **Mat_emit, Int4 **Ins_emit, stpb_typ *Stpb, a_type A)
{
	rpts = 1;
	stpb = Stpb;	
	bits = stpb->Bits();
	f_open = 5*bits;
        f_ox = (Int4) (floor) (5.5*bits + 0.5);
	lambda = 0.; u = 0.; chi = 0.;
	Slambda = 0.; Su = 0.;
        wthreshold = 0; uxparameter = 0; uthreshold = 0;
        xparameter = 0; gthreshold = 0; vthreshold = 0;
	mat_emit = Mat_emit;
	ins_emit = Ins_emit;
	nule = NULE;
	shm_len = stpb->Length();
	consensus = NULL;
	desc = Desc;
	name = Name;
	AB = A;
}

void shm_typ::init(char *Name, char *Desc, Int4 *NULE, Int4 **Mat_emit, Int4 **Ins_emit, stpb_typ *Stpb, a_type A,
	double Lambda, double U, double SLambda, double SU, Int4 Wthreshold, Int4 Uxparameter,
        Int4 Uthreshold, Int4 Xparameter, Int4 Gthreshold, Int4 Vthreshold, e_type Cons)
{
	rpts = 1;
	stpb = Stpb;	
	bits = stpb->Bits();
	f_open = 5*bits;
        f_ox = (Int4) (floor) (5.5*bits + 0.5);
	lambda = Lambda; u = U; chi = 0.;
	Slambda = SLambda; Su = SU;
        wthreshold = Wthreshold; 
	uxparameter = Uxparameter;
	uthreshold = Uthreshold;
	xparameter = Xparameter;
	gthreshold = Gthreshold;
	vthreshold = Vthreshold;
	mat_emit = Mat_emit;
	ins_emit = Ins_emit;
	nule = NULE;
	shm_len = stpb->Length();
	consensus = Cons;
	name = Name;
	desc = Desc;
	AB = A;
}

void shm_typ::readShmB(FILE *fp, a_type A)
{
	Int4 i;
	NEW(name,11,char);
	fread(name,sizeof(char),10,fp);
	fread(&shm_len,sizeof(Int4),1,fp);
	fread(&bits,sizeof(Int4),1,fp);
	fread(&lambda,sizeof(double),1,fp);
	fread(&u,sizeof(double),1,fp);
	fread(&Slambda,sizeof(double),1,fp);
	fread(&Su,sizeof(double),1,fp);

        fread(&wthreshold,sizeof(Int4),1,fp);
        fread(&uxparameter,sizeof(Int4),1,fp);
        fread(&uthreshold,sizeof(Int4),1,fp);
        fread(&xparameter,sizeof(Int4),1,fp);
        fread(&gthreshold,sizeof(Int4),1,fp);
        fread(&vthreshold,sizeof(Int4),1,fp);

	unsigned char *ptr_cons;
	MEW(ptr_cons,shm_len+2,unsigned char);
	fread(ptr_cons,sizeof(unsigned char),shm_len+1,fp);
	consensus = MkSeq(name,shm_len,ptr_cons);
	free(name); free(ptr_cons);
	Int4 *m2m,*m2i,*m2d,*i2i,*i2m,*d2d,*d2m,*b2m,*m2e;

	AB = A;
	Int4 nalpha = nAlpha(A);

	NEWP(mat_emit,shm_len+2,Int4);
	NEWP(ins_emit,shm_len+2,Int4);

        MEW(m2m,shm_len+2,Int4);
        MEW(m2i,shm_len+2,Int4);
        MEW(m2d,shm_len+2,Int4); 
        MEW(i2i,shm_len+2,Int4);   
        MEW(i2m,shm_len+2,Int4);
        MEW(d2d,shm_len+2,Int4);
        MEW(d2m,shm_len+2,Int4);
        MEW(b2m,shm_len+2,Int4);
        MEW(m2e,shm_len+2,Int4);

        fread(m2m,sizeof(Int4),shm_len+1,fp);
        fread(m2i,sizeof(Int4),shm_len+1,fp);
        fread(m2d,sizeof(Int4),shm_len+1,fp);
        fread(i2i,sizeof(Int4),shm_len+1,fp);
        fread(i2m,sizeof(Int4),shm_len+1,fp);
        fread(d2d,sizeof(Int4),shm_len+1,fp);
        fread(d2m,sizeof(Int4),shm_len+1,fp);
        fread(b2m,sizeof(Int4),shm_len+1,fp);
        fread(m2e,sizeof(Int4),shm_len+1,fp);

	Int4 nb,nn,ec,ej,ct,cc,jb,jj;
        fread(&nb,sizeof(Int4),1,fp);
        fread(&nn,sizeof(Int4),1,fp);
        fread(&ec,sizeof(Int4),1,fp);
        fread(&ej,sizeof(Int4),1,fp);
        fread(&ct,sizeof(Int4),1,fp);
        fread(&cc,sizeof(Int4),1,fp);
        fread(&jb,sizeof(Int4),1,fp);
        fread(&jj,sizeof(Int4),1,fp);

	Int4 bm,bd;
        fread(&bm,sizeof(Int4),1,fp);
        fread(&bd,sizeof(Int4),1,fp);

	Int4 gg,gf;
        fread(&gg,sizeof(Int4),1,fp);
        fread(&gf,sizeof(Int4),1,fp);

	MEW(nule,nalpha+1,Int4);
	fread(nule,sizeof(Int4),nalpha+1,fp);

	Int4 *tmat_emit,*tins_emit;
	MEW(tmat_emit,shm_len*(nalpha+1),Int4);
	MEW(tins_emit,shm_len*(nalpha+1),Int4);
	fread(tmat_emit,sizeof(Int4),shm_len*(nalpha+1),fp);
	fread(tins_emit,sizeof(Int4),shm_len*(nalpha+1),fp);

	for(i=1;i<=shm_len;i++){
		mat_emit[i] = tmat_emit;
		ins_emit[i] = tins_emit;
		tmat_emit += (nalpha+1);
		tins_emit += (nalpha+1);
	}

        stpb = new stpb_typ(shm_len,m2m,m2i,m2d,i2i,i2m,d2d,d2m,b2m,m2e,
		nb,nn,ec,ej,ct,cc,jb,jj,bm,bd,gg,gf,bits);
}

void shm_typ::Put(FILE *fp)
{
	stpb->Put(fp);
}

char *shm_typ::loc_Align(e_type E, Int4 *Score, Int4 *StartQ, Int4 *StartS, Int4 *Oper_len)
{ return loc_Viterbi(E,Score,StartQ,StartS,Oper_len); }

char *shm_typ::loc_AlignWOTB(e_type E)
{ return loc_ViterbiWOTB(E); }

char	*shm_typ::Align(e_type E, Int4 *Score, Int4 *Start, Int4 *Oper_len) 
{ return Viterbi(E,Score,Start,Oper_len); }

void shm_typ::Free( )
{
	Int4 i;
	if(mat_emit != NULL) {
		for(i=0;i<=shm_len;i++) free(mat_emit[i]);
		free(mat_emit);
	}

	if(ins_emit != NULL) {
		for(i=0;i<=shm_len;i++) free(ins_emit[i]);
		free(ins_emit);
	}
	if(name != NULL) free(name);
	if(desc != NULL) free(desc);
	if(nule != NULL) free(nule);
	if(consensus != NULL) NilSeq(consensus);
	delete stpb; 
}

void shm_typ::loc_inner_loop(Int4 j, Int4 *matj, Int4 *matjm1, Int4 *insj,
	Int4 *insjm1, Int4 *delj, Int4 *deljm1, Int4 m2m, Int4 d2m, Int4 i2m, 
	Int4  m2d, Int4  d2d, Int4  m2i, Int4 i2i, Int4 seq_len, unsigned char *seq, 
	Int4 *scorej, Int4 *ins_emition)
{
   register Int4 i,im1,M,I,D;
        
   for(im1=0,i=1;i<=seq_len;im1++,i++){
        
        // Find optimum path to match(j,i) state.
	if(matjm1[im1] <= 0) M = 0; else M = matjm1[im1] + m2m;
	if(deljm1[im1] <= 0) D = 0; else D = deljm1[im1] + d2m;
	if(insjm1[im1] <= 0) I = 0; else I = insjm1[im1] + i2m;
        if(M >= D && M >= I){ matj[i] = M + scorej[seq[i]]; traceM[j][i] = 'm'; }
        else if(D >= I){ matj[i] = D + scorej[seq[i]]; traceM[j][i] = 'd'; }
        else { matj[i] = I + scorej[seq[i]]; traceM[j][i] = 'i'; }
	if(matj[i] <= 0) { matj[i] = 0; traceM[j][i] = 'B';}

        // Find optimum path to delete(j,i) state.
	if(matjm1[i] <= 0) M = 0; else M = matjm1[i] + m2d;
        if(deljm1[i] <= 0) D = 0; else D = deljm1[i] + d2d;
        if(M >= D){ delj[i] = M; traceD[j][i] = 'm'; }
        else { delj[i] = D; traceD[j][i] = 'd'; }
	if(delj[i] <= 0) { delj[i] = 0; traceD[j][i] = 'B';}

        // Find optimum path to insert(j,i) state.
	if(matj[im1] <= 0 ) M = 0; else M = matj[im1] + m2i; 
	if(insj[im1] <= 0) I = 0; else I = insj[im1] + i2i;
        if(M >= I){ insj[i] = M + ins_emition[seq[i]]; traceI[j][i] = 'm'; }
        else { insj[i] = I + ins_emition[seq[i]]; traceI[j][i] = 'i'; }
	if(insj[i] <= 0) { insj[i] = 0; traceI[j][i] = 'B';}
   }
}

void shm_typ::inner_loop(Int4 j, Int4 *matj, Int4 *matjm1, Int4 *insj, 
	Int4 *insjm1, Int4 *delj, Int4 *deljm1, Int4 m2m, Int4 d2m, 
	Int4 i2m, Int4  m2d, Int4  d2d, Int4  m2i, Int4 i2i, 
	Int4 seq_len, unsigned char *seq, Int4 *scorej, Int4 *ins_emition)
// about 80% of compute time is spent here...
{
   register Int4 i,im1,M,I,D;

   for(im1=0,i=1;i<=seq_len;im1++,i++){

	// Find optimum path to match(j,i) state.
	M = matjm1[im1] + m2m; D = deljm1[im1] + d2m; I = insjm1[im1] + i2m;  // = to M transitions
        if(M >= D && M >= I){ matj[i]= M + scorej[seq[i]]; traceM[j][i]='m'; }	// m2m 
        else if(D >= I){ matj[i] = D + scorej[seq[i]]; traceM[j][i]='d'; } 	// d2m
        else { matj[i] = I + scorej[seq[i]]; traceM[j][i]='i'; }		// i2m

	// Find optimum path to delete(j,i) state.			 
	M = matjm1[i] + m2d; D = deljm1[i] + d2d;		// == to D transitions
        if(M > D){ delj[i]=M; traceD[j][i]='m'; } 		// m2d
	else { delj[i] = D; traceD[j][i]='d'; }			// d2d

	// Find optimum path to insert(j,i) state.
	M = matj[im1] + m2i; I = insj[im1] + i2i;		// == to I transitions
        if(M > I){ insj[i]= M + ins_emition[seq[i]]; traceI[j][i]='m'; }  	// m2i
	else { insj[i]= I + ins_emition[seq[i]]; traceI[j][i]='i'; }		// i2i
   }
}

char *shm_typ::loc_Viterbi(e_type E, Int4 *Score, Int4 *StartQ, Int4 *StartS, Int4 *Oper_len)
{
        Int4    i,j,seq_len=LenSeq(E);
		
        NEWP(MAT,shm_len+2,Int4); NEWP(DEL,shm_len+2,Int4); NEWP(INS,shm_len+2,Int4);
        NEWP(traceM,shm_len+2,char); NEWP(traceD,shm_len+2,char); 
        NEWP(traceI,shm_len+2,char);
        for(i=0;i<=shm_len;i++){
                MEW(MAT[i],seq_len+1,Int4); MEW(DEL[i],seq_len+1,Int4);
                MEW(INS[i],seq_len+1,Int4);
                NEW(traceM[i],seq_len+2,char); NEW(traceD[i],seq_len+2,char); 
                NEW(traceI[i],seq_len+2,char);
        }
        traceM[0][0]=traceI[0][0]=traceD[0][0]='B';
        Int4 *d2d = stpb->DelToDel();
        for(j=0;j<=shm_len;j++) {
                MAT[j][0]=0; traceM[j][0] = 'B';
                INS[j][0]=0; traceI[j][0] = 'B';
                DEL[j][0]=0; traceD[j][0] = 'B';
        }
        Int4 i2i0 = stpb->InsToIns(0);
        for(i=1;i<=seq_len;i++) {
           DEL[0][i] = 0; traceD[0][i] = 'B';
           MAT[0][i] = 0; traceM[0][i] = 'B';
           INS[0][i] = 0; traceI[0][i] = 'B';
        }
        for(Int4 J=1,Jm1=0;J<=shm_len;J++,Jm1++){
           loc_inner_loop(J,MAT[J],MAT[Jm1],INS[J],INS[Jm1],DEL[J],DEL[Jm1],
                stpb->MatToMat(Jm1),stpb->DelToMat(Jm1),stpb->InsToMat(Jm1),
                stpb->MatToDel(Jm1),stpb->DelToDel(Jm1),
                stpb->MatToIns(J),stpb->InsToIns(J),
                seq_len,SeqPtr(E),mat_emit[J],ins_emit[1]);

        }
        char *operation=loc_get_traceback(seq_len,Oper_len,StartQ,StartS,Score);
        for(i=0;i<=shm_len;i++){free(MAT[i]);free(DEL[i]);free(INS[i]);}
        for(i=0;i<=shm_len;i++){free(traceM[i]); free(traceI[i]); free(traceD[i]); }
        free(traceM); free(traceI); free(traceD); free(MAT);free(DEL);free(INS);
        return operation;
}

char *shm_typ::loc_ViterbiWOTB(e_type E)
{
        Int4    i,j,seq_len=LenSeq(E);
		
        NEWP(MAT,shm_len+2,Int4); NEWP(DEL,shm_len+2,Int4); NEWP(INS,shm_len+2,Int4);
        NEWP(traceM,shm_len+2,char); NEWP(traceD,shm_len+2,char); 
        NEWP(traceI,shm_len+2,char);
        for(i=0;i<=shm_len;i++){
                MEW(MAT[i],seq_len+1,Int4); MEW(DEL[i],seq_len+1,Int4);
                MEW(INS[i],seq_len+1,Int4);
                NEW(traceM[i],seq_len+2,char); NEW(traceD[i],seq_len+2,char); 
                NEW(traceI[i],seq_len+2,char);
        }
        traceM[0][0]=traceI[0][0]=traceD[0][0]='B';
        Int4 *d2d = stpb->DelToDel();
        for(j=0;j<=shm_len;j++) {
                MAT[j][0]=0; traceM[j][0] = 'B';
                INS[j][0]=0; traceI[j][0] = 'B';
                DEL[j][0]=0; traceD[j][0] = 'B';
        }
        Int4 i2i0 = stpb->InsToIns(0);
        for(i=1;i<=seq_len;i++) {
           DEL[0][i] = 0; traceD[0][i] = 'B';
           MAT[0][i] = 0; traceM[0][i] = 'B';
           INS[0][i] = 0; traceI[0][i] = 'B';
        }
        for(Int4 J=1,Jm1=0;J<=shm_len;J++,Jm1++){
           loc_inner_loop(J,MAT[J],MAT[Jm1],INS[J],INS[Jm1],DEL[J],DEL[Jm1],
                stpb->MatToMat(Jm1),stpb->DelToMat(Jm1),stpb->InsToMat(Jm1),
                stpb->MatToDel(Jm1),stpb->DelToDel(Jm1),
                stpb->MatToIns(J),stpb->InsToIns(J),
                seq_len,SeqPtr(E),mat_emit[J],ins_emit[1]);

        }
//        char *operation=loc_get_traceback(seq_len,Oper_len,StartQ,StartS,Score);
        for(i=0;i<=shm_len;i++){free(MAT[i]);free(DEL[i]);free(INS[i]);}
        for(i=0;i<=shm_len;i++){free(traceM[i]); free(traceI[i]); free(traceD[i]); }
        free(traceM); free(traceI); free(traceD); free(MAT);free(DEL);free(INS);
        return NULL;
}

void	shm_typ::MkConsensus( )
{
	unsigned char *seq;

	if(consensus) {
		// PutSeq(stderr,consensus,AB); 
		return; 
	}
	NEW(seq,shm_len+3,unsigned char);
        for(Int4 J=1;J<=shm_len;J++){
	   Int4 max_score=0;
	   unsigned char max_r=0;
           for(unsigned char r=1;r<=nAlpha(AB);r++){
		if(max_score < mat_emit[J][r]) {
			max_score=mat_emit[J][r];
			max_r=r;
		}
	   } seq[J]=max_r;
	}
	consensus=MkSeq("concensus seq for hmm",shm_len,seq);
	// PutSeq(stderr,consensus,AB);
}

char *shm_typ::Viterbi(e_type E, Int4 *Score, Int4 *Start, Int4 *Oper_len)
{
        Int4    i,j,im1,jm1,seq_len=LenSeq(E);

        NEWP(MAT,shm_len+2,Int4); NEWP(DEL,shm_len+2,Int4); NEWP(INS,shm_len+2,Int4);
	NEWP(traceM,shm_len+2,char); NEWP(traceD,shm_len+2,char); 
	NEWP(traceI,shm_len+2,char);
        for(i=0;i<=shm_len;i++){
                MEW(MAT[i],seq_len+1,Int4); MEW(DEL[i],seq_len+1,Int4);
                MEW(INS[i],seq_len+1,Int4);
		NEW(traceM[i],seq_len+2,char); NEW(traceD[i],seq_len+2,char); 
		NEW(traceI[i],seq_len+2,char);
        }
        DEL[0][0]=stpb->MatToDel(0);	// A kludge; these should be done at DEL[1][0] 
	INS[0][0]=stpb->MatToIns(0);	// set to zero in tpb_typ for now.
	MAT[0][0]=stpb->MatToMat(0);     // set to zero in tpb_typ for now.
	traceM[0][0]=traceI[0][0]=traceD[0][0]='B';
	Int4 *d2d=stpb->DelToDel();
        for(jm1=0,j=1;j<=shm_len;jm1++,j++) {
		MAT[j][0]=BigNegative; traceM[j][0] = ' ';	// impossible;
		INS[j][0]=BigNegative; traceI[j][0] = ' ';	// impossible;
		DEL[j][0]=DEL[jm1][0] + d2d[jm1]; traceD[j][0] = 'd';
        }
	Int4 i2i0 = stpb->InsToIns(0);
	for(im1=0,i=1;i<=seq_len;im1++,i++) {
	   DEL[0][i] = BigNegative; traceD[0][i] = ' ';  // impossible
	   MAT[0][i] = BigNegative; traceM[0][i] = ' ';  // impossible
	   INS[0][i] = INS[0][im1] + i2i0; traceI[0][i] = 'i';
	}
        for(Int4 J=1,Jm1=0;J<=shm_len;J++,Jm1++){
           inner_loop(J,MAT[J],MAT[Jm1],INS[J],INS[Jm1],DEL[J],DEL[Jm1],
		stpb->MatToMat(Jm1),stpb->DelToMat(Jm1),stpb->InsToMat(Jm1),
		stpb->MatToDel(Jm1),stpb->DelToDel(Jm1),
		stpb->MatToIns(J),stpb->InsToIns(J),
		seq_len,SeqPtr(E),mat_emit[J],ins_emit[J]);
        } j=shm_len; i=seq_len;
	if(MAT[j][i] > DEL[j][i] && MAT[j][i] > INS[j][i]) *Score=MAT[j][i];
	else if(DEL[j][i] > INS[j][i]) *Score = DEL[j][i]; else *Score = INS[j][i];
        char *operation=get_traceback(seq_len,Oper_len,Start);
        for(i=0;i<=shm_len;i++){free(MAT[i]);free(DEL[i]);free(INS[i]);}
        for(i=0;i<=shm_len;i++){free(traceM[i]); free(traceI[i]); free(traceD[i]); }
	free(traceM); free(traceI); free(traceD); free(MAT);free(DEL);free(INS);
        return operation;
}

char *shm_typ::loc_get_traceback(Int4 seq_len, Int4 *oper_len, Int4 *startQ, Int4 *startS, Int4 *score)
{
        Int4    i,j,k,max_i=seq_len,max_j=shm_len,max,*MAT_j;
        char    state,*operation,*back_operation;

        for(max=0,j=1;j<=shm_len;j++){
                MAT_j = MAT[j];
                for(i=1;i<=seq_len;i++){
                        if(MAT_j[i] > max) { max = MAT_j[i]; max_j = j; max_i = i; } 
                }
        }
	*startS=max_i+1; *startQ=max_j+1;
        NEW(back_operation,seq_len+shm_len+3,char);
        *score = MAT[max_j][max_i];
        for(k=0,j=max_j,i=max_i,state='E';state != 'X'; ){
           switch(state){
                case 'B': // begin state (exit traceback).
                        back_operation[k++]='E'; state='X'; break;
                case 'E': // end state.  
                           back_operation[k++]='E'; 
                           state = 'm'; break;
                case 'm': case 'M': // previously sampled matched state.
                           state = traceM[j][i]; 
			   if (state != 'B') {
			     *startS -= 1; *startQ -= 1;
                             back_operation[k++]='m';
			     i--; j--; 
			   } break;
                case 'd': case 'D': // previously sampled deletion
                           state = traceD[j][i];
			   if (state != 'B') {
			     *startQ -= 1;
                             back_operation[k++]='d';
			     j--; 
			   } break;
                case 'i':  case 'I': // previously sampled insertion.
                       	    state = traceI[j][i]; 
			    if (state != 'B') {
			      *startS -= 1;
                              back_operation[k++]='I';
			      i--; 
			    } break;
                default: 
                     fprintf(stderr,"seq_len=%d\n",seq_len);
                     fprintf(stderr,"state = %c; i=%d; j = %d; k = %d\n", state,i,j,k);
                     back_operation[k]=0; std::cerr << back_operation; std::cerr << std::endl;
                     print_error("this should not happen");
           }
        } *oper_len = k-2;
	for(i=k-2;i>1;i--) {
		if(back_operation[i] != 'm') {
			*oper_len -= 1;
			if(back_operation[i] == 'I') *startS += 1;
			else *startQ += 1;
			back_operation[i] = 'E'; 
		} else break;
	}
        NEW(operation,k+4,char);
        for(i=0;i<=*oper_len+1;i++) operation[i]=back_operation[*oper_len+1-i];
        free(back_operation);
        return operation;
}


char *shm_typ::get_traceback(Int4 seq_len, Int4 *oper_len, Int4 *start)
{
        Int4    i,j,k;
        char    state,*operation, *back_operation;

        NEW(back_operation,seq_len+shm_len+3,char);
        for(k=0,j=shm_len,i=seq_len,state='E';state != 'X'; ){
           switch(state){
                case 'B': // begin state (exit traceback).
                        back_operation[k++]='E'; *start=i; state='X'; break;
                case 'E': // end state.  
                        back_operation[k++]='E'; 
			if(MAT[j][i] > DEL[j][i] && MAT[j][i] > INS[j][i]) state = 'm';
			else if(DEL[j][i] > INS[j][i]){ state = 'd'; }
			else state = 'i';
                        break;
                case 'm': case 'M': // previously sampled matched state.
                        if(j <= 1) { back_operation[k++]='M'; state='B'; }
			else {
                           back_operation[k++]='m';
			   state = traceM[j][i]; j--; i--; 
			} break;
                case 'd': case 'D': // previously sampled deletion
                        if(j <= 1){ back_operation[k++]='D'; state='B'; i++; }
			else {
                             back_operation[k++]='d';
			     state = traceD[j][i]; j--;
			} break;
                case 'i':  case 'I': // previously sampled insertion.
                        back_operation[k++]='I';
			state = traceI[j][i]; i--;
                        break;                   
                default: 
		     fprintf(stderr,"seq_len=%d\n",seq_len);
		     fprintf(stderr,"state = %c; i=%d; j = %d; k = %d\n", state,i,j,k);
		     back_operation[k]=0; std::cerr << back_operation; std::cerr << std::endl;
		     print_error("this should not happen");
           }
        } *oper_len = k;
        NEW(operation,k+4,char);
        for(i=0;i<k;i++) operation[i]=back_operation[k-i-1];
        free(back_operation);
#if 1	// Temporary fix...
	for(i=*oper_len-1; i > 0; i--){
	   if(!(operation[i] == 'i' || operation[i] == 'E')){
		operation[i+1] = 'E'; operation[i+2] = 0; 
		*oper_len = strlen(operation); break;
	   }
	} 
#endif 
        return operation;
}

shm_typ *shm_typ::ReverseSHM()
{
	shm_typ *reversed_shm;
	stpb_typ *reversed_stpb = stpb->ReverseStpb();
	Int4 **reversed_mat_emit;
	NEWP(reversed_mat_emit,shm_len+3,Int4);
	Int4 **reversed_ins_emit;
	NEWP(reversed_ins_emit,shm_len+3,Int4);
	for(Int4 i=0;i<=shm_len+1;i++) { 
		NEW(reversed_mat_emit[i],nAlpha(AB)+3,Int4);
		NEW(reversed_ins_emit[i],nAlpha(AB)+3,Int4);
		for(Int4 j=0;j<=nAlpha(AB);j++) {
			reversed_mat_emit[i][j] = mat_emit[i][j];
			reversed_ins_emit[i][j] = ins_emit[i][j];
		}
	}
	reversed_shm = new shm_typ((char *) "reversed", (char *) "none", nule, reversed_mat_emit,
			reversed_ins_emit,reversed_stpb,AB);
	return reversed_shm;
}

shm_typ *shm_typ::ConcatenateSHM(Int4 rpts)
{
	shm_typ *concat_shm;
	Int4 i,j,k;
	stpb_typ *concat_stpb = stpb->ConcatenateStpb(rpts);
	Int4 **concat_mat_emit, *tconcat_mat_emit, *mat_emit_j;
	Int4 **concat_ins_emit, *tconcat_ins_emit, *ins_emit_j;
	NEWP(concat_mat_emit,rpts*shm_len+3,Int4);
	NEWP(concat_ins_emit,rpts*shm_len+3,Int4);
	for(i=0;i<rpts*shm_len+3;i++) {
		NEW(concat_mat_emit[i],nAlpha(AB)+3,Int4);
		NEW(concat_ins_emit[i],nAlpha(AB)+3,Int4);
	}

	for(i=0;i<rpts;i++) {
		for(j=1;j<=shm_len;j++) {
			tconcat_mat_emit = concat_mat_emit[i*shm_len+j];
			mat_emit_j = mat_emit[j];
			for(k=0;k<=nAlpha(AB);k++) {
				tconcat_mat_emit[k] = mat_emit_j[k];
			}
		}
	}

	for(i=0;i<rpts;i++) {
		for(j=1;j<=shm_len-1;j++) {
			tconcat_ins_emit = concat_ins_emit[i*shm_len+j];
			ins_emit_j = ins_emit[j];
			for(k=0;k<=nAlpha(AB);k++) {
				tconcat_ins_emit[k] = ins_emit_j[k];
			}
		}
	}

	for(i=1;i<rpts;i++) {
		for(k=0;k<=nAlpha(AB);k++) {
//			concat_ins_emit[i*shm_len][k] = nule[k];
			concat_ins_emit[i*shm_len][k] = 0;

		}
	}
	concat_shm = new shm_typ((char *) "concat", (char *) "none", nule, concat_mat_emit,
			concat_ins_emit, concat_stpb, AB);
	return concat_shm;
}

shm_typ *shm_typ::TempSHM(double temp, double temp_trans)
{
	shm_typ 	*temp_shm;
        stpb_typ 	*temp_stpb = stpb->TempStpb(temp_trans);
	Int4 		i,j;
	double 		sum;
        Int4 		**temp_mat_emit,**temp_ins_emit;
        NEWP(temp_mat_emit,shm_len+3,Int4);
        NEWP(temp_ins_emit,shm_len+3,Int4);
	for(i=0;i<=shm_len+1;i++) {
		NEW(temp_mat_emit[i],nAlpha(AB)+1,Int4);
		NEW(temp_ins_emit[i],nAlpha(AB)+1,Int4);
	}
	double		**mat_emit_prob = Mat_emit_Prob(), *mat_emit_prob_i;
	double		**ins_emit_prob = Ins_emit_Prob(), *ins_emit_prob_i;
	double 		*mat_emit_p, *ins_emit_p;
	NEW(mat_emit_p,nAlpha(AB)+1,double);
	NEW(ins_emit_p,nAlpha(AB)+1,double);

        for(i=1;i<=shm_len;i++) {
		mat_emit_prob_i = mat_emit_prob[i];
                for(sum = 0.,j=1;j<=nAlpha(AB);j++) {
                        sum += (mat_emit_p[j] = pow(mat_emit_prob_i[j],1./temp));
                }
                for(j=1;j<=nAlpha(AB);j++) {
                        temp_mat_emit[i][j] = ProbToBits(mat_emit_p[j]/sum,nule[j],bits);
                }
		ins_emit_prob_i = ins_emit_prob[i];
                for(sum = 0.,j=1;j<=nAlpha(AB);j++) {
                        sum += (ins_emit_p[j] = pow(ins_emit_prob_i[j],1./temp));
                }
                for(j=1;j<=nAlpha(AB);j++) {
                        temp_ins_emit[i][j] = ProbToBits(ins_emit_p[j]/sum,nule[j],bits);
                }
        }
        for(i=1;i<=shm_len;i++) { free(mat_emit_prob[i]); free(ins_emit_prob[i]); }
        free(mat_emit_prob); free(ins_emit_prob);
        free(mat_emit_p);free(ins_emit_p);
        Int4 *NULE;
        NEW(NULE,nAlpha(AB)+1,Int4);
        for(j=1;j<=nAlpha(AB);j++) {
                NULE[j] = nule[j];
        }
	temp_shm = new shm_typ(NULL, NULL, NULE, temp_mat_emit,
			temp_ins_emit,temp_stpb,AB);
	return temp_shm;
}

Int4 shm_typ::WordScore(unsigned char *word, Int4 word_len, Int4 position)
{
	Int4 j,k,score=0;
	assert(position >=1 && position <= shm_len-word_len+1);
	for(k=1,j=position;j<=position+word_len-1;j++,k++) {
		score += mat_emit[j][word[k]];
	}
	return score;
}

Int4 shm_typ::WordHits(unsigned char *word, Int4 word_len, Int4 threshold, unsigned short *position)
//returns number of hits with score above threshold and positions of those hits
{
	unsigned short i;
	Int4 j,tem_score,p,k;
	assert(word_len < shm_len);
	for(k=1,i=1;i<=shm_len-word_len+1;i++) {
		tem_score = 0;
		for(j=1;j<=word_len;j++){
			p = word[j];
			tem_score += mat_emit[i+j-1][p];
		}
		if(tem_score >= threshold) {
			position[k] = i;
			k++;
		}
	}
	return (k-1);
}

double **shm_typ::Mat_emit_Prob()
{
	Int4 	i,j,nA = nAlpha(AB);
	double  **mat_emit_prob,*mat_emit_prob_i,sum;
	Int4	*mat_emit_i;
	NEWP(mat_emit_prob,shm_len+2,double);
	for(i=1;i<=shm_len;i++) {
		NEW(mat_emit_prob[i],nA+1,double);
	}
	for(i=1;i<=shm_len;i++) {
		mat_emit_prob_i = mat_emit_prob[i];
		mat_emit_i = mat_emit[i];
		for(sum=0.,j=1;j<=nA;j++) {
			sum += (mat_emit_prob_i[j] = BitsToProb(mat_emit_i[j],BitsToProb(nule[j],1./nA,bits),bits));
		}
		for(j=1;j<=nA;j++) {
			mat_emit_prob_i[j] /= sum;
		}
	}
	return mat_emit_prob;
}

double **shm_typ::Ins_emit_Prob()
{
	Int4 	i,j,nA = nAlpha(AB);
	double  **ins_emit_prob,*ins_emit_prob_i,sum;
	Int4	*ins_emit_i;
	NEWP(ins_emit_prob,shm_len+2,double);
	for(i=1;i<=shm_len;i++) {
		NEW(ins_emit_prob[i],nA+1,double);
	}
	for(i=1;i<=shm_len-1;i++) {
		ins_emit_prob_i = ins_emit_prob[i];
		ins_emit_i = ins_emit[i];
		for(sum=0.,j=1;j<=nA;j++) {
			sum += (ins_emit_prob_i[j] = BitsToProb(ins_emit_i[j],BitsToProb(nule[j],1./nA,bits),bits));
		}
		for(j=1;j<=nA;j++) {
			ins_emit_prob_i[j] /= sum;
		}
	}
	return ins_emit_prob;
}

double *shm_typ::Mat_emit_Prob(Int4 pos)
{
	assert(pos>0 && pos<=shm_len);
	Int4 	j;
	Int4 	*mat_emit_pos = mat_emit[pos];
	double  *mat_emit_prob;
	NEW(mat_emit_prob,nAlpha(AB)+1,double);
	for(j=1;j<=nAlpha(AB);j++) {
		mat_emit_prob[j] = BitsToProb(nule[j],1./20.,bits)*exp2((double) mat_emit_pos[j]/(double) bits);
	}
	return mat_emit_prob;
}

double *shm_typ::Ins_emit_Prob(Int4 pos)
{
	assert(pos>0 && pos<=shm_len);
	Int4 	j;
	Int4    *ins_emit_pos = ins_emit[pos];
	double  *ins_emit_prob;
	NEW(ins_emit_prob,nAlpha(AB)+1,double);
	for(j=1;j<=nAlpha(AB);j++) {
		ins_emit_prob[j] = BitsToProb(nule[j],1./20.,bits)*exp2((double) ins_emit_pos[j]/(double) bits);
	}
	return ins_emit_prob;
}

void shm_typ::EmitSeqSHM(FILE *fp, Int4 nseq, double emiss_tmpr, double trans_tmpr)
{
	e_type *E = EmitSeqSHM(nseq, emiss_tmpr, trans_tmpr);
	for(Int4 i=1;i<=nseq;i++) {
		PutSeq(fp,E[i],AB);
	}
}

e_type *shm_typ::EmitSeqSHM(Int4 nseq, double emiss_tmpr, double trans_tmpr)
{
// FIND STPB
        Int4    	i,j,nA=nAlpha(AB);
        double          *m2mp,*m2ip,*m2dp,*i2mp,*i2ip,*d2mp,*d2dp,*b2mp,*m2ep;
        double          nbp,nnp,ecp,ejp,ctp,ccp,jbp,jjp;
        double          sum=0.,p1 = 350./351;
        NEW(m2mp,shm_len+2,double);NEW(m2ip,shm_len+2,double);NEW(m2dp,shm_len+2,double);
        NEW(i2mp,shm_len+2,double);NEW(i2ip,shm_len+2,double);NEW(d2mp,shm_len+2,double);
        NEW(d2dp,shm_len+2,double);NEW(b2mp,shm_len+2,double);NEW(m2ep,shm_len+2,double);

        Int4 *m2m,*m2i,*m2d,*i2i,*i2m,*d2d,*d2m,*b2m,*m2e;
        m2m = MatToMat(); m2i = MatToIns(); m2d = MatToDel();
        i2i = InsToIns(); i2m = InsToMat(); d2d = DelToDel();
        d2m = DelToMat(); b2m = BegToMat(); m2e = MatToEnd();
        Int4 nb,nn,ec,ej,ct,cc,jb,jj;
        nb = NB(); nn = NN();
        ec = EC(); ej = EJ();
        ct = CT(); cc = CC();
        jb = JB(); jj = JJ();
           
        for(i=1;i<=shm_len-1;i++) {
                m2mp[i] = pow(BitsToProb(m2m[i],p1,bits),1./trans_tmpr);
                m2ip[i] = pow(BitsToProb(m2i[i],p1,bits),1./trans_tmpr);
                m2dp[i] = pow(BitsToProb(m2d[i],1.0,bits),1./trans_tmpr);
                i2mp[i] = pow(BitsToProb(i2m[i],p1,bits),1./trans_tmpr);
                i2ip[i] = pow(BitsToProb(i2i[i],p1,bits),1./trans_tmpr);
                if (i == 1) { d2mp[i] = d2dp[i] = 0; }
                else {
                        d2mp[i] = pow(BitsToProb(d2m[i],p1,bits),1./trans_tmpr);
                        d2dp[i] = pow(BitsToProb(d2d[i],1.0,bits),1./trans_tmpr);
                }
                b2mp[i] = pow(BitsToProb(b2m[i],p1,bits),1./trans_tmpr);
                m2ep[i] = pow(BitsToProb(m2e[i],1.0,bits),1./trans_tmpr);
        }
        d2dp[shm_len-1] = 0.; 
        b2mp[shm_len] = pow(BitsToProb(b2m[shm_len],p1,bits),1./trans_tmpr);
        m2ep[shm_len] = 1.;
        d2mp[shm_len] = 1.;

        for(i=1;i<=shm_len-1;i++) {
                sum = m2mp[i] + m2ip[i] + m2dp[i] + m2ep[i];
                m2mp[i] /= sum; m2ip[i] /= sum; m2dp[i] /= sum; m2ep[i] /= sum;
                sum = i2mp[i] + i2ip[i];
                i2mp[i] /= sum; i2ip[i] /= sum;
        }

        for(i=2;i<=shm_len-2;i++) {
                sum = d2mp[i] + d2dp[i];
                d2mp[i] /=sum; d2dp[i] /= sum;
        }

        for(sum = 0.,i=1;i<=shm_len;i++) {
                sum += b2mp[i];
        }

        for(i=1;i<=shm_len;i++) {
                b2mp[i] /= sum;
        }

        nbp = BitsToProb(nb,1.0,bits);
        nnp = BitsToProb(nn,p1,bits);
        sum = nbp + nnp;
        nbp /= sum; nnp /= sum;
        ecp = BitsToProb(ec,1.0,bits);
        ejp = BitsToProb(ej,1.0,bits);
        sum = ecp + ejp;
        ecp /= sum; ejp /= sum;
        ctp = BitsToProb(ct,1.0-p1,bits);
        ccp = BitsToProb(cc,p1,bits);
        sum = ctp + ccp;
        ctp /= sum; ccp /= sum;
        jbp = BitsToProb(jb,1.0,bits);
        jjp = BitsToProb(jj,p1,bits);
        sum = jbp + jjp;
        jbp /= sum; jjp /= sum;

printf("nb=%f nn=%f ec=%f ej=%f jb=%f jj=%f ct=%f cc=%f\n",nbp,nnp,ecp,ejp,jbp,jjp,ctp,ccp);
//FIND SHM

        double            **temp_mat_emit,**temp_ins_emit;
        NEWP(temp_mat_emit,shm_len+3,double);
        NEWP(temp_ins_emit,shm_len+3,double);
        for(i=0;i<=shm_len+1;i++) {
                NEW(temp_mat_emit[i],nA+1,double);
                NEW(temp_ins_emit[i],nA+1,double);
        }

        double          **mat_emit_prob = Mat_emit_Prob(), *mat_emit_prob_i;
        double          **ins_emit_prob = Ins_emit_Prob(), *ins_emit_prob_i;
        double          *mat_emit_p, *ins_emit_p;
        NEW(mat_emit_p,nA+1,double);
        NEW(ins_emit_p,nA+1,double);

        for(i=1;i<=shm_len;i++) {
                mat_emit_prob_i = mat_emit_prob[i];
                for(sum = 0.,j=1;j<=nA;j++) {
                        sum += (mat_emit_p[j] = pow(mat_emit_prob_i[j],1./emiss_tmpr));
                }
                for(j=1;j<=nA;j++) {
                        temp_mat_emit[i][j] = mat_emit_p[j]/sum;
                }
                ins_emit_prob_i = ins_emit_prob[i];
                for(sum = 0.,j=1;j<=nA;j++) {
                        sum += (ins_emit_p[j] = pow(ins_emit_prob_i[j],1./emiss_tmpr));
                }
                for(j=1;j<=nA;j++) {
                        temp_ins_emit[i][j] = ins_emit_p[j]/sum;
                }
        }
        for(i=1;i<=shm_len;i++) { free(mat_emit_prob[i]); free(ins_emit_prob[i]); }
        free(mat_emit_prob); free(ins_emit_prob);
        free(mat_emit_p);free(ins_emit_p);

// EMIT SEQ

        Int4            k,pos;
        char            state;
        e_type          *E;
        NEW(E,nseq+1,e_type);
        unsigned char   *seq;
        NEW(seq,shm_len+100000,unsigned char);
        
        for(j=1;j<=nseq;j++) {
                state = 'n'; i = 1; k = 1;
                while(state != 't') {
                        switch(state) {
                                case 'n': if(Sample2(nnp) == 1) {
                                                seq[k++] = Sample(Pfam_back,nA);
                                                state = 'n';
                                          } else state = 'b';
                                          break;
                                case 'b': i = Sample(b2mp,shm_len);
                                          state = 'm'; 
                                          break;
                                case 'm': seq[k++] = Sample(temp_mat_emit[i],nA);
//printf("match=%c\n",AlphaChar(seq[k-1],AB));
                                          if(i == shm_len) { state = 'e'; i = 1; }
                                          else if(i == shm_len - 1) {
                                                pos = Sample3(m2mp[i],m2ip[i]);
                                                if (pos == 1) { state = 'm'; i++; }
                                                else if(pos == 2) state = 'i';
                                                else { state = 'e'; i = 1; }
                                          } else {
                                                pos = Sample4(m2mp[i],m2ip[i],m2dp[i]);
                                                if (pos == 1) { state = 'm'; i++; }
                                                else if(pos == 2) state = 'i';
                                                else if(pos == 3) { state = 'd'; i++; }
                                                else { state = 'e'; i = 1; }
                                          }
                                          break;
                                case 'd': pos = Sample2(d2mp[i]);
                                          if((pos == 1) || (i == shm_len - 1)) state = 'm';
                                          else state = 'd';
                                          i++; 
                                          break;
                                case 'i': seq[k++] = Sample(temp_ins_emit[i],nA);
//printf("inserted=%c\n",AlphaChar(seq[k-1],AB));
                                          pos = Sample2(i2mp[i]);
                                          if(pos == 1) { state = 'm'; i++; }
                                          else state = 'i';
                                          break;
                                case 'e': pos = Sample2(ejp);
                                          if(pos == 1) state = 'j';
                                          else state = 'c';
                                          i = 1;
                                          break;
                                case 'j': pos = Sample2(jjp);
                                          if(pos == 1) { 
                                                seq[k++] = Sample(Pfam_back,nA);
                                                state = 'j';
                                          } else state = 'b';
                                          break;
                                case 'c': pos = Sample2(ccp);
                                          if(pos == 1) {
                                                seq[k++] = Sample(Pfam_back,nA);
                                                state = 'c';
                                          } else state = 't';
                                          break;
                                default: print_error("unknown state");
                        }
                }
                E[j] = MkSeq((char *) "simulated ",k-1,seq);
        }
        for(i=0;i<=shm_len+1;i++) { free(temp_mat_emit[i]);free(temp_ins_emit[i]); }
        free(temp_mat_emit); free(temp_ins_emit);
        free(m2mp);free(m2ip);free(m2dp);free(i2mp);free(i2ip);
        free(d2mp);free(d2dp);free(b2mp);free(m2ep);
        return E;
}


#if 0
e_type *shm_typ::EmitSeqSHM(Int4 nseq, BooLean adjust_ej)
{
	Int4 		i,j,k,pos,nA=nAlpha(AB);
	char 		state;
	e_type 		*E;
	NEW(E,nseq+1,e_type);
	unsigned char 	*seq;
	NEW(seq,shm_len+100000,unsigned char);
	double	**mat_emit_prob = Mat_emit_Prob();
	double	**ins_emit_prob = Ins_emit_Prob();
	double	*m2mp,*m2ip,*m2dp,*i2mp,*i2ip,*d2mp,*d2dp,*b2mp,*m2ep;
	double nbp,nnp,ecp,ejp,ctp,ccp,jbp,jjp;
	double p1 = 350./351;
	NEW(m2mp,shm_len+2,double);
	NEW(m2ip,shm_len+2,double);
	NEW(m2dp,shm_len+2,double);
	NEW(i2mp,shm_len+2,double);
	NEW(i2ip,shm_len+2,double);
	NEW(d2mp,shm_len+2,double);
	NEW(d2dp,shm_len+2,double);
	NEW(b2mp,shm_len+2,double);
	NEW(m2ep,shm_len+2,double);

        Int4 *m2m,*m2i,*m2d,*i2i,*i2m,*d2d,*d2m,*b2m,*m2e;
        m2m = MatToMat(); m2i = MatToIns(); m2d = MatToDel();
        i2i = InsToIns(); i2m = InsToMat(); d2d = DelToDel();
        d2m = DelToMat(); b2m = BegToMat(); m2e = MatToEnd();

	for(i=1;i<=shm_len-1;i++) {
		m2mp[i] = BitsToProb(m2m[i],p1,bits);
		m2ip[i] = BitsToProb(m2i[i],p1,bits);
		m2dp[i] = BitsToProb(m2d[i],1.0,bits);
		i2mp[i] = BitsToProb(i2m[i],p1,bits);
		i2ip[i] = BitsToProb(i2i[i],p1,bits);
		d2mp[i] = BitsToProb(d2m[i],p1,bits);
		d2dp[i] = BitsToProb(d2d[i],1.0,bits);
		b2mp[i] = BitsToProb(b2m[i],1.0,bits);
		m2ep[i] = BitsToProb(m2e[i],1.0,bits);
	}
	double sum;
        m2dp[shm_len-1] = 0.;
        b2mp[shm_len] = BitsToProb(b2m[shm_len],1.0,bits);
        
        for(i=1;i<=shm_len-1;i++) {
//printf("m2mp[%d]=%f m2ip[%d]=%f m2dp[%d]=%f m2ep[%d]=%f \n",i,m2mp[i],i,m2ip[i],i,m2dp[i],i,m2ep[i]);
                sum = m2mp[i] + m2ip[i] + m2dp[i] + m2ep[i]; 
                m2mp[i] /= sum; m2ip[i] /= sum; m2dp[i] /= sum; m2ep[i] /= sum;
        }
        m2ep[shm_len] = 1.;
//printf("m2ep[%d]=%f \n",shm_len,m2ep[shm_len]);
        for(i=2;i<=shm_len-2;i++) {
//printf("d2mp[%d]=%f d2dp[%d]=%f\n",i,d2mp[i],i,d2dp[i]);
                sum = d2mp[i] + d2dp[i];
                d2mp[i] /= sum; d2dp[i] /= sum;
        }
        d2dp[shm_len-1] = 0.; d2mp[shm_len-1] = 1.;
//printf("d2mp[%d]=%f d2dp[%d]=%f\n",shm_len-1,d2mp[shm_len-1],shm_len-1,d2dp[shm_len-1]);
        
        for(i=1;i<=shm_len-1;i++) {
//printf("i2mp[%d]=%f i2ip[%d]=%f\n",i,i2mp[i],i,i2ip[i]);
                sum = i2mp[i] + i2ip[i];
                i2mp[i] /= sum; i2ip[i] /= sum;
        }
        
        for(sum=0.,i=1;i<=shm_len;i++) {
                sum += b2mp[i];
        }
        
        for(i=1;i<=shm_len;i++) {
//printf("b2mp[%d]=%f\n",i,b2mp[i]);
                b2mp[i] /= sum;
        }

	Int4 nb,nn,ec,ej,ct,cc,jb,jj;
	nb = NB(); nn = NN(); 
	ec = EC(); ej = EJ(); 
	ct = CT(); cc = CC();
	jb = JB(); jj = JJ();

        nbp = BitsToProb(nb,1.0,bits);
        nnp = BitsToProb(nn,p1,bits);
        sum = nbp + nnp;
        nbp /= sum; nnp /= sum;
	if(!adjust_ej || (shm_len >= 135)) {
        	ecp = BitsToProb(ec,1.0,bits);
        	ejp = BitsToProb(ej,1.0,bits);
        	sum = ecp + ejp;
        	ecp /= sum; ejp /= sum;
	} else {
		ejp = 1 - 1./(135./shm_len + 1 );
		ecp = 1. - ejp;
	}

        ctp = BitsToProb(ct,1.0-p1,bits);
        ccp = BitsToProb(cc,p1,bits);
        sum = ctp + ccp;
        ctp /= sum; ccp /= sum;  
        jbp = BitsToProb(jb,1.0,bits);
        jjp = BitsToProb(jj,p1,bits);
        sum = jbp + jjp;
        jbp /= sum; jjp /= sum;
        
	char nm[20];

	for(j=1;j<=nseq;j++) {
		state = 'n'; i = 1; k = 1;
		while(state != 't') {
			switch(state) {
				case 'n': if(Sample2(nnp) == 1) {
						seq[k++] = Sample(Pfam_back,nA);
						state = 'n';
					  } else state = 'b';
					  break;
				case 'b': i = Sample(b2mp,shm_len);
					  state = 'm'; 
					  break;
				case 'm': seq[k++] = Sample(mat_emit_prob[i],nA);
					  if(i == shm_len) { state = 'e'; i = 1; }
					  else if(i == shm_len - 1) {
                                                pos = Sample3(m2mp[i],m2ip[i]);
                                                if (pos == 1) { state = 'm'; i++; }
                                                else if(pos == 2) state = 'i';
                                                else { state = 'e'; i = 1; }
					  } else {
					  	pos = Sample4(m2mp[i],m2ip[i],m2dp[i]);
					  	if (pos == 1) { state = 'm'; i++; }
					  	else if(pos == 2) state = 'i';
					  	else if(pos == 3) { state = 'd'; i++; }
					  	else { state = 'e'; i = 1; }
					  }
					  break;
				case 'd': pos = Sample2(d2mp[i]);
					  if((pos == 1) || (i == shm_len - 1)) state = 'm';
					  else state = 'd';
					  i++; 
					  break;
				case 'i': seq[k++] = Sample(ins_emit_prob[i],nA);
					  pos = Sample2(i2mp[i]);
					  if(pos == 1) { state = 'm'; i++; }
					  else state = 'i';
					  break;
				case 'e': pos = Sample2(ejp);
					  if(pos == 1) state = 'j';	
					  else state = 'c';
					  i = 1;
					  break;
				case 'j': pos = Sample2(jjp);
					  if(pos == 1) { 
					  	seq[k++] = Sample(Pfam_back,nA);
						state = 'j';
					  } else state = 'b';
					  break;
				case 'c': pos = Sample2(ccp);
					  if(pos == 1) {
						seq[k++] = Sample(Pfam_back,nA);
						state = 'c';
					  } else state = 't';
					  break;
				default: print_error("unknown state");
			}
		}	
		E[j] = MkSeq((char *) "simulated ",k-1,seq);
	}
	for(i=1;i<=shm_len;i++) { free(mat_emit_prob[i]);free(ins_emit_prob[i]); }
	free(mat_emit_prob); free(ins_emit_prob);
	free(m2mp);free(m2ip);free(m2dp);free(i2mp);free(i2ip);
	free(d2mp);free(d2dp);free(b2mp);free(m2ep);

	return E;
}
#endif


Int4 shm_typ::ScoreLocalSAP(e_type sE, sap_typ sap)
{
	dsp_typ 	dsp = sap->segs;
	Int4 		*starts = dsp->starts,*lens = dsp->lens;
	Int4 		len,numseg = dsp->numseg;
	Int4		q, s;
	Int4		 pos = 0, score = 0;
	unsigned char	*ptr = SeqPtr(sE);

        Int4 *m2m,*m2i,*m2d,*i2i,*i2m,*d2d,*d2m;
        m2m = stpb->MatToMat();
        m2i = stpb->MatToIns();
        m2d = stpb->MatToDel();
        i2i = stpb->InsToIns();                 
        i2m = stpb->InsToMat();
        d2d = stpb->DelToDel();
        d2m = stpb->DelToMat();

	Int4 start_seg,end_seg,k;

	for(k=0;k<=numseg-1;k++) {
		q = starts[2*k];
		s = starts[2*k+1];
		if(q != -1 && s != -1) {
			start_seg = k;
			break;
		} else continue;
	}

        for(k=numseg-1;k>=0;k--) {
                q = starts[2*k];  
                s = starts[2*k+1];
                if(q != -1 && s != -1) {
                        end_seg = k;
                        break;
                } else continue;
        }

	Int4 i, pos_q = starts[2*start_seg]+1, pos_s = starts[2*start_seg+1]+1;

	for(k=start_seg;k<=end_seg;k++) {
		len = lens[k];
		q = starts[2*k];
		s = starts[2*k+1];
		if(q != -1 && s != -1) {
			for(i=1;i<=len-1;i++) {
				score += mat_emit[pos_q][ptr[pos_s]] + m2m[pos_q];
				pos_s++; pos_q++;
			}
			score += mat_emit[pos_q++][ptr[pos_s++]];
		} else if(q == -1 && s != -1) {
			score += m2i[pos_q-1];
			for(i=1;i<=len-1;i++) {
				score += ins_emit[pos_q-1][ptr[pos_s]] + i2i[pos_q-1];
				pos_s++;
			}
			score += ins_emit[pos_q-1][ptr[pos_s]] + i2m[pos_q-1];
			pos_s++;
		} else {
			score += m2d[pos_q-1];
			for(i=1;i<=len-1;i++) {
				score += d2d[pos_q];
				pos_q++;
			}
			score += d2m[pos_q];
			pos_q++;
		}
	}

	return score;
}

shm_typ *shm_typ::RenormalizeSHM()
{
	shm_typ *rn_shm;
	Int4 i,j;
	stpb_typ *rn_stpb = stpb->RenormalizeStpb();
	Int4 **rn_mat_emit, *mat_emit_i, *rn_mat_emit_i;
	NEWP(rn_mat_emit,shm_len+3,Int4);
	Int4 **rn_ins_emit, *ins_emit_i, *rn_ins_emit_i;
	NEWP(rn_ins_emit,shm_len+3,Int4);
	for(i=0;i<=shm_len;i++) { 
		MEW(rn_mat_emit[i],nAlpha(AB)+3,Int4);
		MEW(rn_ins_emit[i],nAlpha(AB)+3,Int4);
	}
	double *nl;
	MEW(nl,nAlpha(AB)+3,double);
	for(j=1;j<=nAlpha(AB);j++) {
		nl[j] = BitsToProb(nule[j],1./nAlpha(AB),bits);
	}

	double *tmp_emiss,sum;
	MEW(tmp_emiss,nAlpha(AB)+2,double);
	for(i=1;i<=shm_len;i++) {
		mat_emit_i = mat_emit[i];
		rn_mat_emit_i = rn_mat_emit[i];
		ins_emit_i = ins_emit[i];
		rn_ins_emit_i = rn_ins_emit[i];
		for(j=0;j<=nAlpha(AB);j++) {
			rn_mat_emit_i[j] = mat_emit_i[j];
			rn_ins_emit_i[j] = ins_emit_i[j];
		}
	}
	free(nl),free(tmp_emiss);
	Int4 *rn_nule;
	MEW(rn_nule,nAlpha(AB)+3,Int4);
	for(j=1;j<=nAlpha(AB);j++) rn_nule[j] = nule[j];

	rn_shm = new shm_typ(NULL, NULL, rn_nule, rn_mat_emit,
			rn_ins_emit,rn_stpb,AB);
	rn_shm->SetSLambda(Slambda);
	rn_shm->SetSU(Su);
	rn_shm->SetLambda(lambda);
	rn_shm->SetU(u);
        rn_shm->SetWthreshold(wthreshold);
        rn_shm->SetUxparameter(uxparameter);
        rn_shm->SetUthreshold(uthreshold);
        rn_shm->SetXparameter(xparameter);
        rn_shm->SetGthreshold(gthreshold);
        rn_shm->SetVthreshold(vthreshold);
	return rn_shm;
}

sap_typ AlignToSAP(e_type qE, e_type sE, char *operation, Int4 start_qE, Int4 start_sE, Int4 score)
{
	if(operation[1]=='E') return NULL;
        Int4 i=1,k,j,s,t,seg;
        short numseg = 0;
        while(operation[i] != 'E'){
                if(tolower(operation[i]) != tolower(operation[i-1])) numseg++;
                i++;
        }
        Int4 *starts, *lens;
        MEW(starts,2*numseg+1,Int4); MEW(lens,numseg+2,Int4);
        i=2;k=1;seg=0;
        while(operation[i] != 'E'){
                if(tolower(operation[i]) != tolower(operation[i-1])) {
                        lens[seg++] = k;
                        k = 1;
                } else { k++; }
                i++;
        }
        lens[seg] = k;
        k=1; s=1; t=1; j=0;
        for(i=0;i<numseg;i++){
                switch(operation[t]){
                        case 'd': 
                        case 'D': starts[j++] = -1;
				  starts[j++] = k+start_qE-2;
                                  k += lens[i];
                                  break;
                        case 'm': 
                        case 'M': starts[j++] = s+start_sE-2;
				  starts[j++] = k+start_qE-2;
                                  k += lens[i]; s += lens[i];                             
                                  break;
                        case 'I': starts[j++] = s+start_sE-2;
				  starts[j++] = -1;
                                  s += lens[i];
                                  break;
                }
                t += lens[i];
        }
        sap_typ sap = MakeGSeqAlign(numseg, SeqI(sE), SeqI(qE), sE, qE, starts, lens);
	sap->score = score;
        return sap;
}
#if 0
sap_typ AlignToSAP2(e_type qE, e_type sE, char *operation, Int4 start_qE, Int4 start_sE, Int4 score)
{
	if(operation[1]=='E') return NULL;
        Int4 i=1,k,j,s,t,seg;
        short numseg = 0;
        while(operation[i] != 'E'){
                if(tolower(operation[i]) != tolower(operation[i-1])) numseg++;
                i++;
        }
        Int4 *starts, *lens;
        MEW(starts,2*numseg+1,Int4);
        MEW(lens,numseg+2,Int4);
        i=2;k=1;seg=0;
        while(operation[i] != 'E'){
                if(tolower(operation[i]) != tolower(operation[i-1])) {
                        lens[seg++] = k;
                        k = 1;
                } else { k++; }
                i++;
        }
        lens[seg] = k;
        k=1; s=1; t=1; j=0;
        for(i=0;i<numseg;i++){
                switch(operation[t]){
                        case 'd': 
                        case 'D': starts[j++] = k+start_qE-2;
                                  starts[j++] = -1;
                                  k += lens[i];
                                  break;
                        case 'm': 
                        case 'M': starts[j++] = k+start_qE-2;
                                  starts[j++] = s+start_sE-2;
                                  k += lens[i]; s += lens[i];                             
                                  break;
                        case 'I': starts[j++] = -1;
                                  starts[j++] = s+start_sE-2;
                                  s += lens[i];
                                  break;
                }
                t += lens[i];
        }
        sap_typ sap = MakeGSeqAlign(numseg, SeqI(qE), SeqI(sE), qE, sE, starts, lens);
	sap->score = score;
        return sap;
}
#endif

double shm_typ::ExpectedCNJScore()
{
	double escore=0.;
	Int4 i;
	for(i=1;i<=nAlpha(AB);i++) {
		escore += BitsToProb(nule[i],1./nAlpha(AB),bits)*nule[i];
	}
	return escore;
}

double shm_typ::ExpectedNScore()
{
	double escore=0.,sum;
	double p1=350./351.;

	Int4 n2n = NN(), n2b = NB();
	double n2np, n2bp;

	n2np = BitsToProb(n2n,p1,bits);
	n2bp = BitsToProb(n2b,1.0,bits);
	sum = n2np + n2bp;
	n2np /= sum; n2bp /= sum;
	double CNJ = ExpectedCNJScore();
//printf("n2np=%f n2bp=%f CNJ=%f n2b/(1.-n2np)=%f
//n2n*n2np/((1.-n2np)*(1.-n2np))=%f\n",n2np,n2bp,CNJ,n2b/(1.-n2np),
//n2n*n2np/((1.-n2np)*(1.-n2np)));
	escore = n2bp*(n2b/(1.-n2np) + n2n*n2np/((1.-n2np)*(1.-n2np)));
	return escore;
}

double shm_typ::ExpectedJScore()
{
	double escore=0.,sum;
	double p1=350./351.;

	Int4 j2b = JB(), j2j = JJ();
	double j2bp,j2jp;

	double CNJ = ExpectedCNJScore();
	j2bp = BitsToProb(j2b,1.0,bits);
	j2jp = BitsToProb(j2j,p1,bits);
	sum = j2bp + j2jp;
	j2bp /= sum; j2jp /= sum;
	escore = j2bp*(j2b/(1.-j2jp) + j2j*j2jp/((1.-j2jp)*(1.-j2jp)));
	return escore;
}

double shm_typ::ExpectedCScore()
{
	double escore=0.;
	double p1=350./351.;

	Int4 c2t = CT(), c2c = CC();
	double c2tp,c2cp;

	double CNJ = ExpectedCNJScore();
	c2tp = BitsToProb(c2t,1.0-p1,bits);
	c2cp = BitsToProb(c2c,p1,bits);
	escore = c2tp*(c2t/(1.-c2cp) + c2c*c2cp/((1.-c2cp)*(1.-c2cp)));
	return escore;
}

double shm_typ::ExpectedMatchStateScore(Int4 pos)
{
	assert(pos>0 && pos<=shm_len);
	double escore=0.;
	Int4 i;
	double *mat_emit_prob = Mat_emit_Prob(pos);
	Int4   *mat_emit_pos = mat_emit[pos];
	for(i=1;i<=nAlpha(AB);i++) {
		escore += mat_emit_prob[i]*mat_emit_pos[i];
	}
	return escore;
}

double shm_typ::ExpectedMatchScore(Int4 pos)
{
	double escore=0.;
	double p1=350./351.;

	Int4 *m2e,*m2m,*m2i,*m2d;
	double m2ep,m2mp,m2ip,m2dp;

	if (pos == shm_len-1){
		m2e = MatToEnd();
		m2m = MatToMat();
		m2i = MatToIns();
		m2ep = BitsToProb(m2e[pos],1.0,bits);
		m2mp = BitsToProb(m2m[pos],p1,bits);
		m2ip = BitsToProb(m2i[pos],p1,bits);
		escore = m2ep*m2e[pos]+m2mp*m2m[pos]+m2ip*m2i[pos];
	} else if(pos == shm_len){
		m2e = MatToEnd();
		m2ep = BitsToProb(m2e[pos],1.0,bits);
		escore = m2ep*m2e[pos];
	} else {
		m2e = MatToEnd();
		m2m = MatToMat();
		m2i = MatToIns();
		m2d = MatToDel();
		m2ep = BitsToProb(m2e[pos],1.0,bits);
		m2mp = BitsToProb(m2m[pos],p1,bits);
		m2ip = BitsToProb(m2i[pos],p1,bits);
		m2dp = BitsToProb(m2d[pos],1.0,bits);
		escore = m2ep*m2e[pos]+m2mp*m2m[pos]+m2ip*m2i[pos]+m2dp*m2d[pos];
	}
	escore += ExpectedMatchStateScore(pos);
	return escore;
}

double shm_typ::ExpectedInsertionStateScore(Int4 pos)
{
	assert(pos>0 && pos<=shm_len);
	double escore=0.;
	Int4 i;
	double *ins_emit_prob = Ins_emit_Prob(pos);
	Int4   *ins_emit_pos = ins_emit[pos];
	for(i=1;i<=nAlpha(AB);i++) {
		escore += ins_emit_prob[i]*ins_emit_pos[i];
	}
	return escore;
}

double shm_typ::ExpectedInsertionScore(Int4 pos)
{
	double escore=0.;
	double p1=350./351.;

	Int4 *i2m,*i2i;
	double i2mp,i2ip;
	double expIns = ExpectedInsertionStateScore(pos);

	i2m = InsToMat();
	i2i = InsToIns();
	i2mp = BitsToProb(i2m[pos],p1,bits);
	i2ip = BitsToProb(i2i[pos],p1,bits);
	escore = i2m[pos] + expIns*(1 + i2ip/i2mp) + i2i[pos]*i2ip*(1 + i2ip/i2mp);
	return escore;
}

double shm_typ::ExpectedDeletionScore(Int4 pos)
{
	assert(pos>1 && pos<shm_len);
	double escore=0.;
	double p1=350./351.;
	Int4 i;

	Int4 *d2m = DelToMat();
	Int4 *d2d = DelToDel();

	double d2mp = BitsToProb(d2m[pos],p1,bits);
	double d2dp = BitsToProb(d2d[pos],1.0,bits);
	if(pos == shm_len-1) d2dp = 0.;

	escore = d2mp*d2m[pos] + d2dp*d2d[pos];

	return escore;
}

double shm_typ::ExpectedBeginScore()
{
	double escore=0.;
	Int4 i;
	Int4 *b2m = BegToMat();

	for(i=1;i<=shm_len;i++) {
		escore += BitsToProb(b2m[i],1.0,bits)*b2m[i];
	}

	return escore;
}

double shm_typ::ExpectedMainModelScore()
{
	double 	escore=0.,sum;
	Int4 	i,im1;
	double	*m2mp,*m2ip,*m2dp,*i2mp,*i2ip,*d2mp,*d2dp,*b2mp,*m2ep;
	double p1 = 350./351;
	NEW(m2mp,shm_len+2,double);
	NEW(m2ip,shm_len+2,double);
	NEW(m2dp,shm_len+2,double);
	NEW(i2mp,shm_len+2,double);
	NEW(i2ip,shm_len+2,double);
	NEW(d2mp,shm_len+2,double);
	NEW(d2dp,shm_len+2,double);
	NEW(b2mp,shm_len+2,double);
	NEW(m2ep,shm_len+2,double);

        Int4 *m2m,*m2i,*m2d,*i2i,*i2m,*d2d,*d2m,*b2m,*m2e;
        m2m = MatToMat(); m2i = MatToIns(); m2d = MatToDel();
        i2i = InsToIns(); i2m = InsToMat(); d2d = DelToDel();
        d2m = DelToMat(); b2m = BegToMat(); m2e = MatToEnd();

	for(i=1;i<=shm_len-1;i++) {
		m2mp[i] = BitsToProb(m2m[i],p1,bits);
		m2ip[i] = BitsToProb(m2i[i],p1,bits);
		m2dp[i] = BitsToProb(m2d[i],1.0,bits);
		i2mp[i] = BitsToProb(i2m[i],p1,bits);
		i2ip[i] = BitsToProb(i2i[i],p1,bits);
		d2mp[i] = BitsToProb(d2m[i],p1,bits);
		d2dp[i] = BitsToProb(d2d[i],1.0,bits);
		b2mp[i] = BitsToProb(b2m[i],1.0,bits);
		m2ep[i] = BitsToProb(m2e[i],1.0,bits);
	}
	m2dp[shm_len-1] = 0.;
	b2mp[shm_len] = BitsToProb(b2m[shm_len],1.0,bits);

	for(i=1;i<=shm_len-1;i++) {
//printf("m2mp[%d]=%f m2ip[%d]=%f m2dp[%d]=%f m2ep[%d]=%f \n",i,m2mp[i],i,m2ip[i],i,m2dp[i],i,m2ep[i]);
		sum = m2mp[i] + m2ip[i] + m2dp[i] + m2ep[i];
		m2mp[i] /= sum; m2ip[i] /= sum; m2dp[i] /= sum; m2ep[i] /= sum;
	}
	m2ep[shm_len] = 1.;
//printf("m2ep[%d]=%f \n",shm_len,m2ep[shm_len]);
	for(i=2;i<=shm_len-2;i++) {
//printf("d2mp[%d]=%f d2dp[%d]=%f\n",i,d2mp[i],i,d2dp[i]);
		sum = d2mp[i] + d2dp[i];
		d2mp[i] /= sum; d2dp[i] /= sum;
	}
	d2dp[shm_len-1] = 0.; d2mp[shm_len-1] = 1.; d2dp[1] = 0.;
//printf("d2mp[%d]=%f d2dp[%d]=%f\n",shm_len-1,d2mp[shm_len-1],shm_len-1,d2dp[shm_len-1]);

	for(i=1;i<=shm_len-1;i++) {
//printf("i2mp[%d]=%f i2ip[%d]=%f\n",i,i2mp[i],i,i2ip[i]);
		sum = i2mp[i] + i2ip[i];
		i2mp[i] /= sum; i2ip[i] /= sum;
	}

	for(sum=0.,i=1;i<=shm_len;i++) {
		sum += b2mp[i];
	}

	for(i=1;i<=shm_len;i++) {
//printf("b2mp[%d]=%f\n",i,b2mp[i]);
		b2mp[i] /= sum;
	}

	double *gM,*gD,*gI,*gtempI;
        NEW(gM,shm_len+3,double);NEW(gI,shm_len+3,double);NEW(gD,shm_len+3,double);
	NEW(gtempI,shm_len+3,double);
	gM[1] = b2mp[1];
	gD[1] = 0.;
	gtempI[1] = m2ip[1]*gM[1]/(1.-i2ip[1]);
	gI[1] = m2ip[1]*gM[1];	
        for(im1=1,i=2;i<=shm_len-1;im1++,i++) {
                gM[i] = b2mp[i] + d2mp[im1]*gD[im1] + i2mp[im1]*gtempI[im1] + m2mp[im1]*gM[im1];
                gD[i] = d2dp[im1]*gD[im1] + m2dp[im1]*gM[im1];
		gtempI[i] = m2ip[i]*gM[i]/(1.-i2ip[i]);
                gI[i] = m2ip[i]*gM[i];
        }

	gM[shm_len] = b2mp[shm_len] + d2mp[shm_len-1]*gD[shm_len-1] + i2mp[shm_len-1]*gtempI[shm_len-1] + 
				m2mp[shm_len-1]*gM[shm_len-1];
	gD[shm_len] = 0.;
	gI[shm_len] = 0.;

	for(i=1;i<=shm_len;i++) {
		escore += ExpectedMatchScore(i)*gM[i]; 
		printf("gM[%d]=%f ExpectedMatchScore[%d]=%f contr=%f\n",
			i,gM[i],i,ExpectedMatchScore(i),ExpectedMatchScore(i)*gM[i]);
	}
	for(i=1;i<shm_len;i++) {
		escore += ExpectedInsertionScore(i)*gI[i]; 
		printf("gI[%d]=%f ExpectedInsertionScore[%d]=%f contr=%f\n",
			i,gI[i],i,ExpectedInsertionScore(i),ExpectedInsertionScore(i)*gI[i]);
	}
	for(i=2;i<shm_len;i++) {
		escore += ExpectedDeletionScore(i)*gD[i]; 
		printf("gD[%d]=%f ExpectedDeletionScore[%d]=%f contr=%f\n",
			i,gD[i],i,ExpectedDeletionScore(i),ExpectedDeletionScore(i)*gD[i]);
	}
	escore += ExpectedBeginScore();
	
	free(gM);free(gI);free(gD);
	return escore;
}

double shm_typ::ExpectedScore()
{
	double escore=0.;
	double sum;
        Int4 nb,nn,ec,ej,ct,cc,jb,jj;
        nb = NB(); nn = NN(); ec = EC();
        ej = EJ(); ct = CT(); cc = CC();
        jb = JB(); jj = JJ();
               
	double p1 = 350./351.;
 
        double nbp = BitsToProb(nb,1.0,bits);
        double nnp = BitsToProb(nn,p1,bits);
	sum = nbp + nnp;
	nbp /= sum; nnp /= sum;
        double ecp = BitsToProb(ec,1.0,bits);
        double ejp = BitsToProb(ej,1.0,bits);
        sum = ecp + ejp;
        ecp /= sum; ejp /= sum;
        double ctp = BitsToProb(ct,1.0-p1,bits);
        double ccp = BitsToProb(cc,p1,bits);
        sum = ctp + ccp;
        ctp /= sum; ccp /= sum;
        double jbp = BitsToProb(jb,1.0,bits);
        double jjp = BitsToProb(jj,p1,bits);
        sum = jbp + jjp;
        jbp /= sum; jjp /= sum;

printf("nb=%d nn=%d ec=%d ej=%d ct=%d cc=%d jb=%d jj=%d\n",nb,nn,ec,ej,ct,cc,jb,jj);
printf("nbp=%f nnp=%f ecp=%f ejp=%f ctp=%f ccp=%f jbp=%f jjp=%f\n",nbp,nnp,ecp,ejp,ctp,ccp,jbp,jjp);

	double mscore = ExpectedMainModelScore();
	double N = ExpectedNScore();	
	double C = ExpectedCScore();
	double J = ExpectedJScore();
	double bscore = N;
	double tscore = C;

	double rscore = (mscore+ec)*ecp + 
			mscore*ecp*ejp*(1./((1.-ejp)*(1.-ejp))+1./(1.-ejp)) + 
			(ej+J)*ecp*ejp/((1.-ejp)*(1.-ejp)) +
			ec*ecp*ejp*1./(1.-ejp);
printf("first=%f sec=%f third=%f fourth=%f\n",(mscore+ec)*ecp,mscore*ecp*ejp*(1./((1.-ejp)*(1.-ejp))+1./(1.-ejp)),
		(ej+J)*ecp*ejp/((1.-ejp)*(1.-ejp)), ec*ecp*ejp*1./(1.-ejp));


printf("rscore=%f mscore=%f bscore=%f tscore=%f C=%f N=%f\n",rscore,mscore,bscore,tscore,C,N);
	escore = rscore + bscore + tscore;
	return escore;
}

Int4 *shm_typ::RunningScore(dsp_typ dsp, Int4 *length, char stype)
// if stype=S return running score at each step including the emiss at that step
// if stype=W without the last emiss
// if stype=B without any transitions (like BLAST)
{
        Int4	*rscore;
        Int4	i,k,len,pos_s,pos_p;
	Int4	s,p,prev_s,prev_p;
	Int4	lng,br=1,score=0,wo_emiss_score=0,wo_trans_score=0;
        Int4 seq_len = LenSeq(dsp->qE); 
        Int4 numseg = dsp->numseg;
        unsigned char *ptr = SeqPtr(dsp->qE);

        Int4 *m2m,*m2i,*m2d,*i2i,*i2m,*d2d,*d2m,*b2m,*m2e;
        m2m = MatToMat(); m2i = MatToIns(); m2d = MatToDel();
        i2i = InsToIns(); i2m = InsToMat(); d2d = DelToDel();
        d2m = DelToMat(); b2m = BegToMat(); m2e = MatToEnd();
        
        Int4 *starts = dsp->starts;
        Int4 *lens = dsp->lens;

        for(lng=0,i=0;i<numseg;i++) { lng += dsp->lens[i]; }
        NEW(rscore,lng+2,Int4);
        *length = lng;

        len = lens[0];
        s = prev_s = starts[0];
        p = prev_p = starts[1];

        if(s != -1 && p != -1) { 
                pos_s = s+1; pos_p = p+1;
		score = mat_emit[pos_p][ptr[pos_s]];
		wo_emiss_score = 0;
		wo_trans_score = score;
		if(stype != 'W') rscore[br++] = score;
		else rscore[br++] = wo_emiss_score;
                pos_s++; pos_p++;
                for(i=2;i<=len;i++) {
                        score += (mat_emit[pos_p][ptr[pos_s]] + m2m[pos_p-1]);
			wo_emiss_score = score - mat_emit[pos_p][ptr[pos_s]];
			wo_trans_score += mat_emit[pos_p][ptr[pos_s]];
			if(stype == 'S') rscore[br++] = score;
			else if(stype == 'W') rscore[br++] = wo_emiss_score;
			else rscore[br++] = wo_trans_score;
                        pos_s++; pos_p++;
                }
        } else if(s == -1 && p != -1) { 
                pos_p = p+1;
                score += m2d[pos_p-1];
		wo_emiss_score = score;
                if(stype == 'S' || stype == 'W') rscore[br++] = score;
                else rscore[br++] = wo_trans_score;
                pos_p++;
                for(i=2;i<=len;i++) {
                                score += d2d[pos_p-1];
				wo_emiss_score = score;
				if(stype == 'S' || stype == 'W') rscore[br++] = score;
				else rscore[br++] = wo_trans_score; 
                                pos_p++;
                }
        } else {
                pos_s = s+1; pos_p = starts[3];
                score += (ins_emit[pos_p][ptr[pos_s]] + m2i[pos_p]);
		wo_emiss_score = score - ins_emit[pos_p][ptr[pos_s]];		
		if(stype == 'S') rscore[br++] = score;
		else if(stype == 'W') rscore[br++] = wo_emiss_score;
		else rscore[br++] = wo_trans_score;
                pos_s++;
                for(i=2;i<=len;i++) {
                        score += (ins_emit[pos_p][ptr[pos_s]] + i2i[pos_p]);
	                wo_emiss_score = score - ins_emit[pos_p][ptr[pos_s]];
        	        if(stype == 'S') rscore[br++] = score;
                	else if(stype == 'W') rscore[br++] = wo_emiss_score;
                	else rscore[br++] = wo_trans_score;
                        pos_s++;
                }
        }

        for(k=1;k<numseg;k++) {
                len = lens[k]; 
                s = starts[2*k];
                p = starts[2*k+1];
                if(p != -1 && s != -1) {
                        if (prev_s == -1) { 
                                score += (mat_emit[pos_p][ptr[pos_s]] + d2m[pos_p-1]);
				wo_emiss_score = score - mat_emit[pos_p][ptr[pos_s]];
				wo_trans_score += mat_emit[pos_p][ptr[pos_s]];
				if(stype == 'S') rscore[br++] = score;
				else if(stype == 'W') rscore[br++] = wo_emiss_score;
				else rscore[br++] = wo_trans_score;
                                pos_p++; pos_s++;
                        } else {
                                score += (mat_emit[pos_p][ptr[pos_s]] + i2m[pos_p-1]);
				wo_emiss_score = score - mat_emit[pos_p][ptr[pos_s]];
				wo_trans_score += mat_emit[pos_p][ptr[pos_s]];
				if(stype == 'S') rscore[br++] = score;
				else if(stype == 'W') rscore[br++] = wo_emiss_score;
				else rscore[br++] = wo_trans_score;
                                pos_p++; pos_s++;
                        }
                        for(i=2;i<=len;i++) {
                                score += (mat_emit[pos_p][ptr[pos_s]] + m2m[pos_p-1]);
				wo_emiss_score = score - mat_emit[pos_p][ptr[pos_s]];
				wo_trans_score += mat_emit[pos_p][ptr[pos_s]];
				if(stype == 'S') rscore[br++] = score;
				else if(stype == 'W') rscore[br++] = wo_emiss_score;
				else rscore[br++] = wo_trans_score;
                                pos_p++; pos_s++;
                        }
                } else if(p == -1 && s != -1) {
                        score += (ins_emit[pos_p-1][ptr[pos_s]] + m2i[pos_p-1]);
			wo_emiss_score = score - ins_emit[pos_p-1][ptr[pos_s]];
			if(stype == 'S') rscore[br++] = score;
			else if(stype == 'W') rscore[br++] = wo_emiss_score;
			else rscore[br++] = wo_trans_score;
                        pos_s++;
                        for(i=2;i<=len;i++) {
                                score += (ins_emit[pos_p-1][ptr[pos_s]] + i2i[pos_p-1]);
				wo_emiss_score = score - ins_emit[pos_p-1][ptr[pos_s]];
				if(stype == 'S') rscore[br++] = score;
				else if(stype == 'W') rscore[br++] = wo_emiss_score;
				else rscore[br++] = wo_trans_score;
                                pos_s++;
                        }
                } else {
                        score += m2d[pos_p-1];
			wo_emiss_score = score;
			if(stype == 'S' || stype == 'W') rscore[br++] = score;
			else rscore[br++] = wo_trans_score;
                        pos_p++;
                        for(i=2;i<=len;i++) {
                                score += d2d[pos_p];
				wo_emiss_score = score;
	                        if(stype == 'S' || stype == 'W') rscore[br++] = score;
        	                else rscore[br++] = wo_trans_score;
                                pos_p++;
                        }
                }
                prev_p = p; prev_s = s;
        }
        return rscore;
}

