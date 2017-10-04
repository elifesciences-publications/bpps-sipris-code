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

#include "stpb_typ.h"

#if 1   //************ Make public stpb_typ routines eventually ************
static  Int4 **Concat(Int4 **smth, Int4 smth_len, Int4 smth_width, Int4 rpts)
{
        Int4 **rpt_smth, *rpt_smth_temp, *smth_i;
        Int4 i,j,k;
        NEWP(rpt_smth,rpts*smth_len+3,Int4);
        for(i=0;i<=rpts*smth_len;i++) NEW(rpt_smth[i],smth_width+3,Int4);
        for(j=0;j<rpts;j++) {
                for(i=1;i<=smth_len;i++) {
                        rpt_smth_temp = rpt_smth[j*smth_len+i];
                        smth_i = smth[i];
                        for(k=0;k<=smth_width;k++) {
                                rpt_smth_temp[k] = smth_i[k];
                        }
                }
        }
        for(i=0;i<=smth_len;i++) free(smth[i]); 
                free(smth);
        return rpt_smth;
}

static Int4 *Concat(Int4 *smth,Int4 smth_len,Int4 rpts)
{
        Int4 i,j;
        Int4 *rpt_smth;
        NEW(rpt_smth,rpts*smth_len+3,Int4);
        rpt_smth[0] = smth[0];
        for(i=1;i<=smth_len;i++) {
                for(j=0;j<rpts;j++) {
                        rpt_smth[j*smth_len+i] = smth[i];
                }
        }
        free(smth);
        return rpt_smth;
}

#endif

stpb_typ::stpb_typ(Int4 Length, double *m2M, double *m2I, double *m2D, double *i2I, 
	double *i2M, double *d2D, double *d2M, double *b2M, double *m2E, double nB, double nN, double eC, 
	double eJ, double cT, double cC, double jB, double jJ, double bM, double bD, double gG, double gF,
	Int4 Bits)
{
printf("nB=%f\n",nB);
        Int4 NN = ProbToBits(nN,1.0,Bits);
        Int4 NB = ProbToBits(nB,1.0,Bits);
        Int4 EJ = ProbToBits(eJ,1.0,Bits);
        Int4 EC = ProbToBits(eC,1.0,Bits);
        Int4 CC = ProbToBits(cC,1.0,Bits);
        Int4 CT = ProbToBits(cT,1.0,Bits);
        Int4 JJ = ProbToBits(jJ,1.0,Bits);
        Int4 JB = ProbToBits(jB,1.0,Bits);
printf("stpb->NB=%d\n",NB);

	Int4 BM = ProbToBits(bM,1.0,Bits);      
	Int4 BD = ProbToBits(bD,1.0,Bits);      
	Int4 GG = ProbToBits(gG,1.0,Bits);      
	Int4 GF = ProbToBits(gF,1.0,Bits);      

	Int4 *M2M,*M2D,*M2I,*I2M,*I2I,*D2M,*D2D,*B2M,*M2E;
	NEW(M2M,Length+1,Int4); NEW(M2I,Length+1,Int4); NEW(M2D,Length+1,Int4);
	NEW(I2M,Length+1,Int4); NEW(I2I,Length+1,Int4); NEW(D2M,Length+1,Int4);
	NEW(D2D,Length+1,Int4); NEW(B2M,Length+1,Int4); NEW(M2E,Length+1,Int4);
        for(Int4 i=1;i<=Length;i++) {
                M2M[i] = ProbToBits(m2M[i],1.0,Bits);
                M2D[i] = ProbToBits(m2D[i],1.0,Bits);
                M2I[i] = ProbToBits(m2I[i],1.0,Bits);
                I2M[i] = ProbToBits(i2M[i],1.0,Bits);
                I2I[i] = ProbToBits(i2I[i],1.0,Bits);
                D2M[i] = ProbToBits(d2M[i],1.0,Bits);
                D2D[i] = ProbToBits(d2D[i],1.0,Bits);
                B2M[i] = ProbToBits(b2M[i],1.0,Bits);
                M2E[i] = ProbToBits(m2E[i],1.0,Bits);
        }
printf("m2E[1]=%f M2E[1]=%d\n",m2E[1],M2E[1]);
	Init(Length,M2M,M2I,M2D,I2I,I2M,D2D,D2M,B2M,M2E,NB,NN,EC,EJ,CT,CC,JB,JJ,BM,BD,
                GG,GF,Bits);
}

stpb_typ::stpb_typ(Int4 Length, Int4 *M2M, Int4 *M2I, Int4 *M2D, Int4 *I2I, 
	Int4 *I2M, Int4 *D2D, Int4 *D2M, double	FreqOfMatchAtPositionOne)
// AFN: added to allow NAR manuscript profiles to be entered...
// Uses default HMMER save file values as much as possible...
//
// B2M & M2E are all really log(prob) not likelihoods!
// So they always need to be negative!
{
	Int4 NB=-8455,NN=-4,EC=-1000,EJ=-1000,CT=-8455,CC=-4,JB=-8455,JJ=-4;
	// XT line parameters in save file... EJ==-1000 implies 50% more than one repeat
	// These correspond to transitions within Plan7 architecture.
	// NOTE that prob == 0 corresponds to -8455 and prob == 1 corresponds to -4.
	// Save file:
	// m->m   m->i   m->d   i->m   i->i   d->m   d->d   b->m   m->e
	//	-4      *    -8455
	//     B->M          B->E  at first position...
	         
	Int4 GG=-4,GF=-8455;	// NULT line = null model parameters.
	Int4 Bits=1000;
       	Int4 BD=BigNegative;	// for global models...
		// as this approaches zero the models becomes more & more local.
	Int4 *B2M, *M2E; // These are all BigNegative for global models.
	Int4 BM;	// 
	// Set BM to be the 
	assert(FreqOfMatchAtPositionOne > 0.0 && FreqOfMatchAtPositionOne <= 1.0);
	BM = (Int4) floor(Bits*log2(FreqOfMatchAtPositionOne));
	// A reasonable value 

	NEW(B2M,Length+3,Int4);
	NEW(M2E,Length+3,Int4);
	// WARNING: this is not freed up! Create flag to tell us to free it...
	for(Int4 i=0; i <= Length; i++) B2M[i]=M2E[i]=BigNegative;

	Init(Length,M2M,M2I,M2D,I2I,I2M,D2D,D2M,B2M,M2E,NB,NN,EC,EJ,CT,CC,JB,JJ,BM,BD,
		GG,GF,Bits);
}

stpb_typ::stpb_typ(Int4 Length, Int4 *M2M, Int4 *M2I, Int4 *M2D, Int4 *I2I, 
	Int4 *I2M, Int4 *D2D, Int4 *D2M, Int4 *B2M, Int4 *M2E, Int4 NB, Int4 NN, Int4 EC, 
	Int4 EJ, Int4 CT, Int4 CC, Int4 JB, Int4 JJ, Int4 BM, Int4 BD, Int4 GG, Int4 GF, 
	Int4 Bits)
{
	Init(Length,M2M,M2I,M2D,I2I,I2M,D2D,D2M,B2M,M2E,NB,NN,EC,EJ,CT,CC,JB,JJ,BM,BD,
		GG,GF,Bits);
}

stpb_typ::stpb_typ(stpb_typ *Stpb,Int4 rpts)
// concatenating stpb - needs work
{
	Int4 *M2M = Concat(Stpb->MatToMat(),Stpb->Length(),rpts);
	Int4 *M2I = Concat(Stpb->MatToIns(),Stpb->Length(),rpts);
	Int4 *M2D = Concat(Stpb->MatToDel(),Stpb->Length(),rpts);
	Int4 *I2I = Concat(Stpb->InsToIns(),Stpb->Length(),rpts);
	Int4 *I2M = Concat(Stpb->InsToMat(),Stpb->Length(),rpts);
	Int4 *D2D = Concat(Stpb->DelToDel(),Stpb->Length(),rpts);
	Int4 *D2M = Concat(Stpb->DelToMat(),Stpb->Length(),rpts);
	Int4 *B2M = Concat(Stpb->BegToMat(),Stpb->Length(),rpts);
	Int4 *M2E = Concat(Stpb->MatToEnd(),Stpb->Length(),rpts);

	Init(rpts*Stpb->Length(),M2M,M2I,M2D,I2I,I2M,D2D,D2M,B2M,M2E,Stpb->NB(),Stpb->NN(),
		Stpb->EC(),Stpb->EJ(),Stpb->CT(),Stpb->CC(),Stpb->JB(),Stpb->JJ(),
			Stpb->BM(),Stpb->BD(),Stpb->GG(),Stpb->GF(),Stpb->Bits());
}

stpb_typ& stpb_typ::operator=(const stpb_typ& stpb)
{
	if (this != &stpb) { 
		Free(); 
		Init(stpb.length,stpb.m2m,stpb.m2i,stpb.m2d,stpb.i2i,stpb.i2m,stpb.d2d,stpb.d2m,
			stpb.b2m,stpb.m2e,stpb.nb,stpb.nn,stpb.ec,stpb.ej,stpb.ct,stpb.cc,stpb.jb,
				stpb.jj,stpb.bm,stpb.bd,stpb.gg,stpb.gf,stpb.bits);
	} return *this;
}

void stpb_typ::Init(Int4 Length, Int4 *M2M, Int4 *M2I, Int4 *M2D, Int4 *I2I, 
        Int4 *I2M, Int4 *D2D, Int4 *D2M, Int4 *B2M, Int4 *M2E, Int4 NB, Int4 NN, Int4 EC, Int4 EJ, 
		Int4 CT, Int4 CC, Int4 JB, Int4 JJ, Int4 BM, Int4 BD, Int4 GG, Int4 GF, Int4 Bits)
{
	length = Length;
	m2m = M2M; m2i = M2I; m2d = M2D;
	i2i = I2I;i2m = I2M; d2d = D2D; d2m = D2M;
	b2m = B2M; m2e = M2E;
	nb = NB; nn = NN; ec = EC; ej = EJ;
	ct = CT; cc = CC; jb = JB; jj = JJ; 
	bm = BM; bd = BD;
	gg = GG; gf = GF; 
	bits = Bits;
}

void	stpb_typ::Free( ) 
{
	if(m2m != NULL) free(m2m);
	if(m2i != NULL) free(m2i);
	if(m2d != NULL) free(m2d);
	if(i2i != NULL) free(i2i);
	if(i2m != NULL) free(i2m);
	if(d2d != NULL) free(d2d);
	if(d2m != NULL) free(d2m);
	if(b2m != NULL) free(b2m);
	if(m2e != NULL) free(m2e);
}

void	stpb_typ::Put(FILE *fp)
{
        fprintf(fp,"POS:  m2m  m2i  m2d  i2i  i2m  d2d  d2m  b2m  m2e\n");
        for(Int4 i=0;i<=length;i++){
            fprintf(fp,"%-3ld: %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld\n",
                       i,-m2m[i],-m2i[i],-m2d[i],-i2i[i],-i2m[i],-d2d[i],-d2m[i],-b2m[i],-m2e[i]);
        }
}

stpb_typ *stpb_typ::ReverseStpb()
{
	register Int4 i;
	register Int4 im1;
	stpb_typ *reversed_stpb;

	double *pm2m,*pm2i,*pm2d,*pi2i,*pi2m,*pd2d,*pd2m;
        double *gM2m,*gM2i,*gM2d,*gI2m,*gD2d,*gD2m;
        double *RgM2m,*RgM2i,*RgM2d,*RgI2m,*RgD2d,*RgD2m;
        Int4 *rM2m,*rM2i,*rM2d,*rI2i,*rI2m,*rD2d,*rD2m;
	Int4 *frM2m,*frM2i,*frM2d,*frI2i,*frI2m,*frD2d,*frD2m,*frB2m,*frM2e;

	NEW(pm2m,length+3,double);NEW(pm2i,length+3,double);
	NEW(pm2d,length+3,double);NEW(pi2i,length+3,double);
	NEW(pi2m,length+3,double);NEW(pd2d,length+3,double);
	NEW(pd2m,length+3,double);


	double sum,a,b,c;
	a = BitsToProb(m2m[0],1.0,bits);
	b = BitsToProb(m2d[0],1.0,bits);
	sum = a+b;
	assert(sum>0);
	pm2m[0] = a/sum; pm2d[0] = b/sum;
	for(i=1;i<=length-1;i++) {
		a = BitsToProb(m2m[i],1.0,bits);
		b = BitsToProb(m2i[i],1.0,bits);
		c = BitsToProb(m2d[i],1.0,bits);
		sum = a+b+c;
		pm2m[i]= a/sum; pm2i[i] = b/sum; pm2d[i] = c/sum;
		a = BitsToProb(i2i[i],1.0,bits);
		b = BitsToProb(i2m[i],1.0,bits);
		sum = a+b;
		pi2i[i] = a/sum; pi2m[i] = b/sum;
		a = BitsToProb(d2d[i],1.0,bits);
		b = BitsToProb(d2m[i],1.0,bits);
		sum = a+b;
		pd2d[i] = a/sum; pd2m[i] = b/sum;
	}
	
	pm2m[length] = 1.;
	pd2m[length] = 1.;

        NEW(gM2m,length+3,double);NEW(gM2i,length+3,double);
        NEW(gM2d,length+3,double);
        NEW(gI2m,length+3,double);NEW(gD2d,length+3,double);
        NEW(gD2m,length+3,double);

        NEW(RgM2m,length+3,double);NEW(RgM2i,length+3,double);
        NEW(RgM2d,length+3,double);
        NEW(RgI2m,length+3,double);NEW(RgD2d,length+3,double);
        NEW(RgD2m,length+3,double);

        NEW(rM2m,length+3,Int4);NEW(rM2i,length+3,Int4);
        NEW(rM2d,length+3,Int4);NEW(rI2i,length+3,Int4);
        NEW(rI2m,length+3,Int4);NEW(rD2d,length+3,Int4);
        NEW(rD2m,length+3,Int4);

	NEW(frM2m,length+3,Int4);NEW(frM2i,length+3,Int4);
	NEW(frM2d,length+3,Int4);NEW(frI2i,length+3,Int4);
	NEW(frI2m,length+3,Int4);NEW(frD2d,length+3,Int4);
	NEW(frD2m,length+3,Int4);

	gM2d[0] = pm2d[0]; gM2m[0] = 1. - gM2d[0];

        for(im1=0,i=1;i<=length-1;im1++,i++) {
                gM2m[i] = pm2m[i]*(gM2m[im1]+gI2m[im1]+gD2m[im1]); 
                gM2i[i] = pm2i[i]*(gM2m[im1]+gI2m[im1]+gD2m[im1]);
                gM2d[i] = pm2d[i]*(gM2m[im1]+gI2m[im1]+gD2m[im1]);
                gI2m[i] = gM2i[i];
                gD2d[i] = pd2d[i]*(gM2d[im1]+gD2d[im1]);
                gD2m[i] = pd2m[i]*(gM2d[im1]+gD2d[im1]);
        }

	gM2m[length] = pm2m[length]*(gM2m[length-1]+gI2m[length-1]+gD2m[length-1]);
	gD2m[length] = pd2m[length]*(gM2d[length-1]+gD2d[length-1]);

	for(i=0;i<=length;i++) {
		a = gM2m[length-i];
		b = gI2m[length-i];
        	c = gD2m[length-i];
		sum = a+b+c;
		RgM2m[i] = a/sum; RgM2i[i] = b/sum; RgM2d[i] = c/sum;

		RgI2m[i] = 1.-pi2i[length-i];

		a = gD2d[length-i];
		b = gM2d[length-i];
		sum = a+b;
		RgD2d[i] = a/sum; RgD2m[i] = b/sum;
	}

        rM2m[0] = (Int4) (floor) (0.5 + bits*log2(RgM2m[0])); 
        rM2d[0] = (Int4) (floor) (0.5 + bits*log2(RgM2d[0])); 

	for(i=1;i<=length-1;i++) {
		rM2m[i] = (Int4) (floor) (0.5 + bits*log2(RgM2m[i]));
		rM2i[i] = (Int4) (floor) (0.5 + bits*log2(RgM2i[i]));
		rM2d[i] = (Int4) (floor) (0.5 + bits*log2(RgM2d[i]));
		rI2i[i] = (Int4) (floor) (0.5 + bits*log2(pi2i[length-i]));
		rI2m[i] = (Int4) (floor) (0.5 + bits*log2(RgI2m[i]));
		rD2d[i] = (Int4) (floor) (0.5 + bits*log2(RgD2d[i]));
		rD2m[i] = (Int4) (floor) (0.5 + bits*log2(RgD2m[i]));
	}

	rM2m[length] = (Int4) (floor) (0.5 + bits*log2(RgM2m[length]));
	rD2m[length] = (Int4) (floor) (0.5 + bits*log2(RgD2m[length]));

//for gpxdrop only!!
	for(i=0;i<=length;i++) {
		frM2m[i] = rM2m[length-i];
		frM2d[i] = rM2d[length-i];
		frM2i[i] = rM2i[length-i];
		frD2m[i] = rD2m[length-i];
		frD2d[i] = rD2d[length-i];
		frI2i[i] = rI2i[length-i];
		frI2m[i] = rI2m[length-i];
		if(frM2d[i] == 0) frM2d[i] = -1;
		if(frM2i[i] == 0) frM2i[i] = -1;
		if(frD2d[i] == 0) frD2d[i] = -1;
		if(frI2i[i] == 0) frI2i[i] = -1;
	}

	NEW(frB2m,length+3,Int4);NEW(frM2e,length+3,Int4);
	for(i=0;i<=length;i++) { frB2m[i] = frM2e[i] = BigNegative; }
	Int4 frNB = BigNegative;Int4 frNN = BigNegative; Int4 frEC = BigNegative;
	Int4 frEJ = BigNegative;Int4 frCT = BigNegative; Int4 frCC = BigNegative;
	Int4 frJB = BigNegative;Int4 frJJ = BigNegative;
	Int4 frBM = BigNegative;Int4 frBD = BigNegative;
	Int4 frGG = BigNegative;Int4 frGF = BigNegative;
	reversed_stpb = new stpb_typ(length,frM2m,frM2i,frM2d,frI2i,frI2m,frD2d,frD2m,
		frB2m,frM2e,frNB,frNN,frEC,frEJ,frCT,frCC,frJB,frJJ,frBM,frBD,frGG,
			frGF,bits);

	free(pm2m);free(pm2i);free(pm2d);
	free(pi2i);free(pi2m);free(pd2d);
	free(pd2m);
        free(gM2m);free(gM2i);free(gM2d);
        free(gI2m);free(gD2d);
        free(gD2m);
        free(RgM2m);free(RgM2i);free(RgM2d);
        free(RgI2m);free(RgD2d);
        free(RgD2m);
	free(rM2m);free(rM2i);free(rM2d);
        free(rI2i);free(rI2m);free(rD2d);
        free(rD2m);
	return reversed_stpb;
}

stpb_typ *stpb_typ::ConcatenateStpb(Int4 rpts)
{
	stpb_typ *concat_stpb;
	Int4 i,j;
	Int4 *cM2m,*cM2i,*cM2d,*cI2i,*cI2m,*cD2d,*cD2m;
	NEW(cM2m,rpts*length+3,Int4);
	NEW(cM2i,rpts*length+3,Int4);
	NEW(cM2d,rpts*length+3,Int4);
	NEW(cI2i,rpts*length+3,Int4);
	NEW(cI2m,rpts*length+3,Int4);
	NEW(cD2d,rpts*length+3,Int4);
	NEW(cD2m,rpts*length+3,Int4);
	for(i=0;i<rpts;i++) {
		for(j=1;j<=length-1;j++) {
			cM2m[i*length+j] = m2m[j];
			cM2d[i*length+j] = m2d[j];
			cM2i[i*length+j] = m2i[j];
			cI2m[i*length+j] = i2m[j];
			cI2i[i*length+j] = i2i[j];
			cD2m[i*length+j] = d2m[j];
			cD2d[i*length+j] = d2d[j];
		}
	}

	double pjb = BitsToProb(jb,1.0,bits);
	double pjj = BitsToProb(jj,1.0,bits);
	double pbm = BitsToProb(bm,1.0,bits);
	double pbd = BitsToProb(bd,1.0,bits);
	double pMI = 1-pjb;

	for(i=1;i<=rpts-1;i++) {		
		cM2m[i*length] = ProbToBits(pjb*pbm,1.0,bits);
		cM2d[i*length] = ProbToBits(pjb*pbd,1.0,bits);
		cM2i[i*length] = ProbToBits(pMI,1.0,bits);
		cI2m[i*length] = jb;
		cI2i[i*length] = jj;
		cD2m[i*length] = bm;
		cD2d[i*length] = bd;
	}


	Int4 *cB2m,*cM2e;
	NEW(cB2m,rpts*length+3,Int4);NEW(cM2e,rpts*length+3,Int4);
	for(i=0;i<=rpts*length;i++) { cB2m[i] = cM2e[i] = BigNegative; }

	concat_stpb = new stpb_typ(rpts*length,cM2m,cM2i,cM2d,cI2i,cI2m,cD2d,cD2m,
		cB2m,cM2e,nb,nn,ec,ej,ct,cc,jb,jj,bm,bd,gg,
			gf,bits);

	return concat_stpb;
}

stpb_typ *stpb_typ::RenormalizeStpb()
//"wing retraction"
{
	stpb_typ *rn_stpb;
	
	Int4 i;
	Int4 *rn_m2m,*rn_m2i,*rn_m2d,*rn_i2i,*rn_i2m,*rn_d2d,*rn_d2m,*rn_b2m,*rn_m2e;

	NEW(rn_m2m,length+3,Int4);NEW(rn_m2i,length+3,Int4);NEW(rn_m2d,length+3,Int4);
	NEW(rn_i2i,length+3,Int4);NEW(rn_i2m,length+3,Int4);NEW(rn_d2d,length+3,Int4);
	NEW(rn_d2m,length+3,Int4);NEW(rn_b2m,length+3,Int4);NEW(rn_m2e,length+3,Int4);

	double *t_m2m,*t_m2i,*t_m2d,*t_i2i,*t_i2m,*t_d2d,*t_d2m,*t_b2m,*t_m2e;

	MEW(t_m2m,length+3,double);MEW(t_m2i,length+3,double);MEW(t_m2d,length+3,double);
	MEW(t_i2i,length+3,double);MEW(t_i2m,length+3,double);MEW(t_d2d,length+3,double);
	MEW(t_d2m,length+3,double);MEW(t_b2m,length+3,double);MEW(t_m2e,length+3,double);


	Int4 rn_nb,rn_nn,rn_ec,rn_ej,rn_ct,rn_cc,rn_jb,rn_jj;
	double t_nb,t_nn,t_ec,t_ej,t_ct,t_cc,t_jb,t_jj;


	for(i=1;i<=length;i++) {
		t_m2m[i] = BitsToProb(m2m[i],1.0,bits);
		t_m2i[i] = BitsToProb(m2i[i],1.0,bits);
		t_m2d[i] = BitsToProb(m2d[i],1.0,bits);
		t_i2m[i] = BitsToProb(i2m[i],1.0,bits);
		t_i2i[i] = BitsToProb(i2i[i],1.0,bits);
		t_d2m[i] = BitsToProb(d2m[i],1.0,bits);
		t_d2d[i] = BitsToProb(d2d[i],1.0,bits);
		t_b2m[i] = BitsToProb(b2m[i],1.0,bits);
		t_m2e[i] = BitsToProb(m2e[i],1.0,bits);
	}

	t_nb = BitsToProb(nb,1.0,bits);
	t_nn = BitsToProb(nn,1.0,bits);
	t_ec = BitsToProb(ec,1.0,bits);
	t_ej = BitsToProb(ej,1.0,bits);
	t_ct = BitsToProb(ct,1.0,bits);
	t_cc = BitsToProb(cc,1.0,bits);
	t_jb = BitsToProb(jb,1.0,bits);
	t_jj = BitsToProb(jj,1.0,bits);

	double accum,tbm,tme,p1=350./351.;
	Int4 k;
	double t_bd = BitsToProb(bd,1.0,bits);

//// Sean Eddy's routine ! /////////
  /* B->M entry transitions. Note how D_1 is folded out.
   * M1 is just B->M1
   * M2 is sum (or max) of B->M2 and B->D1->M2
   * M_k is sum (or max) of B->M_k and B->D1...D_k-1->M_k
   * These have to be done in log space, else you'll get
   * underflow errors; and we also have to watch for log(0).
   * A little sloppier than it probably has to be; historically,
   * doing in this in log space was in response to a bug report.
   */
  accum = t_bd > 0.0 ? log(t_bd) : -9999.;
  for (k = 1; k <= length; k++)
    {
      tbm = t_b2m[k] > 0. ? log(t_b2m[k]) : -9999.;   /* B->M_k part */
     
      /* B->D1...D_k-1->M_k part we get from accum*/
      if (k > 1 && accum > -9999.)
        {
          if (t_d2m[k-1] > 0.0)
            {
              tbm =  MAXIMUM(double,tbm, accum + log(t_d2m[k-1]));
            }
             
          accum = (t_d2d[k-1] > 0.0) ? accum + log(t_d2d[k-1]) : -9999.;
        }
                                /* Convert from log_e to scaled integer log_2 odds. */
      if (tbm > -9999.)
        rn_b2m[k] = (int) floor(0.5 + bits * 1.44269504 * (tbm - log(p1)));
      else
        rn_b2m[k] = -BigNegative;
//printf("hmm->bsc[%d]=%d\n",k,hmm->bsc[k]);
    }


  /* M->E exit transitions. Note how D_M is folded out.
   * M_M is 1 by definition
   * M_M-1 is sum of M_M-1->E and M_M-1->D_M->E, where D_M->E is 1 by definition
   * M_k is sum of M_k->E and M_k->D_k+1...D_M->E
   * Must be done in log space to avoid underflow errors.
   * A little sloppier than it probably has to be; historically,
   * doing in this in log space was in response to a bug report.
   */
  rn_m2e[length] = 0;
  accum = 0.;
  for (k = length-1; k >= 1; k--)
    {
      tme = t_m2e[k] > 0. ? log(t_m2e[k]) : -9999.;
      if (accum > -9999.)
        {
          if (t_m2d[k] > 0.0)
            {
              tme = MAXIMUM(double, tme, accum + log(t_m2d[k]));
            }
          accum = (t_d2d[k] > 0.0) ? accum + log(t_d2d[k]) : -9999.;
        } 
                                /* convert from log_e to scaled integer log odds. */
      rn_m2e[k] = (tme > -9999.) ? (int) floor(0.5 + bits * 1.44269504 * tme) : -BigNegative;
    }
                                
                                /* special transitions */
	rn_nn = ProbToBits(t_nn,p1,bits);
	rn_nb = ProbToBits(t_nb,1.0,bits);
	rn_ej = ProbToBits(t_ej,1.0,bits);
	rn_ec = ProbToBits(t_ec,1.0,bits);
	rn_cc = ProbToBits(t_cc,p1,bits);
	rn_ct = ProbToBits(t_ct,1.-p1,bits);
	rn_jj = ProbToBits(t_jj,p1,bits);
	rn_jb = ProbToBits(t_jb,1.0,bits);

	for(i=1;i<=length;i++) {
		rn_m2m[i] = ProbToBits(t_m2m[i],p1,bits);
		rn_m2d[i] = ProbToBits(t_m2d[i],1.0,bits);
		rn_m2i[i] = ProbToBits(t_m2i[i],p1,bits);
		rn_i2m[i] = ProbToBits(t_i2m[i],p1,bits);
		rn_i2i[i] = ProbToBits(t_i2i[i],p1,bits);
		rn_d2m[i] = ProbToBits(t_d2m[i],p1,bits);
		rn_d2d[i] = ProbToBits(t_d2d[i],1.0,bits);
	}


	rn_stpb = new stpb_typ(length,rn_m2m,rn_m2i,rn_m2d,rn_i2i,rn_i2m,rn_d2d,rn_d2m,rn_b2m,rn_m2e,
		rn_nb,rn_nn,rn_ec,rn_ej,rn_ct,rn_cc,rn_jb,rn_jj,bm,bd,gg,gf,bits);

	free(t_m2m); free(t_m2d); free(t_m2i); free(t_i2m); free(t_i2i);
	free(t_d2m); free(t_d2d); free(t_b2m); free(t_m2e);

	return rn_stpb;
}
#if 0
stpb_typ *stpb_typ::TempStpb(double temp, Int4 s)
// s is a scaling factor, i.e.
// to bring probability p to converge to q we use formula:
// (1+s)/(t+s)*p + (1 - (1+s)/(t+s))*q
{
        stpb_typ        *temp_stpb;

	Int4 *temp_m2m,*temp_m2i,*temp_m2d,*temp_i2i,*temp_i2m,*temp_d2d,*temp_d2m,*temp_b2m,*temp_m2e;
	Int4 temp_nb,temp_nn,temp_ec,temp_ej,temp_ct,temp_cc,temp_jb,temp_jj;

	NEW(temp_m2m,length+3,Int4);NEW(temp_m2i,length+3,Int4);NEW(temp_m2d,length+3,Int4);
	NEW(temp_i2i,length+3,Int4);NEW(temp_i2m,length+3,Int4);NEW(temp_d2d,length+3,Int4);
	NEW(temp_d2m,length+3,Int4);NEW(temp_b2m,length+3,Int4);NEW(temp_m2e,length+3,Int4);

        Int4            i;
        double          *m2mp,*m2ip,*m2dp,*i2mp,*i2ip,*d2mp,*d2dp,*b2mp,*m2ep;
        double          nbp,nnp,ecp,ejp,ctp,ccp,jbp,jjp;
        double          sum,p1 = 350./351;
        NEW(m2mp,length+2,double);NEW(m2ip,length+2,double);NEW(m2dp,length+2,double);
        NEW(i2mp,length+2,double);NEW(i2ip,length+2,double);NEW(d2mp,length+2,double);
        NEW(d2dp,length+2,double);NEW(b2mp,length+2,double);NEW(m2ep,length+2,double);

	double backm = 1./4;
	double backd = 1./2;
	double backi = 1./2;
	double backb = 1./length;
        for(i=1;i<=length-1;i++) {
		if (i == length - 1) backm = 1./3;  
                m2mp[i] = (1.+s)/(temp+s)*BitsToProb(m2m[i],p1,bits) + (1.-(1.+s)/(temp+s))*backm;
                m2ip[i] = (1.+s)/(temp+s)*BitsToProb(m2i[i],p1,bits) + (1.-(1.+s)/(temp+s))*backm;
                m2dp[i] = (1.+s)/(temp+s)*BitsToProb(m2d[i],1.0,bits) + (1.-(1.+s)/(temp+s))*backm;
                i2mp[i] = (1.+s)/(temp+s)*BitsToProb(i2m[i],p1,bits) + (1.-(1.+s)/(temp+s))*backi;
                i2ip[i] = (1.+s)/(temp+s)*BitsToProb(i2i[i],p1,bits) + (1.-(1.+s)/(temp+s))*backi;
		if (i == 1) { d2mp[i] = d2dp[i] = 0; }
		else {
                	d2mp[i] = (1.+s)/(temp+s)*BitsToProb(d2m[i],p1,bits) + (1.-(1.+s)/(temp+s))*backd;
                	d2dp[i] = (1.+s)/(temp+s)*BitsToProb(d2d[i],1.0,bits) + (1.-(1.+s)/(temp+s))*backd;
		}
                b2mp[i] = (1.+s)/(temp+s)*BitsToProb(b2m[i],p1,bits) + (1.-(1.+s)/(temp+s))*backb;
                m2ep[i] = (1.+s)/(temp+s)*BitsToProb(m2e[i],1.0,bits) + (1.-(1.+s)/(temp+s))*backm;
        }
	d2dp[length-1] = 0.; 
	b2mp[length] = (1.+s)/(temp+s)*BitsToProb(b2m[length],p1,bits) + (1.-(1.+s)/(temp+s))*backb;
	m2ep[length] = 1.;
	d2mp[length] = 1.;

	for(i=1;i<=length-1;i++) {
		sum = m2mp[i] + m2ip[i] + m2dp[i] + m2ep[i];
		m2mp[i] /= sum; m2ip[i] /= sum; m2dp[i] /= sum; m2ep[i] /= sum;
		temp_m2m[i] = ProbToBits(m2mp[i],p1,bits);
		temp_m2i[i] = ProbToBits(m2ip[i],p1,bits);
		temp_m2d[i] = ProbToBits(m2dp[i],1.0,bits);
		temp_m2e[i] = ProbToBits(m2ep[i],1.0,bits);
		sum = i2mp[i] + i2ip[i];
		i2mp[i] /= sum; i2ip[i] /= sum;
		temp_i2m[i] = ProbToBits(i2mp[i],p1,bits);
		temp_i2i[i] = ProbToBits(i2ip[i],p1,bits);
	}
	temp_m2e[length] = ProbToBits(1.0,1.0,bits);

	for(i=2;i<=length-2;i++) {
		sum = d2mp[i] + d2dp[i];
		d2mp[i] /=sum; d2dp[i] /= sum;
		temp_d2m[i] = ProbToBits(d2mp[i],p1,bits); 
		temp_d2d[i] = ProbToBits(d2dp[i],1.0,bits);
	}
	temp_d2m[length-1] = ProbToBits(1.0,1.0,bits);;

	for(sum = 0.,i=1;i<=length;i++) {
		sum += b2mp[i];
	}

	for(i=1;i<=length;i++) {
		b2mp[i] /= sum;
		temp_b2m[i] = ProbToBits(b2mp[i],p1,bits);
	}

	double backs = 1./2;

        nbp = (1.+s)/(temp+s)*BitsToProb(nb,1.0,bits) + (1.-(1.+s)/(temp+s))*backs;
        nnp = (1.+s)/(temp+s)*BitsToProb(nn,p1,bits) + (1.-(1.+s)/(temp+s))*backs;
	sum = nbp + nnp;
	nbp /= sum; nnp /= sum;
        ecp = (1.+s)/(temp+s)*BitsToProb(ec,1.0,bits) + (1.-(1.+s)/(temp+s))*backs;
        ejp = (1.+s)/(temp+s)*BitsToProb(ej,1.0,bits) + (1.-(1.+s)/(temp+s))*backs;
	sum = ecp + ejp;
        ecp /= sum; ejp /= sum;
        ctp = (1.+s)/(temp+s)*BitsToProb(ct,1.0-p1,bits) + (1.-(1.+s)/(temp+s))*backs;
        ccp = (1.+s)/(temp+s)*BitsToProb(cc,p1,bits) + (1.-(1.+s)/(temp+s))*backs;
	sum = ctp + ccp;
        ctp /= sum; ccp /= sum;
        jbp = (1.+s)/(temp+s)*BitsToProb(jb,1.0,bits) + (1.-(1.+s)/(temp+s))*backs;
        jjp = (1.+s)/(temp+s)*BitsToProb(jj,p1,bits) + (1.-(1.+s)/(temp+s))*backs;
        sum = jbp + jjp;
        jbp /= sum; jjp /= sum;

	temp_nb = ProbToBits(nbp,1.0,bits);
	temp_nn = ProbToBits(nnp,p1,bits);
	temp_ec = ProbToBits(ecp,1.0,bits);
	temp_ej = ProbToBits(ejp,1.0,bits);
	temp_ct = ProbToBits(ctp,1.0-p1,bits);
	temp_cc = ProbToBits(ccp,p1,bits);
	temp_jb = ProbToBits(jbp,1.0,bits);
	temp_jj = ProbToBits(jjp,p1,bits);

	temp_stpb = new stpb_typ(length,temp_m2m,temp_m2i,temp_m2d,temp_i2i,temp_i2m,temp_d2d,temp_d2m,
			temp_b2m,temp_m2e,temp_nb,temp_nn,temp_ec,temp_ej,temp_ct,temp_cc,temp_jb,
			temp_jj,bm,bd,gg,gf,bits);

	free(m2mp); free(m2dp); free(m2ip); free(i2mp); free(i2ip);
	free(d2mp); free(d2dp); free(b2mp); free(m2ep);

	return temp_stpb;
}
#endif

stpb_typ *stpb_typ::TempStpb(double temp)
//temperature stpb
{
        stpb_typ        *temp_stpb;

	Int4 *temp_m2m,*temp_m2i,*temp_m2d,*temp_i2i,*temp_i2m,*temp_d2d,*temp_d2m,*temp_b2m,*temp_m2e;
	Int4 temp_nb,temp_nn,temp_ec,temp_ej,temp_ct,temp_cc,temp_jb,temp_jj;

	NEW(temp_m2m,length+3,Int4);NEW(temp_m2i,length+3,Int4);NEW(temp_m2d,length+3,Int4);
	NEW(temp_i2i,length+3,Int4);NEW(temp_i2m,length+3,Int4);NEW(temp_d2d,length+3,Int4);
	NEW(temp_d2m,length+3,Int4);NEW(temp_b2m,length+3,Int4);NEW(temp_m2e,length+3,Int4);

        Int4            i;
        double          *m2mp,*m2ip,*m2dp,*i2mp,*i2ip,*d2mp,*d2dp,*b2mp,*m2ep;
        double          nbp,nnp,ecp,ejp,ctp,ccp,jbp,jjp;
        double          sum,p1 = 350./351;
        NEW(m2mp,length+2,double);NEW(m2ip,length+2,double);NEW(m2dp,length+2,double);
        NEW(i2mp,length+2,double);NEW(i2ip,length+2,double);NEW(d2mp,length+2,double);
        NEW(d2dp,length+2,double);NEW(b2mp,length+2,double);NEW(m2ep,length+2,double);

        for(i=1;i<=length-1;i++) {
                m2mp[i] = pow(BitsToProb(m2m[i],p1,bits),1./temp);
                m2ip[i] = pow(BitsToProb(m2i[i],p1,bits),1./temp);
                m2dp[i] = pow(BitsToProb(m2d[i],1.0,bits),1./temp);
                i2mp[i] = pow(BitsToProb(i2m[i],p1,bits),1./temp);
                i2ip[i] = pow(BitsToProb(i2i[i],p1,bits),1./temp);
		if (i == 1) { d2mp[i] = d2dp[i] = 0; }
		else {
                	d2mp[i] = pow(BitsToProb(d2m[i],p1,bits),1./temp);
                	d2dp[i] = pow(BitsToProb(d2d[i],1.0,bits),1./temp);
		}
                b2mp[i] = pow(BitsToProb(b2m[i],p1,bits),1./temp);
                m2ep[i] = pow(BitsToProb(m2e[i],1.0,bits),1./temp);
        }
	d2dp[length-1] = 0.; 
	b2mp[length] = pow(BitsToProb(b2m[length],p1,bits),1./temp);
	m2ep[length] = 1.;
	d2mp[length] = 1.;

	for(i=1;i<=length-1;i++) {
		sum = m2mp[i] + m2ip[i] + m2dp[i] + m2ep[i];
		m2mp[i] /= sum; m2ip[i] /= sum; m2dp[i] /= sum; m2ep[i] /= sum;
		temp_m2m[i] = ProbToBits(m2mp[i],p1,bits);
		temp_m2i[i] = ProbToBits(m2ip[i],p1,bits);
		temp_m2d[i] = ProbToBits(m2dp[i],1.0,bits);
		temp_m2e[i] = ProbToBits(m2ep[i],1.0,bits);
		sum = i2mp[i] + i2ip[i];
		i2mp[i] /= sum; i2ip[i] /= sum;
		temp_i2m[i] = ProbToBits(i2mp[i],p1,bits);
		temp_i2i[i] = ProbToBits(i2ip[i],p1,bits);
	}
	temp_m2e[length] = ProbToBits(1.0,1.0,bits);

	for(i=2;i<=length-2;i++) {
		sum = d2mp[i] + d2dp[i];
		d2mp[i] /=sum; d2dp[i] /= sum;
		temp_d2m[i] = ProbToBits(d2mp[i],p1,bits); 
		temp_d2d[i] = ProbToBits(d2dp[i],1.0,bits);
	}
	temp_d2m[length-1] = ProbToBits(1.0,1.0,bits);;

	for(sum = 0.,i=1;i<=length;i++) {
		sum += b2mp[i];
	}

	for(i=1;i<=length;i++) {
		b2mp[i] /= sum;
		temp_b2m[i] = ProbToBits(b2mp[i],p1,bits);
	}
#if 0
        nbp = pow(BitsToProb(nb,1.0,bits),1./temp);
        nnp = pow(BitsToProb(nn,p1,bits),1./temp);
	sum = nbp + nnp;
	nbp /= sum; nnp /= sum;
        ecp = pow(BitsToProb(ec,1.0,bits),1./temp);
        ejp = pow(BitsToProb(ej,1.0,bits),1./temp);
	sum = ecp + ejp;
        ecp /= sum; ejp /= sum;
        ctp = pow(BitsToProb(ct,1.0-p1,bits),1./temp);
        ccp = pow(BitsToProb(cc,p1,bits),1./temp);
	sum = ctp + ccp;
        ctp /= sum; ccp /= sum;
        jbp = pow(BitsToProb(jb,1.0,bits),1./temp);
        jjp = pow(BitsToProb(jj,p1,bits),1./temp);
        sum = jbp + jjp;
        jbp /= sum; jjp /= sum;

	temp_nb = ProbToBits(nbp,1.0,bits);
	temp_nn = ProbToBits(nnp,p1,bits);
	temp_ec = ProbToBits(ecp,1.0,bits);
	temp_ej = ProbToBits(ejp,1.0,bits);
	temp_ct = ProbToBits(ctp,1.0-p1,bits);
	temp_cc = ProbToBits(ccp,p1,bits);
	temp_jb = ProbToBits(jbp,1.0,bits);
	temp_jj = ProbToBits(jjp,p1,bits);

	temp_stpb = new stpb_typ(length,temp_m2m,temp_m2i,temp_m2d,temp_i2i,temp_i2m,temp_d2d,temp_d2m,
			temp_b2m,temp_m2e,temp_nb,temp_nn,temp_ec,temp_ej,temp_ct,temp_cc,temp_jb,
			temp_jj,bm,bd,gg,gf,bits);
#endif

	temp_stpb = new stpb_typ(length,temp_m2m,temp_m2i,temp_m2d,temp_i2i,temp_i2m,temp_d2d,temp_d2m,
			temp_b2m,temp_m2e,nb,nn,ec,ej,ct,cc,jb,jj,bm,bd,gg,gf,bits);

	free(m2mp); free(m2dp); free(m2ip); free(i2mp); free(i2ip);
	free(d2mp); free(d2dp); free(b2mp); free(m2ep);

	return temp_stpb;
}

