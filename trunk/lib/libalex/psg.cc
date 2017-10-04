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

#include "psg.h"

psg_typ::psg_typ() { assert(!"psg_typ() illegal constructor"); }

psg_typ::psg_typ(cma_typ cma, Int4 Pernats)
{
	Mkpsg_typ(cma,'h','h',0,Pernats);
}

psg_typ::psg_typ(cma_typ cma, char w, char p, double *m, Int4 Pernats)
// m = position specific constant (total pseudocounts);
{ Mkpsg_typ(cma, w, p, m, Pernats); }

void psg_typ::Mkpsg_typ(cma_typ cma, char w, char pc, double *m, Int4 Pernats)
{
	Int4		i, l, t, n, posit, st, p, col;
	Int4 		sum = 0, pos[4];
	unsigned char	res;

	Int4 len=TotalLenCMSA(cma);
	NEW(h_mult,len+2,double);
	if(m == 0) for(i=1;i <= len+1;i++) h_mult[i] = 5.;
	else for(i=1;i <= len+1;i++) h_mult[i] = m[i];

	weightMeth = w, psdcntMeth = pc; pernats = Pernats;

	A = AlphabetCMSA(cma);
	SeqSet = DataCMSA(cma);
	nSeq = NumSeqsCMSA(cma);
	nBlks = nBlksCMSA(cma);

	NEWP(blPos,nBlks+1,Int4);
	for(i=0;i<=nBlks;i++) NEW(blPos[i],nSeq+1,Int4);	

	NEW(blLen,nBlks+1,Int4);
	for(i=1;i<=nBlks;i++) {
		blLen[i] = LengthCMSA(i,cma);
		sum += blLen[i];
	}
	nCol = sum;

	NEWP(rawCnts,nCol+1,Int4);
	for(t=0;t<=nCol;t++){ NEW(rawCnts[t],nAlpha(A)+1,Int4); }
	
	for(n=1;n<=nSeq;n++){
		for(t=1;t<=nBlks;t++){
			st = PosSiteCMSA(t,n,pos,cma);
			if(st != 0) blPos[t][n] = pos[1];
			else blPos[t][n] = 0;
		}
	}
	for(n=1;n<=nSeq;n++){
		col = 0;
		for(t=1;t<=nBlks;t++){
			p = blPos[t][n];
			if(p != 0){
				for(l=1;l<=blLen[t];l++){
					col += 1;	
					posit = p+l-1;
					res=SeqP(n,posit,SeqSet);
					rawCnts[col][res] += 1;
				}
			}
		}
	}
	
	NEW(nResTyp,nCol+1,Int4);
	for(t=1;t<=nCol;t++){
		for(i=1;i<=nAlpha(A);i++){
			if (rawCnts[t][i] != 0)
				nResTyp[t] += 1;		
		}		
	}
	weightedCnts = NULL;
	psdCnts = NULL;

}
psg_typ::psg_typ(cma_typ cma, char w, char p, double *m, Int4 Pernats, char meth)
{
	if (meth == 'r'){ //without weighting
		Mkpsg_typ(cma, w, p, m, Pernats);		
	}
	else if (meth == 'w'){ //with weighting
		Mkpsg_typ(cma, w, p, m, Pernats);
		weightedCnts = WeightedCounts();
		psdCnts = PseudoCounts();
	}
	else print_error("psg_typ(): Unknown method");
}

psg_typ::~psg_typ() { Free(); }

void psg_typ::Free()
{
	Int4 t;

	if(rawCnts) {
		for(t=0;t<=nCol;t++){ free(rawCnts[t]); } 
		free(rawCnts); 
	}
	if(nResTyp) free(nResTyp);
	if(blPos) {
		for(t=0;t<=nBlks;t++){ free(blPos[t]); } 
		free(blPos);
	}
	if(h_mult) free(h_mult);
	if(weightedCnts) {
		for(t=0;t<=nCol;t++) free(weightedCnts[t]); 
		free(weightedCnts); 
	}
	if(psdCnts) {
		for(t=0;t<=nCol;t++) free(psdCnts[t]); 
		free(psdCnts);
	}
	if(blLen) free(blLen);

}

smx_typ *psg_typ::SampleBlcksMatrices(Int4 nRpts, Int4 wt)
{
	static Int4 Seed = 0;

        if (Seed == 0) Seed = Random();

	double 		**counts;
	if(weightedCnts) counts = weightedCnts;
	else print_error("SampleBlcksMatrices(): have to use weighting");
        smx_typ *M;
        Int4            i,k,a,b,r,col,o_l,totalws_l;
        double          o,*observed,*cnts,*pseudo,obs,ave,totalws,totalw;
	double 		t1,t2,sll,slh,shl,shh,s,*totalW, *totalWS, *totalPS;
	double 		random_num1, random_num2;
	NEW(totalWS,nCol+2,double);
	NEW(totalPS,nCol+2,double);
	NEW(totalW,nCol+2,double);
	for(i=1;i<=nCol;i++) {
		cnts = counts[i];
		pseudo = psdCnts[i];
		for(a=1;a<=nAlpha(A);a++){
			totalWS[i] += cnts[a];
			totalPS[i] += pseudo[a];
		}
		totalW[i] = totalWS[i] + totalPS[i];
	}

 	if(wt <= 0 || wt > 24){
		fprintf(stderr,"wt = %d\n",wt);
		assert(wt > 0 && wt < 20);        
	}

        NEW(M,nRpts*nBlks+1,smx_typ);
        for(b=1;b<=nBlks;b++){
            M[b] = MkSMatrix(2.0,blLen[b],blosum62freq,A);
            for(r=1; r < nRpts;r++){ M[r*nBlks+b] = M[b]; }
        }
        for(col=b=1;b<=nBlks;b++){
           for(i=1;i<=blLen[b];i++){
	     observed = counts[col];
	     pseudo = psdCnts[col];
	     totalws = totalWS[col];
	     totalws_l = floor(totalws);
	     totalw = totalW[col];
             for(ave=0,a=1; a<=nAlpha(A); a++){
		o = observed[a];
		o_l=floor(o);
		random_num1 = SampleUniformProb();
		random_num2 = SampleUniformProb();
		if(random_num1 > o - o_l) o_l++;
		if(random_num2 > totalws - totalws_l) totalws_l++;
        	if(o_l > 0 && totalws_l > 0) {
                	obs = (double)totalws_l*betadev(o_l*wt,totalws_l*wt,&Seed);
                	if (obs > 0.0) {
                        	s = ((obs+pseudo[a])/(totalw+(obs-o_l)))/blosum62freq[a];
                	} else s = ((o_l+pseudo[a])/totalw)/blosum62freq[a];
        	} else if(o_l) s = 1.0/blosum62freq[a]; else s = (pseudo[a]/totalw)/blosum62freq[a];
		SetSMatrix(a,i,(Int4) floor((pernats*log(s)+0.5)),M[b]);
                ave += s*blosum62freq[a]; // expected background score.
             } col++;
             s = (Int4)floor(ave);  
             SetSMatrix(0,i,s,M[b]);
           }
        }
	free(totalW); free(totalWS); free(totalPS);
        return M;
}


double *psg_typ::HenikoffWeights()
{
	Int4 		n, p, posit, t, l, col=0;
	double 		maxw=0., *w, sum;
	unsigned char 	res;

	NEW(w,nSeq+1,double);
	
	for(t=1;t<=nBlks;t++){
		for(l=1;l<=blLen[t];l++){
			col += 1;
			for(n=1;n<=nSeq;n++){
				p = blPos[t][n];
				if(p !=0){
					posit = p+l-1;
					res = SeqP(n,posit,SeqSet);
					w[n] += 1./(rawCnts[col][res]*nResTyp[col]);
				}
			}
        	}
	}

	for(n=1;n<=nSeq;n++){
		if(w[n] > maxw) {
			maxw = w[n];
		}
	}

	for(n=1;n<=nSeq;n++) { w[n] /= maxw; }

	return w;
}

double **psg_typ::WeightedCounts()
{
	if (weightMeth == 'h') return HenikoffWeightedCounts();
		else print_error("unknown weighting method");
}

double **psg_typ::HenikoffWeightedCounts()
{
	double 		**wc;
	Int4		t, l, n, p, posit, col;
	unsigned char	res;

	NEWP(wc,nCol+1,double);
	for(t=0;t<=nCol;t++)
		NEW(wc[t],nAlpha(A)+1,double);

	double *w = HenikoffWeights();

        for(n=1;n<=nSeq;n++){
                col = 0; 
                for(t=1;t<=nBlks;t++){
                        p = blPos[t][n];
			if(p != 0){
                        	for(l=1;l<=blLen[t];l++){
                                	col += 1;
                                	posit = p+l-1;
                                	res = SeqP(n,posit,SeqSet);
                                	wc[col][res] += w[n];
                        	}
			}
                }
        }

	free(w);
	return wc;
}

double **psg_typ::PseudoCounts()
{
	if (psdcntMeth == 'h') return HenikoffPseudoCounts();
	if (psdcntMeth == 'o') return OnePseudoCount();
	else print_error("unknown pseudocounts method");
}

double **psg_typ::OnePseudoCount()
{
	Int4 a,t;
	double **b;
	NEWP(b,nCol+1,double);
	for(t=0;t<=nCol;t++)
		NEW(b[t],nAlpha(A)+1,double);
	for(t=1;t<=nCol;t++){
		for(a=1;a<=nAlpha(A);a++)
			b[t][a]=1.;
	}
	return b;
}

double **psg_typ::HenikoffPseudoCounts()
{
	double	**b;
	Int4 	a, r, t, i;
	double 	*N, *Q, sum, B;
	double **wc = WeightedCounts();

	NEWP(b,nCol+1,double);
	for(t=0;t<=nCol;t++) NEW(b[t],nAlpha(A)+1,double);	
	NEW(N,nCol+1,double); NEW(Q,nCol+1,double);
	for(t=1;t<=nCol;t++){
		for(a=1;a<=nAlpha(A);a++){ N[t] += wc[t][a]; }
	}	
	for(a=1;a<=nAlpha(A);a++){
		for(r=1;r<=nAlpha(A);r++){ Q[a] += blosum62P[a][r]; }
	}
	for(t=1;t<=nCol;t++){
		B=h_mult[t]*nResTyp[t];
		for(a=1;a<=nAlpha(A);a++){
			for(sum=0,i=1;i<=nAlpha(A);i++){
				sum+=(wc[t][i]*blosum62P[i][a])/(N[t]*Q[i]);
			} b[t][a]=B*sum;
		}
	} free(N); free(Q);

	for(t=0;t<=nCol;t++) free(wc[t]); free(wc);
	return b;
}
double **psg_typ::TargetFreq()
{
	Int4 	t, a;
	double 	**p;
	double 	N, B;

	NEWP(p,nCol+1,double);
        for(t=1;t<=nCol;t++) NEW(p[t],nAlpha(A)+1,double);

	double **wc = WeightedCounts();
	double **b = PseudoCounts(); 

	for(t=1;t<=nCol;t++){
		N = 0.;
		for(a=1;a<=nAlpha(A);a++) N += wc[t][a];
		if (psdcntMeth == 'h') B = h_mult[t]*nResTyp[t];
		else if (psdcntMeth == 'o') B = 20.;
		else print_error("unknown pseudocounts method");
		for(a=1;a<=nAlpha(A);a++)
			p[t][a] = (wc[t][a]+b[t][a])/(N+B);	
	}
	for(t=0;t<=nCol;t++) {free(wc[t]); free(b[t]);}
	free(wc); free(b);
	return p;
}

Int4 **psg_typ::LogOddMatrix()
{
	Int4 **mtrx,i,a,t;
	double **p = TargetFreq();
	NEWP(mtrx,nCol+1,Int4);
	for(i=1;i<=nCol;i++)
		NEW(mtrx[i],nAlpha(A)+1,Int4);
	for(t=1;t<=nCol;t++){
		for(a=1;a<=nAlpha(A);a++)
			mtrx[t][a] = pernats*log(p[t][a]/blosum62freq[a]); 
	}
        for(t=1;t<=nCol;t++) free(p[t]); free(p);
	return mtrx;
}

Int4 **psg_typ::LogOddMatrix(Int4 PerBits)
	// return log-odds matrix in PerBits units.
{
	Int4 **mtrx,i,a,t;
	double **p = TargetFreq();
	NEWP(mtrx,nCol+1,Int4);
	for(i=1;i<=nCol;i++)
		NEW(mtrx[i],nAlpha(A)+1,Int4);
	for(t=1;t<=nCol;t++){
		for(a=1;a<=nAlpha(A);a++)
			mtrx[t][a] = PerBits*1.44269504*log(p[t][a]/0.05); 
	}
        for(t=1;t<=nCol;t++) free(p[t]); free(p);
	return mtrx;
}

smx_typ *psg_typ::BlcksMatrices(Int4 nRpts)
{
	smx_typ *M;
        Int4            r,i,t,k,a,b,s,score;
	double 		tmp,ave;

	NEW(M,nRpts*nBlks+1,smx_typ);
        for(b=1;b<=nBlks;b++){
	    M[b] = MkSMatrix(2.0,blLen[b],blosum62freq,A);
	    for(r=1; r < nRpts;r++){ M[r*nBlks+b] = M[b]; }
	}
	double **p = TargetFreq();
	for(s=b=1;b<=nBlks;b++){
	   for(i=1;i<=blLen[b];i++){
	     for(ave=0,a=1; a<=nAlpha(A); a++){
		tmp = pernats*log(p[s][a]/blosum62freq[a]);
		if(tmp >= 0.) score = (Int4)(tmp+0.5);
		else score = (Int4)(tmp-0.5);
		ave += tmp*blosum62freq[a]; // expected background score.
		SetSMatrix(a,i,score,M[b]);
	     } s++;
	     score = (Int4)floor(ave);
	     SetSMatrix(0,i,score,M[b]);
	   }
	}
	for(t=0;t<=nCol;t++) free(p[t]); free(p); 
	return M;
}

double *psg_typ::RelEntropy(Int4 rpts)
{
        double  *info;
        Int4    n,t, a;
        double  sum;

        NEW(info,rpts*nCol+1,double);
        double **q = TargetFreq();

        for(t=1;t<=nCol;t++){
                for(a=1;a<=nAlpha(A);a++){
                        info[t] += q[t][a]*log(q[t][a]/blosum62freq[a]);
                }
        }
	for(n=0;n<rpts;n++){
		for(t=1;t<=nCol;t++){
			info[n*nCol+t] = info[t];
		}
	}
        for(t=0;t<=nCol;t++) free(q[t]);
        free(q);
        return info;
}

Int4 *psg_typ::RelEnPosit(double cutoff, Int4 rpts, Int4 *length)
{
        Int4    i, k=1;
        Int4    *pos;

        NEW(pos,rpts*nCol+1,Int4);
        double  *info = RelEntropy(rpts);

        for(i=1;i<=rpts*nCol;i++){
                if(info[i] > cutoff){
                        pos[k++] = i;
                }
        }
        *length = k-1;

        free(info);
        return pos;
}
