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

#include "marg_prob.h"

#if 0	//*** Moved to my_ncbi.h
sap_typ TrimSAP(sap_typ sap, Int4 trimleft, Int4 trimright)
//trims trimleft residues in subj seq from sap at the begining and
//trimright at the end
{ 
	e_type qE = QuerySeqGSAP(sap);
	e_type sE = SubjectSeqGSAP(sap);
	dsp_typ dsp = sap->segs;
	Int4 numseg = dsp->numseg;
	Int4 *starts= dsp->starts;
	Int4 *lens = dsp->lens;
	Int4 len_s = starts[2*numseg-1]+lens[numseg-1]-1;
	assert(trimleft >=0 && trimright >= 0 && trimleft + trimright <= len_s);
	if(trimleft + trimright == len_s) return NULL;
	Int4 newNumseg, block_left, block_right, offset_left=0, offset_right=0, i, j;
	Int4 x,y;
	if(trimleft == 0) {
		block_left = 0; offset_left = 0;
	} else {	
		for(i=0;i<numseg;i++) {
			if ((x=starts[2*i+1]+lens[i]) >= (y=starts[1]+trimleft)) {
				break;
			}
		}
		if((x==y) && (i != numseg-1) && starts[2*i] != -1 && starts[2*i+1] != -1) { 
			i += 2; block_left = i; offset_left = 0; 
		}
		else if (starts[2*i] == -1 || starts[2*i+1] == -1) {
			i++; block_left = i; offset_left = 0;
		} else {
			block_left = i;
			for(j=1; j<=lens[i];j++) {
				if(starts[2*i+1]+j >= starts[1]+trimleft) {
					offset_left = j; break;
				}
			}
		}
	}
	if(trimright == 0) {
		block_right = numseg-1; offset_right = 0;
	} else {
		for(i=0;i<=numseg-1;i++) {
			if ((x=starts[2*i+1]+lens[i]) >= (y=starts[2*numseg-1]+lens[numseg-1]-trimright)) {
				break;
			}
		}
		if(starts[2*i] == -1 || starts[2*i+1] == -1){
			i--; block_right = i;  offset_right = 0;
		} else {
			block_right = i; 
			for(j=0;j<lens[i];j++) {
				if(starts[2*i+1] + lens[i] - j <= starts[2*numseg-1] + lens[numseg-1] - trimright) {
					offset_right = j; break;
				}
			}
		}
	}
	newNumseg = block_right - block_left + 1;
	sap_typ newSAP = GSeqAlignNew();
        dsp_typ newDSP = GDenseSegNew();
        newDSP->dim = 2; newDSP->numseg = newNumseg;
        newDSP->subject_id = SeqI(sE); newDSP->sE=sE;
        newDSP->query_id=0; newDSP->qE=qE;
	Int4 *newLens,*newStarts;
	NEW(newLens,newNumseg,Int4);
	NEW(newStarts,2*newNumseg,Int4);
	if(block_left == block_right) { 
		newLens[0] = lens[block_left] - offset_left - offset_right;
	}
	else { 
		newLens[0] = lens[block_left] - offset_left; 
		newLens[newNumseg-1] = lens[block_right] - offset_right; 
		for(i=1;i<newNumseg-1;i++) { 
			newLens[i] = lens[block_left+i]; 
		}
	}
	newStarts[0] = starts[2*block_left] + offset_left;
	newStarts[1] = starts[2*block_left+1] + offset_left;
	for(i=1;i<=newNumseg-1;i++) { 
		newStarts[2*i] = starts[2*(block_left+i)]; 
		newStarts[2*i+1] = starts[2*(block_left+i)+1];
	}
        newDSP->starts = newStarts; newDSP->lens = newLens;
        newSAP->segs = newDSP; newSAP->next = NULL;
        newSAP->dim =2;
	return newSAP;
}
#endif

sap_typ	MargProbTrimSAP(sap_typ sap, Int4 **matrix, double PerNats, double misalncut,
	Int4 gapopen, Int4 gapextend, Int4 window, a_type AB)
{
	Int4 trim[2];
	if(MarginalProbTrim(sap,matrix,PerNats,misalncut,gapopen, gapextend,window,trim,AB)){
		return TrimSAP(sap,trim[0],trim[1]);
	} else return NULL;
}

BooLean MarginalProbTrim(sap_typ sap, Int4 **matrix, double PerNats, double misalncut,
	Int4 gapopen, Int4 gapextend, Int4 window, Int4 trim[2], a_type AB)
//returns TRUE if alignment is trimmed; FALSE - otherwise
//matrix starts from 1 
//trim[0]=how many residues from subject seq are trimmed from the begining
//trim[1]=how many residues from subject seq are trimmed from the end
{
	assert(window>0);
	assert(misalncut >= 0.0 && misalncut <= 1.0);
	Int4 i,t,j,pos_q,pos_s;
	e_type qE = QuerySeqGSAP(sap);
	e_type sE = SubjectSeqGSAP(sap);
	Int4 length = LenSeq(qE);
	dsp_typ dsp = sap->segs;
	Int4 numseg = dsp->numseg;
	Int4 *starts= dsp->starts;
	Int4 *lens = dsp->lens;
	Int4 start_Q1,start_S1,end_Q1,end_S1,start_Q2,start_S2,end_Q2,end_S2;
	start_Q1 = starts[0] + 1; 
	start_S1 = starts[1] + 1;
	Int4 e_S = starts[2*numseg-1] + lens[numseg-1];
	Int4 sap_subj_len = e_S - start_S1 + 1;
	Int4 mwindow = MINIMUM(Int4,window,sap_subj_len);
	double **marg;
	for(i=0;i<numseg;i++) {
		if(starts[2*i+1] == -1) continue;	
		if(starts[2*i+1] + lens[i] >= starts[1] + mwindow) break;
	}
	if(starts[2*i] == -1) {
		end_S1 = starts[2*i+1];
		end_Q1 = starts[2*(i+1)];
	} else {
		end_S1 = starts[1] + mwindow;
		end_Q1 = starts[2*i] + starts[1] + mwindow - starts[2*i+1];
	}
	Int4 len_prof = end_Q1-start_Q1+1;
	Int4 len_seq = end_S1-start_S1+1;
	unsigned char *ptr = SeqPtr(sE) + start_S1 - 1;
	marg = LocalMarginalProb(len_seq,ptr,len_prof,matrix+start_Q1-1,PerNats,gapopen,gapextend,AB);

	BooLean ST = TRUE; i=0; t=0;
	Int4 x=1; Int4 y=1;
	while (i<numseg && ST) {
		if(starts[2*i] != -1 && starts[2*i+1] != -1) {
			for(j=1;j<=lens[i];j++) {
				if(marg[x++][y++] >= misalncut) {
					ST = FALSE; break;
				} 
				t++;
				if(t >= mwindow) { ST = FALSE; break; }
			} x--; y--;
		} else if (starts[2*i] == -1) {
			t += lens[i]; y += lens[i]; 
			if(t >= mwindow) { t-= lens[i]; break; }
		} else {
			x += lens[i];
			if(t >= mwindow) break;
		} i++;
	}
	trim[0] = t;
	for(i=1;i<=len_prof;i++) free(marg[i]); free(marg);

	end_S2 = starts[2*numseg-1] + lens[numseg-1];
	end_Q2 = starts[2*numseg-2] + lens[numseg-1];
	for(i=0;i<numseg;i++) {
		if(starts[2*i+1] == -1) continue;
		if(starts[2*i+1] + lens[i] - 1 >= end_S2 - mwindow) break;
	}
	if (starts[2*i] == -1) {
		start_S2 = starts[2*(i+1)+1];
		start_Q2 = starts[2*(i+1)];		
	} else {
		start_S2 = end_S2 - mwindow + 1;
		start_Q2 = starts[2*i] + end_S2 - mwindow - starts[2*i+1] + 1;
	}
	len_prof = end_Q2-start_Q2+1;
	len_seq = end_S2-start_S2+1;
	ptr = SeqPtr(sE) + start_S2 - 1;
	marg = LocalMarginalProb(len_seq,ptr,len_prof,matrix+start_Q2-1,PerNats,gapopen,gapextend,AB);

	x=len_prof; y = len_seq;
	ST = TRUE; i = numseg-1; t=0;
	while (i>=0 && ST) {
		if(starts[2*i] != -1 && starts[2*i+1] != -1) {
			for(j=1;j<=lens[i];j++) {
				if(marg[x--][y--] >= misalncut) {
					ST = FALSE; break;
				} 
				t++; 
				if(t >= mwindow) { ST = FALSE; break; }
			} x++; y++;
		} else if (starts[2*i] == -1) {
			t += lens[i]; y -= lens[i];
			if(t >= mwindow) { t -= lens[i]; break; }
		} else {
			x -= lens[i];
			if(t >= mwindow) break;
		}
		i--;
	}
	trim[1] = t;
	for(i=1;i<=len_prof;i++) free(marg[i]); free(marg);
	if(trim[0] > 0 || trim[1] > 0) return TRUE;
	else return FALSE;
}


double  **LocalMarginalProb(Int4 seq_len, unsigned char *seq, Int4 length, Int4 **matrix, double PerNats,
	Int4 gapopen, Int4 gapextend, a_type AB) 
//2*log2(odd) = C*ln(odd);  C=2/ln(2); so for blast PerNats shold be 2/ln(2);
{
        Int4            i,j;
	Int4		*matrix_i;
        double **pmat_emit_odd;
        NEWP(pmat_emit_odd,length+2,double);
        for(i=0;i<=length+1;i++){
                NEW(pmat_emit_odd[i],nAlpha(AB)+2,double);
        }
        double *pmat_emit_odd_i,*WtFreq_i;
        for(i=1;i<=length;i++){
		matrix_i = matrix[i];
                pmat_emit_odd_i = pmat_emit_odd[i];
                for(j=1;j<=nAlpha(AB);j++){
                	pmat_emit_odd_i[j] = exp((double) matrix_i[j]/PerNats);
                }
                pmat_emit_odd_i[0] = 0.5;
        }
	double **TMAT = MargProbCore(seq_len,seq,length,pmat_emit_odd,PerNats,gapopen,gapextend,AB);
	for(i=0;i<=length+1; i++) free(pmat_emit_odd[i]); free(pmat_emit_odd);
        return TMAT;
}

double **MargProbCore(Int4 seq_len, unsigned char *seq, Int4 length, double **pmat_emit_odd, double PerNats,
                Int4 gapopen, Int4 gapextend, a_type AB)
{
	double **TMAT;
	Int4 	i,im1,j,jm1;
        unsigned char   *rseq;
	NEW(rseq,seq_len+1,unsigned char);
	for(i=1;i<=seq_len;i++) rseq[i] = seq[seq_len-i+1];
        double **SMAT,**BMAT,**SDEL,**BDEL,**SINS,**BINS;
        NEWP(SMAT,length+2,double); NEWP(BMAT,length+2,double); NEWP(TMAT,length+2,double);
        NEWP(SDEL,length+2,double); NEWP(BDEL,length+2,double);
        NEWP(SINS,length+2,double); NEWP(BINS,length+2,double);
        for(i=0;i<=length+1;i++){
                NEW(SMAT[i],seq_len+1,double); NEW(BMAT[i],seq_len+1,double); 
                NEW(TMAT[i],seq_len+1,double);
                NEW(SDEL[i],seq_len+1,double); NEW(BDEL[i],seq_len+1,double);
                NEW(SINS[i],seq_len+1,double); NEW(BINS[i],seq_len+1,double);
        }
        double m2d,m2i,d2d,i2i;
        m2d = m2i = exp((double) -gapopen/PerNats);
        d2d = i2i = exp((double) -gapextend/PerNats);
	double **rmat_emit_odd;
        NEWP(rmat_emit_odd,length+2,double);
        for(i=1;i<=length;i++){ rmat_emit_odd[i] = pmat_emit_odd[length-i+1]; }

        SDEL[1][1] = 0; SMAT[1][1] = pmat_emit_odd[1][seq[1]];
        for(jm1=1,j=2;j<=length;jm1++,j++) {
                SINS[j][1] = 0; SDEL[j][1] = 0;
                SMAT[j][1] = pmat_emit_odd[j][seq[1]];
        }
        SINS[1][1] = 0; SMAT[1][2] = pmat_emit_odd[1][seq[2]];
        SINS[1][2] = pmat_emit_odd[1][seq[1]] * m2i; 
        SDEL[1][2] = 0;
        for(im1=2,i=3;i<=seq_len;im1++,i++) {
                SDEL[1][i] = 0; SMAT[1][i] = pmat_emit_odd[1][seq[i]];
                SINS[1][i] = (SINS[1][im1] * i2i) + (SMAT[1][im1] * m2i); 
        }

        BDEL[1][1] = 0; BMAT[1][1] = rmat_emit_odd[1][rseq[1]];
        for(jm1=1,j=2;j<=length;jm1++,j++) {
                BINS[j][1] = 0; BDEL[j][1] = 0;
                BMAT[j][1] = rmat_emit_odd[j][rseq[1]];
        }                        
        BINS[1][1] = 0; BMAT[1][2] = rmat_emit_odd[1][rseq[2]];
        BINS[1][2] = rmat_emit_odd[1][rseq[1]] * m2i;
        BDEL[1][2] = 0;
        for(im1=2,i=3;i<=seq_len;im1++,i++) {
                BDEL[1][i] = 0; BMAT[1][i] = rmat_emit_odd[1][rseq[i]];
                BINS[1][i] = (BINS[1][im1] * i2i) + (BMAT[1][im1] * m2i);
        }

        double *smatj,*sdelj,*sinsj,*bmatj,*bdelj,*binsj;
        double *smatjm1,*sdeljm1,*sinsjm1,*bmatjm1,*bdeljm1,*binsjm1;
        double *pmat_emit_oddj,*rmat_emit_oddj,pmat_emit_oddjres,rmat_emit_oddjres;
        unsigned char res,rres;
        for(jm1=1,j=2;j<=length;j++,jm1++){
                smatj=SMAT[j];sdelj=SDEL[j];sinsj=SINS[j];
                bmatj=BMAT[j];bdelj=BDEL[j];binsj=BINS[j];
                smatjm1=SMAT[jm1];sdeljm1=SDEL[jm1];sinsjm1=SINS[jm1];
                bmatjm1=BMAT[jm1];bdeljm1=BDEL[jm1];binsjm1=BINS[jm1];
                pmat_emit_oddj=pmat_emit_odd[j];
                rmat_emit_oddj=rmat_emit_odd[j];
                for(im1=1,i=2;i<=seq_len;im1++,i++){
                        res=seq[i]; rres=rseq[i];
                        pmat_emit_oddjres=pmat_emit_oddj[res]; 
                        rmat_emit_oddjres=rmat_emit_oddj[rres];
                        smatj[i] += (smatjm1[im1] * pmat_emit_oddjres);
                        smatj[i] += (sdeljm1[i] * pmat_emit_oddjres); 
                        smatj[i] += (sinsjm1[im1] * pmat_emit_oddjres);
                        smatj[i] += pmat_emit_oddjres;
                        sdelj[i] += (smatjm1[im1] * m2d);
                        sdelj[i] += (sdeljm1[i] * d2d);
                        sinsj[i] += (smatj[im1] * m2i);
                        sinsj[i] += (sinsj[im1] * i2i);
                        bmatj[i] += (bmatjm1[im1] * rmat_emit_oddjres);
                        bmatj[i] += (bdeljm1[i] * rmat_emit_oddjres);
                        bmatj[i] += (binsjm1[im1] * rmat_emit_oddjres);
                        bmatj[i] += rmat_emit_oddjres;
                        bdelj[i] += (bmatjm1[im1] * m2d);
                        bdelj[i] += (bdeljm1[i] * d2d);
                        binsj[i] += (bmatj[im1] * m2i);
                        binsj[i] += (binsj[im1] * i2i);
                }
        }
        double *TMAT_j,*SMAT_j,*BMAT_l,*pmat_emit_odd_j;
        for(j=1;j<=length;j++){
                TMAT_j = TMAT[j]; SMAT_j = SMAT[j]; BMAT_l = BMAT[length-j+1];
                pmat_emit_odd_j = pmat_emit_odd[j];
                for(i=1;i<=seq_len;i++){
                        res=seq[i];
                        TMAT_j[i]=BMAT_l[seq_len-i+1]*(SMAT_j[i]/pmat_emit_odd_j[res]);
                }
        }
        double all = 0.;
        for(j=1;j<=length;j++){
                smatj=SMAT[j];
                for(i=1;i<=seq_len;i++){
                        all += smatj[i]; 
                }
        }
        for(j=1;j<=length;j++){
                TMAT_j = TMAT[j];
                for(i=1;i<=seq_len;i++){
                        TMAT_j[i] /= all;
//printf("TMAT[%d][%d] = %lf\n",j,i,TMAT_j[i]);
                }
        }
        for(i=0;i<=length+1;i++) { 
                free(SMAT[i]);free(SDEL[i]);free(SINS[i]);
                free(BMAT[i]);free(BDEL[i]);free(BINS[i]);
        }
        free(SMAT);free(SDEL);free(SINS);free(BMAT);free(BDEL);free(BINS);
        free(rmat_emit_odd); free(rseq);
        if(!(all > 0.0)){
          for(i=0;i<=length+1;i++) free(TMAT[i]); free(TMAT);
          return 0;
        } free(TMAT[0]); TMAT[0]=0;
        free(TMAT[length+1]); TMAT[length+1]=0;
	return TMAT;
}

double  **LocalMarginalProb2(Int4 seq_len, unsigned char *seq, char *operation, Int4 length, 
	double **WtFreq, double *fract_seq_aln, double *bfreq, 
		double PerNats, Int4 gapopen, Int4 gapextend, a_type AB) 
{
        Int4            i,j;
        double **pmat_emit_odd;
        NEWP(pmat_emit_odd,length+2,double);
        for(i=0;i<=length+1;i++){
                NEW(pmat_emit_odd[i],nAlpha(AB)+2,double);
        }
        double *pmat_emit_odd_i,*WtFreq_i;
        for(i=1;i<=length;i++){
                pmat_emit_odd_i = pmat_emit_odd[i];
		WtFreq_i = WtFreq[i];
		if(fract_seq_aln[i] > 0.0) {
                	for(j=1;j<=nAlpha(AB);j++){
                        	pmat_emit_odd_i[j] = WtFreq_i[j]/bfreq[j];
                	}
		} else {
			for(j=1;j<=nAlpha(AB);j++){
				pmat_emit_odd_i[j] = 1.0;
			}
		}
		pmat_emit_odd_i[0] = 0.5;
        }
	double **TMAT = MargProbCore(seq_len,seq,length,pmat_emit_odd,PerNats,gapopen,gapextend,AB);
        for(i=0;i<=length+1; i++) free(pmat_emit_odd[i]); free(pmat_emit_odd);
        return TMAT;
}

double *ScoreMargOperArr(double **TMAT, char *operation, Int4 Start, Int4 Oper_len)
{
	double *score;
	NEW(score,Oper_len,double);
	Int4 i=1,j=1;
	Int4 posProf=1, posSeq=Start;
	while(operation[i] != 'E'){
	    switch(operation[i++]){
		case 'M':
		case 'm': score[j++] = TMAT[posProf++][posSeq++]; break;
		case 'D':
		case 'd': j++; posProf++; break;
		case 'I':
		case 'i': j++; posSeq++; break;
		default : print_error("ScoreMargOperArr(): error in the operation array");
	    }
	}
	return score;
}

double  **MarginalProb2(e_type sE, Int4 length, Int4 **mat_emit_odd, Int4 PerNats, Int4 gapopen,
		Int4 gapextend, a_type AB)
{
        double          **TMAT;
        Int4            i,im1,j,jm1,seq_len=LenSeq(sE);
        unsigned char   *seq=SeqPtr(sE);
        e_type          rE=CopySeq(sE);


        rE = ReverseSeq(rE);
        unsigned char   *rseq=XSeqPtr(rE);
        double **SMAT,**BMAT,**SDEL,**BDEL,**SINS,**BINS;
        NEWP(SMAT,length+2,double); NEWP(BMAT,length+2,double); NEWP(TMAT,length+2,double);
        NEWP(SDEL,length+2,double); NEWP(BDEL,length+2,double);
        NEWP(SINS,length+2,double); NEWP(BINS,length+2,double);
        for(i=0;i<=length+1;i++){
                NEW(SMAT[i],seq_len+1,double); NEW(BMAT[i],seq_len+1,double); 
                NEW(TMAT[i],seq_len+1,double);
                NEW(SDEL[i],seq_len+1,double); NEW(BDEL[i],seq_len+1,double);
                NEW(SINS[i],seq_len+1,double); NEW(BINS[i],seq_len+1,double);
        }

        double m2d,m2i,d2d,i2i;

        m2d = m2i = exp((double) -gapopen/PerNats);
        d2d = i2i = exp((double) -gapextend/PerNats);

        double **pmat_emit_odd,**rmat_emit_odd;

        NEWP(pmat_emit_odd,length+2,double);
        NEWP(rmat_emit_odd,length+2,double);
        for(i=0;i<=length+1;i++){
                NEW(pmat_emit_odd[i],nAlpha(AB)+2,double);
        }
        double *pmat_emit_odd_i;
        Int4 *mat_emit_odd_i;
        for(i=1;i<=length;i++){
                pmat_emit_odd_i = pmat_emit_odd[i];
                mat_emit_odd_i = mat_emit_odd[i];
                for(j=1;j<=nAlpha(AB);j++){
                        pmat_emit_odd_i[j]=exp((double) mat_emit_odd_i[j]/PerNats);
                }
        }

        for(i=1;i<=length;i++){ rmat_emit_odd[i] = pmat_emit_odd[length-i+1]; }

        SDEL[1][1] = m2d; SMAT[1][1] = pmat_emit_odd[1][seq[1]];
        for(jm1=1,j=2;j<=length;jm1++,j++) {
                SINS[j][1] = 0; SDEL[j][1] = SDEL[jm1][1] * d2d;
                SMAT[j][1] = SDEL[jm1][1] * pmat_emit_odd[j][seq[1]];
        }
        SINS[1][1] = 0; SMAT[1][2] = m2i * pmat_emit_odd[1][seq[2]];
        SINS[1][2] = pmat_emit_odd[1][seq[1]] * m2i; 
        SDEL[1][2] = 0;
        for(im1=2,i=3;i<=seq_len;im1++,i++) {
                SDEL[1][i] = 0; SMAT[1][i] = m2i * pow(i2i,i-2) * pmat_emit_odd[1][seq[i]];
                SINS[1][i] = (SINS[1][im1] * i2i) + (SMAT[1][im1] * m2i); 
        }

        BDEL[1][1] = m2d; BMAT[1][1] = rmat_emit_odd[1][rseq[1]];
        for(jm1=1,j=2;j<=length;jm1++,j++) {
                BINS[j][1] = 0; BDEL[j][1] = BDEL[jm1][1] * d2d;
                BMAT[j][1] = BDEL[jm1][1] * rmat_emit_odd[j][rseq[1]];
        }                        
        BINS[1][1] = 0; BMAT[1][2] = m2i * rmat_emit_odd[1][rseq[2]];
        BINS[1][2] = rmat_emit_odd[1][rseq[1]] * m2i;
        BDEL[1][2] = 0;
        for(im1=2,i=3;i<=seq_len;im1++,i++) {
                BDEL[1][i] = 0; BMAT[1][i] = m2i * pow(i2i,i-2) * rmat_emit_odd[1][rseq[i]];
                BINS[1][i] = (BINS[1][im1] * i2i)  + (BMAT[1][im1] * m2i);
        }

        double *smatj,*sdelj,*sinsj,*bmatj,*bdelj,*binsj;
        double *smatjm1,*sdeljm1,*sinsjm1,*bmatjm1,*bdeljm1,*binsjm1;
        double *pmat_emit_oddj,*rmat_emit_oddj,pmat_emit_oddjres,rmat_emit_oddjres;
        unsigned char res,rres;
        for(jm1=1,j=2;j<=length;j++,jm1++){
                smatj=SMAT[j];sdelj=SDEL[j];sinsj=SINS[j];
                bmatj=BMAT[j];bdelj=BDEL[j];binsj=BINS[j];
                smatjm1=SMAT[jm1];sdeljm1=SDEL[jm1];sinsjm1=SINS[jm1];
                bmatjm1=BMAT[jm1];bdeljm1=BDEL[jm1];binsjm1=BINS[jm1];
                pmat_emit_oddj=pmat_emit_odd[j];
                rmat_emit_oddj=rmat_emit_odd[j];
                for(im1=1,i=2;i<=seq_len;im1++,i++){
                        res=seq[i]; rres=rseq[i];
                        pmat_emit_oddjres=pmat_emit_oddj[res]; 
                        rmat_emit_oddjres=rmat_emit_oddj[rres];
                        smatj[i] += (smatjm1[im1] * pmat_emit_oddjres);
                        smatj[i] += (sdeljm1[i] * pmat_emit_oddjres); 
                        smatj[i] += (sinsjm1[im1] * pmat_emit_oddjres);
                        sdelj[i] += (smatjm1[im1] * m2d);
                        sdelj[i] += (sdeljm1[i] * d2d);
                        sinsj[i] += (smatj[im1] * m2i);
                        sinsj[i] += (sinsj[im1] * i2i);
                        bmatj[i] += (bmatjm1[im1] * rmat_emit_oddjres);
                        bmatj[i] += (bdeljm1[i] * rmat_emit_oddjres);
                        bmatj[i] += (binsjm1[im1] * rmat_emit_oddjres);
                        bdelj[i] += (bmatjm1[im1] * m2d);
                        bdelj[i] += (bdeljm1[i] * d2d);
                        binsj[i] += (bmatj[im1] * m2i);
                        binsj[i] += (binsj[im1] * i2i);
                }
        }
        double *TMAT_j,*SMAT_j,*BMAT_l,*pmat_emit_odd_j;
        for(j=1;j<=length;j++){
                TMAT_j = TMAT[j]; SMAT_j = SMAT[j]; BMAT_l = BMAT[length-j+1];
                pmat_emit_odd_j = pmat_emit_odd[j];
                for(i=1;i<=seq_len;i++){
                        res=seq[i];
                        TMAT_j[i]=BMAT_l[seq_len-i+1]*(SMAT_j[i]/pmat_emit_odd_j[res]);
                }
        }

        double sdel=0;
        double sins=SINS[length][seq_len];
        // fprintf(stderr,"SMAT=%g  SDEL=%g  SINS=%g\n",SMAT[length][seq_len],sdel,sins);

        double last_cell=SMAT[length][seq_len]+sdel+sins;
        for(j=1;j<=length;j++){
                TMAT_j = TMAT[j];
                for(i=1;i<=seq_len;i++){
                        TMAT_j[i]/=last_cell;
#if 0
if(TMAT[j][i] > 0.01) printf("TMAT[%d][%d]=%g\n",j,i,TMAT[j][i]);
//printf("TMAT[%d][%d]=%g SMAT[%d][%d]=%gBMAT[%d][%d]=%g\n",
 j,i,TMAT[j][i],j,i,SMAT[j][i],length-j+1,seq_len-i+1,BMAT[length-j+1][seq_len-i+1]);
#endif

//printf("TMAT[%d][%d]=%g\n",j,i,TMAT[j][i]);
                }
        }

        for(i=0;i<=length+1;i++) { 
                free(SMAT[i]);free(SDEL[i]);free(SINS[i]);
                free(BMAT[i]);free(BDEL[i]);free(BINS[i]);
        }
        free(SMAT);free(SDEL);free(SINS);free(BMAT);free(BDEL);free(BINS);
        for(i=0;i<=length+1;i++) { free(pmat_emit_odd[i]); } 
        free(pmat_emit_odd); free(rmat_emit_odd); 
	NilSeq(rE);
        if(!(last_cell > 0.0)){
          for(i=0;i<=length+1;i++) free(TMAT[i]); free(TMAT);
          return 0;
        } free(TMAT[0]); TMAT[0]=0; // AFN: fixes memory leak.
        free(TMAT[length+1]); TMAT[length+1]=0; // AFN: fixes memory leak.
        return TMAT;
}
