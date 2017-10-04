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

#include "dp.h"

dp_typ::dp_typ() 
{
}

dp_typ::dp_typ(e_type E1, smx_typ M1, Int4 *inso, Int4 *inse, Int4 *delo, Int4 *dele)
{
	Mkdp_typ(E1, M1, inso, inse, delo, dele);
}

void dp_typ::Mkdp_typ(e_type E1, smx_typ M1, Int4 *inso, Int4 *inse, Int4 *delo, Int4 *dele)
{
        Int4    i,j,k,l,im1,jm1,r,s,s0,s1,v,temp=0;
	E = E1;
	M = M1;
	seq_len = LenSeq(E);
	unsigned char *seq = SeqPtr(E);	
        prof_len=LenSMatrix(M); 
        alph_len = nAlpha(SMatrixA(M));

	NEW(io,prof_len+2,Int4);
	NEW(ie,prof_len+2,Int4);
	NEW(od,prof_len+2,Int4);
	NEW(de,prof_len+2,Int4);

	for(i=1; i<=prof_len; i++){
		io[i] = inso[i]; ie[i] = inse[i];
		od[i] = delo[i]; de[i] = dele[i];
	}
	for(i=1; i<=prof_len; i++){
		io[i-1] = temp;
		temp = (io[i]+io[i+1])/2;
	}
	temp=0;
	for(i=1; i<=prof_len; i++){
		ie[i-1] = temp;
		temp = (ie[i]+ie[i+1])/2;
	}
	
        NEWP(MAT,seq_len+1,Int4); NEWP(DEL,seq_len+1,Int4);
        NEWP(INS,seq_len+1,Int4);
        for(i=0;i<=seq_len;i++){
                NEW(MAT[i],prof_len+1,Int4); NEW(DEL[i],prof_len+1,Int4);
                NEW(INS[i],prof_len+1,Int4);
        }
        NEWP(match,alph_len+1,Int4);
        for(i=0;i<=alph_len;i++) NEW(match[i],prof_len+1,Int4);

        for(s=1;s<=prof_len;s++){
                for(r=1;r<=alph_len;r++){
                        match[r][s]=ValSMatrix(s,r,M);
                }
        }
        DEL[0][1]=INS[0][1]=MAT[0][1]=od[1]+de[1];
        for(jm1=1,j=2;j<=prof_len;jm1++,j++) {
                MAT[0][j]=DEL[0][jm1]+od[j]+de[j]; 
                DEL[0][j]=DEL[0][jm1]+od[j]+de[j];
                INS[0][j]=DEL[0][jm1]+od[j]+de[j];
        }
        for(i=1;i<=seq_len;i++) {DEL[i][0] = od[1];}

	NEW(last_column,seq_len+2,Int4);
	NEW(max_last,seq_len+2,Int4);
	NEW(end_state,seq_len+2,char);
	NEW(loc_max,seq_len+2,Int4);

        for(im1=0,i=1;i<=seq_len;im1++,i++){
            register Int4 J,Jm1;
            register Int4 *mati= MAT[i];
            register Int4 *matim1= MAT[im1];
            register Int4 *insi= INS[i];
            register Int4 *insim1= INS[im1];
            register Int4 *deli= DEL[i];
            r = seq[i];
            for(J=1,Jm1=0;J<=prof_len;J++,Jm1++){ // INNER LOOP!!!
                    mati[J] = match[r][J]+MAXIMUM(Int4,
                        MAXIMUM(Int4,matim1[Jm1],DEL[im1][Jm1]),insim1[Jm1]);
                    s0 = od[J]+de[J]+mati[Jm1];
                    if(s0 > (s1=de[J]+deli[Jm1])) deli[J]=s0; else deli[J]=s1;
                    s0 = io[J]+ie[J]+matim1[J];
                    if(s0 > (s1=ie[J]+insim1[J])) insi[J]=s0; else insi[J]=s1;
            }
	    if(mati[prof_len] >= deli[prof_len]){
		last_column[i] = mati[prof_len];
		end_state[i] = 'M';
	    }
	    else{ 
		last_column[i] = deli[prof_len];
		end_state[i] = 'D';
	    }
        }
        dh_type H;
        H = dheap(seq_len+1,4);
        for(i=1;i<=seq_len;i++) insrtHeap(i,(keytyp) last_column[i],H);
        for(i=seq_len;i>=1;i--){
		max_last[i] = delminHeap(H);
	}
        Nildheap(H);
	for(i=1;i<=seq_len;i++){
		if((last_column[i-1] <= last_column[i]) && (last_column[i] >= last_column[i+1]))
			loc_max[i] = last_column[i];
		else loc_max[i] = -999;
	}
}

dp_typ::~dp_typ()
{
	Free();
}

void dp_typ::Free()
{
	Int4 i,j;
	if(MAT != NULL){
		for(i=0;i<=seq_len;i++){
			free(MAT[i]);
		}
		free(MAT); 
	}
        if(DEL != NULL){
                for(i=0;i<=seq_len;i++){  
                        free(DEL[i]);
                }
                free(DEL);    
        }
        if(INS != NULL){
                for(i=0;i<=seq_len;i++){  
                        free(INS[i]);
                }
                free(INS);    
        }
	if(match != NULL){
                for(i=0;i<=alph_len;i++){
                        free(match[i]);
                }
                free(match);  
        }
	if(last_column != NULL) free(last_column);
	if(max_last != NULL) free(max_last);
	if(end_state != NULL) free(end_state);
	if(loc_max != NULL) free(loc_max);
	if(io != NULL) free(io);
	if(ie != NULL) free(ie);
	if(od != NULL) free(od);
	if(de != NULL) free(de);
}
char *dp_typ::TraceBackWRTpos(Int4 end_of_aln, Int4 *alignscore, Int4 *start_of_aln, Int4 *oper_len)
{
	Int4	i,j,k;
	char 	state, *operation, *back_operation;

        NEW(operation,seq_len+prof_len+3,char);
        NEW(back_operation,seq_len+prof_len+3,char);
	*alignscore = last_column[end_of_aln];
	i = end_of_aln, j = prof_len;
        for(state='E',k=0;state != 'X'; ){
           switch(state){
                case 'B': // begin state (exit traceback).
                        back_operation[k++]='E'; *start_of_aln=i; state='X'; 
			break;
                case 'E': // end state.  
                        back_operation[k++]='E';
			if (end_state[end_of_aln]=='M') state='M';
			else if(end_state[end_of_aln]=='D') state='D';
			else print_error("TraceBackWRTpos(): unknown end_state");
                        break;
                case 'M': // previously sampled matched state.
                        if(j <= 1) { back_operation[k++]='M'; state='B'; break; }
                        else { // else we are within a block & use affine gaps.
                           j--; i--; back_operation[k++]='m';
                           if(MAT[i][j]>INS[i][j] && MAT[i][j]>DEL[i][j]) state='M';
                           else {if(INS[i][j]>DEL[i][j]) state='I';
                                 else state='D';
                           }
                        }
                        break;                          
                case 'D': // previously sampled deletion
                          if(j <= 1){ back_operation[k++]='D'; state='B'; i++; break; }
                          else { j--; back_operation[k++]='d';
                             if(DEL[i][j]+de[j]>=MAT[i][j]+od[j]+de[j]) state='D';
                             else state='M';
                          }
                          break;
                case 'I': // previously sampled insertion.
                          back_operation[k++]='I'; i--; 
                          if (INS[i][j]+ie[j]>MAT[i][j]+io[j]+ie[j]) state='I';
                          else state='M';
                          break;
                default:  print_error("TraceBackWRTpos(): this should not happen");
           }
        }
        *oper_len = k-2;                           
        for(i=0;i<k;i++) operation[i]=back_operation[k-i-1];
        free(back_operation);
        return operation;
}

Int4 *dp_typ::LastColumn()
{
	return last_column;
}

Int4 *dp_typ::MaxLast() 
{       
        return max_last;
}

Int4 *dp_typ::LocMax()
{
        return loc_max;
}

void dp_typ::PrintGaps()
{
	Int4 i;
	printf("POS: ins_open ins_extd del_open del_extnd\n");
	for(i=1;i<=prof_len;i++){
		printf("%-3d: %8d %8d %8d %9d\n",i,-io[i],-ie[i],-od[i],-de[i]);
	}
}

char **MaximizeSum(e_type E, smx_typ M1, Int4 *io, Int4 *ie, Int4 *od, Int4 *de, Int4 start1, Int4 end1,
	Int4 start2, Int4 end2, Int4 tflank, Int4 *fstart1, Int4 *fend1, Int4 *fstart2, Int4 *fend2,
		Int4 *foper_len1, Int4 *foper_len2, Int4 *falignscore1, Int4 *falignscore2)
{
	char 	**operation;
	char	*temp1, *temp2;
	Int4	i, max_score=-9999, flank;
	Int4 	alignscore1, alignscore2;
	Int4	start_of_aln1, start_of_aln2;
	Int4	oper_len1, oper_len2;
	Int4	end_aln1, end_aln2;
	Int4 	*pos_seq1, *pos_seq2, *pos_mas1, *pos_mas2;
	Int4	*nio, *nie, *nod, *nde;
	Int4 	len = LenSMatrix(M1);

	smx_typ M2 = ReverseSMatrix(M1);

	if(start2+tflank >= end2 || end1-tflank <= start1){
			if(end2-start2-1 < end1-start1-1)
				flank = end2-start2-1;
			else flank = end1-start1-1;
	}
	else flank = tflank;

	e_type E1 = MkSubSeq(start1,start2+flank,E);
	e_type E2 = ReverseSeq(MkSubSeq(end1-flank,end2,E));

	NEWP(operation,3,char);
//	NEW(operation[1],LenSeq(E1)+len+2,char);
//	NEW(operation[2],LenSeq(E2)+len+2,char);

	NEW(pos_seq1,LenSeq(E)+1,Int4);
	NEW(pos_seq2,LenSeq(E)+1,Int4);
	NEW(pos_mas1,LenSeq(E1)+1,Int4);
	NEW(pos_mas2,LenSeq(E2)+1,Int4); 

	for(i=1;i<=LenSeq(E);i++) pos_seq1[i] = i-start1+1;
	for(i=1;i<=LenSeq(E);i++) pos_seq2[i] = end2-i+1;

	for(i=1;i<=LenSeq(E1);i++) pos_mas1[i] = start1+i-1;
	for(i=1;i<=LenSeq(E2);i++) pos_mas2[i] = end2-i+1;

	nio = ReversePen(io,len); nie = ReversePen(ie,len); 
	nod = ReversePen(od,len); nde = ReversePen(de,len);

	dp_typ d1(E1,M1,io,ie,od,de);
	dp_typ d2(E2,M2,nio,nie,nod,nde);

	for(end_aln1 = pos_seq1[end1-flank]; end_aln1 <= pos_seq1[start2+flank]; end_aln1++){ 
		temp1 = d1.TraceBackWRTpos(end_aln1,&alignscore1,&start_of_aln1,&oper_len1);
		for(end_aln2 = pos_seq2[pos_mas1[end_aln1]]-1; end_aln2 >= pos_seq2[start2+flank]; end_aln2--){
			temp2 = d2.TraceBackWRTpos(end_aln2,&alignscore2,&start_of_aln2,&oper_len2);
			if((alignscore1+alignscore2) > max_score){
				max_score = alignscore1+alignscore2;
				*fstart1 = pos_mas1[start_of_aln1]; *fstart2 = pos_mas2[end_aln2];
				*fend1 = pos_mas1[end_aln1]; *fend2 = pos_mas2[start_of_aln2];
				*foper_len1 = oper_len1; *foper_len2 = oper_len2;
				*falignscore1 = alignscore1; *falignscore2 = alignscore2;
			}
			free(temp2);
		}
		free(temp1);
	}
	operation[1]=d1.TraceBackWRTpos(pos_seq1[*fend1],&alignscore1,&start_of_aln1,&oper_len1);
	temp2=d2.TraceBackWRTpos(pos_seq2[*fstart2],&alignscore2,&start_of_aln2,&oper_len2);
	operation[2] = ReverseOperArray(temp2, *foper_len2);
	free(temp2);
	free(pos_seq1); free(pos_seq2); free(pos_mas1); free(pos_mas2);
	free(nio); free(nie); free(nod); free(nde);
	NilSeq(E1); NilSeq(E2);
	NilSMatrix(M2);
	return operation;
}

char **DoAlignm(e_type E, smx_typ M, Int4 *io, Int4 *ie, Int4 *od, Int4 *de, Int4 flank, 
	Int4 cutoff, Int4 *nstarts, Int4 *nends, Int4 *nalignscores, Int4 *noper_lens, Int4 *n_rpts)
{
	Int4	*last_column, *max_last, *loc_max;
	Int4	i=1, k=2, j, m, signal;
	Int4 	seq_len = LenSeq(E);
	Int4 	prof_len = LenSMatrix(M);
	Int4 	alignscore, start, oper_len;
	Int4 	fstart_j,fstart_i,fend_j,fend_i;
	Int4	foper_len_j,foper_len_i,falignscore_j,falignscore_i;
	Int4	*starts, *ends, *alignscores, *oper_lens;
	char 	**operation, **toperation, **noperation, *temp_oper;
	dp_typ 	d(E,M,io,ie,od,de);

	if (cutoff < 1) print_error("DoAlignm(): low cutoff");
	if (flank > prof_len/4) print_error("DoAlignm(): flank too big");

	NEWP(operation,seq_len+1,char); 
	NEWP(noperation,seq_len+1,char);

	NEW(starts,seq_len+1,Int4);
	NEW(ends,seq_len+1,Int4);
	NEW(alignscores,seq_len+1,Int4);
	NEW(oper_lens,seq_len+1,Int4);

	last_column = d.LastColumn();
	max_last = d.MaxLast();
	loc_max = d.LocMax();

	operation[1] = d.TraceBackWRTpos(max_last[1],&alignscore,&start,&oper_len);
	starts[1] = start, ends[1] = max_last[1]; oper_lens[1] = oper_len;
	alignscores[1] = alignscore;
	i = 2;
	while(i <= seq_len){
		m = max_last[i];
		if(loc_max[m] >= cutoff){
			signal = 1;
			temp_oper = d.TraceBackWRTpos(m,&alignscore,&start,&oper_len);
			operation[k] = temp_oper; alignscores[k] = alignscore;
			starts[k] = start, ends[k] = m; oper_lens[k] = oper_len;
			for(j=1; j<=k-1; j++){
				if((start >= starts[j] && start <= ends[j]-flank) 
					|| (starts[j] >= start && starts[j] <= m-flank)
						|| (start >= starts[j] && m <= ends[j])
							|| (starts[j] >= start && ends[j] <= m)){
					free(temp_oper); signal = 0; break;
				}
			}
			if (signal){
				k++;
			}
		}
		i++;
	}
	*n_rpts = k-1;

	for(i=1;i<=*n_rpts;i++){
		for(j=1;j<i;j++){
			if((starts[i] >= ends[j]-flank && starts[i] <= ends[j])){
				toperation = MaximizeSum(E,M,io,ie,od,de,starts[j],ends[j],starts[i],ends[i],
					flank,&fstart_j,&fend_j,&fstart_i,&fend_i,&foper_len_j,&foper_len_i,
						&falignscore_j,&falignscore_i);
				operation[j] = toperation[1]; operation[i] = toperation[2];
				ends[j] = fend_j; ends[i] = fend_i;
				starts[j] = fstart_j; starts[i] = fstart_i;
				oper_lens[j] = foper_len_j; oper_lens[i] = foper_len_i;
				alignscores[j] = falignscore_j; alignscores[i] = falignscore_i;
			}
			else if((starts[j] >= ends[i]-flank && starts[j] <= ends[i])){
                                toperation = MaximizeSum(E,M,io,ie,od,de,starts[i],ends[i],starts[j],ends[j],
                                        flank,&fstart_i,&fend_i,&fstart_j,&fend_j,&foper_len_i,&foper_len_j,
                                                &falignscore_i,&falignscore_j);
                                operation[i] = toperation[1]; operation[j] = toperation[2];
                                ends[i] = fend_i; ends[j] = fend_j;
                                starts[i] = fstart_i; starts[j] = fstart_j;
                                oper_lens[i] = foper_len_i; oper_lens[j] = foper_len_j;
                                alignscores[i] = falignscore_i; alignscores[j] = falignscore_j;
                        }
		}
	}

        dh_type H;
	Int4 ky;
        H = dheap(*n_rpts+1,4);
        for(i=1;i<=*n_rpts;i++) insrtHeap(i,(keytyp) starts[i],H);
        for(i=1;i<=*n_rpts;i++){
		ky = delminHeap(H);
                noperation[i] = operation[ky];
		nstarts[i] = starts[ky];
		nends[i] = ends[ky];
		nalignscores[i] = alignscores[ky];
		noper_lens[i] = oper_lens[ky];
        }
        Nildheap(H);

	for(i = *n_rpts+1; i <= seq_len; i++) {free(noperation[i]);}
	free(operation);
	free(starts); free(ends); free(alignscores); free(oper_lens);
	return noperation;
}
