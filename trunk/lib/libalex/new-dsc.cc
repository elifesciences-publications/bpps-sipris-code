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

#include "new-dsc.h"

dsc_typ	MkDSC(FILE *fp)
{
	Int4 i;
	dsc_typ	D;
	NEW(D,1,dsc_type);
	D->error_line=0;
	D->sequence_length=0;
	D->hom_len = 0;
        D->filter_level = 1;    /* Default level of smoothing, as used in paper */
        D->clean = 1;              /* Default remove areas of bad alignment */
        D->clean_length = 40;      /* Default size of area used for cleaning */
        D->clean_percent = 20;     /* Default percent identity for cleaning */
	D->NumHomSeq = Max_hom_seqs;	// AFN fix...
	NEW(D->att_array,Max_length+2,dsc_att_vector);
	NEWP(D->sequence,Max_hom_seqs+1,char);
		for(i=0;i<=Max_hom_seqs;i++) NEW(D->sequence[i],Max_length,char);
	NEWP(D->names,Max_hom_seqs+1,char);
                for(i=0;i<=Max_hom_seqs;i++) NEW(D->names[i],500,char);
	read_clustalw_sequence(fp, D);
	return D;
}

dsc_typ MkDSC(Int4 block, Int4 left_flank, Int4 right_flank, cma_typ cma)
{
	dsc_typ 	D;
	Int4		i,j,k,s,t,len,l,pos[3],read_length=0,true_site,bestSq;
	Int4		sq;
	double		prob,best_prob;
	e_type 		fakeE,trueE;
	unsigned char 	*ptrF, *ptrT;
	gss_typ 	*gss = gssCMSA(cma);
	a_type		A=AlphabetCMSA(cma);
	char		c,*tmp;

	NEW(D,1,dsc_type);
	D->hom_len = NumSeqsCMSA(cma);
	D->NumHomSeq = D->hom_len;
	D->error_line=0; D->filter_level = 1;
	D->clean = 1; D->clean_length = 40; D->clean_percent = 20;
	NEW(D->att_array,Max_length+2,dsc_att_vector);
	D->names=NULL;
        NEWP(D->sequence,NumSeqsCMSA(cma)+2,char);
                for(i=0;i<=NumSeqsCMSA(cma);i++) NEW(D->sequence[i],Max_length,char);
	read_length = LengthCMSA(block,cma);
	// find best sequence and put it first.
	best_prob=GetProbCMSA(block,1,cma); bestSq=1;
	for(i=2; i <= NumSeqsCMSA(cma); i++){
	    prob = GetProbCMSA(block,i,cma);
	    if(prob > best_prob) { best_prob=prob; bestSq=i; }
	} //bestSq=1;  // for old method...
	for(sq=0,i=1; i <= NumSeqsCMSA(cma); i++,sq++){
		ptrF = SeqPtr(gss->FakeSeq(i)); 
		trueE = gss->TrueSeq(i); ptrT = SeqPtr(trueE);
	        len = LengthCMSA(block,cma);
                PosSiteCMSA(block, i, pos, cma);
		s=pos[1]; true_site = gss->TruePos(i,s);
		for(j=0; j < left_flank; j++){
		   k=true_site-left_flank+j;
		   if(k<1){ if(i != bestSq) c='.'; else c='X'; } 
		   else c=AlphaChar(ptrT[k],A);
		   D->sequence[sq][j]=c;
		}
        	for(l=1; l <= len; l++,j++) {
                   s=pos[1]+l-1;
		   // NEW: don't allow deletions in first sequence...
		   if(IsDeletedCMSA(i,s,cma)){ if(i != bestSq) c='.'; else c='X'; } 
		   else c=AlphaChar(ptrF[s],A);
		   D->sequence[sq][j]=c;
		}
                s=pos[1]+l-2; true_site = gss->TruePos(i,s);
                for(t=0; t < right_flank; t++,j++){
                   k=true_site+t+1;
                   if(k>LenSeq(trueE)){ if(i != bestSq) c='.'; else c='X';} 
		   else c=AlphaChar(ptrT[k],A);
		   D->sequence[sq][j]=c;
                }
	} bestSq--;
	// std::cerr << D->sequence[bestSq]; std::cerr << "\nbest prob = "; std::cerr << best_prob; std::cerr << "; bestSq = ";
	// std::cerr << bestSq+1; std::cerr << "\n\n";
	if(bestSq){ tmp=D->sequence[0]; D->sequence[0]=D->sequence[bestSq]; D->sequence[bestSq]=tmp; }
	D->sequence_length = edit_data(read_length+left_flank+right_flank,D);
	return D;
}

dsc_typ MkDSCWithConsensus(Int4 block, Int4 left_flank, Int4 right_flank, cma_typ cma)
{
	dsc_typ 	D;
	Int4		i,j,k,s,t,len,l,pos[3],read_length=0,true_site;
	e_type 		fakeE,trueE;
	unsigned char 	*ptrF, *ptrT;
	gss_typ 	*gss = gssCMSA(cma);
	char		r,R;
	Int4		*number,best;
	a_type		A=AlphabetCMSA(cma);

	// NEW: Use concensus sequence as first.
	NEW(D,1,dsc_type);
	// char    *ConsensusSeqCMSA(Int4 t, cma_typ cma);
	// D->hom_len = NumSeqsCMSA(cma);		// OLD

	D->hom_len = NumSeqsCMSA(cma)+1;	// add one for concensus sequence.
	D->NumHomSeq = D->hom_len;
	D->error_line=0; D->filter_level = 1; D->clean = 1;
	D->clean_length = 40; D->clean_percent = 20;
	NEW(D->att_array,Max_length+2,dsc_att_vector);
	D->names=NULL;
	NEW(number, nAlpha(A) + 3, Int4);
        NEWP(D->sequence,NumSeqsCMSA(cma)+2,char);
        for(i=0;i <= NumSeqsCMSA(cma);i++) NEW(D->sequence[i],Max_length,char);
	read_length = LengthCMSA(block,cma);
	for(j=0; j < left_flank; j++){
	   for(r = 0; r <= nAlpha(A); r++) number[r]=0;
	   for(i=1; i <= NumSeqsCMSA(cma); i++){
		ptrF = SeqPtr(gss->FakeSeq(i)); ptrT = SeqPtr(gss->TrueSeq(i));
		PosSiteCMSA(block, i, pos, cma);
		true_site = gss->TruePos(i,pos[1]);
		k=true_site-left_flank+j;
		if(k<1){ D->sequence[i][j]='.'; r=0; }
		else { r = ptrT[k]; D->sequence[i][j]=AlphaChar(r,A); }
		number[r]++;
	   }
	   for(best=-1,r=1; r <= nAlpha(A); r++) if(best < number[r]) { best=number[r]; R = r; }
	   // std::cerr << AlphaChar(R,A);
	   D->sequence[0][j] = AlphaChar(R,A);
	}
	char *tmp_str = ConsensusSeqCMSA(block, cma);
        for(l=0; l < LengthCMSA(block,cma); l++,j++) {
	   // for(r = 0; r <= nAlpha(A); r++) number[r]=0;
	   for(i=1; i <= NumSeqsCMSA(cma); i++){ 
		ptrF = SeqPtr(gss->FakeSeq(i)); ptrT = SeqPtr(gss->TrueSeq(i));
		PosSiteCMSA(block, i, pos, cma);
		true_site = gss->TruePos(i,pos[1]);
                s=pos[1]+l;
		// NEW: don't allow deletions in first sequence...
		if(IsDeletedCMSA(i,s,cma) && i > 1){ D->sequence[i][j]='.'; r=0; }
		else { r = ptrF[s]; D->sequence[i][j]=AlphaChar(r,A); }
		number[r]++;
	   }
	   // std::cerr << tmp_str[l]; 
	   //D->sequence[0][j] = tmp_str[l];
	   for(best=-1,r=1; r <= nAlpha(A); r++) if(best < number[r]) { best=number[r]; R = r; }
	   // std::cerr << AlphaChar(R,A);
	   D->sequence[0][j] = AlphaChar(R,A);
	} free(tmp_str);
        for(t=0; t < right_flank; t++,j++){
	   for(r = 0; r <= nAlpha(A); r++) number[r]=0;
	   for(i=1; i <= NumSeqsCMSA(cma); i++){ 
		ptrF = SeqPtr(gss->FakeSeq(i)); ptrT = SeqPtr(gss->TrueSeq(i));
		PosSiteCMSA(block, i, pos, cma);
		true_site = gss->TruePos(i,pos[1]);
	        s=pos[1]+l-2;
		true_site = gss->TruePos(i,s);
               	k=true_site+t+1;
               	if(k>LenSeq(gss->TrueSeq(i))){ D->sequence[i][j]='.'; r=0; }
                else { r = ptrT[k]; D->sequence[i][j]=AlphaChar(r,A); }
		number[r]++;
	   }
	   for(best=-1,r=1; r <= nAlpha(A); r++) if(best < number[r]) { best=number[r]; R = r; }
	   // std::cerr << AlphaChar(R,A);
	   D->sequence[0][j] = AlphaChar(R,A);
	}
	// std::cerr << "\n";
	free(number);
	D->sequence_length = edit_data(read_length+left_flank+right_flank,D);
	return D;
}

void    NilDSC(dsc_typ D)
{
	int i;

	free(D->att_array);
	for(i=0;i<=D->NumHomSeq;i++) free(D->sequence[i]); free(D->sequence);
	if(D->names!=NULL){
		for(i=0;i<=D->NumHomSeq;i++) free(D->names[i]); free(D->names);
	}
	free(D);
}

double **ComputeDSC(FILE *outf, double *predicted_accuracy, FILE *fp)
{
	double  **prob;
        Int4    i;   

	dsc_typ D = MkDSC(fp);
        dsc_att_vector *att_array = D->att_array;
                                
        NEWP(prob,4,double);
	for(i=0;i<4;i++) NEW(prob[i],D->sequence_length+2,double);
        PredictSeqDSC(D);
        *predicted_accuracy = pred_acc_DSC(D);
        for (i = 0; i < D->sequence_length+2; i++) {         /* write in column format */
                prob[1][i+1] = att_array[i].prob_a;
                prob[2][i+1] = att_array[i].prob_b;
                prob[3][i+1] = att_array[i].prob_c;
        }                       
	if (outf != NULL){
		for(i=0;i<D->sequence_length;i++)
			fprintf(outf, "%c  p1=%3.3f p2=%3.3f p3=%3.3f\n",D->sequence[0][i],
				prob[1][i+1],prob[2][i+1],prob[3][i+1]);
                                
        }                       
	NilDSC(D);

        return prob;
}                       

double **ComputeDSC(FILE *outf, double *predicted_accuracy, Int4 block, 
		Int4 left_f, Int4 right_f, cma_typ cma)
{
        double  	**tprob,**prob;
        Int4    	i,j,s,pos[3],len,max_left,max_right,left_flank,right_flank;
	e_type 		trueE;
	gss_typ  	*gss = gssCMSA(cma);

        trueE = gss->TrueSeq(1);
        len=LengthCMSA(block,cma);
        PosSiteCMSA(block, 1, pos, cma);
        s=pos[1];
        max_left=gss->TruePos(1,s)-1;
	if(max_left < 0) max_left=0;
        if(left_f>max_left) left_flank=max_left;
	else left_flank=left_f;
        s=pos[1]+len-1;
        max_right=LenSeq(trueE)-gss->TruePos(1,s);
        if(right_f>max_right) right_flank=max_right;
	else right_flank=right_f;

	dsc_typ D = MkDSC(block, left_flank, right_flank, cma);
        dsc_att_vector *att_array = D->att_array;

        NEWP(tprob,4,double);
        for(i=0;i<4;i++) NEW(tprob[i],D->sequence_length+2,double);
 
	NEWP(prob,4,double);
        for(i=1;i<4;i++)
                NEW(prob[i],D->sequence_length-left_flank-right_flank+2,double);
// assert(left_flank== 0 && right_flank ==0);

        PredictSeqDSC(D);
        *predicted_accuracy = pred_acc_DSC(D);   
        for (i = 0; i < D->sequence_length+2; i++) {         /* write in column format */
                tprob[1][i] = att_array[i].prob_a;
                tprob[2][i] = att_array[i].prob_b;
                tprob[3][i] = att_array[i].prob_c;
        }

        for(j=1,i=left_flank;i<D->sequence_length-right_flank;j++,i++){
		prob[1][j] = tprob[1][i];
		prob[2][j] = tprob[2][i];
		prob[3][j] = tprob[3][i];
// fprintf(stderr,"prob(%d) = %.2f %.2f %.2f\n",i+1,prob[1][i+1],prob[2][i+1],prob[3][i+1]);
	}
	if (outf != NULL){
		for(i=0;i<D->sequence_length;i++)
			fprintf(outf, "%c  p1=%3.3f p2=%3.3f p3=%3.3f\n",D->sequence[0][i],
				tprob[1][i],tprob[2][i],tprob[3][i]);
	
	}
	for(i=0;i<4;i++) free(tprob[i]); free(tprob);
	NilDSC(D);
        return prob;
}               

/*************** INSERT predict.c here *****************/

/* Predict accuracy of prediction */
double pred_acc_DSC(dsc_typ D)
	{
	dsc_att_vector *att_array = D->att_array;
	
	double diff,mdiff;
	double predicted_accuracy;
	double res_ratio_v, res_ratio_e;
	int i;

	double total_diff = 0.0;
	double pred_const = 0.2302;
	double cdiff = 0.750;
	double crv = 0.94;
	double cre = 0.72;

	for (i = 0; i <D->sequence_length; i++) {
		diff = pcd(att_array[i].prob_a, att_array[i].prob_b, att_array[i].prob_c);
		total_diff = total_diff + diff;		/* add up diffs betweens to 2 probs  */
		}
	mdiff = total_diff / D->sequence_length;

	res_ratio_v = res_ratio('v',D);		/* proportion of residues valine */
	res_ratio_e = res_ratio('e',D);		/* proportion of residues glutamic acid  */

	predicted_accuracy = pred_const + (cdiff * mdiff) + (crv * res_ratio_v) + (cre * res_ratio_e);

	return predicted_accuracy;

	}




/* decide prediction diff of highest and second highest prob*/
double pcd(double tot_ian, double tot_ibn, double tot_icn)
	{
        double diff;

        if ((tot_icn >= tot_ian) && (tot_icn >= tot_ibn)) {
                if (tot_ian >= tot_ibn)
                       diff = tot_icn - tot_ian;       /* order c a b */
                else
                        diff = tot_icn - tot_ibn;       /* order c b a */
                }
        if ((tot_ian > tot_icn) && (tot_ian >= tot_ibn)) {
                if (tot_icn >= tot_ibn)
                        diff = tot_ian - tot_icn;       /* order a c b */
                else
                        diff = tot_ian - tot_ibn;       /* order a b c */
                }
        if ((tot_ibn > tot_ian) && (tot_ibn > tot_icn)) {
                if (tot_ian >= tot_icn)
                        diff = tot_ibn - tot_ian;       /* order b a c */
                else
                        diff = tot_ibn - tot_icn;       /* order b c a */
                }

        return diff;

        }



/* Proportion of residue type in chain */
double res_ratio(char residue_type, dsc_typ D)
        {
        int i,j;
        double ratio;
        double total = 0.0;
        double count = 0.0;
        char res1;
        int res2;

        for (i = 0; i < D->hom_len; i++)        /* each homologous chain */
                for (j = 0; j < D->sequence_length; j++) { /* each position */
                        res1 = D->sequence[i][j];
                        res2 = tolower(res1);

                        if ((res1 != ':') && (res1 != '.'))  {  /* not deletion */
                                total = total + 1.0;
                                if (res2 == residue_type)
                                count = count + 1.0;
                                }
                        }
        ratio = count / total;
        return(ratio);

        }

int PredictSeqDSC(dsc_typ D)
/* predict sequence: non-homologous */
{
        double ratio_a, ratio_b;
        double  res_ratio_h, res_ratio_e, res_ratio_q, res_ratio_d, res_ratio_r;
        int i;

        form_res_atts_DSC(D);
	smooth_all_DSC(D);
        discrim1_DSC(&ratio_a, &ratio_b, D);

        res_ratio_h = res_ratio('h',D);
        res_ratio_e = res_ratio('e',D);
        res_ratio_q = res_ratio('q',D);
        res_ratio_d = res_ratio('d',D);
        res_ratio_r = res_ratio('r',D);

        discrim2_DSC(ratio_a, ratio_b, res_ratio_h, res_ratio_e, res_ratio_q,
res_ratio_d, res_ratio_r, D);

        estimate_probs_DSC(D);

        for (i = 0; i < D->filter_level; i++) filter_DSC(D);
        remove_isolated_DSC(D);
        return 0;
}


/*************** ADD SMOOTHING ROUTINES HERE *******************/

/* Copy info into working array, add boundary coditions */
int copy_info(int flag, int info_type,
double  work_array[Max_length], dsc_typ D)
	{
dsc_att_vector *att_array = D->att_array;
        int i;

	for (i = 0; i < D->sequence_length; i++) {
		if (flag == 0) 	
			switch (info_type) {
			case 0: work_array[i] = att_array[i].infoa; break;
			case 1: work_array[i] = att_array[i].infob; break;
			case 2: work_array[i] = att_array[i].infoc; break;
			case 3: work_array[i] = att_array[i].edge_dist; break;
			case 4: work_array[i] = att_array[i].deletion; break;
                        case 5: work_array[i] = att_array[i].insertion; break;
			case 6: work_array[i] = att_array[i].hydro_a; break;
                        case 7: work_array[i] = att_array[i].hydro_b; break;
                        case 8: work_array[i] = att_array[i].cons_a; break;
			case 9: work_array[i] = att_array[i].cons_b; break; 
		      }
				
		if (flag == 1)
			switch (info_type) {	
			case 0: att_array[i].s_infoa = work_array[i]; break;
			case 1: att_array[i].s_infob = work_array[i]; break;
			case 2: att_array[i].s_infoc = work_array[i]; break;
			case 3: att_array[i].s_edge_dist = work_array[i]; break; 
			case 4: att_array[i].s_deletion = work_array[i]; break;
			case 5: att_array[i].s_insertion = work_array[i]; break;
			case 6: att_array[i].s_hydro_a = work_array[i]; break;
			case 7: att_array[i].s_hydro_b = work_array[i]; break;
			case 8: att_array[i].s_cons_a = work_array[i]; break;
			case 9: att_array[i].s_cons_b = work_array[i]; break; 
			}
		}
	return(0);
	}


/* Calculate median of 3 in sequence */
int median3(double work_array1[Max_length+1], double
work_array2[Max_length+1], dsc_typ D)
	{
	int i;

	work_array2[0] = work_array1[0];
	work_array2[D->sequence_length - 1] = work_array1[D->sequence_length - 1];

	for (i = 1; i <D->sequence_length - 1; i++) {
		work_array2[i] = m3(i,work_array1);		/* place median in array */
		}
		
	return(0);

	}

/* get median of 3 */
double m3(int i, double work_array[Max_length+1])
	{
	double w[Max_median];
	double n1,m,p1,median;

	clear_w(w);
	n1 = work_array[i-1];
	insert(n1,w); 			/* add n1 in sorted position in list */
	m = work_array[i];
	insert(m,w); 		
	p1 = work_array[i+1];
	insert(p1,w); 		

	median = medi(3,w);		/* get median */

	return(median);
	}

/* Calculate median of 5 in sequence */
int median5(double work_array1[Max_length+1], double
work_array2[Max_length+1], dsc_typ D)
	{
	double w[Max_median];
	double n1,n2,m,p1,p2;
	int i;

	work_array2[0] = work_array1[0];
	work_array2[D->sequence_length - 1] = work_array1[D->sequence_length - 1];

	work_array2[1] = m3(1,work_array1);					/* calculate median of 3 at edge */
	work_array2[D->sequence_length - 2] = m3(D->sequence_length - 2,work_array1);

	for (i = 2; i <D->sequence_length - 2; i++) {
		clear_w(w);
		n2 = work_array1[i-2];
		insert(n2,w); 	
		n1 = work_array1[i-1];
		insert(n1,w); 	
		m = work_array1[i];
		insert(m,w); 		
		p1 = work_array1[i+1];
		insert(p1,w); 		
		p2 = work_array1[i+2];
		insert(p2,w); 		

		work_array2[i] = medi(5,w);		/* place median in array */
		}
		
	return(0);

	}



/* Calculate median of 2 in sequence */
int median2(double work_array1[Max_length+1], double
work_array2[Max_length+1], dsc_typ D)
	{
	double w[Max_median];
	int i;

	work_array2[0] = work_array1[0];
	work_array2[D->sequence_length - 1] = work_array1[D->sequence_length];	/* assume 4 smooth
done */

	for (i = 1; i <D->sequence_length - 1; i++) {
		clear_w(w);

		work_array2[i] = m2(i,work_array1);	/* place median in array */
		}
		
	return(0);

	}


/* get median of 2 */
double m2(int i, double work_array[Max_length+1])
	{
	double m,p1,median;

	m = work_array[i];
	p1 = work_array[i+1];

	median  = (m + p1) / 2.0;		

	return(median);
	}


/* Calculate median of 2 in sequence */
int median4(double work_array1[Max_length+1], double
work_array2[Max_length+1], dsc_typ D)
{
	double w[Max_median];
	double n1,n2,m,p1;
	int i;

	work_array2[0] = work_array1[0];
	work_array2[D->sequence_length] = work_array1[D->sequence_length - 1];	/* increase length by
1 */

	work_array2[1] = m2(0,work_array1);		/* calculate median of 2 at edge */
	work_array2[D->sequence_length - 1] = m2(D->sequence_length - 2,work_array1);

	for (i = 2; i <D->sequence_length - 1; i++) {
		clear_w(w);
		n2 = work_array1[i-2];
		insert(n2,w); 	
		n1 = work_array1[i-1];
		insert(n1,w); 	
		m = work_array1[i];
		insert(m,w); 		
		p1 = work_array1[i+1];
		insert(p1,w); 		
		work_array2[i] = medi(4,w);		/* place median in array */
	}
		
	return(0);
}

int clear_w(double w[Max_median])
{
	int l;

	for (l = 0; l < Max_median; l++) 
		w[l] = 100.0;	/* clear w N.B. asumes all values less than 99 */
	return(0);
}

/* Put p1 in ordered position in array w */
int insert(double p1,double w[Max_median])
	{
	double buffer1, buffer2;
	int i = 0;

	while (w[i] < p1) 	/* find possition for insertion */
		i++;
	
	buffer1 = w[i];		/* swap */
	w[i] = p1;

	while (w[i] < 99.0) {		/* move rest up 1 position */
		i++;
		buffer2 = w[i];
		w[i] = buffer1;
		buffer1 = buffer2;
		}
		
	w[i] = buffer1;		/* put last one in */
	return(0);

	}

/* Calculate median of sequence of length m_length in array w */
double medi(int m_length, double w[Max_median])
	{
	int i, c1, c2;
	int even = 1;
	double cm1, cm2, median;

	for (i = 1; i <= m_length; i++) 	/* decide parity */
		if (even == 1) 
			even = 0;
		else
			even = 1;

	if (even == 1) {			/* average two in centre */
		c1 = m_length / 2;
		c2 = c1 - 1;
		cm1 = w[c1];
		cm2 = w[c2];
		median = (cm1 + cm2) / 2.0;
		}
	else {					/* simple meadian */
		c1 = (m_length - 1) / 2;
		median = w[c1];
		}

	return(median);

	}

/* Deal with 0 and n endpoints */
int endpoints(double work_array[Max_length], dsc_typ D) /* ??? Minitab dpesn't always use
endpoints ??? */
	{
	double w[Max_median];
	double e0,e1,e2;

	e1 = work_array[1];
	e2 = work_array[2];
	e0 = (3.0 * e1) - (2.0 * e2);

	clear_w(w);
	insert(e0,w); 			/* add e0 in sorted position in list */
	insert(e1,w);
	insert(e2,w);
	work_array[0] = medi(3,w);		/* place median in array */

	e1 = work_array[D->sequence_length - 2];
	e2 = work_array[D->sequence_length - 3];
	e0 = (3.0 * e1) - (2.0 * e2);

	clear_w(w);
	insert(e0,w); 			/* add e0 in sorted position in list */
	insert(e1,w);
	insert(e2,w);
	work_array[D->sequence_length - 1] = medi(3,w);		/* place median in array */

	return(0);

	}



/* Calculate Hanning smoothing */
int hanning(double work_array1[Max_length + 1], double work_array2[Max_length + 1], dsc_typ D)
	{
	double n1,m,p1;
	int i;

	work_array2[0] = work_array1[0];
	work_array2[D->sequence_length - 1] = work_array1[D->sequence_length - 1];

	for(i = 1; i < D->sequence_length  - 1; i++) {
		n1 = work_array1[i-1];
		m = work_array1[i];	
		p1 = work_array1[i+1];
	 	work_array2[i] = (0.25 * n1) + (0.5 * m) + (0.25 * p1);
		}

	return(0);

	}

/* calculate rough */
int rough(double work_array1[Max_length + 1], double work_array2[Max_length + 1], double
work_array3[Max_length + 1], dsc_typ D)
	{
	int i;
	for(i = 0; i <D->sequence_length; i++) 
		work_array3[i] = work_array1[i] - work_array2[i];	/* rough = data - smooth */

	return(0);
	}


/* add smoothed values and smoothed rough */
int add(double work_array1[Max_length + 1], double work_array2[Max_length + 1], double
work_array3[Max_length + 1], dsc_typ D)
	{
	int i;

	for(i = 0; i < D->sequence_length; i++) 
		work_array3[i] = work_array1[i] + work_array2[i];

	return(0);
	}


/*************** ADD SMOOTHING ROUTINES HERE *******************/
int smooth_all_DSC(dsc_typ D)
/* Smooth all the attributes using 4253EH,twice */
{
dsc_att_vector *att_array = D->att_array;
        double work_array1[Max_length + 1];
        double work_array2[Max_length + 1];
        double work_array3[Max_length + 1];
        int info_type;

        for (info_type = 0; info_type <10; info_type++) {
                copy_info(0,info_type,work_array1,D);
                median4(work_array1,work_array2,D);      
                median2(work_array2,work_array1,D);
                median5(work_array1,work_array2,D);
                median3(work_array2,work_array1,D);
                endpoints(work_array1,D);
                hanning(work_array1,work_array3,D);
                copy_info(0,info_type,work_array2,D);     
                rough(work_array2,work_array3,work_array1,D);    
                median4(work_array1,work_array2,D);             
                median2(work_array2,work_array1,D);
                median5(work_array1,work_array2,D);
                median3(work_array2,work_array1,D);
                endpoints(work_array1,D);
                hanning(work_array1,work_array2,D);
                add(work_array3,work_array2,work_array1,D);     
                copy_info(1,info_type,work_array1,D); 
                };

        return(0);
}
 
/*********************insert attributes.c*******************/

int form_res_atts_DSC(dsc_typ D)
{
	dsc_att_vector *att_array = D->att_array;
	int i;
	double infoa,infob,infoc;
	const double alpha = 3.6;
	const double beta = 2.0;

	for (i = 0; i < D->sequence_length; i++) {
		PredictGorDSC(i,&infoa, &infob, &infoc, D);
		att_array[i].infoa = infoa; att_array[i].infob = infob; att_array[i].infoc = infoc; 		
		att_array[i].edge_dist = calc_edge_dist(i,D);
		att_array[i].deletion = predict_del(i, D);
		att_array[i].insertion = predict_ins(i,D);
		att_array[i].hydro_a = predict_h_moment(i,alpha, D);
		att_array[i].hydro_b = predict_h_moment(i,beta, D);
		att_array[i].cons_a = predict_c_moment(i, alpha,D);
		att_array[i].cons_b = predict_c_moment(i, beta, D);
	} return(0);
}


double calc_edge_dist(int i,dsc_typ D)
	{
	double edge_dist = 5.0;

	if (i <= 4 ) {
		edge_dist = (double) i + 1.0;
		}

	if ((D->sequence_length - i) <= 5 ) {
		edge_dist = (double) D->sequence_length - i;
		}

	return edge_dist;

	}


double predict_del(int pos1, dsc_typ D)
	{
	double deletion = 0.0;
	int i = 0;

	while ((i < D->hom_len) && (deletion == 0.0)) { /* each homologous chain until found*/
		if (D->sequence[i][pos1] == '.') {	/* deletion character */
			deletion = 1.0;
			}
		i++;
		}

	return deletion;

	}

	
double predict_ins(int pos1,dsc_typ D)
	{
	double insertion = 0.0;
	int i = 0;
	int inschar1,inschar2,inschar3;
	int pos2, pos3;

	while ((i < D->hom_len) && (insertion == 0.0)) { /* each homologous chain until found*/

		inschar1 = (int) D->sequence[i][pos1];
		if ((inschar1 >= 65) && (inschar1 <= 90)) {	/* Capital letter */
			insertion = 1.0;
			}
		else {						/* check neighbours */
			pos2 = pos1 - 1;
			pos3 = pos1 + 1;

			if ((pos2 >= 0) && (pos3 <= D->sequence_length - 1)) { /* in bounds */
				inschar2 = (int) D->sequence[i][pos2];
				inschar3 = (int) D->sequence[i][pos3];

				if ((inschar2 >= 65) && (inschar2 <= 90)) {	
					insertion = 0.5;
					}

				if ((inschar3 >= 65) && (inschar3 <= 90)) {	
					insertion = 0.5;
					}
				}
			}
	
		i++;
		}

	return insertion;

	}


double predict_h_moment(int pos, double rot_angle, dsc_typ D)
	{
	char nc;
	int n1,j;
	double pi, pi2, a;
	double sin_no, cos_no;
	double sin_tot, cos_tot;
	double stot, tot;
	double i;

	int missing = 0;
	double total = 0.0;
	double const_length = 7.0;

	pi = acos(-1.0);
	pi2 = pi * 2.0;
	a = pi2 / rot_angle;	/* angles in radians */

	for (j = 0; j < D->hom_len; j++) {
		n1 = pos - 3;

		nc = get_res1(j, pos, D);	/* character at pos */

		if ((nc != ':') && (nc != '.'))	{	/* not a deletion or end position at i position*/
			sin_tot = 0.0;
			cos_tot = 0.0;
			for (i = 1.0; i < const_length + 1.0; i = i + 1.0) {
				get_res_eis(j, n1, i, a, &sin_no, &cos_no, D);	/*
contribution from -3 to +3*/
				sin_tot = sin_tot + sin_no;
				cos_tot = cos_tot + cos_no;
				n1++;
				}
			stot = (sin_tot * sin_tot) + (cos_tot * cos_tot);
			tot = sqrt(stot);
			total = total + tot;
			}
		else 
			missing++;
		}

	total = total / (D->hom_len - missing);

	return total;
	}

	
double predict_c_moment(int pos, double rot_angle, dsc_typ D)
	{
	int n1;
	double pi, pi2, a;
	double sin_no, cos_no;
	double sin_tot = 0.0;
	double cos_tot = 0.0;
	double stot, tot;
	double i, const_length;
	const_length = 7.0;

	n1 = pos - 3;

	pi = acos(-1.0);
	pi2 = pi * 2;
	a = pi2 / rot_angle;	/* angles in radians */

	for (i = 1.0; i < const_length + 1.0; i = i + 1.0) {
		get_res_con(n1, i, a, &sin_no, &cos_no, D);	/* contribution from 0 - 3
*/
		sin_tot = sin_tot + sin_no;
		cos_tot = cos_tot + cos_no;
		n1++;
		}

	stot = (sin_tot * sin_tot) + (cos_tot * cos_tot);

	tot = sqrt(stot);

	return tot;
	}



int get_res_eis(int hom, int pos, double n, double a, double *sin2, 
double *cos2, dsc_typ D)
	{
	char nc;
	double nf, x, s, c;

	nc = get_res1(hom, pos,D);      /*character at pos  */

	nf = get_eisenberg(nc);		/* get Eisenberg homology value */
	x = n * a;

	s = sin(x);
	*sin2 = nf * s;

	c = cos(x);
	*cos2 = nf * c;
	
	return 0;

	}


int get_res_con(int pos, double n, double a, 
	double *sin2, double *cos2, dsc_typ D)
	{
	double nf, x, s, c;

	nf = get_cons(pos, D);		/* get conservation at position  */

	x = n * a;

	s = sin(x);
	*sin2 = nf * s;

	c = cos(x);
	*cos2 = nf * c;
	
	return 0;

	}


char get_res1(int hom, int n1, dsc_typ D)
	{
	char nc1;

	if ((n1 >= 0) && (n1 <= D->sequence_length - 1))  /* in bounds */
		nc1 = D->sequence[hom][n1];
	else  
		nc1 = 'x'; 

	return nc1;
	}

/* calculate entropy (normalised) at position */
double get_cons(int n1, dsc_typ D)
	{
	int sqn;
	char nc1;
	int nc2;
	int no_seqs;
	int a = 0; int c = 0; int d = 0; int e = 0; int f = 0; int g = 0; int h = 0;
	int i = 0; int k = 0; int l = 0; int m = 0; int n = 0; int p = 0; int q = 0;
	int r = 0; int s = 0; int t = 0; int v = 0; int w = 0; int y = 0; int x = 0;
	double entropy_total = 0.0;
	double entropy;

	for (sqn = 0; sqn < D->hom_len; sqn++) {
		if ((n1 < 0) || (n1 > D->sequence_length)) 
			nc1 = 'x';
		else
			nc1 = D->sequence[sqn][n1];

		nc2 = tolower(nc1);

		switch (nc2) {	
		case 'a': a++ ; break;
		case 'c': c++ ; break;
		case 'd': d++ ; break;
		case 'e': e++ ; break;
		case 'f': f++ ; break;
		case 'g': g++ ; break;
		case 'h': h++ ; break;
		case 'i': i++ ; break;
		case 'k': k++ ; break;
		case 'l': l++ ; break;
		case 'm': m++ ; break;
		case 'n': n++ ; break;
		case 'p': p++ ; break;
		case 'q': q++ ; break;
		case 'r': r++ ; break;
		case 's': s++ ; break;
		case 't': t++ ; break;
		case 'v': v++ ; break;
		case 'w': w++ ; break;
		case 'y': y++ ; break;
		default: x++ ; }
		}

	no_seqs = D->hom_len - x;

	entropy = entrop(a, no_seqs); entropy_total = entropy_total + entropy;
	entropy = entrop(c, no_seqs); entropy_total = entropy_total + entropy;
	entropy = entrop(d, no_seqs); entropy_total = entropy_total + entropy;
	entropy = entrop(e, no_seqs); entropy_total = entropy_total + entropy;
	entropy = entrop(f, no_seqs); entropy_total = entropy_total + entropy;
	entropy = entrop(g, no_seqs); entropy_total = entropy_total + entropy;
	entropy = entrop(h, no_seqs); entropy_total = entropy_total + entropy;
	entropy = entrop(i, no_seqs); entropy_total = entropy_total + entropy;
	entropy = entrop(k, no_seqs); entropy_total = entropy_total + entropy;
	entropy = entrop(l, no_seqs); entropy_total = entropy_total + entropy;
	entropy = entrop(m, no_seqs); entropy_total = entropy_total + entropy;
	entropy = entrop(n, no_seqs); entropy_total = entropy_total + entropy;
	entropy = entrop(p, no_seqs); entropy_total = entropy_total + entropy;
	entropy = entrop(q, no_seqs); entropy_total = entropy_total + entropy;
	entropy = entrop(r, no_seqs); entropy_total = entropy_total + entropy;
	entropy = entrop(s, no_seqs); entropy_total = entropy_total + entropy;
	entropy = entrop(t, no_seqs); entropy_total = entropy_total + entropy;
	entropy = entrop(v, no_seqs); entropy_total = entropy_total + entropy;
	entropy = entrop(w, no_seqs); entropy_total = entropy_total + entropy;
	entropy = entrop(y, no_seqs); entropy_total = entropy_total + entropy;
		
	
	return entropy_total;
	}



/* Eisenberg homology nos */
double get_eisenberg(char nc0)
	{
	double rno;
	int nc1;

	nc1 = tolower(nc0);
	
	switch (nc1) {	
 	case 'a': rno = 0.62; break;
	case 'c': rno = 0.29; break;
	case 'd': rno = -0.90; break;
	case 'e': rno = -0.74; break;
	case 'f': rno = 1.2; break;
	case 'g': rno = 0.48; break;
	case 'h': rno = -0.40; break;
	case 'i': rno = 1.4; break;
	case 'k': rno = -1.5; break;
	case 'l': rno = 1.1; break;
	case 'm': rno = 0.64; break;
	case 'n': rno = -0.78; break;
	case 'p': rno = 0.12; break;
	case 'q': rno = -0.85; break;
	case 'r': rno = -2.5; break;
	case 's': rno = -0.18; break;
	case 't': rno = -0.05; break;
	case 'v': rno = 1.1; break;
	case 'w': rno = 0.81; break;
	case 'y': rno = 0.26; break;
	default: rno = 0.0;
	}

	return rno;
	}


double entrop(int res_no, int no_seqs)
	{
	double entrop;
	double p_res;
	double logp;
	double c;
	
	if (res_no == 0) 
		entrop = 0.0;
	else {
		c = log(2.0);
		c = 1.0 / c;	/* const to convert to log2 */	

		p_res = (double) res_no / no_seqs;
		logp = log(p_res);
		logp = c * logp;

		entrop = logp * p_res * -1.0;
		}

	return entrop;
	}


/*****************************insert discrim.c*****************/

/* First level discrimnation using residue attributes */
int discrim1_DSC(double *ratio_a, double *ratio_b, dsc_typ D)
	{
	int i;
	double discr_a, discr_b, discr_c;
	double na, nb,nc;
	char prediction;

	na = nb = nc = 0.0;

	for (i = 0; i <D->sequence_length; i++) {
		discr_a = discrima1(i, D);
		discr_b = discrimb1(i, D);
		discr_c = discrimc1(i, D);
		prediction = pc(discr_a, discr_b, discr_c);
		switch (prediction) {	
			case 'a': na = na + 1.0; break;
			case 'b': nb = nb + 1.0; break;
			case 'c': nc = nc + 1.0; break;
			}
		}

	 *ratio_a = na /D->sequence_length; 	/* predicted proportion of residues a-helix */
	 *ratio_b = nb /D->sequence_length;	/* predicted proportion of residues b-strand */
	return(0);

	}


/* decide prediction class */
char pc(double tot_ian, double tot_ibn, double tot_icn)
{
        char c;

        if ((tot_icn >= tot_ian) && (tot_icn >= tot_ibn)) {
                c = 'c';
              }
        if ((tot_ian > tot_icn) && (tot_ian >= tot_ibn)) {
                c = 'a';
              }
        if ((tot_ibn > tot_ian) && (tot_ibn > tot_icn)) {
                c = 'b';
              }

        return c;

      }



/* Second level discrimnation using sequence attributes */
int discrim2_DSC(double ratio_a, double ratio_b, double res_ratio_h, 
double res_ratio_e, double res_ratio_q, double res_ratio_d, double res_ratio_r, dsc_typ D)
{
	dsc_att_vector *att_array = D->att_array;
	int i;

	for (i = 0; i < D->sequence_length; i++) {
		att_array[i].prob_a = discrima2(i, ratio_a, ratio_b, res_ratio_h,
				res_ratio_e, res_ratio_q, res_ratio_d, res_ratio_r, D);
		att_array[i].prob_b = discrimb2(i, ratio_a, ratio_b, res_ratio_h,
				res_ratio_e, res_ratio_q, res_ratio_d, res_ratio_r, D);
		att_array[i].prob_c = discrimc2(i, ratio_a, ratio_b, res_ratio_h,
				res_ratio_e, res_ratio_q, res_ratio_d, res_ratio_r, D);
	 	att_array[i].prediction = pc(att_array[i].prob_a, att_array[i].prob_b, 
				att_array[i].prob_c);
	} return(0);

}


double discrima1(int i, dsc_typ D)
	{
 dsc_att_vector *att_array = D->att_array;
	double discrim;
	
	double constant = -50.904;
	double cga = -8.429;
	double cgb = -11.077;
	double cgc = -15.417;
	double cedge = 10.529;
	double chydro_a = -1.278;
	double chydro_b = -0.445;
	double cdel = -1.717;
	double cins = -0.961;
	double ccons_a = -0.277;
	double ccons_b = -0.132;
	double c_sga = -12.523;
	double c_sgb = -5.847;
	double c_sgc = -10.187;
	double c_sedge = 3.960;
	double c_shydro_a = 4.334;
	double c_shydro_b = 1.986;
	double c_sdel = 3.350;
	double c_sins = -2.787;
	double c_scons_a = 2.098;
	double c_scons_b = 0.635;
	
	discrim = 	constant + 
			(cga * att_array[i].infoa) + 
			(cgb * att_array[i].infob) + 
			(cgc * att_array[i].infoc) + 
			(cedge * att_array[i].edge_dist) + 
			(cdel * att_array[i].deletion) + 
			(cins * att_array[i].insertion) + 
			(chydro_a * att_array[i].hydro_a) + 
			(chydro_b * att_array[i].hydro_b) + 
			(ccons_a * att_array[i].cons_a) + 
			(ccons_b * att_array[i].cons_b) + 
			(c_sga * att_array[i].s_infoa) + 
			(c_sgb * att_array[i].s_infob) + 
			(c_sgc * att_array[i].s_infoc) + 
			(c_sedge * att_array[i].s_edge_dist) + 
			(c_sdel * att_array[i].s_deletion) + 
			(c_sins * att_array[i].s_insertion) + 
			(c_shydro_a * att_array[i].s_hydro_a) + 
			(c_shydro_b * att_array[i].s_hydro_b) + 
			(c_scons_a * att_array[i].s_cons_a) + 
			(c_scons_b * att_array[i].s_cons_b); 
	return discrim;
}


double discrimb1(int i, dsc_typ D)
	{
dsc_att_vector *att_array = D->att_array;
	double discrim;
	
	double constant = -49.951;
	double cga = -9.294;
	double cgb = -11.094;
	double cgc = -15.898;
	double cedge = 11.469;
	double chydro_a = -1.219;
	double chydro_b = -0.518;
	double cdel = -1.522;
	double cins = -1.302;
	double ccons_a = -0.285;
	double ccons_b = -0.112;
	double c_sga = -12.733;
	double c_sgb = -5.073;
	double c_sgc = -10.309;
	double c_sedge = 2.942;
	double c_shydro_a = 3.495;
	double c_shydro_b = 2.532;
	double c_sdel = 2.709;
	double c_sins = -1.957;
	double c_scons_a = 1.633;
	double c_scons_b = 0.936;
	
	discrim = 	constant + 
			(cga * att_array[i].infoa) + 
			(cgb * att_array[i].infob) + 
			(cgc * att_array[i].infoc) + 
			(cedge * att_array[i].edge_dist) + 
			(cdel * att_array[i].deletion) + 
			(cins * att_array[i].insertion) + 
			(chydro_a * att_array[i].hydro_a) + 
			(chydro_b * att_array[i].hydro_b) + 
			(ccons_a * att_array[i].cons_a) + 
			(ccons_b * att_array[i].cons_b) + 
			(c_sga * att_array[i].s_infoa) + 
			(c_sgb * att_array[i].s_infob) + 
			(c_sgc * att_array[i].s_infoc) + 
			(c_sedge * att_array[i].s_edge_dist) + 
			(c_sdel * att_array[i].s_deletion) + 
			(c_sins * att_array[i].s_insertion) + 
			(c_shydro_a * att_array[i].s_hydro_a) + 
			(c_shydro_b * att_array[i].s_hydro_b) + 
			(c_scons_a * att_array[i].s_cons_a) + 
			(c_scons_b * att_array[i].s_cons_b);
	
	return discrim;
	}


double discrimc1(int i, dsc_typ D)
	{
dsc_att_vector *att_array = D->att_array;

	double discrim;
	
	double constant = -46.430;
	double cga = -8.927;
	double cgb = -10.945;
	double cgc = -14.721;
	double cedge = 10.679;
	double chydro_a = -1.101;
	double chydro_b = -0.474;
	double cdel = -0.349;
	double cins = 0.502;
	double ccons_a = -0.227;
	double ccons_b = -0.129;
	double c_sga = -12.909;
	double c_sgb = -6.137;
	double c_sgc = -10.401;
	double c_sedge = 3.109;
	double c_shydro_a = 3.697;
	double c_shydro_b = 2.272;
	double c_sdel = 2.557;
	double c_sins = -1.765;
	double c_scons_a = 1.762;
	double c_scons_b = 0.691;

	discrim = 	constant + 
			(cga * att_array[i].infoa) + 
			(cgb * att_array[i].infob) + 
			(cgc * att_array[i].infoc) + 
			(cedge * att_array[i].edge_dist) + 
			(cdel * att_array[i].deletion) + 
			(cins * att_array[i].insertion) + 
			(chydro_a * att_array[i].hydro_a) + 
			(chydro_b * att_array[i].hydro_b) + 
			(ccons_a * att_array[i].cons_a) + 
			(ccons_b * att_array[i].cons_b) + 
			(c_sga * att_array[i].s_infoa) + 
			(c_sgb * att_array[i].s_infob) + 
			(c_sgc * att_array[i].s_infoc) + 
			(c_sedge * att_array[i].s_edge_dist) + 
			(c_sdel * att_array[i].s_deletion) + 
			(c_sins * att_array[i].s_insertion) + 
			(c_shydro_a * att_array[i].s_hydro_a) + 
			(c_shydro_b * att_array[i].s_hydro_b) + 
			(c_scons_a * att_array[i].s_cons_a) + 
			(c_scons_b * att_array[i].s_cons_b);
	
	return discrim;
	}


double discrima2(int i, double ratio_a, double ratio_b, double res_ratio_h, 
		double res_ratio_e, double res_ratio_q, double
		res_ratio_d, double res_ratio_r, dsc_typ D)
{
	dsc_att_vector *att_array = D->att_array;
	double discrim;

	double constant = -86.08;
	double cga = -8.92;
	double cgb = -12.56;
	double cgc = -17.30;
	double cedge = 19.18;
	double chydro_a = -1.31;
	double chydro_b = -0.30;
	double cdel = -2.09;
	double cins = -1.69;
	double ccons_a = -0.42;
	double ccons_b = -0.14;
	double c_sga = -14.64;
	double c_sgb = -3.98;
	double c_sgc = -8.69;
	double c_sedge = -6.00;
	double c_shydro_a = 4.22;
	double c_shydro_b = 1.42;
	double c_sdel = 8.45;
	double c_sins = 0.37;
	double c_scons_a = 3.40;
	double c_scons_b = 0.96;
	double cratio_a = 27.16;
	double cratio_b = 105.98;
	double cratio_rh = 146.37;
	double cratio_re = 162.70;
	double cratio_rq = 339.57;
	double cratio_rd = 322.84;
	double cratio_rr = 98.78;

	discrim = 	constant + 
			(cga * att_array[i].infoa) + 
			(cgb * att_array[i].infob) + 
			(cgc * att_array[i].infoc) + 
			(cedge * att_array[i].edge_dist) + 
			(cdel * att_array[i].deletion) + 
			(cins * att_array[i].insertion) + 
			(chydro_a * att_array[i].hydro_a) + 
			(chydro_b * att_array[i].hydro_b) + 
			(ccons_a * att_array[i].cons_a) + 
			(ccons_b * att_array[i].cons_b) + 
			(c_sga * att_array[i].s_infoa) + 
			(c_sgb * att_array[i].s_infob) + 
			(c_sgc * att_array[i].s_infoc) + 
			(c_sedge * att_array[i].s_edge_dist) + 
			(c_sdel * att_array[i].s_deletion) + 
			(c_sins * att_array[i].s_insertion) + 
			(c_shydro_a * att_array[i].s_hydro_a) + 
			(c_shydro_b * att_array[i].s_hydro_b) + 
			(c_scons_a * att_array[i].s_cons_a) + 
			(c_scons_b * att_array[i].s_cons_b) + 
			(cratio_a * ratio_a) + 
			(cratio_b * ratio_b) + 
			(cratio_rh * res_ratio_h) + 
			(cratio_re * res_ratio_e) + 
			(cratio_rq * res_ratio_q) + 
			(cratio_rd * res_ratio_d) + 
			(cratio_rr * res_ratio_r); 
	
	return discrim;
}

double discrimb2(int i, double ratio_a, double ratio_b, double res_ratio_h, 
		double res_ratio_e, double res_ratio_q, double res_ratio_d, 
		double res_ratio_r, dsc_typ D)
{
dsc_att_vector *att_array = D->att_array;
	double discrim;

	double constant = -88.92;
	double cga = -10.10;
	double cgb = -12.66;
	double cgc = -18.11;
	double cedge = 20.40;
	double chydro_a = -1.26;
	double chydro_b = -0.36;
	double cdel = -1.96;
	double cins = -1.94;
	double ccons_a = -0.45;
	double ccons_b = -0.11;
	double c_sga = -14.90;
	double c_sgb = -3.50;
	double c_sgc = -9.00;
	double c_sedge = -7.36;
	double c_shydro_a = 3.39;
	double c_shydro_b = 1.93;
	double c_sdel = 8.08;
	double c_sins = 0.94;
	double c_scons_a = 3.02;
	double c_scons_b = 1.19;
	double cratio_a = 24.97;
	double cratio_b = 117.07;
	double cratio_rh = 141.50;
	double cratio_re = 181.76;
	double cratio_rq = 349.61;
	double cratio_rd = 333.41;
	double cratio_rr = 110.87;
	
	discrim = 	constant + 
			(cga * att_array[i].infoa) + 
			(cgb * att_array[i].infob) + 
			(cgc * att_array[i].infoc) + 
			(cedge * att_array[i].edge_dist) + 
			(cdel * att_array[i].deletion) + 
			(cins * att_array[i].insertion) + 
			(chydro_a * att_array[i].hydro_a) + 
			(chydro_b * att_array[i].hydro_b) + 
			(ccons_a * att_array[i].cons_a) + 
			(ccons_b * att_array[i].cons_b) + 
			(c_sga * att_array[i].s_infoa) + 
			(c_sgb * att_array[i].s_infob) + 
			(c_sgc * att_array[i].s_infoc) + 
			(c_sedge * att_array[i].s_edge_dist) + 
			(c_sdel * att_array[i].s_deletion) + 
			(c_sins * att_array[i].s_insertion) + 
			(c_shydro_a * att_array[i].s_hydro_a) + 
			(c_shydro_b * att_array[i].s_hydro_b) + 
			(c_scons_a * att_array[i].s_cons_a) + 
			(c_scons_b * att_array[i].s_cons_b) + 
			(cratio_a * ratio_a) + 
			(cratio_b * ratio_b) + 
			(cratio_rh * res_ratio_h) + 
			(cratio_re * res_ratio_e) + 
			(cratio_rq * res_ratio_q) + 
			(cratio_rd * res_ratio_d) + 
			(cratio_rr * res_ratio_r); 
	
	return discrim;
}

double discrimc2(int i, double ratio_a, double ratio_b, double res_ratio_h, 
		double res_ratio_e, double res_ratio_q, double res_ratio_d, 
		double res_ratio_r, dsc_typ D)
{
	dsc_att_vector *att_array = D->att_array;
        double discrim;

        double constant = -84.23;
        double cga = -9.54;
        double cgb = -12.47;
        double cgc = -16.75;
        double cedge = 19.63;
        double chydro_a = -1.14;
        double chydro_b = -0.32;
        double cdel = -0.76;
        double cins = -0.20;
        double ccons_a = -0.38;
        double ccons_b = -0.13;
        double c_sga = -15.06;
        double c_sgb = -4.36;
        double c_sgc = -8.94;
        double c_sedge = -7.20;
        double c_shydro_a = 3.58;
        double c_shydro_b = 1.68;
        double c_sdel = 7.87;
        double c_sins = 1.35;
        double c_scons_a = 3.12;
        double c_scons_b = 0.99;
        double cratio_a = 26.67;
        double cratio_b = 112.69;
        double cratio_rh = 147.84;
        double cratio_re = 175.70;
        double cratio_rq = 346.37;
        double cratio_rd = 329.59;
        double cratio_rr = 106.72;

        discrim =       constant +
                        (cga * att_array[i].infoa) +
                        (cgb * att_array[i].infob) +
                        (cgc * att_array[i].infoc) +
                        (cedge * att_array[i].edge_dist) +
                        (cdel * att_array[i].deletion) +
                        (cins * att_array[i].insertion) +
                        (chydro_a * att_array[i].hydro_a) +
                        (chydro_b * att_array[i].hydro_b) +
                        (ccons_a * att_array[i].cons_a) +
                        (ccons_b * att_array[i].cons_b) +
                        (c_sga * att_array[i].s_infoa) +
                        (c_sgb * att_array[i].s_infob) +
                        (c_sgc * att_array[i].s_infoc) +
                        (c_sedge * att_array[i].s_edge_dist) +
                        (c_sdel * att_array[i].s_deletion) +
                        (c_sins * att_array[i].s_insertion) +
                        (c_shydro_a * att_array[i].s_hydro_a) +
                        (c_shydro_b * att_array[i].s_hydro_b) +
                        (c_scons_a * att_array[i].s_cons_a) +
                        (c_scons_b * att_array[i].s_cons_b) +
                        (cratio_a * ratio_a) +
                        (cratio_b * ratio_b) +
                        (cratio_rh * res_ratio_h) +
                        (cratio_re * res_ratio_e) +
                        (cratio_rq * res_ratio_q) +
                        (cratio_rd * res_ratio_d) +
                        (cratio_rr * res_ratio_r);

        return discrim;
}

/**********************insert estimate_probs.c********************/
/* Convert discrimination score to probability of class */
int estimate_probs_DSC(dsc_typ D)
{
double  probs [96][3] = {
   {0.443,0.344,0.213},
   {0.416,0.348,0.236},
   {0.510,0.314,0.176},
   {0.524,0.381,0.095},
   {0.603,0.167,0.230},
   {0.542,0.301,0.157},
   {0.734,0.233,0.033},
   {0.667,0.267,0.067},
   {0.686,0.149,0.165},
   {0.851,0.085,0.064},
   {0.857,0.048,0.095},
   {0.980,0.010,0.010},
   {0.910,0.036,0.054},
   {0.914,0.052,0.034},
   {0.980,0.010,0.010},
   {0.980,0.010,0.010},
   {0.422,0.216,0.362},
   {0.519,0.156,0.325},
   {0.484,0.118,0.398},
   {0.517,0.046,0.437},
   {0.621,0.113,0.266},
   {0.662,0.080,0.259},
   {0.626,0.096,0.278},
   {0.589,0.036,0.375},
   {0.745,0.119,0.136},
   {0.790,0.067,0.143},
   {0.751,0.065,0.184},
   {0.723,0.020,0.257},
   {0.893,0.032,0.075},
   {0.945,0.014,0.041},
   {0.905,0.015,0.080},
   {0.877,0.004,0.119},
   {0.292,0.491,0.217},
   {0.323,0.516,0.161},
   {0.259,0.611,0.130},
   {0.404,0.586,0.010},
   {0.238,0.561,0.201},
   {0.276,0.512,0.212},
   {0.309,0.681,0.010},
   {0.066,0.924,0.010},
   {0.126,0.680,0.194},
   {0.098,0.771,0.131},
   {0.097,0.806,0.097},
   {0.010,0.980,0.010},
   {0.046,0.867,0.087},
   {0.055,0.931,0.014},
   {0.029,0.961,0.010},
   {0.010,0.933,0.057},
   {0.198,0.418,0.384},
   {0.092,0.511,0.397},
   {0.070,0.466,0.464},
   {0.056,0.476,0.468}, 
   {0.163,0.589,0.248},
   {0.102,0.624,0.274},
   {0.087,0.608,0.307},
   {0.045,0.613,0.342},
   {0.144,0.648,0.208},
   {0.094,0.650,0.256},
   {0.026,0.696,0.278},
   {0.026,0.661,0.313},
   {0.027,0.888,0.085},
   {0.025,0.888,0.087},
   {0.013,0.888,0.099},
   {0.020,0.831,0.149},
   {0.302,0.299,0.399},
   {0.393,0.150,0.457},
   {0.319,0.128,0.553},  
   {0.366,0.042,0.592},
   {0.249,0.234,0.517},
   {0.239,0.180,0.581},
   {0.210,0.090,0.700},
   {0.255,0.051,0.694},
   {0.189,0.134,0.677},
   {0.189,0.138,0.673},
   {0.207,0.122,0.671},
   {0.185,0.021,0.794},
   {0.110,0.069,0.821},
   {0.070,0.040,0.890},
   {0.091,0.027,0.882},
   {0.087,0.020,0.893},
   {0.275,0.330,0.395}, 
   {0.130,0.429,0.441},
   {0.181,0.348,0.471},
   {0.067,0.365,0.568},
   {0.194,0.267,0.539},
   {0.145,0.192,0.665},
   {0.119,0.260,0.621},
   {0.047,0.240,0.713},
   {0.091,0.215,0.694},
   {0.101,0.175,0.725},
   {0.086,0.160,0.754},
   {0.059,0.209,0.732},
   {0.088,0.085,0.827},
   {0.062,0.088,0.850},
   {0.062,0.121,0.817},
   {0.032,0.100,0.868},};

	dsc_att_vector *att_array = D->att_array;
	int i,order,d1,d2,position;
	double diff1,diff2;
	for (i = 0; i < D->sequence_length; i++) { 
					 /* array filled with actual values for discrim not probs */
		if ((att_array[i].prob_a >= att_array[i].prob_b) 
				&& (att_array[i].prob_b > att_array[i].prob_c)) {	
						/* on equality choose c then a */
			order = 0;		/* ordering of predictions, e.g. here abc */
			diff1 = att_array[i].prob_a - att_array[i].prob_b; /* difference between most prob and second most prob */
			diff2 = att_array[i].prob_b - att_array[i].prob_c; /* difference between second most prob and least prob */
			}

		if ((att_array[i].prob_a > att_array[i].prob_c) 
				&& (att_array[i].prob_c >= att_array[i].prob_b)) {
			order = 1;
			diff1 = att_array[i].prob_a - att_array[i].prob_c;     
			diff2 = att_array[i].prob_c - att_array[i].prob_b; 
		}

		if ((att_array[i].prob_b > att_array[i].prob_a) 
				&& (att_array[i].prob_a > att_array[i].prob_c)) {
			order = 2;
			diff1 = att_array[i].prob_b - att_array[i].prob_a;     
			diff2 = att_array[i].prob_a - att_array[i].prob_c; 
		}

		if ((att_array[i].prob_b > att_array[i].prob_c) 
				&& (att_array[i].prob_c >= att_array[i].prob_a)) {
			order = 3;
			diff1 = att_array[i].prob_b - att_array[i].prob_c;     
			diff2 = att_array[i].prob_c - att_array[i].prob_a; 
		}

		if ((att_array[i].prob_c >= att_array[i].prob_a) 
				&& (att_array[i].prob_a >= att_array[i].prob_b)) {
			order = 4;
			diff1 = att_array[i].prob_c - att_array[i].prob_a;     
			diff2 = att_array[i].prob_a - att_array[i].prob_b; 
		}

		if ((att_array[i].prob_c >= att_array[i].prob_b) 
				&& (att_array[i].prob_b > att_array[i].prob_a)) {
			order = 5;
			diff1 = att_array[i].prob_c - att_array[i].prob_b;     
			diff2 = att_array[i].prob_b - att_array[i].prob_a; 
		}
		d1 = get_discrete(diff1);
		d2 = get_discrete(diff2);
		position = (16 * order) + (4 * d1) + d2;
		att_array[i].prob_a = probs[position][0];	/* put empirical prob in array */
		att_array[i].prob_b = probs[position][1];
		att_array[i].prob_c = probs[position][2];
	} return(0);
}

int get_discrete(double diff)
{
	if ((diff >= 0.0) && (diff < 0.5)) return 0;
	else if ((diff >= 0.5) && (diff < 1.0)) return 1;
	else if ((diff >= 1.0) && (diff < 1.5)) return 2;
	else if (diff >= 1.5) return 3;
	else return 1000;	/* error */
}

/****************insert filter.c*****************/

/* Filter predictions */
int filter_DSC(dsc_typ D)
{
	dsc_att_vector *att_array = D->att_array;
	char work[Max_length];
	int i;

	for (i = 0; i< D->sequence_length; i++) 
		if (att_array[i].prediction == 'c') 	// no change if predict c.
			work[i] = 'c';
		else if ((att_array[i].prediction == 'b') && (1 == f1b(i, D))) {
			work[i] = 'c';			// conditions to change b to c.
			att_array[i].prob_c = 0.5;	// exchange probs to make sensible.
			att_array[i].prob_b = 0.4;
			att_array[i].prob_a = 0.1;
		} else if ((att_array[i].prediction == 'a') 
			&& ((1 == f1a(i, D)) || (1 == f2a(i, D)))) { 
			work[i] = 'b'; 			// conditions to change a to b.
			att_array[i].prob_b = 0.5;	// exchange probs to make sensible.
			att_array[i].prob_a = 0.4;
			att_array[i].prob_c = 0.1;
		} else if ((att_array[i].prediction == 'a') 
				&&  ((1 == f3a(i, D)) || (1 == f4a(i, D)) || (1 == f5a(i, D)) 
				|| (1 == f6a(i, D)) ||  (1 == f7a(i, D)) || (1 == f8a(i, D)) 
				|| (1 == f9a(i, D)) || (1 == f10a(i, D)))) {
			work[i] = 'c';			// conditions to change a to c.
			att_array[i].prob_c = 0.5;	// exchange probs to make sensible.
			att_array[i].prob_a = 0.4;
			att_array[i].prob_b = 0.1;
		} else work[i] = att_array[i].prediction;
	        for (i = 0; i<D->sequence_length; i++) att_array[i].prediction =  work[i];

	return(0);
}

int f1b(int i, dsc_typ D)
{
	char cn1,cn2,cn3,cp1,cp2,cp3;

	neighbour(i, &cn3, &cn2, &cn1, &cp1, &cp2, &cp3, D);

	if ((cn3 != 'a') && (cn2 != 'a') && (cn3 == 'c') && (cp2 != 'b')) return(1);
	else return(0);

}

int f1a(int i, dsc_typ D)
{
	char cn1,cn2,cn3,cp1,cp2,cp3;

	neighbour(i, &cn3, &cn2, &cn1, &cp1, &cp2, &cp3, D);

	if ((cn3 != 'a') && (cp1 == 'b')) return(1);
	else return(0);

}

int f2a(int i, dsc_typ D)
{
	char cn1,cn2,cn3,cp1,cp2,cp3;

	neighbour(i, &cn3, &cn2, &cn1, &cp1, &cp2, &cp3, D);

	if ((cn3 == 'c') && (cn2 == 'b') && (cn1 == 'b') && (cp1 == 'a') && (cp3 == 'a')) return(1);
	else return(0);

}

int f3a(int i, dsc_typ D)
{
	char cn1,cn2,cn3,cp1,cp2,cp3;

	neighbour(i, &cn3, &cn2, &cn1, &cp1, &cp2, &cp3, D);

	if ((cn3 != 'a') && (cp1 == 'c')) return(1);
	else return(0);
}

int f4a(int i, dsc_typ D)
{
	char cn1,cn2,cn3,cp1,cp2,cp3;

	neighbour(i, &cn3, &cn2, &cn1, &cp1, &cp2, &cp3, D);

	if ((cn3 == 'a') && (cp1 == 'c') && (cp3 != 'c')) return(1);
	else return(0);

}

int f5a(int i, dsc_typ D)
{
	char cn1,cn2,cn3,cp1,cp2,cp3;

	neighbour(i, &cn3, &cn2, &cn1, &cp1, &cp2, &cp3, D);

	if ((cn2 != 'a') && (cn1 != 'a') && (cp1 == 'a') && (cp2 == 'c') && (cp3 != 'a')) return(1);
	else return(0);

}

int f6a(int i, dsc_typ D)
{
	char cn1,cn2,cn3,cp1,cp2,cp3;

	neighbour(i, &cn3, &cn2, &cn1, &cp1, &cp2, &cp3, D);

	if ((cn3 != 'a') && (cn2 == 'c') && (cn1 != 'c') && (cp1 == 'a') && (cp2 == 'c') && (cp3 != 'a')) return(1);
	else return(0);

}

int f7a(int i, dsc_typ D)
{
	char cn1,cn2,cn3,cp1,cp2,cp3;

	neighbour(i, &cn3, &cn2, &cn1, &cp1, &cp2, &cp3, D);

	if ((cn3 != 'a') && (cn2 == 'c') && (cn1 == 'c') && (cp1 == 'a') && (cp2 != 'b') && (cp3 != 'a')) return(1);
	else return(0);

}

int f8a(int i, dsc_typ D)
{
	char cn1,cn2,cn3,cp1,cp2,cp3;

	neighbour(i, &cn3, &cn2, &cn1, &cp1, &cp2, &cp3, D);

	if ((cn3 == 'a') && (cn2 == 'c') && (cp1 == 'a') && (cp2 == 'a') && (cp3 != 'a')) return(1);
	else return(0);

}

int f9a(int i, dsc_typ D)
{
	char cn1,cn2,cn3,cp1,cp2,cp3;

	neighbour(i, &cn3, &cn2, &cn1, &cp1, &cp2, &cp3, D);

	if ((cn2 == 'c') && (cp1 == 'a') && (cp2 == 'b') && (cp3 != 'a')) return(1);
	else return(0);

}

int f10a(int i, dsc_typ D)
{
	char cn1,cn2,cn3,cp1,cp2,cp3;

	neighbour(i, &cn3, &cn2, &cn1, &cp1, &cp2, &cp3, D);

	if ((cn3 == 'c') && (cp1 == 'a') && (cp2 != 'a') && (cp3 == 'a')) return(1);
	else return(0);

}

/* Get neighbourhood residues i-3 to i+3 */
int neighbour(int i, char *cn1, char *cn2, char *cn3, char *cp1, char *cp2,
char *cp3, dsc_typ D)
{
	dsc_att_vector *att_array = D->att_array;
	int j;

	j = i - 3;
	if (j > 0) *cn3 = att_array[j].prediction; else *cn3 = 'c';
	j = i - 2;
        if (j > 0) *cn2 = att_array[j].prediction; else *cn2 = 'c';
	j = i - 1;
        if (j > 0) *cn1 = att_array[j].prediction; else *cn1 = 'c';
	j = i + 3;
	if (j > 0) *cp3 = att_array[j].prediction; else *cp3 = 'c';
	j = i + 2;
        if (j > 0) *cp2 = att_array[j].prediction; else *cp2 = 'c';
	j = i + 1;
        if (j > 0) *cp1 = att_array[j].prediction; else *cp1 = 'c';
	return(0);
}

/* remove singlet predictions to make more pretty */
int remove_isolated_DSC(dsc_typ D)
{
dsc_att_vector *att_array = D->att_array;
	char work[Max_length];
	int i;

	for (i = 0; i< D->sequence_length; i++) 
		work[i] = att_array[i].prediction;	/* copy prediction into working */

	if ((att_array[1].prediction == 'a') && (att_array[2].prediction == 'c')) {	/* clean beginning */
		work[1] = 'c';
		att_array[i].prob_c = 0.5;		/* exchange probs to make sensible*/
		att_array[i].prob_a = 0.4;
		att_array[i].prob_b = 0.1;
		}
	if ((att_array[1].prediction == 'b') && (att_array[2].prediction == 'c')) {
		work[1] = 'c';
		att_array[i].prob_c = 0.5;		/* exchange probs to make sensible*/
		att_array[i].prob_b = 0.4;
		att_array[i].prob_a = 0.1;
		}

	if ((att_array[D->sequence_length - 1].prediction == 'c') &&
(att_array[D->sequence_length].prediction == 'a')) { /* clean end  */
		work[D->sequence_length] = 'c';
		att_array[i].prob_c = 0.5;		/* exchange probs to make sensible*/
		att_array[i].prob_a = 0.4;
		att_array[i].prob_b = 0.1;
		}
	if ((att_array[D->sequence_length - 1].prediction == 'c') &&
(att_array[D->sequence_length].prediction == 'b')) {
		work[D->sequence_length] = 'c';
		att_array[i].prob_c = 0.5;		/* exchange probs to make sensible*/
		att_array[i].prob_b = 0.4;
		att_array[i].prob_a = 0.1;
		}

	for (i = 1; i< D->sequence_length - 2; i++) {	/* remove isolated */
		if ((att_array[i-1].prediction == 'c') && (att_array[i].prediction == 'a') && (att_array[i+1].prediction == 'c')) {
			work[i] = 'c';
			att_array[i].prob_c = 0.5;		/* exchange probs to make sensible*/
			att_array[i].prob_a = 0.4;
			att_array[i].prob_b = 0.1;
			}
		if ((att_array[i-1].prediction == 'c') && (att_array[i].prediction == 'b') && (att_array[i+1].prediction == 'c')) {
			work[i] = 'c';
			att_array[i].prob_c = 0.5;		/* exchange probs to make sensible*/
			att_array[i].prob_b = 0.4;
			att_array[i].prob_a = 0.1;
			}
		if ((att_array[i-1].prediction == 'a') && (att_array[i].prediction == 'b') && (att_array[i+1].prediction == 'a')) {
			work[i] = 'a';
			att_array[i].prob_c = 0.5;		/* exchange probs to make sensible*/
			att_array[i].prob_b = 0.4;
			att_array[i].prob_a = 0.1;
			}
		if ((att_array[i-1].prediction == 'b') && (att_array[i].prediction == 'a') && (att_array[i+1].prediction == 'b')) {
			work[i] = 'b';
			att_array[i].prob_c = 0.5;		/* exchange probs to make sensible*/
			att_array[i].prob_a = 0.4;
			att_array[i].prob_b = 0.1;
			}

		if ((att_array[i-1].prediction == 'c') && (att_array[i].prediction == 'a') && (att_array[i+1].prediction == 'b')) {
			work[i] = 'c';
			att_array[i].prob_c = 0.5;		/* exchange probs to make sensible*/
			att_array[i].prob_a = 0.4;
			att_array[i].prob_b = 0.1;
			}
		if ((att_array[i-1].prediction == 'c') && (att_array[i].prediction == 'b') && (att_array[i+1].prediction == 'a')) {
			work[i] = 'c';
			att_array[i].prob_c = 0.5;		/* exchange probs to make sensible*/
			att_array[i].prob_b = 0.4;
			att_array[i].prob_a = 0.1;
			}
		if ((att_array[i-1].prediction == 'a') && (att_array[i].prediction == 'b') && (att_array[i+1].prediction == 'c')) {
			work[i] = 'c';
			att_array[i].prob_c = 0.5;		/* exchange probs to make sensible*/
			att_array[i].prob_b = 0.4;
			att_array[i].prob_a = 0.1;
			}
		if ((att_array[i-1].prediction == 'b') && (att_array[i].prediction == 'a') && (att_array[i+1].prediction == 'c')) {
			work[i] = 'c';
			att_array[i].prob_c = 0.5;		/* exchange probs to make sensible*/
			att_array[i].prob_a = 0.4;
			att_array[i].prob_b = 0.1;
			}
		}

	for (i = 0; i< D->sequence_length; i++) att_array[i].prediction =  work[i];
	return(0);
}


/********************insert gor.c**********************/

int ia[21][17] = {       /* gor info for a-helices */
{8, 13, 18, 21, 27, 33, 42, 47, 60, 60, 50, 44, 37, 31, 24, 18, 14},
{-25, -29, -32, -32, -28, -30, -33, -37, -38, -33, -31, -34, -44, -40,
-36, -37, -41},
{13, 14, 15, 14, 14, 15, 10, 10, -17, -17, -31, -37, -23, -19, -18, -13,
-8},
{17, 20, 21, 27, 32, 32, 37, 51, 62, 54, 29, 13, 17, 14, 9, 13, 15},
{-3, -5, -3, -4, -4, -3, -3, -10, 0, 3, 4, 5, -8, -7, 0, -1, -6},
{-3, -8, -12, -20, -33, -45, -64, -81, -101, -67, -40, -28, -20, -19,
-15, -13, -13},
{-9, -12, -8, -9, -13, -17, -26, -22, -19, -10, -9, -7, -4, -4, -5, -5,
0},
{-1, -4, -6, -5, -2, 0, 0, -6, 1, 1, 3, 9, -5, -2, 6, 7, 8},
{-10, -10, -10, -8, -9, -9, -2, 10, 24, 24, 26, 30, 36, 32, 23, 22, 22},
{4, 5, 7, 10, 13, 19, 30, 31, 45, 44, 48, 50, 39, 34, 35, 25, 17},
{13, 12, 15, 16, 17, 22, 32, 30, 39, 47, 52, 49, 37, 26, 27, 17, 10},
{-2, -2, -3, -6, -11, -11, -18, -20, -36, -20, -17, -16, -6, -6, -8, -6,
-4},
{-7, -8, -12, -19, -23, -28, -33, -53, -87, -190, -140, -99, -56, -37,
-32, -21, -10},
{7, 10, 14, 16, 21, 24, 28, 35, 46, 48, 44, 31, 23, 21, 13, 12, 12},
{-2, -2, 0, 3, 4, 10, 14, 23, 32, 33, 33, 35, 32, 25, 20, 18, 17},
{-3, -1, -4, -7, -10, -14, -19, -15, -33, -30, -39, -41, -34, -31, -33,
-30, -25},
{-2, -3, -4, -4, -6, -13, -21, -27, -45, -47, -43, -45, -39, -35, -32,
-28, -23},
{-11, -12, -11, -14, -12, -15, -18, -25, -20, -24, -20, -18, -27, -22,
-14, -9, -7},
{-2, -6, -10, -15, -10, -5, 1, 0, 7, 0, -7, -12, -16, -13, -10, -12,
-20},
{-10, -11, -13, -13, -11, -11, -13, -16, -2, -2, -2, -3, -8, -5, -2, -6,
-10},
{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        };


int ib[21][17] = {       /* gor info for b-strands  */
{0, -5, -12, -19, -24, -32, -33, -33, -34, -29, -28, -23, -22, -25, -21,
-22, -19},
{-1, 6, -1, -5, -9, -11, 0, 18, 44, 45, 27, 18, 9, 6, 11, 2, 2},
{0, 6, 9, 6, 2, -9, -35, -65, -77, -46, -13, 17, 26, 25, 21, 13, 8}, 
{-12, -14, -12, -14, -22, -30, -46, -47, -46, -45, -36, -21, -14, -7,
-6, -8, -7},
{-15, -21, -26, -27, -22, -12, 17, 36, 45, 37, 19, 1, -14, -20, -23,
-15, -12},
{8, 23, 34, 39, 42, 39, 26, -19, -50, -32, 3, 24, 42, 42, 39, 31, 27},
{-2, 9, 16, 13, 11, 11, 4, 2, 3, -4, 4, 9, 9, 7, 6, 7, 4},
{-13, -19, -22, -19, -7, 11, 38, 62, 77, 63, 33, 4, -19, -29, -28, -21,
-14},
{18, 18, 20, 20, 14, 6, -8, -13, -26, -36, -32, -22, -15, -11, -9, -6,
-6},
{-10, -19, -27, -33, -34, -28, -8, 8, 17, 15, -5, -28, -41, -42, -38,
-34, -28},
{-21, -21, -25, -29, -23, -13, -4, 9, 7, 9, -5, -24, -30, -31, -26, -21,
-10},
{9, 14, 20, 17, 14, 6, -17, -50, -58, -23, -1, 14, 25, 25, 19, 15, 9},
{7, 14, 16, 23, 18, -1, -27, -68, -108, -74, -23, 2, 12, 20, 23, 20,
15},
{3, 2, -7, -6, -15, -20, -29, -30, -31, -43, -40, -27, -21, -25, -20,
-11, 1},
{-6, -7, -2, -3, -1, -1, -4, -6, -15, -26, -29, -29, -25, -25, -22, -17,
-18},
{14, 14, 15, 16, 12, 2, -10, -17, -16, -2, 15, 24, 30, 31, 28, 27, 22},
{6, 4, 9, 13, 14, 18, 19, 27, 31, 27, 22, 20, 17, 23, 26, 24, 24},
{-4, -11, -15, -9, 8, 32, 56, 79, 92, 77, 48, 21, -1, -5, -10, -2, -1},
{-6, -8, -25, -27, -20, 0, 15, 39, 42, 26, 16, 5, -3, -10, -4, -4, 1},
{-1, -8, -9, -6, -5, 5, 16, 34, 43, 32, 15, -2, -6, -4, -6, -9, -8},
{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
};

int ic[21][17] = {       /* gor info for coil  */
{-8, -9, -9, -8, -11, -12, -21, -26, -40, -42, -32, -29, -22, -14, -9,
-4, -1},
{23, 21, 28, 31, 30, 33, 29, 19, -1, -6, 7, 17, 31, 29, 23, 30, 33},
{-13, -18, -21, -18, -14, -8, 13, 27, 56, 42, 35, 20, 1, -1, 0, 2, 2},
{-8, -10, -12, -16, -17, -12, -8, -22, -35, -26, -6, 2, -8, -9, -4, -7,
-9},
{13, 18, 19, 21, 18, 11, -10, -18, -36, -32, -18, -6, 16, 19, 15, 11,
13},
{-2, -9, -15, -12, -3, 9, 33, 76, 107, 74, 31, 7, -14, -16, -17, -12,
-9},
{10, 4, -4, -2, 4, 7, 20, 18, 15, 11, 5, 0, -3, -1, 1, -1, -3},
{10, 16, 20, 16, 7, -8, -30, -46, -70, -54, -28, -11, 17, 20, 12, 7, 2},
{-4, -4, -5, -7, -2, 4, 7, 0, -6, -1, -5, -14, -25, -24, -16, -17, -17},
{3, 8, 11, 11, 9, 0, -24, -37, -59, -57, -45, -32, -13, -8, -11, -3, 2},
{2, 2, 2, 3, -2, -12, -29, -37, -44, -55, -50, -33, -17, -6, -9, -3,
-2},
{-5, -8, -12, -7, -1, 6, 27, 47, 64, 32, 16, 4, -13, -13, -7, -6, -3},
{1, -2, -1, 0, 8, 25, 46, 82, 126, 171, 112, 73, 38, 17, 11, 4, -1},
{-9, -11, -9, -12, -10, -10, -9, -15, -26, -21, -19, -13, -9, -4, 1, -4,
-12},
{6, 6, 1, 0, -3, -9, -11, -18, -21, -15, -14, -15, -15, -8, -5, -6, -4},
{-7, -9, -7, -5, 1, 11, 23, 24, 38, 28, 22, 17, 8, 4, 8, 6, 6},
{-2, -1, -2, -5, -5, -1, 5, 4, 15, 19, 20, 24, 21, 13, 9, 6, 3},
{13, 18, 20, 18, 5, -10, -29, -46, -65, -45, -20, 2, 25, 23, 19, 9, 7},
{6, 10, 24, 31, 22, 4, -13, -31, -40, -20, -6, 7, 16, 19, 12, 14, 17},
{10, 15, 18, 15, 13, 7, 0, -11, -32, -23, -9, 4, 11, 7, 6, 12, 14},
{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
};

/* predict gor secondary structure nos for position */
int PredictGorDSC(int pos1, double *infa, double *infb, double *infc, dsc_typ D)
{
        char residue;
        int i, j, pos2;
        int ian, ibn, icn;
        int tot_ian, tot_ibn, tot_icn;
        int missing;
        double tot_jan = 0.0;
        double tot_jbn = 0.0;
        double tot_jcn = 0.0;

        for (j = -8; j < 9; j++) {      /* each position */
                pos2 = pos1 + j;
                if ((pos2 > -1) && (pos2 < D->sequence_length)) {  /* within range */
                        missing = 0;
                        tot_ian = tot_ibn = tot_icn = 0;
                        for (i = 0; i < D->hom_len; i++) {      /* each homologous chain */
                                residue = D->sequence[i][pos2];    /* get residue type */

                                if ((residue == '.') || (residue == ':')) { /* missing residue, deletion, or end */
                                        missing++;
                                        ian = ibn = icn = 0;
                                        }
                                else 
                                       acces(residue,j,&ian,&ibn,&icn);        /* get gor value */

                                tot_ian = tot_ian + ian;
                                tot_ibn = tot_ibn + ibn;
                                tot_icn = tot_icn + icn;
                        }

                        tot_jan = tot_jan + (tot_ian / (D->hom_len - missing));  /* divide by no. of actual homologs */
                        tot_jbn = tot_jbn + (tot_ibn / (D->hom_len - missing));
                        tot_jcn = tot_jcn + (tot_icn / (D->hom_len - missing));
                }
        }
        *infa = tot_jan / 100.0;
        *infb = tot_jbn / 100.0;
        *infc = tot_jcn / 100.0;
        return 0;
}


/* access GOR info */
int acces(char c, int pos, int *ian, int *ibn, int *icn)
        {

        int residue, seq_pos;

        residue = aa_no(c);
        seq_pos = pos_to_pos(pos);

        if ((residue != 999) && (seq_pos != 999)) {
                *ian = ia[residue][seq_pos];
                *ibn = ib[residue][seq_pos];
                *icn = ic[residue][seq_pos];
                }
        else    {
                *ian = *ibn = *icn = 0;  /* out of range */
                }

        return 0;

        }


/* convert sequence position to array row */
int pos_to_pos(int pos1)
        {
        int pos2;

        switch (pos1) {
        case -8: pos2 = 0; break;
        case -7: pos2 = 1; break;
        case -6: pos2 = 2; break;
        case -5: pos2 = 3; break;
        case -4: pos2 = 4; break;
        case -3: pos2 = 5; break;
        case -2: pos2 = 6; break;
        case -1: pos2 = 7; break;
        case  0: pos2 = 8; break;
        case  1: pos2 = 9; break;
        case  2: pos2 = 10; break;
        case  3: pos2 = 11; break;
        case  4: pos2 = 12; break;
        case  5: pos2 = 13; break;
        case  6: pos2 = 14; break;
        case  7: pos2 = 15; break;
        case  8: pos2 = 16; break;
        default: pos2 = 999;
        }

        return pos2;
        }


/* convert residue type to array vcolumn */
int aa_no(char c)
        {
        int cno, rno;

        cno  = tolower(c);

        switch (cno) {
        case 'a': rno = 0; break;
        case 'c': rno = 1; break;
        case 'd': rno = 2; break;
        case 'e': rno = 3; break;
        case 'f': rno = 4; break;
        case 'g': rno = 5; break;
        case 'h': rno = 6; break;
        case 'i': rno = 7; break;
        case 'k': rno = 8; break;
        case 'l': rno = 9; break;
        case 'm': rno = 10; break;
        case 'n': rno = 11; break;
        case 'p': rno = 12; break;
        case 'q': rno = 13; break;
        case 'r': rno = 14; break;
        case 's': rno = 15; break;
        case 't': rno = 16; break;
        case 'v': rno = 17; break;
        case 'w': rno = 18; break;
        case 'y': rno = 19; break;
        case 'z': rno = 20; break;
        case 'x': rno = 20; break;
        case 'b': rno = 20; break;
        case ':': rno = 20; break;
        case '.': rno = 20; break;
        default: rno = 999;
        }

        return rno;
        }


/************************insert edit_input.c here**********************************/

/* remove insertions relative to  sequence to be predicted and edit ends*/
int edit_data(int read_length, dsc_typ D)
	{
	int seq_length;

	upper_lower(read_length, D);		/* sort out case info */

	edit_begin(read_length, D);		/* make all seqs start the same as predicted */

	e_beginnings(D);	/* distinguish insertions at beginning from standard insertions */

	e_ends(D);		/* distinguish insertions at end from standard insertions */

	edit_probe_inserts(read_length, D);	/* find insertions and mark */

	remove_stars(D);				/* clean up insertions */

	seq_length = get_seq_length(D);

	if (1 == D->clean)
		clean_up(seq_length, D);		/* remove bad domain alignments */

	return(seq_length);
	}
	

/* convert into lower case form */
int upper_lower(int read_length, dsc_typ D)
	{
	int i,j;
	int upper,lower;

	upper = lower = 0;

	for (i = 0; i <D->hom_len; i++)	{	/* go through all hom seqs */
		for (j = 0; j < read_length; j++)	{	/* go through all res */
				if (isupper(D->sequence[i][j]))	/* found upper case char */
					upper = 1;
				if (islower(D->sequence[i][j]))	/* found upper case char */
					lower = 1;
				D->sequence[i][j] = (char) tolower(D->sequence[i][j]);	/* convert to lower case */
				}
				
		}

	return(0);
	}


/* remove insertions relative to  sequence */
int edit_probe_inserts(int read_length, dsc_typ D)
	{

	int i,insert_length;

	i = 0;

	while (':' == D->sequence[0][i]) 	/* find beginning of sequence to be predicted */
		i++;

	while (i <= read_length) {
		if ('.' == D->sequence[0][i]) {	/* insertion in seq to be predicted */
			insert_length = calc_insert_lengths(i, D);	
			edit_seqs(i,insert_length, D);	/* process homologous seqs */
			i = i + insert_length;		/* jump over insertion */
			}
		else
			i++;			
		}
	return(0);
	}
	


int calc_insert_lengths(int insert_position, dsc_typ D)
	{
	int insert_length,i;

	insert_length = 0;
	i = insert_position;

	while('.' == D->sequence[0][i]) {	/* get length of insert and edit*/
		D->sequence[0][i] = '*';	/* star for removal */
		i++;
		insert_length++;
		}

	return(insert_length);
	}


/* Go through homologous seqs and mark insertions */
int edit_seqs(int insert_position, int insert_length, dsc_typ D)
	{
	int i,j;


	for (i = 1; i <D->hom_len; i++) 	/* do all but probe */
		for (j = insert_position; j < insert_position + insert_length; j++) {
			if (('.' != D->sequence[i][j]) && (':' != D->sequence[i][j])) {
				D->sequence[i][insert_position - 1] = (char) toupper(D->sequence[i][insert_position - 1]);
				}
			D->sequence[i][j] = '*';	/* star for removal */
			}
	return(0);
	}



/* Make all sequences start where the sequence to predict starts */
int edit_begin(int read_length, dsc_typ D)
	{
	int i,homology_no,start_length,pos1,pos2;

	for (i = 0; i <D->hom_len; i++)
		D->sequence[i][read_length] = '!';	/* mark max ends of sequences */

	start_length = 0;
	while ('.' == D->sequence[0][start_length])
		start_length++;		/* get length of trailing '.' */


	for (homology_no = 0; homology_no < D->hom_len; homology_no++) {
		pos1 = 0;
		pos2 = start_length;
		do {
			D->sequence[homology_no][pos1] = D->sequence[homology_no][pos2];	/* move all up */
			pos1++;
			pos2++; 
			}
		while ('!' != D->sequence[homology_no][pos1]) ;
		}

	return(0);
	}


/* Removes all residues marked with a star */
int remove_stars(dsc_typ D)
	{
	int i,j,k;

	for (i = 0; i <D->hom_len; i++)	{	/* go through all seqs and remove marked by star */
		j = 0;
		while ('!' != D->sequence[i][j]) {
			if ('*' == D->sequence[i][j]) {	/* star */
				k = j;
				do {
					D->sequence[i][k] = D->sequence[i][k+1];	/* move all up */
					k++; 
					}
				while ('!' != D->sequence[i][k]);
				}
			else
				j++;			/* next position */
			}
		}

	return(0);
	}


/* converts first '.' in sequences to ':' */
int e_beginnings(dsc_typ D)
	{
	int i,j;

	for (i = 1; i <D->hom_len; i++)	{	/* go through all hom seqs and remove beginning deletions */
		j = 0;
		while ('.' == D->sequence[i][j]) {
			D->sequence[i][j] = ':';	/* convert */
			j++;
		}
	}
		
	return(0);
	}


/* converts final '.' in sequences to ':' */
int e_ends(dsc_typ D)
	{
	int i,j;

	for (i = 1; i <D->hom_len; i++)	{	/* go through all hom  seqs and remove end deletions */
		j = 0;
		while ('!' != D->sequence[i][j])	/* last residue */
			j++;

		j--;			/* move back and edit */
		while ('.' == D->sequence[i][j]) {
			D->sequence[i][j] = ':';
			j--;
		}
	}
		

	return(0);
	}


int get_seq_length(dsc_typ D)
{
	int length;

	length = 0;
	while ('!' != D->sequence[0][length]) length++;
	return(length);
}


/* Searches for and masks areas of poor alignment */
int clean_up(int seq_length, dsc_typ D)	
	{
	int temp[Max_length];
	double r1,r2;
	int even,l,half;
	int i,j,k;
	int identity;

	even = 1;
	l = D->clean_length;
	while (l > 0) {
		l--;
		if (1 == even)
			even = 0;
		else 
			even = 1;
		}

	if (0 == even)
		half = ((D->clean_length - 1) / 2);
	else
		half = D->clean_length / 2;

	for (i = 1; i < D->hom_len; i++) {				/* go through all hom seqs */

		for (j = half; j <seq_length - half; j++) { 		/* go through all positions */
			temp[j] = 0;					/* initialise pos positions for centre of masking */
			identity = 0;
			for (k = j - half; k < j + half; k++) {		/* go through all residues - not efficient*/
				if ((D->sequence[0][k] == D->sequence[i][k]) || (':' == D->sequence[i][k]))/* sequence identity */
					identity++;
				}

			r1 = ((double) identity / (double) D->clean_length);
			r2 = ((double) D->clean_percent / 100.0);
			if (r1 < r2) 					/* bad alignment as too few identities */
				temp[j] = 1;			/* set for masking */
									

			}

		for (j = half; j < seq_length - half; j++) { 		/* go through positions and mask those marked */
			if (1 == temp[j])				/* marked position */
				for (k = j - half; k < j + half; k++) 		/* go through all residues - not efficient*/
					D->sequence[i][k] = ':';			/* mask position */
			}
		}

	return(0);
	}


/* read in file of Clustal W format */
int read_clustalw_sequence(FILE *fp, dsc_typ D)
	{
	char first[LINE_LENGTH];

	int test,read_length,block,end;

	test = read_length = block = end = 0;

	test = find_clustalw_start(fp,stdout,D);	/* output to stdout not used */

	if (test == 1)
		return(1);						/* error in format */
	else {
		test = read_clustalw_block(fp,stdout,first,block,&read_length,&end,D);/* output to stdout not used */

		if (0 != test) {
			fprintf(stderr,"ERROR: no sequences found \n");
			return(1);
			}
		else {
			block++;
			while ((0 == test) && (0 == end)) {
				test = read_clustalw_block(fp,stdout,first,block,&read_length,&end, D); /* output to stdout not used */
				block++;
				}

			if (0 != test) 
					return(1);
			else {
				D->sequence_length = edit_data(read_length,D);
				return(0);						/* set to 1 for testing */
				}
			}
		}
	}
		
	


/* get start of file, assumes lines < LINE_LENGTH */
int find_clustalw_start(FILE *fp_in, FILE *fp_out, dsc_typ D)
	{
	char temp[LINE_LENGTH];
	char *pt1,*pt2,*test;
	int found_msf, error;

	found_msf = error = 0;

	while ((found_msf == 0) && (error == 0)) {
		test = fgets(temp,Max_length,fp_in);
		D->error_line++;

		if(test == NULL) {
			error = 1;		/* end of file */
			fprintf(stderr,"ERROR: failed to find start of file marked by CLUSTAL and W");
			}
		else {
			pt1 = strstr(temp,"CLUSTAL");
			pt2 = strstr(temp,"W");
			}
		if ((pt1 != NULL) && (pt2 != NULL))
			found_msf = 1;
		}

	return(0);

	}


/* Reads first block of aligned sequnces in Clustal W format, reads names into array for use in checking other blocks */
int read_clustalw_block(FILE *fp_in, FILE *fp_out, char first[LINE_LENGTH], int block, 
		int *out_sequence_length, int *file_end, dsc_typ D)
	{
	char temp[LINE_LENGTH];
	char name[MAX_NAME_LENGTH];

	char *test;
	int start,block_end,error;
	int old_sequence_length,new_sequence_length;
	int homologous_sequences;


	start = block_end = error = old_sequence_length = homologous_sequences = 0;

	while ((0 == *file_end) && (0 == block_end) && (0 == error)) {

		test = fgets(temp,Max_length,fp_in);
		D->error_line++;

		if (test == NULL)		/* end of file found */
			*file_end = 1;

		new_sequence_length = read_clustalw_line(temp, name, block, homologous_sequences, *out_sequence_length, D);


		if (new_sequence_length < 0)
			error = 1;
		else if (0 == new_sequence_length) {		/* empty line */
			if (1 == start)
				block_end = 1;
			}
		else if (new_sequence_length > 0) {
			if (0 == start) {
				strcpy(first,temp);
				start = 1;
				}
			if (0 == homologous_sequences) {
				old_sequence_length = new_sequence_length;
				}
			else if (old_sequence_length != new_sequence_length) {
				error = 1;
				fprintf(stderr,"ERROR: different string lengths, line %d, expected = %d, actual = %d\n",
					D->error_line, old_sequence_length, new_sequence_length);
				}

			if ((0 == error) && (0 == block)) {
				D->names[homologous_sequences] = (char*) malloc(strlen(name) + 1);	/* copy string */
				strcpy(D->names[homologous_sequences],name);
				}

			if (0 == error) 	
				homologous_sequences++;		/* sequence read ok */
			}

		}

	if ((0 == error) && (0 == block))			/* save no. of homologous sewqunces */
		D->hom_len = homologous_sequences;
	else if ((0 == error) && (block > 0) && (D->hom_len != homologous_sequences) && (1 == start)) { 
		fprintf(stderr,"ERROR: different number of homologous sequences read, line %d, expected = %d, actual = %d\n",
D->error_line, D->hom_len, homologous_sequences);
		}

	if (0 == error)						/* number of residues read */
		*out_sequence_length = *out_sequence_length + old_sequence_length;

	return(error);
	}



/* Read single line of CLUSTAL W line */
int read_clustalw_line(char temp[LINE_LENGTH],char name[MAX_NAME_LENGTH], int block, 
	int homologous_sequences, int out_sequence_length, dsc_typ D)
	{
	char *pt;
	int error,name_found,start;
	int i,j,h;
	int new_sequence_length;

	error = name_found = start = 0;
	i = 0;

	while (0 == name_found) {
		if (0 != isalnum(temp[i]))		/* character found */
			name_found = 1;	
		else if (('\n' == temp[i]) || (i < LINE_LENGTH))	/* end of line */
			name_found = 2;
		else
			i++;
		}


	if (1 == name_found) {		/* possible name found for line */

		j = 0;

		if (i < LINE_LENGTH) { 		/* copy name from temp into name array and check ok */
			do {
				name[j] = temp[i];
				i++;
				j++;
				}
			while ((' ' != temp[i]) && (i < LINE_LENGTH));
			name[j++] = '\0';       /* add end */
			}

		if (block > 0)  {	/* check name has been seen before */
			h = 0;
			while ((h < D->hom_len) && (0 == start)) {
				pt = strstr(name,D->names[h]);	/* check if known name present */
				if (pt != NULL) 	/* name seen before */
					start = 1;
				else
					h++;		/* check next name */
				}
			if (0 == start) {
				error = 1;
				fprintf(stderr,"ERROR: new name in block, line %d, name = %s\n",D->error_line,name);
				}
			}

		if (0 == error) {			
			while (' ' == temp[i])	/* read over blancks */
				i++;

			if ('\n' == temp[i]) {
				error = 1; 	/* could not find sequence */
				fprintf(stderr,"ERROR: could not read CLUSTAL W sequence, line %d\n",D->error_line);
				}
			}

		if (0 == error) {	
			new_sequence_length = read_clustalw_seq(temp,homologous_sequences,out_sequence_length,i,D);/* Read sequence */
			return(new_sequence_length);
			}
		else
			return(-1);	/* error in line */

		}

	else
		return(0);	/* empty line */

	}



/* Read in sequence of residues in clustal w format */
int read_clustalw_seq(char temp[LINE_LENGTH], int seq_no, int out_sequence_length, int i, dsc_typ D)
{
	int end;
	int j,k;

	end = 0;
	j = i;
	k = 0;

	while (0 == end) {
		if (' ' == temp[j])
			j++;
		else if ((0 != isalpha(temp[j])) || ('-' == temp[j])) {
			if ('-' == temp[j])
				D->sequence[seq_no][out_sequence_length + k] = '.';
			else
				D->sequence[seq_no][out_sequence_length + k] = temp[j];
			j++;
			k++;
		}
		else 
			end = 1;
	}

	return(k);
}


