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

#include "sseq.h"

ssq_typ::ssq_typ() {}

ssq_typ::ssq_typ(ptm_typ PM, double pernats, Int4 max_gap_len)
{ make_ssq(PM, pernats, max_gap_len);}

void ssq_typ::make_ssq(ptm_typ PM, double pernats, Int4 max_gap_len)
{
	Int4 i,j,r,s,g,n=0,total=0;
	double sum=0.;	
	a_type	A;

	MaxGap=max_gap_len;
        nmbr_of_blocks=NumModelsPrtnModel(PM);

        NEW(block_lengths, nmbr_of_blocks+1, Int4);
	NEW(start_prof, nmbr_of_blocks+2, Int4);	

	smx_typ *M;
	NEW(M,nmbr_of_blocks+1,smx_typ);

        for(i=1;i<=nmbr_of_blocks;i++){
		start_prof[i]=n+1;
		M[i]=SMatrixPrtnModel(i,PM);
                block_lengths[i] = LenSMatrix(M[i]);
                n+=block_lengths[i];
        }
        prof_len=n; 
        A = SMatrixA(M[1]);
	alph_len=nAlpha(A);

        NEWP(gpen,nmbr_of_blocks+1,double);
        for(i=0;i<=nmbr_of_blocks;i++) NEW(gpen[i],MaxGap+1,double);

        Int4    gapscore;
        for(s=1;s < nmbr_of_blocks;s++){
           for(sum=0.,g=0;g<=MaxGap;g++){
                gapscore=GapScorePrtnModel(s,g,PM);
                if(gapscore <= SHRT_MIN) gpen[s][g]=0;
                else gpen[s][g]=exp((double)gapscore/pernats);
		sum+=gpen[s][g];
           }
	   for(g=0;g<=MaxGap;g++)
		gpen[s][g]/=sum;
        }
	
        NEWP(match,alph_len+1,double);
        for(i=0;i<=alph_len;i++) NEW(match[i],prof_len+1,double);

        for(i=1;i<=nmbr_of_blocks;i++){
           for(s=1;s<=block_lengths[i];s++){
                for(r=1;r<=alph_len;r++){
                   match[r][total+s]=exp((double)ValSMatrix(s,r,M[i])/pernats);		    
		}
           }
           total+=block_lengths[i];
        }
	for(j=1;j<=prof_len;j++){
		for(sum=0.,r=1;r<=alph_len;r++) sum+=match[r][j];
		for(r=1;r<=alph_len;r++) match[r][j]/=sum;
	}
	free(M);
}

e_type ssq_typ::sample_ssq(Int4 I)
{
	Int4	i,j,k,t,g,s,s_len;
	e_type	E;
	unsigned char	*sequence;
	double		rand_num,rand_num1,sum;
	
	s_len=(nmbr_of_blocks-1)*MaxGap+prof_len+30;
	NEW(sequence,s_len,unsigned char);

	static double backfreq[21] = {  0.0, /*  X
            C     G     A     S     T     N     D     E     Q     K */
        0.025,0.074,0.074,0.057,0.051,0.045,0.054,0.054,0.034,0.058,
        /*  R     H     W     Y     F     V     I     L     M     P */
        0.052,0.026,0.013,0.032,0.047,0.073,0.068,0.099,0.025,0.039 };

	for(sum=0.0,k=1;k<=alph_len;k++) sum+=backfreq[k];
	for(i=1;i<=10;i++){
		rand_num=sum*SampleUniformProb();
		for(k=1;k<=alph_len;k++){
			if ((rand_num-=backfreq[k])<=0.0) {sequence[i]=k; break;}
		}
	}
	s=2;
	for(j=1;j<=prof_len;j++){
		    rand_num=SampleUniformProb();
		    for(k=1;k<=alph_len;k++)
			if ((rand_num-=match[k][j])<=0.0) {sequence[i++]=k;break;}
            	    if (j==start_prof[s]-1){
		    	rand_num=SampleUniformProb();
		  	for(g=0;g<=MaxGap;g++){
		            if ((rand_num-=gpen[s-1][g])>0.0){
				rand_num1=SampleUniformProb();
				for(k=1;k<=alph_len;k++)
		        		if ((rand_num1-=backfreq[k])<=0.0) {sequence[i++]=k;break;}
		     	    } else break;
		  	} s++;
		     }
	}
	for(t=0;t<=9;t++){
	   rand_num=SampleUniformProb();
           for(k=1;k<=alph_len;k++)
                 if ((rand_num-=backfreq[k])<=0.0) {sequence[i++]=k; break;}
	}
	E = StringToSeq2(sequence, i-1, "Simulated", I); free(sequence);
	return E;
}

e_type ssq_typ::sample_ssq(Int4 I,Int4 block,Int4 column,Int4 inslen, double *left_relen, double *right_relen)
{
        Int4 i,j,jp,k,r,t,g,s,s_len;
        e_type E;
        unsigned char *sequence;
        double rand_num,rand_num1,sum;
        
        s_len=(nmbr_of_blocks-1)*MaxGap+prof_len+100;
        NEW(sequence,s_len,unsigned char);
        
        static double backfreq[21] = {  0.0, /*  X
            C     G     A     S     T     N     D     E     Q     K */
        0.025,0.074,0.074,0.057,0.051,0.045,0.054,0.054,0.034,0.058,
        /*  R     H     W     Y     F     V     I     L     M     P */
        0.052,0.026,0.013,0.032,0.047,0.073,0.068,0.099,0.025,0.039 };
                
        for(sum=0.0,k=1;k<=alph_len;k++) sum+=backfreq[k];
        for(i=1;i<=10;i++){
                rand_num=sum*SampleUniformProb();
                for(k=1;k<=alph_len;k++){
                        if ((rand_num-=backfreq[k])<=0.0) {sequence[i]=k; break;}
                }
        }
        s=2;
	double temp_relen;
	Int4 position=start_prof[block]+column-1;
	Int4 end_block=start_prof[block]+block_lengths[block]-1;
	if (inslen==0){
                for(temp_relen=0.,j=start_prof[block];j<position;j++){
                        for(r=1;r<=alph_len;r++)
                                temp_relen+=match[r][j]*log(match[r][j]/backfreq[r]);
                }
                *left_relen=temp_relen;
                                        
                for(temp_relen=0.,j=position+1;j<=end_block;j++){
                        for(r=1;r<=alph_len;r++)
                                temp_relen+=match[r][j]*log(match[r][j]/backfreq[r]);
                }
                *right_relen=temp_relen;
                
	        for(j=1;j<position;j++){
       	             	rand_num=SampleUniformProb();
       	             	for(k=1;k<=alph_len;k++)
       	                	if ((rand_num-=match[k][j])<=0.0) {sequence[i++]=k;break;}
       	         	if (j==start_prof[s]-1){
       	           		rand_num=SampleUniformProb();
       	          	 	for(g=0;g<=MaxGap;g++){
       	        	    	  	if ((rand_num-=gpen[s-1][g])>0.0){
       	                			rand_num1=SampleUniformProb();
       	       		        	  	for(k=1;k<=alph_len;k++)
       	                	 			if ((rand_num1-=backfreq[k])<=0.0) {sequence[i++]=k;break;}
					}
					else break;
				}					
       	         		s++;
			}
       		}
		for(j=position;j<prof_len;j++){
                        rand_num=SampleUniformProb();
                        for(k=1;k<=alph_len;k++)
                        	if ((rand_num-=match[k][j+1])<=0.0) {sequence[i++]=k;break;}
                        if (j==start_prof[s]-1){
                                rand_num=SampleUniformProb();
                                for(g=0;g<=MaxGap;g++){
                                        if ((rand_num-=gpen[s-1][g])>0.0){
                                                rand_num1=SampleUniformProb();
                                                for(k=1;k<=alph_len;k++)
                                                        if ((rand_num1-=backfreq[k])<=0.0) {sequence[i++]=k;break;}
                                        } else break;
				}
                        	s++;
			}
		}
	}
	else {
                for(temp_relen=0.,j=start_prof[block];j<=position;j++){
                        for(r=1;r<=alph_len;r++)
                                temp_relen+=match[r][j]*log(match[r][j]/backfreq[r]);
                }
                *left_relen=temp_relen;
                                
                for(temp_relen=0.,j=position+1;j<=end_block;j++){
                        for(r=1;r<=alph_len;r++)
                                temp_relen+=match[r][j]*log(match[r][j]/backfreq[r]);
                }
                *right_relen=temp_relen;
                             
                for(j=1;j<=position;j++){
                        rand_num=SampleUniformProb();
                        for(k=1;k<=alph_len;k++)
                                if ((rand_num-=match[k][j])<=0.0) {sequence[i++]=k;break;}
                        if (j==start_prof[s]-1){
                                rand_num=SampleUniformProb();
                                for(g=0;g<=MaxGap;g++){
                                        if ((rand_num-=gpen[s-1][g])>0.0){
                                                rand_num1=SampleUniformProb();
                                                for(k=1;k<=alph_len;k++)
                                                        if ((rand_num1-=backfreq[k])<=0.0) {sequence[i++]=k;break;}
                                        } else break;
				}
                        	s++;
			}
                }
		for(t=1;t<=inslen;t++){
                	rand_num=SampleUniformProb();
                	for(k=1;k<=alph_len;k++)
                        	if ((rand_num-=backfreq[k])<=0.0) {sequence[i++]=k;break;}
                }
			
                for(j=position+1;j<prof_len;j++){
                        rand_num=SampleUniformProb();
                        for(k=1;k<=alph_len;k++)
                                if ((rand_num-=match[k][j])<=0.0) {sequence[i++]=k;break;}
                        if (j==start_prof[s]-1){
                                rand_num=SampleUniformProb();
                                for(g=0;g<=MaxGap;g++){
                                        if ((rand_num-=gpen[s-1][g])>0.0){
                                                rand_num1=SampleUniformProb();
                                                for(k=1;k<=alph_len;k++)
                                	                        if ((rand_num1-=backfreq[k])<=0.0) {sequence[i++]=k;break;}
                                        } else break;
				}
                        	s++;
			}
                }

	}

        for(t=0;t<=9;t++){
           	rand_num=SampleUniformProb();
           	for(k=1;k<=alph_len;k++)
                 	if ((rand_num-=backfreq[k])<=0.0) {sequence[i++]=k; break;}
        }
//        a_type A = MkAlpha(AMINO_ACIDS,PROT_BLOSUM62);
//for (j=1;j<=tem_len;j++) printf("%c", AlphaChar(tem_sequence[j],A));
//printf("\n");
//for (j=1;j<=sq_len;j++) printf("%c", AlphaChar(sequence[j],A));
	char defline[200];
	sprintf(defline,"Simulated left_info=%.2f; right_info=%.2f",*left_relen,*right_relen);
        E = StringToSeq2(sequence, i-1, defline, I); free(sequence);
        return E;


}

ssq_typ::~ssq_typ() {Free();}

void ssq_typ::Free()
{
	Int4 i;

        for(i=0;i<=alph_len;i++) free(match[i]);
	free(match);

        for(i=0;i<=nmbr_of_blocks;i++)
		free(gpen[i]);
	free(gpen);
	free(block_lengths); free(start_prof);

}
