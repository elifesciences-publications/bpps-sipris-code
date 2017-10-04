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

#include "sdp.h"

char *AlignBlockWithGps(Int4 seq_len, unsigned char *seq, Int4 nblks, smx_typ *M, 
	Int4 **gpen, Int4 *io, Int4 *ie, Int4 *od, Int4 *de, Int4 *alignscore, Int4 *start, 
	Int4 *oper_len)
{
	Int4	pen,i,j,k,l,im1,jm1,i0,n,r,s,t,prof_len,total,g,s1,s0;
	Int4	maxMAT,maxDEL;
	Int4	*start_prof;
	Int4	**MAT,**DEL,**INS,**score;
	a_type	A=SMatrixA(M[1]);

	MEW(start_prof, nblks+2, Int4);	
        for(n=0,i=1;i<=nblks;i++){
		start_prof[i]=n+1; n+=LenSMatrix(M[i]);
        } prof_len=n; start_prof[nblks+1]=0;

        NEWP(MAT,prof_len+1,Int4); NEWP(DEL,prof_len+2,Int4);
	NEWP(INS,prof_len+1,Int4); NEWP(score,prof_len+2,Int4);
	NEW(MAT[0],seq_len+1,Int4); NEW(DEL[0],seq_len+2,Int4);
	NEW(INS[0],seq_len+1,Int4); NEW(score[0],nAlpha(A)+2,Int4);
        for(i=1;i<=prof_len;i++){
                MEW(MAT[i],seq_len+1,Int4); MEW(DEL[i],seq_len+1,Int4);
		MEW(INS[i],seq_len+1,Int4); NEW(score[i],nAlpha(A)+2,Int4);
        }

        for(total=0,i=1;i<=nblks;i++){
           for(s=1;s<=LenSMatrix(M[i]);s++){
                for(r=0;r<=nAlpha(A);r++){
                   score[total+s][r]=ValSMatrix(s,r,M[i]);
		}
           } total+=LenSMatrix(M[i]);
        }
	DEL[0][0]=INS[0][0]=MAT[0][0]=0;
	DEL[1][0]=INS[1][0]=MAT[1][0]=od[1]+de[1];
	// gpen[s-1][0] is gap function for block s-1 with gap length = 0.
        for(s=2,jm1=1,j=2;j<=prof_len;jm1++,j++) {
                if(j!=start_prof[s]){
			MAT[j][0]=INS[j][0]=DEL[j][0]=DEL[jm1][0]+de[j];
                } else {
			MAT[j][0]=DEL[jm1][0]+od[j]+de[j]+gpen[s-1][0]; 
			DEL[j][0]=DEL[jm1][0]+od[j]+de[j]+gpen[s-1][0];
			INS[j][0]=DEL[jm1][0]+od[j]+de[j]+gpen[s-1][0];
			s++; 
		}
        }

	for(i=1;i<=seq_len;i++) {DEL[0][i] = od[1];  MAT[0][i] = 0; INS[0][i] = 0;}

	Int4 J,Jm1;	
        for(s=2,J=1,Jm1=0;J<=prof_len;J++,Jm1++){
	    	register Int4 *matj= MAT[J];
	    	register Int4 *matjm1= MAT[Jm1];
	    	register Int4 *insj= INS[J];
	    	register Int4 *insjm1= INS[Jm1];
	    	register Int4 *delj= DEL[J];
	    	register Int4 *deljm1= DEL[Jm1];
	      	register Int4  delo = od[J];
	    	register Int4  dele = de[J];
	    	register Int4  inso = io[J];
	    	register Int4  inse = ie[J];	  
	    	register Int4  *scorej = score[J];  
	   if(J!=start_prof[s]){
	        for(im1=0,i=1;i<=seq_len;im1++,i++){
        	        r = seq[i];
		    	if(matjm1[im1]>=deljm1[im1] && matjm1[im1]>=insjm1[im1])
				matj[i] = scorej[r]+matjm1[im1];
		    	else if(deljm1[im1]>=insjm1[im1])
				matj[i] = scorej[r]+deljm1[im1];
		    	else matj[i] = scorej[r]+insjm1[im1];

		    	s0 = delo+dele+matjm1[i];
		    	if(s0 > (s1=dele+deljm1[i])) delj[i]=s0; else delj[i]=s1;
		    	s0 = inso+inse+matj[im1];
		    	if(s0 > (s1=inse+insj[im1])) insj[i]=s0; else insj[i]=s1;
		}
	   } else{			
	    	register Int4 *gpen_s = gpen[s-1];
                for(im1=0,i=1;i<=seq_len;im1++,i++){
        	        r = seq[i];
		   	maxMAT=SHRT_MIN;
                   	for(i0=im1,g=0;g < i;g++,i0--){
				if((pen=gpen_s[g]) <= SHRT_MIN) break;
				if(maxMAT < (s1=pen+matjm1[i0])) maxMAT=s1;
				if(maxMAT < (s1=pen+deljm1[i0])) maxMAT=s1;
		   	}
		   	maxDEL=SHRT_MIN;
		   	for(i0=i,g=0;g <= i;g++,i0--){
				if((pen=gpen_s[g]) <= SHRT_MIN) break;
				if(maxDEL < (s1=pen+matjm1[i0])) maxDEL=s1;
				if(maxDEL < (s1=pen+deljm1[i0])) maxDEL=s1;
		   	}
                   	matj[i] = scorej[r]+maxMAT;
                   	delj[i] = delo+dele+maxDEL;
		   	s0 = inso+inse+matj[im1];
		   	if(s0 > (s1=inse+insj[im1])) insj[i]=s0; else insj[i]=s1;
		} s++;
	   }
	}

        Int4    bl_nmbr,best_i;
        Int4    *gapfunct,gapend,i_gap,gapstart;
        Int4	mxm, gap;
        char    state,*operation, *back_operation;

        NEW(operation,seq_len+prof_len+3,char);
        NEW(back_operation,seq_len+prof_len+3,char);
        j=prof_len; i=seq_len;

        for(state='E',k=0;state != 'X'; ){
           switch(state){
                case 'B': // begin state (exit traceback).
			back_operation[k++]='E'; *start=i; state='X'; break;
                case 'E': // end state.  
			back_operation[k++]='E';
			bl_nmbr=nblks; s=start_prof[bl_nmbr];
			mxm=MAT[j][i];state='M';i0=i;
			for(i=1;i<=seq_len;i++){
			   if(MAT[j][i]>mxm){ state='M'; mxm=MAT[j][i];i0=i;}
			   if(DEL[j][i]>mxm){ state='D'; mxm=DEL[j][i];i0=i;}
			}
			i=i0;*alignscore=mxm;
			break;
                case 'M': // previously sampled Matched state.
                        if(j <= 1) { back_operation[k++]='M'; state='B'; break; }
                        if(j==s){ j--; back_operation[k++]='M'; state='G'; gapstart=i-1; }
                        else { // else we are within a block & use affine gaps.
			   j--; i--; back_operation[k++]='m';
                           if(MAT[j][i]>INS[j][i] && MAT[j][i]>DEL[j][i]) state='M';
                           else {if(INS[j][i]>DEL[j][i]) state='I';
                           	 else state='D';
			   }
			}
                        break;
                          
                case 'D': // previously sampled deletion
                          if(j<=1){ back_operation[k++]='D'; state='B'; i++; break; }
                          if(j==s){ j--; back_operation[k++]='D'; state='G'; gapstart=i; }
                          else { j--; back_operation[k++]='d';
                             if(DEL[j][i]+de[j]>=MAT[j][i]+od[j]+de[j]) state='D';
                             else state='M';
                          }
                          break;

                case 'I': // previously sampled insertion.
                          back_operation[k++]='I'; i--; 
                          if (INS[j][i]+ie[j]>MAT[j][i]+io[j]+ie[j]) state='I';
                          else state='M'; 
                          break;
                   					
                case 'G':  // use gap function (between blocks).
                          bl_nmbr--; s=start_prof[bl_nmbr];
                          gapfunct=gpen[bl_nmbr];
                          gapend = i_gap = gapstart;
                          mxm=gapfunct[0]+MAT[j][i_gap];
                          for(g=0; g <= gapend; g++,i_gap--){
                             gap=gapfunct[g]+MAT[j][i_gap];
                             if(gapfunct[g]+MAT[j][i_gap]>=mxm){best_i=i_gap; state='M';mxm=gap;}
                             gap=gapfunct[g]+DEL[j][i_gap];
                             if(gapfunct[g]+DEL[j][i_gap]>=mxm){best_i=i_gap; state='D';mxm=gap;}
                          }   
                          g = (gapstart-best_i);
                          for(l=0;l < g; l++) { back_operation[k++]='i';}
                          i =best_i; 
                          break;

                default:  print_error("this should not happen");
           }
        }

	*oper_len = k-2;                           
        for(i=0;i<k;i++) operation[i]=back_operation[k-i-1];
        free(back_operation);free(start_prof);
        for(i=0;i<=prof_len;i++){free(MAT[i]);free(DEL[i]);free(INS[i]);free(score[i]);}
        free(MAT);free(DEL);free(INS); free(score);
	
        return operation;
}
