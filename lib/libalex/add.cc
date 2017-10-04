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


char    *ComplexAlnSMatrix(idpn_typ *idpn, Int4 seq_len, unsigned char *seq,
        Int4 Rpts, smx_typ *M, Int4 *Score, Int4 *Start, Int4 *Oper_len)
{
        Int4    pen,i,j,k,l,im1,jm1,i0,n,r,s,t,prof_len;
        Int4    maxMAT,maxDEL,alph_len,total,g,s1,s0;
        Int4    *block_lengths, *start_prof;
        Int4    **MAT,**DEL,**INS,**score;
        Int4    *alignscore=Score; 
        Int4    *start=Start;
        Int4    *oper_len=Oper_len;
        Int4    nblks = Rpts*idpn->nBlks();
         
        assert(Rpts <= idpn->MaxRpts());
        MEW(block_lengths, nblks+1, Int4);
        MEW(start_prof, nblks+2, Int4);

        for(n=0,i=1;i<=nblks;i++){
                start_prof[i]=n+1;
                block_lengths[i] = LenSMatrix(M[i]);
                n+=block_lengths[i];
        } start_prof[nblks+1]=0;
        prof_len=n; 

        Int4    *grease_pen=0;
        grease_pen=idpn->SeqPenalties(seq,seq_len); // returns 0 if ixal==0 && ixah==0;
        alph_len = nAlpha(SMatrixA(M[1]));
        NEWP(MAT,prof_len+1,Int4); NEWP(DEL,prof_len+1,Int4);
        NEWP(INS,prof_len+1,Int4);
        for(i=0;i<=prof_len;i++){
                MEW(MAT[i],seq_len+1,Int4); MEW(DEL[i],seq_len+1,Int4);
                MEW(INS[i],seq_len+1,Int4);
        }
        NEWP(score,prof_len+1,Int4);
        for(i=0;i<=prof_len;i++) NEW(score[i],alph_len+1,Int4);
        for(total=0,i=1;i<=nblks;i++){
           for(s=1;s<=block_lengths[i];s++){
                for(r=0;r<=alph_len;r++){
                   score[total+s][r]=ValSMatrix(s,r,M[i]);
                }
           } total+=block_lengths[i];
        }

        Int4 *io=idpn->InsOpen();
        Int4 *ie=idpn->InsExtend();
        Int4 *od=idpn->DelOpen();
        Int4 *de=idpn->DelExtend();
        Int4    **gpen=idpn->GapFunct();

        DEL[0][0]=INS[0][0]=MAT[0][0]=0;
        DEL[1][0]=INS[1][0]=MAT[1][0]=od[1]+de[1];
        for(s=2,jm1=1,j=2;j<=prof_len;jm1++,j++) {
                if(j!=start_prof[s]){ MAT[j][0]=INS[j][0]=DEL[j][0]=DEL[jm1][0]+de[j];
                } else {
                   if(gpen[s-1]){
                        MAT[j][0]=DEL[jm1][0]+od[j]+de[j]+gpen[s-1][0]; 
                        DEL[j][0]=DEL[jm1][0]+od[j]+de[j]+gpen[s-1][0];
                        INS[j][0]=DEL[jm1][0]+od[j]+de[j]+gpen[s-1][0];
                   } else {
                        MAT[j][0]=DEL[jm1][0]+od[j]+de[j]; 
                        DEL[j][0]=DEL[jm1][0]+od[j]+de[j];
                        INS[j][0]=DEL[jm1][0]+od[j]+de[j];
                   } s++; 
                }
        } for(i=1;i<=seq_len;i++) {DEL[0][i] = od[1];  MAT[0][i] = 0; INS[0][i] = 0;}

        register Int4 J,Jm1;
        for(s=2,J=1,Jm1=0;J<=prof_len;J++,Jm1++){
           Int4 End = start_prof[s-1]+block_lengths[s-1]-1;
           if(J!=start_prof[s] || gpen[s-1]==0){
              if((grease_pen && J!=End)){
                complex_affine_aln_smx(MAT[J],MAT[Jm1],INS[J],INS[Jm1],DEL[J],DEL[Jm1],
                  od[J],de[J],io[J],ie[J],seq_len,seq,score[J],grease_pen);
              } else { // don't use grease penalities between blocks...
                complex_affine_aln_smx(MAT[J],MAT[Jm1],INS[J],INS[Jm1],DEL[J],DEL[Jm1],
                  od[J],de[J],io[J],ie[J],seq_len,seq,score[J]);
              }
           } else { // impose gap penalties between blocks when at first position...
              complex_gap_align_smx(MAT[J],MAT[Jm1],INS[J],DEL[J],DEL[Jm1],
                od[J]+de[J],io[J],ie[J],seq_len,seq,score[J],gpen[s-1],grease_pen);
           }
           if(J==start_prof[s]) s++;
        }

        Int4    bl_nmbr,best_i,inse;
        Int4    *gapfunct,gapend,i_gap,gapstart;
        Int4    mxm, gap;
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
                        } i=i0;*alignscore=mxm;
                        break;
                case 'M': // previously sampled matched state.
                        if(j <= 1) { back_operation[k++]='M'; state='B'; break; }
                        if(j==s){
                           bl_nmbr--; s=start_prof[bl_nmbr];
                             j--; back_operation[k++]='M'; 
                             state='G'; gapstart=i-1;  
                        } else { // else we are within a block & use affine gaps.
                           j--; i--; back_operation[k++]='m';
                           if(MAT[j][i]>INS[j][i] && MAT[j][i]>DEL[j][i]) state='M';
                           else {if(INS[j][i]>DEL[j][i]) state='I';
                                 else state='D';
                           }
                        } break;
                case 'D': // previously sampled deletion
                        if(j<=1){ back_operation[k++]='D'; state='B'; i++; break; }
                        if(j==s){
                            j--; back_operation[k++]='D'; 
                            bl_nmbr--; s=start_prof[bl_nmbr];
                            state='G'; gapstart=i; 
                        } else { j--; back_operation[k++]='d';
                             if(DEL[j][i]+de[j]>=MAT[j][i]+od[j]+de[j]) state='D';
                             else state='M';
                        } break;
                case 'I': // previously sampled insertion.
                          back_operation[k++]='I'; i--;
                          if(grease_pen){
                                if(j==start_prof[bl_nmbr]+block_lengths[bl_nmbr]-1) {inse=0;}
                                else { inse = ie[j] + grease_pen[i]; }
                          } else inse = ie[j];
                          if (INS[j][i]+inse > MAT[j][i]+io[j]+inse) state='I';
                          else state='M'; 
                          break;                   
                case 'G':  // use gap function (between blocks).
                  if(gpen[bl_nmbr]){
                      gapfunct=gpen[bl_nmbr];
                      gapend = i_gap = gapstart;
                      mxm=gapfunct[0]+MAT[j][i_gap];
                      for(g=0; g <= gapend; g++,i_gap--){
                         if(gapfunct[g] <= SHRT_MIN) break; // fixes bug???
                         gap=gapfunct[g]+MAT[j][i_gap];
                         if(gapfunct[g]+MAT[j][i_gap]>=mxm) {best_i=i_gap; state='M';mxm=gap;}
                         gap=gapfunct[g]+DEL[j][i_gap];
                         if(gapfunct[g]+DEL[j][i_gap]>=mxm) {best_i=i_gap; state='D';mxm=gap;}
                      } g = (gapstart-best_i);
                  } else {
                      gapend = i_gap = gapstart;
                      mxm=MAT[j][i_gap];
                      for(g=0; g <= gapend; g++,i_gap--){
                         gap=MAT[j][i_gap];
                         if(MAT[j][i_gap]>=mxm) {best_i=i_gap; state='M';mxm=gap;}
                         gap=DEL[j][i_gap];
                         if(DEL[j][i_gap]>=mxm) {best_i=i_gap; state='D';mxm=gap;}
                      } g=(gapstart-best_i);
                  }
                  for(l=0;l < g; l++) { back_operation[k++]='i';}
                  i =best_i; 
                  break;
                default:  print_error("this should not happen");
           }
        } *oper_len = k-2;
        for(i=0;i<k;i++) operation[i]=back_operation[k-i-1];
        free(back_operation);free(block_lengths);free(start_prof);
        for(i=0;i<=prof_len;i++){free(MAT[i]);free(DEL[i]);free(INS[i]);}
        free(MAT);free(DEL);free(INS);
        for(i=0;i<=prof_len;i++) free(score[i]);
        free(score);
        if(grease_pen) free(grease_pen);
// std::cerr << operation; std::cerr << std::endl;
        return operation;
}

