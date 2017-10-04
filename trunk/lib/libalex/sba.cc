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

#include "sba.h"

static const Int4 dim_weights = 20;

sba_typ::sba_typ() { init(); }

sba_typ::sba_typ(Int4 sequence_len, unsigned char *sequence, Int4 *io, 
		Int4 *ie, Int4 *od, Int4 *de, Int4 nblocks, smx_typ *M, 
			Int4 **gapfnct,double pernats)
{ initialize(sequence_len,sequence,io,ie,od,de,nblocks,M,gapfnct,pernats); }


void    sba_typ::initialize(Int4 sequence_len, unsigned char *sequence, Int4 *io, 
		Int4 *ie, Int4 *od, Int4 *de, Int4 nblocks, smx_typ *M, 
			Int4 **gapfnct,double pernats)
{
        Int4 i,j,n=0,r,s,g,total=0;

	nmbr_of_blocks=nblocks;
        seq_len=sequence_len;
        seq=sequence;

        NEW(block_lengths, nblocks+1, Int4);
        NEW(start_prof, nblocks+2, Int4);
	NEW(maxgap,nblocks+1,Int4);
        for(i=1;i<=nblocks;i++){
                start_prof[i]=n+1;
                block_lengths[i] = LenSMatrix(M[i]);
                n+=block_lengths[i];
        }
        prof_len=n; 
        a_type A = SMatrixA(M[1]);
        alph_len=nAlpha(A);

        NEWP(gpen,nblocks+1,double);
        for(i=0;i<=nblocks;i++) NEW(gpen[i],seq_len+1,double);

        Int4    gapscore;
        for(s=1;s < nblocks;s++){
           for(g=0;g<=seq_len;g++){
                gapscore=gapfnct[s][g];
                if(gapscore <= SHRT_MIN) {
			gpen[s][g]=0.0;
			maxgap[s]=g;
			break;
		} else gpen[s][g]=exp((double)gapscore/pernats);
           }
        }

        NEWP(MAT,seq_len+1,double);
        NEWP(INSO,seq_len+1,double);NEWP(INSE,seq_len+1,double);
        NEWP(DELO,seq_len+1,double);NEWP(DELE,seq_len+1,double);
        for(i=0;i<=seq_len;i++){
                NEW(MAT[i],prof_len+1,double);
                NEW(INSO[i],prof_len+1,double);
                NEW(INSE[i],prof_len+1,double);
                NEW(DELO[i],prof_len+1,double);
                NEW(DELE[i],prof_len+1,double);
        }

        NEWP(match,alph_len+1,double);
        for(i=0;i<=alph_len;i++) NEW(match[i],prof_len+1,double);

        NEW(inso,prof_len+1,double);NEW(inse,prof_len+1,double);
        NEW(delo,prof_len+1,double);NEW(dele,prof_len+1,double);

        for(j=1;j<=prof_len;j++){
                inso[j]=exp((double)io[j]/pernats);
                inse[j]=exp((double)ie[j]/pernats);
                delo[j]=exp((double)od[j]/pernats);
                dele[j]=exp((double)de[j]/pernats);
        }

        for(i=1;i<=nblocks;i++){
           for(s=1;s<=block_lengths[i];s++){
                for(r=0;r<=alph_len;r++)
                   match[r][total+s]=exp((double)ValSMatrix(s,r,M[i])/pernats);
           }
           total+=block_lengths[i];
        }
        store_i = NULL; storeState = NULL; storeLike = NULL; 
        Init( );
}

static void     put_dp_matrix(double **D, Int4 len, Int4 prolen)
{
        Int4    r,m,v,start,end;

  for(start=0,end=14; start <= prolen; start+=15,end+=15) {
        fprintf(stderr,"     |");
        for(m=start; m<=end && m <= prolen; m++) { fprintf(stderr,"%3d  |",m); }
        fprintf(stderr,"\n");
        for(r=0; r<=len; r++) {
           fprintf(stderr,"%5d|",r);
           for(m=start; m<=end && m <= prolen; m++) {
		if(D[r][m] == 0) fprintf(stderr,"  -  ");
		else {
			v = (Int4) log10(D[r][m]);
                	fprintf(stderr,"%5d",v);
		}
                fprintf(stderr,"|");
           }
           fprintf(stderr,"\n");
        }
        fprintf(stderr,"\n");
  }
}

BooLean	sba_typ::Init( )
{
	Int4	i,j,k,s,g,r,im1,jm1,i0,end;
	double	sumMAT,sumDEL,vM,vDO,vDE;
        Int4            *Store_i;
        char            *StoreState;
        double          *StoreLike;

        for(i=0;i<=seq_len;i++) MAT[i][0]=1.0; 
        DELO[0][1]=delo[1];  DELE[0][2]=delo[1]*dele[2];
        for(j=3;j<=block_lengths[1];j++) DELE[0][j]=DELE[0][j-1]*dele[j];
	s=2;
	for(j=block_lengths[1]+1;j<=prof_len;j++){
	   if(j==start_prof[s]){
		DELO[0][j]=DELE[0][j-1]*gpen[s-1][0]*delo[j];
			s++;
	   } else if (j==start_prof[s]+1) DELE[0][j]=DELO[0][j-1]*dele[j];
	   else DELE[0][j]=DELO[0][j-1]*dele[j];
	}
        
        for(im1=0,i=1;i<=seq_len;im1++,i++){
            r = seq[i];
	    s=2;
            for(j=1,jm1=0;j<=prof_len;j++,jm1++){
		if(j!=start_prof[s]){
		    MAT[i][j]=match[r][j]*(MAT[im1][jm1]+INSO[im1][jm1]+
                               	INSE[im1][jm1]+DELO[im1][jm1]+DELE[im1][jm1]);
                    INSO[i][j]=inso[j]*MAT[im1][j];
                    INSE[i][j]=inse[j]*(INSO[im1][j]+INSE[im1][j]);
                    DELO[i][j]=delo[j]*MAT[i][jm1];
                    DELE[i][j]=dele[j]*(DELO[i][jm1]+DELE[i][jm1]);
		} else{
		   sumMAT=0.;sumDEL=0.;
		  // end = MINIMUM(Int4,i,maxgap[s-1]);
                   for(i0=im1,g=0;g < i;g++,i0--)
			sumMAT+= gpen[s-1][g]*(MAT[i0][jm1]
				+DELO[i0][jm1]+DELE[i0][jm1]);
		   for(i0=i,g=0;g <= i;g++,i0--)
		      	sumDEL+=gpen[s-1][g]*(MAT[i0][jm1]
				+DELO[i0][jm1]+DELE[i0][jm1]);
                   MAT[i][j]=match[r][j]*sumMAT;
                   DELO[i][j]=delo[j]*sumDEL; DELE[i][j]=0;
                   INSO[i][j]=inso[j]*MAT[im1][j];
		   INSE[i][j]=inse[j] *(INSO[im1][j]+INSE[im1][j]);
		   s++;
		}
	    }
	}
#if 0
	fprintf(stderr,"====================== MAT: ====================\n");
	put_dp_matrix(MAT, seq_len, prof_len);
	fprintf(stderr,"====================== INSO: ====================\n");
	put_dp_matrix(INSO, seq_len, prof_len);
	fprintf(stderr,"====================== INSE: ====================\n");
	put_dp_matrix(INSE, seq_len, prof_len);
	fprintf(stderr,"====================== DELO: ====================\n");
	put_dp_matrix(DELO, seq_len, prof_len);	
	fprintf(stderr,"====================== DELE: ====================\n");
	put_dp_matrix(DELE, seq_len, prof_len);
#endif
	dh_type	Heap=dheap(3*seq_len,3);
	j=prof_len; totalLike=0.0;
	NEW(store_i,3*seq_len + 3, Int4);
	NEW(storeState,3*seq_len + 3, char);
	NEW(storeLike,3*seq_len + 3, double);
	NEW(Store_i,3*seq_len + 3, Int4);
	NEW(StoreState,3*seq_len + 3, char);
	NEW(StoreLike,3*seq_len + 3, double);
	Int4 item;

	for(item=1,i=1; i <= seq_len; i++)
	    totalLike+=MAT[i][j]+DELO[i][j]+DELE[i][j];
	for(item=1,i=1; i <= seq_len; i++) {
	    vM=MAT[i][j]/totalLike; vDO=DELO[i][j]/totalLike; vDE=DELE[i][j]/totalLike;
	    insrtHeap(item, -vM, Heap); 
	    Store_i[item] = i; StoreState[item] = 'M'; StoreLike[item]=vM;
	    item++;
	    insrtHeap(item, -vDO, Heap); 
	    Store_i[item] = i; StoreState[item] = 'D'; StoreLike[item]=vDO;
	    item++;
	    insrtHeap(item, -vDE, Heap); 
	    Store_i[item] = i; StoreState[item] = 'd'; StoreLike[item]=vDO;
	    item++;
	}
	Int4	item2;
	for(item2=1; (item=delminHeap(Heap)) != NULL; item2++){
	    store_i[item2] = Store_i[item];; 
	    storeState[item2] = StoreState[item]; 
	    storeLike[item2]=StoreLike[item];
	}
	num_stored=item2-1;
	free(Store_i); free(StoreState); free(StoreLike);
	Nildheap(Heap);
	if(totalLike >= HUGE_VAL) print_error("SBA( ) overflow error");
	return 1;
}

sba_typ::~sba_typ( ) { Free(); }

void sba_typ::Free( )
{
	Int4 i;

	if(weight!=NULL) {
	  for(i=0; i<=dim_weights; i++) free(weight[i]); free(weight);
	}
        for(i=0;i<=nmbr_of_blocks;i++) free(gpen[i]);

        free(gpen);free(block_lengths);free(start_prof);free(maxgap);
                                        
        for(i=0;i<=seq_len;i++){
                free(MAT[i]);free(INSO[i]);
                free(INSE[i]);free(DELO[i]); free(DELE[i]);
        }
	free(store_i); free(storeState); free(storeLike);
        free(MAT);free(INSO);free(INSE);free(DELO);free(DELE);
        for(i=0;i<=alph_len;i++) free(match[i]);
        free(match);
        free(inso);free(inse);
        free(delo);free(dele);
}

char    *sba_typ::TraceBack(Int4 *alignscore, Int4 *trace_length, Int4 *start,
                        double pernats, double *RE)
{
	Int4	N,i,j,k,l,i0,s,t,r,n,oper_len,bl_nmbr, temp;
	double	total, rand_num, score, sum;
	char	state;
	double	*GAP,*gapfunct;
	Int4	gapend,i_gap,g,gapstart;
	char	*operation, *back_operation;
	double	p,q;

	*RE=0.0;
	NEW(GAP,seq_len+1,double);
        NEW(operation,seq_len+prof_len+3,char);
        NEW(back_operation,seq_len+prof_len+3,char);
	j=prof_len; i=seq_len;
	for(score=1.0,state='E',k=0 ; state != 'X'; ){
	   switch(state){
                case 'B': // begin state (exit traceback).
		   back_operation[k++]='E'; *start=i; state='X'; break;
                case 'E': // end state.  
		   back_operation[k++]='E';
		   bl_nmbr=nmbr_of_blocks; s=start_prof[bl_nmbr];
		   gapfunct=gpen[bl_nmbr];

	q=1.0/(double)(3*seq_len);
	for(total=0.,n=1;n<=seq_len;n++){
		p=MAT[n][j]/totalLike;
		if(p > 0.0) total += p*log(p/q);
		p=DELO[n][j]/totalLike;
	       	if(p > 0.0) total += p*log(p/q);
		p=DELE[n][j]/totalLike;
                if(p > 0.0) total += p*log(p/q);
	} *RE+=total;

		   rand_num=totalLike*SampleUniformProb();
		   do {
            	     if((rand_num-=MAT[i][j]) <= 0.0){ state='M'; break; }
		     else if((rand_num-=DELO[i][j]) <= 0.0){ state='D'; break; }
		     else if((rand_num-=DELE[i][j]) <= 0.0){ state='d'; break; }
		   } while(i-- > 0);
		   break;
                case 'M': // previously sampled matched state.
			r=seq[i]; score*=match[r][j]; 

			if(j <= 1) { back_operation[k++]='M'; state='B'; break; }
			if(j==s){ j--; back_operation[k++]='M'; state='G'; gapstart=i-1; }
			else { // else we are within a block & use affine gaps.
			   j--; back_operation[k++]='m'; i--; 
                           rand_num = MAT[i][j] + INSO[i][j]+
				INSE[i][j]+DELO[i][j]+DELE[i][j];

	q=1.0/5.; total=0;
	p=MAT[i][j]/rand_num;
	if (p>0.0) total += p*log(p/q);
	p=INSO[i][j]/rand_num;
        if (p>0.0) total += p*log(p/q);
	p=INSE[i][j]/rand_num;
        if (p>0.0) total += p*log(p/q);
	p=DELO[i][j]/rand_num;
        if (p>0.0) total += p*log(p/q);
	p=DELE[i][j]/rand_num;
        if (p>0.0) total += p*log(p/q);
	*RE+=total;

                           rand_num *= SampleUniformProb();
                           if((rand_num-=INSO[i][j])<=0.0) {state='I';}
                           else if((rand_num-=INSE[i][j])<=0.0) {state='i';}
                           else if((rand_num-=DELO[i][j])<=0.0) {state='D';}
                           else if((rand_num-=DELE[i][j])<=0.0) {state='d';}
                           else { state='M'; }
			} break;
                          
		case 'D': // previously sampled deletion opening.
			  score*=delo[j];
                          if(j<=1){ back_operation[k++]='D'; state='B'; i++; break; }
                          if(j==s){ j--; back_operation[k++]='D'; state='G'; gapstart=i; }
                          else { j--; back_operation[k++]='d'; state='M'; }
			  break;

                case 'd': // previously sampled deletion extension.
			  score*=dele[j];

			  if(j <= 1) print_error("this should not happen");
			  if(j == s) print_error("this should not happen");
			  j--; back_operation[k++]='d'; 
                          rand_num=DELO[i][j]+DELE[i][j];

	q=1.0/2.; total=0.,
	p=DELO[i][j]/rand_num;
	if (p>0.0) total += p*log(p/q);
	p=DELE[i][j]/rand_num;
        if (p>0.0) total += p*log(p/q);
	*RE+=total;

                          rand_num*=SampleUniformProb();
                          if ((rand_num-=DELO[i][j])<0.0) state='D';
                          else state='d';
                          break;
                      
                case 'I': // previously sampled insertion opening.
			  score*=inso[j]; i--; 
			  back_operation[k++]='I'; state='M'; 
			  break;
                   
                case 'i': // previously sampled insertion extension.
			  score*=inse[j]; i--;
			  back_operation[k++]='I'; 

                          rand_num=INSO[i][j]+INSE[i][j];

	q=1.0/2.0; total=0.,
        p=INSO[i][j]/rand_num;
        if (p>0.0) total += p*log(p/q);
        p=INSE[i][j]/rand_num;
        if (p>0.0) total += p*log(p/q);
        *RE+=total;

                          rand_num*=SampleUniformProb();
                          if((rand_num-=INSO[i][j]) < 0.0) state='I';
                          else state='i'; break;
		case 'G':  // use gap function from match (between blocks).
			  bl_nmbr--; s=start_prof[bl_nmbr];
		   	  gapfunct=gpen[bl_nmbr];
			  gapend = i_gap = gapstart;
			  for(N=0,sum=0.,g=0; g <= gapend; g++,i_gap--){
			    GAP[g]=gapfunct[g]*(MAT[i_gap][j]
				+DELO[i_gap][j]+DELE[i_gap][j]);
			    sum+=GAP[g];
			    N++;
			  } 
// Add relative entropy calc. here...
	i_gap = gapstart;
	  q = 1.0/(double)(3*N);
	  for(total=0.,g=0; g <= gapend; g++,i_gap--){
	  p = gapfunct[g]*MAT[i_gap][j]/sum;
	  if(p > 0.0) total += p*log(p/q);
	  p = gapfunct[g]*DELO[i_gap][j]/sum;
	  if(p > 0.0) total += p*log(p/q);
	  p = gapfunct[g]*DELE[i_gap][j]/sum;
	  if(p > 0.0) total += p*log(p/q);
	} *RE += total;
// end RE...
			  
			  rand_num=sum*SampleUniformProb();
			  i_gap = gapstart;
			  for(g=0;g <= gapend; g++,i_gap--){
	                          if ((rand_num-=GAP[g])<=0.0){
        	                  	for(l=0;l < g; l++) {back_operation[k++]='i';}
                	                score*=gapfunct[g];
                                	break;
                          	  }
			  } i =i_gap; 
                          sum=MAT[i][j]+DELO[i][j]+DELE[i][j];
                          rand_num=sum*SampleUniformProb();
                          if ((rand_num-=DELO[i][j])<0.0) state='D';   
                          else if ((rand_num-=DELE[i][j])<0.0) state='d';
                          else state='M';
                          break;

                default:  print_error("this should not happen");
	   }
	}
        *alignscore=(Int4) (floor(pernats*log(score) +0.5));
        for(i=0;i < k;i++) operation[i]=back_operation[k-i-1];
        *trace_length=k;            
	*RE = *RE/(double)(k-2);
	
        free(GAP);free(back_operation);
        return operation;
}

char	*sba_typ::TraceBack(Int4 *alignscore, Int4 *trace_length, Int4 *start, 
	double pernats)
{
	Int4	i,j,k,l,i0,s,t,r,oper_len,bl_nmbr, temp;
	double	rand_num, score, sum;
	char	state;
	double	*GAP,*gapfunct;
	Int4	gapend,i_gap,g,gapstart;
	char	*operation, *back_operation;

	NEW(GAP,seq_len+1,double);
        NEW(operation,seq_len+prof_len+3,char);
        NEW(back_operation,seq_len+prof_len+3,char);
	j=prof_len;i=seq_len;
	for(score=1.0,state='E',k=0 ; state != 'X'; ){
	   switch(state){
                case 'B': // begin state (exit traceback).
		   back_operation[k++]='E'; *start=i; state='X'; break;
                case 'E': // end state.  
		   back_operation[k++]='E';
		   bl_nmbr=nmbr_of_blocks; s=start_prof[bl_nmbr];
		   gapfunct=gpen[bl_nmbr];
		   rand_num=totalLike*SampleUniformProb();
		   do {
            	     if((rand_num-=MAT[i][j]) <= 0.0){ state='M'; break; }
		     else if((rand_num-=DELO[i][j]) <= 0.0){ state='D'; break; }
		     else if((rand_num-=DELE[i][j]) <= 0.0){ state='d'; break; }
		   } while(i-- > 0);
		   break;
                case 'M': // previously sampled matched state.
			r=seq[i]; score*=match[r][j]; 
//fprintf(stderr,"match = %d\n",(Int4) (floor(pernats*log(match[r][j])+0.5)));

			if(j <= 1) { back_operation[k++]='M'; state='B'; break; }
			if(j==s){ j--; back_operation[k++]='M'; state='G'; gapstart=i-1; }
			else { // else we are within a block & use affine gaps.
			   j--; back_operation[k++]='m'; i--; 
                           rand_num = MAT[i][j]+INSO[i][j]+
				INSE[i][j]+DELO[i][j]+DELE[i][j];
                           rand_num *= SampleUniformProb();
                           if((rand_num-=INSO[i][j])<=0.0) {state='I';}
                           else if((rand_num-=INSE[i][j])<=0.0) {state='i';}
                           else if((rand_num-=DELO[i][j])<=0.0) {state='D';}
                           else if((rand_num-=DELE[i][j])<=0.0) {state='d';}
                           else { state='M'; }
			} break;
                          
		case 'D': // previously sampled deletion opening.
			  score*=delo[j];
//fprintf(stderr,"delo= %d\n",(Int4) (floor(pernats*log(delo[j])+0.5)));
                          if(j<=1){ back_operation[k++]='D'; state='B'; i++; break; }
                          if(j==s){ j--; back_operation[k++]='D'; state='G'; gapstart=i; }
                          else { j--; back_operation[k++]='d'; state='M'; }
			  break;

                case 'd': // previously sampled deletion extension.
			  score*=dele[j];
//fprintf(stderr,"dele = %d\n",(Int4) (floor(pernats*log(dele[j])+0.5)));

			  if(j <= 1) print_error("this should not happen");
			  if(j == s) print_error("this should not happen");
			  j--; back_operation[k++]='d'; 
                          rand_num=DELO[i][j]+DELE[i][j];
                          rand_num*=SampleUniformProb();
                          if ((rand_num-=DELO[i][j])<0.0) state='D';
                          else state='d';
                          break;
                      
                case 'I': // previously sampled insertion opening.
			  score*=inso[j]; i--; 
//fprintf(stderr,"inso = %d\n",(Int4) (floor(pernats*log(inso[j])+0.5)));
			  back_operation[k++]='I'; state='M'; 
			  break;
                   
                case 'i': // previously sampled insertion extension.
			  score*=inse[j]; i--;
//fprintf(stderr,"inse = %d\n",(Int4) (floor(pernats*log(inse[j])+0.5)));
			  back_operation[k++]='I'; 

                          rand_num=INSO[i][j]+INSE[i][j];
                          rand_num*=SampleUniformProb();
                          if((rand_num-=INSO[i][j]) < 0.0) state='I';
                          else state='i'; break;
		case 'G':  // use gap function from match (between blocks).
			  bl_nmbr--; s=start_prof[bl_nmbr];
		   	  gapfunct=gpen[bl_nmbr];
			  gapend = i_gap = gapstart;
			  for(sum=0.,g=0; g <= gapend; g++,i_gap--){
			    GAP[g]=gapfunct[g]*(MAT[i_gap][j]
				+DELO[i_gap][j]+DELE[i_gap][j]);
			    sum+=GAP[g];
			  } 
			  rand_num=sum*SampleUniformProb();
			  i_gap = gapstart;
			  for(g=0;g <= gapend; g++,i_gap--){
	                          if ((rand_num-=GAP[g])<=0.0){
        	                  	for(l=0;l < g; l++) {back_operation[k++]='i';}
                	                score*=gapfunct[g];
                                	break;
                          	  }
			  } i =i_gap; 
//for( ;g > 0 ; g--) fprintf(stderr,"state = i; i=%d; j=%d\n",i+g-1,j);
                          sum=MAT[i][j]+DELO[i][j]+DELE[i][j];
                          rand_num=sum*SampleUniformProb();
                          if ((rand_num-=DELO[i][j])<0.0) state='D';   
                          else if ((rand_num-=DELE[i][j])<0.0) state='d';
                          else state='M';
                          break;

                default:  print_error("this should not happen");
	   }
#if 0
           fprintf(stderr,"state = %c; i=%d; j=%d; (s=%d)",state,i,j,s);
           fprintf(stderr,"score= %d; ", (Int4) (floor(pernats*log(score) +0.5)));
           fprintf(stderr,"back_operation = %c\n", back_operation[k-1]);
#endif
	}
        *alignscore=(Int4) (floor(pernats*log(score) +0.5));
//	fprintf(stderr,"score = %d\n%s\n",*alignscore,back_operation);
                          
        for(i=0;i < k;i++) operation[i]=back_operation[k-i-1];
//	fprintf(stderr,"operation = %s\n",operation);
        *trace_length=k;            
        free(GAP);free(back_operation);
        return operation;
}
char    *sba_typ::TraceBack(Int4 *alignscore, Int4 *trace_length, Int4 *start, 
        double pernats, double tm)
{
        Int4    i,j,k,l,n,i0,s,t,r,oper_len,bl_nmbr, temp;
        double  rand_num, score, sum;
        char    state;
        double  *GAP,*gapfunct,maxLike;
        Int4    gapend,i_gap,g,gapstart;
        char    *operation, *back_operation;
        double  m,ino,ine,dlo,dle,sumX;
        double  *storePowLike, sm;

        NEW(GAP,seq_len+1,double);
        NEW(operation,seq_len+prof_len+3,char);
        NEW(back_operation,seq_len+prof_len+3,char);
        j=prof_len;i=seq_len;
        for(score=1.0,state='E',k=0 ; state != 'X'; ){
           switch(state){
                case 'B': // begin state (exit traceback).
                   back_operation[k++]='E'; *start=i; state='X'; break;
                case 'E': // end state.  
                   back_operation[k++]='E';
                   bl_nmbr=nmbr_of_blocks; s=start_prof[bl_nmbr];
                   gapfunct=gpen[bl_nmbr];
                   NEW(storePowLike,3*seq_len + 3, double);

                   sumX = maxLike=storePowLike[1]=pow(storeLike[1],tm);
                   for(n=2; n <=num_stored; n++){
                        storePowLike[n]=pow(storeLike[n],tm);
                        sumX += storePowLike[n];
                        if((storePowLike[n]/maxLike) < 0.001) { break; }
                   }

                   sumX*=SampleUniformProb();
                   for(l=1;l<=n;l++){
                        if((sumX-=storePowLike[l])<=0.0) {
                                state = storeState[l]; i=store_i[l];
break;
                        }
                   }
                   free(storePowLike); break;
                case 'M': // previously sampled matched state.
                        r=seq[i]; score*=match[r][j]; 
//fprintf(stderr,"match = %d\n",(Int4) (floor(pernats*log(match[r][j])+0.5)));

                        if(j <= 1) { back_operation[k++]='M'; state='B'; break; }
                        if(j==s){ j--; back_operation[k++]='M'; state='G';
gapstart=i-1; }
                        else { // else we are within a block & use affine gaps.
                           j--; back_operation[k++]='m'; i--; 
			   sm=MAT[i][j]+INSO[i][j]+INSE[i][j]+DELO[i][j]+DELE[i][j];
                           m=pow(MAT[i][j]/sm,tm);
                           ino=pow(INSO[i][j]/sm,tm);
                           ine=pow(INSE[i][j]/sm,tm);
                           dlo=pow(DELO[i][j]/sm,tm);
                           dle=pow(DELE[i][j]/sm,tm);
                           rand_num = m+ino+ine+dlo+dle;
                           rand_num *= SampleUniformProb();
                           if((rand_num-=ino)<=0.0) {state='I';}
                           else if((rand_num-=ine)<=0.0) {state='i';}
                           else if((rand_num-=dlo)<=0.0) {state='D';}
                           else if((rand_num-=dle)<=0.0) {state='d';}
                           else { state='M'; }
                        } break;
                          
                case 'D': // previously sampled deletion opening.
                          score*=delo[j];
//fprintf(stderr,"delo= %d\n",(Int4) (floor(pernats*log(delo[j])+0.5)));
                          if(j<=1){ back_operation[k++]='D'; state='B'; i++;
break; }
                          if(j==s){ j--; back_operation[k++]='D'; state='G';
gapstart=i; }
                          else { j--; back_operation[k++]='d'; state='M'; }
                          break;

                case 'd': // previously sampled deletion extension.
                          score*=dele[j];
//fprintf(stderr,"dele = %d\n",(Int4) (floor(pernats*log(dele[j])+0.5)));

                          if(j <= 1) print_error("this should not happen");
                          if(j == s) print_error("this should not happen");
                          j--; back_operation[k++]='d';
			  sm=DELO[i][j]+DELE[i][j]; 
                          dlo=pow(DELO[i][j]/sm,tm);
                          dle=pow(DELE[i][j]/sm,tm); 
                          rand_num=dlo+dle;
                          rand_num*=SampleUniformProb();
                          if ((rand_num-=dlo)<0.0) state='D';
                          else state='d';
                          break;
                      
                case 'I': // previously sampled insertion opening.
                          score*=inso[j]; i--; 
//fprintf(stderr,"inso = %d\n",(Int4) (floor(pernats*log(inso[j])+0.5)));
                          back_operation[k++]='I'; state='M'; 
                          break;
                   
                case 'i': // previously sampled insertion extension.
                          score*=inse[j]; i--;
//fprintf(stderr,"inse = %d\n",(Int4) (floor(pernats*log(inse[j])+0.5)));
                          back_operation[k++]='I';
			  sm=INSO[i][j]+INSE[i][j]; 
                          ino=pow(INSO[i][j]/sm,tm);
                          ine=pow(INSE[i][j]/sm,tm);
                          rand_num = ino+ine;
                          rand_num*=SampleUniformProb();
                          if((rand_num-=ino) < 0.0) state='I';
                          else state='i'; break;
                case 'G':  // use gap function from match (between blocks).
                          bl_nmbr--; s=start_prof[bl_nmbr];
                          gapfunct=gpen[bl_nmbr];
                          gapend = i_gap = gapstart; sm=0.;
			  for(rand_num=0.,g=0; g <= gapend; g++,i_gap--){
			  	sm+=gapfunct[g]*(MAT[i_gap][j]+DELO[i_gap][j]+DELE[i_gap][j]);
			  }			  
			  gapend = i_gap = gapstart;
                          for(rand_num=0.,g=0; g <= gapend; g++,i_gap--){
                            GAP[g]=pow(gapfunct[g]*(MAT[i_gap][j]
                                +DELO[i_gap][j]+DELE[i_gap][j])/sm,tm);
                            rand_num+=GAP[g];
                          } 
                          rand_num*=SampleUniformProb();
                          i_gap = gapstart;
                          for(g=0;g <= gapend; g++,i_gap--){
                                  if ((rand_num-=GAP[g])<=0.0){
                                        for(l=0;l < g; l++)
{back_operation[k++]='i';}
                                        score*=gapfunct[g];
                                        break;
                                  }
                          } i =i_gap; 
//for( ;g > 0 ; g--) fprintf(stderr,"state = i; i=%d; j=%d\n",i+g-1,j);
			  sm=MAT[i][j]+DELO[i][j]+DELE[i][j];
                          m=pow(MAT[i][j]/sm,tm);
                          dlo=pow(DELO[i][j]/sm,tm);
                          dle=pow(DELE[i][j]/sm,tm);
                          rand_num = m+dlo+dle;
                          rand_num*=SampleUniformProb();
                          if ((rand_num-=dlo)<0.0) state='D';   
                          else if ((rand_num-=dle)<0.0) state='d';
                          else state='M';
                          break;

                default:  print_error("this should not happen");
           }
#if 0
           fprintf(stderr,"state = %c; i=%d; j=%d; (s=%d)",state,i,j,s);
           fprintf(stderr,"score= %d; ", (Int4) (floor(pernats*log(score) +0.5)));
           fprintf(stderr,"back_operation = %c\n", back_operation[k-1]);
#endif
        }
        *alignscore=(Int4) (floor(pernats*log(score) +0.5));
//      fprintf(stderr,"score = %d\n%s\n",*alignscore,back_operation);
                          
        for(i=0;i < k;i++) operation[i]=back_operation[k-i-1];
//      fprintf(stderr,"operation = %s\n",operation);
        *trace_length=k;            
        free(GAP);free(back_operation);
        return operation;


}

char	*sba_typ::EMTraceBack(Int4 *alignscore, Int4 *trace_length, Int4 *start, 
	double pernats)
{
	Int4	i,j,k,l,i0,s,t,r,oper_len,bl_nmbr, temp;
	Int4	best_i;
	double	rand_num, score, sum;
	char	state;
	double	*GAP,*gapfunct;
	Int4	gapend,i_gap,g,gapstart;
	char	*operation, *back_operation;

	NEW(GAP,seq_len+1,double);
        NEW(operation,seq_len+prof_len+3,char);
        NEW(back_operation,seq_len+prof_len+3,char);
	j=prof_len; i=seq_len;
	for(score=1.0,state='E',k=0 ; state != 'X'; ){
	   switch(state){
                case 'B': // begin state (exit traceback).
		   back_operation[k++]='E'; *start=i; state='X'; break;
                case 'E': // end state.  
		   back_operation[k++]='E';
		   bl_nmbr=nmbr_of_blocks; s=start_prof[bl_nmbr];
		   gapfunct=gpen[bl_nmbr];
		   rand_num=0;
		   for(i=1;i<=seq_len;i++){
            	     if(rand_num<MAT[i][j]){ state='M'; rand_num=MAT[i][j];i0=i;}
		     if(rand_num<DELO[i][j]){ state='D'; rand_num=DELO[i][j];i0=i;}
		     if(rand_num<DELE[i][j]){ state='d';rand_num=DELE[i][j];i0=i;}
		   }
		   i=i0;
		   break;
                case 'M': // previously sampled matched state.
			r=seq[i]; score*=match[r][j]; 
//fprintf(stderr,"match = %d\n",(Int4) (floor(pernats*log(match[r][j])+0.5)));

			if(j <= 1) { back_operation[k++]='M'; state='B'; break; }
			if(j==s){ j--; back_operation[k++]='M'; state='G'; gapstart=i-1; }
			else { // else we are within a block & use affine gaps.
			   j--; back_operation[k++]='m'; i--; 
			   state='I';rand_num=INSO[i][j];
                           if(rand_num<INSE[i][j]) {state='i';rand_num=INSE[i][j];}
                           if(rand_num<DELO[i][j]) {state='D';rand_num=DELO[i][j];}
                           if(rand_num<DELE[i][j]) {state='d';rand_num=DELE[i][j];}
                           if(rand_num<MAT[i][j])  { state='M'; rand_num=MAT[i][j];}}
			break;
                          
		case 'D': // previously sampled deletion opening.
			  score*=delo[j];
//fprintf(stderr,"delo= %d\n",(Int4) (floor(pernats*log(delo[j])+0.5)));
                          if(j<=1){ back_operation[k++]='D'; state='B'; i++; break; }
                          if(j==s){ j--; back_operation[k++]='D'; state='G'; gapstart=i; }
                          else { j--; back_operation[k++]='d'; state='M'; }
			  break;

                case 'd': // previously sampled deletion extension.
			  score*=dele[j];
//fprintf(stderr,"dele = %d\n",(Int4) (floor(pernats*log(dele[j])+0.5)));

			  if(j <= 1) print_error("this should not happen");
			  if(j == s) print_error("this should not happen");
			  j--; back_operation[k++]='d'; 
                          if (DELE[i][j]<DELO[i][j]) state='D';
                          else state='d';
                          break;
                      
                case 'I': // previously sampled insertion opening.
			  score*=inso[j]; i--; 
//fprintf(stderr,"inso = %d\n",(Int4) (floor(pernats*log(inso[j])+0.5)));
			  back_operation[k++]='I'; state='M'; 
			  break;
                   
                case 'i': // previously sampled insertion extension.
			  score*=inse[j]; i--;
//fprintf(stderr,"inse = %d\n",(Int4) (floor(pernats*log(inse[j])+0.5)));
			  back_operation[k++]='I'; 

                          if(INSE[i][j]<INSO[i][j]) state='I';
                          else state='i'; break;
                                           
		case 'G':  // use gap function (between blocks).
			  bl_nmbr--; s=start_prof[bl_nmbr];
		   	  gapfunct=gpen[bl_nmbr];
			  gapend = i_gap = gapstart;
			  for(sum=0.,g=0; g <= gapend; g++,i_gap--){
			    GAP[g]=gapfunct[g]*(MAT[i_gap][j]
				+DELO[i_gap][j]+DELE[i_gap][j]);
			    sum+=GAP[g];
			  } 
			  rand_num=0;
			  i_gap = gapstart;
			  best_i=gapstart;
			  for(g=0;g <= gapend; g++,i_gap--){
	                          if (rand_num<GAP[g]){rand_num=GAP[g];best_i=i_gap;}
			  }
			  g = (gapstart-best_i);
        	          for(l=0;l < g; l++) { back_operation[k++]='i';}
                	  score*=gapfunct[g];
			  i =best_i; 
//for( ;g > 0 ; g--) fprintf(stderr,"state = i; i=%d; j=%d\n",i+g-1,j);
                          rand_num=0;
                          if (rand_num<DELO[i][j]) {state='D'; rand_num=DELO[i][j];}
                          if (rand_num<DELE[i][j]) {state='d'; rand_num=DELE[i][j];}
                          if (rand_num<MAT[i][j])  {state='M'; }
                          break;

                default:  print_error("this should not happen");
	   }
#if 0
           fprintf(stderr,"state = %c; i=%d; j=%d; (s=%d)",state,i,j,s);
           fprintf(stderr,"score= %d; ", (Int4) (floor(pernats*log(score) +0.5)));
           fprintf(stderr,"back_operation = %c\n", back_operation[k-1]);
#endif
	}
        *alignscore=(Int4) (floor(pernats*log(score) +0.5));
//	fprintf(stderr,"score = %d\n%s\n",*alignscore,back_operation);
                          
        for(i=0;i < k;i++) operation[i]=back_operation[k-i-1];
//	fprintf(stderr,"operation = %s\n",operation);
        *trace_length=k;            
        free(GAP);free(back_operation);
        return operation;
}

void	sba_typ::init()
{ 
	weight=NULL;
}

void	sba_typ::alloc_weights()
//	ins 0..20; del 0..20.
{
	NEWP(weight,dim_weights+2,double);
	for(Int4 i=0; i<=dim_weights; i++) NEW(weight[i],dim_weights+2,double);
}

#if 0
double	sba_typ::indel_weight_sba(Int4 ins, Int4 del, double alignscore)
{
	double	tmp=alignscore;

	print_error("not yet implemented");
	if(weight == NULL) alloc_weights();
	assert(prof_len >= del);
	tmp-=lgamma(prof_len-nmbr_of_blocks+ins+1);
	tmp+=lgamma(ins+1);
	tmp+=lgamma(prof_len-nmbr_of_blocks+1);
	tmp-=lgamma(prof_len+1);
	if(tmp >= HUGE_VAL) print_error("indel_weight_sba() overflow problem");
	tmp+=lgamma(del+1);
	tmp+=lgamma(prof_len-del+1);
	tmp-=lgamma(seq_len-ins+del-prof_len+nmbr_of_blocks+1);
	if(tmp >= HUGE_VAL) print_error("indel_weight_sba() overflow problem");
	tmp+=lgamma(seq_len-ins+del-prof_len+1);
	tmp+=lgamma(nmbr_of_blocks+1);
	if(tmp >= HUGE_VAL) print_error("indel_weight_sba() overflow problem");
	return tmp;
}
#endif

double	sba_typ::calc_indel_weight_sba(Int4 ins, Int4 del)
{
	double	tmp=0.0;

	assert(prof_len >= del);
	tmp-=lgamma(prof_len-nmbr_of_blocks+ins+1);
	tmp+=lgamma(ins+1);
	tmp+=lgamma(prof_len-nmbr_of_blocks+1);
	tmp-=lgamma(prof_len+1);
	tmp+=lgamma(del+1);
	tmp+=lgamma(prof_len-del+1);
	if(tmp >= HUGE_VAL) print_error("indel_weight_sba() overflow problem");
	return tmp;
}

double	sba_typ::indel_weight_sba(Int4 ins, Int4 del, double score)
{
	double	tmp;

	if(weight == NULL) alloc_weights();
	if(ins <= dim_weights && del <= dim_weights){
	   if(weight[ins][del] > 0.0) tmp=weight[ins][del];
	   else tmp=weight[ins][del]=calc_indel_weight_sba(ins,del);
	} else tmp=calc_indel_weight_sba(ins,del);
	return score+tmp;
}

char    *sba_typ::TraceBack(Int4 *start, Int4 *alnscore, Int4 num_samples,
	double pernats)
{
	Int4	*score,trace_len,*strt;
	char	*operation,**seq_set;
	double	rand_num,total,*like,wt,tmp;
	Int4	i,j,ins=0, del=0,rtn;
// h_type	H=Histogram("number of indels",-2000,5000,1.0);

// for(i=1; i<=2000; i++) { std::cerr << lgamma(i); std::cerr << std::endl;}

	NEWP(seq_set,num_samples+2,char);
	NEW(like,num_samples+2,double);
	NEW(strt,num_samples+2,Int4);
	NEW(score,num_samples+2,Int4);
	for(total=0.0,i=1;i<=num_samples;i++){
		seq_set[i]=TraceBack(&score[i], &trace_len, &strt[i], pernats);
		for(del=0,ins=0,j=1;j<=trace_len;j++){
			if ((seq_set[i][j]=='d')||(seq_set[i][j]=='D')) del++;
			if (seq_set[i][j]=='I') ins++;
		}
// IncdHist((ins+del),H);
		tmp=score[i]/pernats;
#if 0
std::cerr << prof_len; std::cerr << std::endl;
std::cerr << nmbr_of_blocks; std::cerr << std::endl;
std::cerr << ins; std::cerr << std::endl;
std::cerr << del; std::cerr << std::endl;
#endif
		wt=indel_weight_sba(ins,del,tmp);
		// wt=indel_weight_sba(ins,del,seq_len,tmp);
		total += like[i] = exp(wt);
	}
// PutHist(stdout,60,H); NilHist(H);
	if(total >= HUGE_VAL) print_error("SampleASeq() overflow problem");
	rand_num=total*SampleUniformProb();
	for(i=1;i<=num_samples;i++){
		if((rand_num-=like[i]) <= 0.0) { rtn=i; break; }
	}
	operation=seq_set[rtn]; *start=strt[rtn]; *alnscore=score[rtn];
	for(j=1;j<=num_samples;j++) if(j!= rtn) free(seq_set[j]); 
	free(seq_set); free(like); free(strt);
	return operation;


}

Int4	sba_typ::CountIndels(char *trace, Int4 trace_len, Int4 *ins, Int4 *del)
{
	Int4 j,i,d; 

	for(d=0,i=0,j=1;j<=trace_len;j++){
		if((trace[j]=='d')||(trace[j]=='D')) d++;
		if(trace[j]=='I') i++;
	} *ins = i; *del = d;
	return i+d;
}


char    *sba_typ::TraceBack(Int4 *start, Int4 *alnscore, Int4 num_samples,
	double pernats, double Temperature)
{
	Int4	*score,trace_len,*strt;
	char	*operation,**seq_set;
	double	rand_num,total,*like,wt,sum;
	Int4	i,j,ins=0, del=0,rtn;
// h_type	H=Histogram("smatrix scores",-2000,5000,5.0);

	if(Temperature==1.0) return TraceBack(start,alnscore,num_samples,pernats);
	NEWP(seq_set,num_samples+2,char); NEW(like,num_samples+2,double);
	NEW(strt,num_samples+2,Int4); NEW(score,num_samples+2,Int4);
	for(total=0.0,i=1;i<=num_samples;i++){
		// seq_set[i]=TraceBack(&score[i], &trace_len, &strt[i], pernats,Temperature);
		seq_set[i]=TraceBack(&score[i], &trace_len, &strt[i], pernats);
		CountIndels(seq_set[i],trace_len,&ins,&del);
		wt=indel_weight_sba(ins,del,score[i]/pernats);
		total += like[i] = exp(wt);
// IncdHist(log(like[i]),H);
	}
// PutHist(stdout,60,H); NilHist(H);
	if(total >= HUGE_VAL) print_error("SampleASeq() overflow problem");
	for(sum=0.0,i=1;i<=num_samples;i++){
		sum += like[i] = pow(like[i]/total,Temperature);
		//sum += like[i] = pow(like[i]/total,1.0);
	}
// std::cerr << sum; std::cerr << std::endl; std::cerr << Temperature; std::cerr << std::endl;
	rand_num=sum*SampleUniformProb();
	for(i=1;i<=num_samples;i++){
		if((rand_num-=like[i]) <= 0.0) { rtn=i; break; }
	}
	operation=seq_set[rtn]; *start=strt[rtn]; *alnscore=score[rtn];
	for(j=1;j<=num_samples;j++) if(j!= rtn) free(seq_set[j]); 
	free(seq_set); free(like); free(strt);
	return operation;
}

