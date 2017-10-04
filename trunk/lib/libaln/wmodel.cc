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

#include "wmodel.h"

wm_type	WModel(Int4 length, double pseudo, double *freq, a_type A)
{ return MkWModel(NULL,length, pseudo, freq, A); }

wm_type	MkWModel(char *null, Int4 length, double pseudo, double *freq, 
	a_type A)
/* create and return a wmodel with totsites = 0 and unfragmented columns. */
{
	wm_type	M;
	Int4	j,r;
	double	total;
	
	NEW(M,1,wmodel_type);
	M->A = A; M->D = NULL; M->length = length;
	M->method = 'm';  /** default is modified Gribskov's  method **/
	NEW(M->null, length+2, char);
	if(null != NULL){ for(j=1; j<=length; j++) M->null[j] = null[j]; }
	else { for(j=1; j<=length; j++) M->null[j] = '*'; }
// if(null!=NULL) printf("%s\n",null+1);
	M->pseudo = MAXIMUM(double,pseudo,0.0);
	M->pseudo = MINIMUM(double,M->pseudo,10.0);
	M->totsites = 0;
	NEW(M->N0,nAlpha(M->A)+2,double);
	NEW(M->temp,nAlpha(M->A)+2,double);
	NEW(M->freq,nAlpha(A)+2,double);
	for(total=0.0, r=0; r<= nAlpha(A); r++) total+= freq[r]; 
	for(r=0; r<= nAlpha(A); r++) M->freq[r] = freq[r]/total;
	NEWP(M->site_freq, M->length+2, double);
	NEWP(M->likelihood, M->length+2, double);
	NEW(M->tmp_val, M->length+2, double);
	for(j=0; j <= M->length; j++){	
	    NEW(M->site_freq[j],nAlpha(M->A) +2, double);
	    NEW(M->likelihood[j], nAlpha(M->A)+2, double);
	}
	M->smx = NULL;
	M->update = TRUE;
	return M;
}

void	InitWModel(wm_type M)
/* Initialize model M to no sites */
{
	Int4	j,b;

	M->totsites = 0.0;
	for( j = 1; j <= M->length; j++){	
	   for(b=0; b <= nAlpha(M->A); b++){ M->site_freq[j][b] = 0.0; }
	}
	M->update = TRUE;
}

wm_type NilWModel(wm_type M)
/* Destroy model M */
{
   Int4 j;
   if(M!=NULL){ 
	free(M->freq); free(M->N0); free(M->temp);
	for(j = 0; j <= M->length; j++){
	    free(M->site_freq[j]); 
	    free(M->likelihood[j]);
	}
	free(M->site_freq); 
	free(M->likelihood); free(M->tmp_val);
	if(M->smx != NULL) NilSMatrix(M->smx);
	if(M->D != NULL) NilDMPriors(M->D);
	free(M->null);
	free(M);
   }
   return (wm_type) NULL;
}

Int4	MaxSegWModel(unsigned char *maxseg, wm_type M)
{
	if(M->update || M->smx == NULL) get_smx_wmodel(M);
	return MaxSegSMatrix(maxseg, M->smx);
}

smx_typ	GetSMatrixWModel(wm_type M)
/* Return the score matrix for model M; score is in half bits */
{
	if(M->update || M->smx == NULL) get_smx_wmodel(M);
	return M->smx;
}

double	ExpectedScoreWModel(BooLean *use, wm_type M)
/* Return the expected log likelihood score for a member sequence to model M */
/* score is in half bits */
{
	Int4	i,j,s,b;
	double	score,factor,q;

	if(M->update || M->smx == NULL) get_smx_wmodel(M);
	factor = M->totsites + M->npseudo;
	for(score=0.0,j=1; j <= M->length; j++){
	    if(use == NULL || use[j]){
              for(b=1; b <= nAlpha(M->A); b++){
if(M->site_freq[j][b] > 0.0){
		q = M->site_freq[j][b]/M->totsites;
#if 0
		q = (M->site_freq[j][b]+M->N0[b])/factor;
#endif
		s=ValSMatrix(j,b,M->smx);
		score += (double)s*q;
}
	      }
	    }
	}
	return score;
}

Int4	SubScoreWModel(unsigned char *seq, Int4 pos, Int4 start, Int4 end, wm_type M)
{
	if(start < 1 || end < start || end > LenWModel(M)) 
			print_error("input error in SubScoreWModel( )");
	if(M->update || M->smx == NULL) get_smx_wmodel(M);
	return SubScoreSMatrix(seq, pos, start, end, M->smx);
}

Int4	ScoreWModel(register unsigned char *seq, register Int4 pos, register wm_type M)
/* Return the log likelihood score for site at pos in seq using model M */
/* score is in half bits */
{
	if(M->update || M->smx == NULL) get_smx_wmodel(M);
	return ScoreSMatrix(seq,pos,M->smx);
}

double	PvalWModel(register unsigned char *seq, register Int4 pos, register wm_type M)
/* Return the p-value for log likelihood score of site at pos in seq using 
   model M score is in half bits. */
{
        register Int4	score;
	if(M->update || M->smx == NULL) get_smx_wmodel(M);
	score = ScoreSMatrix(seq,pos,M->smx);
	/** fprintf(stderr,"method = %c; score = %d\n",M->method,score);/***/
        return SMatrixProb(score, M->smx);
}

double  InfoColWModel(Int4 col, wm_type M)
/* Report the information content of column "observed" relative to model M. */
{
        Int4    b;
        double  total,info,p,q,*site_freq=M->site_freq[col];

        for(total=0.0, b = 1; b <= nAlpha(M->A); b++){
                total += (double)site_freq[b] + M->N0[b];
        }
        for(info=0.0,b = 1; b <= nAlpha(M->A); b++){
                p = (site_freq[b]+M->N0[b])/total;
                if(p > 0.0){
                        /** q = M->freq[b]; /***/
                        q = 0.05; /***/
                        info += p*log(p/q); /* info in nats */
                }
        }
        return info;       /* info in nats */
	/** return 1.442695*info;       /* info in bits */
}

void	get_smx_wmodel(register wm_type M)
{
        Int4	i,score,c,ave,x,R;
	double	s,factor; /*factor only needed for regular model*/
	double	*prob,v,B,b,L,m=5.0,n,N;
	double	weight,aveinfo; 
	double	*freq = M->freq; /** blosum62freq **/
	double  NewBlosumL[21][21];
        BooLean debug=FALSE;


	if(M->update) update_wmodel_freqN(M);
	if(M->smx == NULL) {
           if(M->method == 'b'){
             for(c=0; debug && c<=nAlpha(M->A); c++){
                fprintf(stderr,
                   "blosum62freq[%c]/freq[%c] = %.4f/%.4f = %.4f\n",
                        AlphaChar(c,M->A), AlphaChar(c,M->A),
                        blosum62freq[c], freq[c],
                        blosum62freq[c]/freq[c]);
             }
             for(c=1; c<=nAlpha(M->A); c++){
                for(x = 1; x<=nAlpha(M->A); x++){
                        v = (double)blosum62L[c][x];
                        // v = (double)OurBlosumL[c][x];
                        v *= blosum62freq[c]/freq[c];
                        v *= blosum62freq[x]/freq[x];
                        NewBlosumL[c][x] = v;
                    // if(debug) fprintf(stderr,"NewBlosumL[%c][%c] = %.4f (%.4f)\n",
                    //     AlphaChar(c,M->A),AlphaChar(x,M->A),v, blosum62L[c][x]);
                }
             }
           }
	   M->smx = MkSMatrix(2.0,M->length,M->freq,M->A);
	   N = M->totsites;
	   switch (M->method) {
		case 'r': case 'd':
			factor=5.; 
	   	// 	factor=2.8853901/(M->maxscore/WMODEL_MAX_SCORE);
		  break;
                case 'b': factor=10.; break;
		case 'h': factor=5.; break;
		case 'f': case 'm': factor=10.; break;
		default: factor = 5.; break;
	   }
#if 0
	   fprintf(stderr,"DEBUG: maxscore = %g; factor = %g\n",
		M->maxscore, factor); 
#endif
           for(i=1; i<=M->length; i++){
	     if(M->null[i] == '!'){   // '!' = active site
            	for(ave=0,c=1; c<=nAlpha(M->A); c++){
#if 1
		     score=(Int4)floor(factor*log(M->likelihood[i][c])+0.5);
		     // score=(Int4)floor(5.0*log(M->likelihood[i][c])+0.5);
#endif
#if 0
               	    s= (double)M->site_freq[i][c]/(double)(M->totsites);
		    if(s <= 0.) s = (0.2/(double) (M->totsites));
		    s /= M->freq[c];
		    score=(Int4)floor(factor*log(s)+0.5);
#endif
		    SetSMatrix(c,i,score,M->smx);
		    ave += score;
		}
		// s = 0.1/(double) (M->totsites);
		// score=(Int4)floor(factor*log(s)+0.5);
                score = (Int4)floor(((double)ave/(double)nAlpha(M->A)) -factor + 0.5);
        	SetSMatrix(0,i,score,M->smx);
             } else if(M->null[i] != '.' && M->null[i] != '^'){
	       if(M->method == 'h' || M->method == 'b') {
                 for(R=0,c=1; c<=nAlpha(M->A); c++){
		    if(M->site_freq[i][c] > 0) R++;
		 }
		 B = m * (double)R;
	       }
               for(ave=0,c=1; c<=nAlpha(M->A); c++){
		if(M->N0[c] != 0){
		 switch (M->method) {
		  case 'r': /****** Regular Model *******/
		  case 'd': /****** Dirichlet Mixtures *******/
		     score=(Int4)floor(factor*log(M->likelihood[i][c])+0.5);
		    break;
                  case 'b': /*** MOD. GRIBSKOV with motif background freq ***/
                    for(s = 0.0,x = 1; x<=nAlpha(M->A); x++){
                      if(M->site_freq[i][x] > 0.0){
                        s += (M->site_freq[i][x]/N)*NewBlosumL[c][x];
                      }
                    } score = (Int4) floor((factor*log(s) + 0.5));
                    break;
		  case 'h': /****** HENIKOFF's METHOD *******/
		    n = M->site_freq[i][c];
		    if(n == 0.0){ L = 0.0; }
		    else { L = (n/N)/freq[c]; }
		    // for(b = 0.0,x = 0; x<=nAlpha(M->A); x++) // changed...
		    for(b=0.0,x=1; x<=nAlpha(M->A); x++)
		    {
			b += (M->site_freq[i][x]/N)* (double)blosum62L[c][x];
			// b += (M->site_freq[i][x]/N)* (double)OurBlosumL[c][x];
		    } s = (N*L + B*b)/(N+B); 
                    score = (Int4) floor((factor*log(s) + 0.5));
		    /*** in fifth nats... ***/
		    break;
		  case 'f': /****** MODIFIED GRIBSKOV with Blosum45 *******/
		  case 'm': /****** MODIFIED GRIBSKOV METHOD *******/
		    for(s = 0.0,x = 0; x<=nAlpha(M->A); x++){
			s += (M->site_freq[i][x]/M->totsites)*(double) blosum62L[c][x]; 
			// s += (M->site_freq[i][x]/M->totsites)*(double)OurBlosumL[c][x];
			// s += (M->site_freq[i][x]/M->totsites)*(double) blosum45L[c][x];
		    }
                    score = (Int4) floor((factor*log(s) + 0.5));
		    /*** score is in tenth nats ***/
		    break;
		  default: /**** GRIBSKOV'S PROFILE METHOD ****/
		    for(s = 0.0,x = 0; x<=nAlpha(M->A); x++){
			  s += M->site_freq[i][x]*blosum62[c][x];
		    }
                    score = (Int4) floor(factor*s + 0.5);
		 }
		 SetSMatrix(c,i,score,M->smx);
		 ave += score;
                } else SetSMatrix(c,i,0,M->smx);
	       }
               // score = (Int4) floor(((double)ave/(double)nAlpha(M->A)) + 0.5);
               score = (Int4) floor(((double)ave/(double)nAlpha(M->A)) -factor + 0.5);
               SetSMatrix(0,i,score,M->smx);
	     } else if(M->method == 'f'){
               for(ave=0,c=1; c<=nAlpha(M->A); c++){
                if(M->N0[c] != 0){
                    for(s = 0.0,x = 0; x<=nAlpha(M->A); x++){
                        s += ((double)M->site_freq[i][x]/(double)M->totsites)
                                        * (double) blosum45L[c][x]; 
                    }
                    score = (Int4) floor((factor*log(s) + 0.5));
                    SetSMatrix(c,i,score,M->smx);
                    ave += score;
                } else SetSMatrix(c,i,0,M->smx);
               }
               // score = (Int4) floor(((double)ave/(double)nAlpha(M->A)) + 0.5);
               score = (Int4) floor(((double)ave/(double)nAlpha(M->A)) -factor + 0.5);
               SetSMatrix(0,i,score,M->smx);
	     } 
           }
	} else wmodel_error("second update not yet implemented");
#if 0
	PutSMatrix(stdout, M->smx); PutWModel(stdout, M);
#endif
}

wm_type	MergeWModels(Int4 N, Int4 *p, Int4 leng, double *freq, wm_type *M)
/**********************************************************************
  Merge N models into one of leng with submodels at positions P.
 **********************************************************************/
{
	wm_type	M2;
	Int4	i,j,k,s,end,total,start,r;
	double	totsites,pseudo;
	char	*null;
	a_type	A;

	if(N < 1) return NULL;
	A = M[1]->A;
	totsites=M[1]->totsites;
	pseudo=M[1]->pseudo;
	for(total=end=0,i=1; i<= N; i++){  /** consistency check **/
		if(end >= p[i]) return NULL;
		if(M[i]->totsites != totsites) return NULL;
		if(M[i]->pseudo != pseudo) return NULL;
		end = p[i] + LenWModel(M[i]) - 1;
		if(end > leng) return NULL;
		total +=LenWModel(M[i]);
	}
	if(leng < total) return NULL;
	NEW(null,leng + 3, char);
	for(s=1,i=1; i<= N; i++){
		while(s < p[i]) null[s++] = '.';
		end = p[i] + LenWModel(M[i]) - 1;
		while(s <= end) null[s++] = '*';
	} while(s <= leng) null[s++] = '.';
	M2 = MkWModel(null, leng, pseudo, freq, A);
	for(i=1; i<= N; i++){
	   end = p[i] + LenWModel(M[i]) - 1;
	   for(s=p[i],j=1; j<=M[i]->length; j++,s++){
	   	if(M[i]->site_freq[j] == NULL)
			print_error("input error in MergeWModels( )");
		for(r=0; r <= nAlpha(A); r++){
			M2->site_freq[s][r] = M[i]->site_freq[j][r];
		}
	   }
	}
	M2->totsites=totsites;
        update_wmodel_freqN(M2);
	free(null);
	return M2;
}

double  **RealScoresWModel(wm_type M)
/** return a 2 dimensional array of real values log-odds scores **/
{
        Int4    i,c,x;
        double  s,factor; /*factor only needed for regular model*/
        double  ave,score,**scores;

        if(M->update) update_wmodel_freqN(M);
        /** printf("method = %c\n",M->method);/***/
        if(M->method == 'r' || M->method == 'a' || M->method == 'd') {
                factor=2.8853901/(M->maxscore/WMODEL_MAX_SCORE);
        }
        NEWP(scores,nAlpha(M->A)+2,double);
        for(c=0; c<=nAlpha(M->A); c++) NEW(scores[c],M->length+2,double);
        for(i=1; i<=M->length; i++){
             if(M->null[i] != '.'){
               for(ave=0,c=1; c<=nAlpha(M->A); c++){
                if(M->N0[c] != 0){
                 switch (M->method) {
                  case 'r': /****** Regular Model *******/
                  case 'd': /****** Dirichlet Mixtures *******/
                  case 'a': /****** asymptotic Model *******/
                     score=factor*log(M->likelihood[i][c]);
                    break;
		  case 'm': /****** MODIFIED GRIBSKOV METHOD *******/
                    for(s = 0.0,x = 1; x<=nAlpha(M->A); x++){
#if 1
                        s += ((double)M->site_freq[i][x]/(double)M->totsites)
                                        * (double) blosum62L[c][x]; 
#endif
#if 0
                        s += ((double)M->site_freq[i][x]/(double)M->totsites)
                                        * (double) OurBlosumL[c][x]; 
#endif
                    }
                    score = log(s);
                    break;
                  default: /**** GRIBSKOV'S PROFILE METHOD ****/
                    for(s = 0.0,x = 1; x<=nAlpha(M->A); x++){
                          s += (double)M->site_freq[i][x]*blosum62[c][x];
                    }
                    score = s;
                 }
                 scores[c][i] = score;
                 ave += score;
                } else scores[c][i] = 0.0;
               }
               /** scores[0][i] = 0.0; /****/
               scores[0][i] = (ave/(double) nAlpha(M->A)); 
             } else for(c=0; c<=nAlpha(M->A); c++) scores[c][i] = 0.0;
        }
        return scores;
}

Int4	CellScoreWModel(Int4 r, Int4 pos, wm_type M)
/* r = residue; pos = position */
{
	if(M->update || M->smx == NULL) get_smx_wmodel(M);
        if(pos > 0 && pos <= M->length && r >= 0 && r <= nAlpha(M->A)){
		return ValSMatrix(pos,r,M->smx); 
	} else wmodel_error("Referenced cell is out of bounds");
}

double	*ObservedWModel(Int4 j, wm_type M)
{
	if(M->update) update_wmodel_freqN(M);
	if(j < 1 || j > M->length) return NULL;
	else return M->site_freq[j];
}

void	Add2WModel(unsigned char *seq, Int4 site, double w, wm_type M)
/* Add the segment at site in seq to model with inteher weight w */
{
	Int4 j;
	for(j=1; j<=M->length; j++,site++)M->site_freq[j][seq[site]]+=w;
	M->totsites+=w; M->update = TRUE;
}

void	update_wmodel_freqN(register wm_type M)
/* Normalize frequencies to avoid overflow */
{
	register Int4	j,b,score,max;
	register double factor;
	register a_type	A=M->A;

	if(M->smx != NULL) { NilSMatrix(M->smx); M->smx = NULL; }
	/*** DETERMINE REGULAR AND NORMALIZED NONSITE FREQUENCIES ***/
	/** printf("method = %c\n",M->method);/***/
	M->npseudo = M->pseudo * sqrt(M->totsites);
	if(M->method == 'd' && M->D == NULL){
		M->D=MkDMPriors(TRUE);
	}
	for(b=1; b <= nAlpha(A); b++){
		M->N0[b]= M->npseudo * M->freq[b];
		M->likelihood[0][b] = M->freq[b];
	}
	/*** DETERMINE NORMALIZED SITE FREQUENCIES ***/
	factor = (double) (M->totsites + M->npseudo);
        for(j=1; j<=M->length; j++){
              M->likelihood[j][0] = 1.0;
	      if(M->method == 'd'){
		CalcDMPriors(M->site_freq[j], M->temp, M->D, M->A);
        	for(b=1; b <= nAlpha(A); b++){
               	    M->likelihood[j][b] = M->temp[b];
               	    M->likelihood[j][b] /= M->freq[b];
		}
	      } else {
        	for(b=1; b <= nAlpha(A); b++){
		  if(M->freq[b]==0.0){
		    M->likelihood[j][b] = 0.0;
		  } else {
               	    M->likelihood[j][b] = (M->site_freq[j][b]+M->N0[b])/factor;
               	    M->likelihood[j][b] /= M->likelihood[0][b];
#if 0
		    fprintf(stderr,"DEBUG: likelihood[%d][%c] = %f\n",
			j,AlphaChar(b,M->A), M->likelihood[j][b]);
#endif
		  }
           	}
	      }
        }
        for(M->maxscore=0.0,j=1; j<=M->length; j++){
                for(max=0,b=1; b<=nAlpha(M->A); b++){
		  if(M->freq[b] > 0.0){
		    score = 
		     (Int4) floor(2.8853901*log(M->likelihood[j][b])+0.5);
		    max = MAXIMUM(Int4,score,max);
		  }
		}
		M->maxscore += max;
	}
	M->update = FALSE;
}

void    wmodel_error(char *s) { fprintf(stderr,"wmodel: %s\n",s); exit(1); }

BooLean	NullSiteWModel(Int4 s,wm_type M)
/* if site s == null in model M return TRUE; else return FALSE. */
{
	if(s < 1 || s > M->length) return TRUE;
	else if(M->null[s] == '.') return TRUE; 
	else return FALSE;
}

double	PutWModel(FILE *fptr, wm_type M)
/* Report the current frequency model. */
{
	Int4	j,b,r;
	double	total,info,tot_info,p,q;
	BooLean flag = TRUE;

	if(M->update) update_wmodel_freqN(M);
	fprintf(fptr,"POS  ");
	for(tot_info=0.0, j=1; j<= M->length; j++){
	  if(M->null[j]!='.'){
	    for(total=0.0, b = 1; b <= nAlpha(M->A); b++){
	      total += M->site_freq[j][b] + M->N0[b];
	      if(flag) fprintf(fptr,"%3c", AlphaChar(b, M->A));
	    }
	    if(flag) { fprintf(fptr,"  Info (bits)\n"); flag = FALSE; }
	    fprintf(fptr,"%4d ",j);
	    for(info=0.0,b = 1; b <= nAlpha(M->A); b++){
	    	p = (M->site_freq[j][b]+M->N0[b])/total;
		if(p > 0.0){
			q = M->freq[b];
			info += p*log(p/q)/log(2.0);
		}
	    	fprintf(fptr,"%3d", (Int4)(100*p+0.5));
	    }
	    fprintf(fptr,"   %1.3lf\n", info);
	    tot_info += info;
	  }
	}
	fprintf(fptr,"non-\nsite ");	
	for(b = 1; b <= nAlpha(M->A); b++){
		r = (Int4) (100.0 * M->N0[b]/M->npseudo);
	    	fprintf(fptr,"%3d", r);
	}
	fprintf(fptr,"\n\tinformation = %g\n",tot_info);
	return tot_info;
}

unsigned char	GetBackGroundWModel(wm_type M)
{
        Int4    i,j,c;
        double  r,totLike;

	if(M->update) update_wmodel_freqN(M);
        r = (double) Random()/(double) RANDOM_MAX; /* 0 <= r <= 1 */
        for(c=1; c <= nAlpha(M->A); c++){
           if((r-=M->freq[c]) <= 0.0) return c; 
        }
	return 0;
}

Int4    GetSegWModel(unsigned char *seg, wm_type M)
/* Get a random segment drawn from the model M. Note: positions
   corresponding to null columns are drawn from the non-site model */
{
        Int4    i,j,c;
        double  r,T;

	if(M->update) update_wmodel_freqN(M);
	for(i=0, j=1; j<= M->length; j++){
            r = (double) Random()/(double) RANDOM_MAX; /* 0 <= r <= 1 */
            for(T=0.,c=1; c <= nAlpha(M->A); c++) T+=M->site_freq[j][c];
            for(c=1; c <= nAlpha(M->A); c++){
                    if((r-=M->site_freq[j][c]/T) <= 0.0){ seg[++i]=c; break; }
            }
        }
        return i;
}

