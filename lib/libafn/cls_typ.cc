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

#include "cls_typ.h"

static double charttwo[21][7] = {
       {0.000,0.000,0.000,0.000,0.000,0.000,0.000}, // 'X'
       {0.918,0.002,0.385,0.440,0.138,0.432,0.079}, // 'C'
       {0.084,0.215,0.432,0.111,0.153,0.367,0.125},
       {1.283,1.364,1.077,2.219,0.490,1.265,0.903},
       {0.332,0.753,0.930,0.424,0.734,0.801,0.518},
       {0.197,0.543,0.647,0.680,0.905,0.643,0.808},
       {1.231,1.683,2.157,0.197,1.653,2.430,2.065},
       {0.068,2.103,1.646,0.182,0.664,1.581,1.401},
       {0.281,3.351,2.998,0.789,4.868,2.735,3.812},
       {0.311,2.290,2.330,0.811,2.596,2.155,2.585},
       {1.233,2.194,1.817,0.611,2.095,1.686,2.027},
       {1.014,1.476,1.771,0.114,1.667,2.006,1.844},
       {0.590,0.646,0.584,0.842,0.307,0.611,0.396},
       {0.066,0.064,0.065,0.747,0.006,0.115,0.014},
       {1.319,0.064,0.081,1.526,0.204,0.118,0.096},
       {0.490,0.075,0.391,0.639,0.125,0.081,0.038},
       {1.525,0.479,0.350,0.887,0.286,0.350,0.362},
       {2.408,0.261,0.345,0.931,0.402,0.440,0.289},
       {2.998,0.269,0.367,3.852,0.510,0.514,0.562},
       {2.161,0.605,0.442,1.441,0.607,0.457,0.570},
       {0.004,0.108,0.018,0.006,0.010,0.004,0.007}};

static double chartone[21][7] = { 
       {0.000,0.000,0.000,0.000,0.000,0.000,0.000},
       {0.824,0.022,0.308,0.152,0.180,0.156,0.044},
       {0.045,0.275,0.578,0.216,0.211,0.426,0.156},
       {1.297,1.551,1.084,2.612,0.377,1.248,0.877},
       {0.382,0.583,1.052,0.419,0.525,0.916,0.628},
       {0.169,0.702,0.955,0.654,0.791,0.843,0.647},
       {0.835,1.475,1.534,0.039,1.722,2.456,2.280},
       {0.030,2.352,2.268,0.237,0.663,1.620,1.448},
       {0.262,3.496,3.108,0.998,5.685,2.494,3.048},
       {0.179,2.114,1.778,0.631,2.550,1.578,2.526},
       {1.375,2.639,1.763,0.191,1.815,1.961,2.795},
       {0.659,1.163,1.210,0.031,1.358,1.937,1.798},
       {0.347,0.275,0.679,0.395,0.294,0.579,0.213},
       {0.240,0.000,0.000,0.456,0.019,0.000,0.000},
       {1.417,0.090,0.122,1.659,0.190,0.130,0.155},
       {0.531,0.076,0.403,0.662,0.189,0.106,0.013},
       {1.665,0.403,0.386,0.949,0.211,0.342,0.360},
       {2.597,0.098,0.345,0.894,0.514,0.471,0.431},
       {3.167,0.297,0.398,3.902,0.585,0.501,0.483},
       {2.240,0.370,0.480,1.409,0.541,0.772,0.663},
       {0.000,0.008,0.000,0.013,0.000,0.000,0.000} };

cls_typ::cls_typ() { init = 0; }

cls_typ::cls_typ(a_type a, Int4 maxlength ) 
{
	char	default_ab[22] = "XCGASTNDEQKRHWYFVILMP";
	Int4	r,i,j;

	init = 0;
	A = a;
	chartoption = 2;

	TheChart = new double*[21];
	for(i=0; i < 21; i++ ){ TheChart[i] = new double[7]; }
	MaxLength=maxlength;
	init_memory(maxlength);
	for(i=0; i < 21; i++){
	  r = AlphaCode(default_ab[i],A);
	  for(j=0; j<7; j++) TheChart[r][j] = charttwo[i][j];
	} init = 1;
}

cls_typ::~cls_typ() {
  Int4  i;

  for ( i = 0 ; i < 21 ; i++ ) { delete [] TheChart[i]; }
  delete [] TheChart;
  delete [] heptnum;
  delete [] temphept;
  delete [] calcnumb;
  delete [] tempcalc;
}

BooLean	cls_typ::MaskCoils(e_type E, Int4 minseg, double cutoff)
{
    assert(LenSeq(E) <= MaxLength);
    assert(cutoff > 0.0 && cutoff <= 1.0);
    assert(minseg > 0);

    float	*f14, *f21, *f28;
    Int4	j,start=0,end=0;
    char	state;

    f14 = getArray(SeqPtr(E),LenSeq(E),'B');
    f21 = getArray(SeqPtr(E),LenSeq(E),'T');
    f28 = getArray(SeqPtr(E),LenSeq(E),'Q');
    for(state='u',start=0,j=1; j < LenSeq(E) ; j++ ) {
       if(f14[j] >= cutoff || f21[j] >= cutoff || f28[j] >= cutoff ) {
	  if(state == 'u'){
		state = 'm'; // masking...
		start= end = j;
	  } else {	// state == masking...
		end++;
	  }
       } else {
	  if(state == 'm'){ 
		if((end-start+1) >= minseg) MaskSeq(start, end, E);
		state = 'u'; 
	  }
       }
    } if(state == 'm'){ if((end-start+1) >= minseg) MaskSeq(start, end, E); }
    delete [] f14; delete [] f21; delete [] f28;
    if(end==0) return FALSE; else return TRUE;
}

float *cls_typ::getArray( unsigned char *seq, Int4 num, char opt )
{
  Int4 i;
  float *f;

  if(init){
    residue = seq;
    res_total = num;

    switch ( opt ) {
    case 'b': window = 14; weightchar = 'n'; break;
    case 'B': window = 14; weightchar = 'y'; break;
    case 't': window = 21; weightchar = 'n'; break;
    case 'T': window = 21; weightchar = 'y'; break;
    case 'q': window = 28; weightchar = 'n'; break;
    case 'Q': window = 28; weightchar = 'y'; break;
    default: print_error("cls_typ::getArray( ) input error");
    }
    f = new float[res_total+1];
    if(res_total < window) {
      for(i=0;i<=res_total;i++){ f[i] = 0.0; }
      return (f);
    }
    if((weightchar=='y') || (weightchar=='Y')) {
      // printf("weights: a,d=2.5 and b,c,e,f,g=1.0\n");
      hept_weight=10; ad_weight=2.5;
    } else {
      //  printf("no weights\n");
      hept_weight=7; ad_weight=1.0;
    }
    Calculate(window); Array_probs(f); return (f);
  } else { assert( init ); return (0); }
}

void cls_typ::init_memory(Int4 maxres)
{
   heptnum  = new Int4[maxres+1];
   temphept = new Int4[maxres+1];
   calcnumb = new double[maxres+1];
   tempcalc = new double[maxres+1];
}

void cls_typ::Calculate(Int4 startwindow)
  /* calculates best scores for each window frame */
// 98% of the time is spent in this routine...
{
  Int4 x,y,extras,heptad_pos,window_pos,res_pos;
  double root_inverse,scores,misc;
  Int4 hept,Window;
  double weight;
  Int4 root;

   // fprintf(stderr,"Calculating...\n");
   Window=startwindow; res_pos=0;
   root=(Window / 7) * hept_weight;
   root_inverse=1.0/(double) root;
   do {
    res_pos=res_pos+1;
    tempcalc[res_pos]=0;
    // go through each residue in each heptad pos.
    for(heptad_pos=0; heptad_pos<=6; heptad_pos++) {
      scores=1.0;
      hept=1;
      // get values at all 21 positions: innermost loop is here...
      for(window_pos=0;window_pos<=Window-1;window_pos++) {
        x = residue[window_pos + res_pos];
        y=(Int4) fmod((double)(window_pos + res_pos + heptad_pos),7.0);
	if(y==0) y=7;
        if(window_pos==0) hept=y;
	if((y==1) || (y==4)) weight=ad_weight; else weight=1.0;
	if(weight != 1.0) misc=pow(TheChart[x][y-1], weight);
	else misc=TheChart[x][y-1];
        scores=scores*(pow(misc,root_inverse));
// fprintf(stderr,"scores = %f; TheChart[%d][%d] = %f\n",scores,x,y-1,TheChart[x][y-1]);
      } /* end of window_pos loop */
      if(scores>tempcalc[res_pos]){
        tempcalc[res_pos]=scores;
        temphept[res_pos]=(Int4) fmod((double) (hept-1),7.0);
      }
    }  /* end of heptad_pos loop*/
   } while (res_pos+Window != res_total+1);
   for (extras=1; extras<=Window-1;extras++){
    tempcalc[res_pos+extras]=tempcalc[res_pos];
    temphept[res_pos+extras]=(Int4) fmod((double) (temphept[res_pos]+extras),7.0);
   }

   // maximize loop.
   res_pos=0;
   do {
    res_pos=res_pos+1;
    calcnumb[res_pos]=tempcalc[res_pos];
    heptnum[res_pos]=temphept[res_pos];
    window_pos=0;
    do {
      window_pos=window_pos+1;
      if(res_pos-window_pos<1) window_pos=Window-1;
      else if (tempcalc[res_pos-window_pos]>calcnumb[res_pos]) {
         calcnumb[res_pos]=tempcalc[res_pos-window_pos];
         heptnum[res_pos]=(Int4) fmod((temphept[res_pos-window_pos]+window_pos),7.0);
      } /* endif */
    } while (window_pos!=Window-1); /* enddo */
   } while (res_pos!=res_total); /* enddo */
}

double cls_typ::Calcprob(double x,double meancc,double stddevcc,double meangl,double stddevgl,double ratio_gl_cc)
{
  double prob1,prob2,prob3,prob4;

  prob1=(0.5) * pow(((x-meancc) / stddevcc),2.0);
  prob2=(0.5) * pow(((x-meangl) / stddevgl),2.0);
  prob3=stddevgl * exp(-prob1);
  prob4=ratio_gl_cc * stddevcc * exp(-prob2);
  return (prob3)/(prob3+prob4);
}

void cls_typ::Array_probs(float *f)
/* A similar and duplicated function from Row_probs  */
/* Row_probs was made for computing 3 windows probs, */ 
/* at the same time, a cost of code size for speed.  */
{
  Int4 res_pos;
  
  //printf("%s\n",nameprot);
  for (res_pos=1;res_pos<=res_total;res_pos++) {
      if (chartoption==2) {
	  if ((weightchar=='y') || (weightchar=='Y')) {
	      switch(window) {
	      case 14:
		f[res_pos]=(Calcprob(calcnumb[res_pos],1.89,0.30,1.04,0.27,20));
		break;
	      case 21:
		f[res_pos]=(Calcprob(calcnumb[res_pos],1.79,0.24,0.92,0.22,25));
		break;
	      case 28:
		f[res_pos]=(Calcprob(calcnumb[res_pos],1.74,0.20,0.86,0.18,30));
		break;
	      } /* end */
	  } else {
	      switch (window) {
	      case 14:
		f[res_pos]=(Calcprob(calcnumb[res_pos],1.82,0.28,0.95,0.26,20));
		break;
	      case 21:
		f[res_pos]=(Calcprob(calcnumb[res_pos],1.74,0.23,0.86,0.21,25));
		break;
	      case 28:
		f[res_pos]=(Calcprob(calcnumb[res_pos],1.69,0.18,0.80,0.18,30));
		break;
	      } /* end */
	    } /* end */
      } else {
	  if ((weightchar=='y') || (weightchar=='Y')) {
	      switch(window) {
	      case 28:
		if (calcnumb[res_pos] > 0)
		  f[res_pos]=(Calcprob(calcnumb[res_pos],1.70,0.24,0.79,0.23,30));
		else { f[res_pos]=0;  }
		break;
	      case 21:
		if (calcnumb[res_pos] > 0)
		  f[res_pos]=(Calcprob(calcnumb[res_pos],1.76,0.28,0.86,0.26,25));
		else { f[res_pos]=0;  }
		break;
	      case 14:
		if (calcnumb[res_pos] > 0)
		  f[res_pos]=(Calcprob(calcnumb[res_pos],1.88,0.34,1.00,0.33,20));
		else { f[res_pos]=0;  }
		break;
	      } /* end */
	  } else {
	      switch(window) {
	      case 28:
		if (calcnumb[res_pos] > 0)
		  f[res_pos]=(Calcprob(calcnumb[res_pos],1.628,0.243,0.770,0.202,30));
		else { f[res_pos]=0;  }
		break;
	      case 21:
		if (calcnumb[res_pos] > 0)
		  f[res_pos]=(Calcprob(calcnumb[res_pos],1.683,0.285,0.828,0.236,25));
		else { f[res_pos]=0;  }
		break;
	      case 14:      
		if (calcnumb[res_pos] > 0)
		  f[res_pos]=(Calcprob(calcnumb[res_pos],1.782,0.328,0.936,0.289,20));
		else { f[res_pos]=0;  }
		break;
	      }
	  } /* end */ 
	}
// fprintf(stderr,"res[%d] = %f\n",res_pos,f[res_pos]);
    }
}

BooLean	cls_typ::Put(FILE *fp, e_type E, double cutoff)
{
    assert(LenSeq(E) <= MaxLength);
    assert(cutoff > 0.0 && cutoff <= 1.0);

    float	*f14, *f21, *f28;
    Int4	j,start=0,end=0;
    BooLean	flag=FALSE;

    char *buffer = new char[LenSeq(E)+5];
    SeqToString(buffer+1,E,A); 
    printf("\n>%s\n", SeqKey(E));
    f14 = getArray(SeqPtr(E), LenSeq(E), 'B' );
    f21 = getArray(SeqPtr(E), LenSeq(E), 'T' );
    f28 = getArray(SeqPtr(E), LenSeq(E), 'Q' );

    // "Tree" output option, coiled region on left, reg. region on right 
    fprintf(fp,"\n");
    j = 1;
    Int4 prevJ = j;
    while(j <= LenSeq(E)) {
	if(f14[j] >= cutoff || f21[j] >= cutoff || f28[j] >= cutoff ) {
	  flag=TRUE; j++;
	  while(j <= LenSeq(E) &&
		(f14[j] >= cutoff || f21[j] >= cutoff || f28[j] >= cutoff)) j++;
	  printOutput(fp,buffer,prevJ,j-1,0,f14,f21,f28);
	} else {
	  while(j <= LenSeq(E) &&
		!( f14[j] >= cutoff || f21[j] >= cutoff || f28[j] >= cutoff)) j++;
	  printOutput(fp,buffer, prevJ, j-1, 1, NULL, NULL, NULL );
	} prevJ = j;
    } delete [] f14; delete [] f21; delete [] f28; delete [] buffer;
    return flag;
}

static const char * THIRTY = "                              ";
static const char * ELEVEN = "           ";

static float MaX ( float one, float two ) {
  return ( one > two?one:two );
}

void cls_typ::printOutput(FILE *fp, char *buffer,Int4 start,Int4 end,Int4 option,
	float *f14, float *f21,float *f28)
{
  short statHolder;
  Int4 i, j;

  if ( !option ) {
    if ( end - start <= 29 ) {
      for ( i = 1 ; i <= ( 29 - ( end - start ) ) ; i++ ) {
	fprintf(fp," ");
      }
      for ( i = start ; i <= end ; i++ ) {
	fprintf(fp, "%c", tolower(buffer[i]) ); 
      }
      fprintf(fp,"%5d-%-5d", start, end);
      for ( i = start ; i <= end ; i++ ) {
	statHolder = (short)floor(MaX(MaX(f28[i],f21[i]),f14[i])*10);
	fprintf(fp,"%hd", statHolder!=10?statHolder:9 );
      }
      fprintf(fp,"\n");
      return;
    } else {
      for ( i = start ; i <= start + 29 ; i++ ) {
	fprintf(fp, "%c", tolower(buffer[i]) );
      }
      fprintf(fp,"%5d-%-5d", start, end);
      for ( i = start ; i <= start + 29 ; i++ ) {
	statHolder = (short)floor(MaX(MaX(f28[i],f21[i]),f14[i])*10);
	fprintf(fp,"%hd", statHolder!=10?statHolder:9 );
      }
      fprintf(fp,"\n");
      for ( i = 30 + start ; i + 29 <= end ; i+=30 ) {
	for ( j = i ; j < i + 30 ; j++ ) {
	  fprintf(fp, "%c", tolower(buffer[j]) );
	}
	fprintf(fp, ELEVEN );
	for ( j = i ; j < i + 30; j++ ) {
	  statHolder = (short)floor(MaX(MaX(f28[j],f21[j]),f14[j])*10);
	  fprintf(fp,"%hd", statHolder!=10?statHolder:9 );
	}
	fprintf(fp,"\n");
      }
      if ( i <= end ) {
	for ( j = 0 ; j < 29 - ( end - i ) ; j++ ) {
	  fprintf(fp," ");
	}
	for ( j = i ; j <= end ; j++ ) {
	  fprintf(fp, "%c", tolower(buffer[j]) );
	}
	fprintf(fp, ELEVEN );
	for ( j = i ; j <= end ; j++ ) {
	  statHolder = (short)floor(MaX(MaX(f28[j],f21[j]),f14[j])*10);
	  fprintf(fp,"%hd", statHolder!=10?statHolder:9 );
	}
      } fprintf(fp,"\n");
    }
  } else {
    fprintf(fp, THIRTY );
    fprintf(fp,"%5d-%-5d", start, end);
    for ( i = start ; i <= end  && i <= 29 + start ; i++ ) { 
      fprintf(fp,"%c", buffer[i]);
    }
    if ( i <= end ) {
      fprintf(fp,"\n");
    }
    for ( i = start + 30 ; i + 29 <= end ; i += 30 ) {
      fprintf(fp, THIRTY );
      fprintf(fp, ELEVEN );

      for ( j = i ; j <= i + 29 ; j++ ) {
	fprintf(fp,"%c", buffer[j]);
      }
      fprintf(fp,"\n");
    }
    if( i <= end ) {
      fprintf(fp, THIRTY ); fprintf(fp, ELEVEN );
      for( ; i <= end ; i++ ) { fprintf(fp,"%c", buffer[i]); }
    } fprintf(fp,"\n");
  }
}


