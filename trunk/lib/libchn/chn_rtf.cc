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

#include "chn_typ.h"

void	chn_typ::CreateRTF( )
//******************* create rich text format objects *************************
//******************* create rich text format objects *************************
{
  Int4 length = LengthCMSA(1,CMA[1]);
  Int4 a,anal;
  if(NumAnalysis <= 1){ // then show consensus alignment only.
        assert(!"This option not yet implemented");
  } else {
   //*********** a==0: display seqs = FG; Std Freq == BG. ************
   rtf[0]=new rtf_typ(FractHighlightTOP,fontsize,length,-cbp_cut,
			Obs[0],NullFreq[0],minbars,
			(Int4)0,AB,(Int4)0,MinResFreqCutoff,sets_mode);
   rtfQ[0]=new rtf_typ(FractHighlightTOP,fontsize,length,-cbp_cut,
			Obs[0], NullFreq[0],minbars,
			keyE,AB,0,MinResFreqCutoff,sets_mode);
   rtf[0]->SetExpPttrns(ExpPatterns);
   rtfQ[0]->SetExpPttrns(ExpPatterns);
   Int4	color_code;
   ColorCode[0]=0;
   if(verbose) fprintf(stderr,"Computing sequence constraints.\n");
   //******************* Begin of rtf_typ Main Loop *************************
   Int4 position0=0;
   if(NumResidues==1) position0=Position0[NumResidues];
   for(a=1,anal=1; anal <= NumberMCMA; anal++){
    if(anal == 1){	//****** a==1: display seqs == FG; next highest level == BG. ******
      if(!use_sfbg){
        rtf[a] = new rtf_typ(FractHighlight,fontsize,length,-cbp_cut,
			Obs[anal-1], WtFrq[anal],
			minbars,0,AB,position0,MinResFreqCutoff,sets_mode);
        rtfQ[a] = new rtf_typ(FractHighlight,fontsize,length,-cbp_cut,Obs[anal-1],
			WtFrq[anal],minbars,keyE,AB,position0,MinResFreqCutoff,sets_mode);
      } else {		// Test alternative background....
        rtf[a] = new rtf_typ(FractHighlight,fontsize,length,-cbp_cut,
			Obs[anal-1], WtFrq[NumberMCMA],
			minbars,0,AB,position0,MinResFreqCutoff,sets_mode);
        rtfQ[a] = new rtf_typ(FractHighlight,fontsize,length,-cbp_cut,Obs[anal-1],
			WtFrq[NumberMCMA],minbars,keyE,AB,position0,
				MinResFreqCutoff,sets_mode);
      }
      rtf[a]->SetExpPttrns(ExpPatterns);
      rtfQ[a]->SetExpPttrns(ExpPatterns);
      BackGrnd[a]=a; ColorCode[a]=1; color_code=2;
      Hist[a] = TRUE; FractSeqAln[a]=0; NumSqMain[a]=0; 
      Freq[a]=WtFrq[anal-1]; SuperAln[a]=FALSE;
      a++;  // increments a to be one larger than anal...
    }
    if(anal < NumberMCMA){	//***** a==x: Actual BPPS FG and BG (next highest level) *****
      if(ModeLPR == 'U'){	// called by mcBPPS routines...
	assert(cmc_res_vals && cmc_res_vals[anal]);
        rtf[a] = new rtf_typ(FractHighlight,fontsize,length,-cbp_cut,
				Obs[anal], WtFrq[anal+1],
				minbars,0,AB,position0,MinResFreqCutoff, sets_mode,
				Obs[anal+1], WtFrq[anal],	// new info: ObsBG, freqFG
				cmc_res_vals[anal]);
        rtfQ[a] = new rtf_typ(FractHighlight,fontsize,length,-cbp_cut,
				Obs[anal], WtFrq[anal+1],
				minbars,keyE,AB,position0,MinResFreqCutoff, sets_mode,
				Obs[anal+1], WtFrq[anal],	// new info: ObsBG, freqFG
				cmc_res_vals[anal]);
      } else if(ModeLPR == 'P'){	// Split these alternative IWtSq here...
	Alpha = Alphas[anal]; A0=A0s[anal]; B0=B0s[anal]; sets_mode = SetMode[anal];
	if(verbose) fprintf(stderr,"%d: Alpha=%.4f; a0=%d; b0=%d; sets_mode=%c\n",anal,Alpha,A0,B0,sets_mode);
	// This is where I need to pass in the Intergerized weights for BPPS from the cma files.
        rtf[a] = new rtf_typ(FractHighlight,fontsize,length,-cbp_cut,
				Obs[anal], WtFrq[anal+1],
				minbars,0,AB,position0,MinResFreqCutoff, sets_mode,
				Obs[anal+1], WtFrq[anal],	// new info: ObsBG, freqFG
				ModeLPR,Alpha,A0,B0);
        rtfQ[a] = new rtf_typ(FractHighlight,fontsize,length,-cbp_cut,
				Obs[anal], WtFrq[anal+1],
				minbars,keyE,AB,position0,MinResFreqCutoff, sets_mode,
				Obs[anal+1], WtFrq[anal],	// new info: ObsBG, freqFG
				ModeLPR,Alpha,A0,B0);
      } else if(ModeLPR == 'B'){
        rtf[a] = new rtf_typ(FractHighlight,fontsize,length,-cbp_cut,
				Obs[anal],WtFrq[anal+1],
				minbars,0,AB,position0,MinResFreqCutoff,sets_mode);
        rtfQ[a] = new rtf_typ(FractHighlight,fontsize,length,-cbp_cut,
				Obs[anal], WtFrq[anal+1],
				minbars,keyE,AB,position0,MinResFreqCutoff,sets_mode);
      } else print_error("ModeLPR invalid");
      rtf[a]->SetExpPttrns(ExpPatterns);
      rtfQ[a]->SetExpPttrns(ExpPatterns);
      BackGrnd[a]=a;	// which Analysis are observed counts coming from? 

      char	*AlignName=NameCMSA(MCMA[anal]);
      if(strcmp(AlignName,"BackGround")!=0){ ColorCode[a]=color_code; color_code++; }
      else ColorCode[a]=0;

      Hist[a] = TRUE; 
      FractSeqAln[a]=swt[anal]->FractSeqAln(); 
      NumSqMain[a]=NumSeqsCMSA(MCMA[anal]);
      Freq[a]=WtFrq[anal];
      SuperAln[a]=TRUE;	
      a++; 
    } else {	// This is the Main Set...or motif model from cma file...
	double **MainNullFreq=0;
        if(BG_CMA && use_bg_cma) MainNullFreq=NullFreqBG; else MainNullFreq=NullFreq[anal];
#if 0
        rtf[a] = new rtf_typ(FractHighlight,fontsize,length,-cbp_cut,Obs[anal],
				MainNullFreq,minbars,0,AB,position0,MinResFreqCutoff,sets_mode);
        rtfQ[a] = new rtf_typ(FractHighlight,fontsize,length,-cbp_cut,Obs[anal],
				MainNullFreq,minbars,keyE,AB,position0,
					MinResFreqCutoff,sets_mode);
#else	// make root node to use different residue pattern sets...(from calling enviroment.
        rtf[a] = new rtf_typ(FractHighlight,fontsize,length,-cbp_cut,Obs[anal],MainNullFreq,
				minbars,0,AB,position0,MinResFreqCutoff,root_sets_mode);
        rtfQ[a] = new rtf_typ(FractHighlight,fontsize,length,-cbp_cut,Obs[anal],
				MainNullFreq,minbars,keyE,AB,position0,
					MinResFreqCutoff,root_sets_mode);
#endif
        rtf[a]->SetExpPttrns(ExpPatterns);
        rtfQ[a]->SetExpPttrns(ExpPatterns);
	BackGrnd[a]=a; // set to higher number if cma input contains background sets...
	// if(use_sfbg){ BackGrnd[1]=a-1; }	// Set Query subfamily to SF background.
        ColorCode[a]=color_code; color_code++;
	Hist[a] = TRUE; 
	FractSeqAln[a]=swt[anal]->FractSeqAln(); 
	NumSqMain[a]=NumSeqsCMSA(MCMA[anal]);
	Freq[a]= WtFrq[anal]; 
	SuperAln[a]=TRUE; 
        a++; 
      if(BG_CMA){ // then create an addition rtf for motifs...
        rtf[a] = new rtf_typ(FractHighlightBG,fontsize,length,-cbp_cut,ObsBG,
				NullFreq[anal],minbars,0,AB,position0,
				MinResFreqCutoff,sets_mode);
        rtfQ[a] = new rtf_typ(FractHighlightBG,fontsize,length,-cbp_cut,ObsBG,
				NullFreq[anal],minbars,keyE,AB,position0,
				MinResFreqCutoff,sets_mode);
        rtf[a]->SetExpPttrns(ExpPatterns);
        rtfQ[a]->SetExpPttrns(ExpPatterns);
	BackGrnd[a]=a;
        // ColorCode[a]=color_code; color_code++;
        ColorCode[a]=0; // black == default.
	Hist[a] = TRUE; 
	FractSeqAln[a]=FractSeqAlnBG; 
	NumSqMain[a]=NumSeqsCMSA(BG_CMA);
	Freq[a]= WtFrqBG; 
	SuperAln[a]=TRUE; 
        a++; 
      }
    }
   }
  } //******************* End of rtf_typ Main Loop *************************
}


