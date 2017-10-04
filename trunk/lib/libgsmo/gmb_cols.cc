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

#include "gmb_typ.h"

#define Debug_FB 0

char	*gmb_typ::FindBlocks(ssx_typ *ssx, char mode, double cutoff, double bild_cut,Int4 Limit)
{
// ssx->InitNDL(0.0);     // initialize model parameters.
	cma_typ	cma=ssx->RtnCMA(),rcma=0,xcma=cma;
	assert(nBlksCMSA(cma) == 1);
	FILE *efp=0; // efp=stderr;
	Int4	blk=1,i,j,*N_INS,*INS_RES,n_ins,n_del,ins_res;
	char	*Report,c;
	double	*BILD,d,D,*deleted,BS,BSps,TotalWtSq=ssx->TotalWtSeq();
	enum location { start, middle, end, stop };
	BooLean	trigger;
	//=============== 1. Allocate memory. ================
	location *SITE; NEW(SITE,LengthCMSA(blk,cma) +3, location);
	NEW(Report,LengthCMSA(blk,cma) +3, char);
	NEW(BILD,LengthCMSA(blk,cma) +3, double);
	NEW(N_INS,LengthCMSA(blk,cma) +3, Int4);
	NEW(INS_RES,LengthCMSA(blk,cma) +3, Int4);
	//=============== 2. Compute data. ================
	double **xBS=BILD_ScoresCMSA(ssx),*BldSc=xBS[1]; free(xBS);
	for(SITE[0]=end, i=1; i <= LengthCMSA(blk,cma); i++){
		Report[i]=' '; 
#if 0
		BILD[i]=ssx->BildScore(blk,i);  
#else
		BILD[i]=ssx->ContextBildScore(blk,i,BldSc);
#endif
		n_ins=NumInsertsCMSA(blk,i,cma,ins_res); 
		N_INS[i]=n_ins; INS_RES[i]=ins_res;
#if Debug_FB
	   if(n_ins > 0){
		  d=(double)n_ins/(double)NumAlnSeqsCMSA(cma);
		fprintf(stderr,"%d. n_ins=%d; ins_res=%d (%.2f vs %.2f)\n",i,n_ins,ins_res,d,cutoff);
	   }
#endif
	        if(mode == 'A'){	// AddColumns() mode.
		  d=(double)n_ins/(double)NumAlnSeqsCMSA(cma);
		  // d=(double)n_ins/(double)NumSeqsCMSA(cma);
		  if(d >= cutoff) trigger=TRUE; else trigger=FALSE;
		} else if(mode =='R'){	// RmColumns() mode; treat this as start of a "block".
		   if(ins_res > 10) trigger=TRUE; else trigger=FALSE; 
		} else {		// RmShortBlk() mode.
		   if(ins_res > 10) trigger=TRUE; else trigger=FALSE; 
		}
		if(trigger){ SITE[i]=end; SITE[i+1]=start; }
                else if(SITE[i-1]==end) SITE[i]=start; 
		else SITE[i]=middle;
	} i=LengthCMSA(blk,cma); SITE[i]=end; SITE[i+1]=stop; Report[0]=Report[i+1]=' '; 
	free(BldSc);
	if(mode== 'A'){		// for AddColumns()
	  //=============== 3. Find columns next to insertions. ================
	  for(i=1; i <= LengthCMSA(blk,cma); i++){
	    if(SITE[i] == middle) continue;
	    switch (SITE[i]) {
		case start: Report[i]='B'; break;		// insert column before..
		case end: Report[i]='A'; break;		// insert column after...
		default: print_error("FindBlocks(): This should not happen");
	    }
	  }
	} else if(mode == 'R'){	// for RmColumns()
	  //=============== 3. Label bad columns at start or end of blocks. ================
	  Int4	*NDEL=0; 
	  if(cutoff > 0){
 	    NEW(NDEL, LengthCMSA(blk,cma)+9, Int4);
	    for(i=1; i <= LengthCMSA(blk,cma); i++) NDEL[i]=NumDeletionsCMSA(blk,i,cma);
	  }
#if 1	// avoid removing all columns...
	  double Cutoff=cutoff, Bild_cut=bild_cut,mfact=0.0,inc=0.05,dfact=1.0;
	  do {
	    Int4 num_cols=0,ndel=0,nbild=0,nboth=0;
            for(i=1; i <= LengthCMSA(blk,cma); i++){
	      d=(double)NDEL[i]/(double)NumAlnSeqsCMSA(cma);
	      if(BILD[i] >= bild_cut && d <= cutoff){ nboth++; num_cols++; }
	      if(BILD[i] >= bild_cut){ nbild++; }
	      if(d <= cutoff){ ndel++; }
	    } if(efp) fprintf(stderr,"%d columns retained; nbild=%d; ndel=%d; NumSeqs=%d (%.3g; %.3g; %d)\n",
				num_cols,nbild,ndel,NumAlnSeqsCMSA(cma),bild_cut,cutoff,Limit);
	    if(num_cols < Limit){	// Limit == 3 by default...
	      if(nbild < ndel){		// fewer BILD-signifcant columns than filled columns. 
		mfact = mfact + inc; // mfact = 0.05, 0.10, 0.15, 0.25...10.0
		d=fabs(bild_cut)*mfact;	bild_cut = Bild_cut - d;
	      } else { 			// fewer filled columns than BILD-signifcant columns. 
		dfact = dfact + inc; // mfact = 0.50, 0.55, 0.60, ... 0.75
		cutoff = Cutoff * dfact;	// e.g., cutoff = 0.50 ... 0.75
	      }
	      if(efp){
		// WriteCMSA("junk.cma",cma);
		fprintf(stderr,"%d high deletion columns retained\n",ndel);
		fprintf(stderr,"%d BILD scored columns retained\n",nbild);
		fprintf(stderr,"Too many columns removed; adjusting BILD parameter (%.3g)\n",bild_cut);
	      }
	    } else break;
	  } while(mfact < 10.0 && cutoff <= 0.75);
#endif
	  for(i=1; i <= LengthCMSA(blk,cma); i++){
	   if(BILD[i] >= bild_cut) continue;
           if(SITE[i] == middle) continue;
           switch (SITE[i]) {
                case start: do { Report[i]='*'; i++; } while(BILD[i] < bild_cut && SITE[i] != stop);
                  break;
                case end: j=i;
                  do { Report[j]='*'; j--; } while(j > 0 && BILD[j] < bild_cut && Report[j] == ' ');
                  break;
                default: print_error("FindBlocks(): This should not happen");
           }
	  }
#if 1	// if fraction of deletions > cutoff then remove the column.
	  if(cutoff > 0){
	    for(i=1; i <= LengthCMSA(blk,cma); i++){
		d=(double)NDEL[i]/(double)NumAlnSeqsCMSA(cma);
		if(d > cutoff) Report[i]='*';
	    }
	  } free(NDEL);
#endif
	} else {		// for RmShortBlk() ...
          //=============== 4. Label short blocks. ================
#if 0
          for(j=0,i=1; i <= LengthCMSA(blk,cma); i++){
                if(SITE[i] == start) j=0;
                if(Report[i]==' '){ j++;  continue; }
                if(Report[i]=='*'){
                   if(j < 3){ while(j > 0){ Report[i-j]='*'; j--; } } j=0;
                }
          }
#else	// New: remove all short blocks...
          for(j=0, i=1; i <= LengthCMSA(blk,cma); i++){
                if(SITE[i] == start) j=1;
                else if(SITE[i] == end){
                   if(j < 3){ while(j >= 0){ Report[i-j]='*'; j--; } } j=0;
		} else { assert(SITE[i] == middle); j++; }
          }
#endif
	}
	//=============== 4. Print out results. ================
	if(0){	// for testing...
	   this->RtnMap();
	   for(i=1; i <= LengthCMSA(blk,cma); i++){
		BS=BILD[i]; BSps=BS/TotalWtSq;
		n_del=NumDeletionsCMSA(blk,i,cma);
		n_ins=N_INS[i]; ins_res=INS_RES[i];
		switch (SITE[i]) {
		  case start: c='S'; break;
		  case middle: c=' '; break;
		  case end: c='E'; break;
		  default: print_error("FindBlocks(): This should not happen");
		} fprintf(stderr,"%c%c%d: bild=%.1f (%.2f npws); ins=%d (%d) del=%d (%.1f%c).\n",
			c,Report[i],i,BS,BSps,n_ins,ins_res,n_del,100.0*FractDeletionsCMSA(blk,i,cma),'%');
		if(ins_res > 10) ssx->PutPenalties(stderr,blk,i);
		else if(cutoff <= 0 && c==' ' && ins_res > 0) ssx->PutPenalties(stderr,blk,i);
	   } fprintf(stderr,"\n\n %s\n",Report+1);
	} free(BILD); free(N_INS); free(INS_RES); free(SITE);  
	return Report;
}

#define debug_AC 0

cma_typ	gmb_typ::AddColumns(double bild_cut, BooLean Widen, set_typ InSet,
			double MinMatchFrq, BooLean EndsOnly)
{
    Int4	blk=1,i,j,x,startA,endA,startB,endB,numcol,maxcols=4;
    cma_typ	cma=SSX->RtnCMA(),rcma=0,xcma=cma,bcma,acma;
    double	d,BS,bBS,aBS,prior_wt=SSX->GetPriorWt();
    BooLean	improved,found;

    str_typ	*head = new str_typ('E'),*Str=head,*tail;
    if(this->AddOp){ 
	x=strlen(this->AddOp); assert(x == (LengthCMSA(blk,cma)+2));
       	for(i=1; i <= LengthCMSA(blk,cma); i++){ Str->Append(AddOp[i]); Str=Str->Next(); }
    } else for(i=1; i <= LengthCMSA(blk,cma); i++){ Str->Append('m'); Str=Str->Next(); }
    Str->Append('E'); tail=Str->Next(); // head->Print(stderr);

    // if(EndsOnly) Widen=TRUE;
SSX->InitNDL(0.0);     // initialize model parameters.
  do {
    improved=FALSE;
    gmb_typ *xgmb= new gmb_typ(aa_per_io,aa_per_do,exp_ie,exp_de,xcma,dms_mode,prior_wt); 
    char    *BeforeOrAfter=this->FindBlocks(xgmb->RtnSSX(),'A',MinMatchFrq,bild_cut); delete xgmb;
    if(Widen){ endA=LengthCMSA(blk,xcma); startA=1; endB=endA; startB=1; }
    else { endA=LengthCMSA(blk,xcma)-1; startA=1; endB=endA+1; startB=startA+1; }
    for(Str=tail->Prev(),i=LengthCMSA(blk,xcma); i >= 0; Str=Str->Prev(),i--)
    // for(Str=tail->Prev(),i=105; i >= 100; Str=Str->Prev(),i--)
    {
        for(found=FALSE,numcol=1; numcol <= maxcols; numcol++){
           j=i+1; BS=-999999999;
#if debug_AC
    fprintf(stderr,"%d. Try adding %d cols: BOA[j=%d]=%c; BOA[i=%d]=%c.\n",
		i,numcol,j,BeforeOrAfter[j],i,BeforeOrAfter[i]);
#endif
	   if(BeforeOrAfter[j] == 'B' || BeforeOrAfter[i] == 'A'){
	     if(i > endA || i < startA){ acma=0; aBS=-99999999; }
	     else {	// Add column AFTER position i.
		if(i == LengthCMSA(blk,xcma)) acma=ExtendBlkCMSA(xcma, blk,0,numcol);
		else acma=InsertColumnsCMSA(xcma,blk,i,numcol);
		xgmb=new gmb_typ(aa_per_io,aa_per_do,exp_ie,exp_de,acma,dms_mode,prior_wt); 
		// xgmb->SetPriorWt(prior_wt);
		aBS=xgmb->ScoreBILD(blk,i+numcol); 
		ssx_typ *xssx=xgmb->RtnSSX(); d=xssx->RtnPercentGaps(blk,i+numcol);
#if debug_AC
		fprintf(stderr,"A %d: BS=%.2f(cut=%.2f)[%d] ",
			i+numcol,aBS,bild_cut,NumAlnSeqsCMSA(acma));
		xssx->PutWtCnts(stderr,blk,i+numcol);
#endif
		delete xgmb;
		if(d >= 50.0){ if(acma) NilCMSA(acma); acma=0; aBS = -9999999.0; }
	     }
	     if(j > endB || j < startB){ bcma=0;  bBS=-9999999; }
	     else {	// add column before position j.
		  if(j == 1) bcma=ExtendBlkCMSA(xcma, blk,numcol,0);
		  else bcma=InsertColumnsCMSA(xcma,blk,j,-numcol);
		  xgmb=new gmb_typ(aa_per_io,aa_per_do,exp_ie,exp_de,bcma,dms_mode,prior_wt); 
		  bBS=xgmb->ScoreBILD(blk,i+1); 
		  ssx_typ *xssx=xgmb->RtnSSX(); d=xssx->RtnPercentGaps(blk,i+1);
#if debug_AC
		  fprintf(stderr,"B %d: BS=%.2f(cut=%.2f)[%d] ",
				i+1,bBS,bild_cut,NumAlnSeqsCMSA(bcma));
		  xssx=xgmb->RtnSSX(); xssx->PutWtCnts(stderr,blk,i+1);
#endif
		  delete xgmb;
		  if(d >= 50.0){ if(bcma) NilCMSA(bcma); bcma=0; bBS = -9999999.0; }
	     }
	     if(aBS > bBS){ if(bcma) NilCMSA(bcma); bcma=0; BS=aBS; rcma=acma; }
	     else { if(acma) NilCMSA(acma); acma=0; BS=bBS; rcma=bcma; }
	     if(rcma==0){ numcol=maxcols; continue; }	// Need not lengthen any further.
	     if(BS < bild_cut){ if(rcma) NilCMSA(rcma); rcma=0;  continue; }
	     else { 
		// if(bcma) for(x=1; x <= numcol; x++){ Str=Str->Prev(); Str->Append('d'); } else
		for(x=1; x <= numcol; x++){ Str->Append('d'); }
		found=TRUE;
	     }
	   }
	   if(found){
#if debug_AC
	      if(bcma) fprintf(stderr,"Added %d columns before position %d\n",numcol,j);
	      if(acma) fprintf(stderr,"Added %d column after position %d\n",numcol,i);
#endif
	      if(xcma != cma) NilCMSA(xcma); xcma=rcma; improved=TRUE; 
	      Str=Str->Next(); i++; break; // break from numcol loop...
	   }
	} if(found) break;		// break from i loop.
    } free(BeforeOrAfter);
   } while(improved);
// head->Print(); 
if(AddOp) free(AddOp); AddOp=head->Return(); AddOp[1]=toupper(AddOp[1]);
// fprintf(stderr,"%s\n",AddOp); // exit(1);
delete head;
   if(xcma != cma) return xcma; else return 0;
}

cma_typ	gmb_typ::RmColumns(double bild_cut,BooLean Shrink,double MinDelFrq)
{
	SSX->InitNDL(0.0);     // initialize model parameters.
	cma_typ	cma=SSX->RtnCMA(),rcma=0,xcma=cma;
	Int4	limit=3,i,j,blk=1,start,end,Len=LengthCMSA(blk,xcma);
	char	c,*ToRm=this->FindBlocks(this->RtnSSX(),'R',MinDelFrq,bild_cut,limit);
#if 0
    str_typ	*head = new str_typ('E'),*Str=head,*tail;
    assert(this->AddOp == 0);	// free up before using this!!
    if(this->AddOp){ 
	x=strlen(this->AddOp); assert(x == (LengthCMSA(blk,cma)+2));
       	for(i=1; i <= LengthCMSA(blk,cma); i++){ Str->Append(AddOp[i]); Str=Str->Next(); }
    } else for(i=1; i <= LengthCMSA(blk,cma); i++){ Str->Append('m'); Str=Str->Next(); }
    Str->Append('E'); tail=Str->Next(); // head->Print(stderr);
#else
    assert(this->AddOp == 0);	// free up before using this!!
#endif
	if(Shrink){ end=LengthCMSA(blk,xcma); start=0; }
	else { end=LengthCMSA(blk,xcma)-1; start=1; }
	for(i=end; i > start; i--){
	   if(ToRm[i] != '*') continue;
	   for(j=i; ToRm[j-1] == '*'; ) j--;
	   if(i == LengthCMSA(blk,xcma)){
// assert(LengthCMSA(blk,xcma) >= ((i-j+1)+limit)); 
	   	rcma=TrimBlkCMSA(xcma,blk,0,i-j+1,limit); i=j;  // Blk=1; RmLeft=0; RmRight=i-j+1.
		if(rcma==0) print_error("gmb_typ::RmColumns(): too many aligned columns removed");
	   } else if(j==1){
// fprintf(stderr,"'%s'\n",ToRm);
// assert(LengthCMSA(blk,xcma) >= (i+limit)); 
		rcma=TrimBlkCMSA(xcma,blk,i,0,limit); i=j;	// Blk=1; RmLeft=i; RmRight=0.
		if(rcma==0) print_error("gmb_typ::RmColumns(): too many columns removed from alignment");
	   } else { rcma=ConvertColsToInsertsCMSA(xcma,1,j,i); i=j; }
	   if(xcma != cma) NilCMSA(xcma); xcma=rcma;
	}
	if(!Shrink){ 	// this->AddOp=rtn...
	  char	*rtn; NEW(rtn, Len+5, char); this->AddOp=rtn; rtn[0]='E';  rtn[1]='M';
	  for(j=2; j <= end && (c=ToRm[j]); j++){
	     if(c=='*') rtn[j]='d'; else rtn[j]='m'; 
	  } rtn[j]='m'; j++; rtn[j]='E'; j++; rtn[j]=0;   // assume extra length added by FindBlocks();
	} free(ToRm);
	assert(NumColumnsCMSA(xcma) >= 3);
	if(xcma != cma) return xcma; else return 0;
}

cma_typ	gmb_typ::RmShortBlks(double bild_cut,BooLean Shrink)
{
// return 0;
SSX->InitNDL(0.0);     // initialize model parameters.
	cma_typ	cma=SSX->RtnCMA(),rcma=0,xcma=cma;
	char	*ToRm=this->FindBlocks(this->RtnSSX(),'S',-1.0,bild_cut);
	Int4	i,j,blk=1,start,end;
	double	BS;
	if(Shrink){ end=LengthCMSA(blk,xcma); start=0; }
	else { end=LengthCMSA(blk,xcma)-1; start=1; }
	for(i=end; i > start; i--){
	   if(ToRm[i] != '*') continue;
	   for(j=i; ToRm[j-1] == '*'; ) j--;
	   if(i == LengthCMSA(blk,xcma)){
	   	rcma=TrimBlkCMSA(xcma,1,0,i-j+1,1); i=j;  // Blk=1; RmLeft=0; RmRight=i-j+1; limit=1.
		if(rcma==0) print_error("TrimBlkCMSA(): too many columns removed from alignment");
	   } else if(j==1){
		rcma=TrimBlkCMSA(xcma,1,i,0,1); i=j;	// Blk=1; RmLeft=i; RmRight=0; limit=1.
		if(rcma==0) print_error("TrimBlkCMSA(): too many columns removed from alignment");
	   } else { rcma=ConvertColsToInsertsCMSA(xcma,1,j,i); i=j; }
	   if(xcma != cma) NilCMSA(xcma); xcma=rcma;
	} free(ToRm);
	if(xcma != cma) return xcma; else return 0;
}

