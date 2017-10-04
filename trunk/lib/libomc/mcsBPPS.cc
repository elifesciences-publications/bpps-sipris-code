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

#include "mcs_typ.h"
#include "blosum62.h"

Int4    RunmcsBPPS(Int4 argc,char *argv[])
{
        Int4	a,iter=0,NumRandom=0;
	BooLean	Converged=FALSE;
        time_t  time1=time(NULL);
	mcs_typ	*cmc=0;
	FILE *fp=0; 
	hsw_typ hsw=0;
	swt_typ *swt=0;
	cma_typ TrueMainCMA=0;

	if(argc < 2){ cmc = new mcs_typ( ); return 0; }
	char	str[1000];
        fp = open_file(argv[1],".cmd","w");
	for(a = 0; a < argc; a++) fprintf(fp,"%s ",argv[a]);
	fprintf(fp,"\n"); fclose(fp);
	// fp=open_file(argv[1],".mma","r");
	sprintf(str,"%s.mma",argv[1]);
        if((fp=fopen(str,"r")) == NULL){        // call other constructor for hpt operations
           cmc = new mcs_typ( ); return 0; // don't actually run program; show usage only.
	}
	a_type  AB=MkAlpha(AMINO_ACIDS,PROT_BLOSUM62);
        TrueMainCMA=ReadCMSA(fp,AB);
        fclose(fp); fp=0;
        if(NumRandom==0) NumRandom=1+(NumSeqsCMSA(TrueMainCMA)/3); // sam

        sprintf(str,"%s.hsw",argv[1]);
        if((fp=fopen(str,"r")) == NULL){        // create file...
                swt = new swt_typ(TrueMainCMA,FALSE);
                hsw=swt->RtnHSW( );
                fp = open_file(argv[1],".hsw","w");
                FWriteHSW(fp,hsw); fclose(fp);
        } else { hsw=FReadHSW(fp,AB,TrueMainCMA); fclose(fp); swt = new swt_typ(hsw); }

	cma_typ rcma,in_mcma=MkMainFileCMSA(TrueMainCMA,NumRandom,rcma);
        hsw_typ HSW=AddRandomHSW(hsw,TrueMainCMA,rcma,in_mcma);
        TotalNilCMSA(rcma);     // destroy the temporary Random sequence alignment.

	cmc = new mcs_typ(TrueMainCMA,in_mcma,HSW,argc,argv);
	fflush(stdout);
        do {
		iter++;
		// break;
		Converged=cmc->Sample( );
		double lpr=cmc->CalcTotalLPR();
		if(lpr <= 0.0)
			print_error("Failed to find significant set assignments for this Hyperpartition");
	} while(!Converged);
	cmc->RestoreBest();
        cmc->Put();     // creates <infile>_grp.chn
	fprintf(stderr,"done printing results\n");
	cmc->PutHyperPartition( );
	if(cmc->IsTreeHpt){
          fp= open_file(argv[1],".contrib","w");
	  cmc->PutMapContributions(fp); fclose(fp);
	}
#if 1   // save sets...
	if(cmc->SaveSets){
           set_typ *sets=cmc->CopyOfSeqSets();
           Int4 NumSets=cmc->RtnNumElmntSetCMA( );
           fp=open_file(argv[1],".sets","w");
           WriteSets(fp,NumSets,sets); fclose(fp);
	}
#endif

        delete cmc;
	if(swt) delete swt; 
	if(TrueMainCMA) TotalNilCMSA(TrueMainCMA); // needed with new random sequence method...
        double runtime=difftime(time(NULL),time1);
        fprintf(stderr, "\ttime(cmcBPPS): %0.1f seconds (%0.2f minutes)\n",runtime,runtime/60.0);
	return 0;
}

