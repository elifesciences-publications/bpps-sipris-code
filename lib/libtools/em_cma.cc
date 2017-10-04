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

#include "em_cma.h"

cma_typ	PurgeCMSA(double cutoff, cma_typ cma)
{
        FILE *fp=tmpfile(); 
	PutGoodCMSA(fp,cutoff,cma);
	rewind(fp);
        cma_typ cma2=ReadCMSA(fp,AlphabetCMSA(cma)); fclose(fp);
	return cma2;
}

cma_typ	PurgeCMSA(Int4 min_rpt, Int4 min_spacing, cma_typ cma)
// remove 
{
        FILE *fp=tmpfile(); 
        PutGoodRptsCMSA(fp,min_rpt,min_spacing,cma);
	rewind(fp);
        cma_typ cma2=ReadCMSA(fp,AlphabetCMSA(cma));
        fclose(fp);
	return cma2;
}

cma_typ	CleanCMSA(cma_typ cma, double N, double C, double c, char *name)
{
	char	str[200];

	sprintf(str,"%s_tmp",name); WriteMtfCMSA(str,cma,NULL);

        char *Argv[6];
        char str0[]="cleancma";
        char str2[50],str3[50],str4[50];
        Argv[0] = str0; Argv[1] = str; Argv[2] = str2; 
        Argv[3] = str3; Argv[4] = str4; Argv[5]=0;

	sprintf(str2,"-N%g",N); sprintf(str3,"-C%g",C); sprintf(str4,"-c%g",c);
        fprintf(stderr,"cleancma %s %s %s %s\n",Argv[1],Argv[2],Argv[3],Argv[4]);

        CleanCMSA(5,Argv,AlphabetCMSA(cma));
        sprintf(str,"%s_tmp.cln.cma",name);
        cma_typ cma2=ReadCMSA2(str,AlphabetCMSA(cma));
	return cma2;
}

cma_typ	EMrunCMSA(Int4 max_rpts,char *psm_arg, char *name, cma_typ cma)
// add idp_typ idp  to RelMapCMSA(cma) arguments.
{

	BooLean flag;
	Int4	score,start,oper_len,i;
	char	str[108];

	// oldmap=RelMapCMSA(cma);
        for(i=0,flag=TRUE; flag; i++){
                std::cerr << "        ************************* EM ********************\n";
                flag=FALSE;
                psm_typ pssm(max_rpts,psm_arg,cma); // recompute this each cycle
                for(Int4 sq=1;sq<=NumSeqsCMSA(cma);sq++){
                 gss_typ *gss=gssCMSA(cma);
                 e_type E=gss->TrueSeq(sq);
                 char *operation = pssm.Align(E,1, &score, &start, &oper_len);
                 // if(ReAlignGSqCMSA(sq,operation,start,&cma)){ // OLD....
		 if(gss->NewAlign(operation,oper_len,start,sq)){
		    if(!ReAlignGSqCMSA(sq,operation,start,&cma)){ 
                        fprintf(stderr,"score=%d sq=%d\n",score,sq);
                        std::cerr << operation; std::cerr << "\n";
			gss->Put(stderr,sq);
			print_error("EMrunCMSA( ): False Positive");
		    }
                    std::cerr << operation; std::cerr << "\n";
                    fprintf(stderr,"score=%d sq=%d\n",score,sq);
                    pssm.ReComputeSMX(cma); flag=TRUE;
#if 0
			map = RelMapCMSA(cma);
			fprintf(stderr,"map %.2f --> %.2f\n", oldmap,map);
			oldmap = map;
#endif
#if 0		// Check gss->NewAlign( ) to make sure it works...
		 } else if(ReAlignGSqCMSA(sq,operation,start,&cma)){
			fprintf(stderr,"score=%d sq=%d\n",score,sq);
		 std::cerr << operation; std::cerr << "\n";
			print_error("EMrunCMSA( ): False Negative");
#endif
                 } free(operation);
                }
                sprintf(str,"%s.em",name); WriteMtfCMSA(str,cma,NULL);
	} return cma;
}

cma_typ	HMM_EMrunCMSA(Int4 max_rpts,char *psm_arg, char *name, cma_typ cma)
// add idp_typ idp  to RelMapCMSA(cma) arguments.
{

	BooLean flag;
	Int4	score,start,oper_len,i,pernats=200;
	char	str[108];
	double	map,oldmap;

	oldmap=0.0; //RelMapCMSA(cma);
	Int4	same=0;
        for(i=0,flag=TRUE; flag && same < 4; i++){
                std::cerr << "        ************************* EM ********************\n";
                flag=FALSE;
                HMM_typ hmm(max_rpts,psm_arg,cma,pernats,0); // recompute this each cycle
#if 0
map=hmm.RelMap(cma);
fprintf(stderr,"map = %.2f\n",map);
exit(1);
#endif
                for(Int4 sq=1;sq<=NumSeqsCMSA(cma);sq++){
                 gss_typ *gss=gssCMSA(cma);
                 e_type E=gss->TrueSeq(sq);
                 char *operation = hmm.Align(E,1, &score, &start, &oper_len);
// std::cerr << operation; std::cerr << "\n";
//fprintf(stderr,"start = %d\n",start);
// hmm.PutAlign(stderr,E,operation,oper_len,start,1);
//fprintf(stderr,"Score = %d\n ======================== \n",score);
		 if(gss->NewAlign(operation,oper_len,start,sq)){
		    if(!ReAlignGSqCMSA(sq,operation,start,&cma)){ 
                        fprintf(stderr,"score=%d sq=%d\n",score,sq);
                        // std::cerr << operation; std::cerr << "\n";
			gss->Put(stderr,sq);
			print_error("EMrunCMSA( ): False Positive");
		    }
                    std::cerr << operation; std::cerr << "\n";
                    fprintf(stderr,"score=%d sq=%d\n",score,sq);
                    hmm.ReComputeSMX(cma); flag=TRUE;
#if 0
			map = RelMapCMSA(cma);
			fprintf(stderr,"map %.2f --> %.2f\n", oldmap,map);
			oldmap = map;
#endif
                 } free(operation);
                }
		map = RelMapCMSA(cma);
		// map = hmm.RelMap(cma);
		fprintf(stderr,"map %.2f --> %.2f\n", oldmap,map);
		if(map == oldmap) same++; else { oldmap=map; same=0; }
                sprintf(str,"%s.em",name); WriteMtfCMSA(str,cma,NULL);
	} return cma;
}

