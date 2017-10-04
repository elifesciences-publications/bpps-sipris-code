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

#include "rst_typ.h"

rst_typ::rst_typ(char Mode,a_type A) 
{ AB=A;  own_AB=FALSE; init_res_sets(Mode); }

rst_typ::rst_typ(char Mode) 
{ AB=MkAlpha(AMINO_ACIDS,GBLAST_BLOSUM62);  own_AB=TRUE; init_res_sets(Mode); }

void    rst_typ::FindEquivalent( )
// find sets that are equal to an earlier set
//  does Set[r][s] == Set[r'][s'] where r > r'?
{
	Int4	r,q,rset,qset;
	sst_typ	rsst;

	NEWP(EqualsPrevious,nAlpha(AB) + 3,BooLean);
	NEW(EqualsPrevious[0],max_residue_set +3, BooLean);
        for(r=1; r <= nAlpha(AB); r++){
	   NEW(EqualsPrevious[r],max_residue_set +3, BooLean);
	   for(rset=1; ResidueSSet[r][rset] != 0; rset++){
	     rsst = ResidueSSet[r][rset];
             for(q=1; q < r; q++){
		for(qset=1; ResidueSSet[q][qset] != 0; qset++){
		    if(rsst == ResidueSSet[q][qset]) EqualsPrevious[r][rset]=TRUE;
		}
	     }
	   }
	}
}

void	rst_typ::ValidityCheck(FILE *fp)
// check for programing errors...
{
	Int4	r,x,s,t,C,c,n,i,j;
	sst_typ **RST=this->LegalResSets( );
	sst_typ	*xsst,*rsst,sst,usst;
	for(r =1; r <= nAlpha(AB); r++) {
	   rsst=RST[r];
	   for(usst=0,s=0,i=1; rsst[i] != 0; i++){ usst=UnionSset(usst,rsst[i]); s++; }
#if 0 // DEBUG
	   fprintf(fp," r=%c: ",AlphaChar(r,AB));
	   PutSST(fp,rsst[0],AB); fprintf(fp," ==? "); PutSST(fp,usst,AB); fprintf(fp,"\n");
#endif
	   if(usst != rsst[0]){	// check the universal set for residue r.
		fprintf(fp," r=%c: ",AlphaChar(r,AB));
	      	PutSST(fp,usst,AB); fprintf(fp," != "); PutSST(fp,rsst[0],AB); fprintf(fp,"\n");
		print_error("FATAL: rst_typ - invalid universal residue set");
	   }
	   if(s != NumResSets[r]){
	      	PutSST(fp,rsst[1],AB);
		print_error("FATAL: rst_typ - residues sets != NumResSets");
	   }
	   for(C=1,s=1; rsst[s] != 0; s++){
	      sst=rsst[s]; n=CardSST(sst,AB);
	      if(n < C){		// Check to make sure sets get larger only.
	      	PutSST(fp,sst,AB);
		fprintf(fp," current size = %d; last size = %d\n",n,C);
		print_error("FATAL: rst_typ - residues set sizes out of order");
	      } C = n;
	      for(t=s+1; rsst[t] != 0; t++){	// check to ensure no redundant sets.
		if(sst == rsst[t]){
	      		PutSST(fp,sst,AB);
		  	fprintf(fp," (r=%c (%d  & %d)\n",AlphaChar(r,AB),s,t);
			print_error("FATAL: rst_typ - redundant residue sets");
		} 
	      }
	      for(x =1; x <= nAlpha(AB); x++){	// Check for symmetry of sets.
		if(x==r) continue;
		if(!MemSset(x,sst)) continue;
	        xsst=RST[x];
	        for(t=1; xsst[t] != 0; t++){ if(xsst[t] == sst) break; }   // should find a match.
		if(xsst[t] == 0){
	      	  PutSST(fp,sst,AB); fprintf(fp,"\n");
	      	  PutSST(fp,xsst[1],AB);
		  fprintf(fp," (r = %c; s = %d; x = %c; t = %d\n",AlphaChar(r,AB),s,AlphaChar(x,AB),t);
	          for(i=0; xsst[i] != 0; i++){ PutSST(fp,xsst[i],AB); fprintf(fp," "); } fprintf(fp,"\n");
		  print_error("FATAL: rst_typ - residue sets asymmetric");
		}
	      }
	   }
	}
}

void    rst_typ::Put(FILE *fp, a_type ab)
{
	if(ab==0) ab=AB;
	PutAlpha(stderr,ab);
        for(Int4 r = 1; r <= nAlpha(ab); r++){
		fprintf(fp,"%c: '%s'",AlphaChar(r,ab),ResSetString[r]);
		for(Int4 rset=1; ResidueSSet[r][rset] != 0; rset++){
		   sst_typ	sst = ResidueSSet[r][rset];
		   if(EqualsPrevious[r][rset]) fprintf(fp,",\""); else fprintf(fp,",'"); 
        	   for(unsigned char c = 1; c <= nAlpha(ab); c++){
			if(MemSset(c,sst)){
				fprintf(fp,"%c",AlphaChar(c,ab));
			}
		   } if(EqualsPrevious[r][rset]) fprintf(fp,"\""); else fprintf(fp,"'"); 
		} fprintf(fp,"\n");
	}
}

void    rst_typ::Init_Res_Sets(const char *ResidueSets[21][res_array_size])
{
        char    rc;
        Int4    rsets,rs,r,r1,r2;
	sst_typ	sst;

        assert(nAlpha(AB) == 20);
	NEWP(ResidueSSet,nAlpha(AB)+3,sst_typ);
	NEW(NumResSets,nAlpha(AB)+3,unsigned char);
	NEWP(ResSetString,nAlpha(AB)+3,char);
        for(r = 0; r <= nAlpha(AB); r++){
	    NEW(ResidueSSet[r],max_residue_set+3,sst_typ); // for r=0; leave as empty sst...
	    ResSetString[r] = 0;
        }
        for(r = 1; r <= nAlpha(AB); r++){
		ResSetString[r] = AllocString(ResidueSets[r][0]);
                r1=AlphaCode(ResidueSets[r][1][0],AB);  // first set.
		assert(NumResSets[r1] == 0); 		// first time this res.
                for(rsets=0; rsets <= max_residue_set && ResidueSets[r][rsets]!=0; rsets++){
                   for(rs=0; (rc=ResidueSets[r][rsets][rs]) != 0; rs++){
                        r2=AlphaCode(rc,AB);
			if(rs==0){
			    if(r1 != r2){
				fprintf(stderr,"r1 = %c; r2 = %c\n",AlphaChar(r1,AB),AlphaChar(r2,AB));
			    } assert(r2==r1);	// format check.
			}
			sst = SsetLet(r2);
                        ResidueSSet[r1][rsets]=UnionSset(ResidueSSet[r1][rsets],sst);
                   }
                } rsets--; 
		NumResSets[r1] = rsets;
		// assert(rsets == NumResSets[r1]);	// numbering check
        } FindEquivalent( );
}

void    rst_typ::init_res_sets(char mode)
{
	if(islower(mode)) mode=toupper(mode); 
	// fprintf(stderr,"sets_mode = %c\n",mode);
	const char *ResidueSets_O[21][res_array_size] = {	// one_only
                {"X","X",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                {"C","C",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                {"G","G",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, // G + A(0.0798).
                {"A","A",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                {"S","S",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                {"T","T",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                {"N","N",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                {"D","D",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                {"E","E",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                {"Q","Q",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                {"K","K",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                {"R","R",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                {"H","H",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                {"W","W",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                {"Y","Y",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                {"F","F",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                {"V","V",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                {"I","I",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                {"L","L",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                {"M","M",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                {"P","P",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}};
	//************************ WARNING: READ THIS! **********************************
	// WARNING: need to include all residues in first string above 
	// and need to start each string with key residue.
	//*******************************************************************************
	const char *ResidueSets_M[21][res_array_size] = {	// Medium (Based on chemical-geometrical mymicry)
		// e.g., Q cat look like part of H if in the right conformation.
                {"X","X",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                {"C","C",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                {"GA","G","GA",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, 
                {"ASG","A","AS","AG",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                {"SATN","S","SA","ST","SN",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                {"TS","T","TS",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                {"NSDHQE","N","NS","ND","NH","NQ","NQH","NDE",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                {"DNEQ","D","DN","DE","DNE","DEQ",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                {"EDQN","E","ED","EQ","EDQ","EDN",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                {"QENKRHD","Q","QE","QK","QR","QH","QN","QED","QKR","QHN",0,0,0,0,0,0,0,0 }, 
                {"KRQ","K","KQ","KR","KQR",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                {"RQK","R","RQ","RK","RQK",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                {"HYNQ","H","HY","HN","HQ","HNQ",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                {"WYF","W","WY","WYF",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                {"YHWF","Y","YH","YW","YF","YHW","YHF","YWF","YHWF",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                {"FYWL","F","FY","FL","FWY",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                {"VILM","V","VI","VL","VIL","VILM",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                {"IVLM","I","IV","IL","IVL","IVLM",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                {"LVIMF","L","LV","LI","LM","LF","LVI","LVIM",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                {"MVIL","M","ML","MVIL",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                {"P","P",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}};
	const char *ResidueSets_L[21][res_array_size] = {	// Large (Based on donor type; e.g., -OH) 
                {"X","X",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                {"C","C",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                {"GA","G","GA",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                {"ASG","A","AS","AG",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                {"SATN","S","SA","ST","SN","STN",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                {"TSN","T","TS","TN","TSN",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                {"NSTDHEQ","N","NS","NT","ND","NQ","NH","NST","NDE","NQE","NHQ","NQED",0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                {"DNEQ","D","DN","DE","DEN","DEQ","DEQN",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                {"EDQN","E","ED","EQ","EDQ","EDN","EQN","EDNQ",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                {"QENDHKR","Q","QE","QN","QK","QR","QH","QKR","QED","QHN","QNE","QNDE",0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                {"KRQ","K","KQ","KR","KQR",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                {"RQK","R","RQ","RK","RQK",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                {"HYNQ","H","HY","HN","HQ","HNQ",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                {"WYF","W","WF","WY","WYF",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                {"YHWF","Y","YH","YW","YF","YFW",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                {"FWY","F","FW","FY","FWY",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                {"VILM","V","VI","VL","VM","VLM","VIM","VIL","VILM",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                {"IVLM","I","IV","IL","IM","IVL","IVM","ILM","IVLM",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                {"LVIM","L","LV","LI","LM","LVI","LVM","LIM","LVIM",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                {"MVIL","M","MV","MI","ML","MVI","MVL","MIL","MVIL",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                {"P","P",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}};
	const char *ResidueSets_G[21][res_array_size] = {	// "Gold standard" (only one based on careful pattern analysis)
                {"X","X",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                {"CASV","C","CA","CS","CV","CVA","CAS",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                {"GASND","G","GA","GS","GD","GN","GAS","GND",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                {"ASTGVCILPE","A","AS","AG","AC","AT","AV","AP","AE","AGS","AST","ATV","ACS","ASP","AVI","ACV","AVL",
			"AVIL",0,0,0,0,0,0,0,0},
                {"SATNGPC","S","SA","ST","SN","SC","SG","SP","STN","SAG","SAT","SCA","SAP",0,0,0,0,0,0,0,0,0,0,0,0,0},
                {"TSNAVI","T","TS","TN","TA","TV","TI","TSN","TSA","TAV","TVI",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                {"NSTDHEQGK","N","NS","NT","ND","NQ","NH","NG","NK","NST","NDE","NHQ","NGD",0,0,0,0,0,0,0,0,0,0,0,0,0},
                {"DNEQG","D","DN","DE","DG","DEN","DEQ","DNG",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                {"EDQNAK","E","ED","EQ","EA","EK","EDQ","EDN",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                {"QENDHKRPM","Q","QE","QN","QK","QR","QH","QP","QM","QKR","QED","QHN",0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                {"KRQNE","K","KQ","KR","KN","KE","KQR",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                {"RQKH","R","RQ","RK","RH","RQK",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                {"HYNQRF","H","HY","HN","HQ","HR","HNQ","HYF",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                {"WYF","W","WF","WY","WYF",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                {"YHWFL","Y","YH","YW","YF","YL","YFW","YFL","YHF",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                {"FWYILVHM","F","FW","FY","FI","FL","FM","FLI","FWY","FYL","FLM","FYH","FVI",
			"FLIV","FVILM",0,0,0,0,0,0,0,0,0,0,0},
                {"VILMFATC","V","VI","VL","VM","VA","VT","VC","VLM","VIM","VIL","VAT","VAI","VCA","VTI","VIF","VAL",
			"VILM","VAIL","VILF","VFILM",0,0,0,0,0},
                {"IVLMFAT","I","IV","IL","IM","IF","IT","ILF","IVL","IVM","ILM","IVA","ITV","IVF",
			"IVLM","IVLF","IVAL","IFVLM",0,0,0,0,0,0,0,0},
                {"LVIMFAY","L","LV","LI","LM","LF","LY","LVI","LVM","LIM","LIF","LYF","LMF","LVA",
			"LVIM","LVIF","LVIA","LFVIM",0,0,0,0,0,0,0,0},
                {"MVILFQ","M","MV","MI","ML","MF","MQ","MVI","MVL","MIL","MLF","MVIL","MFVIL",0,0,0,0,0,0,0,0,0,0,0,0,0},
                {"PASQ","P","PA","PS","PQ","PAS",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}};
#if 1	// Roughly based on blosum45 matrix...
	const char *ResidueSets_R[21][res_array_size] = {	// Root (for distant structural homologs) 
                {"X","X",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                {"CASVLIM","C","CA","CS","CV","CVA","CAS","CAVILM",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                {"GASND","G","GA","GS","GD","GN","GAS","GND","GSN",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                {"ASTGVCILMPE","A","AS","AG","AC","AT","AV","AP","AE","AGS","AST","ATV","ACS","ASP","AVI","ACV","AVL",
			"AVIL","ACVILM",0,0,0,0,0,0,0},
                {"SATNGPCQ","S","SA","ST","SN","SC","SG","SP","STN","SAG","SAT","SCA","SAP","SNQ","SGN",0,0,0,0,0,0,0,0,0,0,0},
                {"TSNAVI","T","TS","TN","TA","TV","TI","TSN","TSA","TAV","TVI",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                {"NSTDHEQGK","N","NS","NT","ND","NQ","NH","NG","NK","NST","NSQ","NDE","NQE","NHQ","NGS","NQK","NGD",
			"NQED","NDEQK",0,0,0,0,0,0,0},
                {"DNEQGRK","D","DN","DE","DG","DEN","DEQ","DNG","DEQN","DNEQK","DRKQE",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                {"EDQNARK","E","ED","EQ","EA","EK","EDQ","EDN","ENQ","EDNQ","EKQR","EDNQK","ERKQD",0,0,0,0,0,0,0,0,0,0,0,0,0},
                {"QENDHKRSPM","Q","QE","QN","QK","QR","QH","QP","QM","QKR","QED","QHN","QSN","QNK","QNE",
			"QNDE","QEKR","QDNEK","QRKED",0,0,0,0,0,0,0},
                {"KERQND","K","KQ","KR","KN","KE","KQR","KNQ","KEQR","KDNEQ","KRQED",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                {"REQKDH","R","RQ","RK","RH","RQK","REKQ","RKQED",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                {"HYNQWRF","H","HY","HN","HQ","HR","HNQ","HYF","HWYF",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                {"WYFH","W","WF","WY","WYF","WHYF",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                {"YHWFMLIV","Y","YH","YW","YF","YL","YFW","YFL","YHF","YFM","YWHF","YFMLIV",0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                {"FWYILVHM","F","FW","FY","FI","FL","FM","FLI","FWY","FYL","FLM","FYH","FVI","FYM",
			"FLIV","FILM","FWHY","FLMIV","FYMLIV",0,0,0,0,0,0,0},
                {"VILMFYATC","V","VI","VL","VM","VA","VT","VC","VLM","VIM","VIL","VAT","VAI","VCA","VTI","VIF","VAL",
			"VILM","VAIL","VILF","VLFMI","VYFMLI","VACILM",0,0,0},
                {"IVLMFYATC","I","IV","IL","IM","IF","IT","ILF","IVL","IVM","ILM","IVA","ITV","IVF",
			"IVLM","IVLF","IVAL","IMLF","ILFMV","IFYMLV","IACVLM",0,0,0,0,0},
                {"LVIMFAYC","L","LV","LI","LM","LF","LY","LVI","LVM","LIM","LIF","LYF","LMF","LVA",
			"LVIM","LVIF","LVIA","LMIF","LFMIV","LFYMIV","LACVIM",0,0,0,0,0},
                {"MYFVILACQ","M","MV","MI","ML","MF","MQ","MVI","MVL","MIL","MLF","MYF",
			"MVIL","MILF","MLVFI","MFYLIV","MACVIL",0,0,0,0,0,0,0,0,0},
                {"PASQ","P","PA","PS","PQ","PAS",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}};
// "ACVILM" "CAVILM" "VACILM" "IACVLM" "LACVIM" "MACVIL"
#endif

        // sst_typ ResidueSSet[21][9];

#if 0	// Test	
        for(r = 1; r <= nAlpha(AB); r++){
	  fprintf(stdout,"%c: ",AlphaChar(r,AB));
          for(r2 = 1; r2 <= nAlpha(AB); r2++){
		if(blosum62[r][r2] > 0.0){  // == real value log-odds.
	  	   fprintf(stdout,"%c(%.4f)",AlphaChar(r2,AB),blosum62[r][r2]);
		}
	  }
	  fprintf(stdout,"\n");
	}
	exit(1);
#endif

   switch(mode){
     case 'O': Init_Res_Sets(ResidueSets_O); break;
     case 'M': Init_Res_Sets(ResidueSets_M); break;
#if 0	// test using 'R' option throughout...
     case 'L': Init_Res_Sets(ResidueSets_R); break;
#else
     case 'L': Init_Res_Sets(ResidueSets_L); break;
#endif
     case 'G': Init_Res_Sets(ResidueSets_G); break;
#if 0	// Test effect of NOT using root residue environment patterns...
     case 'R': Init_Res_Sets(ResidueSets_L); break;
#else
     case 'R': Init_Res_Sets(ResidueSets_R); break;
#endif
     case 0: Init_Res_Sets(ResidueSets_L); break;
     default: 
	fprintf(stderr,"fatal error in rst_typ: invalid mode = '%c'(%d)\n",mode,(int)mode);
	print_error("init_res_sets( ) input error"); 
      break;
   } ValidityCheck(stderr); 
}

