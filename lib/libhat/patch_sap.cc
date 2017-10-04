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

#include "patch_sap.h"

Int4	*LongestPathGSAP(Int4 num_saps, sap_typ *SAP,UInt4 *QS,
		UInt4 *QE,UInt4 *SS,UInt4 *SE)
// Find the longest path through the sap_typ's 
/*********************************************************************
Case 1: Overlap in query.

        qs1           qe1
        : ......... qs2:            qe2
   -----########====####=====##########---------------- Query
       /       |   /\ /\     |    |\   \
      /        |  /  X  \    |    | \   \
     /         | /  / \  \   |    |  \   \
    /          |/  /   \  \  |    |   \   \
   ########==######-----###==######====#####-------	Subject
   ss1 ........ se1     ss2 ............ se2                      

Case 2: Overlap in Subject.
 or

   s1            e1     s2................e2                      
   ########==######-----###==######====#####-------	Query
    \          |\  \   /  /  |    |   /   /
     \         | \  \ /  /   |    |  /   /
      \        |  \  X  /    |    | /   /
       \       |   \/ \/     |    |/   /
   -----########====####=====##########--------------   Subject
        :           s2 :             e2
        s1            e1

  o-->o    (colinear in query and subject sequences...)
	   weight = length of target node (i.e., gapped HSP).

 *********************************************************************/
{
	wdg_typ 	G;
	UInt4	qs1,qe1,ss1,se1;
	UInt4	qs2,qe2,ss2,se2;
	Int4		n,n1,n2,nodes,edges,s1,s2;
	Int4		begin,end,*path,*tmp_path,*dist;
	Int4		score1,score2,e,v,r,s;

	nodes = num_saps + 2; edges=nodes*nodes/2;
	G = MkWdgraph(nodes,edges);
	begin=num_saps + 1; end = num_saps + 2;
	for(s1=1; s1 < num_saps; s1++){
	  qs1=QS[s1]; ss1=SS[s1]; qe1=QE[s1]; se1=SE[s1];
	  score1=-ScoreGSAP(SAP[s1]);
#if 0	// DEBUG
PutSeqID(stderr,QuerySeqGSAP(SAP[s1])); fprintf(stderr,"\n");
PutSeqID(stderr,SubjectSeqGSAP(SAP[s1])); fprintf(stderr,"\n");
#endif
	  for(s2=s1+1; s2 <= num_saps; s2++){
	     qs2=QS[s2]; ss2=SS[s2]; qe2=QE[s2]; se2=SE[s2];
	     score2=-ScoreGSAP(SAP[s2]);
	     if(qs1 < qs2){		
		if(ss1 < ss2) JoinWdgraph(s1,s2,score2,G);	// HSPs colinear
		else if(ss1 >= ss2) JoinWdgraph(s2,s1,score1,G);	// repeat...
		else print_error("This Shouldn't happen\n");
	     } else if(qs1 > qs2){		
		if(ss1 > ss2) JoinWdgraph(s2,s1,score1,G);	// HSPs colinear
		else if(ss1 <= ss2) JoinWdgraph(s1,s2,score2,G);	// repeat...
		else print_error("This Shouldn't happen\n");
	     } else if(qs1 == qs2){				// possible repeat...
		if(ss1 > ss2) JoinWdgraph(s2,s1,score1,G);	// repeat...
		else if(ss1 <= ss2) JoinWdgraph(s1,s2,score2,G);	// repeat...
		else print_error("This Shouldn't happen\n");
	     } else print_error("This Shouldn't happen\n");
	  }
	}
	Int4	first=0;
	for(s1=1; s1 <= num_saps; s1++){
	  qs1=QS[s1]; ss1=SS[s1]; qe1=QE[s1]; se1=SE[s1];
	  score1=-ScoreGSAP(SAP[s1]);
	  if(!FirstInWdgraph(s1,G)){ assert(first == 0); first=s1; }
	  JoinWdgraph(begin,s1,score1,G);
	  JoinWdgraph(s1,end,-10,G);
	}
	// PutWdgraph(stderr, G);
	NEW(path,nodes+4,Int4); NEW(dist,nodes+4,Int4); NEW(tmp_path,nodes+4,Int4);
	TopoScanWdigraph(G, begin, path, dist);
        // fprintf(stderr,"\nstart ");
	for(s1=1,n=path[end]; n!=begin; n=path[n]){
                // fprintf(stderr,"-> %d",n); 
		tmp_path[s1]=n; s1++;
        } // fprintf(stderr," (dist = %d)\n\n",-dist[end]);
	for(s2=1; s2 < s1; s2++) path[s2]=tmp_path[s1-s2];
	assert(s2-1 == nodes-2);
	NilWdgraph(G); free(dist); free(tmp_path);
	return path;
}

Int4    RawInfoGSAP(Int4 *QS, Int4 *QE, Int4 *SS, Int4 *SE, sap_typ sap)
{
        Int4	s,q,index;

        assert(sap->segs);
        dsp_typ dsp = sap->segs;
        *QS = dsp->starts[0]; *SS = dsp->starts[1];
        for(s=q=index=0; index < dsp->numseg; index++){
                if(dsp->starts[2*index] != -1) q+=dsp->lens[index];
                if(dsp->starts[2*index+1] != -1) s+=dsp->lens[index];
        } *QE = *QS+q-1; *SE = *SS+s-1;
        return s;
}


Int4	*RunningScoreSAP_L(sap_typ sap,Int4 **mx,a_type A)
// compute the running score to the left and return running score array...
//                 start <--- end
{
        Int4    S,I,J,r,i,len,score=0,*rtn_score;
        Int4    ins_open=11,ins_extend=1; // fix this later!!
	e_type	qE,sE;

        assert(sap->segs);
        dsp_typ	dsp=sap->segs;
        qE = dsp->qE; sE = dsp->sE;
	NEW(rtn_score,LenSeq(qE)+3,Int4);
        for(J=score=0,I=dsp->numseg-1; I >= 0; I--){
           len= dsp->lens[I];
	   J=dsp->starts[2*I];
           S=dsp->starts[2*I+1];
           if(J != -1 && S != -1){   // matching region
                for(i=0,J=J+len-1,S=S+len; i < len; S--,i++,J--){
                   r=ResSeq(S,sE);
                   if(mx) score += mx[J][r];
                   else score += valAlphaR(r,ResSeq(J+1,qE),A);
	 	   rtn_score[J]=score;
		   // assert(J < LenSeq(qE));
                }
           } else if(J == -1){        // deletion in query sequence.
		score -= (ins_extend*len + ins_open); 
           } else if(S == -1){        // insertion in query.
		score -= ins_open;
		for(J=J+len-1,i=0; i < len; J--,i++){
		   // assert(J < LenSeq(qE));
		   score -= ins_extend; rtn_score[J]=99999999;
		}		// prohibit trimming within gaps...
           }
	} return rtn_score;
}

Int4	*RunningScoreSAP_R(sap_typ sap,Int4 **mx,a_type A)
// compute the running score to the right and return in array
{
        Int4    I,J,r,i,S,len,score=0,*rtn_score;
	e_type	qE,sE;
        Int4    ins_open=11,ins_extend=1; // fix this later!!

        assert(sap->segs);
        dsp_typ	dsp=sap->segs;
        qE = dsp->qE; sE = dsp->sE;
	NEW(rtn_score,LenSeq(qE)+3,Int4);
        for(score=I=0; I < dsp->numseg; I++){
	   J=dsp->starts[2*I];
           S=dsp->starts[2*I+1];
           len= dsp->lens[I];
           if(J != -1 && S != -1){   // matching region
              for(i=0; i < len; S++,J++,i++){
                r=ResSeq(S+1,sE);
                if(mx) score += mx[J][r];
                else score += valAlphaR(r,ResSeq(J+1,qE),A);
		rtn_score[J]=score;
		// assert(J < LenSeq(qE));
	      }
           } else if(J == -1){        // insertion in subject sequence.
		score -= (ins_extend*len + ins_open); 
           } else if(S == -1){        // insertion in query.
		for(score -= ins_open,i=0; i < len; J++,i++){
		   // assert(J < LenSeq(qE));
		   score -= ins_extend; rtn_score[J]=99999999;
		}		// prohibit trimming within gaps...
           }
	} return rtn_score;
}


Int4	QueryToSubjectSiteSAP(sap_typ sap,Int4 Qsite)
// return the site in the subject sequence that corresponds to the 
// input site in the query. 
// Return -1 if site is outside of the alignment.
{
        Int4    I,J,r,i,S,SS,JJ,len,score=0,*rtn_score;
	e_type	qE,sE;

        assert(sap->segs);
        dsp_typ	dsp=sap->segs;
        qE = dsp->qE; sE = dsp->sE;
	JJ=dsp->starts[0]; SS=dsp->starts[1];
        for(I=0; I < dsp->numseg; I++){
	   J=dsp->starts[2*I]; S=dsp->starts[2*I+1]; len= dsp->lens[I];
           if(J != -1 && S != -1){   // matching region
              for(i=0; i < len; SS++,i++,JJ++) if(JJ==Qsite) return SS;
           } else if(J == -1){        // insertion in subject.
		for(i=0; i < len; i++) SS++;
           } else if(S == -1){        // insertion in query.
		for(i=0; i < len; i++,JJ++) if(JJ==Qsite) return SS; 
           }
	} return -1;
}

BooLean	TrimOverlappingSAPs(sap_typ tail_sap, sap_typ head_sap,Int4 crosspoint)
/***************************************************************************
                            crosspoint
              qs1           |   qe1				tail_sap
    ----------###########==#######----------------------------- query
    ----------#######==###########----------------------------- subject
              ss1___________    se1				
                        qs2 L____________qe2
    --------------------###################-------------------- query
    --------------------###################-------------------- subject
                        ss2 |            se2			head_sap
                            s_xpt
 ***************************************************************************/
{
	Int4	trimleft,trimright,s_xpt;
	Int4	qs1,qe1,ss1,se1,qs2,qe2,ss2,se2;
	sap_typ tmp_sap;
	FILE	*fp=0;

	assert(QuerySeqGSAP(tail_sap) == QuerySeqGSAP(head_sap));
	assert(SubjectSeqGSAP(tail_sap) == SubjectSeqGSAP(head_sap));
        Int4 len1=RawInfoGSAP(&qs1,&qe1,&ss1,&se1,tail_sap);
        Int4 len2=RawInfoGSAP(&qs2,&qe2,&ss2,&se2,head_sap);
	assert(crosspoint >= qs2 && crosspoint <= qe1+1);
	if(crosspoint <= qe1){		// then need to trim tail_sap
	  s_xpt=QueryToSubjectSiteSAP(tail_sap,crosspoint);
	  if(fp) fprintf(fp,"crosspoint %d in query == %d in subject\n", crosspoint+1,s_xpt+1);
	  if(s_xpt < 0) print_error("TrimOverlappingSAPs( ) error");
	  assert(s_xpt >= ss1 && s_xpt <= se1);
	  if(s_xpt <= se1){
	     trimright = se1-s_xpt+1;
	     assert(trimright < len1);
	     if(fp) fprintf(fp,"trimming %d residues (%d..%d) from right\n", trimright,s_xpt,se1);
	     tmp_sap = TrimSAP(tail_sap,0,trimright);
	     if(tmp_sap==0) return FALSE;
	     GDenseSegFree(tail_sap->segs);
	     tail_sap->segs=tmp_sap->segs; tmp_sap->segs=0;
	     GSeqAlignFree(tmp_sap);
	  }
	}
	if(crosspoint > qs2){		// then need to trim head_sap
	  s_xpt=QueryToSubjectSiteSAP(head_sap,crosspoint);
	  if(fp) fprintf(fp,"crosspoint %d in query == %d in subject\n",
				crosspoint+1,s_xpt+1);
	  if(s_xpt < 0) print_error("TrimOverlappingSAPs( ) error");
	  assert(s_xpt >= ss2 && s_xpt <= se2);
	  if(s_xpt > ss2){
	    trimleft = s_xpt-ss2;
	    assert(trimleft < len2);
	    if(fp) fprintf(fp,"trimming %d residues (%d..%d) from left\n",
			trimleft,ss2,s_xpt);
	    tmp_sap = TrimSAP(head_sap,trimleft,0);
	    if(tmp_sap==0) return FALSE;
	    GDenseSegFree(head_sap->segs);
	    head_sap->segs= tmp_sap->segs; tmp_sap->segs=0;
	    GSeqAlignFree(tmp_sap);
	  }
	} return TRUE;
}

Int4	*EvolvingRmOverlapsGSAP(Int4 overlap_cutoff, sap_typ *SAP,Int4 *path,
					Int4 num,Int4 **mx,a_type A)
// remove overlapping regions
/*********************************************************************
Colinear adjacent HSP overlap... (non-adjacent,non-colinear == likely repeats....)

   segments overlap if(s1 <= s2 && s2 <= e1 || s2 <= s1 && s1 <= e2)
 
Case 1: Overlap in query.

        qs1            qe1
        :          qs2 :            qe2
    ----########====####=====##########--------------   Query
       /       |   /\ /\     |    |\   \
      /        |  /  X  \    |    | \   \
     /         | /  / \  \   |    |  \   \
    /          |/  /   \  \  |    |   \   \
   ########==######-----###==######====#####-------	Subject
   ss1          se1     ss2 ............. se2                      

Case 2: Overlap in Subject.
 or

   qs1          qe1     qs2               qe2                      
   ########==######-----###==######====#####-------	Query
    \          |\  \   /  /  |    |   /   /
     \         | \  \ /  /   |    |  /   /
      \        |  \  X  /    |    | /   /
       \       |   \/ \/     |    |/   /
   -----########====####=====##########---------------  Subject
        :          ss2 :            se2
        ss1           se1

 *********************************************************************/
{
	sap_typ	sap,sap1,sap2,sap_rtn=0;
	Int4	qs1,qe1,qs2,qe2,ss1,ss2,se1,se2;
	dsp_typ	dsp1,dsp2;
        Int4	r,i,j,len1,len2,*rpt;
#if 0
		  FILE *fp=stderr; // debug...
#else
		  FILE *fp=0; 
#endif

	// 3. Find overlapping ends.
	Int4	num_ol=0;
	NEW(rpt,num+2,Int4);	// return repeat assignments...
	for(i=1; i<=num; i++) rpt[i]=1;
	r=1;	// repeat # 1.
        for(Int4 V=1; V < num; V++){
	   i=path[V]; rpt[V]=r;
	   sap1=SAP[i]; dsp1=sap1->segs; 
           len1=RawInfoGSAP(&qs1,&qe1,&ss1,&se1,sap1);
           j=path[V+1];	// get adjacent sap...
	   sap2=SAP[j]; dsp2=sap2->segs; 
           len2=RawInfoGSAP(&qs2,&qe2,&ss2,&se2,sap2);

           if(ss1<=ss2 && ss2<=se1){ 		// Case 2: Overlap in Subject:
		if(fp) fprintf(fp,"Case 2a: \n");
		r++; rpt[V+1]=r; // treat this as a repeat for now...
	   } else if(ss2<=ss1 && ss1<=se2){
	     if(se1 <= se2){	// sap1 contained in sap2 -> probably compositional bias.
		r++; rpt[V+1]=r; // treat this as a repeat for now...
	     } else if(qs1 <= qs2 && qe1 <= qe2){	// && se1 > se2
		r++; rpt[V+1]=r; // Tough one: but treat this as a repeat for now...
				// will likely be deleted when selecting over full domain.
		// Case 2a: overlap in both subject and query.
/***************************************************************************
// This was found for protein 84501246, a hypothetical Helicase.
    probably want to trim back and merge together at some point, but happens so seldom
	it doesn't seem worth the effort right now.

   (ss2<=ss1 && ss1<=se2) && (se1 > se2) && (qs1 <= qs2  && qe1 <= qe2)

           qs1            qe1
    -------##################-------------------- query
            \_:_            :\___
              : \           :    \
    -------------##################-------------------- subject1
              :  ss1        :   se1 	
              :             :
              qs2           :qe2
    ----------##################------------------------------ query
    ----------##################------------------------------ subject2
              ss2            se2

>gi|84501246|ref|ZP_00999451.1| {|392(530)|<Proteobacteria(B)>}hypothetical protein OB2597_12813 [Oceanicola batsensis HTCC2597]

        Score = 44 bits (100), Expect = 1.4e-08
        Identities = 28/127(22.0)

 QUERY: 42   LVERLARADVICATLSGAGSdrniKFDTVIIDEATQATEPSCLIPLVK-GKVILVGDHKQ 100
             L      + +  A            FD V+IDEA+Q T  + L  +++  +V++VGD KQ
 SBJCT: 1271 LKPCWMMSPLAVAQYLHGQF----AFDLVVIDEASQMTPENALGAIMRaRQVVIVGDTKQ 1326

 QUERY: 101  LPPLTVQYRMHPEISEFPS-NEFYPGDIGVITPYNGQVALLRQYLQEKGGLKELYSVDGL 159
             LPP     ++  +  E           + +       V  LR + + +       S   +
 SBJCT: 1327 LPPTNFFSKVLDDADEDEDvRTDSESILDMANTTFTPVRQLRWHYRSRHPNLIAISNKMV 1386

 QUERY: 160  QEKKKGVEVATVD 172
              + +  +  A  D
 SBJCT: 1387 YDGQLTIFPAARD 1399

        Score = 42 bits (95), Expect = 5.3e-08
        Identities = 11/59(18.6)

 QUERY: 156  VDGLQEKKKGVEVATVDGFQGREKDVIILSCVRSLSDPRRLNvALTRAKrGLIIVGNSKT 215
             +D                      D++++     ++    L  A+ RA+  ++IVG++K
 SBJCT: 1269 IDLKPCWMMSPLAVAQYLHGQFAFDLVVIDEASQMTPENALG-AIMRAR-QVVIVGDTKQ 1326

 QUERY: 216  L 216
             L
 SBJCT: 1327 L 1327

 ***************************************************************************/
	     } else if(qs1 >= qs2 && qe1 >= qe2){	// && (ss2<=ss1 && ss1<=se2)
/***************************************************************************
	(qs1 >= qs2 && qe1 >= qe2) && (ss2<=ss1 && ss1<=se2)

                qs1            qe1
    ------------##################-------------------- query
                :\             :  \_
                : \            :    \
    ---------------##################-------------------- subject1
                :  ss1         :   se1 	
                :              :
              qs2              :qe2
    ----------##################------------------------------ query
    ----------##################------------------------------ subject2
              ss2            se2

>gi|121594442|ref|YP_986338.1| {<Proteobacteria(B)>}chromosome segregation protein SMC [Acidovorax sp. JS42] >gi|120606522|gb|ABM42262.1| chromosome segregation protein SMC [Acidovorax sp. JS42]

ss2 <= ss1 <= se2 --> 655 <= 707 <= 871
se1 = 935
qs2 =314; qs1 =837; qe2 = 521
qe1 = 1057
        Score = 61 bits (147), Expect = 4.4e-13
        Identities = 42/220(19.1)

 QUERY: 838  RRRLDELEQKIQELRDDIEELQEKLKSLEREIKELKKELEELDQRLEEVEERIEELEKEL 897
             RR   E + +  EL+ +   L +  +      +++  +L E++ +L +++ER    E
 SBJCT: 708  RREASEAQGRAHELQVETLRLSQLAEQTRARSEQIGADLSEVEAQLGDLQERRVTAEARF 767

 QUERY: 898  EELKEEKEELEKELEELEKELEKLQKKIEKLMKKREALLKK-------IAELEAKIRELE 950
             EEL  +  + ++   +L   + + ++++ +  +++ AL ++          L+A+  EL
 SBJCT: 768  EELDMQLADSQERHAQLGDRVIEAERRLAECREQQRALERQaqeaqfsQRSLQARRAELA 827

 QUERY: 951  RKIRDRGALPQRD-DEAQYLCEELGVLPEEAFEKYKETSSKQLFKKlEKVNEELKKYSHK 1009
             R I       +   DE Q   +EL  L + A +   + +     ++ EK     +
 SBJCT: 828  RTIDTAATQARALhDERQRAQDELSRLSDAAAQGGLQRALNLKIER-EKDLAAQRSQYDD 886

 QUERY: 1010 LEVNKKALDQY-NSFTEQREELTKRREELDRSRKSIEELIEVLDQRKDEA 1058
             L    +A D+         E L  R  E     ++    +E       +A
 SBJCT: 887  LTARLRASDELrMQNERTLEPLRARITEFQLKEQAARLGLEQYSTLLADA 936

        Score = 55 bits (131), Expect = 3.3e-11
        Identities = 40/208(19.2)

 QUERY: 315  KDLQDELEGAEEQRERLEKELKELKEEIEEKEAELEELKPEYESLKEEEEELKAELAELE 374
                Q  L    +  E LEKEL+      EE    L   +  Y    +     + E +E +
 SBJCT: 656  DSEQSGLLARAQDIEHLEKELRAQALIAEEARTALVRAESLYADAAQRLVAARREASEAQ 715

 QUERY: 375  QK---RKELYAKQGRRFKSKEERDKWLKNEIKSLNKEINDKKEELAKLQEDIEDLENEI- 430
              +    +    +  +  +    R + +  ++  +  ++ D +E     +   E+L+ ++
 SBJCT: 716  GRaheLQVETLRLSQLAEQTRARSEQIGADLSEVEAQLGDLQERRVTAEARFEELDMQLa 775

 QUERY: 431  -----EANLEKEIEELDQDLEELKDELEELEKKLEELKKKRDELQDERKELWREENKLQS 485
                   A L   + E ++ L E +++   LE++ +E +  +  LQ  R EL R  +   +
 SBJCT: 776  dsqerHAQLGDRVIEAERRLAECREQQRALERQAQEAQFSQRSLQARRAELARTIDTAAT 835

 QUERY: 486  ELANLKEELERAEAQLRAMMGLARGRAAVERLKRNLP 522
             +   L +E +RA+ +L  +   A        L   +
 SBJCT: 836  QARALHDERQRAQDELSRLSDAAAQGGLQRALNLKIE 872

 ***************************************************************************/
		r++; rpt[V+1]=r; // treat this as a repeat for now...
	     } else {
		PutSeq(stderr,QuerySeqGSAP(sap1),A); PutSeq(stderr,SubjectSeqGSAP(sap1),A);
		PutSeq(stderr,QuerySeqGSAP(sap2),A); PutSeq(stderr,SubjectSeqGSAP(sap2),A);
		fprintf(stderr,"ss2 <= ss1 <= se2 --> %d <= %d <= %d\n",ss2,ss1,se2);
		fprintf(stderr,"se1 = %d\n",se1);
		fprintf(stderr,"qs2 =%d; qs1 =%d; qe2 = %d\n",qs2,qs1,qe2);
		fprintf(stderr,"qe1 = %d\n",qe1);
		PutGSeqAlign(stderr,sap1,60,A); PutGSeqAlign(stderr,sap2,60,A);
#if 1
		r++; rpt[V+1]=r; // Tough one: but treat this as a repeat for now...
#else
		print_error("EvolvingRmOverlapsGSAP: This should not happen...(1)");
#endif
	     }
	   } else if(qs1<=qs2 && qs2<=qe1){		// Case 1: Overlap in Query:
		Int4	len_ol=qe1-qs2+1;
		Int4	cost,least_cost,crosspoint=0,*tailscore,*headscore;
		if(ss1 >= ss2){
			PutGSeqAlign(stderr,sap1,60,A); PutGSeqAlign(stderr,sap2,60,A);
#if 1
			r++; rpt[V+1]=r; // Tough one: but treat this as a repeat for now...
#else
			print_error("EvolvingRmOverlapsGSAP: This should not happen...(2)");
#endif
		} else {		// ss1 < ss2
		  tailscore=RunningScoreSAP_L(sap1,mx,A); // to left (from right)
		  headscore=RunningScoreSAP_R(sap2,mx,A); // to right (from left)
		  e_type qE=QuerySeqGSAP(sap1);
		  e_type sE=SubjectSeqGSAP(sap1);
/***************************************************************************
              qs               qe				tail_sap=sap1
    ----------###################------------------------------ query
              ______________* <--xpt  				
                            L______________
    ---------------------##################-------------------- query
    ---------------------##################-------------------- subject
                         qs              qe 			head_sap=sap2

  Want to minimize the lost score points from removing the overlap.
 ***************************************************************************/
		  least_cost=INT4_MAX;
		  Int4 ix,ixm1;
		  char rqe1=AlphaChar(ResSeq(qe1+1,qE),A);
		  char rqs2=AlphaChar(ResSeq(qs2+1,qE),A);
		  for(ixm1=0,ix=1; ix <= LenSeq(qE); ixm1++,ix++){
			if(tailscore[ixm1] && headscore[ix]){
			   cost = tailscore[ix] + headscore[ixm1]; // regions discarded
			   if(cost < least_cost){ least_cost = cost; crosspoint=ix; }
			   if(fp)fprintf(fp,
				"sq1: %c%d..%c%d (%d) sq2: %c%d..%c%d (%d)-> cost = %d\n",
				AlphaChar(ResSeq(ix+1,qE),A),ix+1,rqe1,qe1+1,tailscore[ix],
				rqs2,qs2+1,AlphaChar(ResSeq(ixm1+1,qE),A),ixm1+1,
				headscore[ixm1],cost);
			}
		  }
		  if(fp) PutSeqInfo(fp,sE);
		  if(fp) fprintf(fp,"\n### minimum overlap cost == %d\n",least_cost);
		  if(least_cost < overlap_cutoff){
		     if(fp){
	 	        fprintf(fp,"********* crosspoint = %d *********\n",crosspoint+1);
	     	        PutGSeqAlign(fp,sap1,60,A); PutGSeqAlign(fp,sap2,60,A);
		     }
		     TrimOverlappingSAPs(sap1, sap2,crosspoint);
		     if(fp){
	     	        PutGSeqAlign(fp,sap1,60,A); PutGSeqAlign(fp,sap2,60,A);
	 	        fprintf(fp,"********* Trimmed at %d *********\n",crosspoint+1);
		     }
		  } else { r++; rpt[V+1]=r; }
		  free(tailscore); free(headscore);
		}
	   } else if(qs2<=qs1 && qs1<=qe2){	// Case 1b: not shown...
		   if(ss2 >= ss1){ r++; rpt[V+1]=r; }	// This is a repeat...
		   else {
			PutGSeqAlign(stderr,sap1,60,A); PutGSeqAlign(stderr,sap2,60,A);
#if 1
			r++; rpt[V+1]=r; // Tough one: but treat this as a repeat for now...
#else
			print_error("EvolvingRmOverlapsGSAP: This should not happen...(3)");
#endif
		   }
	   } else {		// no overlap to delete...but look for permutations
		if(ss1 < ss2 && qs1 > qs2 ) { r++; rpt[V+1]=r; }
		else if(ss1 > ss2 && qs1 < qs2 ) { r++; rpt[V+1]=r; }
		// fprintf(stderr,"DEBUG S\n");
	   }
	} return rpt;
}

char    *OperationFromGSAP(sap_typ sap,UInt4 *QS, UInt4 *SS,
		UInt4 *QE,UInt4 *SE)
{
        Int4	i,s,qs,ss,length,index,total;
	UInt4	lenQ,lenS;
	char	*operation;

        assert(sap && sap->segs);
	dsp_typ	dsp= sap->segs; 
        for(total=index=0; index<dsp->numseg; index++){
		total += dsp->lens[index];
	} NEW(operation,total+5,char); 
	operation[0]='E';
        for(lenQ=lenS=0,s=1,index=0; index<dsp->numseg; index++){
		qs = dsp->starts[2*index];
		ss = dsp->starts[2*index+1];
                length = dsp->lens[index];
		if(qs == -1){			// insertion in subject
			lenS+=length;
			for(i=0; i < length; i++){ operation[s] = 'i'; s++; }
		} else if(ss == -1){		// deletion in subject
			lenQ+=length;
			for(i=0; i < length; i++){ operation[s] = 'd'; s++; }
		} else {			// match states
			lenS+=length; lenQ+=length;
			for(i=0; i < length; i++){ operation[s] = 'm'; s++; }
		}
        } operation[s]='E';
	*QS = dsp->starts[0]; *SS = dsp->starts[1]; 
	*QE = *QS + lenQ - 1; *SE = *SS + lenS - 1;
	return operation; 
}

#define	MAX_ALLOWED_GHSPs	1000
#define DEBUG_PatchSAPs		0
#define DEBUG_PatchSAP_2	0

sap_typ	PatchSAPs(Int4 overlap_cutoff, Int4 **mtx,sap_typ raw_sap,Int4 num_seqs,
			e_type queryE, a_type AB)
{ return PatchSAPs(overlap_cutoff,mtx,raw_sap,queryE,AB); }

sap_typ	PatchSAPs(Int4 overlap_cutoff,Int4 **mtx,sap_typ raw_sap,e_type queryE,a_type AB)
// WARNING: Assumes that HSPs for a given sequence are clustered together!
// Other assumptions apply....
/*****************************************************************************
******************************************************************************/
{
	sap_typ	gsap,rtn_sap,last,head,tail,*SAP,sap,tmp_sap,tmp_sap0;
	e_type	qE,sE;
	Int4	r,h,hsp,sq,len,width=60,num_hsps,score,i,j,qs,qe,ss,se;
	double	evalue,bit_score;
	UInt4	start,end,*NumSAPs,QS,SS,QE,SE,ii;
	char    *operation=0;

	assert(mtx);
	// 1. Find all HSPs for each sequence...
	sap_typ	SAP_GRP[MAX_ALLOWED_GHSPs];
	UInt4	Qstart[MAX_ALLOWED_GHSPs];
	char		*Operation[MAX_ALLOWED_GHSPs];
	UInt4	Qend[MAX_ALLOWED_GHSPs];
	UInt4	Sstart[MAX_ALLOWED_GHSPs];
	UInt4	Send[MAX_ALLOWED_GHSPs];

	// A. Remove overlaps between gapped HSPs for each sequence.
	// WARNING: ASSUMES THAT HEAD IS THE QUERY SELF-ALIGNMENT!
	sE=SubjectSeqGSAP(raw_sap); qE=QuerySeqGSAP(raw_sap); 
	assert(sE == qE);
	if(raw_sap->next){
		tmp_sap=raw_sap->next; sE=SubjectSeqGSAP(raw_sap); 
		// assert(sE != qE);
	} else return raw_sap;

#if DEBUG_PatchSAPs
	fprintf(stderr,"/*************************************************/\n");
	fprintf(stderr,"/*********** Raw output start ************/\n");
	for(tmp_sap=raw_sap ; tmp_sap; tmp_sap=tmp_sap->next){
		PutGSeqAlign(stderr, tmp_sap, width, AB);
	}
	fprintf(stderr,"/************ Raw output end ************/\n");
#endif

	for(last=0,head=raw_sap ; head; head=tail){

	   // A. find the head and tail sap for current sequence...
	   Int4 sid = SeqI(SubjectSeqGSAP(head));
#if DEBUG_PatchSAP_2
	   fprintf(stderr,"/************* Find head and tail sap ************/\n");
	   fprintf(stderr,"\n>"); PutSeqInfo(stderr,SubjectSeqGSAP(head));
	   fprintf(stderr,"\n-------------- seq_id=%d --------------\n",sid);
	   PutGSeqAlign(stderr, head, width, AB);
#endif
	   for(tail=head->next; tail && sid == SeqI(SubjectSeqGSAP(tail)); tail=tail->next){
#if DEBUG_PatchSAP_2
	   	fprintf(stderr,"\n-------------- seq_id=%d --------------\n",
			SeqI(SubjectSeqGSAP(tail)));
		PutGSeqAlign(stderr, tail, width, AB);
#endif
	   } 
#if DEBUG_PatchSAP_2
	fprintf(stderr,"/************* Find head and tail sap ************/\n");
	fprintf(stderr,"/*************************************************/\n");
#endif

	   if(tail==head->next){		// only one hsp for this seq.
		if(last){ last->next = head; }	// attach to end of new list.
		else { rtn_sap=head; }
		last=head; head=tail; continue;
	   }			// NOTE: 'head' == 'last' sap

	   // B. Store relevant information for HSPs.
	   num_hsps=0;
	   for(hsp=1,gsap=head; gsap && gsap != tail; gsap=gsap->next,hsp++){
		if(hsp >= MAX_ALLOWED_GHSPs) print_error("PatchSAPs(): too many hsps!");
		RawInfoGSAP(&qs,&qe,&ss,&se,gsap);
		SAP_GRP[hsp]=gsap; num_hsps++;
		Qstart[hsp]=qs; Sstart[hsp]=ss; Qend[hsp]=qe; Send[hsp]=se;
	   }

	   Int4 *path=LongestPathGSAP(num_hsps,SAP_GRP,Qstart,Qend,Sstart,Send);
	   Int4 *rpt=EvolvingRmOverlapsGSAP(overlap_cutoff,SAP_GRP,path,num_hsps,mtx,AB);

	   // C. Sort HSPs by their order in the alignment and group by repeats
	   sap_typ	**SAP_RPT;
	   Int4		num_rpts,*Nr;

	   NEWP(SAP_RPT,num_hsps+2,sap_typ);
	   NEW(Nr,num_hsps+2,Int4);
	// printf("/************ Sort HSPs by order in seq. ************/\n");
	   for(num_rpts=r=0,i=h=1; path[h]; h++,i++){ 
		hsp=path[h]; 
		if(rpt[h] > r){
			r++; i=1; num_rpts++; assert(r == rpt[h]); 
	   		NEW(SAP_RPT[r],num_hsps+2,sap_typ);  // INSURE BUG FIX
	   		// NEW(SAP_RPT[r],num_hsps,sap_typ);  // INSURE BUG
		} Nr[r]++; SAP_RPT[r][i]=SAP_GRP[hsp];
#if DEBUG_PatchSAP_2
		fprintf(stderr,"hsp[%d]->repeat %d: %d (total: %d)\n",hsp,r,i,Nr[r]);
		PutGSeqAlign(stderr, SAP_GRP[hsp], width, AB);
#endif
	   } assert(h-1 == num_hsps);
	// printf("/*************************************************/\n");

    // E. for each repeat group with > 1 sap get operation arrays...
    for(tmp_sap=0,r=1; r<= num_rpts; r++){
	if(Nr[r] > 1){  // then merge HSPs into one large HSP using operation information.
#if DEBUG_PatchSAPs
	  fprintf(stderr,"/*************************************************/\n");
	  fprintf(stderr,"/************ RmOverlapsGSAP() output ************/\n");
	  for(i=1; i<= Nr[r]; i++){
		PutGSeqAlign(stderr,SAP_RPT[r][i], width, AB);
	  }
	  fprintf(stderr,"/************ RmOverlapsGSAP() output ************/\n");
#endif
	  for(hsp=1; hsp<= Nr[r]; hsp++){
	     gsap=SAP_RPT[r][hsp];
	     Operation[hsp]=OperationFromGSAP(SAP_RPT[r][hsp],&QS,&SS,&QE,&SE);
	     Qstart[hsp]=QS; Sstart[hsp]=SS; Qend[hsp]=QE; Send[hsp]=SE;
#if DEBUG_PatchSAPs
	     fprintf(stderr,"Q=%d-%d; S=%d-%d: %s\n",QS,QE,SS,SE,Operation[hsp]);
#endif
	  }

	  gsap=SAP_RPT[r][1];
	  sE=SubjectSeqGSAP(gsap); qE=QuerySeqGSAP(gsap); 
	  NEW(Operation[0],LenSeq(sE)+LenSeq(qE)+5,char);
	  evalue=gsap->evalue; score=gsap->score; bit_score=gsap->bit_score;
	  qs=Qstart[1]; ss=Sstart[1]; qe=Qend[1]; se=Send[1];
	  Operation[0][0]='E'; ii=1;
	  operation=Operation[1];
	  for(i=1; operation[i] != 'E'; i++){ Operation[0][ii]=operation[i]; ii++; }
	  for(h=2; h <=Nr[r]; h++){		// Create a new operation array
		gsap=SAP_RPT[r][h];
		if(gsap->evalue < evalue){
		   evalue=gsap->evalue; score=gsap->score; bit_score=gsap->bit_score;
		}
	   	QS=Qstart[h]; SS=Sstart[h]; QE=Qend[h]; SE=Send[h];
		for(qe++ ; qe < QS; qe++){ Operation[0][ii]='d'; ii++; } // subject deletions
		for(se++ ; se < SS; se++){ Operation[0][ii]='i'; ii++; } // subject insertions

		char    *operation=Operation[h];
		for(i=1; operation[i] != 'E'; i++){ Operation[0][ii]=operation[i]; ii++; }
		qs=QS; ss=SS; qe=QE; se=SE;
	  } Operation[0][ii]='E';
	  QS=Qstart[1]; SS=Sstart[1]; 
#if DEBUG_PatchSAPs
	  fprintf(stderr,"NEW: Q=%d; S=%d: %s\n",QS,SS,Operation[0]);
#endif
	  assert(strlen(Operation[0]) <= (LenSeq(qE)-QS+LenSeq(sE)-SS+2));
	  tmp_sap0=ToGSeqAlign(strlen(Operation[0]),Operation[0],qE,sE,QS,SS);
	  // sap_typ newsap=AlignToGSAP(qE,sE,operation,QS,SS); // Alexandar's 
	  tmp_sap0->score=score;
	  tmp_sap0->bit_score=bit_score;
	  tmp_sap0->evalue=evalue;
	  tmp_sap0->next=0;
#if DEBUG_PatchSAPs
      	  PutGSeqAlign(stderr,tmp_sap0,width,AB);
#endif
	  for(hsp=1; hsp<= Nr[r]; hsp++) free(Operation[hsp]);
	} else {	// i.e., Nr[r] == 1
	  tmp_sap0=SAP_RPT[r][1];
	  tmp_sap0->next=0;
	}
	if(tmp_sap){ tmp_sap0->next=tmp_sap; tmp_sap=tmp_sap0; }
	else tmp_sap=tmp_sap0;
	
      } // end of for(r=1; r<= num_rpts; r++)...

      if(last) last->next = tmp_sap;	// attach new sap to end of current list
      else rtn_sap=tmp_sap;		// or assign as head for return.
	   
      last=tmp_sap;	// assign last to new sap.
      while(last->next) last=last->next;
     } return rtn_sap;
}

double	FractionOverlapSAP(sap_typ sap1, sap_typ sap2)
// Do sap1 and sap2 belong to the same sequence and do regions overlap?
{
        Int4            qs,qe,ss1,ss2,se1,se2;
        Int4            len,len1,len2;
        dsp_typ    dsp1,dsp2;

        dsp1 = sap1->segs; dsp2 = sap2->segs;
        if(dsp1->subject_id != dsp2->subject_id) return 0.0;
        len1=InfoGSAP(&qs,&qe,&ss1,&se1,sap1);
        len2=InfoGSAP(&qs,&qe,&ss2,&se2,sap2);
	len=MINIMUM(Int4,len1,len2);  // longest overlap possible.
        if(ss1 <= ss2){		// start of first is before start of second.
		return (double)(se1-ss2+1)/(double)(len);
					// e.g., se1=100 & ss2 = 150 --> -49 overlap. 
					// e.g., se1=150 & ss2 = 100 --> +51 overlap
	} else if(ss2 < ss1){		//  start of second is before start of first.
		return (double)(se2-ss1+1)/(double)(len);
					// e.g., se2 = 150 & ss1 = 100 --> +51 overlap.
	} else return 0.0;
}

BooLean	AtLeastOneLabeledSAPS(sap_typ sap)
{
        while(sap != NULL){
		if(LabelGSAP(sap)=='U') return TRUE;
		sap=sap->next;
        } return FALSE;
}

void	RemoveRejectsSAP(sap_typ *HEAD, Int4 num_sap,a_type AB)
{
	Int4	s;
	sap_typ	tail,head,gsap,tmp_sap;
	e_type	sE,qE;

	for(s=1; s <= num_sap; s++){
	  if(HEAD[s]){
#if 0
fprintf(stderr,"************* DEBUG1 *****************\n");
PutGSeqAlignList(stderr, HEAD[s], 60, AB);
fprintf(stderr,"************* DEBUG *****************\n");
#endif
	    qE=QuerySeqGSAP(HEAD[s]);
	    tail=head=MakeSelfAlignedGSAP(qE); tail->next=0; 
	    // SetLabelGSAP('U',HEAD[s]); // keep HEAD[s] unlabeled...
	    SetLabelGSAP('U',head);
	    for(gsap=HEAD[s]; gsap!=NULL; ){
		if(LabelGSAP(gsap) == 'U'){
			{ tail->next=gsap; tail=gsap; }
			gsap=gsap->next;
		} else { tail->next=0; tmp_sap=gsap->next; GSeqAlignFree(gsap); gsap=tmp_sap; }
	    } HEAD[s]=head;
#if 0
if(head && head->next){
  fprintf(stderr,"************* File %d: *****************\n",s);
  PutGSeqAlignList(stderr, head->next, 60, AB);
  fprintf(stderr,"****************************************\n");
}
#endif
	  }
	}
}

Int4	LabelBestSAP(sap_typ *sap,UInt4 start,UInt4 N,a_type AB)
// Assume all have same sE; simple procedure to remove overlapping hits:
// remove poorest scoring hits until there are no more overlaps.
// This was added because overlapping pairs of EF-hand domains were being clustered 
// together into a single disjoint set and only the best was saved (with old routine
// LabelBestSAP0() below).
// ---[___]--[___]---			sap1 (Ef-hand pair)
//        ---[___]--[___]---			sap2
//               ---[___]--[___]---			sap3
//                      ---[___]--[___]---			sap4
{
	Int4	i,j,set_i,set_j,end,I,J;
	e_type	sE,qE;
	Int4	bestscore,score,best,score_I, score_J;

	assert(N>0);
	if(N==1){
	   // fprintf(stderr,"*************** single hit %d(%d) ****************\n",N,start);
	   SetLabelGSAP('U',sap[start]); return 0; // use this one.
	}
	// fprintf(stderr,"*************** multiple hit %d(%d) ****************\n",N,start);
	end=start+N-1;
	sE=SubjectSeqGSAP(sap[start]);
	for(i=start; i <= end; i++){ 
		assert(sE==SubjectSeqGSAP(sap[i]));
		SetLabelGSAP('U',sap[i]);   // label all candidate saps
	}
	// Int4	len_I,len_J,s,e;
	for(i=1,I=start; I < end; i++,I++) {
	    if(LabelGSAP(sap[I]) == 'R') continue;
	    score_I = ScoreGSAP(sap[I]);
	    // len_I=GetStartEndGSAP(sap[I],&s,&e);
	    for(j=i+1,J=I+1; J <= end; j++,J++) {
		if(LabelGSAP(sap[J]) == 'R') continue;
		// qE=QuerySeqGSAP(sap0);
		double overlap=FractionOverlapSAP(sap[I],sap[J]);
		// fprintf(stderr,"******************* overlap=%g ********************\n",overlap);
		if(overlap >= 0.33){	// 33% overlap?
	    		score_J = ScoreGSAP(sap[J]);
#if 0
			fprintf(stderr,"******************* profiles %d,%d *********************\n",
					I,J);
			PutSeqInfo2(stderr,SubjectSeqGSAP(sap[I]));
			PutGSeqAlign(stderr, sap[I], 60, AB);
			PutSeqInfo2(stderr,SubjectSeqGSAP(sap[J]));
			PutGSeqAlign(stderr, sap[J], 60, AB);
			fprintf(stderr,"*********************************************\n");
#endif
			if(score_I > score_J){
				SetLabelGSAP('R',sap[J]);
			} else if(score_I <= score_J){
				SetLabelGSAP('R',sap[I]);
			}
			// set_i=linkDSets(set_i,set_j,sets);
			// are i and j in the same set?
		}
	    }
	} return 1;
}

Int4	LabelBestSAP0(sap_typ *sap,UInt4 start,UInt4 N,a_type AB)
// Assume all have same sE;
{
	Int4	i,j,set_i,set_j,end,I,J;
	e_type	sE,qE;

	assert(N>0);
	if(N==1){
	   // fprintf(stderr,"*************** single hit %d(%d) ****************\n",N,start);
	   SetLabelGSAP('U',sap[start]); return 0; // use this one.
	}
	// fprintf(stderr,"*************** multiple hit %d(%d) ****************\n",N,start);
	end=start+N-1;
	sE=SubjectSeqGSAP(sap[start]);
	for(i=start; i <= end; i++){ 
		assert(sE==SubjectSeqGSAP(sap[i]));
		SetLabelGSAP('R',sap[i]);   // reject.
	}
	ds_type sets;
	sets = DSets(N);
	for(i=1,I=start; I < end; i++,I++) {
	    set_i=findDSets(i,sets); 
	    for(j=i+1,J=I+1; J <= end; j++,J++) {
	        set_j=findDSets(j,sets); 
		if(set_i==set_j) continue;
		// qE=QuerySeqGSAP(sap0);
		double overlap=FractionOverlapSAP(sap[I],sap[J]);
		// fprintf(stderr,"******************* overlap=%g ********************\n",overlap);
		if(overlap >= 0.33){	// 33% overlap?
#if 0
			fprintf(stderr,"******************* profiles %d,%d *********************\n",
					I,J);
			PutSeqInfo2(stderr,SubjectSeqGSAP(sap[I]));
			PutGSeqAlign(stderr, sap[I], 60, AB);
			PutSeqInfo2(stderr,SubjectSeqGSAP(sap[J]));
			PutGSeqAlign(stderr, sap[J], 60, AB);
			fprintf(stderr,"*********************************************\n");
#endif
			set_i=linkDSets(set_i,set_j,sets);
			// are i and j in the same set?
		}
	    }
	}
#if 1
	PutDSets(stderr,sets);
#endif
	Int4	bestscore,score,best;
	for(i=1,I=start; I <= end; i++,I++) {
	    set_i=findDSets(i,sets);
	    if(set_i==i){	// canonical member.
	      // fprintf(stderr,"************* Set %d: *****************\n",i);
	      bestscore=ScoreGSAP(sap[I]); best=I;
	      // fprintf(stderr,"bestscore=%d; best=%d\n",bestscore,best);
	      for(j=1,J=start; J <= end; j++,J++) {
		if(j==i) continue;
		set_j=findDSets(j,sets);
		if(set_i == set_j){
		   score=ScoreGSAP(sap[J]);
	           // fprintf(stderr,"score(%d)=%d\n",j,score);
		   if(score > bestscore){
			bestscore=score; best=J;
	           	// fprintf(stderr,"bestscore=%d; best=%d\n",bestscore,j);
		   }
		}
	      }
	      SetLabelGSAP('U',sap[best]);	// use this one.
	      fprintf(stderr,"============> Labeled sap %d to 'U'\n",best);
	    }
	} NilDSets(sets);
	return 1;
}

Int4	FindBestOverlappingSAPs(sap_typ *head, Int4 num_sap_lists,Int4 num_seqs,a_type AB)
{ return FindBestOverlappingSAPs(head, num_sap_lists,AB); }

Int4	FindBestOverlappingSAPs(sap_typ *head, Int4 num_sap_lists,a_type AB)
{
	sap_typ sap0,sap=0,*array;
	e_type	sE,qE;
	double	key;
	UInt4	len,start,end,id,jd,N;
	Int4	JJ,i,j,k,I;

	// 1. Set qE's ident = I: assumes all seqs on list have same qE.
	Int4	hpsz=0;
	for(I=1; I <= num_sap_lists; I++){ 
		// qE=QuerySeqGSAP(head[I]); EqSeqI(I,qE); 
		hpsz+=LengthListGSAP(head[I]);
	}

    // 2. Make a SAP array using seqid & start of gapped HSP as ordering key.
    sph_typ sH=MakeSapHeap(hpsz+2);
    for(JJ=0,I=1; I <= num_sap_lists; I++){ 
	for(sap0=head[I]; sap0; sap0=sap0->next){
		sE=SubjectSeqGSAP(sap0);
		len=GetStartEndGSAP(sap0, &start, &end);
		key=(double) SeqI(sE) + ((double)start/(double)(LenSeq(sE)*2));
#if 0
		qE=QuerySeqGSAP(sap0);
		Int4 min_len = LenSeq(qE)/2;
		if(len >= min_len) {	// keep only if covers half the domain or more.
			assert(InsertSapHeap(sap0,key,sH));
		} else JJ++;
#else		// got rid of half domain option for curated_srch!!!  (10/21/05)
		assert(InsertSapHeap(sap0,key,sH));
#endif
#if 0
fprintf(stderr,"************* DEBUG: SeqI=%d(%g) *****************\n",I,key);
if(sap0)PutGSeqAlign(stderr,sap0, 60, AB);
fprintf(stderr,"************* DEBUG *****************\n");
#endif
	}
    }
    NEW(array,hpsz+3,sap_typ);
    for(i=0; sap0=DelMinSapHeap(&key,&I,sH); ){
#if 0
fprintf(stderr,"************* DEBUG: SeqI=%d(%g) *****************\n",I,key);
if(sap0)PutGSeqAlign(stderr,sap0, 60, AB);
fprintf(stderr,"************* DEBUG *****************\n");
#endif
	i++; array[i]=sap0; 
    } 
    assert((i+JJ)==hpsz);

    // 4. Make an array of SAPs in order of HSP occurances...
    for(start=N=jd=0,i=1; array[i]; i++){
	sE=SubjectSeqGSAP(array[i]);
	id=SeqI(sE);
	if(id == jd){	// 
		N++; 
	} else {
		// Process saps array for another sequence...
		if(start > 0) LabelBestSAP(array,start,N,AB);
		N=1; start=i; jd=id; 
	}
    } if(start > 0) LabelBestSAP(array,start,N,AB);
    NilSapHeap(sH); free(array);
}

Int4	LabelBestSAPs(sap_typ *HEAD, Int4 num_sap, Int4 num_seqs, 
	BooLean get_all_hits,a_type AB, UInt4 *NUM_EACH_CMA)
// merge the SAPs into one list
// Move this to my_ncbi because of private information.
// ScoreGSAP(sap); EvalueGSAP(sap); NextGSAP(sap); 
// SubjectSeqGSAP(sap); SeqI(E); SubjectSeqGSAP(sap);
// sap_typ ConcatenateGSAP(sap_typ head, sap_typ sap);
{
	sap_typ	head,tmp_sap,*BEST_SAP,tail,sap,gsap;
	Int4	len,score,sq,s,*BEST_CMA,*BEST_SCORE,min_len;
	UInt4	J,start,end,n;
	e_type	qE,sE;

	NEW(BEST_SAP,num_seqs+3,sap_typ);
	NEW(BEST_CMA,num_seqs+3,Int4);
	NEW(BEST_SCORE,num_seqs+3,Int4);
	for(s=1; s <= num_sap; s++){
	    n=0;
	    for(gsap=HEAD[s]; gsap!=NULL; gsap=gsap->next){
		SetLabelGSAP('R',gsap);	// reject.
		sE=SubjectSeqGSAP(gsap); sq=SeqI(sE); score=ScoreGSAP(gsap);
		if(sq > num_seqs) print_error("sequence id error in LabelBestSAPs()");
		len = GetStartEndGSAP(gsap, &start,&end);
#if 0	// debug...
		if(start==52 && end == 125){
			PutGSeqAlignList(stderr, gsap, 40, AB);
			fprintf(stderr,"len = %d; score = %d; best = %d\n",
					len,score,BEST_SCORE[sq]);
		}
#endif
		qE=QuerySeqGSAP(gsap);
		min_len = LenSeq(qE)/2;
		if(len < min_len) continue;
		if(score > BEST_SCORE[sq]){ 
		   BEST_SCORE[sq]=score; BEST_CMA[sq]=s; BEST_SAP[sq]=gsap; 
		} 
#if 0	// DEBUG...
		fprintf(stderr,"len = %d; score = %d; best = %d\n",
					len,score,BEST_SCORE[sq]);
#endif
	    } 
	}
	if(get_all_hits){
	   fprintf(stderr,"DEBUG\n");
	   for(sq=1; sq <= num_seqs; sq++){
	     if(BEST_SAP[sq]){	// a best sap exists for sq.
	     	s=BEST_CMA[sq];
	        for(gsap=HEAD[s]; gsap!=NULL; gsap=gsap->next){
		   sE=SubjectSeqGSAP(gsap); 
		   if(SeqI(sE) == sq){ 
		   	qE=QuerySeqGSAP(gsap);
		   	min_len = LenSeq(qE)/2;
		   	len = GetStartEndGSAP(gsap, &start,&end);
		   	if(len >= min_len) {
			    SetLabelGSAP('U',gsap); NUM_EACH_CMA[s]++; 
			}
		   }
		}
	     }
	   }
	} else {
	  for(sq=1; sq <= num_seqs; sq++){
		if(BEST_SAP[sq]){ 
		    qE=QuerySeqGSAP(BEST_SAP[sq]);
		    min_len = LenSeq(qE)/2;
		    len = GetStartEndGSAP(BEST_SAP[sq], &start,&end);
		    if(len >= min_len) {
			SetLabelGSAP('U',BEST_SAP[sq]); // Use this sap...
			NUM_EACH_CMA[BEST_CMA[sq]]++;
		    }
		}
	  }
	} // At this point all saps are appropriately labeled... Now delete unusable saps.

	RemoveRejectsSAP(HEAD, num_sap,AB);
	free(BEST_SAP); free(BEST_SCORE);
	return 0;
}

