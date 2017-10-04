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

#include "psc_typ.h"

#define MIN_angleDHA 90

Int4	psc_typ::ResToHeteroHBonds(FILE *fp,res_typ Res, float dmax, char color,
		unsigned short file_id, pdb_typ P)
// Identify all H-bonds to heteratoms.
{
	Int4	C2,n,a,a2,res;
	atm_typ	H,Donor,Accept;
	BooLean	na;
	rhb_typ *rhb=0,*rhb0=0;

	for(n=0,a = 1; a <= ResidueAtomNumber(Res); a++){
	     if(!ResidueAtomHydrogens(a,Res)) continue;	// No hydrogens in res atom; can't be a donor.
	     Donor = AtomResidue(a,Res);
	     res=ResAtom(Donor);	// residue number.
	     // For all shared hydrogen atoms...
	     for(Int4 h=1; (H=ResidueAtomHydrogen(a,h,Res)); h++){
	       for(C2 = 1; C2 <= P->nchains; C2++){
		if(C2 == ResidueChain(Res)) continue; 
		for(a2 = 1; a2 <= P->natoms[C2]; a2++){
	          Accept = P->atom[C2][a2];
		  if(HydrogenAtom(Accept)) continue;
		  if(!IsHeteroAtom(Accept)) continue;
		  if(IsWaterAtom(Accept)) continue;
		  if(Accept == Donor) continue;
		  double d_HA = DistanceAtoms(H,Accept);
		  if(d_HA > (double) dmax) continue; 
		  double d_DA = DistanceAtoms(Donor,Accept);
		  if(d_DA <= d_HA) continue;
		  double angle_DHA = CalcAngle(Donor,H,Accept);
		  if(angle_DHA >= MIN_angleDHA){
		    if(CarbonAtom(Donor)) NumWeakHBonds++; else NumModHBonds++;
		    if(rhb0 == 0){
			rhb0 = new rhb_typ(H,Donor,Accept,d_DA,d_HA,angle_DHA);
			rhb = rhb0; assert(rhb); // assert(rhb->next==0);
		    } else {
			// assert(rhb->next==0);
			rhb->next = new rhb_typ(H,Donor,Accept,d_DA,d_HA,angle_DHA);
			rhb = rhb->next;
		    } n++;
		  } 
		}
	       }
	     }
	}
	// use Residue as acceptor...
	if(rhb0) {
		if(fp) rhb0->PutAll(fp,color,file_id,P->A,P->nA);
		delete rhb0;
	}
	return n;
}

Int4	psc_typ::HeteroToResHBonds(FILE *fp,res_typ ResA,float dmax, char color,
			unsigned short file_id, res_typ **ResALL_I,Int4 *num_resALL_I, pdb_typ P)
/** return the number of residues in chain C2 within distance dmax of atoms in 
residue res0 of chain C **/
{
	Int4	n,a,a2,res,i,j,C;
	atm_typ	H,Donor,Accept;
	BooLean	na;
	rhb_typ *rhb=0,*rhb0=0;

    for(C=1; C <= P->nchains; C++){
      res_typ *ResD=ResALL_I[C];
      for(n=0,i=1; i <= num_resALL_I[C]; i++){
#if 0	// DEBUG
	Donor = AtomResidue(1,ResD[i]);
	fprintf(stderr,"============== Residue %d: ===============\n",ResAtom(Donor));
	PrintResidueAtoms(stderr,ResD[i]);
	continue;
#endif
        for(a=1; a <= ResidueAtomNumber(ResD[i]); a++){
	   Donor = AtomResidue(a,ResD[i]);
	   if(AminoAcidAtom(Donor)) continue; // only want HeteroAtom donors.
	   // print_error("Found heteroresidue");
	   if(!ResidueAtomHydrogens(a,ResD[i])) continue;	// No hydrogens in res atom; can't be a donor.
	   // For all shared hydrogen atoms...
	   for(Int4 h=1; (H=ResidueAtomHydrogen(a,h,ResD[i])); h++){
	     for(a2=1; a2 <= ResidueAtomNumber(ResA); a2++){
	          Accept = AtomResidue(a2,ResA);
		  if(HydrogenAtom(Accept)) continue;
		  if(IsHeteroAtom(Accept)) continue;
		  if(Accept == Donor) continue;
#if 0
		  double d_HA = DistanceAtoms(H,Accept);
		  double angle_DHA = CalcAngle(Donor,H,Accept);
		  double d_DA = DistanceAtoms(Donor,Accept);
		  if(d_HA <= (double) dmax && d_DA > d_HA && angle_DHA >= MIN_angleDHA){
#else
		  double d_HA = DistanceAtoms(H,Accept);
		  if(d_HA > (double) dmax) continue; 
		  double d_DA = DistanceAtoms(Donor,Accept);
		  if(d_DA <= d_HA) continue;
		  double angle_DHA = CalcAngle(Donor,H,Accept);
		  if(angle_DHA >= MIN_angleDHA){
#endif
		    if(CarbonAtom(Donor)) NumWeakHBonds++; else NumModHBonds++;
		    if(rhb0 == 0){
			rhb0 = new rhb_typ(H,Donor,Accept,d_DA,d_HA,angle_DHA);
			rhb = rhb0; assert(rhb); // assert(rhb->next==0);
		    } else {
			// assert(rhb->next==0);
			rhb->next = new rhb_typ(H,Donor,Accept,d_DA,d_HA,angle_DHA);
			rhb = rhb->next;
		    } n++;
		  } 
		}
	     }
	}
      }
    }
    // use Residue as acceptor...
    if(rhb0) { if(fp) rhb0->PutAll(fp,color,file_id,P->A,P->nA); delete rhb0; }
    return n;
}

// This is new version that does not depend on a given row...
void	psc_typ::FindBestPatternContacts2(char mode,FILE *logfp)
//**************************** look at structural features of patterns **********************
{
	Int4	x,y,i,j,n,os,k1,k2,cmaSeqID,r,N,I,C,R;
	atm_typ	a1,a2;
	double	d,m;
	Int4	file_id=1; // for vsi output file...
	Int4	posI,posJ;
	// FILE	*fp=stderr;
	FILE	*fp=0;
	
        // h_type HGH=Histogram("Contacts between co-conserved residues",0,500,5.0);
        h_type HGH=Histogram("Contact score between co-conserved residues",0,500,5.0);
        h_type HGK=Histogram("Contacts with kth pattern residue",0,500,10.0);
        if(dmax == 0.0) dmax = 4.5;
        if(HA_dmax == 0.0) HA_dmax = 2.5;

	// Int4	*K1Hits; NEW(K1Hits,ptrn->MaxPttrnPos + 3,Int4);

	Int4	**K1Hits; NEWP(K1Hits,nAlpha(AB) + 3,Int4);
	for(r=1; r <= nAlpha(AB); r++) NEW(K1Hits[r],ptrn->MaxPttrnPos + 3,Int4);

	// char	**K1Pttrn; NEWP(K1Pttrn,ptrn->MaxPttrnPos + 3,char);
	// Int4	*K1Pos; NEW(K1Pos,MaxPttrnPos + 3,Int4);

// move the below to init();
	NEWPP(PatternResHbonds,mpdb->NumPDB+3,Int4);
	NEWPP(WeakHBonds,mpdb->NumPDB+3,Int4);
	NEWPP(ModHBonds,mpdb->NumPDB+3,Int4);
	for(I=1; I <=mpdb->NumPDB; I++){
	     NEWP(PatternResHbonds[I],nChainsPDB(mpdb->pdb[I])+3,Int4);
	     NEWP(WeakHBonds[I],nChainsPDB(mpdb->pdb[I])+3,Int4);
	     NEWP(ModHBonds[I],nChainsPDB(mpdb->pdb[I])+3,Int4);
             for(C=1; C <= nChainsPDB(mpdb->pdb[I]); C++){
		NEW(PatternResHbonds[I][C],MAX_NUMBER_INTERNAL_RPTS+3,Int4);
		NEW(WeakHBonds[I][C],MAX_NUMBER_INTERNAL_RPTS+3,Int4);
		NEW(ModHBonds[I][C],MAX_NUMBER_INTERNAL_RPTS+3,Int4);
	     }
	}

    N=ptrn->MaxPttrnPos+3;
    Int4  K_hetero=ptrn->MaxPttrnPos+1;	// hetero node...
    Int4    hp_size=100;
    rih_typ rih(hp_size);	// 
    for(x=ptrn->NumPttrns; x > 0; x--){			// search over all patterns.
      Int4 Col=ptrn->PttrnCategory[x];
      for(posI=1; posI <= ptrn->NumPttrnRes[x]; posI++){	// search over all pattern positions.
        k1 = ptrn->PosPttrn[x][posI]; 
	assert(k1 <= ptrn->MaxPttrnPos);
	for(y=ptrn->NumPttrns; y >= x; y--) {		// search over all pairs of patterns. 
           for(posJ=1; posJ <= ptrn->NumPttrnRes[y]; posJ++){	// search over all pairs of pattern positions.
	     if(x == y && posJ <= posI) continue;	//  do not compare twice (once each way).
             k2 = ptrn->PosPttrn[y][posJ];
	     assert(k2 <= ptrn->MaxPttrnPos);
	     Int4 TotalHits=0,Hits,Total=0,Success=0,Fail=0;
	     for(I=1; I <=mpdb->NumPDB; I++){		// look over all pdb structures.
               for(C=1; C <= nChainsPDB(mpdb->pdb[I]); C++){
        	 for(R=1; R <= NumRpts[I][C]; R++){
		  Int4 Row = RelevantSeq[I][C][R];
		  if(p2c0->hpt->Cell(Row, Col) != '+') continue;	// foreground partition?
		  // if(RelevantSeq[I][C][R] != cmaFileID) continue; 
	  	  if(Col2pdbSeq[I][C][R]==0) continue;	// no sequence.
		  if(!matches[I][C][R][k1]) continue;	// skip non-matching residues.
	          if(!matches[I][C][R][k2]) continue;	// skip non-matching residues.
		  char chain = ChainCharPDB(C,mpdb->pdb[I]);
	          i = Col2pdbSeq[I][C][R][k1];
	          j = Col2pdbSeq[I][C][R][k2];
	          if(i == 0 || j == 0) continue;	// if no residue at these positions (i.e., deletions).
#if 1		  // look for distant hits only!
		  if(abs(i-j) < 4) continue;		// skip sequence adjacent positions.
#endif
	    	  char color='M';
		  os= OffSetSeq(mpdb->pdbSeq[I][C]);
#if 1	// see whether residue is in pattern: is this redundant with matches[][][] above?!?!?
		  Int4 r1=ResSeq(i,mpdb->pdbSeq[I][C]);	  
		  char res = AlphaChar(r1,AB);
		  if(strchr(ptrn->PttrnRes[x][posI],res) ==0) continue;  // residue not in pattern.
		  Int4 r2=ResSeq(j,mpdb->pdbSeq[I][C]);
		  res = AlphaChar(r2,AB);
		  if(strchr(ptrn->PttrnRes[y][posJ],res) ==0) continue;  // residue not in pattern.
#endif
		  Total++;
#if 0	// debug
		  NumWeakHBonds=NumModHBonds=0;
	    	  Hits = FindResResHBonds(stderr,file_id, i+os,j+os,HA_dmax,dmax,C,chain,color,
				mpdb->ResALL[I],mpdb->num_resALL[I],mpdb->pdb[I]);
#endif
		  NumWeakHBonds=NumModHBonds=0;
	    	  Hits = FindResResHBonds(fp,file_id, i+os,j+os,HA_dmax,dmax,C,chain,color,
				mpdb->ResALL[I],mpdb->num_resALL[I],mpdb->pdb[I]);
		  WeakHBonds[I][C][R] += NumWeakHBonds; ModHBonds[I][C][R] += NumModHBonds;
		  if(Hits > 0) Success++; else Fail++;
		  PatternResHbonds[I][C][R]+=Hits;
		  TotalHits+=Hits;
#if 1	// see how many hits each type of residue at each position has.
		  K1Hits[r1][k1] += Hits; 
		  K1Hits[r2][k2] += Hits;
#endif
		}
               }
	     }
	     if(TotalHits > 0){
		// K1Hits[k1] += TotalHits; K1Hits[k2] += TotalHits;
		// if(K1Pttrn[k1] == 0){ K1Pttrn[k1] = ptrn->PttrnRes[x][posI]; }
		// if(K1Pttrn[k2] == 0) K1Pttrn[k2] = ptrn->PttrnRes[y][posJ];
		// double Score=100* (double) TotalHits/ (double) Total;
		// double d=(double) (Success+1.0) / (double)(Fail+1.0);
		double d=(double) (Success) / (double)(Total);	// penalize for poor sucess rate.
		double Score=(double) TotalHits*d;	// strong bonds count more...
		IncdHist(Score,HGH); 	// all hits for k1 vs k2.
#if 1
		char	Lineage=' ';
		Int4	*P;
		if(p2c0->hpt->IsTree(P)){
		    Int4 z,p,rowX,rowY;
#if 0
		    rowX = ptrn->NumPttrns - ptrn->PttrnCategory[x] + 1; // reversed order for *.pttrns file.
		    rowY = ptrn->NumPttrns - ptrn->PttrnCategory[y] + 1; // reversed order for *.pttrns file.
		    for(p=P[rowX]; p != 0; p = P[p]){ if(p == rowY){ Lineage='L'; break; } } 
		    if(Lineage == ' '){
		       for(p=P[rowY]; p != 0; p = P[p]){ if(p == rowX){ Lineage='L'; break; } }
#else
		    rowX = ptrn->PttrnCategory[x]; rowY=ptrn->PttrnCategory[y];
		    for(p=P[rowX]; p != 0; p = P[p]){ if(p == rowY){ Lineage='L'; break; } } 
		    if(Lineage == ' '){
		       for(p=P[rowY]; p != 0; p = P[p]){ if(p == rowX){ Lineage='L'; break; } }
#endif
		    } free(P);
		}
		char	*Iname=this->GetRowName(x),*Jname=this->GetRowName(y);
	      	rpi_typ *rpi = new rpi_typ(Iname,Jname,k1,k2,Score,Total,ptrn->PttrnRes[x][posI],
				ptrn->PttrnRes[y][posJ],0,Lineage);
#elif 1
		char	*Iname=this->GetRowName(x),*Jname=this->GetRowName(y);
	      	rpi_typ *rpi = new rpi_typ(Iname,Jname,k1,k2,Score,Total,ptrn->PttrnRes[x][posI],ptrn->PttrnRes[y][posJ]);
#else
	      	rpi_typ *rpi = new rpi_typ(k1,k2,Score,Total,ptrn->PttrnRes[x][posI],ptrn->PttrnRes[y][posJ]);
#endif		
		if(rih.Insert(rpi,-Score) == 0) delete rpi;;
	     // fprintf(stdout,"k1=%d; k2=%d; TotalHits=%d\n",k1,k2,TotalHits);
	     }
	   }
	}
	for(I=1; I <=mpdb->NumPDB; I++){
           for(C=1; C <= nChainsPDB(mpdb->pdb[I]); C++){
             for(R=1; R <= NumRpts[I][C]; R++){
	      Int4 Row = RelevantSeq[I][C][R];
	      if(p2c0->hpt->Cell(Row, Col) != '+') continue;
	      // if(RelevantSeq[I][C][R] != cmaFileID) continue; 
	      if(Col2pdbSeq[I][C][R]==0) continue;
	      if(!matches[I][C][R][k1]) continue;	// skip non-matching residues.
	      i = Col2pdbSeq[I][C][R][k1]; os= OffSetSeq(mpdb->pdbSeq[I][C]);
	      Int4 res=i+os;
	      for(r = 1; r <= mpdb->num_resALL[I][C]; r++){
		res_typ Res=mpdb->ResALL[I][C][r];
		if(ResidueID(Res) == res){
	    	   char color='M';
		   NumWeakHBonds=NumModHBonds=0;
	   	   Int4 hitsHR=HeteroToResHBonds(fp,Res,dmax,color,file_id,mpdb->ResALL[I],
					mpdb->num_resALL[I],mpdb->pdb[I]);
		   PatternResHbonds[I][C][R]+=hitsHR;
		   Int4 hitsRH=ResToHeteroHBonds(fp,Res,dmax,color,file_id,mpdb->pdb[I]);
		   PatternResHbonds[I][C][R]+=hitsRH;
		   WeakHBonds[I][C][R] += NumWeakHBonds; ModHBonds[I][C][R] += NumModHBonds;
    		   // if(hitsHR > 0 || hitsRH > 0){ grf.AddEdge(k1,K_hetero); }
		   // K1Hits[k1] += hitsRH + hitsHR;
		   // if(K1Pttrn[k1] == 0) K1Pttrn[k1] = ptrn->PttrnRes[x][posI];
		}
	      }
	    }
	   }
	}
	// fprintf(stdout,"%s%d: %d HeteroToRes; %d ResToHetero\n",ptrn->PttrnRes[x][posI],k1,Hits1,Hits2);
      }
    } // end of for(x=ptrn->NumPttrns; x > 0; x--)... loop
    //*************************************************************************
    // Report the number of Pattern residue Hbonds for each pdb file and chain.
#if 0
    for(I=1; I <=mpdb->NumPDB; I++){
           for(C=1; C <= nChainsPDB(mpdb->pdb[I]); C++){
	        if(Col2pdbSeq[I][C]==0) continue;
	        if(C==1) fprintf(stdout,"pdb[%d]: %s\n",I,mpdb->pdb_file[I]);
		fprintf(stdout,"   chain %d: %d hits\n",C,PatternResHbonds[I][C]);
	   }
    } fprintf(stdout,"\n");
#else	// Output to file.
    assert(cmafilename);
  //================== for hierview program... ======================
  //================== for hierview program... ======================
  //================== for hierview program... ======================
  if(verbose==FALSE){
    FILE *sfp=0; 
    BooLean	opened=FALSE;
    Int4 lastI=0,ii,vsi_number=0,best_j=0,best_score=-999,score,best_row=0,best_set=0,best_vsi=0;
    for(ii=0,i=1; i <= NumPDB_Sets; i++){
	if(CardSet(RelevantSet[i]) == 0) continue;
	assert(PDB_SetI[i]);
	vsi_number++; assert(cmafilename);
	char Str[200],LastStr[200];  LastStr[0]=0;
	for(j=1; PDB_SetI[i][j]; j++){
	  I=PDB_SetI[i][j]; C=PDB_SetC[i][j];
	  char *pstr=strrchr(mpdb->pdb_file[I],'/');
	  if(pstr==NULL) pstr=mpdb->pdb_file[I]; else pstr++;
	  sprintf(Str,"%s",pstr);
	  // fprintf(stderr,"Str=%s; LastStr=%s\n",Str,LastStr);
	  if(ii==0){
	     if(logfp==0) sfp=open_file(cmafilename,"_pdb.log","w"); 
	     else opened=TRUE;
	     sprintf(LastStr,"%s",Str);
	  } if(LastStr[0] == 0){
	     sprintf(LastStr,"%s",Str); best_score=-999; best_j=0;
	  } if(strcmp(Str,LastStr) != 0){	// new pdb file found...
	     if(logfp) fprintf(logfp,"    ChnSARP %s_pdb.VSI %d %d\n",cmafilename,best_j,vsi_number);
	     else fprintf(sfp,"    ChnSARP %s_pdb.VSI %d %d\n",cmafilename,best_j,vsi_number);
	     sprintf(LastStr,"%s",Str); best_score=-999; best_j=0;
	  } ii++;
	  for(R=1; R <= NumRpts[I][C]; R++){
	   Int4 Row = RelevantSeq[I][C][R];
	   if(Col2pdbSeq[I][C]==0) continue;		// this probably can't happen...
	   char *pstr=strrchr(mpdb->pdb_file[I],'/');
	   if(pstr==NULL) pstr=mpdb->pdb_file[I]; else pstr++;
	   if(logfp) fprintf(logfp,"  File%d=%s;",j,pstr); 
	   else fprintf(sfp,"  File%d=%s;",j,pstr); 
	   score=(4*ModHBonds[I][C][R] + WeakHBonds[I][C][R]);
	   if(logfp){
		if(opened) fprintf(logfp," Rpt: %d; chain: %c; score= %d (w:%d;m:%d)\n",
		   R,ChainCharPDB(C,mpdb->pdb[I]), score, WeakHBonds[I][C][R],ModHBonds[I][C][R]);
	   } else if(sfp) fprintf(sfp," Rpt: %d; chain: %c; score= %d (w:%d;m:%d)\n",
		   R,ChainCharPDB(C,mpdb->pdb[I]), score, WeakHBonds[I][C][R],ModHBonds[I][C][R]);
	   if(score > best_score){
		best_score=score; best_j=j; best_row=Row; best_set=i; best_vsi=vsi_number;
	   }
	  }
	}
	if(best_j >= 0){	// need to output previous best...
	    if(logfp) fprintf(logfp,"    ChnSARP %s_pdb.VSI %d %d\n",cmafilename,best_j,vsi_number);
	    else fprintf(sfp,"    ChnSARP %s_pdb.VSI %d %d\n",cmafilename,best_j,vsi_number);
	     LastStr[0]=0; best_score=0; best_j=0; best_vsi=0;
	}
    } 
    if(sfp) fclose(sfp);
    //=================== end =======================
    //=================== end =======================
  } else {	// verbose==TRUE;

    FILE *sfp=open_file(cmafilename,"_pdb.eval","w");
    Int4 hpsz = 10000,*best_ii,*best_jj,*best_rr,*best_ss;
    mh_type mH=Mheap(hpsz,3);	
    h_type HGV=Histogram("structural interaction scores.",0,10000,50);
    NEW(best_ii,hpsz+5,Int4); NEW(best_jj,hpsz+5,Int4); NEW(best_rr,hpsz+5,Int4);
    NEW(best_ss,hpsz+5,Int4);
    Int4 lastI=0,ii;
    Int4 vsi_number=0;
    FILE *bfp=0;
    if(mode=='P') bfp=open_file(cmafilename,"_sed","w");
    for(ii=0,i=1; i <= NumPDB_Sets; i++){
	if(CardSet(RelevantSet[i]) == 0) continue;
	assert(PDB_SetI[i]);
	vsi_number++;
	// if(!MemberSet(cmaFileID,RelevantSet[i])) continue;
        assert(cmafilename); ii++;
        // fprintf(sfp,"========== %s_pdb%d.vsi ==========\n",cmafilename,i);
        // fprintf(sfp,"========== %s_pdb%d.vsi ==========\n",cmafilename,ii);
        fprintf(sfp,"========== %s_pdb%d.vsi (Set %d) ==========\n",
			cmafilename,vsi_number,i);
	// fprintf(sfp,"========== vsi file %d ==========\n",i);
    	Int4 best_j=0,best_score=0,score,best_row=0,best_set=0;
	for(j=1; PDB_SetI[i][j]; j++){
	  I=PDB_SetI[i][j]; C=PDB_SetC[i][j];
	  for(R=1; R <= NumRpts[I][C]; R++){
	   Int4 Row = RelevantSeq[I][C][R];
	   // if(p2c0->hpt->Cell(Row, Col) != '+') continue;
	   // if(RelevantSeq[I][C][R] != cmaFileID) continue; 
	   if(Col2pdbSeq[I][C]==0) continue;		// this probably can't happen...
	   char *pstr=strrchr(mpdb->pdb_file[I],'/');
	   if(pstr==NULL) pstr=mpdb->pdb_file[I]; else pstr++;
	   fprintf(sfp,"  File%d=%s;",j,pstr); 
	   score=(4*ModHBonds[I][C][R] + WeakHBonds[I][C][R]);
	   IncdHist(score,HGV);
	   fprintf(sfp," Rpt: %d; chain: %c; score= %d (w:%d;m:%d)\n",
			R,ChainCharPDB(C,mpdb->pdb[I]), score,
			WeakHBonds[I][C][R],ModHBonds[I][C][R]);
			// PatternResHbonds[I][C][R];
#if 1
	   if(bfp){ fprintf(bfp,"    sed_pml %s_pdb.VSI %d %d\n",cmafilename,j,i); }
#endif
	   if(score > best_score){ best_score=score; best_j=j; best_row=Row; best_set=i; }
	  }
	} fprintf(sfp,"\n");
	if(best_score > 0){
		fprintf(sfp,"  best(%d): ChnVSI %s_pdb%d.vsi %d\n\n",
					best_score,cmafilename,vsi_number,best_j);
		Int4 x=InsertMheap((keytyp)-best_score,mH);     // inserts as minkey
		if(x > 0){ 
			// best_ii[x] = ii; best_jj[x]=best_j; 
			// best_ii[x] = i; 
			best_ii[x] = vsi_number; 
			best_jj[x]=best_j; 
			best_rr[x]=best_row;
			best_ss[x]=best_set;
		}
	}
    } fclose(sfp);
    if(bfp) fclose(bfp);
    Int4 rank=0;
    bfp=open_file(cmafilename,"_pdb.best","w");
    while(!EmptyMheap(mH)){
	double	score=-MinKeyMheap(mH);
	rank++; x=DelMinMheap(mH);
	Int4 Row=best_rr[x];
	char str[50]; StrSeqID(str,30,FullSeq[best_ss[x]]);
	fprintf(bfp,"%3d: Score = %.0f; Row = %d (\"%s\"); pdb set=\"%s\".\n",
			rank,score,Row,p2c0->hpt->ElmntSetName(Row),str);
#if 0
	fprintf(bfp,"    ChnVSI %s_pdb%d.vsi %d\n", cmafilename,best_ii[x],best_jj[x]);
#else
	fprintf(bfp,"    ChnSARP %s_pdb.VSI %d %d\n", cmafilename,best_jj[x],best_ii[x]);
#endif
    }
    PutHist(bfp,60,HGV); NilHist(HGV);
    fclose(bfp); NilMheap(mH); free(best_ii); free(best_jj); free(best_rr); free(best_ss);
    if(!verbose) return;
#endif
#if 0	// Key statistics
     For each pattern residue:
	Variance of C_alpha distances at pattern positions (deviation from mean).
	Numbers of predicted strong H-bonds involving Pattern residues.
	Numbers of predicted weak H-bonds (e.g., CH and CH-pi) involving Pattern residues.
	Matches and mis-matches.
	Missing contacts in certain groups 

	I may want to compare over all proteins...Input multiple subgroups!?
	Get better estimate of structural variability and of deviant residues.
	Input all *mma alignment files and pdb files.
	Input all patterns and hpt? or else patterns for each subgroup.
#endif
    //*************************************************************************
#if 0
    fprintf(stdout,"\nBest positions: ");
    for(i=1; i <= ptrn->MaxPttrnPos; i++){
	if(BestPos[i]){
		fprintf(stdout,"%s%d,",K1Pttrn[i],i);
	}
    } fprintf(stdout,"\n\n");
    free(BestPos);
#endif

    //*************************************************************************
    grf_typ grf(N);	// add an extra node for HeteroMolecules.
    // output best pairwise residue-to-residue H-bonds.
    FILE *ofp=0; if(verbose) ofp=open_file(cmafilename,".interact","w");
    if(ofp) PutHist(ofp,60,HGH); // NilHist(HGH);
    rpi_typ *rpi;
    double	key;
    Int4	Item,Id=0;
    while((rpi=rih.DelMin(&key,&Item)) != NULL){
    		grf.AddEdge(rpi->ColumnI(),rpi->ColumnJ());
		Id++; if(ofp){ fprintf(ofp,"%d: ",Id); rpi->Put(ofp); } delete rpi; 
    } if(ofp){ fprintf(ofp,"\n\n"); fclose(ofp); }

    // report the cliques of interacting pattern residues.
    Int4 NumClust,MinClique=3;
    double pcut=0.001;
    set_typ *NodeSet=0; // Not used...
    unsigned short TrueN=K_hetero; // 
    FILE *gfp=0; if(verbose) gfp=open_file(cmafilename,".grph","w");
    if(gfp){ grf.Put(gfp); fflush(gfp); }
#if 1	// this was core dumping due to array overread; may not yet be fixed.
    vst_typ **clique=grf.Bron_Kerbosch(MinClique,100,&NumClust,pcut,NodeSet,TrueN);
    BooLean *BestPos; NEW(BestPos,ptrn->MaxPttrnPos + 3,BooLean);
    for(i=1; i <= NumClust; i++){
       for(j=0; j < clique[i]->Size(); j++){
	Int4    k1=clique[i]->Vertex(j);
	BestPos[k1]=TRUE;
       } if(gfp) clique[i]->Put(gfp); delete clique[i]; 
    } free(clique); if(gfp) fclose(gfp); free(BestPos);
#endif

    //*************************************************************************
	// output individual key pattern residues (those with the most interactions)).
	if(verbose) ofp=open_file(cmafilename,".key_res","w");
	Int4	HpSz=150,*ResX,*PosX,X;
	mH = Mheap(HpSz,4);
	NEW(ResX,HpSz +3, Int4); NEW(PosX,HpSz +3, Int4);
	
	for(i=1; i<=ptrn->MaxPttrnPos; i++){
	    for(r=1; r <= nAlpha(AB); r++){
		if(K1Hits[r][i] > 0){
		   X = InsertMheap(-(keytyp)K1Hits[r][i],mH);
		   IncdHist(K1Hits[r][i],HGK); 	// all hits for k1.
		   ResX[X]=r; PosX[X]=i; 
		}
	    }
	}
	if(ofp) PutHist(ofp,60,HGK); // NilHist(HGK);
	Item=1;
	while(!EmptyMheap(mH)){
		X=DelMinMheap(mH);
		i=PosX[X]; r=ResX[X];
		Int4 x,p;
		if(ofp) fprintf(ofp,"%d: %c%d: %d hits\n",Item++,AlphaChar(r,AB),i,K1Hits[r][i]);
		// fprintf(ofp,"%d: column %d: %d hits\n",Item++,i,K1Hits[i]);
	} if(ofp){ fprintf(ofp,"\n\n");  fclose(ofp); }
	NilMheap(mH); free(ResX); free(PosX);
   } // end of if(!verbose) .... else ...
	NilHist(HGK); NilHist(HGH);
	for(r=1; r <= nAlpha(AB); r++) free(K1Hits[r]); free(K1Hits);
}

Int4	psc_typ::FindHBonds(FILE *fp,res_typ Res,float dmax,res_typ Res2,
		char color, unsigned short file_id)
/** return the number of residues in chain C2 within distance dmax of atoms in 
residue res0 of chain C **/
{
	Int4	n,a,a2,res,res2;
	atm_typ	H,Donor,Accept;
	BooLean	na;
	rhb_typ *rhb=0,*rhb0=0;

#if 0
	FILE	*efp;
	BooLean	comments=FALSE;
	if(efp==stderr) comments=TRUE; else comments=FALSE;
	comments=TRUE; // waiting to debug this option;
	Int4 gly = AlphaCode('G',AB);
#endif
	for(n=0,a = 1; a <= ResidueAtomNumber(Res); a++){
	     if(!ResidueAtomHydrogens(a,Res)) continue;
	     Donor = AtomResidue(a,Res);
	     if(!AminoAcidAtom(Donor)) continue;
	     res=ResAtom(Donor);
	     for(Int4 h=1; (H=ResidueAtomHydrogen(a,h,Res)); h++){
	        for(a2 = 1; a2 <= ResidueAtomNumber(Res2); a2++){
	     	  Accept = AtomResidue(a2,Res2);
		  res2 = ResAtom(Accept);
		  if(HydrogenAtom(Accept)) continue;
	     	  if(!AminoAcidAtom(Accept)) continue;
                  // Skip backbone-to-backbone interactions and hetero backbone interactions (called elsewhere)...
		  if(!SideAtom(Donor) && !SideAtom(Accept)) continue;
		  if(ChainResidue(Res) == ChainResidue(Res2) && abs(res2 - res) < 2) continue;
		  // if(Accept == Donor) continue;   // this is not possible.

		  double d_HA = DistanceAtoms(H,Accept);
		  if(d_HA > (double) dmax) continue; 
		  double d_DA = DistanceAtoms(Donor,Accept);
		  if(d_DA <= d_HA) continue;
		  double angle_DHA = CalcAngle(Donor,H,Accept);
		  // if(d_HA <= 2.5 && d_DA > d_HA && angle_DHA >= 90){
		  // if(d_HA <= (double) dmax && d_DA > d_HA && angle_DHA >= MIN_angleDHA)
		  if(angle_DHA >= MIN_angleDHA){
		    if(CarbonAtom(Donor)) NumWeakHBonds++; else NumModHBonds++;
		    if(rhb0 == 0){
			rhb0 = new rhb_typ(H, Donor, Accept, d_DA, d_HA, angle_DHA);
			rhb = rhb0; assert(rhb);
		    } else {
			rhb->next = new rhb_typ(H, Donor, Accept, d_DA, d_HA, angle_DHA);
			rhb = rhb->next;
		    } n++;
		  } 
		}
	     }
	}
	if(rhb0) {
		if(fp) rhb0->PutAll(fp,color,file_id,AB,nAB);
		delete rhb0;
	}
	return n;
}

Int4	psc_typ::FindResResHBonds(FILE *fptr,Int4 file_id,Int4 ResI,Int4 ResJ,float HA_dmax,float dmax,
		Int4 C,char chain,char color,res_typ **ResALL_I,Int4 *num_resALL_I,pdb_typ P)
{
	Int4 TotalHits=0,hits;
	Int4 C2;
	Int4 r,i,j,num_resA,num_resD;
	if(chain) C=GetChainNumberPDB(P,chain); // C = query residue chain.
	if(C==0){ 
		     fprintf(stderr,"chain = %c\n",chain);
		     fprintf(stderr,"chain not found!");
		     return 0;
	}
	if(chain==0) chain='?';
	res_typ *ResCA,*ResOA,*ResCD,*ResOD;
	Int4	num_resC,num_resO;
	num_resC=num_resALL_I[C]; ResCA=ResCD=ResALL_I[C];

	for(C2=1; C2 <= nChainsPDB(P); C2++){ 
	      ResOD=ResOA=ResALL_I[C2];	// other than query residue.
	      num_resO=num_resALL_I[C2];
	      
	      for(i=1; i <= num_resC; i++){
		if(ResidueID(ResCD[i]) == ResI){	// --> ResCA[i] == res0
	         for(j=1; j <= num_resO; j++){
		   // PrintResidueAtoms(fptr,ResA[i]);
		   if(ResidueID(ResOA[j]) == ResJ){
			TotalHits += FindHBonds(fptr,ResCD[i],HA_dmax,ResOA[j],color,file_id);
			TotalHits += FindHBonds(fptr,ResOA[j],HA_dmax,ResCD[i],color,file_id);
		   }
		 }
		}
	      }
	      // fprintf(fptr,"//=========== pi-bonds: ===========\n");
	      for(i=1; i <= num_resC; i++){
		if(ResidueID(ResCD[i]) == ResI){	// --> ResCA[i] == res0
	         for(j=1; j <= num_resO; j++){
		   // PrintResidueAtoms(fptr,ResA[i]);
		   if(ResidueID(ResOA[j]) == ResJ){
		     hits = FindAromaticHbondsPDB(fptr,ResCD[i],ResOA[j],C,C2,dmax,color,P);
		     NumWeakHBonds+=hits; TotalHits += hits;
		     hits = FindAromaticHbondsPDB(fptr,ResOA[j],ResCD[i],C2,C,dmax,color,P);
		     NumWeakHBonds+=hits; TotalHits += hits;
		     TotalHits += FindOtherPiHbondsPDB(fptr,ResCD[i],ResOA[j],C,C2,dmax,color,P);
		     TotalHits += FindOtherPiHbondsPDB(fptr,ResOA[j],ResCD[i],C2,C,dmax,color,P);
		   }
		 }
		}
	      }
	} return (TotalHits);
} 

Int4	psc_typ::FindVDWsContacts(FILE *fptr,Int4 file_id,Int4 ResI,Int4 ResJ,float HA_dmax,float dmax,
		Int4 C,char chain,char color,res_typ **ResALL_I,Int4 *num_resALL_I,pdb_typ P)
{
	Int4 TotalHits=0,TotalContacts=0;
	Int4 C2;
	Int4 r,i,j,num_resA,num_resD;
	if(chain) C=GetChainNumberPDB(P,chain); // C = query residue chain.
	if(C==0){ 
		     fprintf(stderr,"chain = %c\n",chain);
		     fprintf(stderr,"chain not found!");
		     return 0;
	}
	if(chain==0) chain='?';
	// Int4 Start,End;
	res_typ *ResCA,*ResOA,*ResCD,*ResOD;
	Int4	num_resC,num_resO;
	num_resC=num_resALL_I[C]; ResCA=ResCD=ResALL_I[C];

	// Move this routine to libpdb.a at some point:
	a_type AB = AminoAcidAlphabetPDB(P);
	Int4	MaxNeighbors=50; // maximum residue neighbors to look at
	res_typ	*Contact;
	double	contact_surface,*Ratio;
	NEW(Contact,MaxNeighbors+3,res_typ);
	NEW(Ratio,MaxNeighbors+3,double);

	 Int4	hpsz=0;
	 dh_type dH = dheap(MaxNeighbors,4);

     	 //********** Check out surface VDW contacts of res0 to nearby residues.
	 // fprintf(fptr,"\n");
	 //  fprintf(fptr,"\n//*********** ResidueToResidueContacts (%d): ***********\n",res0);
	 for(j=1; j <= num_resALL_I[C]; j++){
	    if(ResidueID(ResCD[j]) == ResI) {
		double total_ratio = 0.0;
		BooLean	SideChainOnly;
		char c_res = GetCharResidue(ResCD[j],AB);
		if(c_res == 'G') SideChainOnly=FALSE;
		else SideChainOnly=TRUE;
		for(C2=1; C2 <= nChainsPDB(P); C2++){ 
		   res_typ *ResC2 = ResALL_I[C2];
		   if(num_resALL_I[C2] > 0 && !isspace(ChainResidue(ResC2[1])))
			   // fprintf(fptr,"//======== Chain %c ==========\n", ChainResidue(ResC2[1]));
	           for(hpsz=0,i=1; i <= num_resALL_I[C2]; i++){
		   	if(ResCD[j] == ResC2[i]) continue;
			Int4 res2=ResidueID(ResC2[i]);
			if(res2 != ResJ) continue;	// modified here...
			double ratio;
	     		contact_surface = ResidueToResidueContact(ResCD[j],ResC2[i],&ratio,SideChainOnly);
			if(contact_surface > 0.0){
				TotalContacts++; hpsz++;
				insrtHeap(hpsz,(keytyp)-contact_surface,dH);
				Ratio[hpsz]=ratio; total_ratio+= ratio;
				Contact[hpsz]=ResC2[i];
			}
		   }
		   Int4 hp_pos=1;
		   while(!emptyHeap(dH)){
			contact_surface = (double) -minkeyHeap(dH);
			assert((i=delminHeap(dH)) != 0);
			char    c2_res = GetCharResidue(Contact[i],AB);
			if(!isspace(ChainResidue(Contact[i]))){
			   if(fptr) fprintf(fptr, "#%c%d%c.X  // %c%d%c %.1f A^2 contact (%1.f%%:%d).\n",
					c2_res,ResidueID(Contact[i]),ChainResidue(Contact[i]),
					c_res,ResidueID(ResCD[j]),ChainResidue(ResCD[j]),
					contact_surface,Ratio[i]*100.0,hp_pos);
			} else {
			   if(fptr) fprintf(fptr,"#%c%d.X  // %c%d %.1f A^2 contact (%1.f%%:%d).\n",
					c2_res,ResidueID(Contact[i]),
					c_res,ResidueID(ResCD[j]),
					contact_surface,Ratio[i]*100.0,hp_pos);
			} hp_pos++;
			// PrintResidueSurface(fptr,ResCD[j]);
		   } hpsz=0;
		}
		Nildheap(dH);
		// fprintf(fptr,"# total percent surface = %.1f\n",total_ratio*100.0);
	       }
	     }
	 return TotalContacts;
} 


