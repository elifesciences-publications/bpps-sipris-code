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

#include "msites.h"

mst_type	MakeMSites(Int4 ntyps, Int4 *len_elem, gss_typ &gss0)
// makes a copy of gapped sequence seq constructor 
{
    mst_type	S;

    if(ntyps > MAX_NO_TYPESITES) msites_error("Too many types.");
    S = new msites_type[1]; 
    S->gss=gss0;	// copy constructor called.
    mk_msites(ntyps, len_elem, S);
    return S;
}

mst_type	MkMSites(Int4 ntyps, Int4 *len_elem, ss_type data)
// creates a gapped sequence set from ss_type sequence set.
{
    mst_type	S;

    if(ntyps > MAX_NO_TYPESITES) msites_error("Too many types.");
    S = new msites_type[1]; 
    // S->gss is initialized by constructor function.
    S->gss.initialize(data,10,5,5.0,9,9);	// initialize to full dataset.
    mk_msites(ntyps, len_elem, S);
    return S;
}

void	mk_msites(Int4 ntyps, Int4 *len_elem, mst_type S)
/*********************************************************************
 Initialize a site structure containing no sites. 
							 nsites: A B C
  seq[1]    = .................................................* 0 0 0
  seq[2]    = ...........................................*	 0 0 0
     :
  seq[nseq] = ...............................................*   0 0 0
**********************************************************************/
{
    Int4	n,m,t,max;
    ss_type	data=MSitesSeqSet(S);

    S->ntyp = ntyps;
    S->nseq = NSeqsSeqSet(data);
    S->len_seq = LengthsSeqSet(data);
    NEWP(S->type, S->nseq+1, char);
    NEW(S->pos, S->nseq+1, ol_type);
    for(max=0,n = 1; n <= S->nseq; n++){
	max = MAXIMUM(Int4,max,S->len_seq[n]);
	NEW(S->type[n],(S->len_seq[n]+2),char);
	S->type[n][(S->len_seq[n]+1)] = ENDTYPESITE;
	S->pos[n] = Olist(S->len_seq[n]+1);
    }
    NEW(S->tmp,max+1,Int4);
    NEW(S->len_elem, ntyps+1, Int4);
    NEW(S->totsites, ntyps+1, Int4);
    NEWP(S->nsites,ntyps+1,Int4);
    NEWPP(S->site_pos, ntyps+1,unsigned short);
    for(t = 1; t <= ntyps; t++) {
	S->totsites[t] = 0;
        S->len_elem[t] = len_elem[t];
        NEW(S->nsites[t], S->nseq+1,Int4);
	NEWP(S->site_pos[t], NSeqsSeqSet(data)+1,unsigned short);
	for(n = 1; n <= NSeqsSeqSet(data); n++) {
		S->nsites[t][n] = 0;
                m = ((Int4) SqLenSeqSet(n,data)/(Int4)len_elem[t]) +3;
                NEW(S->site_pos[t][n],m+2,unsigned short);
	}
	if(t==1) S->maxinc = S->len_elem[t] - 1;
	else S->maxinc = MINIMUM(Int4,S->maxinc,(S->len_elem[t]-1));
    }
    S->maxinc = MAXIMUM(Int4,1,S->maxinc);
}

void	NilMSites(mst_type S)
/* destroy sites structure S. */
{
    Int4	t,n;

   free(S->len_elem); free(S->len_seq); free(S->tmp); free(S->totsites); 
   for(n = 1; n <= S->nseq; n++) { NilOlist(S->pos[n]); free(S->type[n]); }
   free(S->type); free(S->pos);
   for(t = 1; t <= S->ntyp; t++) {
   	for(n = 1; n <= S->nseq; n++) free(S->site_pos[t][n]);
	free(S->site_pos[t]); free(S->nsites[t]);
   }
   free(S->site_pos); free(S->nsites);
   delete []S; // replaces free(S);
}

void	VacateMSites(Int4 n, mst_type S)
// remove all sites from nth sequence.
{
   Int4		i,t,s,k;

   for(t=1; t<=nTypeMSites(S); t++){
	k = nMSites(t,n,S);
	for(i=1;i<=k;i++){ s=S->site_pos[t][n][i]; VacateMSite(t,n,s,S);}
   }
}

void	VacateMSite(Int4 t, Int4 n, Int4 site, mst_type S)
/*********************** vacate site *******************
 Remove site in sequence n of type t by vacating all positions. 
     site 
      |
      A   a   a   a   ...   a   a      	<- sequence n
    |+0 |+1 |+2 |+3 | ... |w-2|w-1|	
      o   o   o   o   ...   o   o           	
        for p = site ... site + w - 1 -> type[n][p] = VACANT 
*******************************************************/
{
	Int4	i,end;
	unsigned short  *site_pos=S->site_pos[t][n];

	if(S->type[n][site] != t){
		PutTypeMSites(stderr, S);
		PutMSites(stderr, t, S, NULL,NULL);
		fprintf(stderr,"ELEMENT %c; seq %d; site %d\n",
			'A' +t -1, n,site);
		msites_error("attempt to remove site where none exists.");
	}
	RmOlist(site,S->pos[n]);
	end = site + S->len_elem[t] -1; 
	for(i=site; i <= end ; i++) S->type[n][i] = VACANT; 
/****************************************************************
            i= 1   2   3   4    (nsites = 4)       
  remove(35, [12, 35, 78, 104, NULL]) 

	-> [12, 104, 78, 104, NULL] -> [12, 104, 78, NULL ] 

  remove(35, [35, NULL]) 

	-> [35, NULL] -> [NULL]

 ****************************************************************/
	for(i=1; i <= S->nsites[t][n]; i++){
		if(site_pos[i] == site){
			site_pos[i] = site_pos[S->nsites[t][n]];
			site_pos[S->nsites[t][n]] = NULL;
			S->nsites[t][n]--;
			S->totsites[t]--;
			return ;
		}
	}
	msites_error("VacateMSite( ) - this should not happen");
}

void	AddMSite(Int4 t, Int4 n, Int4 site, mst_type S)
/*************************** add site ********************************
 Add type t site in sequence n.
     site 
      |
      o   o   o   o   ...   o   o           	
    |+0 |+1 |+2 |+3 | ... |w-2|w-1|	
      A   a   a   a   ...   a   a      	<- sequence n
	type[n][site] = t  ('A')
        for p = site + 1 ... site + w - 1 -> type[n][p] = BLOCKED (t) 
***********************************************************************/
{
	Int4	p,end;

	if(OccupiedMSite(t, n, site, S)){
		PutTypeMSites(stderr, S);
		PutMSites(stderr, t, S, NULL,NULL);
		fprintf(stderr,
			"ELEMENT %c(length=%d); seq %d; site %d\n",
			'A' +t -1, S->len_elem[t], n,site);
		S->gss.Put(stderr,n);
		S->gss.PutFA(stderr);
		fprintf(stderr,
			"sequence length = %d; numseqs = %d\n", 
			SeqLenMSites(n,S),S->nseq);
		msites_error("attempt to add site where one exists.");
	}
	InsertOlist(site, S->pos[n]);
	S->type[n][site] = t;
	end = site + S->len_elem[t] - 1; 
	for(p=site+1; p <= end ; p++) S->type[n][p] = BLOCKED(t); 
	S->nsites[t][n]++;
	S->totsites[t]++;
	S->site_pos[t][n][S->nsites[t][n]] = (unsigned short) site;
}

BooLean OccupiedMSite(register Int4 t, register Int4 n, register Int4 site, 
	register mst_type S)
/* determine if a site is blocked. */
{
	register Int4	p,end;

	end = site + S->len_elem[t] - 1; 
	for(p=site; p<end ; p+=S->maxinc){if(S->type[n][p])return TRUE;}
	if(S->type[n][end]) return TRUE;
	return FALSE;
}

Int4	PosTMSites(Int4 t, Int4 n, Int4 *pos, mst_type S)
/* modifies array pos to contain the positions of type t sites in seq. n */
{
	Int4	s,site,i;

	GetOlist(S->tmp, S->pos[n]); 
	for(i=s=1; (site=S->tmp[s]) != 0; s++) {
		if(S->type[n][site]==t) pos[i++] = S->tmp[s];
	}
	pos[i] = 0;
	return (i-1);
}

Int4	PosMSites(Int4 n, Int4 *pos, mst_type S)
/* modifies array pos to contain the positions of sites in sequence n */
{
	Int4	s;
	GetOlist(S->tmp, S->pos[n]); 
	for(s=1; S->tmp[s] != 0; s++) pos[s] = S->tmp[s];
	pos[s] = 0;
	return (s-1);
}

Int4	GapBetweenMSites(Int4 n, Int4 t, mst_type S)
// return the gap length between sites
{
	if(t < 1 || t >= S->ntyp) return 0;
	else return (MSitePos(t+1,n,1,S) - EndMSitePos(t,n,1,S) - 1);
}

/********************************** output ******************************/
void	PutMSites(FILE *fptr,Int4 t,mst_type S,double **site_prob, BooLean *off)
{
	ss_type	P=MSitesSeqSet(S);
	Int4	i,n,s,e,end,N,start,length;
	BooLean	are_sites = FALSE;
	e_type	E;
	char	r,c;
        Int4	flank=10; 
   
	fprintf(fptr,"\n\n");
	length = MSiteLen(t,S);
	for(N=0, n=1; n<= NSeqsSeqSet(P); n++){
	   bubble_sort_msites(S->site_pos[t][n]);
	   E=SeqSetE(n,P);
	   if(nMSites(t,n,S) > 0) N++; 
	   for(s=1; s<= nMSites(t,n,S); s++){
		are_sites = TRUE;
	   	fprintf(fptr,"%2d-%-2d ",n,s);
		start= S->site_pos[t][n][s]; 
		fprintf(fptr,"%4d  ",start);
		e = start + length - 1;
		end = e + flank;
		for(i=start-flank; i <= end; i++){
		   if(i < 1 || i > (Int4) LenSeq(E)) fprintf(fptr," ");
		   else{
			r = ResSeq(i,E);
			if(OpenPos(n,i,S)) c = AlphaCharLow(r,SeqSetA(P));
			else c = AlphaChar(r,SeqSetA(P));
			if(i == e) fprintf(fptr,"%c ", c);
			else if(i == start) fprintf(fptr," %c", c);
			else fprintf(fptr,"%c", c);
		   }
		}
		fprintf(fptr," %4d",e);
		if(site_prob != NULL) {
			fprintf(fptr," (%0.2f)",site_prob[n][start]);
			/** PutSeqID(fptr,E); /*** TEMP ***/
		} 
		fprintf(fptr,"\n");
	   }
	}
	if(are_sites){
	   if(off != NULL){
	   	fprintf(fptr,"sites:%17c",' ');
	   	for(s=1; s <= S->len_elem[t]; s++){
			if(off[s])fprintf(fptr," "); 
			else fprintf(fptr,"*"); 
	  	 }
	   }
	   fprintf(fptr,"\n%*s", 23, "");
	   for(i=1; i<= MSiteLen(t,S); i++) {
		if(i%5==0) fprintf(fptr,"%5d",i);
	   }
	}
	fprintf(fptr,"\n\t(%d sites in %d sequences)\n\n",S->totsites[t],N);
}

void    PutTypeMSites(FILE *fptr, mst_type S)
/*  	e.g.    3: D(24)-A(35)-C(45)-[59]  */
{
	Int4 i,n,t,tot;

	fprintf(fptr,"\n");
	for(n=1; n<=S->nseq; n++){
	   for(tot=0,t=1; t<=S->ntyp; t++) tot += S->nsites[t][n];
	   if(tot > 0){
	      fprintf(fptr,"%3d: ",n);
	      for(i=1; i<=S->len_seq[n]; i++){
		if(S->type[n][i] > 0){
		   t = S->type[n][i];
	   	   if(tot == 0) fprintf(fptr,"-"); tot = 0;
		   if(S->ntyp <= 26){
			fprintf(fptr,"%c(%d)",('A'+t-1),i);
		   } else { fprintf(fptr,"%d(%d)",t,i); }
		}
	      }
	      fprintf(fptr,"-[%d]\n",S->len_seq[n]);
	   }
	}
	fprintf(fptr,"\n");
}

void	PutScanMSites(FILE *fptr, Int4 t, mst_type S, BooLean *off)
{ PutScanMSitesProb(fptr, t, S, NULL, off, 0.0); }

void	PutScanMSitesProb(FILE *fptr, Int4 t, mst_type S, double **prob, BooLean *off,
	double cutoff)
/** Create a Scan file (*.sn) eliminating poor matches from the aligned block **/
{
        ss_type P=MSitesSeqSet(S);
        Int4    i,n,s,e,end,N,start,length;
        e_type  E;
        char    r,c;

        if(off != NULL){
                for(s=1; s <= S->len_elem[t]; s++){
                        if(off[s]==TRUE)fprintf(fptr,".");
                        else if(off[s]==FALSE)fprintf(fptr,"*");
                        else fprintf(fptr,"^");
                 }
        } else for(s=1; s <= S->len_elem[t]; s++) fprintf(fptr,"*");
        fprintf(fptr,"\n");
        length = MSiteLen(t,S);
        for(N=0, n=1; n<= NSeqsSeqSet(P); n++){
           bubble_sort_msites(S->site_pos[t][n]);
           E=SeqSetE(n,P);
           if(nMSites(t,n,S) > 0) N++;
           for(s=1; s<= nMSites(t,n,S); s++){
             start= S->site_pos[t][n][s];
             if(prob == NULL || prob[n][start] >= cutoff) {
                e = start + length - 1;
                end = e;
                for(i=start; i <= end; i++){
                   if(i < 1 || i > (Int4) LenSeq(E)) fprintf(fptr," ");
                   else{
                        r = ResSeq(i,E);
                        c = AlphaChar(r,SeqSetA(P));
                        fprintf(fptr,"%c", c);
                   }
                }
                fprintf(fptr,"\n");
             }
           }
        }
        fprintf(fptr,"\n");
}

/********************************** private ******************************/
void	print_msites(unsigned short *L)
{
	Int4	i;
	L++;
	fprintf(stderr,"\n");
	for(i=0; L[i]!=NULL; i++) fprintf(stderr," %d",L[i]);
	fprintf(stderr," (null)\n");
}

Int4     bubble_sort_msites(unsigned short *L)
{
        unsigned short i,j,t,n;
	Int4 temp=0;

	L++;
        for(n=i=0; L[i]!=NULL; i++) n++;
        for(i=1; i < n; i++){
                for(j=i;j > 0 && (L[j-1] > L[j]);j--){
                        t=L[j]; L[j]=L[j-1]; L[j-1]=t;
                }
        }
	return temp;
}

void	msites_error(char *s){fprintf(stderr,"sites error: %s\n",s);exit(1);}

