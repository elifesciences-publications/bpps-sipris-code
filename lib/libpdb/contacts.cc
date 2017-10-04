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

#include "contacts.h"
#include "energies.h"

ct_type  MkContacts(char *infile, Int4 dmax, a_type A)
{ return MkContactsF(infile, dmax, A, NULL); }

ct_type  MkContactsF(char *infile, Int4 dmax, a_type A, char *ss)
/** filter out contacts in secondary structure (ss) coil regions **/
{
	ct_type	X;
        Int4	k,length,N,i,r,r1,r2,d,v,s1,s2;
	char	str[200];
        FILE    *fptr;

	NEW(X,1,contacts_type);
	dmax = dmax - 4;  	/** 1 = 0-5, 2 = 5-6, ..., 6 = 9-10 **/
	if(dmax < 1) contacts_error("dmax input error in MkContacts( )");
	X->dmax = dmax;
	X->A = A;
	fptr = open_file(infile,".rrc","r");
	for(N=0,X->Nrr=0; fgets(str,150,fptr) != NULL; ) {
	    if(sscanf(str,"%d %d %d",&r1,&r2,&d)==3) {
	      if(d <= dmax){
		if(ss == NULL || (ss[r1]!='c' && ss[r2]!='c')){
		    X->Nrr++;
		    if(N < r1) N = r1;
		    if(N < r2) N = r2;
		}
	      }
	    } else contacts_error(".rrc input file error");
	}
	fclose(fptr);

	fptr = open_file(infile,".rpc","r");
	for(X->Nrp=0; fgets(str,150,fptr) != NULL; ) {
	    if(sscanf(str,"%d %d %d",&r1,&r2,&d)==3) {
	      if(d <= dmax){
		if(ss == NULL || (ss[r1]!='c' && ss[r2]!='c')){
		    X->Nrp++;
		    if(N < r1) N = r1;
		    if(N < r2) N = r2;
		}
	      }
	    } else contacts_error(".rpc input file error");
	}
	fclose(fptr);

	length = X->len = N;
	NEW(X->contact,N+2,BooLean);
	fptr = open_file(infile,".rrc","r");
	NEW(X->rrd,X->Nrr+2,Int4);
	NEW(X->rr1,X->Nrr+2,Int4);
	NEW(X->rr2,X->Nrr+2,Int4);
	for(N=1; fgets(str,150,fptr) != NULL; ){
	    if(sscanf(str,"%d %d %d",&r1,&r2,&d)==3){
	      if(d <= dmax){
		if(ss == NULL || (ss[r1]!='c' && ss[r2]!='c')){
		   if(N > X->Nrr) print_error(".rrc input file error!");
		   X->rr1[N] = r1; X->rr2[N] = r2;
	    	   X->rrd[N] = d; N++;
		   X->contact[r1] = TRUE; X->contact[r2] = TRUE;
		}
	      }
	    }
	}
	fclose(fptr);
	
	NEW(X->rpd,X->Nrp+2,Int4);
	NEW(X->rp1,X->Nrp+2,Int4);
	NEW(X->rp2,X->Nrp+2,Int4);
	fptr = open_file(infile,".rpc","r");
	for(N=1; fgets(str,150,fptr) != NULL; ){
	    if(sscanf(str,"%d %d %d",&r1,&r2,&d)==3) {
	      if(d <= dmax){
		if(ss == NULL || (ss[r1]!='c' && ss[r2]!='c')){
		   if(N > X->Nrp) print_error(".rpc input file error!");
		   X->rp1[N] = r1; X->rp2[N] = r2;
	    	   X->rpd[N] = d; N++;
		   X->contact[r1] = TRUE; 
		   /**X->contact[r2]=TRUE; /** Don't contact backbone **/
		}
	      }
	    }
	}
	fclose(fptr);

	for(X->vnm=0,i=1;i<=(Int4)X->len;i++) if(X->contact[i]) X->vnm++;
	NEWPP(X->enm,23,double);
	for(r2=1; r2<=21; r2++){
	   NEWP(X->enm[r2],23,double);
	   NEW(X->enm[r2][0],23,double);
	   for(r1=1; r1<=21; r1++){
	      NEW(X->enm[r2][r1],7,double);
	      for(d=1; d<=6; d++){
		 X->enm[r2][r1][d] = contact_energies[d][r1][r2];
	      }
	   }
	}
        NEWP(X->mcc,length+2,Int4); NEWP(X->mcd,length+2,Int4);
        NEWP(X->mcp,length+2,double); 
	NEW(X->mcn,length+2,Int4); NEW(X->nrp,length+2,Int4);
        for(i=1;i<=X->Nrr;i++){
           s1=X->rr1[i]; X->mcn[s1]++;
	   s2=X->rr2[i]; X->mcn[s2]++;
        }
        for(i=1;i<=X->len;i++){
         if(k=X->mcn[i]){     /** i.e., there are some contacts at i. **/
		X->mcn[i]=0; 
		NEW(X->mcc[i],k+2,Int4);
		NEW(X->mcd[i],k+2,Int4);
         }
         NEW(X->mcp[i],25,double);
        }
        for(i=1;i<=X->Nrr;i++){
           s1=X->rr1[i]; s2=X->rr2[i];
	   X->mcn[s1]++;
           v=X->mcn[s1]; X->mcc[s1][v]=s2; X->mcd[s1][v]=X->rrd[i];
           X->mcn[s2]++;
           v=X->mcn[s2]; X->mcc[s2][v]=s1; X->mcd[s2][v]=X->rrd[i];
        }
        for(i=1;i<=X->Nrp;i++){
           s1=X->rp1[i]; d=X->rpd[i];
	   X->nrp[s1]++;
           for(r=1; r <= nAlpha(X->A); r++){
                X->mcp[s1][r]+=X->enm[backbone_peptide][r][d];
	   }
        }
	return X;
}

void    PutSeqRRContacts(FILE *fptr, e_type  E, ct_type W)
{
        Int4	i,r1,r2,d;
	char	c1,c2,mtrx[21][21];

	for(r1 = 0; r1 <= nAlpha(W->A); r1++){
	   for(r2 = 0; r2 <= nAlpha(W->A); r2++) mtrx[r1][r2] = 0;
	}
        for(i = 1; i <= W->Nrr; i++){
                r1 = ResSeq(W->rr1[i],E);
                r2 = ResSeq(W->rr2[i],E);
                d = W->rrd[i];
		if(r1 <= r2) mtrx[r1][r2]++;
		else mtrx[r2][r1]++;
        }
	for(r1 = 1; r1 <= nAlpha(W->A); r1++){
	   for(r2 = r1; r2 <= nAlpha(W->A); r2++){
        	fprintf(fptr,"%c-%c %-4d\n", AlphaChar(r1,W->A), 
			AlphaChar(r2,W->A), mtrx[r1][r2]);
	   }
	}
}

void    PutSeqContactsShort(FILE *fptr, e_type  E, ct_type W)
{
        Int4     i,r1,r2,d,s,s1,s2;
        double  e;

        for(i = 1; i <= W->Nrr; i++){
                r1 = ResSeq(W->rr1[i],E);
                r2 = ResSeq(W->rr2[i],E);
                d = W->rrd[i];
                fprintf(fptr,"%g\n", contact_energies[d][r1][r2]);
        }
}

void    PutSeqContacts(FILE *fptr, e_type  E, ct_type W)
{
        Int4     i,r1,r2,d,s,s1,s2;
        double  e;
        const char     *str[8] = {"error"," 0-5 ",">5-6 ",">6-7 ",">7-8 ",
                        ">8-9 ",">9-10","error"};

        fprintf(fptr,"Sequence:\n");
        PutSeq(fptr,E,W->A);
        fprintf(fptr,"Residue/Residue contacts:\n");
        for(i = 1; i <= W->Nrr; i++){
                r1 = ResSeq(W->rr1[i],E);
                r2 = ResSeq(W->rr2[i],E);
                d = W->rrd[i];
                fprintf(fptr,"%4d: %4d %c-%c %-4d (%s Angstroms) E = %g\n",
                        i,W->rr1[i],AlphaChar(r1,W->A),
                        AlphaChar(r2,W->A),W->rr2[i],str[d],
                        contact_energies[d][r1][r2]);
        }
        fprintf(fptr,"\nResidue/Backbone contacts:\n");
        for(i = 1; i <= W->Nrp; i++){
                r1 = ResSeq(W->rp1[i],E);
                r2 = ResSeq(W->rp2[i],E);
                d = W->rpd[i];
                fprintf(fptr,"%4d: %4d %c-b(%c) %-4d (%s Angstroms) E = %g\n",
                        i,W->rp1[i],AlphaChar(r1,W->A),
                        AlphaChar(r2,W->A),W->rp2[i],str[d],
                        contact_energies[d][r1][backbone_peptide]);
        }
        fprintf(fptr,"\nContacts energies for each site\n");
        for(s=1; s <= W->len; s++){
          if(W->contact[s]){
            r1 = ResSeq(s,E);
            fprintf(fptr,"%3d(%c): (peptide contacts = %g)\n",
                s,AlphaChar(r1,W->A),W->mcp[s][r1]);
            for(i=1;i<=W->mcn[s];i++){
                d = W->mcd[s][i];
                s2=W->mcc[s][i];
                r2 = ResSeq(s2,E);
                e = W->enm[r2][r1][d]; 
                fprintf(fptr,"\t%2d: %1d %3d(%c) e = %g\n",
                        i,d,s2,AlphaChar(r2,W->A),e);
            }
            fprintf(fptr,"\n");
          }
        }
        fprintf(fptr,"\n\n");
}

Int4	**ResResContacts(ct_type W)
{
	Int4	s,i,**contacts;

	NEWP(contacts,W->Nrr + 3, Int4);
	for(i = 1; i <= W->Nrr; i++){
		NEW(contacts[i],3,Int4);
		contacts[i][0]=W->rrd[i]+4;
		contacts[i][1] = W->rr1[i];
		contacts[i][2] = W->rr2[i];
	}
	return contacts;
}

void    PutContacts(FILE *fptr, ct_type W)
{
	Int4	i,r1,r2,d,s,s1,s2;
        double  e;
	const char	*str[8] = {"error"," 0-5 ",">5-6 ",">6-7 ",">7-8 ",
			">8-9 ",">9-10","error"};

	fprintf(fptr,"Residue/Residue contacts:\n");
	for(i = 1; i <= W->Nrr; i++){
		d = W->rrd[i];
		fprintf(fptr,"%4d: %4d-%-4d (%s Angstroms)\n",
			i,W->rr1[i], W->rr2[i],str[d]);
	}
	fprintf(fptr,"\nResidue/Backbone contacts:\n");
	for(i = 1; i <= W->Nrp; i++){
		d = W->rpd[i];
		fprintf(fptr,"%4d: %4d-b(%-4d) (%s Angstroms)\n",
			i,W->rp1[i], W->rp2[i],str[d]);
	}
}

double	NegEnergyDiffContacts(Int4 s, Int4 New, Int4 old, e_type E, void *X)
{ return -EnergyDiffContacts(s, New, old, E, X); }

double	NegTotalEnergyContacts(register e_type E, register ct_type C)
{ return - TotalEnergyContacts(E, C); }

double	EnergyDiffContacts(Int4 s, Int4 New, Int4 old, e_type E, void *X)
/** find the energy difference when replacing old by New at site s **/
{
	register unsigned char *seq = SeqPtr(E);
        register Int4   i,*s2,*d;
	ct_type		C = (ct_type) X;
        register double **enew=C->enm[New],**eold=C->enm[old];
	register double	score;

	if(s > C->len) return 0.0;
        s2=C->mcc[s]; d=C->mcd[s];
	score = C->mcp[s][New] - C->mcp[s][old];
        for(i=C->mcn[s]; i > 0;i--){	/** new - old residue  **/
                score += (enew[seq[s2[i]]][d[i]]
                		- eold[seq[s2[i]]][d[i]]); 
        }
        return score;
}

double	PutResEnergyContacts(FILE *fptr, e_type E, ct_type W)
/** WARNING: THIS DOESN'T LOOK RIGHT FOR BACKBONE CONTACTS!!! **/
{
        Int4 i,d,r1,r2;
        double *e,energy,*re;

	NEW(e,W->len+3,double);
	NEW(re,nAlpha(W->A)+3,double);
        for(energy=0,i=1; i<=W->Nrr; i++){
                r1 = ResSeq(W->rr1[i],E);
                r2 = ResSeq(W->rr2[i],E);
                d = W->rrd[i];
                energy += contact_energies[d][r1][r2];
		e[W->rr1[i]] += contact_energies[d][r1][r2];
		e[W->rr2[i]] += contact_energies[d][r1][r2];
		re[r1] += contact_energies[d][r1][r2];
		re[r2] += contact_energies[d][r1][r2];
        }
        for(i = 1; i <= W->Nrp; i++){
                r1 = ResSeq(W->rp1[i],E);
                r2 = ResSeq(W->rp2[i],E);
                d = W->rpd[i];
                energy += contact_energies[d][r1][backbone_peptide];
		e[W->rp1[i]] += contact_energies[d][r1][r2];
		e[W->rp2[i]] += contact_energies[d][r1][r2];
		re[r1] += contact_energies[d][r1][r2];
        }
	PutSeqID(fptr,E);
	fprintf(fptr,"\nres Energy\n");
	for(r1=1; r1<= nAlpha(W->A); r1++){
		fprintf(fptr,"%c: %g\n",AlphaChar(r1,W->A),re[r1]);
		fprintf(fptr,"\tpos: Energy\n");
		for(i=1; i<= W->len; i++){
                	r2 = ResSeq(i,E);
			if(r2 == r1){
				fprintf(fptr,"\t%3d: %g\n", i,e[i]);
			}
		}
	}
	fprintf(fptr,"total = %g\n\n",energy);
	free(e); free(re);
        return energy;
}

Int4    NumContacts(Int4 s, ct_type C) { return C->mcn[s]; }

double	ResEnergyContacts(Int4 s, register unsigned char *seq, ct_type W)
/** [needs comment] **/
{
        register Int4   i,*s2,*d, r;
        register double ***e=W->enm,E;

	if(s > W->len) return 0.0;
	r=seq[s];
        s2=W->mcc[s]; d=W->mcd[s];
	E = W->mcp[s][r];
        for(i=W->mcn[s]; i > 0;i--){	/** new - old residue  **/
	    E  += e[r][seq[s2[i]]][d[i]];
        }
        return E;
}

double	TotalEnergyContacts(register e_type E, register ct_type W)
{
        register Int4 i,d,r1,r2;
        register double energy;

        for(energy=0,i=1; i<=W->Nrr; i++){
                r1 = ResSeq(W->rr1[i],E);
                r2 = ResSeq(W->rr2[i],E);
                d = W->rrd[i];
                energy += contact_energies[d][r1][r2];
        }
        for(i = 1; i <= W->Nrp; i++){
                r1 = ResSeq(W->rp1[i],E);
                r2 = ResSeq(W->rp2[i],E);
                d = W->rpd[i];
                energy += contact_energies[d][r1][backbone_peptide];
        }
        return energy;
}

void    NilContacts(ct_type W)
{
	Int4 i,r1,r2,d;

        for(i=1;i<=W->len;i++){
	 if(W->contact[i]){ /** i.e., there are some contacts at i. **/
		free(W->mcc[i]); free(W->mcd[i]);
         }
         free(W->mcp[i]);
        }
        free(W->mcc); free(W->mcd); free(W->mcp); free(W->mcn);
	free(W->nrp);
	for(r2=1; r2<=21; r2++){
	   free(W->enm[r2][0]);
	   for(r1=1; r1<=21; r1++) free(W->enm[r2][r1]);
	   free(W->enm[r2]);
	}
	free(W->enm);
	free(W->rpd); free(W->rp1); free(W->rp2);
	free(W->rrd); free(W->rr1); free(W->rr2);
	free(W->contact);
	free(W);
}

/******************************** private *********************************/
void	contacts_error(char *s){fprintf(stderr,"Contacts: %s\n",s);exit(1);}

