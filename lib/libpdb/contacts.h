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

/****************** contacts.h ***************/
#if !defined(CONTACTS)
#define CONTACTS
#include "sequence.h"
#include "afnio.h"
#include <math.h>
#include "alphabet.h"
#include "residues.h"
#include "histogram.h"
#include "probability.h"
#include "random.h"

/*** #define ThreadingAminoAcids "XARNDCQEGHILKMFPSTWYV" /** use residues.h **/
#define backbone_peptide 21

/************************** contacts datatype *****************************/
typedef struct {
	Int4		len;	/** number of residues in structure **/
	a_type		A;	/** residue alphabet **/
	BooLean		*contact;  /** contact[s] = site s has a contact **/
	Int4		dmax;	/** maximum contact energy distance **/
	Int4		vnm;	/** number of residues contacted **/
	Int4		*rrd;	/** rrd[Nrr] = 1-6 **/
	Int4		*rr1;	/** rrd[Nrr] = residue 1 **/
	Int4		*rr2;	/** rrd[Nrr] = residue 2 **/
	Int4		Nrr;
	Int4		*rpd;	/** rpd[Nrp] = 1-6 **/
	Int4		*rp1;	/** rpd[Nrp] = residue 1 **/
	Int4		*rp2;	/** rpd[Nrp] = residue 2 **/
	Int4		Nrp;
	double		***enm;	/** energy matrix **/
	double		**mcp;	/* mcp[s][r] sum of peptide contact energy */
	Int4		**mcd;  /* mcd[s][1..n] contact distances for site s */
	Int4		**mcc;  /* mcc[s][1..n] other sites contacting s */
	Int4		*mcn;   /* mcn[s]  total # of rrc contacts with s */
	Int4		*nrp;   /* nrp[s]  total # of rpc contacts with s */
} contacts_type;

typedef	contacts_type	*ct_type;

/******************************** private **********************************/
void    contacts_error(char *s);

/******************************** PUBLIC ***********************************/
ct_type  MkContacts(char *infile, Int4 dmax, a_type A);
ct_type  MkContactsF(char *infile, Int4 dmax, a_type A, char *ss);
void    PutContacts(FILE *fptr, ct_type W);
void    PutSeqContacts(FILE *fptr, e_type  E, ct_type W);
void    PutSeqContactsShort(FILE *fptr, e_type  E, ct_type W);
double  PutResEnergyContacts(FILE *fptr, e_type E, ct_type W);
void    PutSeqRRContacts(FILE *fptr, e_type  E, ct_type W);
Int4	NumContacts(Int4 s, ct_type W);
double  ResEnergyContacts(Int4 s, register unsigned char *seq, ct_type W);
Int4    **ResResContacts(ct_type W);
void	NilContacts(ct_type W);
double  EnergyDiffContacts(Int4 s, Int4 New, Int4 old, e_type E, void *X);
double  TotalEnergyContacts(register e_type E, register ct_type W);
double  NegEnergyDiffContacts(Int4 s, Int4 New, Int4 old, e_type E, void *X);
double  NegTotalEnergyContacts(register e_type E, register ct_type C);
/******************************** PUBLIC ***********************************/
#define ContactsA(W)		((W)->A)
#define NumRRContacts(W)	((W)->Nrr)
#define nRRContactsSite(s,W)	((W)->mcn[(s)])
#define nRPContactsSite(s,W)	((W)->nrp[(s)])
#define HasContactFast(s,W)	((W)->contact[(s)])
#define HasContact(s,W)		(((W)->len >= (s)) && (W)->contact[(s)])
#define LenContacts(W)		((W)->len)
#define dmaxContacts(W)		((W)->dmax+4)

#endif

