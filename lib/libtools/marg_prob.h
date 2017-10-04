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

#include "cmsa.h"
#include "my_ncbi.h"

sap_typ TrimSAP(sap_typ sap, Int4 trimleft, Int4 trimright);
sap_typ	MargProbTrimSAP(sap_typ sap, Int4 **matrix, double PerNats, double misalncut,
	Int4 gapopen, Int4 gapextend, Int4 window, a_type AB);
BooLean MarginalProbTrim(sap_typ sap, Int4 **matrix, double PerNats, double misalncut,
	Int4 gapopen, Int4 gapextend, Int4 window, Int4 trim[2], a_type AB);
double  **LocalMarginalProb(Int4 seq_len, unsigned char *seq, Int4 length, Int4 **matrix, double PerNats,
	Int4 gapopen, Int4 gapextend, a_type AB) ;
double **MargProbCore(Int4 seq_len, unsigned char *seq, Int4 length, double **pmat_emit_odd, double PerNats,
                Int4 gapopen, Int4 gapextend, a_type AB);
double  **LocalMarginalProb2(Int4 seq_len, unsigned char *seq, char *operation, Int4 length, 
	double **WtFreq, double *fract_seq_aln, double *bfreq, 
		double PerNats, Int4 gapopen, Int4 gapextend, a_type AB);
double  **MarginalProb2(e_type sE, char *operation, Int4 Start, Int4 Oper_len, 
        double *score, Int4 length, Int4 **mat_emit_odd, Int4 PerNats, Int4 gapopen, 
	        Int4 gapextend, a_type AB);
double *ScoreMargOperArr(double **TMAT, char *operation, Int4 Start, Int4 Oper_len);
