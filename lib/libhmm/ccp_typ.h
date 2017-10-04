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

  /* Program COILS version 2.1 */
  /* written by J.Lupas, 22 JUN 1993 */
  /* edited on 15 JAN 1994 */
  /* revised and corrected 4 APR 1994 */
  /* incorporates the option to use either the old MTK chart or the new MTIDK
     chart AND the choice of weighting positions a & d more (2.5 times) */
  /* 4 output options:
        - probabilities for window sizes 14, 21 & 28 in columns,
        - probabilities for window sizes 14, 21 & 28 in rows,
        - scores only for a user-input window size, or
        - the probabilities above a user-input cutoff.
     Requests input and output files from user. */
  /* transferred to c++ by Larry Harvie, MX%"harvie@owl.WPI.EDU" */ 

/***********************
ccp_typ.[ch]: 08/16/99 
 Modified by Trevor Carlson <tcarlson@andrew.cmu.edu> with assistance from
             Dr. Andrew Neuwald <neuwald@cshl.org>

ccp_typ(): 
     Default ctor is not operational.

ccp_typ( a_type, Int4 ):
     ctor:  supply the a_type used in creating the e_type, and the maximum
                length of any sequence to be submitted;

float *getArray( unsigned char *, Int4, char*, char):
     Accessor: unsigned char *, to be used with an e_type:     SeqPtr(E)
               Int4           , length of this sequence
	       char *         , discription of the sequence(NULL Terminated)
	       char           , option character

     The options are as follows:
               'b' or 'B' - Binary: Heptad window is 7*2 in length (14)
	       't' or 'T' - Tertiary:      window is 7*3 in length (21)
	       'q' or 'Q' - Quartenary:    window is 7*4 in length (28)
     **NOTE: The length of the sequence must be greater than or equal to the window size.

	       Capital letters indicate a higher weight for amino acids of type a and d.
	       Lower case letters will give an equal weight for all amino acids.

     Returns a float array of length long+1, with the first probability at location 1.
     **NOTE: It is the callers responsibility to destroy the float array.

 ************************/

#ifndef __ccp_typ_h__
#define __ccp_typ_h__

#include "stdinc.h"
#include "alphabet.h"
#include "residues.h"
#include "sequence.h"

class ccp_typ {
public:

  ccp_typ    ();
  ccp_typ    (a_type,Int4);
  ~ccp_typ   ();

  float         *getArray( unsigned char *, Int4, char );
  BooLean	MaskCoils(e_type E,Int4 minseg,double cutoff);
  BooLean	Put(FILE *fp, e_type E, double cutoff);
private:
  char          weightchar;
  unsigned char *residue;
  Int4          window;
  Int4          init;
  Int4          res_total;        /* total no. of residues in protein  */
  Int4          chartoption;
  Int4          hept_weight;
  Int4		MaxLength;
  Int4          *heptnum;
  Int4          *temphept;
  // float         *f;                /* Pointer to array of probabilities */
  double         ad_weight;
  double        *calcnumb;
  double        *tempcalc;
  double       	**TheChart;
  a_type        A;

  void          init_memory(Int4);
  void          Calculate(Int4);
  void          Array_probs(float*);
  double        Calcprob(double,double,double,double,double,double);
  void		printOutput(FILE*,char*,Int4,Int4,Int4,float*,float*,float*);
 };

#endif

