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

/* chn_vsi.h - CHAIN analysis 'View Structural Interaction' program. */
#if !defined (VSI_PDB)
#define VSI_PDB

#include "pdb.h"
#include "vsi_typ.h"
#include "ras_typ.h"
#include "histogram.h"
#include "set_typ.h"

vsi_typ PutAutoFormatVSI(FILE *fp,Int4 file_id, vsi_typ HEAD_NODE,
                char KEY_CHAIN,float HA_dmax,float dmax, BooLean only_if_trace,
                pdb_typ P,char mode=0);
Int4    FindHBondsPDB2(FILE *fp,res_typ Res,Int4 C,float dmax,
                unsigned short file_id, set_typ **skip, char **Color,pdb_typ P,
		char mode=0);
void    PutContactsPDB2(FILE *fptr, Int4 file_id, float HA_dmax, float dmax,
                set_typ **skip,char **Color,pdb_typ P, char mode=0);

#endif

