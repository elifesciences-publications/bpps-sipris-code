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

static Int4 heat_pos[]  =
{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34};

Int4 heat_rpts = 19; 

static char heat_oper[] = "E\
Dmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm\
MmmmmmmmmmmmmmmIImmmmmmmmmmmmmmmmmmmiiiiiiiiiiiiiiiiiii\
MmmmmmmmmmmmmmmIIIIdddddmmmmmmmmmmmmmmiiiiiii\
DdmmmmmmmmmmmmmIImmmmmmmmmmmmmmmmmmmiiiiiiii\
MmmmmmmmmmmmmmmIImmmmmmmmmmmmmmmmmmmiiiiiiiii\
Mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmiiiiiiii\
Mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmiiiiiiiiiiiiiiiiiiiiiiiiiiii\
MmmmmmmmmmmmmmmIIIIIIIIIddmmmmmmmmmmmmmmmmmiiii\
Mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmiiiiiiii\
Mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmiiiiiiiii\
Mmmmmmmmmmmmmmdmmmmmmmmmmmmmmmmmmmiiiiiiiiiiiiiiiiiiiii\
MmmmmmmmmmmmmmmIIIIIddmmmmmmmmmmmmmmmmmiiiiiiiiiii\
MmmmmmmmmmmmmmmIIIIIIIIIIImmmmmmmmmmmmmmmmmmmiiiiiiii\
MmmmmmmmmmmmmmmIImmmmmmmmmmmmmmmmmmmiiiiiii\
MmmmmmmmmmmmmmmImmmmmmmmmmmmmmmmmmmiiiiiii\
MmmmmmmmmmmmmmmIImmmmmmmmmmmmmmmmmmmiiiiiii\
MmmmmmmmmmmmmmmIIIIIIIIIIImmmmmmmmmmmmmmmmmmmiiiiiiiiiiiiiiiii\
MmmmmmmmmmmmmmmIImmmmmmmmmmmmmmmmmmmiiiiiiiiii\
MmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmE";


static Int4 mouse_pos[]  =
{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34};

Int4 mouse_rpts = 10; 

static char mouse_oper[] = "Eiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii\
Mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmiiiiiiiii\
MmmmmmmmmmmmmmImmmmmmmmmmmmmmmmmmmdiiiiiiiii\
Mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmdiiiiiiiii\
MmmmmmmmmmmmmIIIIImmmmmmmmmmmmmmmmmmmmmiiiiiiiii\
Dmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmdiiiiiiiii\
Mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmdiiiiiiiii\
Mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmdiiiiiiiii\
Mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmdiiiiiiiiii\
Mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmiiiiiiiiiiiiii\
MmmmmmmmmmmmmmmmmmmmmmmmmmmmmmddddE";


static Int4 yeast_pos[]  =
{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34};

Int4 yeast_rpts = 10;

static char yeast_oper[] = "E\
MmmdddmmmmmmmmmIIImmmmmmmmmmmmmmmmmmmiiiiiiiii\
MmmmmmmmmmmmmmImmmmmmmmmmmmmmmmmmmdiiiiiiiii\
Mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmdiiiiiiiii\
MmmmmmmmmmmmdImmmmmmmmmmmmmmmmmmmmmiiiiiiiii\
Dmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmdiiiiiiiii\
Mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmdiiiiiiiii\
Mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmdiiiiiiiii\
Mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmiiiiiiiiiii\
Mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmiiiiiiiiiiiiiiiiiii\
MmmmmmmmmmmmmmmmmmmmmmmmmmmmmmddddE";


static Int4 arm_pos[]  =
{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34};

Int4 arm_rpts = 12;

static char arm_oper[] = "E\
Mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmiiiiiiii\
MmmmmmmmmmmmmmImmmmmmmmmmmmmmmmmmmmiiiiiii\
Mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmiiiiiiii\
Mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmiiiiiiii\
MmmmmmmmmmmmmmmImmmmmmmmmmmmmmmmmmmiiiiiii\
Mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmiiiii\
Mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmiiiiiiii\
MmmmmmmmmmmmmmmIImmmmmmmmmmmmmmmmmmmiiiiiiiiiii\
MmmmmmmmmmmmmmmImmmmmmmmmmmmmmmmmmmiiiiiii\
MmmmmmmmmmmmmmIIIIIIIIIIIIIIIIIIIIIImmmmmmmmmmmmmmmmmmmmiiiiiii\
Mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmiiiiiii\
MmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmE";

