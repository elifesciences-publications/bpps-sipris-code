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

/****************************************************************
Amino Acid 'A' solvent accessability: ( '=' is 32 counts. )

<   5.00 : 1586     |=================================================
    5.00 : 1099     |==================================
>= 28.00 : 1297     |========================================
   total : 3982    mean = 21.0892 stdev = 23.5987 range = 0 .. 102.173

Amino Acid 'C' solvent accessability: ( '=' is 8 counts. )

<   5.00 : 366      |=============================================
    5.00 : 218      |===========================
>= 28.00 : 99       |============
   total : 683     mean = 11.2024 stdev = 15.134 range = 0 .. 88.2605

Amino Acid 'D' solvent accessability: ( '=' is 35 counts. )

<   5.00 : 347      |=========
    5.00 : 768      |=====================
>= 28.00 : 1705     |================================================
   total : 2820    mean = 34.874 stdev = 22.6047 range = 0 .. 97.7439

Amino Acid 'E' solvent accessability: ( '=' is 37 counts. )

<   5.00 : 238      |======
    5.00 : 660      |=================
>= 28.00 : 1847     |=================================================
   total : 2745    mean = 38.5547 stdev = 22.1733 range = 0 .. 107.064

Amino Acid 'F' solvent accessability: ( '=' is 19 counts. )

<   5.00 : 924      |================================================
    5.00 : 645      |=================================
>= 28.00 : 278      |==============
   total : 1847    mean = 12.4171 stdev = 16.9031 range = 0 .. 93.3904

Amino Acid 'G' solvent accessability: ( '=' is 28 counts. )

<   5.00 : 1194     |==========================================
    5.00 : 1139     |========================================
>= 28.00 : 1353     |================================================
   total : 3686    mean = 23.9412 stdev = 23.9643 range = 0 .. 96.6216

Amino Acid 'H' solvent accessability: ( '=' is 9 counts. )

<   5.00 : 237      |==========================
    5.00 : 407      |=============================================
>= 28.00 : 375      |=========================================
   total : 1019    mean = 24.1886 stdev = 21.3222 range = 0 .. 92.0618

Amino Acid 'I' solvent accessability: ( '=' is 28 counts. )

<   5.00 : 1382     |=================================================
    5.00 : 690      |========================
>= 28.00 : 345      |============
   total : 2417    mean = 11.0991 stdev = 16.0551 range = 0 .. 99.6408

Amino Acid 'K' solvent accessability: ( '=' is 43 counts. )

<   5.00 : 89       |==
    5.00 : 603      |==============
>= 28.00 : 2108     |=================================================
   total : 2800    mean = 42.5456 stdev = 19.9517 range = 0 .. 100.92

Amino Acid 'L' solvent accessability: ( '=' is 40 counts. )

<   5.00 : 1990     |=================================================
    5.00 : 1239     |==============================
>= 28.00 : 632      |===============
   total : 3861    mean = 12.5608 stdev = 17.2373 range = 0 .. 101.676

Amino Acid 'M' solvent accessability: ( '=' is 10 counts. )

<   5.00 : 471      |===============================================
    5.00 : 298      |=============================
>= 28.00 : 158      |===============
   total : 927     mean = 13.4547 stdev = 17.9971 range = 0 .. 93.0394

Amino Acid 'N' solvent accessability: ( '=' is 23 counts. )

<   5.00 : 359      |===============
    5.00 : 670      |=============================
>= 28.00 : 1117     |================================================
   total : 2146    mean = 31.0847 stdev = 22.8935 range = 0 .. 94.531

Amino Acid 'P' solvent accessability: ( '=' is 21 counts. )

<   5.00 : 407      |===================
    5.00 : 640      |==============================
>= 28.00 : 1038     |=================================================
   total : 2085    mean = 30.2343 stdev = 23.1098 range = 0 .. 93.7472

Amino Acid 'Q' solvent accessability: ( '=' is 20 counts. )

<   5.00 : 234      |===========
    5.00 : 481      |========================
>= 28.00 : 978      |================================================
   total : 1693    mean = 33.0613 stdev = 21.6731 range = 0 .. 97.6207

Amino Acid 'R' solvent accessability: ( '=' is 22 counts. )

<   5.00 : 192      |========
    5.00 : 728      |=================================
>= 28.00 : 1056     |================================================
   total : 1976    mean = 31.8264 stdev = 20.722 range = 0 .. 93.1175

Amino Acid 'S' solvent accessability: ( '=' is 25 counts. )

<   5.00 : 699      |===========================
    5.00 : 796      |===============================
>= 28.00 : 1218     |================================================
   total : 2713    mean = 27.4923 stdev = 23.8391 range = 0 .. 94.8263

Amino Acid 'T' solvent accessability: ( '=' is 25 counts. )

<   5.00 : 708      |============================
    5.00 : 824      |================================
>= 28.00 : 1200     |================================================
   total : 2732    mean = 26.0043 stdev = 22.2283 range = 0 .. 100.432

Amino Acid 'V' solvent accessability: ( '=' is 33 counts. )

<   5.00 : 1622     |=================================================
    5.00 : 918      |===========================
>= 28.00 : 548      |================
   total : 3088    mean = 12.6802 stdev = 17.4367 range = 0 .. 93.4012

Amino Acid 'W' solvent accessability: ( '=' is 6 counts. )

<   5.00 : 291      |================================================
    5.00 : 275      |=============================================
>= 28.00 : 103      |=================
   total : 669     mean = 13.1457 stdev = 15.9236 range = 0 .. 84.4212

Amino Acid 'Y' solvent accessability: ( '=' is 16 counts. )

<   5.00 : 548      |==================================
    5.00 : 752      |===============================================
>= 28.00 : 383      |=======================
   total : 1683    mean = 17.5715 stdev = 17.7842 range = 0 .. 96.7243

Overall Amino Acid solvent accessability: ( '=' is 357 counts. )

<   5.00 : 13884    |======================================
    5.00 : 13850    |======================================
>= 28.00 : 17838    |=================================================
   total : 45572   mean = 24.3641 stdev = 23.0368 range = 0 .. 107.064
number of bins = 3 (s:5,e:28,i23)
Mutual Information I(A,B) = 0.105884 nats
 ***********************************************************************/
/* exposed.h = exposed surface likelihood ratios */
#if !defined (EXPOSED_SURFACE)
#define EXPOSED_SURFACE

static char exposed_character[5] = ".BIS";

static double percent_exposed[4] = {0, 5, 28, 999999};

static double exposed_info[21][4] = { {0, 0,0,0}, /** X **/
{ 0,  0.565, 0.049, -0.993 }, /** C **/
{ 0,  0.061, 0.017, -0.064 }, /** G **/
{ 0, 0.268, -0.096, -0.184 }, /** A **/
{ 0, -0.168, -0.035, 0.137 }, /** S **/
{ 0, -0.162, -0.008, 0.115 }, /** T **/
{ 0,  -0.599, 0.027, 0.285 }, /** N **/
{ 0, -0.907, -0.110, 0.435 }, /** D **/
{ 0, -1.257, -0.234, 0.542 }, /** E **/
{ 0, -0.790, -0.067, 0.389 }, /** Q **/
{ 0, -2.260, -0.344, 0.654 }, /** K **/
{ 0,  -1.143, 0.192, 0.311 }, /** R **/
{ 0, -0.270, 0.273, -0.062 }, /** H **/
{ 0,  0.356, 0.302, -0.933 }, /** W **/
{ 0,  0.066, 0.385, -0.542 }, /** Y **/ 
{ 0,  0.496, 0.139, -0.956 }, /** F **/
{ 0, 0.545, -0.022, -0.791 }, /** V **/
{ 0, 0.630, -0.063, -1.009 }, /** I **/
{ 0,  0.526, 0.054, -0.872 }, /** L **/
{ 0,  0.511, 0.056, -0.831 }, /** M **/
{ 0,  -0.445, 0.010, 0.240 }  /** P **/
};

/** XCGASTNDEQKRHWYFVILMP **/
static double maximum_exposed[21] = { 0,  /** X **/
 143.933, /** 'C' **/ 
  96.318, /** 'G' **/ 
 119.993, /** 'A' **/ 
 147.224, /** 'S' **/ 
 165.783, /** 'T' **/ 
 194.498, /** 'N' **/ 
 190.149, /** 'D' **/ 
 217.879, /** 'E' **/ 
 219.730, /** 'Q' **/ 
 228.293, /** 'K' **/ 
 270.874, /** 'R' **/ 
 208.623, /** 'H' **/ 
 274.155, /** 'W' **/ 
 249.471, /** 'Y' **/ 
 232.435, /** 'F' **/ 
 166.971, /** 'V' **/ 
 192.113, /** 'I' **/ 
 194.673, /** 'L' **/ 
 208.215, /** 'M' **/ 
 159.369  /** 'P' **/ 
 };

#endif

