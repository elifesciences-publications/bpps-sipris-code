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

#if !defined (BLOSUM45H)
#define BLOSUM45H

static float blosum45L[21][21] = {
/*    C           G           A           S           T       
      N           D           E           Q           K       
      R           H           W           Y           F       
      V           I           L           M           P       */
{1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.},
{1.,   17.0902,    0.5568,    0.8004,    0.7969,    0.8218, 
   0.6802,    0.5328,    0.5454,    0.4862,    0.5467, 
   0.4718,    0.4906,    0.3336,    0.4889,    0.6017, 
   0.7148,    0.5428,    0.6727,    0.6037,    0.4115, },
{1.,   0.5568,    5.0708,    1.0803,    1.0580,    0.6928, 
   1.0586,    0.7404,    0.5764,    0.6870,    0.6779, 
   0.5702,    0.6621,    0.5913,    0.5492,    0.4803, 
   0.4792,    0.4163,    0.4501,    0.5847,    0.7022, },
{1.,   0.8004,    1.0803,    2.9504,    1.3004,    1.0006, 
   0.7891,    0.6889,    0.8252,    0.8667,    0.7862, 
   0.6997,    0.6541,    0.5655,    0.6387,    0.5874, 
   1.0100,    0.7468,    0.7124,    0.8214,    0.7086, },
{1.,   0.7969,    1.0580,    1.3004,    2.7820,    1.4725, 
   1.1997,    0.9294,    0.9121,    1.0919,    0.8900, 
   0.7991,    0.8544,    0.4284,    0.7059,    0.6095, 
   0.7278,    0.6184,    0.5557,    0.6603,    0.7498, },
{1.,   0.8218,    0.6928,    1.0006,    1.4725,    3.1388, 
   1.1147,    0.8761,    0.8329,    0.7806,    0.8845, 
   0.7151,    0.7062,    0.4541,    0.7434,    0.7161, 
   1.0402,    0.8475,    0.7807,    0.8605,    0.8562, },
{1.,   0.6802,    1.0586,    0.7891,    1.1997,    1.1147, 
   4.4781,    1.5017,    0.8932,    1.0325,    1.1188, 
   0.8978,    1.2201,    0.3897,    0.6058,    0.5724, 
   0.5516,    0.5642,    0.4837,    0.6817,    0.6404, },
{1.,   0.5328,    0.7404,    0.6889,    0.9294,    0.8761, 
   1.5017,    5.3558,    1.6426,    0.9577,    0.9417, 
   0.7712,    0.9762,    0.3734,    0.6449,    0.4310, 
   0.4936,    0.4403,    0.4625,    0.4944,    0.7240, },
{1.,   0.5454,    0.5764,    0.8252,    0.9121,    0.8329, 
   0.8932,    1.6426,    3.8734,    1.5309,    1.2769, 
   1.0107,    0.9617,    0.5193,    0.6166,    0.4978, 
   0.5552,    0.4853,    0.5708,    0.6148,    0.9107, },
{1.,   0.4862,    0.6870,    0.8667,    1.0919,    0.7806, 
   1.0325,    0.9577,    1.5309,    4.4073,    1.3301, 
   1.3291,    1.1505,    0.6452,    0.8292,    0.4436, 
   0.5473,    0.5779,    0.6420,    0.9404,    0.7155, },
{1.,   0.5467,    0.6779,    0.7862,    0.8900,    0.8845, 
   1.1188,    0.9417,    1.2769,    1.3301,    3.3272, 
   1.9426,    0.8903,    0.5616,    0.7371,    0.5292, 
   0.5916,    0.5323,    0.5536,    0.7380,    0.7812, },
{1.,   0.4718,    0.5702,    0.6997,    0.7991,    0.7151, 
   0.8978,    0.7712,    1.0107,    1.3291,    1.9426, 
   4.7469,    0.9731,    0.5802,    0.8074,    0.5896, 
   0.5777,    0.4884,    0.6012,    0.7759,    0.5816, },
{1.,   0.4906,    0.6621,    0.6541,    0.8544,    0.7062, 
   1.2201,    0.9762,    0.9617,    1.1505,    0.8903, 
   0.9731,    9.5123,    0.4518,    1.4720,    0.6791, 
   0.4567,    0.4533,    0.6698,    0.9181,    0.6612, },
{1.,   0.3336,    0.5913,    0.5655,    0.4284,    0.4541, 
   0.3897,    0.3734,    0.5193,    0.6452,    0.5616, 
   0.5802,    0.4518,    29.7022,    1.8010,    1.3549, 
   0.4733,    0.5649,    0.6709,    0.6343,    0.5250, },
{1.,   0.4889,    0.5492,    0.6387,    0.7059,    0.7434, 
   0.6058,    0.6449,    0.6166,    0.8292,    0.7371, 
   0.8074,    1.4720,    1.8010,    5.7533,    2.1845, 
   0.8093,    0.9065,    0.9651,    1.0231,    0.4788, },
{1.,   0.6017,    0.4803,    0.5874,    0.6095,    0.7161, 
   0.5724,    0.4310,    0.4978,    0.4436,    0.5292, 
   0.5896,    0.6791,    1.3549,    2.1845,    5.7482, 
   0.9525,    1.0638,    1.3030,    1.0625,    0.4512, },
{1.,   0.7148,    0.4792,    1.0100,    0.7278,    1.0402, 
   0.5516,    0.4936,    0.5552,    0.5473,    0.5916, 
   0.5777,    0.4567,    0.4733,    0.8093,    0.9525, 
   2.8707,    2.1760,    1.3337,    1.2358,    0.5400, },
{1.,   0.5428,    0.4163,    0.7468,    0.6184,    0.8475, 
   0.5642,    0.4403,    0.4853,    0.5779,    0.5323, 
   0.4884,    0.4533,    0.5649,    0.9065,    1.0638, 
   2.1760,    3.2326,    1.5962,    1.4553,    0.6096, },
{1.,   0.6727,    0.4501,    0.7124,    0.5557,    0.7807, 
   0.4837,    0.4625,    0.5708,    0.6420,    0.5536, 
   0.6012,    0.6698,    0.6709,    0.9651,    1.3030, 
   1.3337,    1.5962,    2.9972,    1.7315,    0.4779, },
{1.,   0.6037,    0.5847,    0.8214,    0.6603,    0.8605, 
   0.6817,    0.4944,    0.6148,    0.9404,    0.7380, 
   0.7759,    0.9181,    0.6343,    1.0231,    1.0625, 
   1.2358,    1.4553,    1.7315,    4.1142,    0.6437, },
{1.,   0.4115,    0.7022,    0.7086,    0.7498,    0.8562, 
   0.6404,    0.7240,    0.9107,    0.7155,    0.7812, 
   0.5816,    0.6612,    0.5250,    0.4788,    0.4512, 
   0.5400,    0.6096,    0.4779,    0.6437,    8.8189, }
};

#endif

