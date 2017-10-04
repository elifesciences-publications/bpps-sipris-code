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

#if !defined (BLOSUM35H)
#define BLOSUM35H

static float blosum35L[21][21] = {
{1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.},
{1.,   13.1930,    0.5644,    0.7570,    0.6166,    0.8482, 
   0.7817,    0.5928,    0.8665,    0.5459,    0.7227, 
   0.6140,    0.5105,    0.4584,    0.4290,    0.5277, 
   0.6527,    0.5297,    0.6982,    0.5434,    0.4751, },
{1.,   0.5644,    3.6138,    1.0698,    1.1659,    0.7013, 
   1.1792,    0.7675,    0.6758,    0.7609,    0.8520, 
   0.6859,    0.7128,    0.8565,    0.7540,    0.5939, 
   0.5877,    0.6053,    0.5968,    0.7974,    0.7632, },
{1.,   0.7570,    1.0698,    2.1886,    1.1788,    1.0434, 
   0.8925,    0.8654,    0.8654,    0.9660,    0.9293, 
   0.8516,    0.7328,    0.7295,    0.8002,    0.7008, 
   1.0554,    0.8580,    0.7494,    0.9890,    0.7542, },
{1.,   0.6166,    1.1659,    1.1788,    2.1723,    1.3547, 
   0.9725,    0.8818,    0.9487,    1.0864,    0.9666, 
   0.8963,    0.8963,    0.6715,    0.7824,    0.8002, 
   0.8223,    0.7674,    0.7248,    0.8070,    0.7117, },
{1.,   0.8482,    0.7013,    1.0434,    1.3547,    2.3475, 
   0.9869,    0.7847,    0.8641,    0.9232,    0.9874, 
   0.7437,    0.7704,    0.7661,    0.7389,    0.8164, 
   1.1685,    0.9155,    0.9486,    0.9432,    0.9956, },
{1.,   0.7817,    1.1792,    0.8925,    0.9725,    0.9869, 
   3.2014,    1.2949,    0.8482,    1.1011,    1.0039, 
   0.8492,    1.1859,    0.6554,    0.6801,    0.7756, 
   0.6506,    0.8704,    0.6539,    0.8704,    0.7284, },
{1.,   0.5928,    0.7675,    0.8654,    0.8818,    0.7847, 
   1.2949,    4.1658,    1.3267,    0.8752,    0.9048, 
   0.8297,    0.9178,    0.6178,    0.6590,    0.5705, 
   0.6762,    0.5775,    0.6806,    0.5730,    0.8156, },
{1.,   0.8665,    0.6758,    0.8654,    0.9487,    0.8641, 
   0.8482,    1.3267,    2.7787,    1.5044,    1.2567, 
   0.8515,    0.8680,    0.9090,    0.7762,    0.6160, 
   0.7102,    0.6081,    0.7994,    0.7151,    0.9482, },
{1.,   0.5459,    0.7609,    0.9660,    1.0864,    0.9232, 
   1.1011,    0.8752,    1.5044,    3.3272,    1.0494, 
   1.3867,    0.9085,    0.7988,    0.9641,    0.5272, 
   0.6228,    0.6600,    0.7211,    0.8634,    0.9465, },
{1.,   0.7227,    0.8520,    0.9293,    0.9666,    0.9874, 
   1.0039,    0.9048,    1.2567,    1.0494,    2.3001, 
   1.3694,    0.6972,    0.9518,    0.8175,    0.8148, 
   0.7280,    0.7230,    0.7041,    0.9870,    0.9514, },
{1.,   0.6140,    0.6859,    0.8516,    0.8963,    0.7437, 
   0.8492,    0.8297,    0.8515,    1.3867,    1.3694, 
   3.7586,    0.8509,    0.9609,    0.9205,    0.8170, 
   0.8269,    0.6405,    0.7273,    0.9876,    0.6851, },
{1.,   0.5105,    0.7128,    0.7328,    0.8963,    0.7704, 
   1.1859,    0.9178,    0.8680,    0.9085,    0.6972, 
   0.8509,    8.5132,    0.5391,    0.9754,    0.6120, 
   0.5366,    0.6372,    0.7433,    1.2364,    0.8285, },
{1.,   0.4584,    0.8565,    0.7295,    0.6715,    0.7661, 
   0.6554,    0.6178,    0.9090,    0.7988,    0.9518, 
   0.9609,    0.5391,    14.9721,    1.7053,    1.1825, 
   0.6621,    0.8140,    0.9666,    1.1309,    0.5200, },
{1.,   0.4290,    0.7540,    0.8002,    0.7824,    0.7389, 
   0.6801,    0.6590,    0.7762,    0.9641,    0.8175, 
   0.9205,    0.9754,    1.7053,    3.9493,    1.8203, 
   0.9613,    1.0273,    1.0804,    1.0806,    0.6299, },
{1.,   0.5277,    0.5939,    0.7008,    0.8002,    0.8164, 
   0.7756,    0.5705,    0.6160,    0.5272,    0.8148, 
   0.8170,    0.6120,    1.1825,    1.8203,    4.3304, 
   1.1093,    1.1275,    1.3384,    0.9217,    0.4732, },
{1.,   0.6527,    0.5877,    1.0554,    0.8223,    1.1685, 
   0.6506,    0.6762,    0.7102,    0.6228,    0.7280, 
   0.8269,    0.5366,    0.6621,    0.9613,    1.1093, 
   2.2422,    1.9082,    1.3261,    1.1289,    0.6424, },
{1.,   0.5297,    0.6053,    0.8580,    0.7674,    0.9155, 
   0.8704,    0.5775,    0.6081,    0.6600,    0.7230, 
   0.6405,    0.6372,    0.8140,    1.0273,    1.1275, 
   1.9082,    2.5235,    1.4800,    1.2356,    0.7830, },
{1.,   0.6982,    0.5968,    0.7494,    0.7248,    0.9486, 
   0.6539,    0.6806,    0.7994,    0.7211,    0.7041, 
   0.7273,    0.7433,    0.9666,    1.0804,    1.3384, 
   1.3261,    1.4800,    2.3203,    1.5684,    0.5949, },
{1.,   0.5434,    0.7974,    0.9890,    0.8070,    0.9432, 
   0.8704,    0.5730,    0.7151,    0.8634,    0.9870, 
   0.9876,    1.2364,    1.1309,    1.0806,    0.9217, 
   1.1289,    1.2356,    1.5684,    2.7643,    0.5709, },
{1.,   0.4751,    0.7632,    0.7542,    0.7117,    0.9956, 
   0.7284,    0.8156,    0.9482,    0.9465,    0.9514, 
   0.6851,    0.8285,    0.5200,    0.6299,    0.4732, 
   0.6424,    0.7830,    0.5949,    0.5709,    6.0202, }
};

#endif
