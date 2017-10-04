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

#include "n2a_typ.h"

static void	PutN2AcodonAA(FILE *fp,Int4 b1,Int4 b2, Int4 b3, Int4 ***table,
			a_type AB, a_type dnaAB)
{
	fprintf(fp,"table[%c][%c][%c]=AlphaCode('%c',AB);\n",
		AlphaChar(b1,dnaAB), AlphaChar(b2,dnaAB),
		AlphaChar(b3,dnaAB),AlphaChar(table[b1][b2][b3],AB));
}

n2a_typ::n2a_typ(Int4 Code, Int4 MinLen, a_type aaAB, a_type naAB)
{ Init(Code,MinLen,aaAB,naAB); }

void	n2a_typ::Init(Int4 Code, Int4 MinLen, a_type aaAB, a_type naAB)
{
	if(nAlpha(aaAB) != 20 || nAlpha(naAB) != 4) print_error("MkN2A( ): inconsistent alphabet");
	for(Int4 rf=1; rf <=NumReadFrames( ); rf++){ num[rf]=0; E[rf]=NULL; }
	code = Code;
	dnaE=NULL;
	AB = aaAB; 
	dnaAB = naAB; 
	min_len=MinLen;
	N=AlphaCode('N',dnaAB);
	A=AlphaCode('A',dnaAB); C=AlphaCode('C',dnaAB);
	G=AlphaCode('G',dnaAB); T=AlphaCode('T',dnaAB);
        Int4	i,j,k;

        NEWPP(table,nAlpha(dnaAB)+1,Int4);
        for(i=0;i<=nAlpha(dnaAB);i++) {
                NEWP(table[i],nAlpha(dnaAB)+1,Int4);
                for(j=0;j<=nAlpha(dnaAB);j++)
                        NEW(table[i][j],nAlpha(dnaAB)+1,Int4);
        }
        for(j=0;j<=nAlpha(dnaAB);j++){
           for(k=0;k<=nAlpha(dnaAB);k++){
               table[N][j][k]=AlphaCode('X',AB);;
               table[j][N][k]=AlphaCode('X',AB);;
           }
        }
	table=set_table( );
}

Int4	***n2a_typ::set_table( )
{

	table[A][A][N]=AlphaCode('X',AB);
	table[A][A][A]=AlphaCode('K',AB); table[A][A][C]=AlphaCode('N',AB);
	table[A][A][G]=AlphaCode('K',AB); table[A][A][T]=AlphaCode('N',AB);

	table[A][C][N]=AlphaCode('T',AB);
	table[A][C][A]=AlphaCode('T',AB); table[A][C][C]=AlphaCode('T',AB);
	table[A][C][G]=AlphaCode('T',AB); table[A][C][T]=AlphaCode('T',AB);

	table[A][G][N]=AlphaCode('X',AB);
	table[A][G][A]=AlphaCode('R',AB); table[A][G][C]=AlphaCode('S',AB);
	table[A][G][G]=AlphaCode('R',AB); table[A][G][T]=AlphaCode('S',AB);

	table[A][T][N]=AlphaCode('X',AB);
	table[A][T][A]=AlphaCode('I',AB); table[A][T][C]=AlphaCode('I',AB);
	table[A][T][G]=AlphaCode('M',AB); table[A][T][T]=AlphaCode('I',AB);

	table[C][A][N]=AlphaCode('X',AB);
	table[C][A][A]=AlphaCode('Q',AB); table[C][A][C]=AlphaCode('H',AB);
	table[C][A][G]=AlphaCode('Q',AB); table[C][A][T]=AlphaCode('H',AB);

	table[C][C][N]=AlphaCode('P',AB);
	table[C][C][A]=AlphaCode('P',AB); table[C][C][C]=AlphaCode('P',AB);
	table[C][C][G]=AlphaCode('P',AB); table[C][C][T]=AlphaCode('P',AB);

	table[C][G][N]=AlphaCode('R',AB);
	table[C][G][A]=AlphaCode('R',AB); table[C][G][C]=AlphaCode('R',AB);
	table[C][G][G]=AlphaCode('R',AB); table[C][G][T]=AlphaCode('R',AB);

	table[C][T][N]=AlphaCode('L',AB);
	table[C][T][A]=AlphaCode('L',AB); table[C][T][C]=AlphaCode('L',AB);
	table[C][T][G]=AlphaCode('L',AB); table[C][T][T]=AlphaCode('L',AB);

	table[G][A][N]=AlphaCode('X',AB);
	table[G][A][A]=AlphaCode('E',AB); table[G][A][C]=AlphaCode('D',AB);
	table[G][A][G]=AlphaCode('E',AB); table[G][A][T]=AlphaCode('D',AB);

	table[G][C][N]=AlphaCode('A',AB);
	table[G][C][A]=AlphaCode('A',AB); table[G][C][C]=AlphaCode('A',AB);
	table[G][C][G]=AlphaCode('A',AB); table[G][C][T]=AlphaCode('A',AB);

	table[G][G][N]=AlphaCode('G',AB);
	table[G][G][A]=AlphaCode('G',AB); table[G][G][C]=AlphaCode('G',AB);
	table[G][G][G]=AlphaCode('G',AB); table[G][G][T]=AlphaCode('G',AB);

	table[G][T][N]=AlphaCode('V',AB);
	table[G][T][A]=AlphaCode('V',AB); table[G][T][C]=AlphaCode('V',AB);
	table[G][T][G]=AlphaCode('V',AB); table[G][T][T]=AlphaCode('V',AB);

	table[T][A][N]=AlphaCode('X',AB);
	table[T][A][A]=AlphaCode('-',AB); table[T][A][C]=AlphaCode('Y',AB);
	table[T][A][G]=AlphaCode('-',AB); table[T][A][T]=AlphaCode('Y',AB);

	table[T][C][N]=AlphaCode('S',AB);
	table[T][C][A]=AlphaCode('S',AB); table[T][C][C]=AlphaCode('S',AB);
	table[T][C][G]=AlphaCode('S',AB); table[T][C][T]=AlphaCode('S',AB);

	table[T][G][N]=AlphaCode('X',AB);
	table[T][G][A]=AlphaCode('-',AB); table[T][G][C]=AlphaCode('C',AB);
	table[T][G][G]=AlphaCode('W',AB); table[T][G][T]=AlphaCode('C',AB);

	table[T][T][N]=AlphaCode('X',AB);
	table[T][T][A]=AlphaCode('L',AB); table[T][T][C]=AlphaCode('F',AB);
	table[T][T][G]=AlphaCode('L',AB); table[T][T][T]=AlphaCode('F',AB);

        switch (code){
                case  1: // The Standard Code:
			break;
                case  2: // The Vertebrate Mitochondrial Code:
			table[A][G][A]=table[A][G][G]=AlphaCode('X',AB); // = Termination!
			table[A][T][A]=AlphaCode('M',AB);
			table[T][G][A]=AlphaCode('W',AB);	// Not termination!
                        break;	// Vertebrata
                case  3: // The Yeast Mitochondrial Code
			table[A][T][A]=AlphaCode('M',AB); // AUA    Met  M 
                 	table[C][T][T]=AlphaCode('T',AB); 
                 	table[C][T][C]=AlphaCode('T',AB); 
                 	table[C][T][A]=AlphaCode('T',AB); 
                 	table[C][T][G]=AlphaCode('T',AB); 
                 	table[C][T][N]=AlphaCode('T',AB); 
			table[T][G][A]=AlphaCode('W',AB);	// Not termination!
			table[C][G][A]=AlphaCode('X',AB);	// absent! == seq. error?
			table[C][G][C]=AlphaCode('X',AB);	// absent! == seq. error?
                        break;	// Saccharomyces cerevisiae, Candida glabrata, Hansenula saturnus, 
				// and Kluyveromyces thermotolerans
                case  4: // The Mold, Protozoan, and Coelenterate Mitochondrial Code 
			 // and the Mycoplasma/Spiroplasma Code
			table[T][G][A]=AlphaCode('W',AB);	// Not termination!
                        break;	// Mycoplasmatales: Mycoplasma, Spiroplasma (Bove et al., 1989);
                case  5: // The Invertebrate Mitochondrial Code:
			table[A][G][N]=table[A][G][A]=table[A][G][G]=AlphaCode('S',AB); // 
			table[A][T][A]=AlphaCode('M',AB);
			table[T][G][A]=AlphaCode('W',AB);	// Not termination!
                        break;	// Nematoda: Ascaris, Caenorhabditis;
                		// Mollusca: Bivalvia (Hoffmann et al., 1992); 
				// Polyplacophora (Boore and Brown, 1994)
                		// Arthropoda/Crustacea: Artemia (Batuecas et al., 1988);
                		// Arthropoda/Insecta: Drosophila 
                case  6: // The Ciliate, Dasycladacean and Hexamita Nuclear Code
			table[T][A][A]=table[T][A][G]=AlphaCode('Q',AB);	// Not termination!
                        break;
                // case  7: case 8: // deleted!
                case  9: // The Echinoderm Mitochondrial Code
			table[A][A][A]=AlphaCode('N',AB);	// 
			table[A][G][A]=table[A][G][G]=table[A][G][N]=AlphaCode('S',AB);
			table[T][G][A]=AlphaCode('W',AB);	// Not termination!
                        break;
                case 10: // The Euplotid Nuclear Code 
			table[T][G][A]=AlphaCode('C',AB);	// Not termination!
                        break;	// Ciliata: Euplotidae.
                case 11: // The Bacterial and Plant Plastid Code
			break;	// same as standard code.
                case 12: 	// The Alternative Yeast Nuclear Code
			table[C][T][G]=AlphaCode('S',AB);	// Not Leu!
                        break;	// Endomycetales (yeasts): Candida albicans, Candida cylindracea, 
				// Candida melibiosica, Candida parapsilosis, 
				// and Candida rugosa (Ohama et al., 1993). 
                case 13: 	// The Ascidian Mitochondrial Code
			table[A][G][A]=table[A][G][G]=AlphaCode('G',AB); // 
			table[A][T][A]=AlphaCode('M',AB);
			table[T][G][A]=AlphaCode('W',AB);	// Not termination!
                        break;	// Ascidiacea (sea squirts): Pyuridae (Durrheim et al. 1993; Yokobori et al., 1993)
                case 14:	//  The Flatworm Mitochondrial Code
			table[A][A][A]=AlphaCode('N',AB);	// 
			table[A][G][A]=table[A][G][G]=table[A][G][N]=AlphaCode('S',AB);
			table[T][A][A]=AlphaCode('Y',AB);	// Not termination!
			table[T][G][A]=AlphaCode('W',AB);	// Not termination!
                        break;	// Platyhelminthes (flatworms) 
                case 15: 	// Blepharisma Nuclear Code 
			table[T][A][G]=AlphaCode('Q',AB);	// Not termination!
                        break; // Ciliata: Blepharisma (Liang and Heckman, 1993) 
                case 16:	// Chlorophycean Mitochondrial Code
			break;
                case 21:	//  Trematode Mitochondrial Code
			break;
                default: print_error("no such conversion table");
        } return table;
}

e_type	*n2a_typ::AllTranslate(Int4 *NumRFs,e_type DnaE)
// return only the single best translated sequence. 
{
    Int4 Num=0;
    if((Num=Translate(DnaE)) > 0){
	e_type *EList=0;
	Int4	rf,i,m;
	NEW(EList,Num+5,e_type);
        for(m=0,rf=1; rf<=NumReadFrames( ); rf++){
           for(i=1; i <= SeqPerReadFrame(rf); i++){
              m++; EList[m] = Sequence(rf,i); 
           }
        }
	*NumRFs=m;
	return EList;
    } else { *NumRFs=0; return NULL; }
}

e_type	n2a_typ::BestTranslate(e_type qE, double *best_score, e_type DnaE)
// return only the single best translated sequence. 
{
    Int4 Num=0;
    if((Num=Translate(DnaE)) > 0){
	e_type *EList=0;
	Int4	rf,i,m;
	NEW(EList,Num+5,e_type);
        for(m=0,rf=1; rf<=NumReadFrames( ); rf++){
           for(i=1; i <= SeqPerReadFrame(rf); i++){
              m++; EList[m] = Sequence(rf,i); 
           }
        }
	assert(Num==m);
	char	name[]="translated";
	ss_type data=Array2SeqSet(EList,Num, name,AB); // WARNING: EList absorbed!!
        double  Ethresh=1.0,Ecutoff=1.0;
        BooLean top_only=FALSE;
	double	*scores;
	Int4	Threshold=11,x_parameter=15.0;

        gpsi_type *gpsi = new gpsi_type(qE,data,Ethresh,Ecutoff,x_parameter,Num,1);
        scores=gpsi->FastSearch(top_only,Num,Threshold,'S');
	double	best=0.0; e_type bE=0;
        for(m = 1; m <= Num; m++){
	    // fprintf(stderr,"%4.0f ",scores[m]);
            if(scores[m] > best){ /** then merge sets **/
		best=scores[m]; bE=SeqSetE(m,data);
            }
        }
	if(bE) bE = CopySeq(bE);
	delete gpsi;
	EList=NilSeqSetRtnSeqs(data);
	free(EList); free(scores);
	*best_score=best;
	// fprintf(stderr,"\n");
	return bE;
    } else return NULL;
}

void	n2a_typ::reset( )
{
    if(dnaE != NULL){
	for(Int4 rf=1; rf <=NumReadFrames( ); rf++){
	     for(Int4 n = 1; n <= num[rf]; n++){
		NilSeq(E[rf][n]);
	     } free(E[rf]);
	} if(dnaE) NilSeq(dnaE);
    } dnaE=NULL;   
}

void	n2a_typ::Free( )
{
	Int4	i,j;
	reset( );
	for(i=0;i<=4;i++){
		for(j=0;j<=4;j++) free(table[i][j]); 
		free(table[i]);
	} free(table);
}

void	n2a_typ::reverse_compl( )
{
        Int4    i, j;
	char a;

        ReverseSeq(dnaE);
	unsigned char *ptr = SeqPtr(dnaE);

        for(i=1;i<=LenSeq(dnaE);i++){
		a = AlphaChar(ptr[i],dnaAB);
                if (a == 'A')  ptr[i] = AlphaCode('T',dnaAB);
                else if (a == 'T') ptr[i] = AlphaCode('A',dnaAB);
		else if (a == 'G') ptr[i] = AlphaCode('C',dnaAB);
		else if (a == 'C') ptr[i] = AlphaCode('G',dnaAB);		
        }
}

e_type	n2a_typ::Sequence(Int4 rf, Int4 i)
{
	if (rf > NumReadFrames() || rf < 1) 
		print_error("SeqN2A( ): wrong number of reading frames");
	else if(i > num[rf] || i <1) 
		print_error("SeqN2A( ): wrong number of seqs in a reading frame");
	return E[rf][i];
}

Int4	n2a_typ::Translate(e_type DnaE)
{
	if(DnaE == NULL) return 0;
	else {
	   reset( );
	   dnaE=CopySeq(DnaE);
	   for(Int4 rf=1;rf<=NumReadFrames( );rf++) 
		NEW(E[rf],LenSeq(dnaE)/min_len+5,e_type);
	   Int4 Total=translate(FALSE,DnaE);
	   Total += translate(TRUE,DnaE);
	   return Total;
	}
}

Int4	n2a_typ::translate(BooLean rev,e_type DnaE)
{
        Int4            start=1, end, rf, j, k, m, n, r, t, x, y, lng, lngth, len[4];
        unsigned char   a, *seq, *P;
        char            *st, str[5000], temp1[50], *temp2, *p, c;
	Int4		TotalTrans;

        NEW(st,5000,char);
        StrSeqInfo(st,DnaE);
        for(j=0;st[j]!=' ' && st[j] != '\n';j++) temp1[j]=st[j]; 
	assert(j > 0);
	if(temp1[j-1] == '|') j--;  // block out '|' at ends...(NCBI change 8/18/04)
	temp1[j]='\0';
// fprintf(stderr,"************ StrSeqInfo = %s **************\n temp1=%s\n",st,temp1);
        temp2=st; temp2+=j;

        if(rev) { reverse_compl( ); c='-'; x=4; y=6;}
        else {c='+'; x=1; y=3;}

        lngth = LenSeq(dnaE);

        NEW(P,lngth/3+5,unsigned char);
        seq = SeqPtr(dnaE);
        r = lngth%3;
        if (r == 0) {len[1] = lngth/3; len[2] = lngth/3-1; len[3] = lngth/3-1;}
        else if (r == 1) {len[1] = lngth/3; len[2] = lngth/3; len[3] = lngth/3-1;}
        else {len[1] = lngth/3; len[2] = lngth/3; len[3] = lngth/3;}

        for(TotalTrans=0,j=1,rf=x; rf<=y; j++, rf++){
                n=1;k=1;t=1;
                while(k<=len[j]){
                        m=3*(k-1)+j;
                        P[t]=table[seq[m]][seq[m+1]][seq[m+2]];
                        if(P[t]==21 || k==len[j]){
                                if(t > min_len){
                                        end=m-1;start=end-(3*t)+4;
                                        sprintf(str,"%s_EST (nt %c%d:%d)%s",
							temp1,c,j,start,temp2);
					// note: (calling function parses out '(nt ...)'
					if(P[t]!=21 && k==len[j]) lng=t;
					else lng=t-1;
// fprintf(stderr,"%s\n rf=%d; n=%d\n",str,rf,n);
                                        E[rf][n]=MkSeq(str,lng,P);
// PutSeq(stderr,E[rf][n],AB);
                                        TotalTrans++; n++;
                                } t=1;k++;
                        } else {k++;t++;}
                } num[rf]=n-1;
        } free(st); free(P);
        return TotalTrans;
}

