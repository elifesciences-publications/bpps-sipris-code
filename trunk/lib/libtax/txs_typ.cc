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

#include "txs_typ.h"

#if 1	// added txc_typ here...

char	**txc_typ::CopyStr(const char **str)
{
	Int4	n,i;
	char	**str2;

	for(n=0; str[n][0] != 0; n++) ;
	NEWP(str2, n+2, char);
	for(i=0; str[i][0] != 0; i++){
		str2[i] = AllocString(str[i]);
	} return str2;
}

void	txc_typ::Free(  )
{
       for(Int4 i=0; i < MAX_TAX_GROUPS; i++){
           if(Group[i]){
              for(Int4 j=0; Group[i][j]; j++){
                 free(Group[i][j]);
              } free(Group[i]);
           }
       }
}

txc_typ::txc_typ(  )
// TO RETREIVE this...use the command:
// grep 'const char' sqtax.l | sed 's/[[].*$//g' | 
// sed 's/^const char[^*]*\*/Group[60]=/' | sed 's/$/;/'
{
const char	*UnKnown[] = {"UnKnown",""};

/******************************* Eukaryotes **********************************/
// Mycetozoa
const char	*E1[]= {"Mycetozoa",""};
const char	*E1a[] = {"Dictyosteliida","Dictyostelia","Dictyostelida",
			"Dictyosteliales","dictyostelid cellular slime molds",""};
const char	*E1b[]= {"Myxogastria","Myxomycetes","Myxogastromycetidae",
				"Stemonitomycetidae",""};
const char	*E1c[]= {"Protostelida",""};

// diplomonads == no mitochondria?? == archezoa??
const char	*E2[] = {"Diplomonadida","diplomonads","Hexamitidae",""}; // e.g., Giardia 
const char	*E3[] = {"Entamoebidae","Entamoebida","Entamoeba",""}; // Entamoeba to be safe.

// Euglenozoa 
const char	*E4a[]= {"Kinetoplastida","Protomonadida","kinetoplasts",""};
const char	*E4b[]= {"Euglenida","Euglenophyta","euglenids",""};

// Oxymonadida
const char	*E5[]= {"Oxymonadida","Pyrsonymphidae",""};

// Parabasalidea
const char	*E6[]= {"Parabasalidea","Parabasalia","parabasalids",
			"parabasalians",""};
const char	*E6a[]= {"Hypermastigida","Holomastigotidae","Trichonymphidae",""};
const char	*E6b[]= {"Trichomonadida","trichomonads",""};

// Heterolobosea
const char	*E7[]= {"Heterolobosea","Lobosa","Acrasida","Schizopyrenida",""};

const char	*E8[]= {"Microsporidia","Microspora","Microsporea","Microsporida",
			"microsporidia",""};

const char	*E9[]= {"Granuloreticulosea","Foraminifera","foraminifers","forams",
			"Reticulomyxa",""};

#if 0	//**********************************************************************
Archaeprotista:
Includes free-living amitochondriate amoeboflagellate.
1. Pelomyxa palustris has no ER, Golgi, mitochondria, chromosomes or centrioles.
	It is multinucleated and has three types of bacterial endosymbionts.
#endif  //**********************************************************************
const char	*E10[]= {"Pelobiontida","pelobionts","Mastigamoebidae","Pelomyxidae",""};

// Includes Chlorarachnion
// Both the cryptomonads and the chlorarachniophytes, both of which retain a
// vestigial endosymbiont nucleus (nucleomorph) between the outer and the inner pair of chloroplast
// envelope membranes 
const char	*AD1[]={"Cercozoa","Cercomonadidae","Euglyphina","Heteromitidae",
			"Thaumatomonadidae","Chlorarachniophyceae","Chlorarachnida",
			"Chloroarachnida",""};

#if 0   //**********************************************************************
Includes the freshwater protozoon Reclinomonas americana, which is a 
heterotrophic flagellate whose mtDNA (69,034 base pairs) contains 
the largest collection of genes (97) so far identified in any mtDNA.
#endif  //**********************************************************************
const char	*AD2[]={"Reclinomonas","core jakobids",""};

/**************************** eukaryote crown group ********************************/

const char	*EC1[]= {"Acanthamoebidae",""};
// const char	*EC2[]= {"Alveolata","alveolates",""};
const char	*EC2a[]= {"Apicomplexa","apicomplexans","apicomplexa",""};
const char	*EC2b[]= {"Ciliophora","Ciliata","ciliates",""};
const char	*EC2c[]= {"Dinophyceae","Dinophyta","Dinophycidae","dinoflagellates",
			   "dinophytes",""};
const char	*EC2d[]= {"Perkinsea","Perkinsiea",""};
const char	*EC3[]= {"Cryptophyta","Cryptophyceae","Cryptomonadales","cryptomonads",""};
const char      *EC4[]= {"Glaucocystophyceae",""};
const char      *EC5[]= {"Haptophyceae","Haptophyta","Prymnesiophyta","Prymnesiophyceae",""};
const char	*EC6[]= {"Rhodophyta","Rhodophyceae","Rhodophycota","red algae",""};
const char	*EC7[]={"Stramenopiles","heterokonts",""};

const char	*EP1[]={"Chlorophyta","Chlorophycota","green algae","algae","Volvocales",""};
const char	*EP2[]={"Streptophyta","Charophyta/Embryophyta group",
			"charophyte/embryophyte group",
			"Embryophyta","land plants","higher plants","plants",""};
const char	*EP2a[]={"Bryata","nonvascular plants",
			"Bryophyta","Musci","Bryopsida","mosses","bryophytes",
 			 "Bryopsida","Polytrichopsida","hair mosses",
			 "Sphagnopsida","Takakiopsida",
			"Anthocerotophyta","hornworts","Anthocerotales",
			  "Anthocerotaceae","Anthoceros",
			"Marchantiophyta","liverworts","Marchantiales","Marchantiaceae",
			""};
const char	*EP2b[]={"Tracheophyta","vascular plants",
			  "Euphyllophyta","Filicophyta","ferns","Spermatophyta","seed plants",
			  "Lycopodiophyta","club mosses","Lycopodiopsida","clubmosses",
			""};
const char	*EP2c[]={"Characeae","stoneworts",""};
const char	*EP2d[]={"Mesostigmatophyceae","Mesostigmatales","Mesostigmataceae",
			  "Mesostigma",""};

/****************************** Metazoans **********************************/

// Metazoa: Bilateria: Acoelomata 
const char	*EM0a[]={"Platyhelminthes","flatworms",""};

/*************************************************************************
Myzostomida: a link between trochozoans and flatworms?
"Myzostomids are obligate symbiotic invertebrates associated with echinoderms..."
Myzostomids have acquired a unique anatomy that obscures their phylogenetic
affinities to other metazoans: they are incompletely segmented, parenchymous, 
acoelomate organisms with chaetae and a trochophore larva. "
Today, they are most often classified within annelids either as an 
aberrant family of polychaetes or as a separate class...  All our 
analyses congruently indicated that myzostomids are not annelids but 
suggested instead that they are more closely related to flatworms than to 
any trochozoan taxon."
Ref:  Proc R Soc Lond B Biol Sci 2000 Jul 22;267(1451):1383-92.
 *************************************************************************/
const char	*EM0b[]={"Myzostomida","Myzostomidae","Myzostoma","Myzostomum",""};

// Metazoa: Bilateria: Pseudocoelomata 
const char      *EM1[]={"Acanthocephala","thorny-headed worms","acanthocephalans",""};
const char      *EM2[]={"Gastrotricha","gastrotrichs",""};
const char	*EM3[]={"Nematoda","Nemata","roundworms",""};
const char      *EM4[]={"Nematomorpha","horsehair worms",""};
const char      *EM5[]={"Rotifera","rotifers",""};

// Metazoa: Bilateria: Coelomata - Deuterostomia:
const char	*EM6[]={"Chordata","chordates",""};
const char	*EM7[]={"Echinodermata","echinoderms",""};
const char	*EM8[]={"Hemichordata","hemichordates",""};

// Metazoa: Bilateria: Coelomata - Protostomia:
const char	*EM9[]={"Annelida","annelid worms",""};
const char	*EM10[]={"Brachiopoda","lampshells","brachiopods","Lophophorata",""};
const char	*EM11[]={"Bryozoa","Ectoprocta","bryozoans",""};
const char	*EM12[]={"Echiura","spoonworms",""};
const char	*EM13[]={"Mollusca","mollusks","molluscs",""};
const char	*EM14[]={"Nemertea","Nemertini","ribbon worms","nemertines","bootlace worms",""};

// Panarthropoda ...
const char	*EM15a[]={"Arthropoda",""};
const char	*EM15b[]={"Onychophora","velvet worms",""};
const char	*EM16[]={"Sipuncula","sipunculids",""};
const char	*EM17[]={"Vestimentifera","Obturata","vestimentiferans",""};
// 
const char	*EM18[]={"Cnidaria","cnidarians",""};
const char	*EM19[]={"Ctenophora","ctenophores",""};

const char	*EM20[]={"Porifera","Parazoa","sponges",""};

// Tardigrada == water bears
const char	*EM21[]={"Tardigrada","water bears","tardigrades",
			"tardigrades","Eutardigrada",""};


/****************************** Fungi **********************************/
// WARNING: Ascomycota includes diverse fungi!! 
// make sure this is below fission and budding yeast!

// Ascomycota (ascomycetes) [ 64877 ] 
const char	*F0[]={"Ascomycota","ascomycetes","sac fungi","mitosporic Ascomycota",
			"unclassified Ascomycota",""};
  // Pezizomycotina [ 13222 ] 
const char	*F0a[]={"Pezizomycotina","Arthoniomycetes","black yeasts","Chaetothyriomycetes",
			"Dothideomycetes","Eurotiomycetes", "Onygenales",
			"Eurotiales", "anamorphic Ascomycota", "Fusarium", 
			"Penicillium","green and blue molds",
			"Lecanoromycetes","Leotiomycetes","powdery mildews",
			"Pezizomycetes", "Sordariomycetes", "unclassified Pezizomycotina",
			"Hypocreales", "Microascales", "Ophiostomatales",
			"Dutch elm disease fungus",
			"Dothideomycetes et Chaetothyriomycetes incertae sedis", 
			"Diaporthales","Halosphaeriales",
              		"Sordariales", "Xylariales", "unclassified Sordariomycetes",
              		"Sordariomycetes incertae sedis","pyrenomycetes", // Neurospora crassa
			""};
  // Saccharomycotina [ 36289 ] 
const char	*F0b[]={"Saccharomycotina","Saccharomycetes","Saccharomycetales",
			"Endomycetales","Hemiascomycetes","budding yeasts",""};
  // Taphrinomycotina [ 15114 ] 
const char	*F0c[]={"Taphrinomycotina","Taphrinomycetes","Neolectomycetes",
			"Pneumocystidomycetes","Archiascomycetes","Schizosaccharomycetes",
			"Schizosaccharomycetales","Schizosaccharomycetaceae",
			"Schizosaccharomycetoideae","fission yeasts",""};

// Basidiomycota (basidiomycetes) [ 2500 ] 
const char	*F1[] = {"Basidiomycota","basidiomycetes",""};
// Chytridiomycota [ 372 ] 
const char	*F2[]={"Chytridiomycota","Mastigomycotina",""};
// Zygomycota [ 635 ] 
const char	*F3[]={"Zygomycota","pin molds",""};
// Fungi incertae sedis [ 52 ] 
const char	*F4[]={"Fungi incertae sedis",""};

/****************************** Archaea **********************************/
const char	*A0[] ={"Archaea","Mendosicutes","Metabacteria","Archaeobacteria",
			"archaebacteria","Procaryotae","Prokaryotae",""};
const char	*A1[] = {"Crenarchaeota","Eocyta",
			"extremely thermophilic archaebacteria",""};
const char	*A1a[] = {"Desulfurococcales","Aeropyrum","Desulfurococcus",
			"Staphylothermus",""};
const char	*A1b[] = {"Pyrodictium",""};

const char	*A2[ ] = {"Euryarchaeota",""};
const char	*A2a[ ] = {"Archaeoglobales",""};
const char	*A2b[ ] = {"Halobacteriales",""};
const char	*A2c[ ] = {"Methanobacteriales",""};
const char	*A2d[ ] = {"Methanococcales",""};
const char	*A2e[ ] = {"Methanomicrobiales",""};
const char	*A2f[ ] = {"Methanosarcinales","Methanosarcina","Methanosarcinaceae",""};
const char	*A2g[] = {"Thermococcales",""};
const char	*A2h[] = {"Thermoplasmales",""};

/****************************** Eubacteria **********************************/
const char	*B1[]={"Aquificales","Oxygen reducing bacteria",
			"Oxygen-reducing bacteria","Hydrogenic oxygen reducing bacteria",
			"Obligately chemolithotrophic hydrogen bacteria",
			"Aquifecales","Oligately chemolithotrophic hydrogen bacteria",""};
const char	*B2[]={"Chlamydiales",""}; 
const char	*B3[]={"Cyanobacteria","Cyanophycota","Oxyphotobacteria",
			"Oxygenic photosynthetic bacteria","blue-green algae",""};
const char	*B4[]={"Cytophagales","BCF group","CFB group","Bacteroidaceae",
			"Cytophaga","Promyxobacterium",
			"Bacteroides-Cytophaga-Flexibacter group",
			"Cytophaga-Flexibacter-Bacteroides phylum",
			"Green sulfur bacteria", ""};
const char	*B5[]={"Firmicutes","Firmacutes","gram-positive bacteria",""};
const char	*B6[]={"Proteobacteria","purple bacteria",
			"purple non-sulfur bacteria",""};
const char	*B7[]={"Spirochaetales",""};	
const char	*B8[]={"Thermotogales",""};	
const char	*B9[]={"Thermus/Deinococcus group","Deinococcaceae","Thermus group",
			"Thermus",""};	
const char	*B10[]={"Deinonema","Deinonema sp",""};	

const char	*B11[]={"Chloroflecales","Chloroflexus/Deinococcus group",
			"Chloroflexus/Deinococcaceae group",
			"Chloroflexaceae/Deinococcaceae group","Green non-sulfur bacteria",""};
const char	*B12[]={"Fibrobacter/Acidobacteria group",
			"Acidobacteria/Holophaga group","Fibrobacter group",""};
const char	*B13[]={"Flexistipes group",""};
const char	*B14[]={"Bacteroidetes","BCF group","Bacteroides-Cytophaga-Flexibacter group",
			"Cytophaga-Flexibacter-Bacteroides phylum","CFB group",""};

/*********************************************************************************/
const char	*EukaryotaSyn[] = {"Eukaryota","Eucarya","Eukarya","Eucaryotae","Eukaryotae",
			"eucaryotes","eukaryotes",""};
const char	*BacteriaSyn[5] = {"Bacteria","Eubacteria","Procaryotae","Prokaryotae",""};
const char	*OtherArchaea[4];
const char	*OtherEukaryotes[4];
const char	*OtherBacteria[4];

const char	**Eukaryota[4];
const char	**Archaea[4];
const char	**Bacteria[4]; 	

// const char	**Group[MAX_TAX_GROUPS]; 
// char		Kingdom[MAX_TAX_GROUPS]; 
// Int4	GroupNum=0;
// UInt4	GroupHits[MAX_TAX_GROUPS];

	Int4	g;

GroupNum=0;
for(g=0; g < MAX_TAX_GROUPS; g++) Group[g]=0;
g=0;
Group[g]=CopyStr(UnKnown); Kingdom[g]='U'; g++;
Group[g]=CopyStr(E1a); Kingdom[g]='P'; g++;
Group[g]=CopyStr(E1b); Kingdom[g]='P'; g++;
Group[g]=CopyStr(E1c); Kingdom[g]='P'; g++;
Group[g]=CopyStr(E1); Kingdom[g]='P'; g++;
Group[g]=CopyStr(E2); Kingdom[g]='P'; g++;
Group[g]=CopyStr(E3); Kingdom[g]='P'; g++;
Group[g]=CopyStr(E4a); Kingdom[g]='P'; g++;
Group[g]=CopyStr(E4b); Kingdom[g]='P'; g++;
Group[g]=CopyStr(E5); Kingdom[g]='P'; g++;
Group[g]=CopyStr(E6); Kingdom[g]='P'; g++;
Group[g]=CopyStr(E6a); Kingdom[g]='P'; g++;
Group[g]=CopyStr(E6b); Kingdom[g]='P'; g++;
Group[g]=CopyStr(E7); Kingdom[g]='P'; g++;
Group[g]=CopyStr(E8); Kingdom[g]='P'; g++;
Group[g]=CopyStr(E9); Kingdom[g]='P'; g++;
Group[g]=CopyStr(E10); Kingdom[g]='P'; g++;
Group[g]=CopyStr(EC1); Kingdom[g]='P'; g++;
Group[g]=CopyStr(EC2a); Kingdom[g]='P'; g++;
Group[g]=CopyStr(EC2b); Kingdom[g]='P'; g++;
Group[g]=CopyStr(EC2c); Kingdom[g]='P'; g++;
Group[g]=CopyStr(EC2d); Kingdom[g]='P'; g++;
Group[g]=CopyStr(EC3); Kingdom[g]='P'; g++;
Group[g]=CopyStr(EC4); Kingdom[g]='P'; g++;
Group[g]=CopyStr(EC5); Kingdom[g]='P'; g++;
Group[g]=CopyStr(EC6); Kingdom[g]='P'; g++;
Group[g]=CopyStr(EC7); Kingdom[g]='P'; g++;
Group[g]=CopyStr(EP1); Kingdom[g]='G'; g++;
Group[g]=CopyStr(EP2a); Kingdom[g]='G'; g++;
Group[g]=CopyStr(EP2b); Kingdom[g]='G'; g++;
Group[g]=CopyStr(EP2c); Kingdom[g]='G'; g++;
Group[g]=CopyStr(EP2d); Kingdom[g]='G'; g++;
Group[g]=CopyStr(EP2); Kingdom[g]='G'; g++;
Group[g]=CopyStr(EM0a); Kingdom[g]='M'; g++;
Group[g]=CopyStr(EM0b); Kingdom[g]='M'; g++;
Group[g]=CopyStr(EM1); Kingdom[g]='M'; g++;
Group[g]=CopyStr(EM2); Kingdom[g]='M'; g++;
Group[g]=CopyStr(EM3); Kingdom[g]='M'; g++;
Group[g]=CopyStr(EM4); Kingdom[g]='M'; g++;
Group[g]=CopyStr(EM5); Kingdom[g]='M'; g++;
Group[g]=CopyStr(EM6); Kingdom[g]='M'; g++;
Group[g]=CopyStr(EM7); Kingdom[g]='M'; g++;
Group[g]=CopyStr(EM8); Kingdom[g]='M'; g++;
Group[g]=CopyStr(EM9); Kingdom[g]='M'; g++;
Group[g]=CopyStr(EM10); Kingdom[g]='M'; g++;
Group[g]=CopyStr(EM11); Kingdom[g]='M'; g++;
Group[g]=CopyStr(EM12); Kingdom[g]='M'; g++;
Group[g]=CopyStr(EM13); Kingdom[g]='M'; g++;
Group[g]=CopyStr(EM14); Kingdom[g]='M'; g++;
Group[g]=CopyStr(EM15a); Kingdom[g]='M'; g++;
Group[g]=CopyStr(EM15b); Kingdom[g]='M'; g++;
Group[g]=CopyStr(EM16); Kingdom[g]='M'; g++;
Group[g]=CopyStr(EM17); Kingdom[g]='M'; g++;
Group[g]=CopyStr(EM18); Kingdom[g]='M'; g++;
Group[g]=CopyStr(EM19); Kingdom[g]='M'; g++;
Group[g]=CopyStr(EM20); Kingdom[g]='M'; g++;
Group[g]=CopyStr(EM21); Kingdom[g]='M'; g++;
Group[g]=CopyStr(AD1); Kingdom[g]='U'; g++;
Group[g]=CopyStr(AD2); Kingdom[g]='U'; g++;
Group[g]=CopyStr(F0a); Kingdom[g]='F'; g++;
Group[g]=CopyStr(F0b); Kingdom[g]='F'; g++;
Group[g]=CopyStr(F0c); Kingdom[g]='F'; g++;
Group[g]=CopyStr(F0); Kingdom[g]='F'; g++;
Group[g]=CopyStr(F1); Kingdom[g]='F'; g++;
Group[g]=CopyStr(F2); Kingdom[g]='F'; g++;
Group[g]=CopyStr(F3); Kingdom[g]='F'; g++;
Group[g]=CopyStr(F4); Kingdom[g]='F'; g++;
Group[g]=CopyStr(A1a); Kingdom[g]='A'; g++;
Group[g]=CopyStr(A1b); Kingdom[g]='A'; g++;
Group[g]=CopyStr(A1); Kingdom[g]='A'; g++;
Group[g]=CopyStr(A2a); Kingdom[g]='A'; g++;
Group[g]=CopyStr(A2b); Kingdom[g]='A'; g++;
Group[g]=CopyStr(A2c); Kingdom[g]='A'; g++;
Group[g]=CopyStr(A2d); Kingdom[g]='A'; g++;
Group[g]=CopyStr(A2e); Kingdom[g]='A'; g++;
Group[g]=CopyStr(A2f); Kingdom[g]='A'; g++;
Group[g]=CopyStr(A2g); Kingdom[g]='A'; g++;
Group[g]=CopyStr(A2h); Kingdom[g]='A'; g++;
Group[g]=CopyStr(A2); Kingdom[g]='A'; g++;
Group[g]=CopyStr(B1); Kingdom[g]='B'; g++;
Group[g]=CopyStr(B2); Kingdom[g]='B'; g++;
Group[g]=CopyStr(B3); Kingdom[g]='B'; g++;
Group[g]=CopyStr(B4); Kingdom[g]='B'; g++;
Group[g]=CopyStr(B5); Kingdom[g]='B'; g++;
Group[g]=CopyStr(B6); Kingdom[g]='B'; g++;
Group[g]=CopyStr(B7); Kingdom[g]='B'; g++;
Group[g]=CopyStr(B8); Kingdom[g]='B'; g++;
Group[g]=CopyStr(B9); Kingdom[g]='B'; g++;
Group[g]=CopyStr(B10); Kingdom[g]='B'; g++;
Group[g]=CopyStr(B11); Kingdom[g]='B'; g++;
Group[g]=CopyStr(B12); Kingdom[g]='B'; g++;
Group[g]=CopyStr(B13); Kingdom[g]='B'; g++;
Group[g]=CopyStr(B14); Kingdom[g]='B'; 
assert(g < MAX_TAX_GROUPS);
GroupNum=g; g++; Group[g]=0; 
// if(g!=GroupNum) { fprintf(stderr,"g = %d != %d = GroupNum\n",g,GroupNum); assert(g==GroupNum); }
	for(g=0; g <=GroupNum; g++) { GroupHits[g]=0; }
}

#endif

txs_typ::txs_typ(char *FileName,a_type A)
{ 
	init(); AB=A;
	filename=AllocString(FileName);
	if(!ReadFile( )){ print_error("txs_typ() input error"); Free(); init(); }
}

txs_typ::txs_typ(char *FileName,Int4 *GrpSize, Int4 NumGrps, 
	const char **GrpName, char *Kingdom, UInt4 NumSeqs, 
	e_type *Seq, a_type A)
{
        Int4    g;

	init(); AB=A; 
	filename=AllocString(FileName);
	total_seqs=NumSeqs;
	num_groups=NumGrps;
	NEW(group_size,num_groups+3,Int4);
	NEW(kingdom,num_groups+3,char);
	NEWP(group_name,num_groups+3,char);
        for(g=1; g <= num_groups; g++){
	     group_size[g]=GrpSize[g]; 
	     group_name[g]=AllocString(GrpName[g]);
	     kingdom[g]=Kingdom[g];
        } Load(Seq);
}

void	txs_typ::RmSubseq( )
// Remove subsequences from groups.
{
	Int4	N,i,j,g,total;
	e_type	*ListE,E1,E2;
	BooLean	swiss,pdb;
	char	Subsq;
	
	for(total=0,g=1; g <= num_groups; g++){
          N = NSeqsSeqSet(group[g]);
	  ListE=NilSeqSetRtnSeqs(group[g]); // destroys group[g].
          for(i=1; i<= N; i++){
           if(ListE[i]){
             E1 = ListE[i];
             swiss = SwissSeq(E1);
             pdb = PdbSeq(E1);
             for(j=i+1; j<= N; j++){
                if(ListE[j]){
                  E2 = ListE[j];
                  if(Subsq=IsSubSeq(E1,E2)){
                    if(Subsq==3){
                        if(pdb){ NilSeq(ListE[j]); ListE[j]=NULL; }
                        else if(PdbSeq(E2)){ NilSeq(ListE[i]); ListE[i]=NULL; }
                        else if(swiss){ NilSeq(ListE[j]); ListE[j]=NULL; }
                        else { NilSeq(ListE[i]); ListE[i]=NULL; break; }
                    } else if(Subsq==1){ NilSeq(ListE[i]); ListE[i]=NULL; }
		    else { NilSeq(ListE[j]); ListE[j]=NULL; }
                  }
                } if(!ListE[i]) break;
             }
           }
          }
          for(j=0,i=1; i<= N; i++){
		if(ListE[i]){ j++; ListE[j]=ListE[i]; }
	  } group_size[g] = j; total+= j;
	  group[g]=Array2SeqSet(ListE,group_size[g],group_name[g],AB);
	} total_seqs=total;
}

void    txs_typ::PutSeqs(FILE *fp)
{
        for(Int4 g=1; g <= num_groups; g++){
	   for(Int4 e=1; e <= NSeqsSeqSet(group[g]); e++){
		PutSeq(fp,SeqSetE(e,group[g]),AB);
	   }
	}
}

void	txs_typ::WriteFile(char *filename0)
{
	FILE	*fp;
	Int4	e,g;

	fp=open_file(filename0,".txs","w");
        fprintf(fp,"<TaxSeqs(%d;%d)={",num_groups,total_seqs);
        for(g=1; g <= num_groups; g++){
                if(g > 1) fprintf(fp,";");
                fprintf(fp,"%s(%c%d)=",group_name[g],kingdom[g],group_size[g]);
                for(e=1; e <= NSeqsSeqSet(group[g]); e++){
		   if(e > 1) fprintf(fp,",");
		fprintf(fp,"%d",LenSeq(SeqSetE(e,group[g])));
		// PutSeq(stderr,SeqSetE(e,group[g]),AB);
                }
        } fprintf(fp,"};\n\n");
        for(g=1; g <= num_groups; g++){
	   for(e=1; e <= NSeqsSeqSet(group[g]); e++){
		PutSeq(fp,SeqSetE(e,group[g]),AB);
		// PutSeq(stderr,SeqSetE(e,group[g]),AB);
	   }
        } fprintf(fp,">.\n"); fclose(fp);
}

void	txs_typ::init( )
{
	filename=0; AB=0;
	num_groups=0; group=0; group_name=0; kingdom=0;
	group_size=0; num_groups=0; total_seqs=0;
}

Int4	txs_typ::ReadFile( )
{
	char	str[200];
	Int4	NumSeqs,g,*size,s,e,sz;
	e_type	*Seq;

	FILE *fp=open_file(filename,"","r");
	if(fscanf(fp,"<TaxSeqs(%d;%d)={",&num_groups,&NumSeqs) != 2){
		return 0;
	}
	if(num_groups < 1 || NumSeqs < 1) return 0;
	total_seqs=NumSeqs;
	NEWP(group_name,num_groups+3,char);
	NEW(group_size,num_groups+3,Int4);
	NEW(kingdom,num_groups+3,char);
	NEW(Seq,NumSeqs+3,e_type);
	NEW(size,NumSeqs+3,Int4);
	for(s=0,g=1; g <=num_groups; g++){
	   if(fscanf(fp,"%[^(]",str) != 1){
		free(Seq); free(size); free(group_name); free(group_size);
		free(kingdom); return 0;
	   } else group_name[g]=AllocString(str);
	   if(fscanf(fp,"%[^=]=",str) != 1){
                free(Seq); free(size); free(group_name); free(group_size);
                free(kingdom); return 0;
           } else {
		if(sscanf(str,"(%d)",&group_size[g]) != 1){
		  if(sscanf(str,"(%c%d)",&kingdom[g],&group_size[g]) != 2){
			free(Seq); free(size); free(group_name); 
			free(group_size); free(kingdom); return 0;
		  }
		} else { kingdom[g]='U'; }
	   }
	   for(e=1; e < group_size[g]; e++){
	      if(fscanf(fp,"%d,",&sz) != 1){ free(Seq); free(size); return 0; }
	      else { s++; size[s]=sz; }
	   }
	   if(g < num_groups){
	      if(fscanf(fp,"%d;",&sz) != 1){ free(Seq); free(size); return 0; }
	      else { s++; size[s]=sz; }
	   } else if(fscanf(fp,"%d",&sz) != 1){ free(Seq); free(size); return 0; }
	   else { s++; size[s]=sz; }
	}
	fscanf(fp,"};\n\n");
	for(s=1; s <= NumSeqs; s++){ 
// fprintf(stderr,"s=%d; size=%d\n",s,size[s]);
		Seq[s] = ReadSeq(fp,s,size[s],AB); 
	}
	free(size); fscanf(fp,">."); fclose(fp);
	Load(Seq); free(Seq); 
	return NumSeqs;
}

void	txs_typ::Load(e_type *Seq)
{
	Int4	g,s,e;
	e_type	*EList;

	NEW(group,num_groups+2,ss_type);
	for(s=0,g=1; g <=num_groups; g++){
#if 0
	fprintf(stderr,"group = %d; group_size = %d (%s)\n",
		g,group_size[g],group_name[g]);
#endif
		NEW(EList,group_size[g]+3,e_type);
		for(e=1; e <= group_size[g]; e++){
		   s++; EList[e]=Seq[s];
		   EqSeqI(e,EList[e]); // this number is needed for gpsiblast.c
		   // PutSeq(stderr,EList[e],AB);
		}
		// WARNING: EList is absorbed into the new object.
		assert(group_size[g] > 0);
		group[g]=Array2SeqSet(EList,group_size[g],group_name[g],AB);
	} assert(s==total_seqs);
}

void	txs_typ::Free( )
{
	Int4 g;
	for(g=1; g <= num_groups; g++){
	    if(group_name && group_name[g]) free(group_name[g]);
	    if(group) NilSeqSet(group[g]);
	}
	if(group_name) free(group_name);
	if(group_size) free(group_size);
	if(kingdom) free(kingdom);
	if(group) free(group);
	if(filename) free(filename);
}

