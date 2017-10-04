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

/******************************** mdl.h - *******************************
   		   Protein Sequence Model Data Type:
 ************************************************************************/
#if !defined(_MDL_)
#define _MDL_
#include "afnio.h"
#include "fmodel.h"
#include "sites.h"
#include "alphabet.h"
#include "random.h"
#include "gss_typ.h"
#include "spouge.h"

class mdl_typ {
public:
                mdl_typ(Int4,Int4 *,Int4,Int4 *,double,BooLean **,a_type);
                mdl_typ(const mdl_typ& mdl);
                ~mdl_typ( ){ Free(); }
        mdl_typ& operator=(const mdl_typ&);     // assignment operator.
	Int4	AddColumn(Int4 t,Int4 *obs,Int4 pos);
	void	Add(Int4 t, unsigned char *seq, Int4 site);
	void	Remove(Int4 t, unsigned char *seq, Int4 site);
	Int4	MvColumn(Int4 t, Int4 *obs, Int4 lemon, Int4 pos);
	Int4	RmColumn(Int4 t, Int4 pos);
	void    Shift(Int4 t,Int4 *obs, BooLean left);
	void    SetPseudo(Int4,double);
	double  LogProbRatio(Int4,Int4 *,double);
	fm_type	Model(Int4 t) { return model[t]; }
	fm_type	*Models( ) { return model; }
	BooLean Enlarge(Int4 t,Int4 maxlen){return EnlargeFModel(maxlen,model[t]);}
	double  Ratio(Int4 t,Int4 *obs,Int4 d){return RatioFModel(obs,d,model[t]);}
	Int4    Observed(Int4 t,Int4 **array){return ObservedFModel(array,model[t]);}

	smx_typ	*SampleSMX(double pernats,Int4 wt);
	char	*GapAlnTrace(Int4,Int4,Int4,unsigned char *,Int4 **,Int4 *);
	char    *gapped_aln_trace(Int4,Int4,Int4,unsigned char *,Int4 **,
		Int4 *,double *,Int4 *);
#if 0
	BooLean	RemoveColumn(Int4 t, Int4 pos) // PRIVATE TO FMODEL???
		{return RmColumnFModel(pos, model[t]);} 

	double	Likelihood(Int4 t,unsigned char*sq,Int4 p)
		{ return LikelihoodFModel(sq,p,model[t]); }
	double	NormLikelihood(Int4 t){ return NormLikelihoodFModel(model[t]); }
	Int4	ChoiceLemon(Int4 t){ return ChoiceLemonFModel(model[t]); }
	Int4	ChoiceOrange(Int4 t){ return ChoiceOrangeFModel(model[t]); }
        void    Put(FILE *fp, Int4 t){ PutFModel(fp,model[t]); }
	Int4	Lemon(Int4 t){ return LemonFModel(model[t]); }
	BooLean NullSite(Int4 t, Int4 s){ return NullSiteFModel(s,model[t]); }
	Int4	Orange(Int4 t){ return OrangeFModel(model[t]); }
	double	Prob(Int4 t, register unsigned char *seq, register Int4 pos, 
		register double p) { return ProbFModel(seq, pos, p, model[t]); }
	smx_typ GetSmatrix(Int4 t, double pernats)
		{ return GetSmatrixFModel(pernats, model[t]); }
	Int4	Length(Int4 t){ return LenFModel(model[t]); }
	BooLean	Contig(Int4 t){ return ContigFModel(model[t]); }
	double	*Pseudo(Int4 t){ return PseudoFModel(model[t]); }
	Int4	nColumns(Int4 t){ return nColsFModel(model[t]); }
	Int4	TotSites(Int4 t){ return TotSitesFModel(model[t]); }
	Int4	MaxLen(Int4 t){ return MaxLenFModel(model[t]); }
        double  DiffLogProbRatio(Int4 , Int4 *, double , Int4 , Int4 );
        double  DiffLogProbRatio(Int4 , Int4 *, double , Int4 *, Int4 , Int4 );
        double  DiffLogProbRatio(Int4,double,double,unsigned char **,unsigned char **);
#endif
	Int4	tNumModels( ){ return NumModels; }
	fm_type	*tmodel( ){ return model; }
	a_type	tAB( ){ return AB; }
	UInt4	*tobservedNS( ){ return observedNS; }
	double	*tfreq() { return freq; }
	UInt4	*tcounts() { return counts; }
private:
	// **** functions: ****
	void    	Free( );
	void    	Check( );
        void    	copy(const mdl_typ&);
	// **** pointers: ****
	// st_type	sites;		// sites pointer...
	a_type		AB;		// alphabet.
	// **** objects: ****
        Int4		NumModels;	// # product multinomial models
        fm_type		*model;		// product multinomial models
	// **** Nonsite models: ****
	UInt4	*observedNS;	// observed counts for nonsites.
	UInt4	*counts;
	double		*freq;		// background frequencies.
	double		pseudo;		// pseudo counts.
	// **** Future additions: ****
	// UInt4	*observedInDel;	// observed indel counts.
	// double	*likelihoodNS;	// nonsite likelihood
	// BooLean	*update;	// update specific models?
	// Int4	**gaps;		// gaps[m][g]: # observed gaps between blocks
	// Int4	**del;		// del[m][s]: deletion at each position.
	// Int4	**ins;		// ins[m][s]: insertions at each position.
        // BooLean      **null;		// column configurations
	// Int4		*counts;
	// Int4		*SeqLens;
	// Int4		NumSeq;
	// gapfunct;
	// unsigned char *wt;	// sequence weight 0..100.
	// blsm62;
};

#endif

