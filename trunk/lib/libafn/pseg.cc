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

/************************************************************************/
/***  (pseg.c)                                                        ***/
/*** #include precomputed ln(fact(n)) array in "lnfac.h"              ***/
/************************************************************************/
#include "pseg.h"
#include "lnfac.h"

static int window = 12;
static int downset, upset;
static double locut = 2.2;
static double hicut = 2.5;
static int period = 0;
static int hilenmin = 0;
static int overlaps = FALSE;
static int maxtrim = 100;

#define LN2	0.69314718055994530941723212145818

static int      thewindow;
static double	*entray=NULL;
static int		aaindex[128];
static unsigned char	aaflag[128];
static char		aachar[ALPHA_SIZE];

/******************************* Public *********************************/
BooLean	ProcessSeqPSeg(int win, double lowcut, double highcut, 
	int mxtrim, e_type E, a_type A)
{
   static char		initialize_pseg=TRUE;
   struct Sequence	*seq;
   struct Segment	*segs;
   BooLean		mask;

   if(initialize_pseg) { initialize_pseg=FALSE; genwininit(); }
   window = win; locut = lowcut; hicut = highcut;
   if(locut > hicut) {
      fprintf(stderr, "Warning: trigger complexity (%.2f) greater than extension (%.2f)\n",locut,hicut);
      hicut = locut;
   }
   downset = (window+1)/2 - 1;
   upset = window - downset;
   maxtrim = mxtrim;
   entropy_init(window);
   seq=CreatePSegSeq(E,A);
   segs = (struct Segment *) NULL;
   segseq(seq, &segs, 0);
   mergesegs(seq, segs);
   // singreport(seq, segs); /** DEBUG **/
   mask=process_pseg(segs, E);
   freesegs(segs); closeseq(seq);
   if(entray!=NULL) free(entray);
   return mask;
}

/******************************* private ********************************/

struct Sequence *seq_phase(struct Sequence *seq, int phase, int Period)
{
   struct Sequence *perseq;
   int len, i, j;

   perseq = (struct Sequence *) malloc(sizeof(struct Sequence));
   len = ((seq->length)/Period) + 1;
   perseq->seq = (char *) malloc((len+1)*sizeof(char));
   perseq->length = 0;
   for (i=0, j=phase; j<seq->length; i++, j+=Period)
     { perseq->seq[i] = seq->seq[j]; perseq->length++; }
   perseq->seq[i] = '\0';
   return(perseq);
}

void    segseq(struct Sequence *seq, struct Segment **segs,int  offset)
{
   struct Segment *seg, *leftsegs;
   struct Sequence *leftseq;
   int		first, last, lowlim;
   int		loi, hii, i;
   int		leftend, rightend, lend, rend;
   double 	*H;

   H = seqent(seq);
   if (H==NULL) return;
   first = downset;
   last = seq->length - upset;
   lowlim = first;
   for (i=first; i<=last; i++) {
      if (H[i]<=locut && H[i]!=-1) {
         loi = findlo(i, lowlim, H);
         hii = findhi(i, last, H);
         leftend = loi - downset;
         rightend = hii + upset - 1;
         trim(openwin(seq, leftend, rightend-leftend+1), &leftend, &rightend);
         if(i+upset-1<leftend){   /* check for trigger window in left trim */
            lend = loi - downset;
            rend = leftend - 1;
            leftseq = openwin(seq, lend, rend-lend+1);
            leftsegs = (struct Segment *) NULL;
            segseq(leftseq, &leftsegs, offset+lend);
            if (leftsegs!=NULL) {
               if (*segs==NULL) *segs = leftsegs;
               else appendseg(*segs, leftsegs);
            }
            closewin(leftseq);
         }
         seg = (struct Segment *) malloc(sizeof(struct Segment));
         seg->begin = leftend + offset;
         seg->end = rightend + offset;
         seg->next = (struct Segment *) NULL;
         if (*segs==NULL) *segs = seg;
         else appendseg(*segs, seg);
         i = MINIMUM(int,hii, rightend+downset);
         lowlim = i + 1;
/***       i = hii;     this ignores the trimmed residues... ***/
      }
   }
   free(H);
}

double *seqent(struct Sequence *seq)
{
   struct Sequence *win;
   double *H;
   int i, first, last;

   if (window>seq->length) { return((double *) NULL); }
   H = (double *) malloc(seq->length*sizeof(double));
   for (i=0; i<seq->length; i++) { H[i] = -1.; }
   win = openwin(seq, 0, window); enton(win);
   first = downset; last = seq->length - upset;
   for(i=first; i<=last; i++) { H[i]=win->entropy; shiftwin1(win); }
   closewin(win);
   return(H);
}

int    findlo(int i, int limit, double *H)
{
   int j;

   for(j=i; j>=limit; j--) { if(H[j]==-1) break; if(H[j]>hicut) break; }
   return(j+1);
}

int    findhi(int i, int limit, double *H)
{
   int j;

   for(j=i; j<=limit; j++) { if(H[j]==-1) break; if(H[j]>hicut) break; }
   return(j-1);
}


void    trim(struct Sequence *seq, int *leftend, int *rightend)
{
   struct Sequence *win;
   double prob, minprob;
// double	realprob;
   int		 shift, len, i;
   int		lend, rend,minlen;
   double 	totseq;


// fprintf(stderr, "%d %d\n", *leftend, *rightend); 

   lend = 0;
   rend = seq->length - 1;
   minlen = 1;
   if((seq->length-maxtrim)>minlen) minlen = seq->length-maxtrim;
   minprob = 1.;
//   fprintf(stderr, "\n");               
   for (len=seq->length; len>minlen; len--) {
//      fprintf(stderr, "%5d ", len);    
      win=openwin(seq, 0, len);

      for(i=0, shift=TRUE; shift; i++) {
         /** was: prob = getprob(win->state, len); - afn ***/
   	 totseq = ((double) len) * 2.9957322735539909; /** = LN20 **/
	 prob = lnass(win->state) + lnperm(win->state, len) - totseq;

//         realprob = exp(prob);             /*-(for tracing the trim)-*/
//         fprintf(stderr, "%2e ", realprob);                       /*-*/
         if (prob<minprob)
           { minprob = prob; lend = i; rend = len + i - 1; }
         shift = shiftwin1(win);
      }
      closewin(win);
//      fprintf(stderr, "\n");                                     /*-*/
     }
//   fprintf(stderr, "%d-%d ", *leftend, *rightend);               /*-*/
   *leftend = *leftend + lend;
   *rightend = *rightend - (seq->length - rend - 1);

//   fprintf(stderr, "%d-%d\n", *leftend, *rightend);              /*-*/

   closewin(seq);
}

double lnperm(int *sv, int tot)
{
   double ans;
   int i;

   for(ans=lnfac[tot], i=0; sv[i]!=0; i++) { ans-=lnfac[sv[i]]; }
   return(ans);
}

double lnass(register int *sv)
{
	double	ans;
	register int	svi, svim1;
	register int	Class, total;
	register int    i;

	ans=lnfac[ALPHA_SIZE];
	if(sv[0]==0) return ans;
	total=ALPHA_SIZE; Class=1; svim1=sv[0];
	for(i=0;; svim1=svi) {
	        if (++i==ALPHA_SIZE) { ans -= lnfac[Class]; break; }
		else if ((svi = *++sv) == svim1) { Class++; continue; }
		else {
			total -= Class;
			ans -= lnfac[Class];
			if (svi == 0) { ans -= lnfac[total]; break; }
			else { Class = 1; continue; }
		}
	}
	return ans;
}

void    mergesegs(struct Sequence *seq, struct Segment *segs)
{
   struct Segment *seg, *nextseg;
   int len;

   if(overlaps) return;
   if(segs==NULL) return;
   if(segs->begin<hilenmin) segs->begin = 0;
   seg = segs;
   nextseg = seg->next;
   while (nextseg!=NULL) {
      if(seg->end>=nextseg->begin && seg->end>=nextseg->end) {
         seg->next = nextseg->next; free(nextseg);
         nextseg = seg->next; continue;
      }
      if(seg->end>=nextseg->begin){              /* overlapping segments */
         seg->end = nextseg->end; seg->next = nextseg->next;
         free(nextseg); nextseg = seg->next; continue;
      }
      len = nextseg->begin - seg->end - 1;
      if(len<hilenmin){                           /* short hient segment */
         seg->end = nextseg->end; seg->next = nextseg->next;
         free(nextseg); nextseg = seg->next; continue;
      }
      seg = nextseg; nextseg = seg->next;
   }
   len = seq->length - seg->end - 1;
   if (len<hilenmin) seg->end = seq->length - 1;
}

struct Segment *per_mergesegs(struct Sequence *seq, struct Segment **persegs)
{
   struct Segment **localsegs;
   struct Segment *firstseg, *segs, *seg;
   int first;
   int phase, savephase;

   segs = (struct Segment *) NULL;
   if (persegs==NULL) return(segs);

   localsegs = (struct Segment **) malloc(period*sizeof(double));
   for (phase=0; phase<period; phase++) localsegs[phase] = persegs[phase];

   while (1) {
      firstseg = (struct Segment *) NULL;
      first = -1;
      savephase = -1;
      for (phase=0; phase<period; phase++) {
         if (localsegs[phase]==NULL) continue;
         if (first==-1 || localsegs[phase]->begin<first) {
            savephase = phase;
            firstseg = localsegs[phase];
            first = firstseg->begin;
         }
      }
      if (firstseg==NULL) break;
      seg = (struct Segment *) malloc(sizeof(struct Segment));
      seg->begin = ((firstseg->begin)*period)+savephase;
      seg->end = ((firstseg->end)*period)+savephase;
      seg->next = (struct Segment *) NULL;
      if (segs==NULL) segs = seg; else appendseg(segs, seg);
      localsegs[savephase] = localsegs[savephase]->next;
   }
   free(localsegs); mergesegs(seq, segs);
   return(segs);
}

void	per_seqprep(struct Sequence *seq, struct Segment **persegs)
{
   char		*proseq;
   struct Segment *segs, *seg;
   int		begin, end, i, phase, pos;

   proseq = seq->seq;
   upper(proseq, seq->length);
   for(phase=0; phase<period; phase++) {
      segs = persegs[phase];
      for(seg=segs; seg!=NULL; seg=seg->next) {
         begin = seg->begin; end = seg->end; 
         for(i=begin; i<=end; i++) {
            pos=phase+period*i;
            if(isalpha(proseq[pos])) proseq[pos]=tolower(proseq[pos]);
         }
      }
   }
   for(i=0; i<seq->length; i++) if(islower(proseq[i])) proseq[i]='x';
}

BooLean	process_pseg(struct Segment *segs, e_type E)
{
	unsigned char	*proseq;
	struct Segment	*seg;
	int	begin, end;
	BooLean	mask=FALSE;

	AddXArraySeq(E);
	proseq = XSeqPtr(E); proseq++;	/** E array starts with 1 not 0 **/
	for (seg=segs; seg!=NULL; seg=seg->next) {
		begin = seg->begin;
		end = seg->end;
		memset(proseq + begin, 0, (end - begin +1));
		mask = TRUE;
	}
	return mask;
}

void	per_process_pseg(struct Sequence *seq, e_type E)
{
   char *proseq, *end;
   int		i;

   AddXArraySeq(E);
   proseq = seq->seq; end = proseq +seq->length;
   for(i=1; proseq<end; ++i, ++proseq) {if(*proseq=='x') EqXSeq(i,0,E);}
}

void    singreport(struct Sequence *seq, struct Segment *segs)
{
	char	*proseq, *proseqmax;
	struct Segment	*seg;
	int	begin, end, i, ctr;

	proseq = seq->seq;
	proseqmax = proseq + seq->length;
	upper(proseq, seq->length);
	for (seg=segs; seg!=NULL; seg=seg->next) {
		begin = seg->begin;
		end = seg->end;
		memset(proseq + begin, 'x', end - begin +1);
	}
	fprintf(stdout, "%s\n", seq->header);
	for (i=0, ctr=0; proseq < proseqmax; ++i, ++ctr, ++proseq) {
		if (ctr==60) {
			putc('\n', stdout);
			ctr = 0;
		}
		putc(*proseq, stdout);
	}
	putc('\n', stdout);
	if (putc('\n', stdout) == EOF) {
		fprintf(stderr, "premature EOF on write\n"); exit(2);
	}
}

#if 0
void	per_singreport(struct Sequence *seq, struct Segment **persegs,
  	struct Segment *psegs)
{
   char *proseq, *proseqmax;
   int i, ctr;

   proseq = seq->seq;
   proseqmax = proseq +seq->length;
   fprintf(stdout, "%s\n", seq->header);

   for (i=0, ctr=0; proseq < proseqmax; ++i, ++ctr, ++proseq) {
      if (ctr==60) { putc('\n', stdout); ctr = 0; }
      putc(*proseq, stdout);
   }
   putc('\n', stdout);
   if (putc('\n', stdout) == EOF) {
      fprintf(stderr, "premature EOF on write\n"); exit(2);
   }
}
#endif

void    appendseg(struct Segment *segs, struct Segment *seg)
{
   struct Segment *temp;

   temp = segs;
   while (1) {
      if (temp->next==NULL) { temp->next = seg; break; }
      else { temp = temp->next; }
   }
}

void    freesegs(struct Segment *segs)
{
   struct Segment *temp;

   while (segs!=NULL) { temp = segs->next; free(segs); segs = temp; }
}

/******************************* pgenwin.c ****************************/
struct Sequence *CreatePSegSeq(e_type E, a_type A)
{
   struct Sequence *seq;
   int		r,i,j;

   seq = (struct Sequence *) malloc(sizeof(struct Sequence));
/*---                       ---[backpointers null at the top]---*/
   seq->parent = (struct Sequence *) NULL;
   seq->root = (struct Sequence *) NULL;
   if(SeqKey(E) != NULL) seq->header = AllocString(SeqKey(E));
   else seq->header = AllocString("random");
   seq->length = LenSeq(E);
   MEW(seq->seq, seq->length+2,char);
   for(i=0,j=1; i < seq->length; i++,j++){
	r = ResSeq(j,E);
	seq->seq[i]= AlphaChar(r,A);
   }
/*--- [set implementation parameters to be initially unconfigured]---*/
   seq->entropy = -2.;
   seq->state = (int *) NULL;
   seq->composition = (int *) NULL;
   return(seq);
}

void	genwininit(void)
{
	char	*cp, *cp0;
	int		i;
	char	c;

	for(i=0; i<sizeof(aaindex)/sizeof(aaindex[0]); ++i){
		aaindex[i] = ALPHA_SIZE;
		aaflag[i] = TRUE;
	}
	for (cp = cp0 = "ACDEFGHIKLMNPQRSTVWY"; (c = *cp) != '\0'; ++cp) {
		i = cp - cp0;
		aaindex[c] = i;
		aaindex[tolower(c)] = i;
		aachar[i] = tolower(c);
		aaflag[c] = FALSE;
		aaflag[tolower(c)] = FALSE;
	}
}
        
void	closeseq(struct Sequence *seq)
{
   if (seq==NULL) return;
   if (seq->header!=NULL)      free(seq->header);
   if (seq->state!=NULL)       free(seq->state);
   if (seq->composition!=NULL) free(seq->composition);
   free(seq->seq); free(seq);
}

struct Sequence *openwin(struct Sequence *parent, int start,
	int length)
{
   struct Sequence *win;

   if (start<0 || length<0 || start+length>parent->length) {
      return((struct Sequence *) NULL);
   }
   win = (struct Sequence *) malloc(sizeof(struct Sequence));
/*---                                   ---[set links, up and down]---*/
   win->parent = parent;
   if (parent->root==NULL) {win->root = parent;}
   else {win->root = parent->root;}
   win->header = (char *) NULL;

/*---                  ---[install the local copy of the sequence]---*/
   win->start = start;
   win->length = length;
   win->seq = parent->seq + start;
/*---     [set local flags to be initially unconfiguerd window]---*/
	win->entropy = -2.;
	win->state = (int *) NULL;
	win->composition = (int *) NULL;
	stateon(win);
	return win;
}

BooLean	shiftwin1(struct Sequence *win)
{
	register int	j, length;
	register int	*comp;

	length = win->length;
	comp = win->composition;
	if ((++win->start + length) > win->parent->length) {
		--win->start; return FALSE;
	}
	if (!aaflag[j = win->seq[0]])
		decrementsv(win->state, comp[aaindex[j]]--);
	j = win->seq[length];
	++win->seq;
	if (!aaflag[j]) incrementsv(win->state, comp[aaindex[j]]++);
	if (win->entropy > -2.)
		win->entropy = entropy(win->state);
	return TRUE;
}

void	closewin(struct Sequence *win)
{
   if (win==NULL) return;
   if (win->state!=NULL)       free(win->state);
   if (win->composition!=NULL) free(win->composition);
   free(win);
}

void	compon(struct Sequence *win)
{
	register int	*comp;
	register int	aa;
	register char	*seq, *seqmax;

	win->composition = comp = (int *) calloc(ALPHA_SIZE*sizeof(*comp), 1);
	seq = win->seq;
	seqmax = seq + win->length;
	while(seq<seqmax) { aa=*seq++; if(!aaflag[aa]) comp[aaindex[aa]]++; }
}

static int state_cmp(int *s1,int *s2) { return (*s2-*s1); }

void	stateon(struct Sequence	*win)
{
	register int	aa, nel, c;

	if (win->composition == NULL) compon(win);
	win->state = (int *) malloc(ALPHA_SIZE_PLUS*sizeof(win->state[0]));
	for (aa = nel = 0; aa < ALPHA_SIZE; ++aa) {
		if ((c = win->composition[aa]) == 0) continue;
		win->state[nel++] = c;
	}
	for (aa = nel; aa < ALPHA_SIZE_PLUS; ++aa)
		win->state[aa] = 0;
	qsort(win->state, nel, sizeof(win->state[0]), 
		(int (* )(const void *  , const void * ))state_cmp);
}

void	enton(struct Sequence *win)
{
   if (win->state==NULL) {stateon(win);}
   win->entropy = entropy(win->state);
}

void entropy_init(int Window)
{
	int	i;
	double	x, xw;

	entray = (double *)malloc((Window+1) * sizeof(*entray));
	xw = Window;
	for (i = 1; i <= Window; ++i) {
		x = i / xw; entray[i] = -x * log(x) / LN2;
	}
	thewindow = Window;
}

double	entropy(register int *sv)
{
	int	*sv0 = sv;
	register double	ent;
	register int	i, total;
	register int	*svmax;
	register double	xtotrecip, xsv;

	for (total = 0; (i = *sv) != 0; ++sv) total += i;
	svmax = sv;
	ent = 0.0;
	if (total == thewindow) {
		for (sv = sv0; sv < svmax; ) { ent += entray[*sv++]; }
		return ent;
	}
	if (total == 0) return 0.;
	xtotrecip = 1./(double)total;
	for (sv = sv0; sv < svmax; ) {
		xsv = *sv++;
		ent += xsv * log(xsv * xtotrecip);
	}
	return -ent * xtotrecip / LN2;
}

void decrementsv(register int *sv, register int Class)
{
	register int	svi;

	while ((svi = *sv++) != 0) {
		if(svi==Class && *sv<Class) {sv[-1]=svi-1; break; }
	}
}

void incrementsv(register int *sv, register int Class)
{ for (;;) { if (*sv++ == Class) { sv[-1]++; break; } } }

void	upper(register char *string, size_t len)
{
	register char	*stringmax, c;

	for (stringmax = string + len; string < stringmax; ++string)
		if (islower(c = *string)) *string = toupper(c);
}

