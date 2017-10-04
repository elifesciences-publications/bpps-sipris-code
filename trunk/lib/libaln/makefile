# makefile for libaln

CC = ${CURRENT_CC}
CFLAGS = $(LIBDIRS) $(INCDIRS) -c

ARCHIVE = ../libafn/libafn.a ../libgsq/libgsq.a
INCDIRS = -I. -I../include -I../libafn -I../libgsq 

LIBS = -lstat -lblst -lafn -lgsq -lstdc++ ${MATH_LIBRARY}
LIBDIRS = -L../libafn/ -L../libgsq/

STDINC = ../include/stdinc.h ../include/sset.h
BIOINC = ../include/residues.h ../include/blosum62.h

OBJS = align.o asset.o cmsa.o cmsa_map.o cmsa_recombine.o colwt.o \
	dfps.o eblocks.o evalue.o fmodel.o fsm.o genetic.o \
	gibbs.o msaheap.o oscan.o \
	pairaln.o pattern.o pheap.o profile.o prtn_model.o \
	scanheap.o scaninfo.o sites.o sma.o smatrix.o psialn.o \
	smooth.o strctprof.o watchgibbs.o wmodel.o \
	align_smx.o sample_smx.o cmsa_put.o cmsa_operations.o \
	cmsa_rtf.o cmsa_sma.o cmsa_sample.o cmsa_io.o sseq.o msites.o \
	swt_typ.o hsw_typ.o mem_typ.o cmsa_emit_sq.o cmsa_recombine_gaps.o \
	mdl.o dom_typ.o cmafull.o pah_typ.o

INC = align.h asset.h cmsa.h colwt.h \
	eblocks.h evalue.h fmodel.h fsm.h genetic.h \
	gibbs.h msaheap.h oscan.h \
	pairaln.h pattern.h pheap.h profile.h prtn_model.h \
	scanheap.h scaninfo.h sites.h sma.h psialn.h \
	smatrix.h smooth.h strctprof.h \
	watchgibbs.h wmodel.h sseq.h msites.h mem_typ.h \
	swt_typ.h hsw_typ.h mdl.h dom_typ.h pah_typ.h $(STDINC) $(ARCHIVE)

libaln.a: $(OBJS) $(INC)
	ar rcs libaln.a $(OBJS) 
#	ranlib libaln.a

.cc.o:
	${CC} $(CFLAGS) $*.cc 
#	${CC} $(CFLAGS) $*.cc ${LIBS}

eblocks.s: eblocks.cc $(INC) $(OBJS) 
	${CC} -O -S eblocks.cc ${LIBS} $(INCDIRS)

junk.s: junk.cc $(INC) $(OBJS)
	${CC} -O -S junk.cc ${LIBS} $(INCDIRS)

# psimatrix needs ncbi blast library when compiled.
# psimatrix.o: psimatrix.h $(STDINC) $(NCBI_INCDIR)/*.h

psialn.o: psialn.h $(STDINC) 

# asset routines
asset.o: $(INC)
dfps.o: $(INC)
pattern.o: $(INC)
pheap.o: $(INC)
evalue.o: $(INC)
align.o: msites.h $(INC)
eblocks.o: $(INC)
profile.o: $(INC)

# pairwise alignment routines:
fsm.o: fsm.h $(INC) ${STDINC}
pairaln.o: pairaln.h $(INC)

# multiple alignment routines:
sites.o: sites.h $(INC)
msites.o: msites.h $(INC)
cmsa_put.o: cmsa.h $(INC)
cmsa_rtf.o: cmsa.h $(INC)
cmsa_sample.o: cmsa.h $(INC)
cmafull.o: cmsa.h $(INC)
cmsa_sma.o: cmsa.h $(INC)
cmsa_io.o: cmsa.h $(INC)
cmsa_operations.o: cmsa.h $(INC)
cmsa.o: cmsa.h $(INC) 
cmsa_recombine.o: cmsa.h $(INC) 
cmsa_emit_sq.o: cmsa.h $(INC) 
cmsa_recombine_gaps.o: cmsa.h $(INC) 
cmsa_map.o: cmsa.h $(INC) 
msaheap.o: cmsa.h fmodel.h sites.h colwt.h $(BIOINC) $(INC) 
colwt.o: colwt.h $(ARCHIVE) $(STDINC)
sma.o: sma.h $(INC)
fmodel.o: fmodel.h $(INC)
mdl.o: mdl.h $(BIOINC) $(INC)
swt_typ.o: swt_typ.h $(INC) $(ARCHIVE)
hsw_typ.o: hsw_typ.h $(INC) $(ARCHIVE)

# scan routines:
smatrix.o: smatrix.h $(ARCHIVE) $(STDINC) $(BIOINC)
align_smx.o: smatrix.h $(ARCHIVE) $(STDINC) $(BIOINC)
sample_smx.o: smatrix.h $(ARCHIVE) $(STDINC) $(BIOINC)
wmodel.o: wmodel.h smatrix.h $(BIOINC) $(INC)
oscan.o: oscan.h wmodel.h sites.h scanheap.h scaninfo.h $(BIOINC) $(INC)
prtn_model.o: prtn_model.h wmodel.h $(BIOINC) $(INC)
scaninfo.o: scaninfo.h $(ARCHIVE) $(BIOINC) $(STDINC)
scanheap.o: scanheap.h scaninfo.h $(ARCHIVE) $(BIOINC) $(STDINC)
smooth.o: smooth.h $(BIOINC) $(INC) $(SPOUGE)

# high level routines
watchgibbs.o: watchgibbs.h $(INC)
dom_typ.o: dom_typ.h $(STDINC) $(ARCHIVE)
pah_typ.o: pah_typ.h $(INC)

# gibbs routines:
gibbs.o: gibbs.h $(INC) 
genetic.o: genetic.h gibbs.h msaheap.h $(ARCHIVE) $(STDINC)
sseq.o: sseq.h $(ARCHIVE) $(INC)
mem_typ.o: mem_typ.h $(INC)

# structural routines:
strctprof.o: strctprof.h $(INC)

clean: 
	@ \rm *.o *.a

