# makefile for libafn

CC = ${CURRENT_CC}
CFLAGS = $(LIBDIRS) $(INCDIRS) -c

INCDIRS = -I. -I../include 

STDINC = ../include/stdinc.h ../include/sset.h
BIOINC = ../include/residues.h ../include/blosum62.h

LIBS = ${MATH_LIBRARY} 
LIBDIRS = -L../libafn -L../libstat


OBJS = afnio.o alphabet.o binomial.o block.o dheap.o digraph.o dms_typ.o dsets.o firm.o \
	guide.o haussler.o histogram.o iarray.o karlin.o list.o mheap.o mlist.o \
	nrutil.o olist.o order.o probability.o pseg.o random.o randwalk.o swaln.o \
	segment.o seqset.o sequence.o set_typ.o sqheap.o stack.o sumstat.o table.o \
	wdigraph.o idp_typ.o cls_typ.o rff_typ.o license.o bstree.o bpps_typ.o \
	bpps_quick.o pps_typ.o rst_typ.o pat_typ.o ptrn_typ.o fdt_typ.o

INC = afnio.h alphabet.h binomial.h block.h dheap.h digraph.h dirichlet.h dms_typ.h dsets.h firm.h \
	guide.h haussler.h histogram.h iarray.h karlin.h list.h lnfac.h mheap.h mlist.h \
	nrutil.h olist.h order.h probability.h pseg.h random.h randwalk.h swaln.h \
	segment.h seqset.h sequence.h set_typ.h sqheap.h stack.h sumstat.h table.h \
	wdigraph.h idp_typ.h cls_typ.h rff_typ.h license.h bstree.h pps_typ.h \
	bpps_typ.h rst_typ.h pat_typ.h ptrn_typ.h pch_typ.h fdt_typ.h $(STDINC)

libafn.a: $(OBJS) $(INC)
	ar rcs libafn.a $(OBJS)
#	ranlib libafn.a

.cc.o :
	${CC} $(CFLAGS) $*.cc

iblock.s: iblock.cc block.h $(STDINC)
	${CC} -O -S $(LIBDIRS) ${INCDIRS} iblock.cc 

lgamma: gamma.h maingamma.o probability.o $(STDINC)
	${CC} -O ${INCDIRS} -o lgamma gamma.o probability.o maingamma.o ${MATH_LIBRARY}

maingamma.o: gamma.o gamma.h $(STDINC)

driver_seqset: $(INC) driver_seqset.cc libafn.a ${STDINC} ${BIOINC}
	${CC} $(LIBDIRS) $(INCDIRS) -o driver_seqset driver_seqset.cc -lafn -lstat ${LIBS}

#driver_seqset: $(INC) driver_seqset.o libafn.a ${STDINC} ${BIOINC}
#	${CC} $(LIBDIRS) $(INCDIRS) -o driver_seqset driver_seqset.o -lafn -lstat ${LIBS}

#basic algorithmic routines
afnio.o: afnio.h $(STDINC)
digraph.o: digraph.h list.h random.h $(STDINC)
wdigraph.o: wdigraph.h list.h random.h $(STDINC)
dheap.o: dheap.h $(STDINC)
mheap.o: mheap.h dheap.h $(STDINC)
dsets.o: dsets.h $(STDINC)
set_typ.o: set_typ.h $(STDINC)
fdt_typ.o: fdt_typ.h $(STDINC)
block.o: block.h $(STDINC)
iarray.o: iarray.h $(STDINC)
list.o: list.h $(STDINC)
mlist.o: mlist.h $(STDINC)
olist.o: olist.h $(STDINC)
order.o: order.h $(STDINC)
stack.o: stack.h $(STDINC)
bstree.o: bstree.h $(STDINC)

# statistical routines
nrutil.o: nrutil.h
license.o: license.h $(STDINC)
random.o: random.h probability.h $(STDINC)
betaprior.o: betaprior.h afnio.h probability.h $(STDINC)
binomial.o: binomial.h random.h afnio.h histogram.h $(STDINC)
firm.o: firm.h afnio.h $(STDINC)
gamma.o: gamma.h afnio.h $(STDINC)
randwalk.o: randwalk.h afnio.h $(STDINC)
haussler.o: haussler.h alphabet.h probability.h dirichlet.h $(INC)
histogram.o: histogram.h $(STDINC)
karlin.o: karlin.h alphabet.h $(STDINC)
sumstat.o: sumstat.h
probability.o: probability.h $(STDINC)
table.o: table.h alphabet.h probability.h $(STDINC)

dms_typ.o: dms_typ.h dms_data.h $(INC) $(BIOINC)
pps_typ.o: pps_typ.h $(INC)
bpps_typ.o: bpps_typ.h $(INC)
bpps_quick.o: bpps_typ.h $(INC)
rst_typ.o: rst_typ.h $(INC)

# sequence routines:
alphabet.o: alphabet.h $(STDINC)
sequence.o: sequence.h afnio.h alphabet.h random.h dheap.h rff_typ.h $(STDINC)
seqset.o: seqset.h afnio.h sequence.h karlin.h rff_typ.h $(INC)
pseg.o: pseg.h afnio.h sequence.h alphabet.h rff_typ.h $(STDINC)
guide.o: guide.h sequence.h alphabet.h random.h dheap.h rff_typ.h $(STDINC)
sqheap.o: sqheap.h sequence.h histogram.h mheap.h ${STDINC}
segment.o: segment.h $(STDINC)
idp_typ.o: idp_typ.h $(STDINC)
cls_typ.o: cls_typ.h alphabet.h sequence.h rff_typ.h $(STDINC)
rff_typ.o: rff_typ.h $(STDINC)
pat_typ.o: pat_typ.h $(INC)
ptrn_typ.o: ptrn_typ.h $(INC)
swaln.o: swaln.h $(BIOINC) $(INC)

clean: 
	@ \rm *.o *.a

