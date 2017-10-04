# makefile for libtools

CC = ${CURRENT_CC}
CFLAGS = $(LIBDIRS) $(INCDIRS) -c

ARCHIVE = ../libafn/libafn.a ../libblst/libblst.a ../libhmm/libhmm.a \
	../libgsq/libgsq.a ../libaln/libaln.a ../libalex/libalex.a
IDIR = ../include 
INCDIRS = -I. -I../include -I../libafn -I../libhmm \
	-I../libblst -I../libgsq -I../libaln  -I../libalex

LIBS = -lstat -laln -lgsq -lafn -lblst -lalex -libhmm -lstdc++ ${MATH_LIBRARY}
LIBDIRS = -L../libafn/ -L../libblst -L../libgsq -L../libaln -L../libalex -L../libhmm

STDINC = ../include/stdinc.h ../include/sset.h
BIOINC = ../include/residues.h ../include/blosum62.h

OBJS = cluster.o goscan.o jackknife.o oscan2msa.o purge.o \
	tblast.o gpsi_typ.o psm_typ.o editcma.o ema_typ.o clncma.o em_cma.o \
	rma_typ.o gpsi2cma.o tpb_typ.o HMM_typ.o hma_typ.o \
	my_posit.o gpsiblst.o cma_chain.o marg_prob.o rtf_typ.o \
	csp_typ.o rtf_put.o csp_bpps.o 

INC = cluster.h goscan.h jackknife.h purge.h tblast.h gpsi_typ.h \
	oscan2msa.h psm_typ.h editcma.h ema_typ.h clncma.h em_cma.h \
	rma_typ.h gpsi2cma.h tpb_typ.h HMM_typ.h hma_typ.h \
	my_posit.h gpsiblst.h cma_chain.h marg_prob.h rtf_typ.h \
	csp_typ.h $(STDINC) $(ARCHIVE)

# removed coils.o coils.h (needs libncib ?)

libtools.a: $(OBJS) $(INC)
	ar rcs libtools.a $(OBJS) 
#	ranlib libtools.a

.cc.o:
	${CC} $(CFLAGS) $*.cc 
#	${CC} $(CFLAGS) $*.cc ${LIBS}

junk.s: junk.cc $(INC) $(OBJS)
	${CC} -O -S junk.cc ${LIBS} $(INCDIRS)

# tools:

tblast.o: tblast.h purge.h $(INC) $(BIOINC)
cluster.o: $(INC)
editcma.o: $(INC)
cma_chain.o: $(INC) 
em_cma.o: $(INC)
clncma.o: $(INC)
ema_typ.o: $(INC)
jackknife.o: $(INC)
purge.o: purge.h $(BIOINC) $(INC)
gpsi_typ.o: gpsi_typ.h $(BIOINC) $(INC)
marg_prob.o: marg_prob.h $(INC)
# rtf_typ.o: rtf_typ.h rst_typ.h $(INC)
rtf_typ.o: rtf_typ.h $(INC)
csp_typ.o: csp_typ.h $(INC)
csp_bpps.o: csp_typ.h $(INC)
rtf_put.o: rtf_typ.h $(INC)
# rst_typ.o: rst_typ.h $(INC)
# swt_typ.o: swt_typ.h $(INC) 
# hsw_typ.o: hsw_typ.h $(INC)
my_posit.o: my_posit.h $(INC)
gpsiblst.o: gpsiblst.h marg_prob.cc $(INC)
rma_typ.o: rma_typ.h $(BIOINC) $(INC)
HMM_typ.o: HMM_typ.h tpb_typ.h $(BIOINC) $(INC)
hma_typ.o: hma_typ.h HMM_typ.h tpb_typ.h $(BIOINC) $(INC)
tpb_typ.o: tpb_typ.h $(BIOINC) $(INC)
gpsi2cma.o: gpsi2cma.h $(BIOINC) $(INC)
goscan.o: goscan.h $(BIOINC) $(INC)
# probe.o: probe.h $(INC)
# gambit.o: gambit.h $(INC)
oscan2msa.o: purge.h $(BIOINC) $(INC)
psm_typ.o: psm_typ.h $(BIOINC) $(INC)
# coils.o: coils.h $(BIOINC) $(INC) 

clean: 
	@ \rm *.o *.a

