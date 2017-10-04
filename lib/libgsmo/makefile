# makefile for libtools

CC = ${CURRENT_CC} 
# -std=c++11
CFLAGS = $(LIBDIRS) $(INCDIRS) -c 

ARCHIVE = ../libafn/libafn.a ../libblst/libblst.a ../libhmm/libhmm.a ../libgsq/libgsq.a \
	../libaln/libaln.a ../libalex/libalex.a ../libtools/libtools.a
IDIR = ../include 
INCDIRS = -I. -I../include -I../libafn -I../libhmm -I../libgsq \
	-I../libblst -I../libgsq  -I../libaln  -I../libalex -I../libtools

LIBS = -lstat -laln -lgsq -lafn -lblst -lalex -libhmm -libtools -lstdc++ ${MATH_LIBRARY}
LIBDIRS = -L../libafn/ -L../libblst -L../libgsq -L../libaln -L../libalex -L../libhmm -L../libtools

STDINC = ../include/stdinc.h ../include/sset.h 
BIOINC = ../include/residues.h ../include/blosum62.h

OBJS = gsm_typ.o gsm_blks.o gsm_smpl.o gsm_init.o gsm_srch.o gmb_typ.o gmb_put.o gmb_smpl.o \
	gmb_sticky.o gmb_subaln.o gmb_cols.o gmb_p2p.o gsm_plus.o gmb_clstr.o gmb_clocl.o gmb_junk.o \
	gmb_debug.o ssx_typ.o ssx_init.o ssx_scores.o ndl_p2p.o jlh_typ.o jlh_put.o jlh_hmm.o jlh_p2p.o \
	ndl_typ.o ndl_hmm.o ndl_smpl.o ndl_debug.o cma_gmb.o str_typ.o subaln_gmb.o \
	ssx_dmp.o gpr_typ.o ppg_typ.o

SSX_INC = ssx_typ.h jlh_typ.h ndl_typ.h p2p_typ.h lgm_typ.h $(STDINC) $(ARCHIVE)
GMB_INC = gmb_typ.h ssx_typ.h jlh_typ.h ndl_typ.h cma_gmb.h lgm_typ.h str_typ.h gpr_typ.h $(STDINC) $(ARCHIVE)
GSM_INC = gsm_typ.h gmb_typ.h ssx_typ.h lgm_typ.h ppg_typ.h $(GMB_INC) $(STDINC) $(ARCHIVE)
JLH_INC = jlh_typ.h p2p_typ.h $(STDINC) $(ARCHIVE)
NDL_INC = ndl_typ.h p2p_typ.h lgm_typ.h $(STDINC) $(ARCHIVE)
DMS_INC = dms_data.h $(STDINC) $(ARCHIVE)


# removed coils.o coils.h (needs libncbi ?)

libgsmo.a: $(OBJS) $(GSM_INC)
	\rm -f libgsmo.a
	ar rcs libgsmo.a $(OBJS) 
#	ranlib libgsmo.a

.cc.o:
	${CC} $(CFLAGS) $*.cc 
#	${CC} $(CFLAGS) $*.cc ${LIBS}

# object files:

str_typ.o: $(STD_INC)
ssx_typ.o: $(SSX_INC)
ssx_init.o: $(SSX_INC)

gpr_typ.o: gpr_typ.h $(STD_INC) $(ARCHIVE)
ppg_typ.o: ppg_typ.h $(STD_INC) $(ARCHIVE)
gsm_typ.o:  $(GSM_INC)
gsm_plus.o:  $(GSM_INC)
gismo_plus.o:  $(GSM_INC)
gsm_smpl.o:  $(GSM_INC)
gsm_blks.o:  $(GSM_INC)
gsm_init.o:  $(GSM_INC)
gsm_srch.o: $(GSM_INC)

gmb_typ.o: $(GMB_INC)
gmb_p2p.o: $(GMB_INC)
ssx_scores.o: $(GMB_INC)
ssx_dmp.o: $(GMB_INC)
gmb_put.o: $(GMB_INC)
gmb_debug.o: $(GMB_INC)
gmb_junk.o: $(GMB_INC)
gmb_smpl.o: $(GMB_INC)
gmb_cols.o: $(GMB_INC)
gmb_sticky.o: $(GMB_INC)
gmb_clstr.o: $(GMB_INC)
gmb_clocl.o: $(GMB_INC)
gmb_subaln.o: $(GMB_INC)

jlh_typ.o: $(JLH_INC)
jlh_put.o: $(JLH_INC)
jlh_hmm.o: $(JLH_INC)
jlh_p2p.o: $(JLH_INC)

ndl_typ.o: $(NDL_INC)
ndl_hmm.o: $(NDL_INC)
ndl_debug.o: $(NDL_INC)
ndl_smpl.o: $(NDL_INC)
ndl_p2p.o: $(NDL_INC)

subaln_gmb.o: $(GMB_INC)
cma_gmb.o: $(GMB_INC)

clean: 
	@ \rm *.o *.a

