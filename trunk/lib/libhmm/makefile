# makefile for libhmm

CC = ${CURRENT_CC}
CFLAGS = $(LIBDIRS) $(INCDIRS) -c

ARCHIVE = ../libafn/libafn.a ../libgsq/libgsq.a ../libaln/libaln.a \
	../libstat/libstat.a ../libblst/libblst.a
IDIR = ../include 
INCDIRS = -I. -I../include -I../include -I../libafn -I../libgsq -I../libaln \
        -I../libstat -I../libblst 
LIBS = -laln -lgsq -lafn -lstat -lblst 
LIBDIRS = -L../libafn -L../libgsq -L../libaln -L../libblst \
	-L../libstat 

STDINC = ../include/stdinc.h ../include/residues.h 

OBJS = hmm_mem.o gpsihmm.o ccp_typ.o tmp_typ.o pda_typ.o sqs_typ.o path2sap.o selexCMSA.o \
	vtb_typ.o shm_typ.o stpb_typ.o lkp_typ.o hmm_math.o sstp_typ.o

INC = hmm_mem.h gpsihmm.h ccp_typ.h tmp_typ.h pda_typ.h sqs_typ.h path2sap.h selexCMSA.h \
	vtb_typ.h shm_typ.h stpb_typ.h lkp_typ.h hmm_math.h sstp_typ.h $(ARCHIVE) $(STDINC)

libhmm.a: $(OBJS) $(INC)
	ar rcs libhmm.a $(OBJS)

.cc.o:
	${CC} $(CFLAGS) $*.cc 
#	${CC} $(CFLAGS) $*.cc ${LIBS}

# routines:
hmm_mem.o: $(INC)
gpsihmm.o: $(INC)
ccp_typ.o: $(INC)
tmp_typ.o: $(INC)
pda_typ.o: $(INC)
sqs_typ.o: $(INC)
sstp_typ.o: $(INC)
path2sap.o: $(INC)
selexCMSA.o: $(INC)
vtb_typ.o: $(INC)
shm_typ.o: $(INC)
stpb_typ.o: $(INC)
hmm_math.o: hmm_math.h $(STDINC) $(ARCHIVE)
lkp_typ.o: $(INC)

clean:
	@ \rm *.o *.a

