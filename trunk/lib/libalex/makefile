# makefile for libalex

CC = ${CURRENT_CC}

CFLAGS = $(LIBDIRS) $(INCDIRS) -c

ARCHIVE = ../libafn/libafn.a ../libblst/libblst.a ../libgsq/libgsq.a ../libaln/libaln.a \
	../libstat/libstat.a 
IDIR = ../include 
INCDIRS = -I. -I../include -I../libafn -I../libblst  -I../libgsq -I../libaln \
        -I../libstat 
LIBS = -laln -lgsq -lblst -lafn -lstat ${MATH_LIBRARY}
LIBDIRS = -L../libafn -L../libblst -L../libgsq -L../libaln -L../libstat 
#STDINC = ../include/stdinc.h  ../include/residues.h ./infopair.h ./gold_stand.h
STDINC = ../include/stdinc.h  ../include/residues.h 

OBJS = dp.o new-dsc.o n2a.o NWprof.o sseq.o sba.o sdp.o s-w.o routines.o psg.o

INC = dp.h new-dsc.h n2a.h NWprof.h sba.h sseq.h sdp.h s-w.h routines.h psg.h $(ARCHIVE) $(STDINC)

libalex.a: $(OBJS) $(INC)
	ar rcs libalex.a $(OBJS)

.cc.o:
	${CC} $(CFLAGS) $*.cc 
#	${CC} $(CFLAGS) $*.cc ${LIBS}

# routines:
dp.o: $(INC)
new-dsc.o: $(INC)
n2a.o: $(INC)
NWprof.o: $(INC)
sseq.o: $(INC)
sba.o: $(INC)
sdp.o: $(INC)
# gor_typ.o: $(INC)
s-w.o: $(INC)
routines.o: $(INC)
psg.o: $(INC)

clean: 
	@ \rm *.o *.a

