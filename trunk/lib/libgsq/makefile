# makefile for libgsq

CC = ${CURRENT_CC}
CFLAGS = $(LIBDIRS) $(INCDIRS) -c

ARCHIVE = ../libafn/libafn.a 
INCDIRS = -I. -I../include -I../libafn 

LIBS = -lstat -lblst -lafn -lstdc++ ${MATH_LIBRARY}
LIBDIRS = -L../libafn/ 

STDINC = ../include/stdinc.h ../include/sset.h
BIOINC = ../include/residues.h ../include/blosum62.h

OBJS = gss_typ.o gsq_typ.o gsq_operate.o gsq_init.o gsq_put.o

INC = gss_typ.h gsq_typ.h $(STDINC) $(ARCHIVE)

libgsq.a: $(OBJS) $(INC)
	ar rcs libgsq.a $(OBJS) 
#	ranlib libgsq.a

.cc.o:
	${CC} $(CFLAGS) $*.cc 
#	${CC} $(CFLAGS) $*.cc ${LIBS}

junk.s: junk.cc $(INC) $(OBJS)
	${CC} -O -S junk.cc ${LIBS} $(INCDIRS)

gsq_typ.o: $(INC)
gss_typ.o: gss_typ.h $(INC)

clean: 
	@ \rm *.o *.a

