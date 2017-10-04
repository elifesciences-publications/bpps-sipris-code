# makefile for scan
CC = ${CURRENT_CC}
#GCC = g++ -static -g -O3
GCC = ${CURRENT_CC}
CFLAGS = $(INCDIRS) $(LIBDIRS)
ARCHIVES = ../libafn/libafn.a ../libpdb/libpdb.a ../libgsq/libgsq.a ../libaln/libaln.a \
	../libblst/libblst.a ../libtools/libtools.a \
	../libstat/libstat.a ../libalex/libalex.a ../libtax/libtax.a 

#CAUTION: make sure that libvsi is not included as both libvsi and libchn archive lex/yacc files.

INCDIRS = -I. -I../include -I../libafn -I../libblst -I../libgsq -I../libaln \
	-I../libpdb -I../libstat -I../libtools -I../libtax \
	-I../libalex 

LIBS = -ltax -ltools -lalex -laln -lgsq -lblst -lpdb -lafn -lstat ${MATH_LIBRARY}
LIBDIRS = -L../libafn -L../libgsq -L../libaln -L../libpdb -L../libblst \
	-L../libstat -L../libalex -L../libtools -L../libtax 

STDINC = ../include/stdinc.h ../include/residues.h

OBJS = lpr_typ.o chn_aln.o chn_pps.o sqd_typ.o ctn_typ.o clique.o chn_typ.o cth_typ.o che_typ.o \
	che_junk.o chn_rtf.o chn_read.o chn_swt.o chn_put.o chn_cha.o hpt_typ.o hpt_sort.o \
	hpt_operators.o hpt_init.o hpt_tree.o btn_typ.o hmm_typ.o tree2hpt.o \
	cnh_typ.o c2h_typ.o lex.yy.o y.tab.o 

INC = lpr_typ.h chn_aln.h chn_pps.h chn_typ.h sqd_typ.h ctn_typ.h clique.h cth_typ.h che_typ.h pcr_typ.h \
	hpt_typ.h tree2hpt.h btn_typ.h hmm_typ.h cnh_typ.h c2h_typ.h bpcp_typ.h \
	$(ARCHIVES) $(STDINC)

libchn.a: $(OBJS) $(INC)
	\rm -f libchn.a
	ar rcs libchn.a $(OBJS)

.cc.o:	$(INC)
	${CC} ${CFLAGS} -c $*.cc

lpr_typ.o: lpr_typ.h $(INC) 
che_typ.o: $(INC)
che_junk.o: $(INC)
chn_aln.o: $(INC)
chn_rtf.o: $(INC)
chn_read.o: $(INC)
chn_swt.o: $(INC)
chn_put.o: $(INC)
chn_cha.o: $(INC)
chn_pps.o: $(INC)
chn_typ.o: $(INC)
c2h_typ.o: $(INC)
cnh_typ.o: $(INC)
hmm_typ.o: $(INC)
clique.o: $(INC)
cth_typ.o: $(INC)
sqd_typ.o: $(INC)
ctn_typ.o: $(INC)
hpt_typ.o: $(INC)
hpt_init.o: $(INC)
hpt_operators.o: $(INC)
hpt_sort.o: $(INC)
hpt_tree.o: $(INC)

# tree2hpt program...
YACC = yacc -vd
OBJS_T2H = lex.yy.o y.tab.o btn_typ.o tree2hpt.o
INC_T2H = tree2hpt.h btn_typ.h 

y.tab.o: y.tab.c $(ARCHIVE)
	${GCC} ${CFLAGS} -c y.tab.c 

y.tab.c: tree2hpt.y ${INC} $(ARCHIVE)
	${YACC} tree2hpt.y 

lex.yy.o: y.tab.o tree2hpt.l y.tab.h ${INC} $(ARCHIVE)
	flex tree2hpt.l
	${GCC} ${CFLAGS} -c lex.yy.c -lfl

btn_typ.o: btn_typ.cc ${INC}
tree2hpt.o: tree2hpt.cc ${INC}
y.tab.h: y.tab.c 

clean: 
	@ \rm lex.yy.c y.tab.c *.o *.a y.tab.h

