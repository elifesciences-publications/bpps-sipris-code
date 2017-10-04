# makefile for libpdb

CC = ${CURRENT_CC}
CFLAGS = $(LIBDIRS) $(INCDIRS) -c

ARCHIVE = ../libafn/libafn.a
IDIR = ../include 
INCDIRS = -I. -I../include -I../libafn 

LIBS = -lstat -lafn ${MATH_LIBRARY}
LIBDIRS = -L../libafn/ 

STDINC = ../include/stdinc.h ../include/sset.h
BIOINC = ../include/residues.h ../include/blosum62.h

OBJS = pdb.o pdb_misc.o atom.o rasmol.o surface.o contacts.o res_typ.o radius.o \
  geometry.o sip_typ.o phi_psi.o ras_typ.o vsi_typ.o vsi_heap.o pml_typ.o \
  vsi_ras.o vsi_pml.o

INC = pdb.h atom.h rasmol.h surface.h exposed.h contacts.h res_typ.h \
 radius.h geometry.h sip_typ.h phi_psi.h ras_typ.h vsi_typ.h pml_typ.h \
 $(STDINC) $(ARCHIVE)

libpdb.a: $(OBJS) $(INC)
	ar rcs libpdb.a $(OBJS)
#	ranlib libpdb.a

.cc.o:
	${CC} $(CFLAGS) $*.cc 
#	${CC} $(CFLAGS) $*.cc ${LIBS}


# basic 3D structure routines
pdb.o: pdb.h $(INC)
pdb_misc.o: pdb.h $(INC)
geometry.o: geometry.h $(INC)
sip_typ.o: sip_typ.h $(INC)
atom.o: atom.h radius.h $(INC)

# structural property routines
surface.o: surface.h exposed.h $(INC)
contacts.o: contacts.h $(INC)

# other structural routines:
#strctprof.o: strctprof.h $(INC)
dsc.o: dsc.h $(INC)
rasmol.o: rasmol.h $(INC)
res_typ.o: res_typ.h $(INC)
radius.o: radius.h $(INC)

pml_typ.o: ${INC} 
ras_typ.o: ${INC} 
vsi_heap.o: ${INC}
vsi_typ.o: ${INC}
vsi_ras.o: ${INC}
vsi_pml.o: ${INC}

phi_psi.o: phi_psi.cc $(INC)

clean: 
	@ \rm *.o *.a

