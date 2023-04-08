#
# SHELL = /bin/sh
# can add -O2 to FFLAGS to get maximum optimization
# use -Bstatic to compile with static library.
#

#CFLAGS = -e -Bstatic 
CFLAGS = -g
FC = gfortran

all: lrna crna

lrna: lrna.o lin.o lmrna1.o lmrna2.o lsort.o formid.o multid.o 
	$(FC) $(CFLAGS) -o lrna lrna.o lin.o lmrna1.o lmrna2.o lsort.o formid.o multid.o  

crna: crna.o circ.o cmrna1.o cmrna2.o csort.o formid.o multid.o 
	$(FC) $(CFLAGS) -o crna crna.o circ.o cmrna1.o cmrna2.o csort.o formid.o multid.o 

formid.o: formid.f
	$(FC) $(CFLAGS) -c -o formid.o formid.f 

multid.o: multid.f
	$(FC) $(CFLAGS) -c -o multid.o multid.f 

lrfd.inc: lin.inc
	cp lin.inc rfd.inc

lrna.o: rna.f lrfd.inc
	$(FC) $(CFLAGS) -c -o lrna.o rna.f 

lin.o: lin.f lrfd.inc
	$(FC) $(CFLAGS) -c -o lin.o lin.f 

lmrna1.o: mrna1.f lrfd.inc
	$(FC) $(CFLAGS) -c -o lmrna1.o mrna1.f 

lmrna2.o: mrna2.f lrfd.inc
	$(FC) $(CFLAGS) -c -o lmrna2.o mrna2.f 

lsort.o: sort.f lrfd.inc
	$(FC) $(CFLAGS) -c -o lsort.o sort.f 

crfd.inc: circ.inc
	cp circ.inc rfd.inc

crna.o: rna.f crfd.inc
	$(FC) $(CFLAGS) -c -o crna.o rna.f 

circ.o: circ.f crfd.inc
	$(FC) $(CFLAGS) -c -o circ.o circ.f

cmrna1.o: mrna1.f crfd.inc
	$(FC) $(CFLAGS) -c -o cmrna1.o mrna1.f

cmrna2.o: mrna2.f crfd.inc
	$(FC) $(CFLAGS) -c -o cmrna2.o mrna2.f

csort.o: sort.f crfd.inc
	$(FC) $(CFLAGS) -c -o csort.o sort.f 
