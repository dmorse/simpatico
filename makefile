include src/compiler.mk
# ==============================================================================
.PHONY: mcMd ddMd mcMd-mpi \
        test-mcMd test-mcMd-mpi test-ddMd clean clean-bin \
        clean-mcMd clean-ddMd clean clean-bin \
        html clean-html veryclean

all:
	-rm -f log
	cd src; make mcMd
	make test-mcMd
	cd src; make mcMd-mpi
	make test-mcMd-mpi
	cd src; make ddMd
	make test-ddMd
	cat log

mcMd:
	cd src; make mcMd
	-rm -f log
	make test-mcMd
	cat log

mcMd-mpi: 
	cd src; make mcMd-mpi
	-rm -f log
	make test-mcMd-mpi
	cat log

ddMd:
	cd src; make ddMd
	-rm -f log
	make test-ddMd
	cat log

# ==============================================================================
test-mcMd:
	./configure -m0
	cd tests/util; make all; ./Test > log;
	echo `grep failed tests/util/log` ", "\
             `grep successful tests/util/log` "in tests/util/log" >> log
	cd tests/inter; make all; ./Test > log;
	echo `grep failed tests/inter/log` ", "\
             `grep successful tests/inter/log` "in tests/inter/log" >> log
	cd tests/mcMd; make all; ./Test > log;
	echo `grep failed tests/mcMd/log` ", "\
             `grep successful tests/mcMd/log` "in tests/mcMd/log" >> log

test-mcMd-mpi:
	./configure -m1
	cd tests/util/mpi; make all; $(MPIRUN) 2 ./Test > log
	echo `grep failed tests/util/mpi/log` ", "\
             `grep successful tests/util/mpi/log` "in tests/util/mpi/log" >> log
	cd tests/util/param/mpi; make all; $(MPIRUN) 2 ./MpiTest > log
	echo `grep failed tests/util/param/mpi/log` ", "\
             `grep successful tests/util/param/mpi/log` "in tests/util/param/mpi/log" >> log

test-ddMd:
	./configure -m1
	cd tests/ddMd; make all; $(MPIRUN) 6 ./Test > log;
	echo `grep failed tests/ddMd/log` ", "\
             `grep successful tests/ddMd/log` "in tests/ddMd/log" >> log
# ==============================================================================
# Clean targets

clean-mcMd:
	cd src; make clean-mcMd
	cd tests/util; make clean;
	cd tests/inter; make clean;
	cd tests/mcMd; make clean;

clean-ddMd:
	cd src; make clean-ddMd
	cd tests/util; make clean;
	cd tests/inter; make clean;
	cd tests/ddMd; make clean;

clean:
	cd src; make clean
	cd tests; make clean
	-rm -f log

clean-bin:
	cd src; make clean-bin
# ==============================================================================
# HTML Documentation
 
html:
	cd doc; make html

clean-html:
	cd doc; make clean

# ==============================================================================
# Remove all generated files

veryclean:
	-rm -f log
	cd html; make clean
	cd tests; make clean
	cd src; make veryclean

# ==============================================================================
