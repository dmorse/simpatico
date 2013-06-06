# ==============================================================================
.PHONY: mcMd ddMd mcMd-mpi \
        test-mcMd test-ddMd clean clean-bin \
        clean-mcMd clean-ddMd clean clean-bin \
        html clean-html veryclean

mcMd:
	cd src; make mcMd
	-rm -f log
	make test-mcMd
	cat log
	rm log

test-mcMd:
	cd tests/util; make all; ./Test > log;
	echo `grep failed tests/util/log` \
             `grep successful tests/util/log` "in tests/util/log" >> log
	cd tests/inter; make all; ./Test > log;
	echo `grep failed tests/inter/log` \
             `grep successful tests/inter/log` "in tests/inter/log" >> log
	cd tests/mcMd; make all; ./Test > log;
	echo `grep failed tests/mcMd/log` \
             `grep successful tests/mcMd/log` "in tests/mcMd/log" >> log

mcMd-mpi: 
	cd src; make mcMd-mpi

ddMd:
	cd src; make ddMd

clean-mcMd:
	cd src; make clean-mcMd

clean-ddMd:
	cd src; make clean-ddMd

clean:
	cd src; make clean

clean-bin:
	cd src; make clean-bin
 
# ==============================================================================
html:
	cd doc; make html

clean-html:
	cd doc; make clean

# ==============================================================================
veryclean:
	cd html; make clean
	cd src; make veryclean

# ==============================================================================
