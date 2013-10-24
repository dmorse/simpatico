include src/compiler.mk
# ==============================================================================
.PHONY: mcMd ddMd mcMd-mpi test \
        clean-serial clean-parallel clean \
        html clean-html 

# ==============================================================================
# Main targets

all:
	cd build/serial; make mcMd
	cd build/parallel; make ddMd
	cd build/parallel; make mcMd-mpi

mcMd:
	cd build/serial; make mcMd

mcMd-mpi: 
	cd build/parallel; make mcMd-mpi

ddMd:
	cd build/parallel; make ddMd

test:
	cd build/serial/util/tests; make all; make run
	cd build/serial/inter/tests; make all; make run
	cd build/serial/mcMd/tests; make all; make run
	cd build/parallel/ddMd/tests; make all; make run
	@cat build/serial/util/tests/count > count
	@cat build/serial/inter/tests/count >> count
	@cat build/serial/mcMd/tests/count >> count
	@cat build/parallel/ddMd/tests/count >> count
	@cat count
	@rm -f count


# ==============================================================================
# Clean targets

clean-serial:
	cd build/serial; make clean

clean-parallel:
	cd build/parallel; make clean

clean:
	cd src; make clean
	cd build/serial; make clean
	cd build/parallel; make clean

veryclean:
	cd html; make clean
	cd build/serial; make veryclean
	cd build/parallel; make veryclean
	cd src; make veryclean

# ==============================================================================
# HTML Documentation
 
html:
	cd doc; make html

clean-html:
	cd doc; make clean

# ==============================================================================
