include src/compiler.mk
# ==============================================================================
.PHONY: mcMd ddMd mcMd-mpi test \
        clean-serial clean-parallel clean \
        html clean-html 

# ==============================================================================
# Main targets

all:
	cd obj/serial; make mcMd
	cd obj/parallel; make ddMd
	cd obj/parallel; make mcMd-mpi

mcMd:
	cd obj/serial; make mcMd

mcMd-mpi: 
	cd obj/parallel; make mcMd-mpi

ddMd:
	cd obj/parallel; make ddMd

test:
	cd obj/serial/util/tests; make all; make run
	cd obj/serial/inter/tests; make all; make run
	cd obj/serial/mcMd/tests; make all; make run
	cd obj/parallel/ddMd/tests; make all; make run
	@cat obj/serial/util/tests/count > count
	@cat obj/serial/inter/tests/count >> count
	@cat obj/serial/mcMd/tests/count >> count
	@cat obj/parallel/ddMd/tests/count >> count
	@cat count
	@rm -f count


# ==============================================================================
# Clean targets

clean-serial:
	cd obj/serial; make clean

clean-parallel:
	cd obj/parallel; make clean

clean:
	cd src; make clean
	cd obj/serial; make clean
	cd obj/parallel; make clean

veryclean:
	cd html; make clean
	cd obj/serial; make veryclean
	cd obj/parallel; make veryclean
	cd src; make veryclean

# ==============================================================================
# HTML Documentation
 
html:
	cd doc; make html

clean-html:
	cd doc; make clean

# ==============================================================================
