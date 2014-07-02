include src/config.mk
# ==============================================================================
.PHONY: all mcMd mcMd-mpi ddMd \
        test-serial test-parallel \
        clean-serial clean-parallel clean clean-bin veryclean \
        html clean-html

# ==============================================================================
# Main build targets

all:
	cd bld/serial; make mcMd
	cd bld/parallel; make ddMd
	cd bld/parallel; make mcMd-mpi

# Build serial mcSim and mdSim MC and MD programs in bld/serial
mcMd:
	cd bld/serial; make mcMd

# Build embarassingly parallel mcSim and mdSim programs in bld/parallel
mcMd-mpi: 
	cd bld/parallel; make mcMd-mpi

# Build parallel mdSim MD program in bld/parallel
ddMd:
	cd bld/parallel; make ddMd

# ==============================================================================
# Test targets

test-serial:
	@cd bld/serial/util/tests; make all; make run
	@cd bld/serial/inter/tests; make all; make run
	@cd bld/serial/mcMd/tests; make all; make run
	@cat bld/serial/util/tests/count > count
	@cat bld/serial/inter/tests/count >> count
	@cat bld/serial/mcMd/tests/count >> count
	@cat count
	@rm -f count

test-parallel:
	cd bld/parallel/ddMd/tests; make all; make run
	@cat bld/parallel/ddMd/tests/count >> count
	@cat count
	@rm -f count


# =========================================================================
# Clean targets

clean-serial:
	cd bld/serial; make clean

clean-parallel:
	cd bld/parallel; make clean

clean:
	cd src; make clean
	cd bld/serial; make clean
	cd bld/parallel; make clean

clean-bin:
	-rm -f $(BIN_DIR)/mcSim*
	-rm -f $(BIN_DIR)/mdSim*
	-rm -f $(BIN_DIR)/ddSim*
 
veryclean:
	cd bld/serial; make veryclean; rm -f makefile configure
	cd bld/serial; rm -f util/makefile inter/makefile mcMd/makefile ddMd/makefile 
	cd bld/parallel; make veryclean; rm -f makefile configure
	cd bld/parallel; rm -f util/makefile inter/makefile mcMd/makefile ddMd/makefile 
	cd doc; make clean
	make clean-bin
	cd src; make veryclean

# =========================================================================
# HTML Documentation
 
html:
	cd doc; make html

clean-html:
	cd doc; make clean

# ==============================================================================
