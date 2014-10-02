include src/config.mk
# ==============================================================================
.PHONY: all mcMd mcMd-mpi ddMd spAn \
        test-serial test-parallel \
        clean-serial clean-parallel clean clean-bin veryclean \
        html clean-html

# ==============================================================================
# Main build targets

all:
	cd bld/serial; $(MAKE) mcMd
	cd bld/parallel; $(MAKE) ddMd
	cd bld/parallel; $(MAKE) mcMd-mpi
	cd bld/serial; $(MAKE) spAn

# Build serial mcSim and mdSim MC and MD programs in bld/serial
mcMd:
	cd bld/serial; $(MAKE) mcMd

# Build embarassingly parallel mcSim and mdSim programs in bld/parallel
mcMd-mpi: 
	cd bld/parallel; $(MAKE) mcMd-mpi

# Build parallel mdSim MD program in bld/parallel
ddMd:
	cd bld/parallel; $(MAKE) ddMd

# Build single-processor analysis program in bld/serial
spAn:
	cd bld/serial; $(MAKE) spAn

# ==============================================================================
# Test targets

test-serial:
	@cd bld/serial/util/tests; $(MAKE) all; $(MAKE) run
	@cd bld/serial/inter/tests; $(MAKE) all; $(MAKE) run
	@cd bld/serial/mcMd/tests; $(MAKE) all; $(MAKE) run
	@cd bld/serial/spAn/tests; $(MAKE) all; $(MAKE) run
	@cat bld/serial/util/tests/count > count
	@cat bld/serial/inter/tests/count >> count
	@cat bld/serial/mcMd/tests/count >> count
	@cat bld/serial/spAn/tests/count >> count
	@cat count
	@rm -f count

test-parallel:
	cd bld/parallel/ddMd/tests; $(MAKE) all; $(MAKE) run
	@cat bld/parallel/ddMd/tests/count >> count
	@cat count
	@rm -f count

# =========================================================================
# Clean targets

clean-serial:
	cd bld/serial; $(MAKE) clean

clean-parallel:
	cd bld/parallel; $(MAKE) clean

clean:
	cd src; $(MAKE) clean
	cd bld/serial; $(MAKE) clean
	cd bld/parallel; $(MAKE) clean

clean-bin:
	-rm -f $(BIN_DIR)/mcSim*
	-rm -f $(BIN_DIR)/mdSim*
	-rm -f $(BIN_DIR)/ddSim*
 
veryclean:
	cd bld/serial; $(MAKE) veryclean; rm -f makefile configure
	cd bld/serial; rm -f util/makefile inter/makefile mcMd/makefile ddMd/makefile spAn/makefile
	cd bld/parallel; $(MAKE) veryclean; rm -f makefile configure
	cd bld/parallel; rm -f util/makefile inter/makefile mcMd/makefile ddMd/makefile  spAn/makefile
	cd doc; $(MAKE) clean
	$(MAKE) clean-bin
	cd src; $(MAKE) veryclean

# =========================================================================
# HTML Documentation
 
html:
	cd doc; $(MAKE) html

clean-html:
	cd doc; $(MAKE) clean

# ==============================================================================
