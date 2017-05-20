include src/config.mk
# ==============================================================================
.PHONY: all mcMd mcMd-mpi ddMd tools \
        test-serial test-parallel \
        clean-serial clean-parallel clean clean-bin veryclean \
        html clean-html

# ==============================================================================
# Main build targets

all:
	cd bld/serial; $(MAKE) tools
	cd bld/serial; $(MAKE) mcMd
	cd bld/parallel; $(MAKE) mcMd-mpi
	cd bld/parallel; $(MAKE) ddMd

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
tools:
	cd bld/serial; $(MAKE) tools

# ==============================================================================
# Test targets

test-serial:
	@cd bld/serial/util/tests; $(MAKE) all; $(MAKE) run
	@cd bld/serial/simp/tests; $(MAKE) all; $(MAKE) run
	@cd bld/serial/mcMd/tests; $(MAKE) all; $(MAKE) run
	@cd bld/serial/tools/tests; $(MAKE) all; $(MAKE) run
	@cat bld/serial/util/tests/count > count
	@cat bld/serial/simp/tests/count >> count
	@cat bld/serial/mcMd/tests/count >> count
	@cat bld/serial/tools/tests/count >> count
	@echo " "
	@echo "Summary"
	@cat count
	@rm -f count

test-parallel:
	cd bld/parallel/ddMd/tests; $(MAKE) all; $(MAKE) run
	@cat bld/parallel/ddMd/tests/count >> count
	@echo " "
	@cat count
	@rm -f count

# =========================================================================
# Clean targets

clean-serial:
	cd bld/serial; $(MAKE) clean

clean-parallel:
	cd bld/parallel; $(MAKE) clean

clean-tests:
	cd src/; $(MAKE) clean-tests
	cd bld/serial; $(MAKE) clean-tests
	cd bld/parallel; $(MAKE) clean-tests

clean:
	cd src; $(MAKE) clean
	cd bld/serial; $(MAKE) clean
	cd bld/parallel; $(MAKE) clean

clean-bin:
	-rm -f $(BIN_DIR)/mcSim*
	-rm -f $(BIN_DIR)/mdSim*
	-rm -f $(BIN_DIR)/ddSim*
	-rm -f $(BIN_DIR)/mdPp*
	-rm -f $(BIN_DIR)/*Maker*
 
veryclean:
	cd bld/serial; $(MAKE) veryclean; rm -f makefile configure config.mk
	cd bld/serial; rm -f util/makefile simp/makefile 
	cd bld/serial; rm -f mcMd/makefile ddMd/makefile tools/makefile
	cd bld/serial; rm -f util/tests/makefile simp/tests/makefile 
	cd bld/serial; rm -f mcMd/tests/makefile ddMd/tests/makefile tools/tests/makefile
	cd bld/parallel; $(MAKE) veryclean; rm -f makefile configure config.mk
	cd bld/parallel; rm -f util/makefile simp/makefile 
	cd bld/parallel; rm -f mcMd/makefile ddMd/makefile tools/makefile
	cd bld/parallel; rm -f util/tests/makefile simp/tests/makefile 
	cd bld/parallel; rm -f mcMd/tests/makefile ddMd/tests/makefile tools/tests/makefile
	$(MAKE) clean-bin
	-rm -f $(BIN_DIR)/makeDep
	-rm -f scripts/python/*.pyc
	cd doc; $(MAKE) clean
	cd src; $(MAKE) veryclean

# =========================================================================
# HTML Documentation
 
html:
	cd doc; $(MAKE) html

clean-html:
	cd doc; $(MAKE) clean

# ==============================================================================
