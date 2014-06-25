include src/config.mk
# ==============================================================================
.PHONY: all mcMd mcMd-mpi ddMd \
        test-serial test-parallel \
        clean-serial clean-parallel clean clean-bin veryclean \
        html clean-html

# ==============================================================================
# Main build targets

all:
	cd obj/serial; make mcMd
	cd obj/parallel; make ddMd
	cd obj/parallel; make mcMd-mpi

# Build serial mcSim and mdSim MC and MD programs in obj/serial
mcMd:
	cd obj/serial; make mcMd

# Build embarassingly parallel mcSim and mdSim programs in obj/parallel
mcMd-mpi: 
	cd obj/parallel; make mcMd-mpi

# Build parallel mdSim MD program in obj/parallel
ddMd:
	cd obj/parallel; make ddMd

# ==============================================================================
# Test targets

test-serial:
	@cd obj/serial/util/tests; make all; make run
	@cd obj/serial/inter/tests; make all; make run
	@cd obj/serial/mcMd/tests; make all; make run
	@cat obj/serial/util/tests/count > count
	@cat obj/serial/inter/tests/count >> count
	@cat obj/serial/mcMd/tests/count >> count
	@cat count
	@rm -f count

test-parallel:
	cd obj/parallel/ddMd/tests; make all; make run
	@cat obj/parallel/ddMd/tests/count >> count
	@cat count
	@rm -f count


# =========================================================================
# Clean targets

clean-serial:
	cd obj/serial; make clean

clean-parallel:
	cd obj/parallel; make clean

clean:
	cd src; make clean
	cd obj/serial; make clean
	cd obj/parallel; make clean

clean-bin:
	-rm -f $(BIN_DIR)/mcSim*
	-rm -f $(BIN_DIR)/mdSim*
	-rm -f $(BIN_DIR)/ddSim*
 
veryclean:
	cd obj/serial; make veryclean; rm -f makefile configure
	cd obj/serial; rm -f util/makefile inter/makefile mcMd/makefile ddMd/makefile 
	cd obj/parallel; make veryclean; rm -f makefile configure
	cd obj/parallel; rm -f util/makefile inter/makefile mcMd/makefile ddMd/makefile 
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
