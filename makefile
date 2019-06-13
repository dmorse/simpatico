include src/config.mk
# ==============================================================================
.PHONY: all mcMd mcMd-mpi ddMd mdPp \
        test-serial test-parallel \
        clean-serial clean-parallel clean clean-bin veryclean \
        html clean-html help

# =========================================================================
# Generate help dialog listing makefile targets

help: ## This help dialog.
	@IFS=$$'\n' ; \
        help_lines=(`fgrep -h "##" $(MAKEFILE_LIST) | fgrep -v fgrep | sed -e 's/\\$$//'`); \
        for help_line in $${help_lines[@]}; do \
           IFS=$$'#' ; \
           help_split=($$help_line) ; \
           help_command=`echo $${help_split[0]} | sed -e 's/^ *//' -e 's/ *$$//'` ; \
           help_info=`echo $${help_split[2]} | sed -e 's/^ *//' -e 's/ *$$//'` ; \
           printf "%-18s %s\n" $$help_command $$help_info ; \
        done

# ==============================================================================
# Main build targets

all:      ## Build all simpatico programs (mdSim, mcSim, ddSim, mdPp)
	cd bld/serial; $(MAKE) mdPp
	cd bld/serial; $(MAKE) mcMd
	cd bld/parallel; $(MAKE) mcMd-mpi
	cd bld/parallel; $(MAKE) ddMd

mcMd:     ## Build serial mcSim and mdSim programs
	cd bld/serial; $(MAKE) mcMd

mcMd-mpi: ## Build multi-processor mcSim and mdSim programs 
	cd bld/parallel; $(MAKE) mcMd-mpi

ddMd:     ## Build ddSim parallel MD program 
	cd bld/parallel; $(MAKE) ddMd

mdPp:    ## Build mdPp MD postprocessor program
	cd bld/serial; $(MAKE) mdPp

# ==============================================================================
# Test targets

test-serial:  ## Run unit tests for serial code (MPI disabled)
	@cd bld/serial/util/tests; $(MAKE) quiet
	@cd bld/serial/simp/tests; $(MAKE) quiet
	@cd bld/serial/mcMd/tests; $(MAKE) quiet
	@cd bld/serial/mdPp/tests; $(MAKE) quiet
	@cat bld/serial/util/tests/count > count
	@cat bld/serial/simp/tests/count >> count
	@cat bld/serial/mcMd/tests/count >> count
	@cat bld/serial/mdPp/tests/count >> count
	@echo " "
	@echo "Summary"
	@cat count
	@rm -f count

test-parallel: ## Run unit tests for parallel code (MPI enabled)
	cd bld/parallel/ddMd/tests; $(MAKE) quiet
	@cat bld/parallel/ddMd/tests/count >> count
	@echo " "
	@cat count
	@rm -f count

# =========================================================================
# Clean targets

clean:          ## Remove files generated by compilation (leave executables)
	cd src; $(MAKE) clean
	cd bld/serial; $(MAKE) clean
	cd bld/parallel; $(MAKE) clean

clean-serial:   ## Clean build directory bld/serial 
	cd bld/serial; $(MAKE) clean

clean-parallel: ## Clean build directory bld/parallel
	cd bld/parallel; $(MAKE) clean

clean-tests:    ## Clean all test directories
	cd src/; $(MAKE) clean-tests
	cd bld/serial; $(MAKE) clean-tests
	cd bld/parallel; $(MAKE) clean-tests

clean-bin:      ## Remove all generated executables 
	-rm -f $(BIN_DIR)/mcSim*
	-rm -f $(BIN_DIR)/mdSim*
	-rm -f $(BIN_DIR)/ddSim*
	-rm -f $(BIN_DIR)/mdPp*
	-rm -f $(BIN_DIR)/*Maker*
 
veryclean:     ## Revert to as-distributed state, pre-setup
	cd bld/serial; $(MAKE) veryclean; rm -f makefile configure config.mk
	cd bld/serial; rm -f util/makefile simp/makefile 
	cd bld/serial; rm -f mcMd/makefile ddMd/makefile mdPp/makefile
	cd bld/serial; rm -f util/tests/makefile simp/tests/makefile 
	cd bld/serial; rm -f mcMd/tests/makefile ddMd/tests/makefile
	cd bld/serial; rm -f mdPp/tests/makefile
	cd bld/parallel; $(MAKE) veryclean; rm -f makefile configure config.mk
	cd bld/parallel; rm -f util/makefile simp/makefile 
	cd bld/parallel; rm -f mcMd/makefile ddMd/makefile mdPp/makefile
	cd bld/parallel; rm -f util/tests/makefile simp/tests/makefile 
	cd bld/parallel; rm -f mcMd/tests/makefile ddMd/tests/makefile 
	cd bld/parallel; rm -f mdPp/tests/makefile
	-rm -f scripts/python/*.pyc
	-rm -f $(BIN_DIR)/*
	cd doc; $(MAKE) clean
	cd src; $(MAKE) veryclean

# =========================================================================
# HTML Documentation
 
html:   ## Build html manual
	cd doc; $(MAKE) html

clean-html:   ## Remove generated html manual files
	cd doc; $(MAKE) clean

