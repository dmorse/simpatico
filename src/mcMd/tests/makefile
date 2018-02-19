BLD_DIR_REL =../..
include $(BLD_DIR_REL)/config.mk
include $(BLD_DIR)/mcMd/config.mk
include $(BLD_DIR)/simp/config.mk
include $(BLD_DIR)/util/config.mk
include $(SRC_DIR)/mcMd/patterns.mk
include $(SRC_DIR)/mcMd/sources.mk
include $(SRC_DIR)/simp/sources.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/mcMd/tests/sources.mk

TEST=$(BLD_DIR)/mcMd/tests/Test
.SECONDARY: $(mcMd_tests_OBJS)

all: $(mcMd_tests_EXES) 

run: $(TEST)
	$(TEST) $(SRC_DIR)/mcMd/tests/
	cd $(SRC_DIR)/mcMd/tests; $(MAKE) clean-outputs

quiet: $(TEST)
	$(TEST) $(SRC_DIR)/mcMd/tests/ > log
	@echo `grep failed log` ", "\
             `grep successful log` "in mcMd/tests/log" > count
	@cat count
	@cd $(SRC_DIR)/mcMd/tests; $(MAKE) clean-outputs

clean:
	rm -f $(mcMd_tests_EXES)
	rm -f $(mcMd_tests_OBJS)
	rm -f $(mcMd_tests_OBJS:.o=.d)
	rm -f log count
ifeq ($(BLD_DIR),$(SRC_DIR))
	cd analyzers; $(MAKE) clean
	cd chemistry; $(MAKE) clean
	cd generators; $(MAKE) clean
	cd mcSimulation; $(MAKE) clean
	cd mdSimulation; $(MAKE) clean
	cd potentials; $(MAKE) clean
	cd mdIntegrators; $(MAKE) clean
	cd neighbor; $(MAKE) clean
	cd simulation; $(MAKE) clean
else
	cd $(SRC_DIR)/mcMd/tests; $(MAKE) clean-outputs
endif

clean-outputs:
	@cd chemistry; $(MAKE) clean-outputs
	@cd mcSimulation; $(MAKE) clean-outputs
	@cd mdSimulation; $(MAKE) clean-outputs

-include $(mcMd_tests_OBJS:.o=.d)
-include $(mcMd_OBJS:.o=.d)
-include $(simp_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)

