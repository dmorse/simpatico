BLD_DIR_REL =../../..
include $(BLD_DIR_REL)/config.mk
include $(BLD_DIR)/util/config.mk
include $(BLD_DIR)/simp/config.mk
include $(SRC_DIR)/simp/patterns.mk
include $(SRC_DIR)/simp/sources.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/simp/tests/interaction/sources.mk

TEST=$(BLD_DIR)/simp/tests/interaction/Test

all: $(simp_tests_interaction_EXES) 

run: $(TEST)
	$(TEST) $(SRC_DIR)/simp/tests/interaction/ > log
	@echo `grep failed log` ", "\
              `grep successful log` "in simp/tests/interaction/log" > count
	@cat count

clean:
	rm -f $(simp_tests_interaction_EXES)
	rm -f $(simp_tests_interaction_OBJS)
	rm -f $(simp_tests_interaction_OBJS:.o=.d) 
	rm -f log count
ifeq ($(BLD_DIR),$(SRC_DIR))
	cd pair; $(MAKE) clean
	cd bond; $(MAKE) clean
	cd angle; $(MAKE) clean
	cd dihedral; $(MAKE) clean
else
	cd $(SRC_DIR)/simp/tests/interaction; $(MAKE) clean-outputs
endif

clean-outputs:
	cd pair; $(MAKE) clean-outputs

-include $(simp_tests_interaction_OBJS:.o=.d)
