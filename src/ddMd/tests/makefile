BLD_DIR_REL =../..
include $(BLD_DIR_REL)/config.mk
include $(BLD_DIR)/ddMd/config.mk
include $(BLD_DIR)/inter/config.mk
include $(BLD_DIR)/util/config.mk
include $(SRC_DIR)/ddMd/patterns.mk
include $(SRC_DIR)/ddMd/sources.mk
include $(SRC_DIR)/inter/sources.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/ddMd/tests/sources.mk

TEST=ddMd/tests/Test

all: $(ddMd_tests_OBJS)

run: $(ddMd_tests_OBJS) $(BLD_DIR)/$(TEST)
	$(MPIRUN) 6 $(BLD_DIR)/$(TEST) $(SRC_DIR)/ddMd/tests/ > log
	@echo `grep failed log` ", "\
              `grep successful log` "in ddMd/tests/log" > count
	@cat count

clean:
	rm -f $(ddMd_tests_OBJS) $(ddMd_tests_OBJS:.o=.d)
	rm -f $(BLD_DIR)/$(TEST)
	cd chemistry; make clean
	cd communicate; make clean
	cd configIos; make clean
	cd neighbor; make clean
	cd potentials; make clean
	cd simulation; make clean
	cd storage; make clean
	make clean-outputs

clean-outputs:
	rm -f log count
	cd chemistry; make clean-outputs
	cd configIos; make clean-outputs
        
-include $(ddMd_tests_OBJS:.o=.d)
-include $(ddMd_OBJS:.o=.d)
-include $(inter_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)

