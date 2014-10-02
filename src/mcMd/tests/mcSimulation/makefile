BLD_DIR_REL =../../..
include $(BLD_DIR_REL)/config.mk
include $(BLD_DIR)/mcMd/config.mk
include $(BLD_DIR)/inter/config.mk
include $(BLD_DIR)/util/config.mk
include $(SRC_DIR)/mcMd/patterns.mk
include $(SRC_DIR)/mcMd/sources.mk
include $(SRC_DIR)/inter/sources.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/mcMd/tests/mcSimulation/sources.mk

TEST=mcMd/tests/mcSimulation/Test

all: $(mcMd_tests_mcSimulation_OBJS)

clean:
	rm -f $(mcMd_tests_mcSimulation_OBJS) 
	rm -f $(mcMd_tests_mcSimulation_OBJS:.o=.d)
	rm -f $(mcMd_tests_mcSimulation_OBJS:.o=)
	$(MAKE) clean-outputs

clean-outputs:
	rm -f out/param
	rm -f *.rst

-include $(mcMd_tests_mcSimulation_OBJS:.o=.d)
-include $(mcMd_OBJS:.o=.d)
-include $(inter_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)

