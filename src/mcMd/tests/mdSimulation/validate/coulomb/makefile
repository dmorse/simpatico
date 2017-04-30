SRC_DIR_REL = ../../../../..

include $(SRC_DIR_REL)/config.mk
include $(BLD_DIR)/util/config.mk
include $(BLD_DIR)/inter/config.mk
include $(BLD_DIR)/mcMd/config.mk
include $(SRC_DIR)/mcMd/patterns.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/inter/sources.mk
include $(SRC_DIR)/mcMd/sources.mk
include $(SRC_DIR)/mcMd/tests/mdSimulation/validate/coulomb/sources.mk

TEST=mcMd/tests/mcSimulation/validate/EwaldTest

all: $(mcMd_tests_mdSimulation_validate_coulomb_OBJS)

clean:
	rm -f $(mcMd_tests_mdSimulation_validate_coulomb_OBJS)
	rm -f $(mcMd_tests_mdSimulation_validate_coulomb_OBJS:.o=.d)
	rm -f $(mcMd_tests_mdSimulation_validate_coulomb_OBJS:.o=)

echo:
	echo $(INCLUDES)
	echo $(LIBS)
	echo $(mcMd_tests_mdSimulation_validate_coulomb_OBJS)

-include $(mcMd_tests_mdSimulation_validate_coulomb_OBJS:.o=.d)
-include $(mdMd_OBJS:o=.d)
-include $(inter_OBJS:o=.d)
-include $(util_OBJS:o=.d)
