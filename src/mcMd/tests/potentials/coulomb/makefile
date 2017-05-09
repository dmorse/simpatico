SRC_DIR_REL = ../../../..

include $(SRC_DIR_REL)/config.mk
include $(BLD_DIR)/util/config.mk
include $(BLD_DIR)/simp/config.mk
include $(BLD_DIR)/mcMd/config.mk
include $(SRC_DIR)/mcMd/patterns.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/simp/sources.mk
include $(SRC_DIR)/mcMd/sources.mk
include $(SRC_DIR)/mcMd/tests/potentials/coulomb/sources.mk

all: $(mcMd_tests_potentials_coulomb_OBJS)

clean:
	rm -f $(mcMd_tests_potentials_coulomb_OBJS)
	rm -f $(mcMd_tests_potentials_coulomb_OBJS:.o=.d)
	rm -f $(mcMd_tests_potentials_coulomb_OBJS:.o=)
	rm -f *.o
	rm -f *.d

echo:
	echo $(INCLUDES)
	echo $(LIBS)
	echo $(mcMd_tests_potentials_coulomb_OBJS)

-include $(mcMd_tests_potentials_coulomb_OBJS:.o=.d)
-include $(mdMd_OBJS:o=.d)
-include $(simp_OBJS:o=.d)
-include $(util_OBJS:o=.d)
