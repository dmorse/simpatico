BLD_DIR_REL = ../../../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/mcMd/patterns.mk
include $(SRC_DIR)/mcMd/tests/potentials/coulomb/sources.mk

EWALD=$(SRC_DIR)/mcMd/tests/potentials/coulomb/EwaldTest
SPME=$(SRC_DIR)/mcMd/tests/potentials/coulomb/SpmeTest

.Secondary: $(mcMd_tests_potentials_coulomb_OBJS)

all: $(mcMd_tests_potentials_coulomb_EXES)

Ewald: $(EWALD)
	$(EWALD)

Spme: $(SPME)
	$(SPME)

clean:
	rm -f $(mcMd_tests_potentials_coulomb_OBJS)
	rm -f $(mcMd_tests_potentials_coulomb_OBJS:.o=.d)
	rm -f $(mcMd_tests_potentials_coulomb_EXES)
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
