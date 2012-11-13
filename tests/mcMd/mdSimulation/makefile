SRC_DIR_REL =../../../src
include $(SRC_DIR_REL)/compiler.mk
include $(SRC_DIR)/util/defines.mk
include $(SRC_DIR)/inter/defines.mk
include $(SRC_DIR)/mcMd/defines.mk
include $(TESTS_DIR)/mcMd/patterns.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/mcMd/sources.mk
include sources.mk

all: $(tests_mcMd_mdSimulation_OBJS) Test

clean:
	rm -f Test
	rm -f $(tests_mcMd_mdSimulation_OBJS) $(tests_mcMd_mdSimulation_OBJS:.o=.d)
	rm -f *.rst *.cfg

clean-deps:
	rm -f $(tests_mcMd_mdSimulation_OBJS:.o=.d)

Test: Test.o $(mcMd_LIB) $(inter_LIB) $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o Test Test.o \
	       $(mcMd_LIB) $(inter_LIB) $(util_LIB)

-include $(tests_mcMd_mdSimulation_OBJS:.o=.d)
-include Test.d

