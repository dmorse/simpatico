OBJ_DIR_REL =../../..
include $(OBJ_DIR_REL)/compiler.mk
include $(OBJ_DIR)/mcMd/defines.mk
include $(OBJ_DIR)/inter/defines.mk
include $(OBJ_DIR)/util/defines.mk
include $(SRC_DIR)/inter/patterns.mk
include $(SRC_DIR)/mcMd/sources.mk
include $(SRC_DIR)/inter/sources.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/mcMd/tests/mcSimulation/sources.mk

all: $(mcMd_tests_mcSimulation_OBJS) Test

clean:
	rm -f $(mcMd_tests_mcSimulation_OBJS) Test $(mcMd_tests_mcSimulation_OBJS:.o=.d)
	rm -f out/param

clean-deps:
	rm -f $(mcMd_tests_mcSimulation_OBJS:.o=.d)

Test: Test.o $(mcMd_LIB) $(inter_LIB) $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o Test Test.o \
	       $(mcMd_LIB) $(inter_LIB) $(util_LIB)

-include $(mcMd_tests_mcSimulation_OBJS:.o=.d)

