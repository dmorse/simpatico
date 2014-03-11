OBJ_DIR_REL=../../..
include $(OBJ_DIR_REL)/config.mk
include $(OBJ_DIR)/ddMd/config.mk
include $(OBJ_DIR)/inter/config.mk
include $(OBJ_DIR)/util/config.mk
include $(OBJ_DIR)/ddMd/patterns.mk
include $(SRC_DIR)/ddMd/sources.mk
include $(SRC_DIR)/inter/sources.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/ddMd/tests/modifiers/sources.mk

all: $(ddMd_tests_modifiers_OBJS)  Test

clean:
	rm -f $(ddMd_tests_modifiers_OBJS) Test $(ddMd_tests_modifiers_OBJS:.o=.d)

clean-deps:
	rm -f $(ddMd_tests_modifiers_OBJS:.o=.d)

#Test: Test.o $(ddMd_LIB) $(inter_LIB) $(util_LIB)
#	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o Test Test.o \
#	       $(ddMd_LIB) $(inter_LIB) $(util_LIB)

-include $(ddMd_tests_modifiers_OBJS:.o=.d)
-include $(ddMd_OBJS:.o=.d)
-include $(inter_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)

