SRC_DIR_REL =..

include $(SRC_DIR_REL)/compiler.mk
include $(SRC_DIR_REL)/util/defines.mk
include $(SRC_DIR_REL)/util/patterns.mk
include $(SRC_DIR_REL)/util/sources.mk

all: $(util_OBJS) $(util_LIB)

clean:
	rm -f $(util_OBJS) $(util_OBJS:.o=.d) $(util_LIB)
	rm -f Test Test.o Test.d

clean-deps:
	rm -f $(util_OBJS:.o=.d)

clean-test:
	rm -f Test Test.o Test.d

Test: $(util_LIB) $(SRC_DIR)/util/Test.o
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o Test $(SRC_DIR)/util/Test.o \
	       $(util_LIB)

-include $(util_OBJS:.o=.d)
-include Test.d

