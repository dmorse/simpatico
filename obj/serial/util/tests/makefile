OBJ_DIR_REL =../..
include $(OBJ_DIR_REL)/compiler.mk
include $(OBJ_DIR)/util/defines.mk
include $(SRC_DIR)/util/patterns.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/util/tests/sources.mk

all: $(util_tests_OBJS) Test

run: $(util_tests_OBJS) Test
	./Test $(SRC_DIR)/util/tests/ > log
	@echo `grep failed log` ", "\
              `grep successful log` "in serial/util/tests/log" > count
	@cat count

clean:
	rm -f $(util_tests_OBJS) Test $(util_tests_OBJS:.o=.d)
	-rm -f dummy log count

clean-deps:
	rm -f $(util_tests_OBJS:.o=.d)

Test: Test.o $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o Test Test.o \
	       $(util_LIB)

-include $(util_tests_OBJS:.o=.d)

