OBJ_DIR_REL =../..
include $(OBJ_DIR_REL)/compiler.mk
include $(OBJ_DIR)/inter/defines.mk
include $(OBJ_DIR)/util/defines.mk
include $(SRC_DIR)/inter/patterns.mk
include $(SRC_DIR)/inter/sources.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/inter/tests/sources.mk

all: $(inter_tests_OBJS) Test

run: $(inter_tests_OBJS) Test
	./Test $(SRC_DIR)/inter/tests/ > log
	@echo `grep failed log` ", "\
              `grep successful log` "in serial/inter/tests/log" > count
	@cat count

clean:
	rm -f $(inter_tests_OBJS) Test $(inter_tests_OBJS:.o=.d)
	rm -f log count

Test: Test.o $(inter_LIB) $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o Test Test.o \
	       $(inter_LIB) $(util_LIB)

-include $(inter_tests_OBJS:.o=.d)

