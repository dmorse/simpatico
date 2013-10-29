OBJ_DIR_REL =../..
include $(OBJ_DIR_REL)/compiler.mk
include $(OBJ_DIR)/ddMd/defines.mk
include $(OBJ_DIR)/inter/defines.mk
include $(OBJ_DIR)/util/defines.mk
include $(SRC_DIR)/ddMd/patterns.mk
include $(SRC_DIR)/ddMd/sources.mk
include $(SRC_DIR)/inter/sources.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/ddMd/tests/sources.mk

all: $(ddMd_tests_OBJS) Test

run: $(ddMd_tests_OBJS) Test
	$(MPIRUN) 6 ./Test $(SRC_DIR)/ddMd/tests/ > log
	@echo `grep failed log` ", "\
              `grep successful log` "in serial/ddMd/tests/log" > count
	@cat count

clean:
	-rm -f $(ddMd_tests_OBJS) Test $(ddMd_tests_OBJS:.o=.d)
	-rm -f log count

Test: Test.o $(ddMd_LIB) $(inter_LIB) $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o Test Test.o \
	       $(ddMd_LIB) $(inter_LIB) $(util_LIB)

-include $(ddMd_tests_OBJS:.o=.d)


