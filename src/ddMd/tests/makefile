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

TEST=ddMd/tests/Test

all: $(ddMd_tests_OBJS) $(OBJ_DIR)/$(TEST)

run: $(ddMd_tests_OBJS) $(OBJ_DIR)/$(TEST)
	$(MPIRUN) 6 ./Test $(SRC_DIR)/ddMd/tests/ > log
	@echo `grep failed log` ", "\
              `grep successful log` "in src/ddMd/tests/log" > count
	@cat count

clean:
	rm -f $(ddMd_tests_OBJS) $(ddMd_tests_OBJS:.o=.d)
	rm -f $(OBJ_DIR)/$(TEST)
	rm -f log count

	# cd communicate; make clean
	#cd configIo; make clean
        
-include $(ddMd_tests_OBJS:.o=.d)

