OBJ_DIR_REL =../..
include $(OBJ_DIR_REL)/config.mk
include $(OBJ_DIR)/mcMd/config.mk
include $(OBJ_DIR)/inter/config.mk
include $(OBJ_DIR)/util/config.mk
include $(SRC_DIR)/mcMd/patterns.mk
include $(SRC_DIR)/mcMd/sources.mk
include $(SRC_DIR)/inter/sources.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/mcMd/tests/sources.mk

TEST=mcMd/tests/Test

all: $(mcMd_tests_OBJS) $(OBJ_DIR)/$(TEST)

run: $(mcMd_tests_OBJS) $(OBJ_DIR)/$(TEST)
	$(OBJ_DIR)/$(TEST) $(SRC_DIR)/mcMd/tests/ > log
	@echo `grep failed log` ", "\
             `grep successful log` "in mcMd/tests/log" > count
	@cat count

clean:
	rm -f $(mcMd_tests_OBJS) $(mcMd_tests_OBJS:.o=.d)
	rm -f $(OBJ_DIR)/$(TEST) 
	rm -f *.cfg log count
	cd chemistry; make clean
	cd mcSimulation; make clean
	cd mdSimulation; make clean
	cd mdIntegrators; make clean
	cd neighbor; make clean
	cd simulation; make clean
	cd species; make clean

-include $(mcMd_tests_OBJS:.o=.d)
-include $(mcMd_OBJS:.o=.d)
-include $(inter_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)

