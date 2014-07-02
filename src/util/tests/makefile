BLD_DIR_REL =../..
include $(BLD_DIR_REL)/config.mk
include $(BLD_DIR)/util/config.mk
include $(SRC_DIR)/util/patterns.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/util/tests/sources.mk

TEST=util/tests/Test

all: $(util_tests_OBJS)

run: $(util_tests_OBJS) $(BLD_DIR)/$(TEST)
	$(BLD_DIR)/$(TEST) $(SRC_DIR)/util/tests/ > log
	@echo `grep failed log` ", "\
              `grep successful log` "in util/tests/log" > count
	@cat count

clean:
	rm -f $(util_tests_OBJS) $(util_tests_OBJS:.o=.d)
	rm -f $(BLD_DIR)/$(TEST) $(BLD_DIR)/$(TEST).d
	rm -f dummy log count
	cd accumulators/unit; make clean
	cd archives; make clean
	cd boundary; make clean
	cd containers; make clean
	cd crystal; make clean
	cd format; make clean
	cd misc; make clean
	cd param; make clean
	cd random; make clean
	cd signal; make clean
	cd space; make clean
	rm -f */*.o */*/*.o 
	rm -f */*.d */*/*.d

-include $(util_tests_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)

