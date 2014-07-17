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
	cd accumulators/unit; $(MAKE) clean
	cd archives; $(MAKE) clean
	cd boundary; $(MAKE) clean
	cd containers; $(MAKE) clean
	cd crystal; $(MAKE) clean
	cd format; $(MAKE) clean
	cd misc; $(MAKE) clean
	cd param; $(MAKE) clean
	cd random; $(MAKE) clean
	cd signal; $(MAKE) clean
	cd space; $(MAKE) clean
	rm -f */*.o */*/*.o 
	rm -f */*.d */*/*.d

-include $(util_tests_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)

