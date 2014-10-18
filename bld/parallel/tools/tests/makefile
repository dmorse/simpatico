BLD_DIR_REL =../..
include $(BLD_DIR_REL)/config.mk
include $(BLD_DIR)/tools/config.mk
include $(BLD_DIR)/inter/config.mk
include $(BLD_DIR)/util/config.mk
include $(SRC_DIR)/tools/patterns.mk
include $(SRC_DIR)/inter/sources.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/tools/sources.mk
include $(SRC_DIR)/tools/tests/sources.mk


TEST=tools/tests/Test

all: $(tools_tests_OBJS)

run: $(tools_tests_OBJS) $(BLD_DIR)/$(TEST)
	$(BLD_DIR)/$(TEST) $(SRC_DIR)/tools/tests/ > log
	@echo `grep failed log` ", "\
              `grep successful log` "in tools/tests/log" > count
	@cat count

clean:
	rm -f $(tools_tests_OBJS) $(tools_tests_OBJS:.o=.d)
	rm -f $(BLD_DIR)/$(TEST)
	$(MAKE) clean-outputs

clean-outputs:
	rm -f log count
        
-include $(tools_tests_OBJS:.o=.d)
-include $(tools_OBJS:.o=.d)
-include $(inter_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)


