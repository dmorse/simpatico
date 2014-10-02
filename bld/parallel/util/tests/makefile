BLD_DIR_REL =../..
include $(BLD_DIR_REL)/config.mk
include $(BLD_DIR)/util/config.mk
include $(SRC_DIR)/util/patterns.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/util/tests/sources.mk

TEST=util/tests/Test

all: $(util_tests_OBJS) $(BLD_DIR)/$(TEST)

run: $(util_tests_OBJS) $(BLD_DIR)/$(TEST)
	$(BLD_DIR)/$(TEST) $(SRC_DIR)/util/tests/ > log
	@echo `grep failed log` ", "\
              `grep successful log` "in parallel/util/tests/log" > count
	@cat count

clean:
	rm -f $(util_tests_OBJS) $(util_tests_OBJS:.o=.d)
	rm -f $(BLD_DIR)/$(TEST) 
	rm -f dummy log count

$(BLD_DIR)/$(TEST): $(BLD_DIR)/$(TEST).o $(util_LIB)

-include $(util_tests_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)

