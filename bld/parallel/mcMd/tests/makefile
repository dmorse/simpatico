BLD_DIR_REL =../..
include $(BLD_DIR_REL)/config.mk
include $(BLD_DIR)/mcMd/config.mk
include $(BLD_DIR)/inter/config.mk
include $(BLD_DIR)/util/config.mk
include $(SRC_DIR)/mcMd/patterns.mk
include $(SRC_DIR)/mcMd/sources.mk
include $(SRC_DIR)/inter/sources.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/mcMd/tests/sources.mk

TEST=mcMd/tests/Test

all: $(mcMd_tests_OBJS) $(BLD_DIR)/$(TEST)

run: $(mcMd_tests_OBJS) $(BLD_DIR)/$(TEST)
	$(BLD_DIR)/$(TEST) $(SRC_DIR)/mcMd/tests/ > log
	@echo `grep failed log` ", "\
             `grep successful log` "in parallel/mcMd/tests/log" > count
	@cat count

clean:
	-rm -f $(mcMd_tests_OBJS) $(BLD_DIR)/$(TEST) $(mcMd_tests_OBJS:.o=.d)
	-rm -f *.cfg log count

-include $(mcMd_tests_OBJS:.o=.d)
-include $(mcMd_OBJS:.o=.d)
-include $(inter_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)

