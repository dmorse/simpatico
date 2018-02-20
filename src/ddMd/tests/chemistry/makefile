BLD_DIR_REL =../../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/ddMd/patterns.mk
include $(SRC_DIR)/ddMd/tests/chemistry/sources.mk

TEST= $(SRC_DIR)/ddMd/tests/chemistry/Test

.Secondary: $(ddMd_tests_chemistry_OBJS) 

all: $(ddMd_tests_chemistry_EXES) 

run: $(TEST)
	$(TEST)

clean:
	rm -f $(ddMd_tests_chemistry_EXES)
	rm -f $(ddMd_tests_chemistry_OBJS) 
	rm -f $(ddMd_tests_chemistry_OBJS:.o=.d)
	make clean-outputs

clean-outputs:
	@rm -f tmp/*

-include $(ddMd_tests_chemistry_OBJS:.o=.d)
-include $(ddMd_OBJS:.o=.d)
-include $(simp_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)

