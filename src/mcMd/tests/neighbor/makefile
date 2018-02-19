BLD_DIR_REL =../../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/mcMd/patterns.mk
include $(SRC_DIR)/mcMd/tests/neighbor/sources.mk

TEST=$(BLD_DIR)/mcMd/tests/neighbor/Test

.Secondary: $(mcMd_tests_neighbor_OBJS)

all: $(mcMd_tests_neighbor_EXES)

run: $(TEST)
	$(TEST)

clean:
	rm -f $(mcMd_tests_neighbor_EXES)
	rm -f $(mcMd_tests_neighbor_OBJS) 
	rm -f $(mcMd_tests_neighbor_OBJS:.o=.d)

-include $(mcMd_tests_neighbor_OBJS:.o=.d)
-include $(mcMd_OBJS:.o=.d)
-include $(simp_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)

