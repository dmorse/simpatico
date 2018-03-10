BLD_DIR_REL =../../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/mcMd/patterns.mk
include $(SRC_DIR)/mcMd/tests/links/sources.mk

TEST=$(BLD_DIR)/mcMd/tests/links/Test

.Secondary: $(mcMd_tests_links_OBJS)

all: $(mcMd_tests_links_EXES)

run: $(TEST)
	$(TEST)

clean:
	rm -f $(mcMd_tests_links_EXES) 
	rm -f $(mcMd_tests_links_OBJS) 
	rm -f $(mcMd_tests_links_OBJS:.o=.d)

clean-deps:
	rm -f $(mcMd_tests_links_OBJS:.o=.d)

-include $(mcMd_tests_links_OBJS:.o=.d)
-include $(mcMd_OBJS:.o=.d)

