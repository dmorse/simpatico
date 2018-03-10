mcMd_tests_links_=mcMd/tests/links/Test.cc

mcMd_tests_links_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mcMd_tests_links_))
mcMd_tests_links_OBJS=\
     $(addprefix $(BLD_DIR)/, $(mcMd_tests_links_:.cc=.o))
mcMd_tests_links_EXES=\
     $(addprefix $(BLD_DIR)/, $(mcMd_tests_links_:.cc=))

