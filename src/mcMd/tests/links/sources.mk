mcMd_tests_links_=mcMd/links/Test.cc

mcMd_tests_links_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mcMd_tests_links_))
mcMd_tests_links_OBJS=\
     $(addprefix $(OBJ_DIR)/, $(mcMd_tests_links_:.cc=.o))

