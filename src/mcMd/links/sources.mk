
mcMd_links_=mcMd/links/LinkMaster.cpp 

mcMd_links_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mcMd_links_))
mcMd_links_OBJS=\
     $(addprefix $(BLD_DIR)/, $(mcMd_links_:.cpp=.o))

