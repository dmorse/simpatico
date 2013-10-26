
mcMd_links_=mcMd/links/LinkMaster.cpp 

mcMd_links_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mcMd_links_))
mcMd_links_OBJS=\
     $(addprefix $(OBJ_DIR)/, $(mcMd_links_:.cpp=.o))

