ddMd_neighbor_= \
    ddMd/neighbor/Cell.cpp \
    ddMd/neighbor/CellList.cpp \
    ddMd/neighbor/PairList.cpp 

ddMd_neighbor_SRCS=\
     $(addprefix $(SRC_DIR)/, $(ddMd_neighbor_))
ddMd_neighbor_OBJS=\
     $(addprefix $(BLD_DIR)/, $(ddMd_neighbor_:.cpp=.o))

