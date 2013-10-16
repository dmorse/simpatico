
mcMd_neighbor_=mcMd/neighbor/Cell.cpp \
    mcMd/neighbor/CellList.cpp \
    mcMd/neighbor/PairList.cpp 

mcMd_neighbor_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mcMd_neighbor_))
mcMd_neighbor_OBJS=\
     $(addprefix $(BLD_DIR)/, $(mcMd_neighbor_:.cpp=.o))

