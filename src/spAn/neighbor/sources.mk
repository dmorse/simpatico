
spAn_neighbor_= \
    spAn/neighbor/Cell.cpp \
    spAn/neighbor/CellList.cpp

# Create lists of source (*.cpp) and object (*.o) files
spAn_neighbor_SRCS=\
    $(addprefix $(SRC_DIR)/, $(spAn_neighbor_))
spAn_neighbor_OBJS=\
    $(addprefix $(BLD_DIR)/, $(spAn_neighbor_:.cpp=.o))

