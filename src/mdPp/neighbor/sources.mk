
mdPp_neighbor_= \
    mdPp/neighbor/Cell.cpp \
    mdPp/neighbor/CellList.cpp

# Create lists of source (*.cpp) and object (*.o) files
mdPp_neighbor_SRCS=\
    $(addprefix $(SRC_DIR)/, $(mdPp_neighbor_))
mdPp_neighbor_OBJS=\
    $(addprefix $(BLD_DIR)/, $(mdPp_neighbor_:.cpp=.o))

