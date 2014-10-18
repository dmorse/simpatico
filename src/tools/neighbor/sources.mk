
tools_neighbor_= \
    tools/neighbor/Cell.cpp \
    tools/neighbor/CellList.cpp

# Create lists of source (*.cpp) and object (*.o) files
tools_neighbor_SRCS=\
    $(addprefix $(SRC_DIR)/, $(tools_neighbor_))
tools_neighbor_OBJS=\
    $(addprefix $(BLD_DIR)/, $(tools_neighbor_:.cpp=.o))

