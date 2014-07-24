
ddMd_sp_neighbor_= \
    ddMd/sp/neighbor/SpCell.cpp \
    ddMd/sp/neighbor/SpCellList.cpp

# Create lists of source (*.cpp) and object (*.o) files
ddMd_sp_neighbor_SRCS=\
    $(addprefix $(SRC_DIR)/, $(ddMd_sp_neighbor_))
ddMd_sp_neighbor_OBJS=\
    $(addprefix $(BLD_DIR)/, $(ddMd_sp_neighbor_:.cpp=.o))

