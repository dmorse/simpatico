mcMd_mcMoves_linear_SRCS=\
    $(SRC_DIR)/mcMd/mcMoves/linear/CfbDoubleRebridgeMove.cpp \
    $(SRC_DIR)/mcMd/mcMoves/linear/CfbEndMove.cpp \
    $(SRC_DIR)/mcMd/mcMoves/linear/CfbRebridgeMove.cpp \
    $(SRC_DIR)/mcMd/mcMoves/linear/CfbReptationMove.cpp \
    $(SRC_DIR)/mcMd/mcMoves/linear/EndSwapMove.cpp 

mcMd_mcMoves_linear_OBJS=$(mcMd_mcMoves_linear_SRCS:.cpp=.o)

