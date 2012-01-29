modules_sliplink_mcMoves_SRCS= \
    $(SRC_DIR)/modules/sliplink/mcMoves/SliplinkerAll.cpp \
    $(SRC_DIR)/modules/sliplink/mcMoves/SliplinkerEnd.cpp \
    $(SRC_DIR)/modules/sliplink/mcMoves/SliplinkMove.cpp \
    $(SRC_DIR)/modules/sliplink/mcMoves/GcSliplinkMove.cpp \
    $(SRC_DIR)/modules/sliplink/mcMoves/SliplinkMcMoveFactory.cpp 

modules_sliplink_mcMoves_OBJS=$(modules_sliplink_mcMoves_SRCS:.cpp=.o)

