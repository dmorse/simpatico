modules_sliplink_mcMoves_=\
    modules/sliplink/mcMoves/SliplinkerAll.cpp \
    modules/sliplink/mcMoves/SliplinkerEnd.cpp \
    modules/sliplink/mcMoves/SliplinkMove.cpp \
    modules/sliplink/mcMoves/GcSliplinkMove.cpp \
    modules/sliplink/mcMoves/SliplinkMcMoveFactory.cpp 

modules_sliplink_mcMoves_SRCS=\
    $(addprefix $(SRC_DIR)/, $(modules_sliplink_mcMoves_))
modules_sliplink_mcMoves_OBJS=\
    $(addprefix $(BLD_DIR)/, $(modules_sliplink_mcMoves_:.cpp=.o))

