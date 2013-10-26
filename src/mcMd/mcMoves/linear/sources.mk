mcMd_mcMoves_linear_=\
    mcMd/mcMoves/linear/CfbDoubleRebridgeMove.cpp \
    mcMd/mcMoves/linear/CfbEndMove.cpp \
    mcMd/mcMoves/linear/CfbRebridgeMove.cpp \
    mcMd/mcMoves/linear/CfbReptationMove.cpp \
    mcMd/mcMoves/linear/EndSwapMove.cpp 

mcMd_mcMoves_linear_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mcMd_mcMoves_linear_))
mcMd_mcMoves_linear_OBJS=\
     $(addprefix $(OBJ_DIR)/, $(mcMd_mcMoves_linear_:.cpp=.o))

