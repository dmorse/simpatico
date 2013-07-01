modules_hoomd_mcMoves_SRCS=$(SRC_DIR)/modules/hoomd/mcMoves/HoomdMcMoveFactory.cpp \
    $(SRC_DIR)/modules/hoomd/mcMoves/HoomdMove.cpp  \
    $(SRC_DIR)/modules/hoomd/mcMoves/HoomdNPTMTKMove.cpp

modules_hoomd_mcMoves_OBJS=$(modules_hoomd_mcMoves_SRCS:.cpp=.o)

