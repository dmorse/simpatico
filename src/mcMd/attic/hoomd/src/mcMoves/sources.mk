modules_hoomd_mcMoves_=\
    modules/hoomd/mcMoves/HoomdMcMoveFactory.cpp \
    modules/hoomd/mcMoves/HoomdMove.cpp  \
    modules/hoomd/mcMoves/HoomdNPTMTKMove.cpp

modules_hoomd_mcMoves_SRCS=\
    $(addprefix $(SRC_DIR)/, $(modules_hoomd_mcMoves_))
modules_hoomd_mcMoves_OBJS=\
    $(addprefix $(BLD_DIR)/, $(modules_hoomd_mcMoves_:.cpp=.o))

#modules_hoomd_mcMoves_NVCC_SRCS=\
#    $(addprefix $(SRC_DIR)/, $(modules_hoomd_mcMoves_NVCC_))
#modules_hoomd_mcMoves_NVCC_OBJS=\
#    $(addprefix $(BLD_DIR)/, $(modules_hoomd_mcMoves_NVCC_:.cpp=.o))

