mcMd_mcMoves_common_=\
    mcMd/mcMoves/common/AtomDisplaceMove.cpp \
    mcMd/mcMoves/common/DpdMove.cpp \
    mcMd/mcMoves/common/HybridMdMove.cpp \
    mcMd/mcMoves/common/HybridNphMdMove.cpp \
    mcMd/mcMoves/common/MdMove.cpp \
    mcMd/mcMoves/common/RigidDisplaceMove.cpp 

mcMd_mcMoves_common_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mcMd_mcMoves_common_))
mcMd_mcMoves_common_OBJS=\
     $(addprefix $(BLD_DIR)/, $(mcMd_mcMoves_common_:.cpp=.o))

