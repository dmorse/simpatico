mcMd_mcMoves_common_SRCS=$(SRC_DIR)/mcMd/mcMoves/common/AtomDisplaceMove.cpp \
    $(SRC_DIR)/mcMd/mcMoves/common/DpdMove.cpp \
    $(SRC_DIR)/mcMd/mcMoves/common/HybridMdMove.cpp \
    $(SRC_DIR)/mcMd/mcMoves/common/HybridNphMdMove.cpp \
    $(SRC_DIR)/mcMd/mcMoves/common/MdMove.cpp \
    $(SRC_DIR)/mcMd/mcMoves/common/RigidDisplaceMove.cpp 

mcMd_mcMoves_common_OBJS=$(mcMd_mcMoves_common_SRCS:.cpp=.o)

