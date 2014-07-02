mcMd_mcMoves_ring_=\
    mcMd/mcMoves/ring/CfbRingRebridgeMove.cpp \
    mcMd/mcMoves/ring/RingOctaRebridgeMove.cpp \
    mcMd/mcMoves/ring/RingTetraRebridgeMove.cpp 

mcMd_mcMoves_ring_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mcMd_mcMoves_ring_))
mcMd_mcMoves_ring_OBJS=\
     $(addprefix $(BLD_DIR)/, $(mcMd_mcMoves_ring_:.cpp=.o))

