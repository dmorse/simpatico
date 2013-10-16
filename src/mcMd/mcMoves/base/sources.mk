mcMd_mcMoves_base_=mcMd/mcMoves/base/CfbEndBase.cpp \
    mcMd/mcMoves/base/CfbRebridgeBase.cpp \
    mcMd/mcMoves/base/GroupRebridgeBase.cpp 

mcMd_mcMoves_base_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mcMd_mcMoves_base_))
mcMd_mcMoves_base_OBJS=\
     $(addprefix $(BLD_DIR)/, $(mcMd_mcMoves_base_:.cpp=.o))

