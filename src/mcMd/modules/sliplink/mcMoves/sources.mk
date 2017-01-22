mcMd_modules_sliplink_mcMoves_=\
  mcMd/modules/sliplink/mcMoves/SliplinkerAll.cpp \
  mcMd/modules/sliplink/mcMoves/SliplinkerEnd.cpp \
  mcMd/modules/sliplink/mcMoves/SliplinkMove.cpp \
  mcMd/modules/sliplink/mcMoves/GcSliplinkMove.cpp \
  mcMd/modules/sliplink/mcMoves/SliplinkMcMoveFactory.cpp 

mcMd_modules_sliplink_mcMoves_SRCS=\
  $(addprefix $(SRC_DIR)/, $(mcMd_modules_sliplink_mcMoves_))
mcMd_modules_sliplink_mcMoves_OBJS=\
  $(addprefix $(BLD_DIR)/, $(mcMd_modules_sliplink_mcMoves_:.cpp=.o))

