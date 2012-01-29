include $(SRC_DIR)/modules/sliplink/mcMoves/sources.mk
include $(SRC_DIR)/modules/sliplink/diagnostics/sources.mk

modules_sliplink_SRCS= \
    $(modules_sliplink_mcMoves_SRCS) \
    $(modules_sliplink_diagnostics_SRCS) \
    $(SRC_DIR)/modules/sliplink/SliplinkMcModule.cpp 

modules_sliplink_OBJS=$(modules_sliplink_SRCS:.cpp=.o)

