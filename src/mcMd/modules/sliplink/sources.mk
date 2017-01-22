include $(SRC_DIR)/mcMd/modules/sliplink/mcMoves/sources.mk
include $(SRC_DIR)/mcMd/modules/sliplink/analyzers/sources.mk

mcMd_modules_sliplink_= \
    mcMd/modules/sliplink/SliplinkMcModule.cpp \
    $(mcMd_modules_sliplink_mcMoves_) \
    $(mcMd_modules_sliplink_analyzers_) \

mcMd_modules_sliplink_SRCS=\
         $(addprefix $(SRC_DIR)/, $(mcMd_modules_sliplink_))
mcMd_modules_sliplink_OBJS=\
         $(addprefix $(BLD_DIR)/, $(mcMd_modules_sliplink_:.cpp=.o))

