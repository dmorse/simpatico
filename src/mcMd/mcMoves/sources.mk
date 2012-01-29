include $(SRC_DIR)/mcMd/mcMoves/common/sources.mk
include $(SRC_DIR)/mcMd/mcMoves/base/sources.mk
include $(SRC_DIR)/mcMd/mcMoves/linear/sources.mk
include $(SRC_DIR)/mcMd/mcMoves/ring/sources.mk
include $(SRC_DIR)/mcMd/mcMoves/semigrand/sources.mk

mcMd_mcMoves_SRCS=$(mcMd_mcMoves_common_SRCS) $(mcMd_mcMoves_base_SRCS) \
    $(mcMd_mcMoves_linear_SRCS) $(mcMd_mcMoves_ring_SRCS) \
    $(mcMd_mcMoves_semigrand_SRCS) $(SRC_DIR)/mcMd/mcMoves/McMove.cpp \
    $(SRC_DIR)/mcMd/mcMoves/McMoveFactory.cpp \
    $(SRC_DIR)/mcMd/mcMoves/McMoveManager.cpp \
    $(SRC_DIR)/mcMd/mcMoves/SystemMove.cpp 

mcMd_mcMoves_OBJS=$(mcMd_mcMoves_SRCS:.cpp=.o)

