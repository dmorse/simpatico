include $(SRC_DIR)/mcMd/mcMoves/common/sources.mk
include $(SRC_DIR)/mcMd/mcMoves/semigrand/sources.mk

mcMd_mcMoves_=\
    mcMd/mcMoves/McMove.cpp \
    mcMd/mcMoves/McMoveFactory.cpp \
    mcMd/mcMoves/McMoveManager.cpp \
    mcMd/mcMoves/SystemMove.cpp \
    $(mcMd_mcMoves_common_) \
    $(mcMd_mcMoves_semigrand_) 

ifdef SIMP_BOND
include $(SRC_DIR)/mcMd/mcMoves/base/sources.mk
include $(SRC_DIR)/mcMd/mcMoves/linear/sources.mk
include $(SRC_DIR)/mcMd/mcMoves/ring/sources.mk
mcMd_mcMoves_+=\
    $(mcMd_mcMoves_base_) \
    $(mcMd_mcMoves_linear_) \
    $(mcMd_mcMoves_ring_) 
endif

mcMd_mcMoves_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mcMd_mcMoves_))
mcMd_mcMoves_OBJS=\
     $(addprefix $(BLD_DIR)/, $(mcMd_mcMoves_:.cpp=.o))

