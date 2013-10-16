mcMd_diagnostics_util_=mcMd/diagnostics/util/PairSelector.cpp 

mcMd_diagnostics_util_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mcMd_diagnostics_util_))
mcMd_diagnostics_util_OBJS=\
     $(addprefix $(BLD_DIR)/, $(mcMd_diagnostics_util_:.cpp=.o))

