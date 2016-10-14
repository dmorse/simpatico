ddMd_analyzers_util_=\
     ddMd/analyzers/util/PairSelector.cpp


ddMd_analyzers_util_SRCS=\
     $(addprefix $(SRC_DIR)/, $(ddMd_analyzers_util_))
ddMd_analyzers_util_OBJS=\
     $(addprefix $(BLD_DIR)/, $(ddMd_analyzers_util_:.cpp=.o))

