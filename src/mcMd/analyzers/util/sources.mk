mcMd_analyzers_util_=mcMd/analyzers/util/PairSelector.cpp 

mcMd_analyzers_util_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mcMd_analyzers_util_))
mcMd_analyzers_util_OBJS=\
     $(addprefix $(OBJ_DIR)/, $(mcMd_analyzers_util_:.cpp=.o))

