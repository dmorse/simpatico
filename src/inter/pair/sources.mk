inter_pair_=\
    inter/pair/DpdPair.cpp \
    inter/pair/LJPair.cpp \
    inter/pair/WcaPair.cpp 

inter_pair_SRCS=\
     $(addprefix $(SRC_DIR)/, $(inter_pair_))
inter_pair_OBJS=\
     $(addprefix $(BLD_DIR)/, $(inter_pair_:.cpp=.o))

