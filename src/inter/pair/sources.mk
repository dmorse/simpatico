inter_pair_SRCS=\
    $(SRC_DIR)/inter/pair/DpdPair.cpp \
    $(SRC_DIR)/inter/pair/LJPair.cpp \
    $(SRC_DIR)/inter/pair/WcaPair.cpp 

inter_pair_OBJS=$(inter_pair_SRCS:.cpp=.o)

