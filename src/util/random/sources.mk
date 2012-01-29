include $(SRC_DIR)/util/random/mersenne/sources.mk

util_random_SRCS=$(util_random_mersenne_SRCS) \
    $(SRC_DIR)/util/random/Ar1Process.cpp \
    $(SRC_DIR)/util/random/Random.cpp 

util_random_OBJS=$(util_random_SRCS:.cpp=.o)

