
util_accumulators_SRCS=$(SRC_DIR)/util/accumulators/Average.cpp \
    $(SRC_DIR)/util/accumulators/AverageStage.cpp \
    $(SRC_DIR)/util/accumulators/Distribution.cpp \
    $(SRC_DIR)/util/accumulators/IntDistribution.cpp \
    $(SRC_DIR)/util/accumulators/RadialDistribution.cpp 

util_accumulators_OBJS=$(util_accumulators_SRCS:.cpp=.o)

