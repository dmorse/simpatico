
ddMd_communicate_SRCS=$(SRC_DIR)/ddMd/communicate/Buffer.cpp \
    $(SRC_DIR)/ddMd/communicate/Collector.cpp \
    $(SRC_DIR)/ddMd/communicate/Exchanger.cpp \
    $(SRC_DIR)/ddMd/communicate/Domain.cpp \
    $(SRC_DIR)/ddMd/communicate/AtomDistributor.cpp \
    $(SRC_DIR)/ddMd/communicate/BondDistributor.cpp 

ddMd_communicate_OBJS=$(ddMd_communicate_SRCS:.cpp=.o)

