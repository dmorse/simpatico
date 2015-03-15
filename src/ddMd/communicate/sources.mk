ddMd_communicate_=\
    ddMd/communicate/Buffer.cpp \
    ddMd/communicate/Domain.cpp \
    ddMd/communicate/AtomDistributor.cpp \
    ddMd/communicate/Exchanger.cpp \
    ddMd/communicate/AtomCollector.cpp \
    ddMd/communicate/Plan.cpp 

ddMd_communicate_SRCS=\
     $(addprefix $(SRC_DIR)/, $(ddMd_communicate_))
ddMd_communicate_OBJS=\
     $(addprefix $(BLD_DIR)/, $(ddMd_communicate_:.cpp=.o))

