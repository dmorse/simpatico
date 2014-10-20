tools_generators_=\
   tools/generators/chainMaker/ChainMaker.cpp \
   tools/generators/atomicMaker/AtomicMaker.cpp \
   tools/generators/rings/RingMaker.cpp 

tools_generators_SRCS=\
     $(addprefix $(SRC_DIR)/, $(tools_generators_))
tools_generators_OBJS=\
     $(addprefix $(BLD_DIR)/, $(tools_generators_:.cpp=.o))

