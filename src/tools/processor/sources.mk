
tools_processor_= \
    tools/processor/Processor.cpp 

# Create lists of source (*.cpp) and object (*.o) files
tools_processor_SRCS=\
     $(addprefix $(SRC_DIR)/, $(tools_processor_))
tools_processor_OBJS=\
     $(addprefix $(BLD_DIR)/, $(tools_processor_:.cpp=.o))

