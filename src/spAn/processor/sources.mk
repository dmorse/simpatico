
spAn_processor_= \
    spAn/processor/Processor.cpp 

# Create lists of source (*.cpp) and object (*.o) files
spAn_processor_SRCS=\
     $(addprefix $(SRC_DIR)/, $(spAn_processor_))
spAn_processor_OBJS=\
     $(addprefix $(BLD_DIR)/, $(spAn_processor_:.cpp=.o))

