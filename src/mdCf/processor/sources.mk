
mdCf_processor_= \
    mdCf/processor/Processor.cpp 

# Create lists of source (*.cpp) and object (*.o) files
mdCf_processor_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mdCf_processor_))
mdCf_processor_OBJS=\
     $(addprefix $(OBJ_DIR)/, $(mdCf_processor_:.cpp=.o))

