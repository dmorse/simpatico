
mdPp_processor_= \
    mdPp/processor/Processor.cpp \
    mdPp/processor/ProcessorAnalyzerFactory.cpp \
    mdPp/processor/ProcessorAnalyzerManager.cpp 

# Create lists of source (*.cpp) and object (*.o) files
mdPp_processor_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mdPp_processor_))
mdPp_processor_OBJS=\
     $(addprefix $(BLD_DIR)/, $(mdPp_processor_:.cpp=.o))

