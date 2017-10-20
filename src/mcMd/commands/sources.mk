
mcMd_commands_= \
    mcMd/commands/Command.cpp \
    mcMd/commands/McCommandFactory.cpp \
    mcMd/commands/MdCommandFactory.cpp \
    mcMd/commands/CommandManager.cpp \
    mcMd/commands/McDeformCommand.cpp 

mcMd_commands_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mcMd_commands_))
mcMd_commands_OBJS=\
     $(addprefix $(BLD_DIR)/, $(mcMd_commands_:.cpp=.o))

