include $(SRC_DIR)/mdPp/configIos/sources.mk
include $(SRC_DIR)/mdPp/analyzers/sources.mk

mdPp_= \
    $(mdPp_configIos_) \
    mdPp/Processor.cpp
    #$(mdPp_analyzers_) $(mdPp_chemistry_)

# Add user source files, if any
# include $(SRC_DIR)/mdPp/user/sources.mk
# mdPp_+= $(mdPp_user_) 

# Create lists of source (*.cpp) and object (*.o) files
mdPp_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mdPp_))
mdPp_OBJS=\
     $(addprefix $(OBJ_DIR)/, $(mdPp_:.cpp=.o))

$(mdPp_LIB): $(mdPp_OBJS)
	$(AR) rcs $(mdPp_LIB) $(mdPp_OBJS)

