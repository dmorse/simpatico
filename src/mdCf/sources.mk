include $(SRC_DIR)/mdCf/chemistry/sources.mk
include $(SRC_DIR)/mdCf/storage/sources.mk
include $(SRC_DIR)/mdCf/configIos/sources.mk
include $(SRC_DIR)/mdCf/analyzers/sources.mk
include $(SRC_DIR)/mdCf/processor/sources.mk

mdCf_= \
    $(mdCf_chemistry_) \
    $(mdCf_storage_) \
    $(mdCf_configIos_) \
    $(mdCf_analyzers_) \
    $(mdCf_processor_) 

# Add user source files, if any
# include $(SRC_DIR)/mdCf/user/sources.mk
# mdCf_+= $(mdCf_user_) 

# Create lists of source (*.cpp) and object (*.o) files
mdCf_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mdCf_))
mdCf_OBJS=\
     $(addprefix $(OBJ_DIR)/, $(mdCf_:.cpp=.o))

$(mdCf_LIB): $(mdCf_OBJS)
	$(AR) rcs $(mdCf_LIB) $(mdCf_OBJS)

