# Include source files lists from subdirectories
include $(SRC_DIR)/simp/interaction/sources.mk
include $(SRC_DIR)/simp/species/sources.mk

# Concatenate source file lists from subdirectories
simp_=\
    $(simp_interaction_) \
    $(simp_species_)

# Create lists of src and object files, with absolute paths
simp_SRCS=\
     $(addprefix $(SRC_DIR)/, $(simp_))
simp_OBJS=\
     $(addprefix $(BLD_DIR)/, $(simp_:.cpp=.o))

# Target to create library file for Simp namespace
$(simp_LIB): $(simp_OBJS)
	$(AR) rcs $(simp_LIB) $(simp_OBJS)

