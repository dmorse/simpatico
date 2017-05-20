include $(SRC_DIR)/simp/interaction/sources.mk

simp_=\
    $(simp_interaction_)

simp_SRCS=\
     $(addprefix $(SRC_DIR)/, $(simp_))
simp_OBJS=\
     $(addprefix $(BLD_DIR)/, $(simp_:.cpp=.o))

$(simp_LIB): $(simp_OBJS)
	$(AR) rcs $(simp_LIB) $(simp_OBJS)

