include $(SRC_DIR)/inter/user/sources.mk

inter_SRCS=\
    $(inter_user_SRCS)

ifndef INTER_NOPAIR
include $(SRC_DIR)/inter/pair/sources.mk
inter_+=$(inter_pair_)
endif

ifdef INTER_BOND
include $(SRC_DIR)/inter/bond/sources.mk
inter_+=$(inter_bond_)
endif

ifdef INTER_ANGLE
include $(SRC_DIR)/inter/angle/sources.mk
inter_+=$(inter_angle_)
endif

ifdef INTER_DIHEDRAL
include $(SRC_DIR)/inter/dihedral/sources.mk
inter_+=$(inter_dihedral_)
endif

ifdef INTER_EXTERNAL
include $(SRC_DIR)/inter/external/sources.mk
inter_+=$(inter_external_)
endif

inter_SRCS=\
     $(addprefix $(SRC_DIR)/, $(inter_))
inter_OBJS=\
     $(addprefix $(BLD_DIR)/, $(inter_:.cpp=.o))


$(inter_LIB): $(inter_OBJS)
	$(AR) rcs $(inter_LIB) $(inter_OBJS)

