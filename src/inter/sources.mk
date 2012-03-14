include $(SRC_DIR)/inter/pair/sources.mk
include $(SRC_DIR)/inter/bond/sources.mk

inter_SRCS=\
    $(inter_pair_SRCS)\
    $(inter_bond_SRCS)

ifdef INTER_ANGLE
include $(SRC_DIR)/inter/angle/sources.mk
inter_SRCS+=$(inter_angle_SRCS)
endif

ifdef INTER_DIHEDRAL
include $(SRC_DIR)/inter/dihedral/sources.mk
inter_SRCS+=$(inter_dihedral_SRCS)
endif

ifdef INTER_EXTERNAL
include $(SRC_DIR)/inter/external/sources.mk
inter_SRCS+=$(inter_external_SRCS)
endif

ifdef INTER_TETHER
include $(SRC_DIR)/inter/tether/sources.mk
inter_SRCS+=$(inter_tether_SRCS)
endif

inter_OBJS=$(inter_SRCS:.cpp=.o)

$(inter_LIB): $(inter_OBJS)
	$(AR) rcs $(inter_LIB) $(inter_OBJS)

