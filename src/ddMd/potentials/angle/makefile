BLD_DIR_REL =../../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/ddMd/patterns.mk

all: $(ddMd_potentials_angle_OBJS)

clean:
	rm -f $(ddMd_potentials_angle_OBJS) $(ddMd_potentials_angle_OBJS:.o=.d)

clean-deps:
	rm -f $(ddMd_potentials_angle_OBJS:.o=.d)

-include $(ddMd_potentials_angle_OBJS:.o=.d)

