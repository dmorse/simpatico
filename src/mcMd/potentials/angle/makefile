BLD_DIR_REL =../../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/mcMd/patterns.mk

all: $(mcMd_potentials_angle_OBJS)

clean:
	rm -f $(mcMd_potentials_angle_OBJS) $(mcMd_potentials_angle_OBJS:.o=.d)

clean-deps:
	rm -f $(mcMd_potentials_angle_OBJS:.o=.d)

-include $(mcMd_potentials_angle_OBJS:.o=.d)

