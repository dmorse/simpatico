BLD_DIR_REL =../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/mcMd/patterns.mk

all: $(mcMd_mdSimulation_OBJS)

clean:
	rm -f $(mcMd_mdSimulation_OBJS) $(mcMd_mdSimulation_OBJS:.o=.d)

clean-deps:
	rm -f $(mcMd_mdSimulation_OBJS:.o=.d)

-include $(mcMd_mdSimulation_OBJS:.o=.d)

