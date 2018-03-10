BLD_DIR_REL =../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/mcMd/patterns.mk

all: $(mcMd_mcSimulation_OBJS)

clean:
	rm -f $(mcMd_mcSimulation_OBJS) $(mcMd_mcSimulation_OBJS:.o=.d)

clean-deps:
	rm -f $(mcMd_mcSimulation_OBJS:.o=.d)

-include $(mcMd_mcSimulation_OBJS:.o=.d)

