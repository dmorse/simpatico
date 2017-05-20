BLD_DIR_REL =../../..
include $(BLD_DIR_REL)/config.mk
include $(BLD_DIR)/ddMd/config.mk
include $(BLD_DIR)/simp/config.mk
include $(BLD_DIR)/util/config.mk
include $(SRC_DIR)/ddMd/patterns.mk
include $(SRC_DIR)/ddMd/sources.mk
include $(SRC_DIR)/simp/sources.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/ddMd/tests/neighbor/sources.mk

all: $(ddMd_tests_neighbor_OBJS)

clean:
	rm -f $(ddMd_tests_neighbor_OBJS) 
	rm -f $(ddMd_tests_neighbor_OBJS:.o=.d)
	rm -f $(ddMd_tests_neighbor_OBJS:.o=)

-include $(ddMd_tests_neighbor_OBJS:.o=.d)
-include $(ddMd_OBJS:.o=.d)
-include $(simp_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)

