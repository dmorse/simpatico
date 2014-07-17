BLD_DIR_REL =../../..
include $(BLD_DIR_REL)/config.mk
include $(BLD_DIR)/ddMd/config.mk
include $(BLD_DIR)/inter/config.mk
include $(BLD_DIR)/util/config.mk
include $(SRC_DIR)/ddMd/patterns.mk
include $(SRC_DIR)/inter/sources.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/ddMd/sources.mk
include $(SRC_DIR)/ddMd/tests/sp/sources.mk

all: $(ddMd_tests_sp_OBJS)

run:
	./Test
	cd configIos; ./Test

clean:
	cd chemistry; $(MAKE) clean
	cd storage; $(MAKE) clean
	cd configIos; $(MAKE) clean
	cd processor; $(MAKE) clean
	rm -f Test Test.d Test.o

-include $(ddMd_tests_sp_OBJS:.o=.d)
-include $(ddMd_sp_OBJS:.o=.d)
-include $(inter_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)

