BLD_DIR_REL =../..
include $(BLD_DIR_REL)/config.mk
include $(BLD_DIR)/spAn/config.mk
include $(BLD_DIR)/inter/config.mk
include $(BLD_DIR)/util/config.mk
include $(SRC_DIR)/spAn/patterns.mk
include $(SRC_DIR)/inter/sources.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/spAn/sources.mk
include $(SRC_DIR)/spAn/tests/sources.mk

all: $(spAn_tests_OBJS)

run:
	./Test
	cd configIos; ./Test

clean:
	cd chemistry; $(MAKE) clean
	cd storage; $(MAKE) clean
	cd configIos; $(MAKE) clean
	cd processor; $(MAKE) clean
	cd neighbor; $(MAKE) clean
	rm -f Test Test.d Test.o

-include $(spAn_tests_OBJS:.o=.d)
-include $(spAn_OBJS:.o=.d)
-include $(inter_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)

