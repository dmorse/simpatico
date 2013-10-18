BLD_DIR_REL =..

include $(BLD_DIR_REL)/compiler.mk
include $(BLD_DIR)/util/defines.mk
include $(BLD_DIR)/inter/defines.mk
include $(SRC_DIR)/inter/patterns.mk
include $(SRC_DIR)/inter/sources.mk

all: $(inter_OBJS) $(inter_LIB)

clean:
	rm -f $(inter_OBJS) $(inter_LIB)
	rm -f $(inter_OBJS:.o=.d)
	cd tests; make clean
	-rm -f */*.o */*/*.o
	-rm -f */*.d */*/*.d

-include $(util_OBJS:.o=.d)
-include $(inter_OBJS:.o=.d)


